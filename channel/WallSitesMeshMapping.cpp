#include <clot/app_namespaces.h>

#include <ibamr/IBFEMethod.h>

#include <ibtk/IBTK_MPI.h>
#include <ibtk/IndexUtilities.h>
#include <ibtk/LEInteractor.h>
#include <ibtk/libmesh_utilities.h>

#include <libmesh/boundary_info.h>
#include <libmesh/enum_xdr_mode.h>
#include <libmesh/explicit_system.h>

// Local includes
#include "WallSitesMeshMapping.h"

namespace clot
{
WallSitesMeshMapping::WallSitesMeshMapping(std::string object_name,
                                           Pointer<Database> input_db,
                                           std::string W_sys_name,
                                           MeshBase* mesh,
                                           const std::string& restart_read_dirname,
                                           unsigned int restart_restore_num)
    : GeneralBoundaryMeshMapping(std::move(object_name), input_db, mesh, restart_read_dirname, restart_restore_num),
      d_W_sys_name(std::move(W_sys_name))
{
    auto var_db = VariableDatabase<NDIM>::getDatabase();
    d_w_var = new CellVariable<NDIM, double>(d_object_name + "::wall_sites");
    d_w_idx = var_db->registerVariableAndContext(d_w_var, var_db->getContext(d_object_name + "::ctx"), 4);
}

WallSitesMeshMapping::~WallSitesMeshMapping()
{
    // Deallocate the wall sites
    Pointer<PatchHierarchy<NDIM>> hierarchy = d_bdry_mesh_partitioners[0]->getPatchHierarchy();
    for (int ln = 0; ln <= hierarchy->getFinestLevelNumber(); ++ln)
    {
        Pointer<PatchLevel<NDIM>> level = hierarchy->getPatchLevel(ln);
        if (level->checkAllocated(d_w_idx)) level->deallocatePatchData(d_w_idx);
    }
}

void
WallSitesMeshMapping::initializeFEData(Pointer<PatchHierarchy<NDIM>> hierarchy)
{
    const bool from_restart = RestartManager::getManager()->isFromRestart();
    if (from_restart) return;
    for (unsigned int part = 0; part < d_bdry_meshes.size(); ++part)
    {
        d_bdry_eq_sys_vec[part]->init();
    }
    d_bdry_mesh_partitioners[0]->setPatchHierarchy(hierarchy);
    updateBoundaryLocation(0.0, false);
}

void
WallSitesMeshMapping::buildBoundaryMesh()
{
    // TODO: We need some way to extract a boundary mesh. For now, because we only have one part, we can extract the
    // part we want.
    auto bdry_mesh = libmesh_make_unique<BoundaryMesh>(*static_cast<BoundaryMesh*>(d_base_meshes[0]));
    d_bdry_meshes.push_back(std::move(bdry_mesh));
}

void
WallSitesMeshMapping::spreadWallSites(int w_idx)
{
    Pointer<PatchHierarchy<NDIM>> hierarchy = d_bdry_mesh_partitioners[0]->getPatchHierarchy();
    if (w_idx == IBTK::invalid_index) w_idx = d_w_idx;
    plog << "Spreading wall sites to the fluid.\n";

    d_bdry_mesh_partitioners[0]->reinitElementMappings();
    System& W_sys = d_bdry_mesh_partitioners[0]->getEquationSystems()->get_system(d_W_sys_name);
    const DofMap& W_dof_map = W_sys.get_dof_map();
    NumericVector<double>* W_vec = d_bdry_mesh_partitioners[0]->buildGhostedSolutionVector(d_W_sys_name);
    System& X_sys = d_bdry_mesh_partitioners[0]->getEquationSystems()->get_system(
        d_bdry_mesh_partitioners[0]->COORDINATES_SYSTEM_NAME);
    const DofMap& X_dof_map = X_sys.get_dof_map();
    NumericVector<double>* X_vec =
        d_bdry_mesh_partitioners[0]->buildGhostedSolutionVector(d_bdry_mesh_partitioners[0]->COORDINATES_SYSTEM_NAME);
    // Multiply by the nodal volume fractions (to convert densities intovalues).
    PetscVector<double>* dX_vec = d_bdry_mesh_partitioners[0]->buildIBGhostedDiagonalL2MassMatrix(d_W_sys_name);
    std::unique_ptr<NumericVector<double>> W_x_dX_vec = W_vec->clone();
    W_x_dX_vec->pointwise_mult(*W_vec, *dX_vec);

    // Extract local form vectors.
    auto W_x_dX_petsc_vec = static_cast<PetscVector<double>*>(W_x_dX_vec.get());
    const double* const W_x_dX_local_soln = W_x_dX_petsc_vec->get_array_read();
    auto X_petsc_vec = static_cast<PetscVector<double>*>(X_vec);
    const double* const X_local_soln = X_petsc_vec->get_array_read();

    for (int ln = 0; ln <= hierarchy->getFinestLevelNumber(); ++ln)
    {
        Pointer<PatchLevel<NDIM>> level = hierarchy->getPatchLevel(ln);
        if (!level->checkAllocated(w_idx)) level->allocatePatchData(w_idx);
        int local_patch_num = 0;
        const std::vector<std::vector<Node*>>& active_patch_node_map =
            d_bdry_mesh_partitioners[0]->getActivePatchNodeMap(ln);
        // Spread from the nodes.
        for (PatchLevel<NDIM>::Iterator p(level); p; p++, ++local_patch_num)
        {
            const Pointer<Patch<NDIM>> patch = level->getPatch(p());
            Pointer<CellData<NDIM, double>> w_data = patch->getPatchData(w_idx);
            w_data->fillAll(0.0);
            // The relevant collection of nodes.
            const std::vector<Node*>& patch_nodes = active_patch_node_map[local_patch_num];
            const size_t num_active_patch_nodes = patch_nodes.size();
            if (!num_active_patch_nodes) continue;

            const Pointer<CartesianPatchGeometry<NDIM>> patch_geom = patch->getPatchGeometry();
            const double* const patch_x_lower = patch_geom->getXLower();
            const double* const patch_x_upper = patch_geom->getXUpper();
            std::array<bool, NDIM> touches_upper_regular_bdry;
            for (unsigned int d = 0; d < NDIM; ++d)
                touches_upper_regular_bdry[d] = patch_geom->getTouchesRegularBoundary(d, 1);

            // Store the values of F_JxW and X at the nodes inside the patch.
            std::vector<double> W_x_dX_node, X_node;
            W_x_dX_node.reserve(num_active_patch_nodes);
            X_node.reserve(NDIM * num_active_patch_nodes);
            std::vector<dof_id_type> F_idxs, X_idxs;
            IBTK::Point X;
            for (unsigned int k = 0; k < num_active_patch_nodes; ++k)
            {
                const Node* const n = patch_nodes[k];
                bool inside_patch = true;
                for (unsigned int d = 0; d < NDIM; ++d)
                {
                    IBTK::get_nodal_dof_indices(X_dof_map, n, d, X_idxs);
                    X[d] = X_local_soln[X_petsc_vec->map_global_to_local_index(X_idxs[0])];
                    inside_patch =
                        inside_patch && (X[d] >= patch_x_lower[d]) &&
                        ((X[d] < patch_x_upper[d]) || (touches_upper_regular_bdry[d] && X[d] <= patch_x_upper[d]));
                }
                if (inside_patch)
                {
                    IBTK::get_nodal_dof_indices(W_dof_map, n, 0, F_idxs);
                    for (const auto& F_idx : F_idxs)
                    {
                        W_x_dX_node.push_back(W_x_dX_local_soln[W_x_dX_petsc_vec->map_global_to_local_index(F_idx)]);
                    }
                    X_node.insert(X_node.end(), &X[0], &X[0] + NDIM);
                }
            }
            TBOX_ASSERT(W_x_dX_node.size() <= num_active_patch_nodes);
            TBOX_ASSERT(X_node.size() <= NDIM * num_active_patch_nodes);

            // Spread values from the nodes to the Cartesian grid patch.
            //
            // \todo Fix this for FE structures with periodic boundaries.  Nodes
            // on periodic boundaries will be "double counted".
            //
            // \todo Add warnings for FE structures with periodic boundaries.
            const Box<NDIM> spread_box = patch->getBox();
            LEInteractor::spread(w_data, W_x_dX_node, 1, X_node, NDIM, patch, spread_box, d_kernel_fcn);
        }
    }

    // Restore local form vectors.
    W_x_dX_petsc_vec->restore_array();
    X_petsc_vec->restore_array();
}
} // namespace clot
