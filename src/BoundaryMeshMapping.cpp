#include <clot/BoundaryMeshMapping.h>
#include <clot/app_namespaces.h>

#include <ibamr/IBFEMethod.h>

#include <ibtk/IBTK_MPI.h>
#include <ibtk/IndexUtilities.h>
#include <ibtk/LEInteractor.h>
#include <ibtk/libmesh_utilities.h>

#include <libmesh/boundary_info.h>
#include <libmesh/enum_xdr_mode.h>
#include <libmesh/explicit_system.h>

namespace clot
{
BoundaryMeshMapping::BoundaryMeshMapping(std::string object_name,
                                         Pointer<Database> input_db,
                                         MeshBase* mesh,
                                         FEDataManager* data_manager,
                                         const std::string& restart_read_dirname,
                                         unsigned int restart_restore_num)
    : GeneralBoundaryMeshMapping(std::move(object_name), input_db, mesh, restart_read_dirname, restart_restore_num),
      d_base_data_manager(data_manager)
{
    // intentionally blank
}

void
BoundaryMeshMapping::buildBoundaryMesh()
{
    // TODO: We need some way to extract a boundary mesh. For now, because we only have one part, we can extract the
    // part we want.
    auto bdry_mesh = libmesh_make_unique<BoundaryMesh>(*static_cast<BoundaryMesh*>(d_base_meshes[0]));

    // Now generate a mapping between nodes
    // TODO: find a more efficient way to do this.
    auto n_it_end = d_base_meshes[0]->nodes_end();
    for (auto it = d_base_meshes[0]->nodes_begin(); it != n_it_end; ++it)
    {
        Node* orig_n = *it;
        auto bdry_n_it_end = bdry_mesh->nodes_end();
        for (auto bdry_n_it = bdry_mesh->nodes_begin(); bdry_n_it != bdry_n_it_end; ++bdry_n_it)
        {
            Node* bdry_n = *bdry_n_it;
            bool found_match = true;
            for (unsigned int d = 0; d < bdry_mesh->spatial_dimension(); ++d)
                found_match = found_match && (*bdry_n)(d) == (*orig_n)(d);
            if (found_match)
            {
                d_bdry_base_node_map[bdry_n] = orig_n;
                // We can quit this loop early.
                break;
            }
        }
    }

    d_bdry_meshes.push_back(std::move(bdry_mesh));
}

void
BoundaryMeshMapping::updateBoundaryLocation(const double time, const unsigned int part, const bool /*end_of_timestep*/)
{
    // Move the boundary mesh to match the volume mesh.
    EquationSystems* eq_sys = d_base_data_manager->getEquationSystems();

    System& X_base_system = eq_sys->get_system(d_base_data_manager->COORDINATES_SYSTEM_NAME);
    const DofMap& X_base_dof_map = X_base_system.get_dof_map();
    NumericVector<double>* X_base_vec = X_base_system.solution.get();

    System& X_bdry_sys = d_bdry_eq_sys_vec[0]->get_system(d_coords_sys_name);
    const DofMap& X_bdry_dof_map = X_bdry_sys.get_dof_map();
    NumericVector<double>* X_bdry_vec = X_bdry_sys.solution.get();

    System& dX_bdry_sys = d_bdry_eq_sys_vec[0]->get_system(d_disp_sys_name);
    NumericVector<double>* dX_bdry_vec = dX_bdry_sys.solution.get();

    auto node_it = d_bdry_meshes[0]->local_nodes_begin();
    auto node_end = d_bdry_meshes[0]->local_nodes_end();
    for (; node_it != node_end; ++node_it)
    {
        Node* node = *node_it;
        Node* base_node = d_bdry_base_node_map.at(node);
        // Grab current position of volumetric mesh.
        std::vector<dof_id_type> X_dof_indices, X_bdry_dof_indices;
        for (int d = 0; d < NDIM; ++d)
        {
            X_base_dof_map.dof_indices(base_node, X_dof_indices, d);
            X_bdry_dof_map.dof_indices(node, X_bdry_dof_indices, d);
            X_bdry_vec->set(X_bdry_dof_indices[0], (*X_base_vec)(X_dof_indices[0]));
            dX_bdry_vec->set(X_bdry_dof_indices[0], (*X_base_vec)(X_dof_indices[0]) - (*node)(d));
        }
    }
    X_bdry_vec->close();
    dX_bdry_vec->close();
    X_bdry_sys.update();
    dX_bdry_sys.update();
    return;
}
} // namespace clot
