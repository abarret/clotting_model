#include "ibamr/app_namespaces.h"
#include <ibamr/IBFEMethod.h>

#include "ibtk/IBTK_MPI.h"
#include "ibtk/IndexUtilities.h"
#include "ibtk/libmesh_utilities.h"
#include <ibtk/LEInteractor.h>

#include "BoundaryMeshMapping.h"

#include "libmesh/enum_xdr_mode.h"
#include "libmesh/explicit_system.h"
#include <libmesh/boundary_info.h>

namespace IBAMR
{
BoundaryMeshMapping::BoundaryMeshMapping(std::string object_name,
                                         Pointer<Database> input_db,
                                         MeshBase* mesh,
                                         FEDataManager* data_manager)
    : d_object_name(std::move(object_name)), d_orig_mesh(mesh), d_orig_data_manager(data_manager), d_input_db(input_db)
{
    auto var_db = VariableDatabase<NDIM>::getDatabase();
    d_w_var = new CellVariable<NDIM, double>("wall_sites");
    d_w_idx = var_db->registerVariableAndContext(d_w_var, var_db->getContext("wall_sites"), 4);

    // TODO: We need some way to extract a boundary mesh. For now, because we only have one part, we can extract the
    // part we want.
    d_bdry_mesh = libmesh_make_unique<BoundaryMesh>(*(static_cast<BoundaryMesh*>(mesh)));
    // Now delete all the nodes of the bdry_mesh that are below the midline
    auto e_it_end = d_bdry_mesh->elements_end();
    for (auto it = d_bdry_mesh->elements_begin(); it != e_it_end; ++it)
    {
        Elem* e = *it;
        // If this element has nodes on the bottom half of the cylinder, delete it.
        for (unsigned int n = 0; n < e->n_nodes(); ++n)
        {
            if (e->node_ref(n)(1) < 0.5)
            {
                d_bdry_mesh->delete_elem(e);
                break;
            }
        }
    }
    auto n_it_end = d_bdry_mesh->nodes_end();
    for (auto it = d_bdry_mesh->nodes_begin(); it != n_it_end; ++it)
    {
        Node* n = *it;
        // If the node is on the bottom half of the cylinder, delete it.
        if ((*n)(1) < 0.5) d_bdry_mesh->delete_node(n);
    }
    d_bdry_mesh->prepare_for_use();
    d_bdry_mesh->print_info();

    // Now generate a mapping between nodes
    // TODO: find a more efficient way to do this.
    n_it_end = d_orig_mesh->nodes_end();
    for (auto it = d_orig_mesh->nodes_begin(); it != n_it_end; ++it)
    {
        Node* orig_n = *it;
        auto bdry_n_it_end = d_bdry_mesh->nodes_end();
        for (auto bdry_n_it = d_bdry_mesh->nodes_begin(); bdry_n_it != bdry_n_it_end; ++bdry_n_it)
        {
            Node* bdry_n = *bdry_n_it;
            bool found_match = true;
            for (int d = 0; d < d_bdry_mesh->mesh_dimension(); ++d)
                found_match = found_match && (*bdry_n)(d) == (*orig_n)(d);
            if (found_match) d_bdry_orig_node_map[bdry_n] = orig_n;
        }
    }

    // Now make equation systems
    d_bdry_eq_sys = std::move(libmesh_make_unique<EquationSystems>(*d_bdry_mesh));
    d_fe_data = std::make_shared<FEData>(d_object_name + "::FEData", *d_bdry_eq_sys, true);

    auto& X_sys = d_bdry_eq_sys->add_system<ExplicitSystem>(d_coords_sys_name);
    for (unsigned int d = 0; d < NDIM; ++d) X_sys.add_variable("X_" + std::to_string(d));
    auto& dX_sys = d_bdry_eq_sys->add_system<ExplicitSystem>(d_disp_sys_name);
    for (unsigned int d = 0; d < NDIM; ++d) dX_sys.add_variable("dX_" + std::to_string(d));
    auto& W_sys = d_bdry_eq_sys->add_system<ExplicitSystem>(d_W_sys_name);
    W_sys.add_variable("W");
    W_sys.assemble_before_solve = false;
    W_sys.assemble();
}

BoundaryMeshMapping::~BoundaryMeshMapping()
{
    // Deallocate the wall sites
    Pointer<PatchHierarchy<NDIM>> hierarchy = d_bdry_data_manager->getPatchHierarchy();
    for (int ln = 0; ln <= hierarchy->getFinestLevelNumber(); ++ln)
    {
        Pointer<PatchLevel<NDIM>> level = hierarchy->getPatchLevel(ln);
        if (level->checkAllocated(d_w_idx)) level->deallocatePatchData(d_w_idx);
    }
}

void
BoundaryMeshMapping::setInitialConditions()
{
    plog << "Setting initial conditions\n";
    // Boundary moves according to what's in the velocity systems
    EquationSystems* eq_sys = d_bdry_data_manager->getEquationSystems();
    auto& W_sys = eq_sys->get_system(d_W_sys_name);
    const DofMap& W_dof_map = W_sys.get_dof_map();
    NumericVector<double>* W_vec = W_sys.solution.get();
    // Loop over nodes
    auto it = d_bdry_mesh->local_nodes_begin();
    const auto it_end = d_bdry_mesh->local_nodes_end();
    for (; it != it_end; ++it)
    {
        const Node* const node = *it;
        std::vector<dof_id_type> W_dofs;
        W_dof_map.dof_indices(node, W_dofs);
        W_vec->set(W_dofs[0], 1.0);
    }
    W_vec->close();
    W_sys.update();

    // Spread wall sites
    spreadWallSites(d_w_idx);
}

void
BoundaryMeshMapping::initializeEquationSystems()
{
    Pointer<PatchHierarchy<NDIM>> hierarchy = d_orig_data_manager->getPatchHierarchy();
    d_bdry_eq_sys->init();
    d_bdry_data_manager = FEDataManager::getManager(d_fe_data,
                                                    d_object_name + "::FEDataManager",
                                                    d_input_db,
                                                    hierarchy->getFinestLevelNumber() + 1,
                                                    d_orig_data_manager->getDefaultInterpSpec(),
                                                    d_orig_data_manager->getDefaultSpreadSpec(),
                                                    FEDataManager::WorkloadSpec());
    d_bdry_data_manager->COORDINATES_SYSTEM_NAME = d_coords_sys_name;
    d_bdry_data_manager->setPatchHierarchy(hierarchy);
    d_input_db.setNull();
    updateBoundaryLocation(0.0, 0.0, true);
    return;
}

void
BoundaryMeshMapping::updateBoundaryLocation(const double t_start, const double t_end, const bool initial_time)
{
    if (initial_time)
    {
        setInitialConditions();
    }
    // Move the boundary mesh to match the volume mesh.
    EquationSystems* eq_sys = d_orig_data_manager->getEquationSystems();

    System& X_system = eq_sys->get_system(d_orig_data_manager->COORDINATES_SYSTEM_NAME);
    const DofMap& X_dof_map = X_system.get_dof_map();
    NumericVector<double>* X_vec;
    X_vec = X_system.solution.get();

    System& X_bdry_sys = d_bdry_eq_sys->get_system(d_coords_sys_name);
    const DofMap& X_bdry_dof_map = X_bdry_sys.get_dof_map();
    NumericVector<double>* X_bdry_vec = X_bdry_sys.solution.get();

    System& dX_bdry_sys = d_bdry_eq_sys->get_system(d_disp_sys_name);
    NumericVector<double>* dX_bdry_vec = dX_bdry_sys.solution.get();

    auto node_it = d_bdry_mesh->local_nodes_begin();
    auto node_end = d_bdry_mesh->local_nodes_end();
    for (; node_it != node_end; ++node_it)
    {
        Node* node = *node_it;
        Node* orig_node = d_bdry_orig_node_map.at(node);
        dof_id_type bdry_node_id = node->id();
        dof_id_type orig_node_id = orig_node->id();
        // Grab current position of volumetric mesh.
        std::vector<dof_id_type> X_dof_indices, X_bdry_dof_indices;
        for (int d = 0; d < NDIM; ++d)
        {
            X_dof_map.dof_indices(d_orig_mesh->node_ptr(orig_node_id), X_dof_indices, d);
            X_bdry_dof_map.dof_indices(node, X_bdry_dof_indices, d);
            X_bdry_vec->set(X_bdry_dof_indices[0], (*X_vec)(X_dof_indices[0]));
            dX_bdry_vec->set(X_bdry_dof_indices[0], (*X_vec)(X_dof_indices[0]) - (*node)(d));
        }
    }
    X_bdry_vec->close();
    dX_bdry_vec->close();
    X_bdry_sys.update();
    dX_bdry_sys.update();

    spreadWallSites(d_w_idx);
    return;
}

void
BoundaryMeshMapping::spreadWallSites(int w_idx)
{
    Pointer<PatchHierarchy<NDIM>> hierarchy = d_orig_data_manager->getPatchHierarchy();
    if (w_idx == IBTK::invalid_index) w_idx = d_w_idx;
    // Double check that wall sites are allocated
    for (int ln = 0; ln <= hierarchy->getFinestLevelNumber(); ++ln)
    {
        Pointer<PatchLevel<NDIM>> level = hierarchy->getPatchLevel(ln);
        if (!level->checkAllocated(w_idx)) level->allocatePatchData(w_idx);
    }
    plog << "Spreading wall sites to the fluid.\n";
    // First zero out w_idx.
    // TODO: Should we always zero out w_idx? Maybe include an option
    for (int ln = 0; ln <= hierarchy->getFinestLevelNumber(); ++ln)
    {
        Pointer<PatchLevel<NDIM>> level = hierarchy->getPatchLevel(ln);
        for (PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            Pointer<Patch<NDIM>> patch = level->getPatch(p());
            Pointer<CellData<NDIM, double>> w_data = patch->getPatchData(w_idx);
            w_data->fillAll(0.0);
        }
    }

    d_bdry_data_manager->reinitElementMappings();
    NumericVector<double>* W_vec = d_bdry_data_manager->buildGhostedSolutionVector(d_W_sys_name);
    NumericVector<double>* X_vec =
        d_bdry_data_manager->buildGhostedSolutionVector(d_bdry_data_manager->COORDINATES_SYSTEM_NAME);
    d_bdry_data_manager->spread(w_idx, *W_vec, *X_vec, d_W_sys_name);
}
} // namespace IBAMR
