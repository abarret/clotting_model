#include "ibamr/app_namespaces.h"
#include <ibamr/IBFEMethod.h>

#include "ibtk/IBTK_MPI.h"
#include "ibtk/IndexUtilities.h"
#include "ibtk/libmesh_utilities.h"
#include <ibtk/LEInteractor.h>

#include "BoundaryMeshMapping.h"

#include "libmesh/enum_xdr_mode.h"
#include "libmesh/explicit_system.h"

namespace IBAMR
{
std::string
get_libmesh_restart_file_name(const std::string& restart_dump_dirname,
                              const std::string& base_filename,
                              unsigned int time_step_number,
                              unsigned int part,
                              const std::string& extension)
{
    std::ostringstream file_name_prefix;
    file_name_prefix << restart_dump_dirname << "/" << base_filename << "_part_" << part << "." << std::setw(6)
                     << std::setfill('0') << std::right << time_step_number << "." << extension;
    return file_name_prefix.str();
}
BoundaryMeshMapping::BoundaryMeshMapping(std::string object_name,
                                         Pointer<Database> input_db,
                                         MeshBase* vol_mesh,
                                         FEDataManager* vol_data_manager,
                                         const std::string& restart_read_dirname,
                                         const unsigned int restart_restore_number)
    : d_object_name(std::move(object_name)), d_vol_meshes({ vol_mesh }), d_vol_data_managers({ vol_data_manager })
{
    commonConstructor(input_db, restart_read_dirname, restart_restore_number);
}

BoundaryMeshMapping::BoundaryMeshMapping(std::string object_name,
                                         Pointer<Database> input_db,
                                         const std::vector<MeshBase*>& vol_meshes,
                                         const std::vector<FEDataManager*>& fe_data_managers,
                                         const std::string& restart_read_dirname,
                                         const unsigned int restart_restore_number)
    : d_object_name(std::move(object_name)), d_vol_meshes(vol_meshes), d_vol_data_managers(fe_data_managers)
{
    commonConstructor(input_db, restart_read_dirname, restart_restore_number);
}

BoundaryMeshMapping::BoundaryMeshMapping(std::string object_name) : d_object_name(std::move(object_name))
{
    // intentionally blank
}

void
BoundaryMeshMapping::commonConstructor(Pointer<Database> input_db,
                                       std::string restart_read_dirname,
                                       unsigned int restart_restore_number)
{
    d_abs_thresh = input_db->getDouble("abs_threshold");
    d_Sw = input_db->getDouble("spring_constant");
    const bool from_restart = RestartManager::getManager()->isFromRestart();
    unsigned int num_parts = d_vol_meshes.size();
    d_bdry_meshes.resize(num_parts);
    d_bdry_eq_sys_vec.resize(num_parts);
    d_fe_data.resize(num_parts);
    d_hierarchy = d_vol_data_managers[0]->getPatchHierarchy();
    d_input_db = input_db;
    for (unsigned int part = 0; part < num_parts; ++part)
    {
        // TODO: We need some way to extract a boundary mesh. For now, because we only have one part, we can extract the
        // part we want.
        d_bdry_meshes[part] =
            libmesh_make_unique<BoundaryMesh>(d_vol_meshes[part]->comm(), d_vol_meshes[part]->mesh_dimension() - 1);
        d_vol_meshes[part]->boundary_info->sync(*d_bdry_meshes[part]);
        // Now delete all the nodes of the bdry_mesh that are below the midline
        auto e_it_end = d_bdry_meshes[part]->elements_end();
        for (auto it = d_bdry_meshes[part]->elements_begin(); it != e_it_end; ++it)
        {
            Elem* e = *it;
            // If this element has nodes on the bottom half of the cylinder, delete it.
            for (unsigned int n = 0; n < e->n_nodes(); ++n)
            {
                if (e->node_ref(n)(1) < 0.0)
                {
                    d_bdry_meshes[part]->delete_elem(e);
                    break;
                }
            }
        }
        auto n_it_end = d_bdry_meshes[part]->nodes_end();
        for (auto it = d_bdry_meshes[part]->nodes_begin(); it != n_it_end; ++it)
        {
            Node* n = *it;
            // If the node is on the bottom half of the cylinder, delete it.
            if ((*n)(1) < 0.0) d_bdry_meshes[part]->delete_node(n);
        }
        d_bdry_meshes[part]->prepare_for_use();
        d_bdry_meshes[part]->print_info();
        d_bdry_eq_sys_vec[part] = std::move(libmesh_make_unique<EquationSystems>(*d_bdry_meshes[part]));
        d_fe_data[part] = std::make_shared<FEData>(
            d_object_name + "::FEData::" + std::to_string(part), *d_bdry_eq_sys_vec[part], true);

        if (from_restart)
        {
            const std::string& file_name = get_libmesh_restart_file_name(
                restart_read_dirname, d_object_name, restart_restore_number, part, d_libmesh_restart_file_extension);
            const XdrMODE xdr_mode = (d_libmesh_restart_file_extension == "xdr" ? DECODE : READ);
            const int read_mode =
                EquationSystems::READ_HEADER | EquationSystems::READ_DATA | EquationSystems::READ_ADDITIONAL_DATA;
            d_bdry_eq_sys_vec[part]->read(file_name, xdr_mode, read_mode, /*partition_agnostic*/ true);
        }
        else
        {
            auto& X_sys = d_bdry_eq_sys_vec[part]->add_system<ExplicitSystem>(d_coords_sys_name);
            for (unsigned int d = 0; d < NDIM; ++d) X_sys.add_variable("X_" + std::to_string(d), FEType());
            auto& dX_sys = d_bdry_eq_sys_vec[part]->add_system<ExplicitSystem>(d_disp_sys_name);
            for (unsigned int d = 0; d < NDIM; ++d) dX_sys.add_variable("dX_" + std::to_string(d), FEType());
            auto& U_sys = d_bdry_eq_sys_vec[part]->add_system<ExplicitSystem>(d_velocity_sys_name);
            for (unsigned int d = 0; d < NDIM; ++d) U_sys.add_variable("U_" + std::to_string(d), FEType());
            auto& Xi_sys = d_bdry_eq_sys_vec[part]->add_system<ExplicitSystem>(d_xi_sys_name);
            Xi_sys.add_variable(d_xi_sys_name);
            auto& Sigma_sys = d_bdry_eq_sys_vec[part]->add_system<ExplicitSystem>(d_sigma_sys_name);
            for (unsigned int d = 0; d < NDIM; ++d) Sigma_sys.add_variable("Sigma_" + std::to_string(d), FEType());
            auto& force_sys = d_bdry_eq_sys_vec[part]->add_system<ExplicitSystem>(d_force_sys_name);
            for (unsigned int d = 0; d < NDIM; ++d) force_sys.add_variable("Force_" + std::to_string(d), FEType());
            X_sys.assemble_before_solve = false;
            X_sys.assemble();
            dX_sys.assemble_before_solve = false;
            dX_sys.assemble();
        }
    }

    // Register forcing term
    auto var_db = VariableDatabase<NDIM>::getDatabase();
    d_f_var = new SideVariable<NDIM, double>(d_object_name + "::F_var");
    d_f_idx = var_db->registerVariableAndContext(
        d_f_var, var_db->getContext(d_object_name + "::context"), IntVector<NDIM>(6));
    return;
}

void
BoundaryMeshMapping::updateBoundaryLocation(const double t_start,
                                            const double t_end,
                                            const bool initial_time)
{
    if (initial_time)
    {
        setInitialConditions();
    }
    else
    {
        for (unsigned int part = 0; part < d_bdry_meshes.size(); ++part) updateBoundaryLocation(t_start, t_end, part);
    }
    return;
}

void
BoundaryMeshMapping::updateBoundaryLocation(const double t_start, const double t_end, const unsigned int part)
{
    const double dt = t_end - t_start;
    // Boundary moves according to what's in the velocity systems
    EquationSystems* eq_sys = d_bdry_eq_sys_vec[part].get();

    auto& X_sys = eq_sys->get_system(d_coords_sys_name);
    const DofMap& X_dof_map = X_sys.get_dof_map();
    NumericVector<double>* X_vec = X_sys.solution.get();
    auto& dX_sys = eq_sys->get_system(d_disp_sys_name);
    const DofMap& dX_dof_map = dX_sys.get_dof_map();
    NumericVector<double>* dX_vec = dX_sys.solution.get();
    auto& U_sys = eq_sys->get_system(d_velocity_sys_name);
    const DofMap& U_dof_map = U_sys.get_dof_map();
    NumericVector<double>* U_vec = U_sys.current_local_solution.get();

    // Loop over nodes
    auto it = d_bdry_meshes[part]->local_nodes_begin();
    const auto it_end = d_bdry_meshes[part]->local_nodes_end();
    for (; it != it_end; ++it)
    {
        const Node* const node = *it;
        std::vector<dof_id_type> X_dofs, dX_dofs, U_dofs;
        X_dof_map.dof_indices(node, X_dofs);
        dX_dof_map.dof_indices(node, dX_dofs);
        U_dof_map.dof_indices(node, U_dofs);
        // Use forward Euler
        for (int d = 0; d < NDIM; ++d)
        {
            X_vec->set(X_dofs[d], (*X_vec)(X_dofs[d]) + dt * (*U_vec)(U_dofs[d]));
            dX_vec->set(dX_dofs[d], (*node)(d) - (*X_vec)(X_dofs[d]));
        }
    }
    X_vec->close();
    dX_vec->close();
    X_sys.update();
    dX_sys.update();
}

void
BoundaryMeshMapping::setInitialConditions()
{
    for (unsigned int part = 0; part < d_bdry_meshes.size(); ++part)
    {
        // Boundary moves according to what's in the velocity systems
        EquationSystems* eq_sys = d_bdry_eq_sys_vec[part].get();

        auto& X_sys = eq_sys->get_system(d_coords_sys_name);
        const DofMap& X_dof_map = X_sys.get_dof_map();
        NumericVector<double>* X_vec = X_sys.solution.get();
        auto& dX_sys = eq_sys->get_system(d_disp_sys_name);
        const DofMap& dX_dof_map = dX_sys.get_dof_map();
        NumericVector<double>* dX_vec = dX_sys.solution.get();

        // Loop over nodes
        auto it = d_bdry_meshes[part]->local_nodes_begin();
        const auto it_end = d_bdry_meshes[part]->local_nodes_end();
        for (; it != it_end; ++it)
        {
            const Node* const node = *it;
            std::vector<dof_id_type> X_dofs, dX_dofs;
            X_dof_map.dof_indices(node, X_dofs);
            dX_dof_map.dof_indices(node, dX_dofs);
            for (int d = 0; d < NDIM; ++d)
            {
                X_vec->set(X_dofs[d], (*node)(d));
                dX_vec->set(dX_dofs[d], 0.0);
            }
        }
        X_vec->close();
        dX_vec->close();
        X_sys.update();
        dX_sys.update();
    }
}

void
BoundaryMeshMapping::initializeEquationSystems()
{
    const bool from_restart = RestartManager::getManager()->isFromRestart();
    if (from_restart) return;
    d_hierarchy = d_vol_data_managers[0]->getPatchHierarchy();
    for (unsigned int part = 0; part < d_bdry_meshes.size(); ++part)
    {
        d_bdry_eq_sys_vec[part]->init();
        d_bdry_data_managers.push_back(
            FEDataManager::getManager(d_fe_data[part],
                                      d_object_name + "::FEDataManager::" + std::to_string(part),
                                      d_input_db,
                                      d_hierarchy->getFinestLevelNumber() + 1,
                                      d_vol_data_managers[part]->getDefaultInterpSpec(),
                                      d_vol_data_managers[part]->getDefaultSpreadSpec(),
                                      FEDataManager::WorkloadSpec()));
        d_bdry_data_managers[part]->COORDINATES_SYSTEM_NAME = d_coords_sys_name;
        d_bdry_data_managers[part]->setPatchHierarchy(d_hierarchy);
    }
    updateBoundaryLocation(0.0, 0.0, true);
    return;
}

void
BoundaryMeshMapping::findVelocity(const double time, const int xi_idx, const int sigma_idx)
{
    for (unsigned int part = 0; part < d_bdry_meshes.size(); ++part)
    {
        d_bdry_data_managers[part]->reinitElementMappings();
        const std::unique_ptr<BoundaryMesh>& bdry_mesh = d_bdry_meshes[part];
        EquationSystems* eq_sys = d_bdry_eq_sys_vec[part].get();
        // First we need to interpolate xi and sigma to the structure. Use the FEDataManager for this
        auto& Xw_sys = eq_sys->get_system<ExplicitSystem>(d_coords_sys_name);
        const DofMap& Xw_dof_map = Xw_sys.get_dof_map();
        NumericVector<double>* Xw_vec = Xw_sys.solution.get();
        auto& xi_sys = eq_sys->get_system<ExplicitSystem>(d_xi_sys_name);
        const DofMap& xi_dof_map = xi_sys.get_dof_map();
        NumericVector<double>* xi_vec = xi_sys.solution.get();
        d_bdry_data_managers[part]->interp(xi_idx, *xi_vec, *Xw_vec, d_xi_sys_name);
        auto& sigma_sys = eq_sys->get_system<ExplicitSystem>(d_sigma_sys_name);
        NumericVector<double>* sigma_vec = sigma_sys.solution.get();
        d_bdry_data_managers[part]->interp(sigma_idx, *sigma_vec, *Xw_vec, d_sigma_sys_name);

        // Then solve for U_w. We need U, X_w, X, Sigma, and Xi to solve for U_w.
        EquationSystems* vol_eq_sys = d_vol_data_managers[part]->getEquationSystems();
        auto& X_sys = vol_eq_sys->get_system<ExplicitSystem>(d_vol_data_managers[part]->COORDINATES_SYSTEM_NAME);
        const DofMap& X_dof_map = X_sys.get_dof_map();
        NumericVector<double>* X_vec = X_sys.current_local_solution.get();
        auto& U_sys = vol_eq_sys->get_system<ExplicitSystem>(IBFEMethod::VELOCITY_SYSTEM_NAME);
        NumericVector<double>* U_vec = U_sys.current_local_solution.get();

        auto& Uw_sys = eq_sys->get_system<ExplicitSystem>(d_velocity_sys_name);
        NumericVector<double>* Uw_vec = Uw_sys.solution.get();

        std::map<dof_id_type, dof_id_type> node_id_map;
        std::map<dof_id_type, unsigned char> side_id_map;
        d_vol_meshes[part]->boundary_info->get_side_and_node_maps(*bdry_mesh, node_id_map, side_id_map);
        auto it = bdry_mesh->local_nodes_begin();
        const auto it_end = bdry_mesh->local_nodes_end();
        for (; it != it_end; ++it)
        {
            const Node* const node = *it;
            dof_id_type bdry_node_id = node->id();
            // Determine volume node. TODO: This is potentially expensive. We should cache our own map.
            auto vol_iter = std::find_if(
                node_id_map.begin(), node_id_map.end(), [bdry_node_id](const std::pair<dof_id_type, dof_id_type>& obj) {
                    return obj.second == bdry_node_id;
                });
            dof_id_type vol_node_id = vol_iter->first;
            // Grab the current location and velocity of the volumetric mesh and location of boundary mesh
            IBTK::VectorNd X, U, Xw;
            std::vector<dof_id_type> dof_indices;
            for (int d = 0; d < NDIM; ++d)
            {
                X_dof_map.dof_indices(d_vol_meshes[part]->node_ptr(vol_node_id), dof_indices, d);
                X[d] = (*X_vec)(dof_indices[0]);
                U[d] = (*U_vec)(dof_indices[0]);
                Xw_dof_map.dof_indices(node, dof_indices, d);
                Xw[d] = (*Xw_vec)(dof_indices[0]);
            }
            // Now grab values of Xi and Sigma
            double xi = 0.0, sigma = 0.0;
            xi_dof_map.dof_indices(node, dof_indices);
            xi = (*xi_vec)(dof_indices[0]);
            sigma = (*sigma_vec)(dof_indices[0]);

            // Now we can find Uw
            if (xi < d_abs_thresh)
            {
                // Basically zero, we set Uw equal to U
                for (int d = 0; d < NDIM; ++d)
                {
                    Xw_dof_map.dof_indices(node, dof_indices, d);
                    Uw_vec->set(dof_indices[d], U[d]);
                }
            }
            else
            {
                // Velocities are different.
                for (int d = 0; d < NDIM; ++d)
                {
                    Xw_dof_map.dof_indices(node, dof_indices, d);
                    double uw = U[d] + (sigma + d_Sw * (X[d] - Xw[d])) / xi;
                    Uw_vec->set(dof_indices[0], uw);
                }
            }
        }
    }
}

void
BoundaryMeshMapping::findForce(const double time, const int f_idx)
{
    for (unsigned int part = 0; part < d_bdry_meshes.size(); ++part)
    {
        const std::unique_ptr<BoundaryMesh>& bdry_mesh = d_bdry_meshes[part];
        EquationSystems* eq_sys = d_bdry_eq_sys_vec[part].get();
        // We need to pull out the boundary systems
        auto& Xw_sys = eq_sys->get_system<ExplicitSystem>(d_coords_sys_name);
        const DofMap& Xw_dof_map = Xw_sys.get_dof_map();
        NumericVector<double>* Xw_vec = Xw_sys.solution.get();
        auto& xi_sys = eq_sys->get_system<ExplicitSystem>(d_xi_sys_name);
        const DofMap& xi_dof_map = xi_sys.get_dof_map();
        NumericVector<double>* xi_vec = xi_sys.solution.get();
        auto& Uw_sys = eq_sys->get_system<ExplicitSystem>(d_velocity_sys_name);
        NumericVector<double>* Uw_vec = Uw_sys.solution.get();

        auto& F_sys = eq_sys->get_system<ExplicitSystem>(d_force_sys_name);
        const DofMap& F_dof_map = F_sys.get_dof_map();
        NumericVector<double>* F_vec = F_sys.solution.get();

        // Pull out volumetric systems
        EquationSystems* vol_eq_sys = d_vol_data_managers[part]->getEquationSystems();
        auto& X_sys = vol_eq_sys->get_system<ExplicitSystem>(d_vol_data_managers[part]->COORDINATES_SYSTEM_NAME);
        const DofMap& X_dof_map = X_sys.get_dof_map();
        NumericVector<double>* X_vec = X_sys.current_local_solution.get();
        auto& U_sys = vol_eq_sys->get_system<ExplicitSystem>(IBFEMethod::VELOCITY_SYSTEM_NAME);
        NumericVector<double>* U_vec = U_sys.current_local_solution.get();

        // Build map between volume and boundary mesh.
        std::map<dof_id_type, dof_id_type> node_id_map;
        std::map<dof_id_type, unsigned char> side_id_map;
        d_vol_meshes[part]->boundary_info->get_side_and_node_maps(*bdry_mesh, node_id_map, side_id_map);
        auto it = bdry_mesh->local_nodes_begin();
        const auto it_end = bdry_mesh->local_nodes_end();
        for (; it != it_end; ++it)
        {
            const Node* const node = *it;
            dof_id_type bdry_node_id = node->id();
            // Determine volume node. TODO: This is potentially expensive. We should cache our own map.
            auto vol_iter = std::find_if(
                node_id_map.begin(), node_id_map.end(), [bdry_node_id](const std::pair<dof_id_type, dof_id_type>& obj) {
                    return obj.second == bdry_node_id;
                });
            dof_id_type vol_node_id = vol_iter->first;
            // Grab the current location and velocity of the volumetric mesh and location of boundary mesh
            IBTK::VectorNd X, U, Xw, Uw;
            std::vector<dof_id_type> dof_indices;
            // Grab position and velocity values.
            for (int d = 0; d < NDIM; ++d)
            {
                X_dof_map.dof_indices(d_vol_meshes[part]->node_ptr(vol_node_id), dof_indices, d);
                X[d] = (*X_vec)(dof_indices[0]);
                U[d] = (*U_vec)(dof_indices[0]);
                Xw_dof_map.dof_indices(node, dof_indices, d);
                Xw[d] = (*Xw_vec)(dof_indices[0]);
                Uw[d] = (*Uw_vec)(dof_indices[0]);
            }
            // Grab value of xi.
            double xi = 0.0;
            xi_dof_map.dof_indices(node, dof_indices);
            xi = (*xi_vec)(dof_indices[0]);

            // Now we can find the force
            if (xi < d_abs_thresh)
            {
                // Basically zero, we set F equal to 0
                for (int d = 0; d < NDIM; ++d)
                {
                    Xw_dof_map.dof_indices(node, dof_indices, d);
                    F_vec->set(dof_indices[d], 0.0);
                }
            }
            else
            {
                // A force is being applied
                for (int d = 0; d < NDIM; ++d)
                {
                    Xw_dof_map.dof_indices(node, dof_indices, d);
                    double Fw = d_Sw * (Xw[d] - X[d]) + xi * (Uw[d] - U[d]);
                    F_vec->set(dof_indices[0], Fw);
                }
            }
        }

        // Now update the forces
        F_vec->close();
        F_sys.update();
        // Now spread the force into d_f_idx:
        // We need to spread from volumetric points.
        // Figure out which points are inside patch
        const std::vector<std::vector<Node*>>& active_bdry_nodes = d_bdry_data_managers[part]->getActivePatchNodeMap();
        // Assume we are on the finest level
        {
            Pointer<PatchLevel<NDIM>> level = d_hierarchy->getPatchLevel(d_hierarchy->getFinestLevelNumber());
            int local_patch_num = 0;
            for (PatchLevel<NDIM>::Iterator p(level); p; p++, ++local_patch_num)
            {
                Pointer<Patch<NDIM>> patch = level->getPatch(p());
                const std::vector<Node*>& patch_nodes = active_bdry_nodes[local_patch_num];
                const size_t num_active_patch_nodes = patch_nodes.size();
                if (patch_nodes.empty()) continue;

                Pointer<CartesianPatchGeometry<NDIM>> pgeom = patch->getPatchGeometry();
                const double* const patch_lower = pgeom->getXLower();
                const double* const patch_upper = pgeom->getXUpper();
                std::array<bool, NDIM> touches_upper_regular_bdry;
                for (int d = 0; d < NDIM; ++d) touches_upper_regular_bdry[d] = pgeom->getTouchesRegularBoundary(d, 1);

                std::vector<double> F_node, X_node;
                X_node.reserve(NDIM * num_active_patch_nodes);
                F_node.reserve(NDIM * num_active_patch_nodes);
                std::vector<dof_id_type> F_idxs, X_idxs;
                IBTK::Point X;
                for (const auto& n : patch_nodes)
                {
                    bool inside_patch = true;
                    // Get corresponding node on volumetric mesh
                    dof_id_type bdry_node_id = n->id();
                    // Determine volume node. TODO: This is potentially expensive. We should cache our own map.
                    auto vol_iter = std::find_if(node_id_map.begin(),
                                                 node_id_map.end(),
                                                 [bdry_node_id](const std::pair<dof_id_type, dof_id_type>& obj) {
                                                     return obj.second == bdry_node_id;
                                                 });
                    dof_id_type vol_node_id = vol_iter->first;
                    const Node* const vol_n = d_vol_meshes[part]->node_ptr(vol_node_id);
                    for (unsigned int d = 0; d < NDIM; ++d)
                    {
                        IBTK::get_nodal_dof_indices(X_dof_map, vol_n, d, X_idxs);
                        X[d] = (*X_vec)(X_idxs[0]);
                        inside_patch =
                            inside_patch && (X[d] >= patch_lower[d]) &&
                            ((X[d] < patch_upper[d]) || (touches_upper_regular_bdry[d] && X[d] <= patch_upper[d]));
                    }
                    if (inside_patch)
                    {
                        for (unsigned int i = 0; i < NDIM; ++i)
                        {
                            IBTK::get_nodal_dof_indices(F_dof_map, n, i, F_idxs);
                            for (const auto& F_idx : F_idxs) F_node.push_back((*F_vec)(F_idx));
                        }
                        X_node.insert(X_node.end(), &X[0], &X[0] + NDIM);
                    }
                }
                TBOX_ASSERT(F_node.size() <= NDIM * num_active_patch_nodes);
                TBOX_ASSERT(X_node.size() <= NDIM * num_active_patch_nodes);

                // Now we spread values
                const Box<NDIM>& spread_box = patch->getBox();
                Pointer<SideData<NDIM, double>> f_data = patch->getPatchData(d_f_idx);
                LEInteractor::spread(f_data,
                                     F_node,
                                     NDIM,
                                     X_node,
                                     NDIM,
                                     patch,
                                     spread_box,
                                     d_bdry_data_managers[part]->getDefaultInterpSpec().kernel_fcn);
            }
        }
        // Now copy data back to f_idx.
        for (int ln = 0; ln <= d_hierarchy->getFinestLevelNumber(); ++ln)
        {
            Pointer<PatchLevel<NDIM>> level = d_hierarchy->getPatchLevel(ln);
            for (PatchLevel<NDIM>::Iterator p(level); p; p++)
            {
                Pointer<Patch<NDIM>> patch = level->getPatch(p());
                Pointer<SideData<NDIM, double>> src_data = patch->getPatchData(d_f_idx);
                Pointer<PatchData<NDIM>> dst_data = patch->getPatchData(f_idx);
                dst_data->copy(*src_data);
            }
        }
    }
}

void
BoundaryMeshMapping::setDataOnPatchHierarchy(int data_idx,
                                             Pointer<hier::Variable<NDIM>> var,
                                             Pointer<PatchHierarchy<NDIM>> hierarchy,
                                             const double data_time,
                                             const bool initial_time,
                                             int coarsest_ln,
                                             int finest_ln)
{
    // We need to fill in the force data
    if (initial_time)
    {
        // Just set it to zero at initial time.
        coarsest_ln = coarsest_ln < 0 ? 0 : coarsest_ln;
        finest_ln = finest_ln < 0 ? hierarchy->getFinestLevelNumber() : finest_ln;
        for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
        {
            Pointer<PatchLevel<NDIM>> level = hierarchy->getPatchLevel(ln);
            setDataOnPatchLevel(data_idx, var, level, data_time, initial_time);
        }
        return;
    }
    // Allocate patch data
    for (int ln = 0; ln <= hierarchy->getFinestLevelNumber(); ++ln)
    {
        Pointer<PatchLevel<NDIM>> level = hierarchy->getPatchLevel(ln);
        if (!level->checkAllocated(d_f_idx)) level->allocatePatchData(d_f_idx);
        for (PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            Pointer<Patch<NDIM>> patch = level->getPatch(p());
            Pointer<SideData<NDIM, double>> f_data = patch->getPatchData(d_f_idx);
            f_data->fillAll(0.0);
        }
    }

    // Find force
    findForce(data_time, data_idx);

    // Deallocate patch data
    for (int ln = 0; ln <= hierarchy->getFinestLevelNumber(); ++ln)
    {
        Pointer<PatchLevel<NDIM>> level = hierarchy->getPatchLevel(ln);
        if (level->checkAllocated(d_f_idx)) level->deallocatePatchData(d_f_idx);
    }
}

void
BoundaryMeshMapping::setDataOnPatch(const int data_idx,
                                    Pointer<hier::Variable<NDIM>> /*var*/,
                                    Pointer<Patch<NDIM>> patch,
                                    const double /*data_time*/,
                                    const bool /*initial_time*/,
                                    Pointer<PatchLevel<NDIM>> /*level*/)
{
    Pointer<CellData<NDIM, double>> cc_data = patch->getPatchData(data_idx);
    Pointer<SideData<NDIM, double>> sc_data = patch->getPatchData(data_idx);
    Pointer<FaceData<NDIM, double>> fc_data = patch->getPatchData(data_idx);
    if (cc_data)
        cc_data->fillAll(0.0);
    else if (sc_data)
        sc_data->fillAll(0.0);
    else if (fc_data)
        fc_data->fillAll(0.0);
    else
        TBOX_ERROR(d_object_name + ": Unknown data centering\n");
}

void
BoundaryMeshMapping::writeFEDataToRestartFile(const std::string& restart_dump_dirname, unsigned int time_step_number)
{
    for (unsigned int part = 0; part < d_bdry_eq_sys_vec.size(); ++part)
    {
        const std::string& file_name = get_libmesh_restart_file_name(
            restart_dump_dirname, d_object_name, time_step_number, part, d_libmesh_restart_file_extension);
        const XdrMODE xdr_mode = (d_libmesh_restart_file_extension == "xdr" ? ENCODE : WRITE);
        const int write_mode = EquationSystems::WRITE_DATA | EquationSystems::WRITE_ADDITIONAL_DATA;
        d_bdry_eq_sys_vec[part]->write(file_name, xdr_mode, write_mode, /*partition_agnostic*/ true);
    }
}
} // namespace IBAMR
