#include "ibamr/app_namespaces.h"
#include <ibamr/IBFEMethod.h>

#include "ibtk/IBTK_MPI.h"
#include "ibtk/IndexUtilities.h"
#include "ibtk/libmesh_utilities.h"

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
            auto& Sigma_sys = d_bdry_eq_sys_vec[part]->add_system<ExplicitSystem>(d_sigma_sys_name);
            for (unsigned int d = 0; d < NDIM; ++d) Sigma_sys.add_variable("Sigma_" + std::to_string(d), FEType());
            auto& force_sys = d_bdry_eq_sys_vec[part]->add_system<ExplicitSystem>(d_force_sys_name);
            for (unsigned int d = 0; d < NDIM; ++d) force_sys.add_variable("Force_" + std::to_string(d), FEType());
            X_sys.assemble_before_solve = false;
            X_sys.assemble();
            dX_sys.assemble_before_solve = false;
            dX_sys.assemble();
        }

        d_bdry_data_managers.push_back(
            FEDataManager::getManager(d_fe_data[part],
                                      d_object_name + "::FEDataManager::" + std::to_string(part),
                                      input_db,
                                      d_hierarchy->getFinestLevelNumber(),
                                      d_vol_data_managers[part]->getDefaultInterpSpec(),
                                      d_vol_data_managers[part]->getDefaultSpreadSpec(),
                                      FEDataManager::WorkloadSpec()));
    }
    return;
}

void
BoundaryMeshMapping::updateBoundaryLocation(const double t_start,
                                            const double t_end,
                                            const bool end_of_timestep,
                                            const bool initial_time)
{
    if (initial_time)
    {
        setInitialConditions();
    }
    else
    {
        for (unsigned int part = 0; part < d_bdry_meshes.size(); ++part)
            updateBoundaryLocation(t_start, t_end, part, end_of_timestep);
    }
    return;
}

void
BoundaryMeshMapping::updateBoundaryLocation(const double t_start,
                                            const double t_end,
                                            const unsigned int part,
                                            const bool end_of_timestep)
{
    const double dt = t_end - t_start;
    // Boundary moves according to what's in the velocity systems
    EquationSystems* eq_sys = d_bdry_eq_sys_vec[part];

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
        if (initial_time)
        {
            X_vec->set(X_dofs[d], (*node)(d));
        }
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
        EquationSystems* eq_sys = d_bdry_eq_sys_vec[part];

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
    for (unsigned int part = 0; part < d_bdry_meshes.size(); ++part)
    {
        d_bdry_eq_sys_vec[part]->init();
    }
    updateBoundaryLocation(0.0, 0.0, false);
    return;
}

void
BoundaryMeshMapping::findVelocity(const double time, const int xi_idx, const int sigma_idx)
{
    for (unsigned int part = 0; part < d_bdry_meshes.size(); ++part)
    {
        const std::unique_ptr<BoundaryMesh>& bdry_mesh = d_bdry_meshes[part];
        EquationSystems* eq_sys = d_bdry_eq_sys_vec[part];
        // First we need to interpolate xi and sigma to the structure. Use the FEDataManager for this
        ExplicitSystem& Xw_sys = eq_sys->get_system(d_coords_sys_name);
        const DofMap& Xw_dof_map = Xw_sys.get_dof_map();
        NumericVector<double>* Xw_vec = Xw_sys.solution.get();
        ExplicitSystem& xi_sys = eq_sys->get_system(d_xi_sys_name);
        const DofMap& xi_dof_map = xi_sys.get_dof_map();
        NumericVector<double>* xi_vec = xi_sys.solution.get();
        d_bdry_data_managers[part]->interp(xi_idx, *xi_vec, *Xw_vec, d_xi_sys_name);
        ExplicitSystem& sigma_sys = eq_sys->get_system(d_sigma_sys_name);
        NumericVector<double>* sigma_vec = sigma_sys.solution.get();
        d_bdry_data_managers[part]->interp(sigma_idx, *sigma_vec, *Xw_vec, d_sigma_sys_name);

        // Then solve for U_w. We need U, X_w, X, Sigma, and Xi to solve for U_w.
        EquationSystems* vol_eq_sys = d_vol_data_managers[part];
        ExplicitSystem& X_sys = vol_eq_sys->get_system(d_vol_data_managers[part]->COORDINATES_SYSTEM_NAME);
        const DofMap& X_dof_map = X_sys.get_dof_map();
        NumericVector<double>* X_vec = X_sys.current_local_solution.get();
        ExplicitSystem& U_sys = vol_eq_sys->get_system(IBFEMethod::VELOCITY_SYSTEM_NAME);
        NumericVector<double>* U_vec = U_sys.current_local_solution.get();

        ExplicitSystem& Uw_sys = eq_sys->get_system(d_velocity_sys_name);
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
        EquationSystems* eq_sys = d_bdry_eq_sys_vec[part];
        // We need to pull out the boundary systems
        ExplicitSystem& Xw_sys = eq_sys->get_system(d_coords_sys_name);
        const DofMap& Xw_dof_map = Xw_sys.get_dof_map();
        NumericVector<double>* Xw_vec = Xw_sys.solution.get();
        ExplicitSystem& xi_sys = eq_sys->get_system(d_xi_sys_name);
        const DofMap& xi_dof_map = xi_sys.get_dof_map();
        NumericVector<double>* xi_vec = xi_sys.solution.get();
        ExplicitSystem& Uw_sys = eq_sys->get_system(d_velocity_sys_name);
        NumericVector<double>* Uw_vec = Uw_sys.solution.get();

        ExplicitSystem& F_sys = eq_sys->get_system(d_force_sys_name);
        NumericVector<double>* F_vec = F_sys.solution.get();

        // Pull out volumetric systems
        EquationSystems* vol_eq_sys = d_vol_data_managers[part];
        ExplicitSystem& X_sys = vol_eq_sys->get_system(d_vol_data_managers[part]->COORDINATES_SYSTEM_NAME);
        const DofMap& X_dof_map = X_sys.get_dof_map();
        NumericVector<double>* X_vec = X_sys.current_local_solution.get();
        ExplicitSystem& U_sys = vol_eq_sys->get_system(IBFEMethod::VELOCITY_SYSTEM_NAME);
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
        // Now spread the force into f_idx:
        d_bdry_data_managers[part]->spread(f_idx, *F_vec, *Xw_vec, d_force_sys_name);
    }
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
