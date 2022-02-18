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

void
BoundaryMeshMapping::commonConstructor(Pointer<Database> input_db,
                                       std::string restart_read_dirname,
                                       unsigned int restart_restore_number)
{
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
            auto& W_sys = d_bdry_eq_sys_vec[part]->add_system<ExplicitSystem>(d_W_sys_name);
            W_sys.add_variable("W", FEType());
            W_sys.assemble_before_solve = false;
            W_sys.assemble();
        }
    }

    return;
}

void
BoundaryMeshMapping::updateBoundaryLocation(const double t_start, const double t_end, const bool initial_time)
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
    FEDataManager* fe_data_manager = d_vol_data_managers[part];
    EquationSystems* eq_sys = fe_data_manager->getEquationSystems();

    System& X_system = eq_sys->get_system(fe_data_manager->COORDINATES_SYSTEM_NAME);
    const DofMap& X_dof_map = X_system.get_dof_map();
    NumericVector<double>* X_vec;
    X_vec = X_system.solution.get();

    System& X_bdry_sys = d_bdry_eq_sys_vec[part]->get_system(d_coords_sys_name);
    const DofMap& X_bdry_dof_map = X_bdry_sys.get_dof_map();
    NumericVector<double>* X_bdry_vec = X_bdry_sys.solution.get();

    System& dX_bdry_sys = d_bdry_eq_sys_vec[part]->get_system(d_disp_sys_name);
    NumericVector<double>* dX_bdry_vec = dX_bdry_sys.solution.get();

    std::map<dof_id_type, dof_id_type> node_id_map;
    std::map<dof_id_type, unsigned char> side_id_map;
    d_vol_meshes[part]->boundary_info->get_side_and_node_maps(*d_bdry_meshes[part], node_id_map, side_id_map);
    auto node_it = d_bdry_meshes[part]->local_nodes_begin();
    auto node_end = d_bdry_meshes[part]->local_nodes_end();
    for (; node_it != node_end; ++node_it)
    {
        Node* node = *node_it;
        dof_id_type bdry_node_id = node->id();
        // TODO: This is potentially expensive. We should cache our own map between bdry nodes and volumetric nodes.
        auto vol_iter = std::find_if(
            node_id_map.begin(), node_id_map.end(), [bdry_node_id](const std::pair<dof_id_type, dof_id_type>& obj) {
                return obj.second == bdry_node_id;
            });
        dof_id_type vol_node_id = vol_iter->first;
        // Grab current position of volumetric mesh.
        std::vector<dof_id_type> X_dof_indices, X_bdry_dof_indices;
        for (int d = 0; d < NDIM; ++d)
        {
            X_dof_map.dof_indices(d_vol_meshes[part]->node_ptr(vol_node_id), X_dof_indices, d);
            X_bdry_dof_map.dof_indices(node, X_bdry_dof_indices, d);
            X_bdry_vec->set(X_bdry_dof_indices[0], (*X_vec)(X_dof_indices[0]));
            dX_bdry_vec->set(X_bdry_dof_indices[0], (*X_vec)(X_dof_indices[0]) - (*node)(d));
        }
    }
    X_bdry_vec->close();
    dX_bdry_vec->close();
    X_bdry_sys.update();
    dX_bdry_sys.update();
    return;
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
    d_input_db.setNull();
    updateBoundaryLocation(0.0, 0.0, true);
    return;
}

void
BoundaryMeshMapping::spreadWallSites(const int w_idx)
{
    for (const auto& bdry_data_manager : d_bdry_data_managers)
    {
        EquationSystems* eq_sys = bdry_data_manager->getEquationSystems();
        std::unique_ptr<NumericVector<double>> W_vec = bdry_data_manager->buildIBGhostedVector(d_W_sys_name);
        std::unique_ptr<NumericVector<double>> X_vec =
            bdry_data_manager->buildIBGhostedVector(bdry_data_manager->COORDINATES_SYSTEM_NAME);
        bdry_data_manager->spread(w_idx, *W_vec, *X_vec, d_W_sys_name);
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
