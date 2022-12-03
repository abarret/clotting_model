// ---------------------------------------------------------------------
//
// Copyright (c) 2014 - 2020 by the IBAMR developers
// All rights reserved.
//
// This file is part of IBAMR.
//
// IBAMR is free software and is distributed under the 3-clause BSD
// license. The full text of the license can be found in the file
// COPYRIGHT at the top level directory of IBAMR.
//
// ---------------------------------------------------------------------

#include <ibamr/config.h>

#include <clot/BondBoundSource.h>
#include <clot/app_namespaces.h>

#include <HierarchyDataOpsManager.h>
#include <SAMRAI_config.h>
#include <math.h>

namespace clot
{
/////////////////////////////// PUBLIC ///////////////////////////////////////
BondBoundSource::BondBoundSource(std::string object_name, const BoundClotParams& clot_params)
    : CartGridFunction(std::move(object_name)), d_clot_params(clot_params)
{
    // intentionally blank
} // BondBoundSource

bool
BondBoundSource::isTimeDependent() const
{
    return true;
} // isTimeDependent

void
BondBoundSource::setDataOnPatchHierarchy(const int data_idx,
                                         Pointer<Variable<NDIM>> var,
                                         Pointer<PatchHierarchy<NDIM>> hierarchy,
                                         const double data_time,
                                         const bool initial_time,
                                         int coarsest_ln,
                                         int finest_ln)
{
    coarsest_ln = (coarsest_ln == -1 ? 0 : coarsest_ln);
    finest_ln = (finest_ln == -1 ? hierarchy->getFinestLevelNumber() : finest_ln);

    // Allocate scratch data
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        Pointer<PatchLevel<NDIM>> level = hierarchy->getPatchLevel(ln);
        level->allocatePatchData(d_phi_b_scr_idx, data_time);
        level->allocatePatchData(d_bond_scr_idx, data_time);
    }

    std::map<Pointer<Variable<NDIM>>, bool> scratch_allocated;
    for (const auto& var_integrator_pair : d_var_integrator_pairs)
    {
        Pointer<Variable<NDIM>> var = var_integrator_pair.first;
        Pointer<HierarchyIntegrator> integrator = var_integrator_pair.second;
        // Check that things have been set correctly
        TBOX_ASSERT(var && integrator);
        // Allocate scratch data when needed.
        auto var_db = VariableDatabase<NDIM>::getDatabase();
        int var_scr_idx = var_db->mapVariableAndContextToIndex(var, integrator->getScratchContext());
        scratch_allocated[var] = integrator->isAllocatedPatchData(var_scr_idx);
        if (!scratch_allocated.at(var)) integrator->allocatePatchData(var_scr_idx, data_time);
        if (!initial_time)
        {
            int var_cur_idx = var_db->mapVariableAndContextToIndex(var, integrator->getCurrentContext());
            int var_new_idx = var_db->mapVariableAndContextToIndex(var, integrator->getNewContext());
            const bool var_new_is_allocated = integrator->isAllocatedPatchData(var_scr_idx);
            auto hier_data_ops_manager = HierarchyDataOpsManager<NDIM>::getManager();
            Pointer<HierarchyDataOpsReal<NDIM, double>> hier_cc_data_ops =
                hier_data_ops_manager->getOperationsDouble(var, hierarchy, /*get_unique*/ true);
            if (integrator->getCurrentCycleNumber() == 0 || !var_new_is_allocated)
                hier_cc_data_ops->copyData(var_scr_idx, var_cur_idx);
            else
                hier_cc_data_ops->linearSum(var_scr_idx, 0.5, var_cur_idx, 0.5, var_new_idx);
            // Do we need to compute ghost cells?
        }
    }

    // Now fill in ghost cells for things that need it
    using ITC = HierarchyGhostCellInterpolation::InterpolationTransactionComponent;
    std::vector<ITC> ghost_cell_comps(2);
    auto var_db = VariableDatabase<NDIM>::getDatabase();

    Pointer<Variable<NDIM>> phi_b_var = d_var_integrator_pairs[0].first;
    Pointer<AdvDiffHierarchyIntegrator> phi_b_integrator = d_var_integrator_pairs[0].second;
    const int phi_b_scr_idx = var_db->mapVariableAndContextToIndex(phi_b_var, phi_b_integrator->getScratchContext());
    ghost_cell_comps[0] = ITC(d_phi_b_scr_idx,
                              phi_b_scr_idx,
                              "CONSERVATIVE_LINEAR_REFINE",
                              false,
                              "NONE",
                              "LINEAR",
                              false,
                              phi_b_integrator->getPhysicalBcCoefs(phi_b_var));

    Pointer<Variable<NDIM>> bond_var = d_var_integrator_pairs[2].first;
    Pointer<AdvDiffHierarchyIntegrator> bond_integrator = d_var_integrator_pairs[2].second;
    const int bond_scr_idx = var_db->mapVariableAndContextToIndex(bond_var, bond_integrator->getScratchContext());
    ghost_cell_comps[1] = ITC(d_bond_scr_idx,
                              bond_scr_idx,
                              "CONSERVATIVE_LINEAR_REFINE",
                              false,
                              "NONE",
                              "LINEAR",
                              false,
                              bond_integrator->getPhysicalBcCoefs(bond_var));

    HierarchyGhostCellInterpolation hier_ghost_cell;
    hier_ghost_cell.initializeOperatorState(ghost_cell_comps, hierarchy, coarsest_ln, finest_ln);
    hier_ghost_cell.fillData(data_time);

    // Compute the source function
    for (int level_num = coarsest_ln; level_num <= finest_ln; ++level_num)
    {
        setDataOnPatchLevel(data_idx, var, hierarchy->getPatchLevel(level_num), data_time, initial_time);
    }

    // Deallocate scratch data when needed.
    for (const auto& var_integrator_pair : d_var_integrator_pairs)
    {
        Pointer<Variable<NDIM>> var = var_integrator_pair.first;
        Pointer<HierarchyIntegrator> integrator = var_integrator_pair.second;
        if (!scratch_allocated.at(var))
        {
            auto var_db = VariableDatabase<NDIM>::getDatabase();
            int var_scr_idx = var_db->mapVariableAndContextToIndex(var, integrator->getScratchContext());
            integrator->deallocatePatchData(var_scr_idx);
        }
    }

    // Deallocate scratch data
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        Pointer<PatchLevel<NDIM>> level = hierarchy->getPatchLevel(ln);
        level->deallocatePatchData(d_phi_b_scr_idx);
        level->deallocatePatchData(d_bond_scr_idx);
    }
    return;
} // setDataOnPatchHierarchy

void
BondBoundSource::setDataOnPatch(const int data_idx,
                                Pointer<Variable<NDIM>> /*var*/,
                                Pointer<Patch<NDIM>> patch,
                                const double /*data_time*/,
                                const bool initial_time,
                                Pointer<PatchLevel<NDIM>> /*patch_level*/)
{
    Pointer<CellData<NDIM, double>> ret_data = patch->getPatchData(data_idx);
    ret_data->fillAll(0.0);
    if (initial_time) return;

    // Patch information
    const Box<NDIM>& patch_box = patch->getBox();
    Pointer<CartesianPatchGeometry<NDIM>> pgeom = patch->getPatchGeometry();
    const double* const dx = pgeom->getDx();

    Pointer<CellData<NDIM, double>> sig_data =
        patch->getPatchData(d_var_integrator_pairs[3].first, d_var_integrator_pairs[3].second->getScratchContext());
    Pointer<CellData<NDIM, double>> phi_a_data =
        patch->getPatchData(d_var_integrator_pairs[1].first, d_var_integrator_pairs[1].second->getScratchContext());
    Pointer<CellData<NDIM, double>> phi_b_data = patch->getPatchData(d_phi_b_scr_idx);
    Pointer<CellData<NDIM, double>> bond_data = patch->getPatchData(d_bond_scr_idx);
    Pointer<CellData<NDIM, double>> w_data = patch->getPatchData(d_w_idx);
    auto psi_fcn = getKernelAndWidth(d_kernel);
    for (CellIterator<NDIM> i(patch_box); i; i++)
    {
        const CellIndex<NDIM>& idx = i();
        const double phi_a = (*phi_a_data)(idx);
        const double phi_b = (*phi_b_data)(idx);
        const double z = (*bond_data)(idx);
        const double w = (*w_data)(idx);
        // Compute sources
        double R2 =
            d_clot_params.Kab * phi_a *
            std::max(
                convolution(
                    1.0, phi_b_data.getPointer(), -2.0, bond_data.getPointer(), psi_fcn.first, psi_fcn.second, idx, dx),
                0.0) /
            d_clot_params.nb_max;
        double R3 = d_clot_params.Kbb * std::pow(phi_b - 2.0 * z, 2.0) / d_clot_params.nb_max;
        double R4 = d_clot_params.Kaw * w * phi_a / d_clot_params.nb_max;

        const double alpha = R2 + R3 + R4;

        // Stress decay
        double trace = 0.0;
        for (int d = 0; d < NDIM; ++d) trace += (*sig_data)(idx, d);
        const double y_brackets = std::sqrt(2.0 * trace / (z * d_clot_params.S0 + 1.0e-8));
        double beta = d_beta_fcn(y_brackets);

        (*ret_data)(idx) = alpha - beta * z;
    }
} // setDataOnPatch
} // namespace clot
