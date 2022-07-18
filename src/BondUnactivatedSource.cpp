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

#include <clot/BondUnactivatedSource.h>
#include <clot/app_namespaces.h>

#include <HierarchyDataOpsManager.h>
#include <SAMRAI_config.h>
#include <math.h>

namespace clot
{
/////////////////////////////// PUBLIC ///////////////////////////////////////
BondUnactivatedSource::BondUnactivatedSource(Pointer<hier::Variable<NDIM>> phi_u_var,
                                             Pointer<hier::Variable<NDIM>> phi_a_var,
                                             Pointer<hier::Variable<NDIM>> z_var,
                                             Pointer<hier::Variable<NDIM>> sig_var,
                                             Pointer<Database> input_db,
                                             Pointer<AdvDiffHierarchyIntegrator> adv_diff_integrator,
                                             Pointer<AdvDiffHierarchyIntegrator> sb_adv_diff_integrator)
    : d_phi_u_var(phi_u_var),
      d_phi_a_var(phi_a_var),
      d_z_var(z_var),
      d_sig_var(sig_var),
      d_adv_diff_integrator(adv_diff_integrator),
      d_sb_adv_diff_integrator(sb_adv_diff_integrator)
{
    // init db vars
    d_a0 = input_db->getDouble("a0");
    d_a0w = input_db->getDouble("a0w");
    d_beta_limit = input_db->getDoubleWithDefault("beta_limit", 300.0);
    d_clot_break_x = input_db->getDoubleWithDefault("clot_break_x", 2.25);
} // BondUnactivatedSource

bool
BondUnactivatedSource::isTimeDependent() const
{
    return true;
} // isTimeDependent

void
BondUnactivatedSource::setDataOnPatchHierarchy(const int data_idx,
                                               Pointer<Variable<NDIM>> var,
                                               Pointer<PatchHierarchy<NDIM>> hierarchy,
                                               const double data_time,
                                               const bool initial_time,
                                               const int coarsest_ln_in,
                                               const int finest_ln_in)
{
    // This implementation atm is identical to ActivatePlatelet Source
    // Loop over variables.
    std::array<std::pair<Pointer<Variable<NDIM>>, Pointer<AdvDiffHierarchyIntegrator>>, 4> var_hier_pairs = {
        std::make_pair(d_phi_u_var, d_sb_adv_diff_integrator),
        std::make_pair(d_z_var, d_sb_adv_diff_integrator),
        std::make_pair(d_phi_a_var, d_sb_adv_diff_integrator),
        std::make_pair(d_sig_var, d_adv_diff_integrator)
    };
    std::map<Pointer<Variable<NDIM>>, bool> scratch_allocated;
    for (const auto var_hier_pair : var_hier_pairs)
    {
        Pointer<Variable<NDIM>> local_var = var_hier_pair.first;
        Pointer<AdvDiffHierarchyIntegrator> integrator = var_hier_pair.second;
        // Allocate scratch data when needed.
        VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
        int var_scr_idx = var_db->mapVariableAndContextToIndex(local_var, integrator->getScratchContext());
        scratch_allocated[local_var] = integrator->isAllocatedPatchData(var_scr_idx);
        if (!scratch_allocated[local_var])
        {
            integrator->allocatePatchData(var_scr_idx, data_time);
        }

        // Communicate ghost-cell data.
        if (!initial_time)
        {
            // Determine values at correct timestep.
            int var_current_idx = var_db->mapVariableAndContextToIndex(local_var, integrator->getCurrentContext());
            int var_new_idx = var_db->mapVariableAndContextToIndex(local_var, integrator->getNewContext());
            const bool var_new_is_allocated = integrator->isAllocatedPatchData(var_new_idx);
            HierarchyDataOpsManager<NDIM>* hier_data_ops_manager = HierarchyDataOpsManager<NDIM>::getManager();
            Pointer<HierarchyDataOpsReal<NDIM, double>> hier_cc_data_ops =
                hier_data_ops_manager->getOperationsDouble(local_var, hierarchy, /*get_unique*/ true);
            if (integrator->getCurrentCycleNumber() == 0 || !var_new_is_allocated)
            {
                // Copy values from step n.
                hier_cc_data_ops->copyData(var_scr_idx, var_current_idx);
            }
            else
            {
                // Values at step n+1 is filled, so we use values at step n+1/2
#if !defined(NDEBUG)
                TBOX_ASSERT(integrator->getCurrentCycleNumber() > 0);
#endif
                hier_cc_data_ops->linearSum(var_scr_idx, 0.5, var_current_idx, 0.5, var_new_idx);
            }

            // Communicate ghost cell values
            using ITC = HierarchyGhostCellInterpolation::InterpolationTransactionComponent;
            ITC ghost_fill_component(var_scr_idx,
                                     "CONSERVATIVE_LINEAR_REFINE",
                                     false,
                                     "CONSERVATIVE_COARSEN",
                                     "LINEAR",
                                     false,
                                     integrator->getPhysicalBcCoefs(local_var));
            HierarchyGhostCellInterpolation ghost_fill_op;
            ghost_fill_op.initializeOperatorState(ghost_fill_component, hierarchy);
            ghost_fill_op.fillData(data_time);
        }
    }

    // Compute the source function
    const int coarsest_ln = (coarsest_ln_in == -1 ? 0 : coarsest_ln_in);
    const int finest_ln = (finest_ln_in == -1 ? hierarchy->getFinestLevelNumber() : finest_ln_in);
    for (int level_num = coarsest_ln; level_num <= finest_ln; ++level_num)
    {
        setDataOnPatchLevel(data_idx, var, hierarchy->getPatchLevel(level_num), data_time, initial_time);
    }

    // Deallocate scratch data when needed.
    for (const auto& var_hier_pair : var_hier_pairs)
    {
        Pointer<Variable<NDIM>> local_var = var_hier_pair.first;
        Pointer<AdvDiffHierarchyIntegrator> integrator = var_hier_pair.second;
        if (!scratch_allocated[local_var])
        {
            auto var_db = VariableDatabase<NDIM>::getDatabase();
            int var_scr_idx = var_db->mapVariableAndContextToIndex(local_var, integrator->getScratchContext());
            d_sb_adv_diff_integrator->deallocatePatchData(var_scr_idx);
        }
    }
    return;
} // setDataOnPatchHierarchy

void
BondUnactivatedSource::setDataOnPatch(const int data_idx,
                                      Pointer<Variable<NDIM>> /*var*/,
                                      Pointer<Patch<NDIM>> patch,
                                      const double /*data_time*/,
                                      const bool initial_time,
                                      Pointer<PatchLevel<NDIM>> /*patch_level*/)
{
    // Computes alpha - beta * z
    Pointer<CellData<NDIM, double>> bond_data = patch->getPatchData(data_idx);
    bond_data->fillAll(0.0);
    if (initial_time) return;
    // grab cell data for variables
    Pointer<CellData<NDIM, double>> phi_a_data =
        patch->getPatchData(d_phi_a_var, d_sb_adv_diff_integrator->getScratchContext());
    Pointer<CellData<NDIM, double>> phi_u_data =
        patch->getPatchData(d_phi_u_var, d_sb_adv_diff_integrator->getScratchContext());
    Pointer<CellData<NDIM, double>> z_data =
        patch->getPatchData(d_z_var, d_sb_adv_diff_integrator->getScratchContext());
    // stress tensor and wall sites
    Pointer<CellData<NDIM, double>> sig_data =
        patch->getPatchData(d_sig_var, d_adv_diff_integrator->getScratchContext());
    Pointer<CellData<NDIM, double>> w_data = patch->getPatchData(d_w_idx);
    const Box<NDIM>& patch_box = patch->getBox();
    Pointer<CartesianPatchGeometry<NDIM>> pgeom = patch->getPatchGeometry();
    const double* const dx = pgeom->getDx();
    const double* const xlow = pgeom->getXLower();
    const hier::Index<NDIM>& idx_low = patch_box.lower();
    // begin cell loop
    for (CellIterator<NDIM> ci(patch_box); ci; ci++)
    {
        // grab the index
        const CellIndex<NDIM>& idx = ci();
        VectorNd x;
        for (int d = 0; d < NDIM; ++d) x[d] = xlow[d] + dx[d] * (static_cast<double>(idx(d) - idx_low(d)) + 0.5);
        // grab the var values w/ the index
        double phi_u = (*phi_u_data)(idx);
        double phi_a = (*phi_a_data)(idx);
        double w = (*w_data)(idx);
        double z = (*z_data)(idx);
        // Compute alpha data
        double alpha = (d_a0 * phi_a * phi_a) + (d_a0w * w * phi_a);
#if (NDIM == 2)
        const double trace = (*sig_data)(idx, 0) + (*sig_data)(idx, 1);
        const double y_brackets = trace / (z + 1.0e-8);
        double beta = d_beta_fcn(y_brackets);
        if (x[0] > d_clot_break_x) beta = d_beta_limit;
#endif
#if (NDIM == 3)
        const double trace = (*sig_data)(idx, 0) + (*sig_data)(idx, 1) + (*sig_data)(idx, 2);
        const double y_brackets = trace / (z + 1.0e-8);
        double beta = d_beta_fcn(y_brackets);
#endif
        (*bond_data)(idx) = alpha - beta * z;
    }
} // setDataOnPatch
} // namespace clot
