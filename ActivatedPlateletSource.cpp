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

#include <ibamr/app_namespaces.h>

#include <HierarchyDataOpsManager.h>
#include <SAMRAI_config.h>
#include <math.h>

// Local includes
#include "ActivatedPlateletSource.h"

/////////////////////////////// PUBLIC ///////////////////////////////////////

ActivatedPlateletSource::ActivatedPlateletSource(Pointer<Variable<NDIM>> pl_n_var,
                                                 Pointer<Variable<NDIM>> c_var,
                                                 Pointer<AdvDiffHierarchyIntegrator> adv_diff_hier_integrator)
    : d_pl_n_var(pl_n_var), d_c_var(c_var), d_adv_diff_hier_integrator(adv_diff_hier_integrator)
{
    // intentionally blank
    d_ct = input_db->getDouble("ct"); // change the get strings to fix whatever the true label is
    d_r0 = input_db->getDouble("r0");
    d_k = input_db->getDouble("k");
    return;
} // ActivatedPlateletSource

bool
ActivatedPlateletSource::isTimeDependent() const
{
    return true;
} // isTimeDependent

void
ActivatedPlateletSource::setDataOnPatchHierarchy(const int data_idx,
                                                 Pointer<Variable<NDIM>> var,
                                                 Pointer<PatchHierarchy<NDIM>> hierarchy,
                                                 const double data_time,
                                                 const bool initial_time,
                                                 const int coarsest_ln_in,
                                                 const int finest_ln_in)
{
    // Loop over variables.
    std::array<Pointer<Variable<NDIM>>, 2> vars = { d_pl_n_var, d_c_var };
    std::map<Pointer<Variable<NDIM>>, bool> scratch_allocated;
    for (const auto local_var : vars)
    {
        // Allocate scratch data when needed.
        VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
        int var_scr_idx =
            var_db->mapVariableAndContextToIndex(local_var, d_adv_diff_hier_integrator->getScratchContext());
        scratch_allocated[local_var] = d_adv_diff_hier_integrator->isAllocatedPatchData(var_scr_idx);
        if (!scratch_allocated[local_var])
        {
            d_adv_diff_hier_integrator->allocatePatchData(var_scr_idx, data_time);
        }

        // Communicate ghost-cell data.
        if (!initial_time)
        {
            // Determine values at correct timestep.
            int var_current_idx =
                var_db->mapVariableAndContextToIndex(local_var, d_adv_diff_hier_integrator->getCurrentContext());
            int var_new_idx =
                var_db->mapVariableAndContextToIndex(local_var, d_adv_diff_hier_integrator->getNewContext());
            const bool var_new_is_allocated = d_adv_diff_hier_integrator->isAllocatedPatchData(var_new_idx);
            HierarchyDataOpsManager<NDIM>* hier_data_ops_manager = HierarchyDataOpsManager<NDIM>::getManager();
            Pointer<HierarchyDataOpsReal<NDIM, double>> hier_cc_data_ops =
                hier_data_ops_manager->getOperationsDouble(local_var, hierarchy, /*get_unique*/ true);
            if (d_adv_diff_hier_integrator->getCurrentCycleNumber() == 0 || !var_new_is_allocated)
            {
                // Copy values from step n.
                hier_cc_data_ops->copyData(var_scr_idx, var_current_idx);
            }
            else
            {
                // Values at step n+1 is filled, so we use values at step n+1/2
#if !defined(NDEBUG)
                TBOX_ASSERT(d_adv_diff_hier_integrator->getCurrentCycleNumber() > 0);
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
                                     d_adv_diff_hier_integrator->getPhysicalBcCoefs(local_var));
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
    for (const auto& local_var : vars)
    {
        if (!scratch_allocated[local_var])
        {
            auto var_db = VariableDatabase<NDIM>::getDatabase();
            int var_scr_idx =
                var_db->mapVariableAndContextToIndex(local_var, d_adv_diff_hier_integrator->getScratchContext());
            d_adv_diff_hier_integrator->deallocatePatchData(var_scr_idx);
        }
    }
    return;
} // setDataOnPatchHierarchy

void
ActivatedPlateletSource::setDataOnPatch(const int data_idx,
                                        Pointer<Variable<NDIM>> /*var*/, // So why do some params have no name?
                                        Pointer<Patch<NDIM>> patch,
                                        const double /*data_time*/,
                                        const bool initial_time,
                                        Pointer<PatchLevel<NDIM>> /*patch_level*/)
{
    Pointer<CellData<NDIM, double>> F_data = patch->getPatchData(data_idx);
    F_data->fillAll(0.0);
    if (initial_time) return;
    Pointer<CellData<NDIM, double>> pl_n_data =
        patch->getPatchData(d_pl_n_var, d_adv_diff_hier_integrator->getScratchContext());
    Pointer<CellData<NDIM, double>> c_data =
        patch->getPatchData(d_c_var, d_adv_diff_hier_integrator->getScratchContext());
    const Box<NDIM>& patch_box = patch->getBox();
    for (CellIterator<NDIM> ci(patch_box); ci; ci++)
    {
        const CellIndex<NDIM>& idx = ci();
        // Compute source data and fill in F_data
        // I assume the values to be filled is that of phi_n, and that H is defined at this point
        double d_c = (*c_data)(idx);
        double d_pl_n = (*pl_n_data)(idx);
        // phi_n+1 = phi_n * exp(-k(R(c)+R_w(w))) (at a specific index)
        // Here, R(c) = R_0 * H(c-c_t), R_0 base reaction & c_t concentration threshold. H is smooth Heaviside function.
        // we expect d_r0 and d_ct to be set correctly by now.
        double Rc = d_r0 * Heaviside(d_c - d_ct);
        // typically R_w(w) will be 0. But how do we get this info when we need it?
        (*F_data)(idx) = d_pl_n * std::exp(-d_k*Rc);
    }
    return;
} // setDataOnPatch


ActivatedPlateletSource::setDataOnPatchLevel(const int data_idx,
                            Pointer<Variable<NDIM>> var /* var */,
                            Pointer<PatchLevel<NDIM>> p_level /* patch level */,
                            const double data_time,
                            const bool initial_time = false)
{
    // Much ado about patch levels, though I don't fully understand their nuances even now
    // I will discuss this with Aaron tomorrow.
} // setDataOnPatchLevel
//////////////////////////////////////////////////////////////////////////////
