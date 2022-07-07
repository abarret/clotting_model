// ---------------------------------------------------------------------
//
// Copyright (c) 2019 - 2020 by the IBAMR developers
// All rights reserved.
//
// This file is part of IBAMR.
//
// IBAMR is free software and is distributed under the 3-clause BSD
// license. The full text of the license can be found in the file
// COPYRIGHT at the top level directory of IBAMR.
//
// ---------------------------------------------------------------------

#include "ibamr/app_namespaces.h" // IWYU pragma: keep

#include "ibtk/ibtk_utilities.h"

#include "CellData.h"
#include "CellIterator.h"
#include "Patch.h"
#include "tbox/Database.h"

// Local headers
#include "clot/CohesionStressUnactivatedRHS.h"

// Namespace
namespace clot
{
CohesionStressUnactivatedRHS::CohesionStressUnactivatedRHS(Pointer<Variable<NDIM>> phi_u_var,
                                                           Pointer<Variable<NDIM>> phi_a_var,
                                                           Pointer<Variable<NDIM>> z_var,
                                                           const std::string& object_name,
                                                           Pointer<Database> input_db,
                                                           Pointer<AdvDiffHierarchyIntegrator> adv_diff_integrator)
    : CFRelaxationOperator(object_name, input_db),
      d_phi_u_var(phi_u_var),
      d_phi_a_var(phi_a_var),
      d_z_var(z_var),
      d_adv_diff_integrator(adv_diff_integrator)
{
    // Get values from inputdb
    d_c4 = input_db->getDouble("c4");
    // a0 Constants
    d_a0 = input_db->getDouble("a0");
    d_a0w = input_db->getDouble("a0w");
    d_beta_limit = input_db->getDoubleWithDefault("beta_limit", 300.0);
    d_clot_break_x = input_db->getDoubleWithDefault("clot_break_x", 2.25);
    return;
} // Constructor

void
CohesionStressUnactivatedRHS::setDataOnPatchHierarchy(const int data_idx,
                                                      Pointer<Variable<NDIM>> var,
                                                      Pointer<PatchHierarchy<NDIM>> hierarchy,
                                                      const double data_time,
                                                      const bool initial_time,
                                                      const int coarsest_ln_in,
                                                      const int finest_ln_in)
{
    // This implementation atm is identical to ActivatePlatelet Source
    // Loop over variables.
    std::array<Pointer<Variable<NDIM>>, 3> vars = { d_phi_u_var, d_z_var, d_phi_a_var };
    std::map<Pointer<Variable<NDIM>>, bool> scratch_allocated;
    for (const auto local_var : vars)
    {
        // Allocate scratch data when needed.
        VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
        int var_scr_idx = var_db->mapVariableAndContextToIndex(local_var, d_adv_diff_integrator->getScratchContext());
        scratch_allocated[local_var] = d_adv_diff_integrator->isAllocatedPatchData(var_scr_idx);
        if (!scratch_allocated[local_var])
        {
            d_adv_diff_integrator->allocatePatchData(var_scr_idx, data_time);
        }

        // Communicate ghost-cell data.
        if (!initial_time)
        {
            // Determine values at correct timestep.
            int var_current_idx =
                var_db->mapVariableAndContextToIndex(local_var, d_adv_diff_integrator->getCurrentContext());
            int var_new_idx = var_db->mapVariableAndContextToIndex(local_var, d_adv_diff_integrator->getNewContext());
            const bool var_new_is_allocated = d_adv_diff_integrator->isAllocatedPatchData(var_new_idx);
            HierarchyDataOpsManager<NDIM>* hier_data_ops_manager = HierarchyDataOpsManager<NDIM>::getManager();
            Pointer<HierarchyDataOpsReal<NDIM, double>> hier_cc_data_ops =
                hier_data_ops_manager->getOperationsDouble(local_var, hierarchy, /*get_unique*/ true);
            if (d_adv_diff_integrator->getCurrentCycleNumber() <= 0 || !var_new_is_allocated)
            {
                // Copy values from step n.
                hier_cc_data_ops->copyData(var_scr_idx, var_current_idx);
            }
            else
            {
                // Values at step n+1 is filled, so we use values at step n+1/2
#if !defined(NDEBUG)
                TBOX_ASSERT(d_adv_diff_integrator->getCurrentCycleNumber() > 0);
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
                                     d_adv_diff_integrator->getPhysicalBcCoefs(local_var));
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
                var_db->mapVariableAndContextToIndex(local_var, d_adv_diff_integrator->getScratchContext());
            d_adv_diff_integrator->deallocatePatchData(var_scr_idx);
        }
    }
    return;
} // setDataOnPatchHierarchy

void
CohesionStressUnactivatedRHS::setDataOnPatch(const int data_idx,
                                             Pointer<Variable<NDIM>> /*var*/,
                                             Pointer<Patch<NDIM>> patch,
                                             const double data_time,
                                             const bool initial_time,
                                             Pointer<PatchLevel<NDIM>> /*patch_level*/)
{
    const Box<NDIM>& patch_box = patch->getBox();
    Pointer<CartesianPatchGeometry<NDIM>> pgeom = patch->getPatchGeometry();
    const double* const dx = pgeom->getDx();
    const double* const xlow = pgeom->getXLower();
    const hier::Index<NDIM>& idx_low = patch_box.lower();
    // SET THE INDICES!!!!
    Pointer<CellData<NDIM, double>> ret_data = patch->getPatchData(data_idx);
    Pointer<CellData<NDIM, double>> in_data = patch->getPatchData(d_W_cc_idx); // this is sigma?
    Pointer<CellData<NDIM, double>> phi_a_data =
        patch->getPatchData(d_phi_a_var, d_adv_diff_integrator->getScratchContext());
    Pointer<CellData<NDIM, double>> phi_u_data =
        patch->getPatchData(d_phi_u_var, d_adv_diff_integrator->getScratchContext());
    Pointer<CellData<NDIM, double>> z_data = patch->getPatchData(d_z_var, d_adv_diff_integrator->getScratchContext());
    // stress tensor and wall sites
    Pointer<CellData<NDIM, double>> w_data = patch->getPatchData(d_w_idx);
    ret_data->fillAll(0.0);
    if (initial_time) return;
    for (CellIterator<NDIM> i(patch_box); i; i++)
    {
        const CellIndex<NDIM>& idx = i();
        VectorNd x;
        for (int d = 0; d < NDIM; ++d) x[d] = xlow[d] + dx[d] * (static_cast<double>(idx(d) - idx_low(d)) + 0.5);
        // Compute the ODE terms for stress here
        // C4 * alpha * I - beta * sigma
        // compute R2 and R4 -> alpha
        // where beta = beta(<y>)
        // where <y> = Tr(sig) / z

        // Index variables
        double phi_a = (*phi_a_data)(idx);
        double phi_u = (*phi_u_data)(idx);
        double w = (*w_data)(idx);
        double z = (*z_data)(idx);
        // Compute the source terms
        const double R2 = d_a0 * phi_a * phi_a;
        const double R4 = d_a0w * w * phi_a;
        const double alpha = R2 + R4;
#if (NDIM == 2)
        const double trace = (*in_data)(idx, 0) + (*in_data)(idx, 1);
        const double y_brackets = trace / (z + 1.0e-8);
        double beta = d_beta_fcn(y_brackets);
        if (x[0] > d_clot_break_x) beta = d_beta_limit;

        (*ret_data)(idx, 0) = d_c4 * alpha - beta * (*in_data)(idx, 0);
        (*ret_data)(idx, 1) = d_c4 * alpha - beta * (*in_data)(idx, 1);
        (*ret_data)(idx, 2) = -beta * (*in_data)(idx, 2);
#endif
#if (NDIM == 3)
        const double trace = (*in_data)(idx, 0) + (*in_data)(idx, 1) + (*in_data)(idx, 2);
        const double y_brackets = trace / (z + 1.0e-8);
        const double beta = d_beta_fcn(y_brackets);

        (*ret_data)(idx, 0) = d_c4 * alpha - beta * (*in_data)(idx, 0);
        (*ret_data)(idx, 1) = d_c4 * alpha - beta * (*in_data)(idx, 1);
        (*ret_data)(idx, 2) = d_c4 * alpha - beta * (*in_data)(idx, 2);
        (*ret_data)(idx, 3) = -beta * (*in_data)(idx, 3);
        (*ret_data)(idx, 4) = -beta * (*in_data)(idx, 4);
        (*ret_data)(idx, 5) = -beta * (*in_data)(idx, 5);
#endif
    }
} // setDataOnPatch
} // namespace clot
