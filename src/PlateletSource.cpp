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

#include <clot/PlateletSource.h>
#include <clot/app_namespaces.h>

#include <HierarchyDataOpsManager.h>
#include <SAMRAI_config.h>
#include <math.h>

namespace clot
{
/////////////////////////////// PUBLIC ///////////////////////////////////////

PlateletSource::PlateletSource(std::string object_name,
                               Pointer<Variable<NDIM>> phi_u_var,
                               Pointer<Variable<NDIM>> phi_a_var,
                               Pointer<Database> input_db,
                               Pointer<AdvDiffHierarchyIntegrator> adv_diff_hier_integrator)
    : CartGridFunction(std::move(object_name)),
      d_phi_u_var(phi_u_var),
      d_phi_a_var(phi_a_var),
      d_adv_diff_hier_integrator(adv_diff_hier_integrator)
{
    // These need to be changed to the relevant parameters
    // a0 Constants
    d_Kua = input_db->getDouble("Kua");
    d_Kuw = input_db->getDouble("Kuw");

    // scratch index
    auto var_db = VariableDatabase<NDIM>::getDatabase();
    d_pl_scr_idx =
        var_db->registerVariableAndContext(d_phi_a_var, var_db->getContext(d_object_name + "::ScrCtx"), 4 /*ghosts*/);
    return;
} // PlateletSource

bool
PlateletSource::isTimeDependent() const
{
    return true;
} // isTimeDependent

void
PlateletSource::setDataOnPatchHierarchy(const int data_idx,
                                        Pointer<Variable<NDIM>> var,
                                        Pointer<PatchHierarchy<NDIM>> hierarchy,
                                        const double data_time,
                                        const bool initial_time,
                                        const int coarsest_ln_in,
                                        const int finest_ln_in)
{
    const int coarsest_ln = (coarsest_ln_in == -1 ? 0 : coarsest_ln_in);
    const int finest_ln = (finest_ln_in == -1 ? hierarchy->getFinestLevelNumber() : finest_ln_in);
    // Loop over variables.
    std::array<Pointer<Variable<NDIM>>, 2> vars = { d_phi_u_var, d_phi_a_var };
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

    // Fill in ghost cells for scratch index
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        Pointer<PatchLevel<NDIM>> level = hierarchy->getPatchLevel(ln);
        level->allocatePatchData(d_pl_scr_idx, data_time);
    }
    auto var_db = VariableDatabase<NDIM>::getDatabase();
    int phi_a_scr_idx =
        var_db->mapVariableAndContextToIndex(d_phi_a_var, d_adv_diff_hier_integrator->getScratchContext());
    using ITC = HierarchyGhostCellInterpolation::InterpolationTransactionComponent;
    ITC ghost_fill_comp(d_pl_scr_idx,
                        phi_a_scr_idx,
                        "CONSERVATIVE_LINEAR_REFINE",
                        false,
                        "CONSERVATIVE_COARSEN",
                        "LINEAR",
                        false,
                        d_adv_diff_hier_integrator->getPhysicalBcCoefs(d_phi_a_var));
    HierarchyGhostCellInterpolation ghost_fill_op;
    ghost_fill_op.initializeOperatorState(ghost_fill_comp, hierarchy);
    ghost_fill_op.fillData(data_time);

    // Compute the source function
    for (int level_num = coarsest_ln; level_num <= finest_ln; ++level_num)
    {
        setDataOnPatchLevel(data_idx, var, hierarchy->getPatchLevel(level_num), data_time, initial_time);
    }

    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        Pointer<PatchLevel<NDIM>> level = hierarchy->getPatchLevel(ln);
        level->deallocatePatchData(d_pl_scr_idx);
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
PlateletSource::setDataOnPatch(const int data_idx,
                               Pointer<Variable<NDIM>> /*var*/,
                               Pointer<Patch<NDIM>> patch,
                               const double /*data_time*/,
                               const bool initial_time,
                               Pointer<PatchLevel<NDIM>> /*patch_level*/)
{
    Pointer<CellData<NDIM, double>> F_data = patch->getPatchData(data_idx);
    F_data->fillAll(0.0);
    if (initial_time) return;
    Pointer<CartesianPatchGeometry<NDIM>> pgeom = patch->getPatchGeometry();
    const double* const dx = pgeom->getDx();
    const double* const xlow = pgeom->getXLower();

    std::array<VectorNd, 2> bdrys;
    bdrys[0](0) = -2.0;
    bdrys[0](1) = 0.0;
    bdrys[1](0) = 6.0;
    bdrys[1](1) = 2.0;
    // It's apparent that when Baaron wrote this, he intended I use them, so I'll ask about this (I think this relates
    // to **)
    Pointer<CellData<NDIM, double>> phi_a_data = patch->getPatchData(d_pl_scr_idx);
    Pointer<CellData<NDIM, double>> phi_u_data =
        patch->getPatchData(d_phi_u_var, d_adv_diff_hier_integrator->getScratchContext());
    Pointer<CellData<NDIM, double>> w_data = patch->getPatchData(d_w_idx);
    const Box<NDIM>& patch_box = patch->getBox();
    const hier::Index<NDIM>& idx_low = patch_box.lower();
    auto psi_fcn = getKernelAndWidth(d_kernel);
    for (CellIterator<NDIM> ci(patch_box); ci; ci++)
    {
        const CellIndex<NDIM>& idx = ci();
        // Compute source data (relaxation term)
        // double phi_a = (*phi_a_data)(idx);
        double phi_u = (*phi_u_data)(idx);
        double w = (*w_data)(idx);
        // convolve phi_a*psi
        // included w_data as the 4 arg since idk how to have an empty "const CellData<NDIM, double>&" object.
        // const double eta_a =
        //    convolution(1.0, phi_a_data.getPointer(), 0.0, nullptr, psi_fcn.first, psi_fcn.second, idx, dx);
        VectorNd x;
        for (unsigned int d = 0; d < NDIM; ++d)
            x[d] = xlow[d] + dx[d] * (static_cast<double>(idx(d) - idx_low(d)) + 0.5);
        const double eta_a = convolution_mask(
            1.0, phi_a_data.getPointer(), 0.0, nullptr, psi_fcn.first, psi_fcn.second, idx, dx, x, bdrys);
        // Compute the f^a_u
        (*F_data)(idx) = d_sign * (d_Kua * phi_u * eta_a + d_Kuw * w * phi_u); // include f^a_u?
    }
    return;
} // setDataOnPatch
} // namespace clot
//////////////////////////////////////////////////////////////////////////////
