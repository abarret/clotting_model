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

#include <clot/BoundExtraStressForcing.h>
#include <clot/app_namespaces.h>

#include <HierarchyDataOpsManager.h>
#include <SAMRAI_config.h>
#include <math.h>

namespace clot
{
/////////////////////////////// PUBLIC ///////////////////////////////////////
BoundExtraStressForcing::BoundExtraStressForcing(std::string object_name, BoundClotParams clot_params)
    : CartGridFunction(std::move(object_name)),
      d_clot_params(std::move(clot_params)),
      d_extra_stress_var(new CellVariable<NDIM, double>(d_object_name + "::Var"))
{
    // Create variable for the extra stress
    // Note it only has depth of one because it's a multiple of the identity matrix.
    auto var_db = VariableDatabase<NDIM>::getDatabase();
    d_extra_stress_idx =
        var_db->registerVariableAndContext(d_extra_stress_var, var_db->getContext(d_object_name + "::ctx"), 1);
} // BoundExtraStressForcing

bool
BoundExtraStressForcing::isTimeDependent() const
{
    return true;
} // isTimeDependent

void
BoundExtraStressForcing::setDataOnPatchHierarchy(const int data_idx,
                                                 Pointer<Variable<NDIM>> var,
                                                 Pointer<PatchHierarchy<NDIM>> hierarchy,
                                                 const double data_time,
                                                 const bool initial_time,
                                                 int coarsest_ln,
                                                 int finest_ln)
{
    coarsest_ln = (coarsest_ln == -1 ? 0 : coarsest_ln);
    finest_ln = (finest_ln == -1 ? hierarchy->getFinestLevelNumber() : finest_ln);
    // Allocate local scratch data
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        Pointer<PatchLevel<NDIM>> level = hierarchy->getPatchLevel(ln);
        level->allocatePatchData(d_extra_stress_idx, data_time);
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
        }
    }

    // Find the multiple of the identity matrix...
    findFactor(hierarchy, coarsest_ln, finest_ln);

    // Now fill in ghost cells
    using ITC = HierarchyGhostCellInterpolation::InterpolationTransactionComponent;
    std::vector<ITC> ghost_cell_comps(1);

    ghost_cell_comps[0] =
        ITC(d_extra_stress_idx, "CONSERVATIVE_LINEAR_REFINE", false, "NONE", "LINEAR", false, nullptr);

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
    // Also deallocate local scratch data
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        Pointer<PatchLevel<NDIM>> level = hierarchy->getPatchLevel(ln);
        level->deallocatePatchData(d_extra_stress_idx);
    }
    return;
} // setDataOnPatchHierarchy

void
BoundExtraStressForcing::setDataOnPatch(const int data_idx,
                                        Pointer<Variable<NDIM>> /*var*/,
                                        Pointer<Patch<NDIM>> patch,
                                        const double /*data_time*/,
                                        const bool initial_time,
                                        Pointer<PatchLevel<NDIM>> /*patch_level*/)
{
    Pointer<SideData<NDIM, double>> F_data = patch->getPatchData(data_idx);
    F_data->fillAll(0.0);
    if (initial_time) return;

    // Patch information
    const Box<NDIM>& patch_box = patch->getBox();
    Pointer<CartesianPatchGeometry<NDIM>> pgeom = patch->getPatchGeometry();
    const double* const dx = pgeom->getDx();

    Pointer<CellData<NDIM, double>> extra_data = patch->getPatchData(d_extra_stress_idx);
    for (int axis = 0; axis < NDIM; ++axis)
    {
        for (SideIterator<NDIM> si(patch_box, axis); si; si++)
        {
            const SideIndex<NDIM> idx = si();
            (*F_data)(idx) = ((*extra_data)(idx.toCell(1)) - (*extra_data)(idx.toCell(0))) / dx[axis];
        }
    }
} // setDataOnPatch

void
BoundExtraStressForcing::findFactor(Pointer<PatchHierarchy<NDIM>> hierarchy, const int coarsest_ln, const int finest_ln)
{
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        Pointer<PatchLevel<NDIM>> level = hierarchy->getPatchLevel(ln);
        for (PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            Pointer<Patch<NDIM>> patch = level->getPatch(p());
            Pointer<CellData<NDIM, double>> bond_data = patch->getPatchData(
                d_var_integrator_pairs[0].first, d_var_integrator_pairs[0].second->getScratchContext());
            Pointer<CellData<NDIM, double>> sig_data = patch->getPatchData(
                d_var_integrator_pairs[1].first, d_var_integrator_pairs[1].second->getScratchContext());
            Pointer<CellData<NDIM, double>> extra_data = patch->getPatchData(d_extra_stress_idx);
            for (CellIterator<NDIM> ci(patch->getBox()); ci; ci++)
            {
                const CellIndex<NDIM>& idx = ci();
                const double z = (*bond_data)(idx);
                double tr = 0.0;
                for (int d = 0; d < NDIM; ++d) tr += (*sig_data)(idx, d);
                (*extra_data)(idx) =
                    z * d_clot_params.S0 * d_clot_params.R0 * std::sqrt(2.0 * tr / (d_clot_params.S0 * z + 1.0e-8));
            }
        }
    }
}
} // namespace clot
