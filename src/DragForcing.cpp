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

#include "clot/DragForcing.h"

#include <ADS/app_namespaces.h>

#include <ibtk/HierarchyGhostCellInterpolation.h>

namespace clot
{
/////////////////////////////// PUBLIC ///////////////////////////////////////

DragForcing::DragForcing(std::string object_name, Pointer<Database> input_db) : CartGridFunction(std::move(object_name))
{
    d_C3 = input_db->getDouble("c3");
    d_vol_pl = input_db->getDouble("vol_pl");
    return;
} // DragForcing

DragForcing::~DragForcing()
{
    // intentionally blank
    return;
} // ~DragForcing

bool
DragForcing::isTimeDependent() const
{
    return true;
} // isTimeDependent

void
DragForcing::setDataOnPatchHierarchy(const int data_idx,
                                     Pointer<Variable<NDIM>> var,
                                     Pointer<PatchHierarchy<NDIM>> hierarchy,
                                     const double data_time,
                                     const bool initial_time,
                                     int coarsest_ln,
                                     int finest_ln)
{
    coarsest_ln = coarsest_ln == -1 ? 0 : coarsest_ln;
    finest_ln = finest_ln == -1 ? hierarchy->getFinestLevelNumber() : finest_ln;

    // Need to fill in scratch indices with ghost data
    using ITC = HierarchyGhostCellInterpolation::InterpolationTransactionComponent;
    std::vector<ITC> ghost_comps(3);

    HierarchyCellDataOpsReal<NDIM, double> hier_cc_data_ops(hierarchy, coarsest_ln, finest_ln);
    auto var_db = VariableDatabase<NDIM>::getDatabase();
    const int phi_scr_idx = var_db->mapVariableAndContextToIndex(d_phi_var, d_phi_integrator->getScratchContext());
    const int phi_cur_idx = var_db->mapVariableAndContextToIndex(d_phi_var, d_phi_integrator->getCurrentContext());
    const int phi_new_idx = var_db->mapVariableAndContextToIndex(d_phi_var, d_phi_integrator->getNewContext());
    const bool phi_new_allocated = d_phi_integrator->isAllocatedPatchData(phi_new_idx, coarsest_ln, finest_ln);
    const bool deallocate_phi_scr = !d_phi_integrator->isAllocatedPatchData(phi_scr_idx, coarsest_ln, finest_ln);
    if (deallocate_phi_scr) d_phi_integrator->allocatePatchData(phi_scr_idx, data_time, coarsest_ln, finest_ln);

    if (phi_new_allocated)
    {
        // Use the average of phi_cur and phi_new
        hier_cc_data_ops.linearSum(phi_scr_idx, 0.5, phi_cur_idx, 0.5, phi_new_idx);
    }
    else
    {
        hier_cc_data_ops.copyData(phi_scr_idx, phi_cur_idx);
    }

    // Now create the ITC
    ghost_comps[0] = ITC(phi_scr_idx,
                         "CONSERVATIVE_LINEAR_REFINE",
                         false,
                         "NONE",
                         "LINEAR",
                         false,
                         d_phi_integrator->getPhysicalBcCoefs(d_phi_var));

    const int ub_scr_idx = var_db->mapVariableAndContextToIndex(d_ub_var, d_ub_integrator->getScratchContext());
    const int ub_cur_idx = var_db->mapVariableAndContextToIndex(d_ub_var, d_ub_integrator->getCurrentContext());
    const int ub_new_idx = var_db->mapVariableAndContextToIndex(d_ub_var, d_ub_integrator->getNewContext());
    const bool sig_ub_allocated = d_ub_integrator->isAllocatedPatchData(ub_new_idx, coarsest_ln, finest_ln);
    const bool deallocate_ub_scr = !d_ub_integrator->isAllocatedPatchData(ub_scr_idx, coarsest_ln, finest_ln);
    if (deallocate_ub_scr) d_ub_integrator->allocatePatchData(ub_scr_idx, data_time, coarsest_ln, finest_ln);

    if (sig_ub_allocated)
    {
        // Use the average of phi_cur and phi_new
        hier_cc_data_ops.linearSum(ub_scr_idx, 0.5, ub_cur_idx, 0.5, ub_new_idx);
    }
    else
    {
        hier_cc_data_ops.copyData(ub_scr_idx, ub_cur_idx);
    }

    // Now create the ITC
    ghost_comps[1] = ITC(ub_scr_idx,
                         "CONSERVATIVE_LINEAR_REFINE",
                         false,
                         "NONE",
                         "LINEAR",
                         false,
                         d_ub_integrator->getPhysicalBcCoefs(d_ub_var));

    HierarchyFaceDataOpsReal<NDIM, double> hier_fc_data_ops(hierarchy, coarsest_ln, finest_ln);
    const int uf_scr_idx = var_db->mapVariableAndContextToIndex(d_uf_var, d_ins_integrator->getScratchContext());
    const int uf_cur_idx = var_db->mapVariableAndContextToIndex(d_uf_var, d_ins_integrator->getCurrentContext());
    const int uf_new_idx = var_db->mapVariableAndContextToIndex(d_uf_var, d_ins_integrator->getNewContext());
    const bool uf_new_allocated = d_ins_integrator->isAllocatedPatchData(uf_new_idx, coarsest_ln, finest_ln);
    const bool deallocate_uf_scr = !d_ins_integrator->isAllocatedPatchData(uf_scr_idx, coarsest_ln, finest_ln);
    if (deallocate_uf_scr) d_ins_integrator->allocatePatchData(uf_scr_idx, data_time, coarsest_ln, finest_ln);

    if (uf_new_allocated)
    {
        // Use the average of phi_cur and phi_new
        hier_fc_data_ops.linearSum(uf_scr_idx, 0.5, uf_cur_idx, 0.5, uf_new_idx);
    }
    else
    {
        hier_fc_data_ops.copyData(uf_scr_idx, uf_cur_idx);
    }

    // Now create the ITC
    ghost_comps[2] = ITC(uf_scr_idx,
                         "LINEAR_REFINE",
                         false,
                         "NONE",
                         "LINEAR",
                         false,
                         d_ins_integrator->getIntermediateVelocityBoundaryConditions());

    // Now that ghost cells are filled, fill data
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        setDataOnPatchLevel(data_idx, var, hierarchy->getPatchLevel(ln), data_time, initial_time);
    }

    // Deallocate patch data if needed
    if (deallocate_phi_scr) d_phi_integrator->deallocatePatchData(phi_scr_idx, coarsest_ln, finest_ln);
    if (deallocate_ub_scr) d_ub_integrator->deallocatePatchData(ub_scr_idx, coarsest_ln, finest_ln);
    if (deallocate_uf_scr) d_ins_integrator->deallocatePatchData(uf_scr_idx, coarsest_ln, finest_ln);
} // setDataOnPatchHierarchy

void
DragForcing::setDataOnPatch(const int data_idx,
                            Pointer<Variable<NDIM>> /*var*/,
                            Pointer<Patch<NDIM>> patch,
                            const double /*data_time*/,
                            const bool initial_time,
                            Pointer<PatchLevel<NDIM>> /*patch_level*/)
{
    Pointer<SideData<NDIM, double>> F_data = patch->getPatchData(data_idx);
    F_data->fillAll(0.0);
    if (initial_time) return;
    Pointer<SideData<NDIM, double>> uf_data = patch->getPatchData(d_uf_var, d_ins_integrator->getScratchContext());
    Pointer<SideData<NDIM, double>> ub_data = patch->getPatchData(d_ub_var, d_ub_integrator->getScratchContext());
    Pointer<CellData<NDIM, double>> phi_data = patch->getPatchData(d_phi_var, d_phi_integrator->getScratchContext());
    const Box<NDIM>& patch_box = patch->getBox();
    for (int axis = 0; axis < NDIM; ++axis)
    {
        for (SideIterator<NDIM> si(patch_box, axis); si; si++)
        {
            const SideIndex<NDIM>& idx = si();
            // We need to be careful on calculating phi near the wall.
            double th_b = d_vol_pl * 0.5 * ((*phi_data)(idx.toCell(1)) + (*phi_data)(idx.toCell(0)));
            double xi = d_C3 * th_b * th_b / std::pow(1.0 - th_b, 3.0);
            (*F_data)(idx) = xi * ((*ub_data)(idx) - (*uf_data)(idx));
        }
    }
    return;
} // setDataOnPatch
} // namespace clot
//////////////////////////////////////////////////////////////////////////////
