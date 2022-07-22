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

#include <clot/BoundPlateletSource.h>
#include <clot/app_namespaces.h>

#include <boost/math/tools/roots.hpp>

#include <HierarchyDataOpsManager.h>
#include <SAMRAI_config.h>
#include <math.h>

namespace clot
{
/////////////////////////////// PUBLIC ///////////////////////////////////////

BoundPlateletSource::BoundPlateletSource(std::string object_name,
                                         Pointer<Database> input_db)
    : CartGridFunction(std::move(object_name))
{
    // Parameters
    d_Kab = input_db->getDouble("kab");
    d_Kaw = input_db->getDouble("kaw");
    d_nb_max = input_db->getDouble("nb_max");
    d_nw_max = input_db->getDouble("nw_max");
    return;
} // BoundPlateletSource

bool
BoundPlateletSource::isTimeDependent() const
{
    return true;
} // isTimeDependent

void
BoundPlateletSource::setDataOnPatchHierarchy(const int data_idx,
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
    std::map<Pointer<Variable<NDIM>>, bool> scratch_allocated;
    for (const auto& var_integrator_pair : d_var_integrator_pairs)
    {
        Pointer<Variable<NDIM>> var = var_integrator_pair.first;
        Pointer<HierarchyIntegrator> integrator = var_integrator_pair.second;
        // Allocate scratch data when needed
        VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
        int var_scr_idx = var_db->mapVariableAndContextToIndex(var, integrator->getScratchContext());
        scratch_allocated[var] = integrator->isAllocatedPatchData(var_scr_idx);
        if (!scratch_allocated.at(var)) integrator->allocatePatchData(var_scr_idx, data_time);
        if (!initial_time)
        {
            // Determine values at correct time
            int var_cur_idx = var_db->mapVariableAndContextToIndex(var, integrator->getCurrentContext());
            int var_new_idx = var_db->mapVariableAndContextToIndex(var, integrator->getNewContext());
            const bool var_new_is_allocated = integrator->isAllocatedPatchData(var_new_idx);
            auto* hier_data_ops_manager = HierarchyDataOpsManager<NDIM>::getManager();
            Pointer<HierarchyDataOpsReal<NDIM, double>> hier_cc_data_ops =
                hier_data_ops_manager->getOperationsDouble(var, hierarchy, /*get_unique*/ true);
            if (integrator->getCurrentCycleNumber() == 0 || !var_new_is_allocated)
                hier_cc_data_ops->copyData(var_scr_idx, var_cur_idx);
            else
                hier_cc_data_ops->linearSum(var_scr_idx, 0.5, var_cur_idx, 0.5, var_new_idx);

            // Do we need to fill ghost cells?
        }
    }

    // Fill ghost cells for values that require it
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        Pointer<PatchLevel<NDIM>> level = hierarchy->getPatchLevel(ln);
        level->allocatePatchData(d_phi_b_scr_idx, data_time);
        level->allocatePatchData(d_bond_scr_idx, data_time);
    }
    using ITC = HierarchyGhostCellInterpolation::InterpolationTransactionComponent;
    std::vector<ITC> ghost_fill_comps(2);

    Pointer<Variable<NDIM>> phi_b_var = d_var_integrator_pairs[0].first;
    Pointer<AdvDiffHierarchyIntegrator> phi_b_integrator = d_var_integrator_pairs[0].second;
    auto var_db = VariableDatabase<NDIM>::getDatabase();
    int phi_b_scr_idx = var_db->mapVariableAndContextToIndex(phi_b_var, phi_b_integrator->getScratchContext());
    ghost_fill_comps[0] = ITC(d_phi_b_scr_idx,
                              phi_b_scr_idx,
                              "CONSERVATIVE_LINEAR_REFINE",
                              false,
                              "NONE",
                              "LINEAR",
                              false,
                              phi_b_integrator->getPhysicalBcCoefs(phi_b_var));

    Pointer<Variable<NDIM>> bond_var = d_var_integrator_pairs[2].first;
    Pointer<AdvDiffHierarchyIntegrator> bond_integrator = d_var_integrator_pairs[2].second;
    int bond_scr_idx = var_db->mapVariableAndContextToIndex(bond_var, bond_integrator->getScratchContext());
    ghost_fill_comps[1] = ITC(d_bond_scr_idx,
                              bond_scr_idx,
                              "CONSERVATIVE_LINEAR_REFINE",
                              false,
                              "NONE",
                              "LINEAR",
                              false,
                              bond_integrator->getPhysicalBcCoefs(bond_var));

    HierarchyGhostCellInterpolation ghost_fill_op;
    ghost_fill_op.initializeOperatorState(ghost_fill_comps, hierarchy);
    ghost_fill_op.fillData(data_time);

    // Compute the source function
    for (int level_num = coarsest_ln; level_num <= finest_ln; ++level_num)
    {
        setDataOnPatchLevel(data_idx, var, hierarchy->getPatchLevel(level_num), data_time, initial_time);
    }

    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        Pointer<PatchLevel<NDIM>> level = hierarchy->getPatchLevel(ln);
        level->deallocatePatchData(d_phi_b_scr_idx);
        level->deallocatePatchData(d_bond_scr_idx);
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
    return;
} // setDataOnPatchHierarchy

void
BoundPlateletSource::setDataOnPatch(const int data_idx,
                                    Pointer<Variable<NDIM>> /*var*/,
                                    Pointer<Patch<NDIM>> patch,
                                    const double data_time,
                                    const bool initial_time,
                                    Pointer<PatchLevel<NDIM>> /*patch_level*/)
{
    Pointer<CellData<NDIM, double>> F_data = patch->getPatchData(data_idx);
    F_data->fillAll(0.0);
    if (initial_time) return;
    Pointer<CartesianPatchGeometry<NDIM>> pgeom = patch->getPatchGeometry();
    const double* const dx = pgeom->getDx();
    Pointer<CellData<NDIM, double>> phi_b_data = patch->getPatchData(d_phi_b_scr_idx);
    Pointer<CellData<NDIM, double>> phi_a_data =
        patch->getPatchData(d_var_integrator_pairs[1].first, d_var_integrator_pairs[1].second->getScratchContext());
    Pointer<CellData<NDIM, double>> bond_data = patch->getPatchData(d_bond_scr_idx);
    Pointer<CellData<NDIM, double>> sig_data =
        patch->getPatchData(d_var_integrator_pairs[3].first, d_var_integrator_pairs[3].second->getScratchContext());
    Pointer<CellData<NDIM, double>> w_data = patch->getPatchData(d_w_idx);
    const Box<NDIM>& patch_box = patch->getBox();
    auto psi_fcn = getKernelAndWidth(d_kernel);
    for (CellIterator<NDIM> ci(patch_box); ci; ci++)
    {
        const CellIndex<NDIM>& idx = ci();
        const double w = (*w_data)(idx);
        const double phi_a = (*phi_a_data)(idx);
        const double phi_b = (*phi_b_data)(idx);
        const double z = (*bond_data)(idx);
        // Unbound activated to bound
        const double R2 = d_Kab * phi_a *
                          convolution(d_nb_max,
                                      phi_b_data.getPointer(),
                                      -2.0,
                                      bond_data.getPointer(),
                                      psi_fcn.first,
                                      psi_fcn.second,
                                      idx,
                                      dx);

        const double R4 = d_Kaw * d_nw_max * w * phi_a;
        const double f_ab = 1.0 / d_nb * R2 + 1.0 / d_nw * R4;

        // Bound to unbound activated
        double trace = 0.0;
        for (int d = 0; d < NDIM; ++d) trace += (*sig_data)(idx, d);
        const double y_brackets = trace / (z + 1.0e-8);
        double beta = d_beta_fcn(y_brackets);

        auto lambda_fcn = [z, phi_b](const double lambda) -> std::pair<double, double> {
            double nb = z / (phi_b + 1.0e-8);
            double fcn = lambda - nb + nb * std::exp(-lambda);
            double fcn_der = 1.0 - nb * std::exp(-lambda);
            return std::make_pair(fcn, fcn_der);
        };
        double lambda = 0.0;
        if (z > 1.0e-8)
            lambda =
                boost::math::tools::newton_raphson_iterate(lambda_fcn, 1.0, 0.0, std::numeric_limits<double>::max(), 5);
        const double f_ba = beta * z * lambda / ((std::exp(lambda) - 1.0 + 1.0e-8));

        (*F_data)(idx) = d_sign * (f_ba - f_ab);
    }

    return;
} // setDataOnPatch
} // namespace clot
//////////////////////////////////////////////////////////////////////////////
