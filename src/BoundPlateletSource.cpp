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

BoundPlateletSource::BoundPlateletSource(std::string object_name, const BoundClotParams& clot_params)
    : CartGridFunction(std::move(object_name)), d_clot_params(clot_params)
{
    // intentionally blank
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
            auto* hier_data_ops_manager = HierarchyDataOpsManager<NDIM>::getManager();
            Pointer<HierarchyDataOpsReal<NDIM, double>> hier_cc_data_ops =
                hier_data_ops_manager->getOperationsDouble(var, hierarchy, /*get_unique*/ true);
            hier_cc_data_ops->copyData(var_scr_idx, var_cur_idx);
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
        const double R2 = d_clot_params.Kab * d_clot_params.nb_max * phi_a *
                          std::max(convolution(d_clot_params.nb_max,
                                               phi_b_data.getPointer(),
                                               -2.0,
                                               bond_data.getPointer(),
                                               psi_fcn.first,
                                               psi_fcn.second,
                                               idx,
                                               dx),
                                   0.0);
        const double R4 = d_clot_params.Kaw * d_clot_params.nw_max * w * d_clot_params.nb_max * phi_a;
        const double f_ab = R2 / d_clot_params.nb + R4 / d_clot_params.nw;

        // Bound to unbound activated
        double trace = 0.0;
        double nb_per_pl = 0.0;
        if (z > 1.0e-12 && phi_b > 3.0e-8) nb_per_pl = 1.0e4 * z / phi_b;
        for (int d = 0; d < NDIM; ++d) trace += (*sig_data)(idx, d);
        double y_brackets = 0.0;
        if (nb_per_pl > 1.0) y_brackets = std::sqrt(2.0 * trace / (z * d_clot_params.S0 + 1.0e-12));
        double beta = d_beta_fcn(y_brackets);

        double P = 1.0;
        if (nb_per_pl > 2.0)
        {
            auto lambda_fcn = [nb_per_pl](const double lambda) -> std::pair<double, double> {
                double den = 1.0 - std::exp(-lambda);
                double fcn = lambda / den - nb_per_pl;
                double fcn_der = 1.0 / den - std::exp(-lambda) * lambda / (den * den);
                return std::make_pair(fcn, fcn_der);
            };
            double lambda = 0.0;
            boost::uintmax_t max_iters = 100;
            lambda = boost::math::tools::newton_raphson_iterate(
                lambda_fcn, 1.0, 1.0, std::numeric_limits<double>::max(), 10, max_iters);
            P = lambda / (std::exp(lambda) - 1.0);
        }
        // Account for unit change.
        const double f_ba = beta * z * P * 1.0e4;

        (*F_data)(idx) = d_sign * (f_ba - f_ab);
    }

    return;
} // setDataOnPatch
} // namespace clot
//////////////////////////////////////////////////////////////////////////////
