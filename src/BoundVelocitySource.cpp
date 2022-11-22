#include <clot/BoundVelocitySource.h>
#include <clot/utility_functions.h>

#include <ADS/app_namespaces.h>

#include <FaceData.h>

namespace clot
{
BoundVelocitySource::BoundVelocitySource(const string& object_name, Pointer<Database> input_db)
    : CartGridFunction(object_name)
{
    d_c3 = input_db->getDouble("c3");
    d_vol_pl = input_db->getDouble("vol_pl");
    d_R0 = input_db->getDouble("r0");
    d_S0 = input_db->getDouble("s0");
    return;
} // UFunction

BoundVelocitySource::~BoundVelocitySource()
{
    // Prevent circular smart pointer dependency.
    d_ins_integrator = nullptr;
    d_sig_integrator = nullptr;
    d_phi_integrator = nullptr;
    d_bond_integrator = nullptr;
    d_ub_integrator = nullptr;
}

void
BoundVelocitySource::setDataOnPatchHierarchy(const int data_idx,
                                             Pointer<Variable<NDIM>> var,
                                             Pointer<PatchHierarchy<NDIM>> hierarchy,
                                             const double data_time,
                                             const bool initial_time,
                                             int coarsest_ln,
                                             int finest_ln)
{
    bool use_new_ctx = IBTK::rel_equal_eps(data_time, d_new_time);
    coarsest_ln = coarsest_ln == -1 ? 0 : coarsest_ln;
    finest_ln = finest_ln == -1 ? hierarchy->getFinestLevelNumber() : finest_ln;

    // Need to fill in scratch indices with ghost data
    using ITC = HierarchyGhostCellInterpolation::InterpolationTransactionComponent;
    std::vector<ITC> ghost_comps(4);

    HierarchyCellDataOpsReal<NDIM, double> hier_cc_data_ops(hierarchy, coarsest_ln, finest_ln);
    auto var_db = VariableDatabase<NDIM>::getDatabase();
    const int phi_scr_idx = var_db->mapVariableAndContextToIndex(d_phi_var, d_phi_integrator->getScratchContext());
    int phi_idx = var_db->mapVariableAndContextToIndex(
        d_phi_var, use_new_ctx ? d_phi_integrator->getNewContext() : d_phi_integrator->getCurrentContext());
    // If index isn't allocated, drop back to current.
    if (!d_phi_integrator->isAllocatedPatchData(phi_idx, coarsest_ln, finest_ln))
    {
        phi_idx = var_db->mapVariableAndContextToIndex(d_phi_var, d_phi_integrator->getCurrentContext());
    }
    TBOX_ASSERT(d_phi_integrator->isAllocatedPatchData(phi_idx, coarsest_ln, finest_ln));
    const bool deallocate_phi_scr = !d_phi_integrator->isAllocatedPatchData(phi_scr_idx, coarsest_ln, finest_ln);
    if (deallocate_phi_scr) d_phi_integrator->allocatePatchData(phi_scr_idx, data_time, coarsest_ln, finest_ln);

    // Now create the ITC
    ghost_comps[0] = ITC(phi_scr_idx,
                         phi_idx,
                         "CONSERVATIVE_LINEAR_REFINE",
                         false,
                         "NONE",
                         "LINEAR",
                         false,
                         d_phi_integrator->getPhysicalBcCoefs(d_phi_var));

    const int sig_scr_idx = var_db->mapVariableAndContextToIndex(d_sig_var, d_sig_integrator->getScratchContext());
    int sig_idx = var_db->mapVariableAndContextToIndex(
        d_sig_var, use_new_ctx ? d_sig_integrator->getNewContext() : d_sig_integrator->getCurrentContext());
    // If index isn't allocated, drop back to current.
    if (!d_sig_integrator->isAllocatedPatchData(sig_idx, coarsest_ln, finest_ln))
    {
        sig_idx = var_db->mapVariableAndContextToIndex(d_sig_var, d_sig_integrator->getCurrentContext());
    }
    TBOX_ASSERT(d_sig_integrator->isAllocatedPatchData(sig_idx, coarsest_ln, finest_ln));
    const bool deallocate_sig_scr = !d_sig_integrator->isAllocatedPatchData(sig_scr_idx, coarsest_ln, finest_ln);
    if (deallocate_sig_scr) d_sig_integrator->allocatePatchData(sig_scr_idx, data_time, coarsest_ln, finest_ln);

    // Now create the ITC
    ghost_comps[1] = ITC(sig_scr_idx,
                         sig_idx,
                         "CONSERVATIVE_LINEAR_REFINE",
                         false,
                         "NONE",
                         "LINEAR",
                         false,
                         d_sig_integrator->getPhysicalBcCoefs(d_sig_var));

    const int bond_scr_idx = var_db->mapVariableAndContextToIndex(d_bond_var, d_bond_integrator->getScratchContext());
    int bond_idx = var_db->mapVariableAndContextToIndex(
        d_bond_var, use_new_ctx ? d_bond_integrator->getNewContext() : d_bond_integrator->getCurrentContext());
    // If index isn't allocated, drop back to current.
    if (!d_bond_integrator->isAllocatedPatchData(bond_idx, coarsest_ln, finest_ln))
    {
        bond_idx = var_db->mapVariableAndContextToIndex(d_bond_var, d_bond_integrator->getCurrentContext());
    }
    TBOX_ASSERT(d_bond_integrator->isAllocatedPatchData(bond_idx, coarsest_ln, finest_ln));
    const bool deallocate_bond_scr = !d_bond_integrator->isAllocatedPatchData(bond_scr_idx, coarsest_ln, finest_ln);
    if (deallocate_bond_scr) d_bond_integrator->allocatePatchData(bond_scr_idx, data_time, coarsest_ln, finest_ln);

    // Now create the ITC
    ghost_comps[2] = ITC(bond_scr_idx,
                         bond_idx,
                         "CONSERVATIVE_LINEAR_REFINE",
                         false,
                         "NONE",
                         "LINEAR",
                         false,
                         d_bond_integrator->getPhysicalBcCoefs(d_bond_var));

    const int ub_scr_idx = var_db->mapVariableAndContextToIndex(d_ub_var, d_ub_integrator->getScratchContext());
    int ub_idx = var_db->mapVariableAndContextToIndex(
        d_ub_var, use_new_ctx ? d_ub_integrator->getNewContext() : d_ub_integrator->getCurrentContext());
    // If index isn't allocated, drop back to current.
    if (!d_ub_integrator->isAllocatedPatchData(ub_idx, coarsest_ln, finest_ln))
    {
        ub_idx = var_db->mapVariableAndContextToIndex(d_ub_var, d_ub_integrator->getCurrentContext());
    }
    TBOX_ASSERT(d_ub_integrator->isAllocatedPatchData(ub_idx, coarsest_ln, finest_ln));
    const bool deallocate_ub_scr = !d_ub_integrator->isAllocatedPatchData(ub_scr_idx, coarsest_ln, finest_ln);
    if (deallocate_ub_scr) d_ub_integrator->allocatePatchData(ub_scr_idx, data_time, coarsest_ln, finest_ln);

    // Now create the ITC
    ghost_comps[3] = ITC(ub_scr_idx,
                         ub_idx,
                         "CONSERVATIVE_LINEAR_REFINE",
                         false,
                         "NONE",
                         "LINEAR",
                         false,
                         d_ub_integrator->getPhysicalBcCoefs(d_ub_var));

    HierarchyGhostCellInterpolation hier_ghost_fill;
    hier_ghost_fill.initializeOperatorState(ghost_comps, hierarchy, coarsest_ln, finest_ln);
    hier_ghost_fill.fillData(data_time);

    // Now that ghost cells are filled, fill data
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        setDataOnPatchLevel(data_idx, var, hierarchy->getPatchLevel(ln), data_time, initial_time);
    }

    // Deallocate patch data if needed
    if (deallocate_phi_scr) d_phi_integrator->deallocatePatchData(phi_scr_idx, coarsest_ln, finest_ln);
    if (deallocate_sig_scr) d_sig_integrator->deallocatePatchData(sig_scr_idx, coarsest_ln, finest_ln);
    if (deallocate_bond_scr) d_bond_integrator->deallocatePatchData(bond_scr_idx, coarsest_ln, finest_ln);
    if (deallocate_ub_scr) d_ub_integrator->deallocatePatchData(ub_scr_idx, coarsest_ln, finest_ln);
}

void
BoundVelocitySource::setDataOnPatch(const int data_idx,
                                    Pointer<Variable<NDIM>> /*var*/,
                                    Pointer<Patch<NDIM>> patch,
                                    const double data_time,
                                    const bool initial_time,
                                    Pointer<PatchLevel<NDIM>> /*level*/)
{
    Pointer<CellData<NDIM, double>> ub_src_data = patch->getPatchData(data_idx);
    if (initial_time)
    {
        ub_src_data->fillAll(0.0);
        return;
    }

    Pointer<SideData<NDIM, double>> uf_data =
        patch->getPatchData(d_uf_var,
                            IBTK::rel_equal_eps(data_time, d_new_time) ? d_ins_integrator->getNewContext() :
                                                                         d_ins_integrator->getCurrentContext());
    Pointer<CellData<NDIM, double>> ub_data = patch->getPatchData(d_ub_var, d_ub_integrator->getScratchContext());
    Pointer<CellData<NDIM, double>> phi_data = patch->getPatchData(d_phi_var, d_phi_integrator->getScratchContext());
    Pointer<CellData<NDIM, double>> sig0_data = patch->getPatchData(d_sig_var, d_sig_integrator->getScratchContext());
    Pointer<CellData<NDIM, double>> bond_data = patch->getPatchData(d_bond_var, d_bond_integrator->getScratchContext());

#if !defined(NDEBUG)
    TBOX_ASSERT(ub_src_data);
    TBOX_ASSERT(uf_data);
    TBOX_ASSERT(phi_data);
    TBOX_ASSERT(sig0_data);
    TBOX_ASSERT(bond_data);
#endif

    // Convert to the correct form of stress
    Pointer<CellData<NDIM, double>> sig_data = convertToStress(*sig0_data, *bond_data, patch, d_S0, d_R0, true);

    const Box<NDIM>& box = patch->getBox();
    Pointer<CartesianPatchGeometry<NDIM>> pgeom = patch->getPatchGeometry();
    const double* const dx = pgeom->getDx();

    for (CellIterator<NDIM> ci(box); ci; ci++)
    {
        const CellIndex<NDIM>& idx = ci();

        // Get drag force
        double th = d_vol_pl * (*phi_data)(idx);
        double xi = d_c3 * th * th / (std::pow(1.0 - th, 3.0) + 1.0e-8);
        VectorNd drag_force;
        for (int d = 0; d < NDIM; ++d)
        {
            SideIndex<NDIM> idx_u(idx, d, 1);
            SideIndex<NDIM> idx_l(idx, d, 0);
            double uf = 0.5 * ((*uf_data)(idx_u) + (*uf_data)(idx_l));
            drag_force[d] = xi * (uf - (*ub_data)(idx, d));
        }

        // Divergence of stress
        IntVector<NDIM> x(1, 0), y(0, 1);
        VectorNd div_sig;
        div_sig[0] = ((*sig_data)(idx + x, 0) - (*sig_data)(idx - x, 0)) / dx[0] +
                     ((*sig_data)(idx + y, 2) - (*sig_data)(idx - y, 2)) / dx[1];
        div_sig[1] = ((*sig_data)(idx + x, 2) - (*sig_data)(idx - x, 2)) / dx[0] +
                     ((*sig_data)(idx + y, 1) - (*sig_data)(idx - y, 1)) / dx[1];

        // Total source value
        for (int d = 0; d < NDIM; ++d) (*ub_src_data)(idx, d) = drag_force[d] + div_sig[d];
    }
    return;
} // setDataOnPatch

//////////////////////////////////////////////////////////////////////////////
} // namespace clot
