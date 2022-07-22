#include <clot/BoundVelocityFunction.h>
#include <clot/utility_functions.h>

#include <ADS/app_namespaces.h>

#include <FaceData.h>

namespace clot
{
BoundVelocityFunction::BoundVelocityFunction(const string& object_name, Pointer<Database> input_db)
    : CartGridFunction(object_name)
{
    d_c3 = input_db->getDouble("c3");
    d_thresh = input_db->getDoubleWithDefault("threshold", d_thresh);
    d_vol_pl = input_db->getDouble("vol_pl");
    d_R0 = input_db->getDouble("r0");
    d_S0 = input_db->getDouble("s0");
    return;
} // UFunction

BoundVelocityFunction::~BoundVelocityFunction()
{
    // Prevent circular smart pointer dependency.
    d_ins_integrator = nullptr;
    d_sig_integrator = nullptr;
    d_phi_integrator = nullptr;
    d_bond_integrator = nullptr;
}

void
BoundVelocityFunction::setDataOnPatchHierarchy(const int data_idx,
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
    std::vector<ITC> ghost_comps(3);

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
}

void
BoundVelocityFunction::setDataOnPatch(const int data_idx,
                                      Pointer<Variable<NDIM>> /*var*/,
                                      Pointer<Patch<NDIM>> patch,
                                      const double data_time,
                                      const bool initial_time,
                                      Pointer<PatchLevel<NDIM>> /*level*/)
{
    Pointer<FaceData<NDIM, double>> ub_data = patch->getPatchData(data_idx);
    if (initial_time)
    {
        ub_data->fillAll(0.0);
        return;
    }

    Pointer<SideData<NDIM, double>> uf_data =
        patch->getPatchData(d_uf_var,
                            IBTK::rel_equal_eps(data_time, d_new_time) ? d_ins_integrator->getNewContext() :
                                                                         d_ins_integrator->getCurrentContext());
    Pointer<CellData<NDIM, double>> phi_data = patch->getPatchData(d_phi_var, d_phi_integrator->getScratchContext());
    Pointer<CellData<NDIM, double>> sig0_data = patch->getPatchData(d_sig_var, d_sig_integrator->getScratchContext());
    Pointer<CellData<NDIM, double>> bond_data = patch->getPatchData(d_bond_var, d_bond_integrator->getScratchContext());

#if !defined(NDEBUG)
    TBOX_ASSERT(ub_data);
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
    for (int axis = 0; axis < NDIM; ++axis)
    {
        for (FaceIterator<NDIM> fi(box, axis); fi; fi++)
        {
            const FaceIndex<NDIM>& idx = fi();
            const SideIndex<NDIM> s_idx(idx.toCell(0), axis, 1);

            // Get the platelet count
            double phi = 0.5 * ((*phi_data)(idx.toCell(0)) + (*phi_data)(idx.toCell(1)));
            // Volume fraction
            double th = d_vol_pl * phi;
            // Drag coefficient
            double xi = d_c3 * th * th / (std::pow(1.0 - th, 3.0) + 1.0e-8);

            if (xi < d_thresh)
            {
                (*ub_data)(idx) = (*uf_data)(s_idx);
            }
            else
            {
                // Now compute the divergence of the stress
                double div_sig = 0.0;
                if (axis == 0)
                {
                    div_sig = ((*sig_data)(idx.toCell(1), 0) - (*sig_data)(idx.toCell(0), 0)) / dx[0];
                    IntVector<NDIM> up(0);
                    up(0) = 1; // Note index is switched because we're looping over FACE indices.
                    double sig_up = 0.5 * ((*sig_data)((idx + up).toCell(1), 2) + (*sig_data)((idx + up).toCell(0), 2));
                    double sig_low =
                        0.5 * ((*sig_data)((idx - up).toCell(1), 2) + (*sig_data)((idx - up).toCell(0), 2));
                    div_sig += (sig_up - sig_low) / (2.0 * dx[1]);
                }
                else if (axis == 1)
                {
                    div_sig = ((*sig_data)(idx.toCell(1), 1) - (*sig_data)(idx.toCell(0), 1)) / dx[1];
                    IntVector<NDIM> up(0);
                    up(1) = 1; // Note index is switched because we're looping over FACE indices.
                    double sig_up = 0.5 * ((*sig_data)((idx + up).toCell(1), 2) + (*sig_data)((idx + up).toCell(0), 2));
                    double sig_low =
                        0.5 * ((*sig_data)((idx - up).toCell(1), 2) + (*sig_data)((idx - up).toCell(0), 2));
                    div_sig += (sig_up - sig_low) / (2.0 * dx[0]);
                }
                else if (axis == 2)
                {
                    TBOX_ERROR("3D NOT SUPPORTED YET!\n\n");
                }
                else
                {
                    TBOX_ERROR("SHOULD NOT GET HERE!\n\n");
                }

                // Now compute the bound velocity field. Note we've handled the case of xi small above.
                (*ub_data)(idx) = (*uf_data)(s_idx) + div_sig / (xi + 1.0e-8);
            }
        }
    }

    return;
} // setDataOnPatch

//////////////////////////////////////////////////////////////////////////////
}
