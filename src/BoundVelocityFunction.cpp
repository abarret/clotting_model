#include <clot/BoundVelocityFunction.h>

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
    return;
} // UFunction

BoundVelocityFunction::~BoundVelocityFunction()
{
    // Prevent circular smart pointer dependency.
    d_ins_integrator = nullptr;
    d_sig_integrator = nullptr;
    d_phi_integrator = nullptr;
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

    const int sig_scr_idx = var_db->mapVariableAndContextToIndex(d_sig_var, d_sig_integrator->getScratchContext());
    const int sig_cur_idx = var_db->mapVariableAndContextToIndex(d_sig_var, d_sig_integrator->getCurrentContext());
    const int sig_new_idx = var_db->mapVariableAndContextToIndex(d_sig_var, d_sig_integrator->getNewContext());
    const bool sig_new_allocated = d_sig_integrator->isAllocatedPatchData(sig_new_idx, coarsest_ln, finest_ln);
    const bool deallocate_sig_scr = !d_sig_integrator->isAllocatedPatchData(sig_scr_idx, coarsest_ln, finest_ln);
    if (deallocate_sig_scr) d_sig_integrator->allocatePatchData(sig_scr_idx, data_time, coarsest_ln, finest_ln);

    if (sig_new_allocated)
    {
        // Use the average of phi_cur and phi_new
        hier_cc_data_ops.linearSum(sig_scr_idx, 0.5, sig_cur_idx, 0.5, sig_new_idx);
    }
    else
    {
        hier_cc_data_ops.copyData(sig_scr_idx, sig_cur_idx);
    }

    // Now create the ITC
    ghost_comps[1] = ITC(sig_scr_idx,
                         "CONSERVATIVE_LINEAR_REFINE",
                         false,
                         "NONE",
                         "LINEAR",
                         false,
                         d_sig_integrator->getPhysicalBcCoefs(d_sig_var));

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
    if (deallocate_sig_scr) d_sig_integrator->deallocatePatchData(sig_scr_idx, coarsest_ln, finest_ln);
    if (deallocate_uf_scr) d_ins_integrator->deallocatePatchData(uf_scr_idx, coarsest_ln, finest_ln);
}

void
BoundVelocityFunction::setDataOnPatch(const int data_idx,
                                      Pointer<Variable<NDIM>> /*var*/,
                                      Pointer<Patch<NDIM>> patch,
                                      const double /*data_time*/,
                                      const bool initial_time,
                                      Pointer<PatchLevel<NDIM>> /*level*/)
{
    Pointer<FaceData<NDIM, double>> ub_data = patch->getPatchData(data_idx);
    Pointer<FaceData<NDIM, double>> uf_data = patch->getPatchData(d_uf_var, d_ins_integrator->getScratchContext());
    Pointer<CellData<NDIM, double>> phi_data = patch->getPatchData(d_phi_var, d_phi_integrator->getScratchContext());
    Pointer<CellData<NDIM, double>> sig_data = patch->getPatchData(d_sig_var, d_sig_integrator->getScratchContext());

#if !defined(NDEBUG)
    TBOX_ASSERT(ub_data);
    TBOX_ASSERT(uf_data);
    TBOX_ASSERT(phi_data);
    TBOX_ASSERT(sig_data);
#endif

    if (initial_time)
    {
        ub_data->fillAll(0.0);
        return;
    }

    const Box<NDIM>& box = patch->getBox();
    Pointer<CartesianPatchGeometry<NDIM>> pgeom = patch->getPatchGeometry();
    const double* const dx = pgeom->getDx();
    for (int axis = 0; axis < NDIM; ++axis)
    {
        for (FaceIterator<NDIM> fi(box, axis); fi; fi++)
        {
            const FaceIndex<NDIM>& idx = fi();

            // Get the platelet count
            double phi = 0.5 * ((*phi_data)(idx.toCell(0)) + (*phi_data)(idx.toCell(1)));
            // Volume fraction
            double th = d_vol_pl * phi;
            // Drag coefficient
            double xi = d_c3 * th * th / std::pow(1.0 - th, 3.0);

            if (xi < d_thresh)
            {
                (*ub_data)(idx) = (*uf_data)(idx);
            }
            else
            {
                // Now compute the divergence of the stress
                double grad_sig = 0.0;
                if (axis == 0)
                {
                    grad_sig = ((*sig_data)(idx.toCell(1)) - (*sig_data)(idx.toCell(0))) / dx[0];
                    IntVector<NDIM> up(0);
                    up(1) = 1;
                    double sig_up = 0.5 * ((*sig_data)((idx + up).toCell(1)) + (*sig_data)((idx + up).toCell(0)));
                    double sig_low = 0.5 * ((*sig_data)((idx - up).toCell(1)) + (*sig_data)((idx - up).toCell(0)));
                    grad_sig += (sig_up - sig_low) / (2.0 * dx[1]);
                }
                else if (axis == 1)
                {
                    grad_sig = ((*sig_data)(idx.toCell(1)) - (*sig_data)(idx.toCell(0))) / dx[1];
                    IntVector<NDIM> up(0);
                    up(0) = 1;
                    double sig_up = 0.5 * ((*sig_data)((idx + up).toCell(1)) + (*sig_data)((idx + up).toCell(0)));
                    double sig_low = 0.5 * ((*sig_data)((idx - up).toCell(1)) + (*sig_data)((idx - up).toCell(0)));
                    grad_sig += (sig_up - sig_low) / (2.0 * dx[0]);
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
                (*ub_data)(idx) = (*uf_data)(idx) + grad_sig / xi;
            }
        }
    }

    return;
} // setDataOnPatch

//////////////////////////////////////////////////////////////////////////////
}
