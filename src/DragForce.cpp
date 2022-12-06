#include <clot/DragForce.h>
#include <clot/utility_functions.h>

#include <ADS/app_namespaces.h>

#include <FaceData.h>

namespace clot
{
DragForce::DragForce(const string& object_name, BoundClotParams clot_params)
    : CartGridFunction(object_name), d_clot_params(std::move(clot_params))
{
    // intentionall blank
    return;
}

DragForce::~DragForce()
{
    // Prevent circular smart pointer dependency.
    d_ins_integrator = nullptr;
    d_phi_integrator = nullptr;
}

void
DragForce::setDataOnPatchHierarchy(const int data_idx,
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
    std::vector<ITC> ghost_comps(2);

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

    const int ub_scr_idx = var_db->mapVariableAndContextToIndex(d_ub_var, d_ub_integrator->getScratchContext());
    int ub_idx = var_db->mapVariableAndContextToIndex(
        d_ub_var, use_new_ctx ? d_phi_integrator->getNewContext() : d_ub_integrator->getCurrentContext());
    // If index isn't allocated, drop back to current.
    if (!d_ub_integrator->isAllocatedPatchData(ub_idx, coarsest_ln, finest_ln))
    {
        ub_idx = var_db->mapVariableAndContextToIndex(d_ub_var, d_ub_integrator->getCurrentContext());
    }
    TBOX_ASSERT(d_phi_integrator->isAllocatedPatchData(ub_idx, coarsest_ln, finest_ln));
    const bool deallocate_ub_scr = !d_phi_integrator->isAllocatedPatchData(ub_scr_idx, coarsest_ln, finest_ln);
    if (deallocate_ub_scr) d_ub_integrator->allocatePatchData(ub_scr_idx, data_time, coarsest_ln, finest_ln);

    // Now create the ITC
    ghost_comps[1] = ITC(ub_scr_idx,
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
    if (deallocate_ub_scr) d_ub_integrator->deallocatePatchData(ub_scr_idx, coarsest_ln, finest_ln);
}

void
DragForce::setDataOnPatch(const int data_idx,
                          Pointer<Variable<NDIM>> /*var*/,
                          Pointer<Patch<NDIM>> patch,
                          const double data_time,
                          const bool initial_time,
                          Pointer<PatchLevel<NDIM>> /*level*/)
{
    Pointer<SideData<NDIM, double>> drag_data = patch->getPatchData(data_idx);
    if (initial_time)
    {
        drag_data->fillAll(0.0);
        return;
    }

    Pointer<SideData<NDIM, double>> uf_data =
        patch->getPatchData(d_uf_var,
                            IBTK::rel_equal_eps(data_time, d_new_time) ? d_ins_integrator->getNewContext() :
                                                                         d_ins_integrator->getCurrentContext());
    Pointer<CellData<NDIM, double>> phi_data = patch->getPatchData(d_phi_var, d_phi_integrator->getScratchContext());
    Pointer<CellData<NDIM, double>> ub_data = patch->getPatchData(d_ub_var, d_ub_integrator->getScratchContext());

#if !defined(NDEBUG)
    TBOX_ASSERT(drag_data);
    TBOX_ASSERT(ub_data);
    TBOX_ASSERT(uf_data);
    TBOX_ASSERT(phi_data);
#endif

    const Box<NDIM>& box = patch->getBox();
    Pointer<CartesianPatchGeometry<NDIM>> pgeom = patch->getPatchGeometry();

    for (int axis = 0; axis < NDIM; ++axis)
    {
        for (SideIterator<NDIM> si(box, axis); si; si++)
        {
            const SideIndex<NDIM>& idx = si();

            // Get the drag coefficient
            double phi = 0.5 * ((*phi_data)(idx.toCell(1)) + (*phi_data)(idx.toCell(0)));
            double th = d_clot_params.vol_pl * phi;
            double xi = d_clot_params.drag_coef * th * th / (std::pow(1.0 - th, 3.0) + 1.0e-12);

            // Now compute drag force
            double ub = 0.5 * ((*ub_data)(idx.toCell(1), axis) + (*ub_data)(idx.toCell(0), axis));
            (*drag_data)(idx) = xi * (ub - (*uf_data)(idx));
        }
    }
    return;
} // setDataOnPatch

//////////////////////////////////////////////////////////////////////////////
} // namespace clot
