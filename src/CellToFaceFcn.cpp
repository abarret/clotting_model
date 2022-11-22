#include <clot/CellToFaceFcn.h>
#include <clot/utility_functions.h>

#include <ADS/app_namespaces.h>

#include <CellDataFactory.h>
#include <FaceData.h>
#include <HierarchyDataOpsManager.h>

namespace clot
{
CellToFaceFcn::CellToFaceFcn(std::string object_name,
                             Pointer<CellVariable<NDIM, double>> Q_var,
                             Pointer<AdvDiffHierarchyIntegrator> integrator)
    : CartGridFunction(std::move(object_name)), d_Q_var(Q_var), d_integrator(integrator)
{
    Pointer<CellDataFactory<NDIM, double>> fac = Q_var->getPatchDataFactory();
    TBOX_ASSERT(fac->getDefaultDepth() == NDIM);
    // intentionally blank
}

void
CellToFaceFcn::setDataOnPatchHierarchy(const int data_idx,
                                       Pointer<Variable<NDIM>> var,
                                       Pointer<PatchHierarchy<NDIM>> hierarchy,
                                       const double data_time,
                                       const bool initial_time,
                                       int coarsest_ln,
                                       int finest_ln)
{
    coarsest_ln = coarsest_ln == IBTK::invalid_level_number ? 0 : coarsest_ln;
    finest_ln = finest_ln == IBTK::invalid_level_number ? hierarchy->getFinestLevelNumber() : finest_ln;

    auto var_db = VariableDatabase<NDIM>::getDatabase();
    const int Q_cur_idx = var_db->mapVariableAndContextToIndex(d_Q_var, d_integrator->getCurrentContext());
    const int Q_scr_idx = var_db->mapVariableAndContextToIndex(d_Q_var, d_integrator->getScratchContext());
    bool scr_allocated = d_integrator->isAllocatedPatchData(Q_scr_idx, coarsest_ln, finest_ln);

    if (!scr_allocated) d_integrator->allocatePatchData(Q_scr_idx, data_time, coarsest_ln, finest_ln);
    if (!initial_time)
    {
        auto* hier_data_ops_manager = HierarchyDataOpsManager<NDIM>::getManager();
        Pointer<HierarchyDataOpsReal<NDIM, double>> hier_cc_data_ops =
            hier_data_ops_manager->getOperationsDouble(d_Q_var, hierarchy, /*get_unique*/ true);
        hier_cc_data_ops->copyData(Q_scr_idx, Q_cur_idx);
    }

    // We need boundary conditions
    using ITC = HierarchyGhostCellInterpolation::InterpolationTransactionComponent;
    std::vector<ITC> ghost_cell_comps(1);
    ghost_cell_comps[0] = ITC(Q_scr_idx,
                              "CONSERVATIVE_LINEAR_REFINE",
                              false,
                              "NONE",
                              "LINEAR",
                              false,
                              d_integrator->getPhysicalBcCoefs(d_Q_var));
    HierarchyGhostCellInterpolation hier_ghost_cell;
    hier_ghost_cell.initializeOperatorState(ghost_cell_comps, hierarchy, coarsest_ln, finest_ln);
    hier_ghost_cell.fillData(data_time);

    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
        setDataOnPatchLevel(data_idx, var, hierarchy->getPatchLevel(ln), data_time, initial_time);

    // Deallocate when necessary
    if (!scr_allocated) d_integrator->deallocatePatchData(Q_scr_idx, coarsest_ln, finest_ln);
}

void
CellToFaceFcn::setDataOnPatch(const int data_idx,
                              Pointer<Variable<NDIM>> /*var*/,
                              Pointer<Patch<NDIM>> patch,
                              const double data_time,
                              const bool initial_time,
                              Pointer<PatchLevel<NDIM>> /*level*/)
{
    Pointer<FaceData<NDIM, double>> Q_fc_data = patch->getPatchData(data_idx);
    Pointer<CellData<NDIM, double>> Q_cc_data = patch->getPatchData(d_Q_var, d_integrator->getScratchContext());
    if (initial_time)
    {
        Q_fc_data->fillAll(0.0);
        return;
    }
    PatchMathOps patch_math_ops;
    patch_math_ops.interp(Q_fc_data, Q_cc_data, patch);
    return;
} // setDataOnPatch

//////////////////////////////////////////////////////////////////////////////
} // namespace clot
