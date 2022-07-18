// ---------------------------------------------------------------------
//
// Copyright (c) 2014 - 2014 by the IBAMR developers
// All rights reserved.
//
// This file is part of IBAMR.
//
// IBAMR is free software and is distributed under the 3-clause BSD
// license. The full text of the license can be found in the file
// COPYRIGHT at the top level directory of IBAMR.
//
// ---------------------------------------------------------------------

#ifndef included_clot_BondBoundSource
#define included_clot_BondBoundSource

#include <clot/utility_functions.h>

#include <ibamr/AdvDiffHierarchyIntegrator.h>

#include <tbox/Database.h>
#include <tbox/Pointer.h>
#include <tbox/Utilities.h>

#include <Box.h>
#include <CartesianPatchGeometry.h>
#include <CellData.h>
#include <CellIndex.h>
#include <CellVariable.h>
#include <Index.h>
#include <IntVector.h>
#include <Patch.h>
#include <PatchGeometry.h>
#include <PatchHierarchy.h>
#include <PatchLevel.h>
#include <Variable.h>

#include <functional>

/////////////////////////////// CLASS DEFINITION /////////////////////////////
namespace clot
{
/*!
 * Class BondBoundSource is a concrete CartGridFunction that implements the growth and decay terms for the bond
 * variable.
 *
 * Note this function should only be used for the model that has unactivated and activated platelets.
 */
class BondBoundSource : public IBTK::CartGridFunction
{
public:
    /*!
     * \brief Class constructor.
     */
    BondBoundSource(std::string object_name, SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> input_db);

    /*!
     * \brief Empty destructor.
     */
    ~BondBoundSource() = default;

    /*!
     * \brief Deleted constructors
     */
    //\{
    BondBoundSource(const BondBoundSource& from) = delete;
    BondBoundSource& operator=(const BondBoundSource& that) = delete;
    //\}

    /*!
     * \brief Indicates whether the concrete BondBoundSource object is
     * time-dependent.
     */
    bool isTimeDependent() const override;

    /*!
     * \brief Evaluate the function on the patch interiors on the specified
     * levels of the patch hierarchy.
     */
    void setDataOnPatchHierarchy(const int data_idx,
                                 SAMRAI::tbox::Pointer<SAMRAI::hier::Variable<NDIM>> var,
                                 SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM>> hierarchy,
                                 const double data_time,
                                 const bool initial_time = false,
                                 const int coarsest_ln = -1,
                                 const int finest_ln = -1) override;

    /*!
     * \brief Evaluate the function on the patch interior.
     */
    void setDataOnPatch(const int data_idx,
                        SAMRAI::tbox::Pointer<SAMRAI::hier::Variable<NDIM>> var,
                        SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM>> patch,
                        const double data_time,
                        const bool initial_time = false,
                        SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM>> patch_level =
                            SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM>>(NULL)) override;

    /*!
     * \brief Methods to set the inputs
     */
    //\{
    void setBoundPlateletData(SAMRAI::tbox::Pointer<SAMRAI::hier::Variable<NDIM>> var,
                              SAMRAI::tbox::Pointer<IBTK::HierarchyIntegrator> integrator)
    {
        d_var_integrator_pairs[0].first = var;
        d_var_integrator_pairs[0].second = integrator;
        // Create extra patch index with large ghost width
        auto var_db = SAMRAI::hier::VariableDatabase<NDIM>::getDatabase();
        d_phi_b_scr_idx =
            var_db->registerVariableAndContext(var, var_db->getContext(d_object_name + "::ScrCtx"), 4 /*ghosts*/);
    }

    void setActivatedPlateletData(SAMRAI::tbox::Pointer<SAMRAI::hier::Variable<NDIM>> var,
                                  SAMRAI::tbox::Pointer<IBTK::HierarchyIntegrator> integrator)
    {
        d_var_integrator_pairs[1].first = var;
        d_var_integrator_pairs[1].second = integrator;
    }

    void setBondData(SAMRAI::tbox::Pointer<SAMRAI::hier::Variable<NDIM>> var,
                     SAMRAI::tbox::Pointer<IBTK::HierarchyIntegrator> integrator)
    {
        d_var_integrator_pairs[2].first = var;
        d_var_integrator_pairs[2].second = integrator;
        // Create extra patch index with large ghost width
        auto var_db = SAMRAI::hier::VariableDatabase<NDIM>::getDatabase();
        d_bond_scr_idx =
            var_db->registerVariableAndContext(var, var_db->getContext(d_object_name + "::ScrCtx"), 4 /*ghosts*/);
    }

    void setStressData(SAMRAI::tbox::Pointer<SAMRAI::hier::Variable<NDIM>> var,
                       SAMRAI::tbox::Pointer<IBTK::HierarchyIntegrator> integrator)
    {
        d_var_integrator_pairs[3].first = var;
        d_var_integrator_pairs[3].second = integrator;
    }

    inline void registerBetaFcn(std::function<double(double, void*)> wrapper, void* beta)
    {
        d_beta_fcn = std::bind(wrapper, std::placeholders::_1, beta);
    }

    inline void setWIdx(const int w_idx)
    {
        d_w_idx = w_idx;
    }

    inline void setKernel(const clot::Kernel kern)
    {
        d_kernel = kern;
    }
    //\}

private:
    double d_Kab = std::numeric_limits<double>::quiet_NaN(), d_Kbb = std::numeric_limits<double>::quiet_NaN(),
           d_Kaw = std::numeric_limits<double>::quiet_NaN(), d_nb_max = std::numeric_limits<double>::quiet_NaN(),
           d_nw_max = std::numeric_limits<double>::quiet_NaN();

    Kernel d_kernel = UNKNOWN_KERNEL;

    std::array<std::pair<SAMRAI::tbox::Pointer<SAMRAI::hier::Variable<NDIM>>,
                         SAMRAI::tbox::Pointer<IBTK::HierarchyIntegrator>>,
               4>
        d_var_integrator_pairs;
    int d_w_idx = IBTK::invalid_index;

    // Scratch indices with extra ghost cells
    int d_phi_b_scr_idx = IBTK::invalid_index, d_bond_scr_idx = IBTK::invalid_index;

    // Beta function pointer
    std::function<double(double)> d_beta_fcn;

    // wall index
    int d_w_idx = IBTK::invalid_index;
};
} // namespace clot

#endif
