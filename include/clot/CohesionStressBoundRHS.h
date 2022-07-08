// ---------------------------------------------------------------------
//
// Copyright (c) 2019 - 2020 by the IBAMR developers
// All rights reserved.
//
// This file is part of IBAMR.
//
// IBAMR is free software and is distributed under the 3-clause BSD
// license. The full text of the license can be found in the file
// COPYRIGHT at the top level directory of IBAMR.
//
// ---------------------------------------------------------------------

/////////////////////////////// INCLUDE GUARD ////////////////////////////////

#ifndef included_clot_CohesionStressBoundRHS
#define included_clot_CohesionStressBoundRHS

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <ibamr/config.h>

#include <clot/utility_functions.h>

#include <ibamr/CFRelaxationOperator.h>

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
#include <string>

/////////////////////////////// CLASS DEFINITION /////////////////////////////

namespace clot
{
/*!
 * \brief Class CohesionStressBoundRHS is a concrete CFRelaxationOperator that computes the relaxation function for the
 * clotting model. This model assumes activated and bound platelets.
 */
class CohesionStressBoundRHS : public IBAMR::CFRelaxationOperator
{
public:
    /*!
     * \brief This constructor reads in the parameters for the model from the input database.
     */
    CohesionStressBoundRHS(std::string object_name, SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> input_db);

    /*!
     * \brief Deleted constructors
     */
    //\{
    CohesionStressBoundRHS() = delete;
    CohesionStressBoundRHS(const CohesionStressBoundRHS& from) = delete;
    CohesionStressBoundRHS& operator=(const CohesionStressBoundRHS& that) = delete;
    //\}

    /*!
     * \brief Empty destructor.
     */
    ~CohesionStressBoundRHS();

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
    //\}

    /*!
     * \name Methods to set patch data.
     */
    //\{

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
    //\}

    inline void registerBetaFcn(std::function<double(double, void*)> wrapper, void* beta)
    {
        // We set beta here
        d_beta_fcn = std::bind(wrapper, std::placeholders::_1, beta);
    }

    inline void setOmegaIdx(const int w_idx)
    {
        d_w_idx = w_idx;
    }

    inline void setKernel(const clot::Kernel kern)
    {
        d_kernel = kern;
    }

private:
    double d_c4 = std::numeric_limits<double>::quiet_NaN(), d_Kab = std::numeric_limits<double>::quiet_NaN(),
           d_Kbb = std::numeric_limits<double>::quiet_NaN(), d_Kaw = std::numeric_limits<double>::quiet_NaN(),
           d_nb_max = std::numeric_limits<double>::quiet_NaN(), d_nw_max = std::numeric_limits<double>::quiet_NaN();

    Kernel d_kernel = UNKNOWN_KERNEL;

    std::array<std::pair<SAMRAI::tbox::Pointer<SAMRAI::hier::Variable<NDIM>>,
                         SAMRAI::tbox::Pointer<IBTK::HierarchyIntegrator>>,
               3>
        d_var_integrator_pairs;
    int d_w_idx = IBTK::invalid_index;

    // Scratch indices with extra ghost cells
    int d_phi_b_scr_idx = IBTK::invalid_index, d_bond_scr_idx = IBTK::invalid_index;

    // Beta function pointer
    std::function<double(double)> d_beta_fcn;

    SAMRAI::tbox::Pointer<IBAMR::AdvDiffHierarchyIntegrator> d_adv_diff_integrator;
};

} // Namespace clot
#endif
