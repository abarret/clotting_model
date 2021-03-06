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

#ifndef included_clot_CohesionStressUnactivatedRHS
#define included_clot_CohesionStressUnactivatedRHS

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
 * \brief Class CohesionStressUnactivatedRHS is a concrete CFRelaxationOperator that computes the relaxation function
 * for the clotting model. This class is for the model that has unactivated and activated platelets.
 */
class CohesionStressUnactivatedRHS : public IBAMR::CFRelaxationOperator
{
public:
    /*!
     * \brief This constructor reads in the parameters for the model from the input database.
     */
    CohesionStressUnactivatedRHS(SAMRAI::tbox::Pointer<SAMRAI::hier::Variable<NDIM>> phi_u_var,
                                 SAMRAI::tbox::Pointer<SAMRAI::hier::Variable<NDIM>> phi_a_var,
                                 SAMRAI::tbox::Pointer<SAMRAI::hier::Variable<NDIM>> z_var,
                                 const std::string& object_name,
                                 SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> input_db,
                                 SAMRAI::tbox::Pointer<IBAMR::AdvDiffHierarchyIntegrator> adv_diff_integrator);

    /*!
     * \brief Default constructor.
     *
     * \note This constructor is not implemented and should not be used.
     */
    CohesionStressUnactivatedRHS() = delete;

    /*!
     * \brief Copy constructor.
     *
     * \note This constructor is not implemented and should not be used.
     */
    CohesionStressUnactivatedRHS(const CohesionStressUnactivatedRHS& from) = delete;

    /*!
     * \brief Assignment operator.
     *
     * \note This operator is not implemented and should not be used.
     */
    CohesionStressUnactivatedRHS& operator=(const CohesionStressUnactivatedRHS& that) = delete;

    /*!
     * \brief Empty destructor.
     */
    ~CohesionStressUnactivatedRHS() = default;

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

private:
    int d_w_idx = IBTK::invalid_index;
    double d_c4 = std::numeric_limits<double>::quiet_NaN();
    double d_a0 = std::numeric_limits<double>::quiet_NaN();
    double d_a0w = std::numeric_limits<double>::quiet_NaN();
    SAMRAI::tbox::Pointer<SAMRAI::hier::Variable<NDIM>> d_phi_u_var, d_phi_a_var, d_z_var;

    // When X is greater than d_clot_break_x, set beta equal to d_beta_limit;
    double d_clot_break_x = std::numeric_limits<double>::quiet_NaN();
    double d_beta_limit = std::numeric_limits<double>::quiet_NaN();

    // Beta function pointer
    std::function<double(double)> d_beta_fcn;

    SAMRAI::tbox::Pointer<IBAMR::AdvDiffHierarchyIntegrator> d_adv_diff_integrator;
};

} // Namespace clot
#endif
