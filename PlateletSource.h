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

#ifndef included_PlateletSource
#define included_PlateletSource

#include <ibamr/AdvDiffHierarchyIntegrator.h>

#include "utility_functions.h"

/////////////////////////////// CLASS DEFINITION /////////////////////////////
namespace IBAMR
{
/*!
 * \brief Class PlateletSource provides a source term for the activated platelet concentration.
 */
class PlateletSource : public IBTK::CartGridFunction
{
public:
    /*!
     * \brief Class constructor.
     */
    PlateletSource(SAMRAI::tbox::Pointer<SAMRAI::hier::Variable<NDIM>> phi_u_var,
                   SAMRAI::tbox::Pointer<SAMRAI::hier::Variable<NDIM>> phi_a_var,
                   SAMRAI::tbox::Pointer<SAMRAI::hier::Variable<NDIM>> w_var,
                   SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> input_db,
                   SAMRAI::tbox::Pointer<AdvDiffHierarchyIntegrator> adv_diff_hier_integrator);

    /*!
     * \brief Empty destructor.
     */
    ~PlateletSource() = default;

    /*!
     * \name Methods to set patch data.
     */
    //\{

    /*!
     * \brief Indicates whether the concrete PlateletSource object is
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

    //\}

    inline void setKernel(const IBAMR::Kernel kern)
    {
        d_kernel = kern;
    }

    inline void setSign(const bool positive)
    {
        d_sign = (positive) ? 1.0 : -1.0;
    }

private:
    PlateletSource() = delete;
    double d_Kua = std::numeric_limits<double>::quiet_NaN();
    double d_Kuw = std::numeric_limits<double>::quiet_NaN();
    double d_w_mx = std::numeric_limits<double>::quiet_NaN();
    double d_sign = 1.0;
    Kernel d_kernel = UNKNOWN_KERNEL;

    PlateletSource(const PlateletSource& from) = delete;

    PlateletSource& operator=(const PlateletSource& that) = delete;

    SAMRAI::tbox::Pointer<SAMRAI::hier::Variable<NDIM>> d_phi_u_var, d_phi_a_var, d_w_var;
    SAMRAI::tbox::Pointer<AdvDiffHierarchyIntegrator> d_adv_diff_hier_integrator;
};
} // namespace IBAMR
//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_PlateletSource
