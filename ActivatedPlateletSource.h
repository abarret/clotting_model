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

#ifndef included_ActivatedPlateletSource
#define included_ActivatedPlateletSource

#include <ibamr/AdvDiffHierarchyIntegrator.h>
#include "utility_functions.h"

/////////////////////////////// CLASS DEFINITION /////////////////////////////
namespace IBAMR
{
/*!
 * \brief Class ActivatedPlateletSource provides a source term for the activated platelet concentration.
 */
class ActivatedPlateletSource : public IBTK::CartGridFunction
{
public:
    /*!
     * \brief Class constructor.
     */
    ActivatedPlateletSource(SAMRAI::tbox::Pointer<SAMRAI::hier::Variable<NDIM>> pl_n_var,
                            SAMRAI::tbox::Pointer<SAMRAI::hier::Variable<NDIM>> c_var,
                            SAMRAI::tbox::Pointer<AdvDiffHierarchyIntegrator> adv_diff_hier_integrator);

    /*!
     * \brief Empty destructor.
     */
    ~ActivatedPlateletSource() = default;

    /*!
     * \name Methods to set patch data.
     */
    //\{

    /*!
     * \brief Indicates whether the concrete ActivatedPlateletSource object is
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

    /*!
     * \brief Evaluate the platlet sources
     */
    void setDataOnPatchLevel(const int data_idx,
                            SAMRAI::tbox::Pointer<SAMRAI::hier::Variable<NDIM>> var,
                            SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM>> p_level,
                            const double data_time,
                            const bool initial_time = false) override;

    inline void setKernel(const IBAMR::Kernel kern)
    {
        d_kernel = kern;
    }


private:
    ActivatedPlateletSource() = delete;
    double d_a0 = std::numeric_limits<double>::quiet_NaN();
    double d_a0w = std::numeric_limits<double>::quiet_NaN();
    double d_w_mx = std::numeric_limits<double>::quiet_NaN();
    IBAMR::Kernel d_kernel = UNKNOWN_KERNEL

    ActivatedPlateletSource(const ActivatedPlateletSource& from) = delete;

    ActivatedPlateletSource& operator=(const ActivatedPlateletSource& that) = delete;

    SAMRAI::tbox::Pointer<SAMRAI::hier::Variable<NDIM>> d_phi_a_var, d_phi_u_var, d_w_var;
    SAMRAI::tbox::Pointer<AdvDiffHierarchyIntegrator> d_adv_diff_hier_integrator;
};
} // namespace IBAMR
//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_ActivatedPlateletSource
