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

    /*!
     * \brief Defines the smooth Heaviside function to be used by this class
     */
    std::function<double(double)> *Heaviside; // correct this if neccesary

private:
    ActivatedPlateletSource() = delete;
    // R(c) constants
    double d_ct = std::numeric_limits<double>::quiet_NaN(); // concentration threshold for activation
    double d_r0 = std::numeric_limits<double>::quiet_NaN(); // base rate of platelet activation
    double d_k = std::numeric_limits<double>::quiet_NaN(); // platelet degredation rate?

    ActivatedPlateletSource(const ActivatedPlateletSource& from) = delete;

    ActivatedPlateletSource& operator=(const ActivatedPlateletSource& that) = delete;

    SAMRAI::tbox::Pointer<SAMRAI::hier::Variable<NDIM>> d_pl_n_var, d_c_var;
    SAMRAI::tbox::Pointer<AdvDiffHierarchyIntegrator> d_adv_diff_hier_integrator;
};
} // namespace IBAMR
//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_ActivatedPlateletSource