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

#ifndef included_clot_BoundPlateletSource
#define included_clot_BoundPlateletSource

#include <clot/utility_functions.h>

#include <ibamr/AdvDiffHierarchyIntegrator.h>

/////////////////////////////// CLASS DEFINITION /////////////////////////////
namespace clot
{
/*!
 * \brief Class BoundPlateletSource provides a source term for the activated and bound platelets.
 *
 * \note This is used in the iteration of the model that has activated and bound platelets only.
 */
class BoundPlateletSource : public IBTK::CartGridFunction
{
public:
    /*!
     * \brief Class constructor.
     */
    BoundPlateletSource(std::string object_name,
                        SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> input_db);
    /*!
     * \brief Deleted constructors.
     */
    //\{
    BoundPlateletSource() = delete;
    BoundPlateletSource(const BoundPlateletSource& from) = delete;
    BoundPlateletSource& operator=(const BoundPlateletSource& that) = delete;
    //\}

    /*!
     * \brief Empty destructor.
     */
    ~BoundPlateletSource() = default;

    inline void setBoundPlateletData(SAMRAI::tbox::Pointer<SAMRAI::hier::Variable<NDIM>> phi_b_var,
                                     SAMRAI::tbox::Pointer<IBAMR::AdvDiffHierarchyIntegrator> phi_b_integrator)
    {
        d_var_integrator_pairs[0].first = phi_b_var;
        d_var_integrator_pairs[0].second = phi_b_integrator;
        // Create scratch index with extra ghost cells
        auto var_db = SAMRAI::hier::VariableDatabase<NDIM>::getDatabase();
        d_phi_b_scr_idx = var_db->registerVariableAndContext(phi_b_var, var_db->getContext(d_object_name + "::CTX"), 4);
    }
    inline void setActivatedPlateletData(SAMRAI::tbox::Pointer<SAMRAI::hier::Variable<NDIM>> phi_a_var,
                                         SAMRAI::tbox::Pointer<IBAMR::AdvDiffHierarchyIntegrator> phi_a_integrator)
    {
        d_var_integrator_pairs[1].first = phi_a_var;
        d_var_integrator_pairs[1].second = phi_a_integrator;
    }
    inline void setBondsData(SAMRAI::tbox::Pointer<SAMRAI::hier::Variable<NDIM>> bond_var,
                             SAMRAI::tbox::Pointer<IBAMR::AdvDiffHierarchyIntegrator> bond_integrator)
    {
        d_var_integrator_pairs[2].first = bond_var;
        d_var_integrator_pairs[2].second = bond_integrator;
        // Create scratch index with extra ghost cells
        auto var_db = SAMRAI::hier::VariableDatabase<NDIM>::getDatabase();
        d_bond_scr_idx = var_db->registerVariableAndContext(bond_var, var_db->getContext(d_object_name + "::CTX"), 4);
    }
    inline void setStressData(SAMRAI::tbox::Pointer<SAMRAI::hier::Variable<NDIM>> sig_var,
                              SAMRAI::tbox::Pointer<IBAMR::AdvDiffHierarchyIntegrator> sig_integrator)
    {
        d_var_integrator_pairs[3].first = sig_var;
        d_var_integrator_pairs[3].second = sig_integrator;
    }

    /*!
     * \name Methods to set patch data.
     */
    //\{

    /*!
     * \brief Indicates whether the concrete BoundPlateletSource object is
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

    inline void setKernel(const clot::Kernel kern)
    {
        d_kernel = kern;
    }

    inline void setSign(const bool positive)
    {
        d_sign = (positive) ? 1.0 : -1.0;
    }

    inline void setWIdx(const int w_idx)
    {
        d_w_idx = w_idx;
    }

    inline void registerBetaFcn(std::function<double(double, void*)> wrapper, void* beta)
    {
        d_beta_fcn = std::bind(wrapper, std::placeholders::_1, beta);
    }

private:
    double d_Kab = std::numeric_limits<double>::quiet_NaN();
    double d_Kaw = std::numeric_limits<double>::quiet_NaN();
    double d_nb_max = std::numeric_limits<double>::quiet_NaN();
    double d_nw_max = std::numeric_limits<double>::quiet_NaN();
    double d_nb = 1.0, d_nw = 1.0;
    double d_sign = 1.0;
    Kernel d_kernel = UNKNOWN_KERNEL;

    std::array<std::pair<SAMRAI::tbox::Pointer<SAMRAI::hier::Variable<NDIM>>,
                         SAMRAI::tbox::Pointer<IBTK::HierarchyIntegrator>>,
               4>
        d_var_integrator_pairs;
    int d_w_idx = IBTK::invalid_index;

    // Note we need our own scratch index with larger ghost width to compute the convolution.
    int d_phi_b_scr_idx = IBTK::invalid_index, d_bond_scr_idx = IBTK::invalid_index;

    // Beta function pointer
    std::function<double(double)> d_beta_fcn;
};
} // namespace clot
//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_clot_BoundPlateletSource
