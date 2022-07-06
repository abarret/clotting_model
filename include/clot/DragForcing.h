// ---------------------------------------------------------------------
//
// Copyright (c) 2014 - 2021 by the IBAMR developers
// All rights reserved.
//
// This file is part of IBAMR.
//
// IBAMR is free software and is distributed under the 3-clause BSD
// license. The full text of the license can be found in the file
// COPYRIGHT at the top level directory of IBAMR.
//
// ---------------------------------------------------------------------

#ifndef included_clot_DragForcing
#define included_clot_DragForcing

/////////////////////////////// INCLUDES /////////////////////////////////////
#include <ibamr/AdvDiffHierarchyIntegrator.h>
#include <ibamr/INSHierarchyIntegrator.h>

#include <ibtk/CartGridFunction.h>

namespace clot
{
/////////////////////////////// CLASS DEFINITION /////////////////////////////

/*!
 * \brief Class DragForcing provides forcing for the momentum equation for the fluid phase. It computes drag between the
 * fluid velocity and the bound velocity.
 */
class DragForcing : public IBTK::CartGridFunction
{
public:
    /*!
     * \brief Class constructor.
     */
    DragForcing(std::string object_name, SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> input_db);

    /*!
     * \brief Deleted constructors
     */
    //\{
    DragForcing() = delete;
    DragForcing(const DragForcing& from) = delete;
    DragForcing& operator=(const DragForcing& that) = delete;
    //\}

    /*!
     * \brief Empty destructor.
     */
    virtual ~DragForcing();

    inline void setBoundVelocityData(SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double>> sig_var,
                                     SAMRAI::tbox::Pointer<IBAMR::AdvDiffHierarchyIntegrator> integrator)
    {
        d_ub_var = sig_var;
        d_ub_integrator = integrator;
        return;
    }

    /*!
     * Set the bound platelet data. Used to compute xi.
     */
    inline void setPlateletData(SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double>> phi_var,
                                SAMRAI::tbox::Pointer<IBAMR::AdvDiffHierarchyIntegrator> integrator)
    {
        d_phi_var = phi_var;
        d_phi_integrator = integrator;
        return;
    }

    /*!
     * Set the fluid velocity data.
     */
    inline void setFluidVelocityData(SAMRAI::tbox::Pointer<SAMRAI::pdat::FaceVariable<NDIM, double>> uf_var,
                                     SAMRAI::tbox::Pointer<IBAMR::INSHierarchyIntegrator> integrator)
    {
        d_uf_var = uf_var;
        d_ins_integrator = integrator;
        return;
    }

    /*!
     * \name Methods to set patch data.
     */
    //\{

    /*!
     * \brief Indicates whether the concrete DragForcing object is
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

private:
    SAMRAI::tbox::Pointer<SAMRAI::hier::Variable<NDIM>> d_uf_var, d_ub_var, d_phi_var;
    SAMRAI::tbox::Pointer<IBAMR::AdvDiffHierarchyIntegrator> d_ub_integrator, d_phi_integrator;
    SAMRAI::tbox::Pointer<IBAMR::INSHierarchyIntegrator> d_ins_integrator;

    double d_C3 = std::numeric_limits<double>::quiet_NaN(), d_vol_pl = std::numeric_limits<double>::quiet_NaN();
};

//////////////////////////////////////////////////////////////////////////////
} // namespace clot
#endif //#ifndef included_DragForcing
