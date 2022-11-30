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

#ifndef included_clot_BoundExtraStressForcing
#define included_clot_BoundExtraStressForcing

#include <clot/ClotParameters.h>
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
 * Class BoundExtraStressForcing is a concrete CartGridFunction that implements the growth and decay terms for the bond
 * variable.
 *
 * Note this function should only be used for the model that has unactivated and activated platelets.
 */
class BoundExtraStressForcing : public IBTK::CartGridFunction
{
public:
    /*!
     * \brief Class constructor.
     */
    BoundExtraStressForcing(std::string object_name, BoundClotParams clot_params);

    /*!
     * \brief Empty destructor.
     */
    ~BoundExtraStressForcing() = default;

    /*!
     * \brief Deleted constructors
     */
    //\{
    BoundExtraStressForcing(const BoundExtraStressForcing& from) = delete;
    BoundExtraStressForcing& operator=(const BoundExtraStressForcing& that) = delete;
    //\}

    /*!
     * \brief Indicates whether the concrete BoundExtraStressForcing object is
     * time-dependent.
     */
    bool isTimeDependent() const override;

    /*!
     * \brief Evaluate the function on the patch interiors on the specified
     * levels of the patch hierarchy.
     */
    void setDataOnPatchHierarchy(int data_idx,
                                 SAMRAI::tbox::Pointer<SAMRAI::hier::Variable<NDIM>> var,
                                 SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM>> hierarchy,
                                 double data_time,
                                 bool initial_time = false,
                                 int coarsest_ln = -1,
                                 int finest_ln = -1) override;

    /*!
     * \brief Evaluate the function on the patch interior.
     */
    void setDataOnPatch(int data_idx,
                        SAMRAI::tbox::Pointer<SAMRAI::hier::Variable<NDIM>> var,
                        SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM>> patch,
                        double data_time,
                        bool initial_time = false,
                        SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM>> patch_level =
                            SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM>>(NULL)) override;

    /*!
     * \brief Methods to set the inputs
     */
    //\{
    void setBondData(SAMRAI::tbox::Pointer<SAMRAI::hier::Variable<NDIM>> var,
                     SAMRAI::tbox::Pointer<IBTK::HierarchyIntegrator> integrator)
    {
        d_var_integrator_pairs[0].first = var;
        d_var_integrator_pairs[0].second = integrator;
    }

    void setStressData(SAMRAI::tbox::Pointer<SAMRAI::hier::Variable<NDIM>> var,
                       SAMRAI::tbox::Pointer<IBTK::HierarchyIntegrator> integrator)
    {
        d_var_integrator_pairs[1].first = var;
        d_var_integrator_pairs[1].second = integrator;
    }
    //\}

private:
    void
    findFactor(SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM>> hierarchy, int coarsest_ln, int finest_ln);

    const BoundClotParams& d_clot_params;

    std::array<std::pair<SAMRAI::tbox::Pointer<SAMRAI::hier::Variable<NDIM>>,
                         SAMRAI::tbox::Pointer<IBTK::HierarchyIntegrator>>,
               2>
        d_var_integrator_pairs;

    // Scratch index for stress
    SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double>> d_extra_stress_var;
    int d_extra_stress_idx = IBTK::invalid_index;
};
} // namespace clot

#endif
