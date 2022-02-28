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

#ifndef included_BondSource
#define included_BondSource

#include <ibamr/AdvDiffHierarchyIntegrator.h>
#include "utility_functions.h"

#include "Box.h"
#include "CartesianPatchGeometry.h"
#include "CellData.h"
#include "CellIndex.h"
#include "CellVariable.h"
#include "Index.h"
#include "IntVector.h"
#include "Patch.h"
#include "PatchGeometry.h"
#include "PatchHierarchy.h"
#include "PatchLevel.h"
#include "Variable.h"
#include "tbox/Database.h"
#include "tbox/Pointer.h"
#include "tbox/Utilities.h"
#include "utility_functions.h"

#include <functional>

/////////////////////////////// CLASS DEFINITION /////////////////////////////
namespace IBAMR
{
    class BondSource : public IBTK::CartGridFunction
    {
        public:
        /*!
         * \brief Class constructor.
         */
        BondSource(SAMRAI::tbox::Pointer<SAMRAI::hier::Variable<NDIM>> phi_u_var,
                    SAMRAI::tbox::Pointer<SAMRAI::hier::Variable<NDIM>> phi_a_var,
                    SAMRAI::tbox::Pointer<SAMRAI::hier::Variable<NDIM>> z_var,
                    SAMRAI::tbox::Pointer<SAMRAI::hier::Variable<NDIM>> sig_var,
                    SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> input_db,
                    SAMRAI::tbox::Pointer<AdvDiffHierarchyIntegrator> adv_diff_hier_integrator);

        /*!
         * \brief Empty destructor.
         */
        ~BondSource() = default;

        /*!
         * \brief Indicates whether the concrete BondSource object is
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

        void registerBetaFcn(std::function<double(double, void*)> wrapper, void* beta);

        inline void setWIdx(const int w_idx)
        {
            d_w_idx = w_idx;
        }
        private:
        // INPUT DB VARS (I assume eq. 8-10 in new model correspond to bond source)
        double d_a0 = std::numeric_limits<double>::quiet_NaN();
        double d_a0w = std::numeric_limits<double>::quiet_NaN();
        double d_w_mx = std::numeric_limits<double>::quiet_NaN();

        // wall index
        int d_w_idx = IBTK::invalid_index;

        // Beta function pointer
        std::function<double(double)> d_beta_fcn;

        // OTHER THINGS
        BondSource(const BondSource& from) = delete;
        BondSource& operator=(const BondSource& that) = delete;
        SAMRAI::tbox::Pointer<SAMRAI::hier::Variable<NDIM>> d_phi_u_var, d_phi_a_var, d_z_var, d_sig_var;
        SAMRAI::tbox::Pointer<AdvDiffHierarchyIntegrator> d_adv_diff_hier_integrator;
    }
}

#endif