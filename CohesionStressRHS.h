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

#ifndef included_CohesionStressRHS
#define included_CohesionStressRHS

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <ibamr/config.h>

#include "ibamr/CFRelaxationOperator.h"
#include "ibamr/ibamr_enums.h"

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
#include <math.h>
#include <string>

/////////////////////////////// CLASS DEFINITION /////////////////////////////

namespace IBAMR
{
/*!
 * \brief Class CFOldroydBRelaxation is a concrete CFRelaxationOperator that computes the relaxation function for an
 * Oldroyd-B fluid model.
 */
class CohesionStressRHS : public IBAMR::CFRelaxationOperator
{
public:
    /*!
     * \brief This constructor reads in the parameters for the model from the input database.
     */
    CohesionStressRHS(const std::string& object_name, SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> input_db);

    /*!
     * \brief Default constructor.
     *
     * \note This constructor is not implemented and should not be used.
     */
    CohesionStressRHS() = delete;

    /*!
     * \brief Copy constructor.
     *
     * \note This constructor is not implemented and should not be used.
     */
    CohesionStressRHS(const CohesionStressRHS& from) = delete;

    /*!
     * \brief Assignment operator.
     *
     * \note This operator is not implemented and should not be used.
     */
    CohesionStressRHS& operator=(const CohesionStressRHS& that) = delete;

    /*!
     * \brief Empty destructor.
     */
    ~CohesionStressRHS() = default;

    /*!
     * \name Methods to set patch data.
     */
    //\{

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

    void registerBetaFcn(std::function<double(double, void*)> wrapper, void* beta);

    inline void setPlateletAIdx(const int phi_a_idx)
    {
        d_phi_a_idx = phi_a_idx;
    }

    inline void setPlateletUIdx(const int phi_u_idx)
    {
        d_phi_u_idx = phi_u_idx;
    }

    inline void setWccIdx(const int w_cc_idx)
    {
        d_W_cc_idx = w_cc_idx;
    }

    inline void setOmegaIdx(const int w_idx)
    {
        d_w_idx = w_idx;
    }

    inline void setZIdx(const int z_idx)
    {
        d_z_idx = z_idx;
    }

    inline void setKernel(const IBAMR::Kernel kern)
    {
        d_kernel = kern;
    }

private:
    int d_phi_a_idx = IBTK::invalid_index, d_z_idx = IBTK::invalid_index, d_w_idx = IBTK::invalid_index,
        d_phi_u_idx = IBTK::invalid_index, d_W_cc_idx = IBTK::invalid_index;
    double d_c4 = std::numeric_limits<double>::quiet_NaN();
    double d_a0 = std::numeric_limits<double>::quiet_NaN();
    double d_a0w = std::numeric_limits<double>::quiet_NaN();
    double d_w_mx = std::numeric_limits<double>::quiet_NaN();
    IBAMR::Kernel d_kernel = UNKNOWN_KERNEL;
    // Beta function pointer
    std::function<double(double)> d_beta_fcn;
};

} // Namespace IBAMR
#endif
