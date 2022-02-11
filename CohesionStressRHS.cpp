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

#include "ibamr/app_namespaces.h" // IWYU pragma: keep

#include "ibtk/ibtk_utilities.h"

#include "CellData.h"
#include "CellIterator.h"
#include "Patch.h"
#include "tbox/Database.h"

// Local headers
#include "CohesionStressRHS.h"

// Namespace
namespace IBAMR
{
CohesionStressRHS::CohesionStressRHS(const std::string& object_name, Pointer<Database> input_db)
    : CFRelaxationOperator(object_name, input_db)
{
    // Get values from inputdb
    d_c4 = input_db->getDouble("c4");
    // a0 Constants
    d_a0 = input_db->detDouble("a0");
    d_a0w = input_db->getDouble("a0w");
    // w constant
    d_w_mx = input_db->getDouble("wmax")
    return;
} // Constructor

void
CohesionStressRHS::setDataOnPatch(const int data_idx,
                                  Pointer<Variable<NDIM>> /*var*/,
                                  Pointer<Patch<NDIM>> patch,
                                  const double /*data_time*/,
                                  const bool initial_time,
                                  Pointer<PatchLevel<NDIM>> /*patch_level*/)
{
    const Box<NDIM>& patch_box = patch->getBox();
    //Pointer<CartesianPatchGeometry<NDIM> > pgeom = patch->getPatchGeometry();
    //const double* const dx = pgeom->getDx();
    // SET THE INDICES!!!!
    Pointer<CellData<NDIM, double>> ret_data = patch->getPatchData(data_idx);
    Pointer<CellData<NDIM, double>> in_data = patch->getPatchData(d_W_cc_idx); //this is sigma?
    Pointer<CellData<NDIM, double>> phi_a_data = patch->getPatchData(d_phi_a_idx);
    Pointer<CellData<NDIM, double>> phi_u_data = patch->getPatchData(d_phi_u_idx);
    Pointer<CellData<NDIM, double>> z_data = patch->getPatchData(d_z_idx);
    Pointer<CellData<NDIM, double>> w_data = patch->getPatchData(d_w_idx);
    auto phi_fcn = IBAMR::getKernelAndWidth(d_kernel);
    ret_data->fillAll(0.0);
    if (initial_time) return;
    for (CellIterator<NDIM> i(patch_box); i; i++)
    {
        const CellIndex<NDIM>& idx = i();
        // Compute the ODE terms for stress here
        // C4 * alpha * I - beta * sigma
        // compute R2 R3 and R4 -> alpha
        // where beta = beta(<y>)
        // where <y> = sqrt(2*Tr(sig) / (S0*z))
        
        // Index variables
        double phi_a = (*phi_a_data)(idx);
        double phi_u = (*phi_u_data)(idx);
        double w = (*w_data)(idx);
        double z = (*z_data)(idx);
        //const double eta_b = IBAMR::convolution(d_n_b_mx, phi_b, -2.0, z_data, phi_fcn.first, phi_fcn.second, idx, dx);
        // Compute the source terms
        const double R2 = d_a0 * phi_a * phi_a;
        const double R4 = d_a0w * (d_w_mx - w) * phi_u;
        const double alpha = R2 + R4;
#if (NDIM == 2)
        const double trace = (*in_data)(idx, 0) + (*in_data)(idx, 1);
        const double y_brackets = std::sqrt(2.0 * trace / (z + 1.0e-12));
        const double beta = d_beta_fcn(y_brackets);

        (*ret_data)(idx, 0) = d_c4 * alpha - beta * (*in_data)(idx, 0);
        (*ret_data)(idx, 1) = d_c4 * alpha - beta * (*in_data)(idx, 1);
        (*ret_data)(idx, 2) = -beta * (*in_data)(idx, 2);
#endif
#if (NDIM == 3)
        const double trace = (*in_data)(idx, 0) + (*in_data)(idx, 1) + (*in_data)(idx, 2);
        const double y_brackets = std::sqrt(2.0 * trace / (z+ 1.0e-12));
        const double beta = d_beta_fcn(y_brackets);

        (*ret_data)(idx, 0) = d_c4 * alpha - beta * (*in_data)(idx, 0);
        (*ret_data)(idx, 1) = d_c4 * alpha - beta * (*in_data)(idx, 1);
        (*ret_data)(idx, 2) = d_c4 * alpha - beta * (*in_data)(idx, 2);
        (*ret_data)(idx, 3) = -beta * (*in_data)(idx, 3);
        (*ret_data)(idx, 4) = -beta * (*in_data)(idx, 4);
        (*ret_data)(idx, 5) = -beta * (*in_data)(idx, 5);
#endif
    }
} // setDataOnPatch

void
CohesionStressRHS::registerBetaFcn(std::function<double(double, void*)> wrapper, void* beta)
{
    // We set beta here
    d_beta_fcn = std::bind(wrapper, std::placeholders::_1, beta);
}

} // namespace IBAMR
