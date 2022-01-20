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
    const double d_a2 = input_db->getDouble("a2"); // TODO DELETE THIS
    const double d_s0 = input_db->getDouble("s0");
    const double d_c4 = input_db->getDouble("c4");
    // K Constants
    const double d_Kab = input_db->getDouble("Kab"); // change the get strings to fix whatever the true label is
    const double d_Kbb = input_db->getDouble("Kbb");
    const double d_Kaw = input_db->getDouble("Kaw");
    // n constants
    const double d_n_b_mx = input_db->getDouble("nbmax")
    const double d_n_w_mx = input_db->getDouble("nwmax")
    // w constant
    const double d_w_mx = input_db->getDouble("wmax")
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
    // SET THE INDICES!!!!
    Pointer<CellData<NDIM, double>> ret_data = patch->getPatchData(data_idx);
    Pointer<CellData<NDIM, double>> in_data = patch->getPatchData(d_W_cc_idx); //this is sigma?
    Pointer<CellData<NDIM, double>> phi_a_data = patch->getPatchData(d_phi_a_idx);
    Pointer<CellData<NDIM, double>> phi_b_data = patch->getPatchData(d_phi_b_idx);
    Pointer<CellData<NDIM, double>> z_data = patch->getPatchData(d_z_idx);
    Pointer<CellData<NDIM, double>> w_data = patch->getPatchData(d_w_idx);
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
        const double d_eta_b = 1.0; // placeholder since idk how to compute it
        // Index variables
        double d_phi_a = (*phi_a_data)(idx);
        double d_phi_b = (*phi_b_data)(idx);
        double d_w = (*w_data)(idx);
        double d_z = (*z_data)(idx);
        // Compute the source terms
        const double d_R2 = d_Kab * d_n_b_mx * d_phi_a * d_eta_b;
        const double d_R3 = d_Kbb * (d_n_b_mx * d_phi_b - 2.0 * d_z) * (d_n_b_mx * d_phi_b - 2.0 * d_z);
        const double d_R4 = d_Kaw * d_n_w_mx * (d_w_mx - d_w) * d_n_b_mx * d_phi_a;
        const double alpha = d_R2 + d_R3 + d_R4;
#if (NDIM == 2)
        const double trace = (*in_data)(idx, 0) + (*in_data)(idx, 1);
        const double d_y_brackets = std::sqrt(2.0 * trace / (d_s0 * d_z));

        (*ret_data)(idx, 0) = d_c4 * alpha - d_beta_fcn(d_y_brackets) * (*in_data)(idx, 0);
        (*ret_data)(idx, 1) = d_c4 * alpha - d_beta_fcn(d_y_brackets) * (*in_data)(idx, 1);
        (*ret_data)(idx, 2) = -d_beta_fcn(d_y_brackets) * (*in_data)(idx, 2);
#endif
#if (NDIM == 3)
        const double trace = (*in_data)(idx, 0) + (*in_data)(idx, 1) + (*in_data)(idx, 2);
        const double d_y_brackets = std::sqrt(2.0 * trace / (d_s0 * d_z));

        (*ret_data)(idx, 0) = d_c4 * alpha - d_beta_fcn(d_y_brackets) * (*in_data)(idx, 0);
        (*ret_data)(idx, 1) = d_c4 * alpha - d_beta_fcn(d_y_brackets) * (*in_data)(idx, 1);
        (*ret_data)(idx, 2) = d_c4 * alpha - d_beta_fcn(d_y_brackets) * (*in_data)(idx, 2);
        (*ret_data)(idx, 3) = -d_beta_fcn(d_y_brackets) * (*in_data)(idx, 3);
        (*ret_data)(idx, 4) = -d_beta_fcn(d_y_brackets) * (*in_data)(idx, 4);
        (*ret_data)(idx, 5) = -d_beta_fcn(d_y_brackets) * (*in_data)(idx, 5);
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
