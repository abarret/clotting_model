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
#include "CohesionStressRHS.h"
#include "Patch.h"
#include "tbox/Database.h"

// Namespace
namespace IBAMR
{
CohesionStressRHS::CohesionStressRHS(const std::string& object_name, Pointer<Database> input_db)
    : CFRelaxationOperator(object_name, input_db)
{
    // Get values from inputdb
    d_a2 = input_db->getDouble("a2");
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
    Pointer<CellData<NDIM, double>> ret_data = patch->getPatchData(data_idx);
    Pointer<CellData<NDIM, double>> in_data = patch->getPatchData(d_W_cc_idx);
    Pointer<CellData<NDIM, double>> phi_data = patch->getPatchData(d_phi_idx);
    Pointer<CellData<NDIM, double>> z_data = patch->getPatchData(d_z_idx);
    ret_data->fillAll(0.0);
    if (initial_time) return;
    //const double l_inv = 1.0 / d_lambda; not used in this calculation
    for (CellIterator<NDIM> i(patch_box); i; i++)
    {
        const CellIndex<NDIM>& idx = i();
	// Compute the ODE terms for stress here
#if (NDIM == 2)
    double trace = (*in_data)(idx, 0) + (*in_data)(idx, 1); 
    (*ret_data)(idx, 0) = d_a2 * (*phi_data)(idx) * (*phi_data)(idx) - d_beta_fcn(trace / (*z_data)(idx)) * (*in_data)(idx, 0);
    (*ret_data)(idx, 1) = d_a2 * (*phi_data)(idx) * (*phi_data)(idx) - d_beta_fcn(trace / (*z_data)(idx)) * (*in_data)(idx, 1);
    (*ret_data)(idx, 2) = -d_beta_fcn(trace / (*z_data)(idx)) * (*in_data)(idx, 2);
#endif
#if (NDIM == 3)
    double trace = (*in_data)(idx, 0) + (*in_data)(idx, 1) + (*in_data)(idx, 2);
    (*ret_data)(idx, 0) = d_a2 * (*phi_data)(idx) * (*phi_data)(idx) - d_beta_fcn(trace / (*z_data)(idx)) * (*in_data)(idx, 0);
    (*ret_data)(idx, 1) = d_a2 * (*phi_data)(idx) * (*phi_data)(idx) - d_beta_fcn(trace / (*z_data)(idx)) * (*in_data)(idx, 1);
    (*ret_data)(idx, 2) = d_a2 * (*phi_data)(idx) * (*phi_data)(idx) - d_beta_fcn(trace / (*z_data)(idx)) * (*in_data)(idx, 2);
    (*ret_data)(idx, 3) = -d_beta_fcn(trace / (*z_data)(idx)) * (*in_data)(idx, 3);
    (*ret_data)(idx, 4) = -d_beta_fcn(trace / (*z_data)(idx)) * (*in_data)(idx, 4);
    (*ret_data)(idx, 5) = -d_beta_fcn(trace / (*z_data)(idx)) * (*in_data)(idx, 5);
#endif
    }
} // setDataOnPatch


void CohesionStressRHS::registerBetaFcn(std::function<double(double, void*)> wrapper, void *beta) {
    // We set beta here
    d_beta_fcn = std::bind(wrapper, std::placeholders::_1, beta);
}

} // namespace IBAMR
