/////////////////////////////// INCLUDE GUARD ////////////////////////////////

#ifndef included_utility_functions
#define included_utility_functions

/////////////////////////////// INCLUDES /////////////////////////////////////
#include <ibamr/config.h>

#include <CellData.h>
#include <CellIndex.h>

#include <functional>
namespace IBAMR
{
/*!
 * Computes the convolution of (alpha*a + beta*b) and the function phi at a given patch index, assuming phi has finite
 * support.
 *
 * Both patch data MUST have enough ghost cells for the width of the phi function. This assumes the phi function is
 * centered at idx and scaled by the grid spacing.
 */
double convolution(double alpha,
                   const SAMRAI::pdat::CellData<NDIM, double>& a_data,
                   double beta,
                   const SAMRAI::pdat::CellData<NDIM, double>& b_data,
                   std::function<double(double)> phi,
                   unsigned int phi_width,
                   const SAMRAI::pdat::CellIndex<NDIM>& idx,
                   const double* const dx);
} // namespace IBAMR

#endif // included_utility_functions
