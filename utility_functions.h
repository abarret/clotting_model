/////////////////////////////// INCLUDE GUARD ////////////////////////////////

#ifndef included_utility_functions
#define included_utility_functions

/////////////////////////////// INCLUDES /////////////////////////////////////
#include <ibamr/config.h>

#include <ibamr/ibamr_enums.h>

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
 * centered at idx and scaled by the grid spacing (i.e. phi(1) evaluates one grid cell away from idx).
 */
double convolution(double alpha,
                   const SAMRAI::pdat::CellData<NDIM, double>& a_data,
                   double beta,
                   const SAMRAI::pdat::CellData<NDIM, double>& b_data,
                   std::function<double(double)> phi,
                   unsigned int phi_width,
                   const SAMRAI::pdat::CellIndex<NDIM>& idx,
                   const double* const dx);

enum Kernel
{
    BSPLINE_2,
    BSPLINE_3,
    BSPLINE_4,
    UNKNOWN_KERNEL = -1
};

template <>
inline Kernel
string_to_enum<Kernel>(const std::string& val)
{
    if (strcasecmp(val.c_str(), "BSPLINE_2")) return BSPLINE_2;
    if (strcasecmp(val.c_str(), "BSPLINE_3")) return BSPLINE_3;
    if (strcasecmp(val.c_str(), "BSPLINE_4")) return BSPLINE_4;
    return UNKNOWN_KERNEL;
}

template <>
inline std::string
enum_to_string<Kernel>(Kernel val)
{
    if (val == BSPLINE_2) return "BSPLINE_2";
    if (val == BSPLINE_3) return "BSPLINE_3";
    if (val == BSPLINE_4) return "BSPLINE_4";
    return "UNKNOWN_KERNEL";
}

/*
 * Function that evaluates the BSpline 2 kernel. Note that this is just a hat function.
 */
double bspline_2(const double x);
/*
 * Function that evaluates the BSpline 3 kernel.
 */
double bspline_3(const double x);
/*
 * Function that evaluates the BSpline 4 kernel.
 */
double bspline_4(const double x);

/*!
 * \brief Function that returns the kernel and the kernel width.
 */
inline std::pair<std::function<double(double)>, int>
getKernelAndWidth(Kernel val)
{
    if (val == BSPLINE_2) return std::make_pair(bspline_2, 2);
    if (val == BSPLINE_3) return std::make_pair(bspline_3, 3);
    if (val == BSPLINE_4) return std::make_pair(bspline_4, 4);
    TBOX_ERROR("Unknown Kernel\n");
    // We've thrown an error, so just return whatever.
    return std::make_pair(bspline_2, -1);
}
} // namespace IBAMR

#include "utility_functions.I"

#endif // included_utility_functions
