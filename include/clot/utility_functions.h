/////////////////////////////// INCLUDE GUARD ////////////////////////////////

#ifndef included_clot_utility_functions
#define included_clot_utility_functions

/////////////////////////////// INCLUDES /////////////////////////////////////
#include <ibamr/config.h>

#include <ibtk/ibtk_utilities.h>

#include <CellData.h>
#include <CellIndex.h>

#include <functional>
#include <vector>

namespace clot
{
/*!
 * Computes the convolution of (alpha*a + beta*b) and the function phi at a given patch index, assuming phi has finite
 * support.
 *
 * Both patch data MUST have enough ghost cells for the width of the phi function. This assumes the phi function is
 * centered at idx and scaled by the grid spacing (i.e. phi(1) evaluates one grid cell away from idx).
 */
double convolution(double alpha,
                   SAMRAI::pdat::CellData<NDIM, double>* a_data,
                   double beta,
                   SAMRAI::pdat::CellData<NDIM, double>* b_data,
                   std::function<double(double)> phi,
                   unsigned int phi_width,
                   const SAMRAI::pdat::CellIndex<NDIM>& idx,
                   const double* const dx);

double convolution_mask(double alpha,
                        SAMRAI::pdat::CellData<NDIM, double>* a_data,
                        double beta,
                        SAMRAI::pdat::CellData<NDIM, double>* b_data,
                        std::function<double(double)> phi,
                        unsigned int phi_width,
                        const SAMRAI::pdat::CellIndex<NDIM>& idx,
                        const double* const dx,
                        const IBTK::VectorNd& x,
                        const std::array<IBTK::VectorNd, 2>& xbds);

/*!
 * Functions to convert the stress to the correct stress in the presence of non-divergence free velocity fields.
 */
//\{
/*!
 * \brief Template type Array must have a valid copy constructor and a member function data() that returns a pointer to
 * the contiguous underlying data.
 */
template <typename Array>
Array convertToStress(const Array& si, double bond, double S0, double R0);
IBTK::MatrixNd convertToStressMatrix(const IBTK::MatrixNd& si, double bond, double S0, double R0);

/*!
 * \brief Modifies the pointer provided to update the stress components.
 */
void convertToStressPtr(double* const si, double bond, double S0, double R0);
/*!
 * \brief Returns the factor that modifies the stress tensor.
 */
double modifiedStressFactor(double tr, double bond, double S0, double R0);

/*!
 * \brief Convert the cell data of the stress tensor on a given patch. An optional flag determines whether to convert
 * ghost cell data too.
 *
 * \note Conversion in the ghost cell region requires both sigma and bond ghost cells be given.
 */
SAMRAI::tbox::Pointer<SAMRAI::pdat::CellData<NDIM, double>>
convertToStress(const SAMRAI::pdat::CellData<NDIM, double>& sig_data,
                const SAMRAI::pdat::CellData<NDIM, double>& bond_data,
                SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM>> patch,
                double S0,
                double R0,
                bool fill_ghosts = false);
//\}

template <typename T>
inline T
string_to_enum(const std::string& /*val*/)
{
    TBOX_ERROR("UNSUPPORTED ENUM TYPE\n");
    return -1;
}

template <typename T>
inline std::string enum_to_string(T /*val*/)
{
    TBOX_ERROR("UNSUPPORTED ENUM TYPE\n");
    return "UNKNOWN";
}

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
double bspline2(const double x);
/*
 * Function that evaluates the BSpline 3 kernel.
 */
double bspline3(const double x);
/*
 * Function that evaluates the BSpline 4 kernel.
 */
double bspline4(const double x);

/*!
 * \brief Function that returns the kernel and the kernel width.
 */
inline std::pair<std::function<double(double)>, int>
getKernelAndWidth(Kernel val)
{
    if (val == BSPLINE_2) return std::make_pair(bspline2, 2);
    if (val == BSPLINE_3) return std::make_pair(bspline3, 3);
    if (val == BSPLINE_4) return std::make_pair(bspline4, 4);
    TBOX_ERROR("Unknown Kernel\n");
    // We've thrown an error, so just return whatever.
    return std::make_pair(bspline2, -1);
}
} // namespace clot

#include "clot/private/utility_functions_inc.h"

#endif // included_clot_utility_functions
