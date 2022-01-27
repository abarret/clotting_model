#include <ibamr/app_namespaces.h>

#include <ibtk/ibtk_utilities.h>

#include <Box.h>

// Local includes
#include "utility_functions.h"

namespace IBAMR
{
double
convolution(double alpha,
            const CellData<NDIM, double>& a_data,
            double beta,
            const CellData<NDIM, double>& b_data,
            std::function<double(double)> phi,
            unsigned int phi_width,
            const CellIndex<NDIM>& idx,
            const double* const dx)
{
    // Box over which integration is performed
    Box<NDIM> int_box(idx, idx);
    int_box.grow(phi_width);
#ifndef NDEBUG
    // Error checking. Need ghost data for patch data.
    TBOX_ASSERT(a_data.getGhostBox() * int_box == int_box);
    TBOX_ASSERT(b_data.getGhostBox() * int_box == int_box);
#endif
    // Now compute convolution
    double convolve = 0.0;
    for (CellIterator<NDIM> ci(int_box); ci; ci++)
    {
        const CellIndex<NDIM>& n_idx = ci();
        double phi_weight = 1.0;
        for (int d = 0; d < NDIM; ++d) phi_weight *= phi(static_cast<double>(n_idx(d) - idx(d))) * dx[d];
        convolve += (alpha * a_data(n_idx) + beta * b_data(n_idx)) * phi_weight;
    }
    return convolve;
}
} // namespace IBAMR
