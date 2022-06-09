#include <clot/app_namespaces.h>
#include <clot/utility_functions.h>

#include <ibtk/ibtk_utilities.h>

#include <Box.h>

namespace clot
{
double
convolution(double alpha,
            CellData<NDIM, double>* a_data,
            double beta,
            CellData<NDIM, double>* b_data,
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
    TBOX_ASSERT(a_data || b_data);
    if (a_data) TBOX_ASSERT(a_data->getGhostBox() * int_box == int_box);
    if (b_data) TBOX_ASSERT(b_data->getGhostBox() * int_box == int_box);
#endif
    // Now compute convolution
    double convolve = 0.0;
    for (CellIterator<NDIM> ci(int_box); ci; ci++)
    {
        const CellIndex<NDIM>& n_idx = ci();
        double phi_weight = 1.0;
        for (int d = 0; d < NDIM; ++d) phi_weight *= phi(static_cast<double>(n_idx(d) - idx(d)));
        double val = 0.0;
        if (a_data) val += alpha * (*a_data)(n_idx);
        if (b_data) val += beta * (*b_data)(n_idx);
        convolve += val * phi_weight;
    }
    return convolve;
}

double
convolution_mask(double alpha,
                 CellData<NDIM, double>* a_data,
                 double beta,
                 CellData<NDIM, double>* b_data,
                 std::function<double(double)> phi,
                 unsigned int phi_width,
                 const CellIndex<NDIM>& idx,
                 const double* const dx,
                 const VectorNd& x,
                 const std::array<VectorNd, 2>& xbds)
{
    double convolve = convolution(alpha, a_data, beta, b_data, phi, phi_width, idx, dx);
    // Now mask the convolution with a heaviside function.
    // Determine if we are outside the box determined by xbds.
    if ((x(0) < xbds[0](0) || x(0) > xbds[1](0)) || (x(1) < xbds[0](1) || x(1) > xbds[1](1)))
        return 0.0;
    else
        return convolve;
}
} // namespace clot
