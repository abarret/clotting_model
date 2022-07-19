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

#if (NDIM == 2)
static int num_stress_comps = 3;
#else
static int num_stress_comps = 6;
#endif

Pointer<CellData<NDIM, double>>
convertToStress(const CellData<NDIM, double>& sig_data,
                const CellData<NDIM, double>& bond_data,
                Pointer<Patch<NDIM>> patch,
                const double S0,
                const double R0,
                bool fill_ghosts)
{
    IntVector<NDIM> gcw =
        fill_ghosts ? std::min(sig_data.getGhostCellWidth().max(), bond_data.getGhostCellWidth().max()) : 0;
    Pointer<CellData<NDIM, double>> ret_data = new CellData<NDIM, double>(patch->getBox(), sig_data.getDepth(), gcw);

    for (CellIterator<NDIM> ci(ret_data->getGhostBox()); ci; ci++)
    {
        const CellIndex<NDIM>& idx = ci();
        double tr = 0.0;
        for (int d = 0; d < NDIM; ++d) tr += sig_data(idx, d);
        tr = modifiedStressFactor(tr, bond_data(idx), S0, R0);
        for (int d = 0; d < NDIM; ++d) (*ret_data)(idx, d) = sig_data(idx, d) - tr;
        for (int d = NDIM; d < num_stress_comps; ++d) (*ret_data)(idx, d) = sig_data(idx, d);
    }
    return ret_data;
}

} // namespace clot
