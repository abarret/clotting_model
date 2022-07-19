#ifndef included_clot_utility_functions_inl
#define included_clot_utility_functions_inl

#include <clot/utility_functions.h>

namespace clot
{
inline double
bspline2(double x)
{
    x = std::abs(x);
    if (x <= 1.0)
        return 1.0 - x;
    else
        return 0.0;
}

inline double
bspline3(double x)
{
    x = std::abs(x);
    const double r = x + 1.5;
    const double r2 = r * r;
    if (x <= 0.5)
        return 0.5 * (-2.0 * r2 + 6.0 * r - 3.0);
    else if (x <= 1.5)
        return 0.5 * (r2 - 6.0 * r + 9.0);
    else
        return 0.0;
}

inline double
bspline4(double x)
{
    x = std::abs(x);
    const double r = x + 2.0;
    const double r2 = r * r;
    const double r3 = r2 * r;
    if (x <= 1.0)
        return (1.0 / 6.0) * (3.0 * r3 - 24.0 * r2 + 60.0 * r - 44.0);
    else if (x <= 2.0)
        return (1.0 / 6.0) * (-r3 + 12.0 * r2 - 48.0 * r + 64.0);
    else
        return 0.0;
}

template <typename Array>
inline Array
convertToStress(const Array& si, const double bond, const double S0, const double R0)
{
    Array ret = si;
    return convertToStressPtr(ret.data(), bond, S0, R0);
}

inline void
convertToStressPtr(double* const si, const double bond, const double S0, const double R0)
{
    if (bond < 1.0e-8) return;
    double tr = 0.0;
    for (int d = 0; d < NDIM; ++d) tr += si[d];
    tr = modifiedStressFactor(tr, bond, S0, R0);
    for (int d = 0; d < NDIM; ++d) si[d] -= tr;
}

inline IBTK::MatrixNd
convertToStressMatrix(const IBTK::MatrixNd& si, const double bond, const double S0, const double R0)
{
    IBTK::MatrixNd ret = si;
    if (bond < 1.0e-8) return ret;
    double tr = modifiedStressFactor(si.trace(), bond, S0, R0);
    ret = si - tr * IBTK::MatrixNd::Ones();
    return ret;
}

inline double
modifiedStressFactor(double tr, const double bond, const double S0, const double R0)
{
    if (bond < 1.0e-8) return tr;
    tr = -bond * S0 * R0 * std::sqrt(2.0 * tr / (S0 * bond + 1.0e-8));
    return tr;
}
} // namespace clot
#endif
