#ifndef included_utility_functions_inl
#define included_utility_functions_inl

#include "utility_functions.h"

namespace IBAMR
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

} // namespace IBAMR
#endif
