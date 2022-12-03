#ifndef included_clot_ClotParameters
#define included_clot_ClotParameters

#include "tbox/Database.h"

#include <limits>

namespace clot
{
struct BoundClotParams
{
public:
    BoundClotParams(SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> db);

    // Diffusion Coefficients
    double phi_a_diff_coef = std::numeric_limits<double>::quiet_NaN(),
           phi_b_diff_coef = std::numeric_limits<double>::quiet_NaN();

    // Momentum coefficients
    double uf_viscosity = std::numeric_limits<double>::quiet_NaN(),
           ub_viscosity = std::numeric_limits<double>::quiet_NaN();
    double drag_coef = std::numeric_limits<double>::quiet_NaN();

    // Platelet parameters
    double vol_pl = std::numeric_limits<double>::quiet_NaN();

    // Stress growth constant
    double stress_growth = std::numeric_limits<double>::quiet_NaN();

    // Bond breaking constants
    double beta_0 = std::numeric_limits<double>::quiet_NaN(), beta_1 = std::numeric_limits<double>::quiet_NaN();
    double R0 = std::numeric_limits<double>::quiet_NaN(), S0 = std::numeric_limits<double>::quiet_NaN();

    // Transition rates
    double Kbb = std::numeric_limits<double>::quiet_NaN(), Kaw = std::numeric_limits<double>::quiet_NaN(),
           Kab = std::numeric_limits<double>::quiet_NaN();

    // Max bonds per platelet or wall site
    double nw_max = std::numeric_limits<double>::quiet_NaN(), nb_max = std::numeric_limits<double>::quiet_NaN();

    // Average bonds per platelet
    double nb = std::numeric_limits<double>::quiet_NaN(), nw = std::numeric_limits<double>::quiet_NaN();
};

} // namespace clot
#endif
