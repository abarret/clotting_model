#include <clot/ClotParameters.h>

namespace clot
{
BoundClotParams::BoundClotParams(SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> db)
{
    phi_a_diff_coef = db->getDouble("phi_a_diff_coef");
    phi_b_diff_coef = db->getDouble("phi_b_diff_coef");
    uf_viscosity = db->getDouble("uf_viscosity");
    ub_viscosity = db->getDouble("ub_viscosity");
    drag_coef = db->getDouble("drag_coefficient");
    vol_pl = db->getDouble("vol_pl");
    stress_growth = db->getDouble("stress_growth");
    beta_0 = db->getDouble("beta_0");
    beta_1 = db->getDouble("beta_1");
    R0 = db->getDouble("r0");
    S0 = db->getDouble("s0");
    Kbb = db->getDouble("kbb");
    Kaw = db->getDouble("kaw");
    Kab = db->getDouble("kab");
    nw_max = db->getDouble("nw_max");
    nb_max = db->getDouble("nb_max");
    nb = db->getDouble("nb");
    nw = db->getDouble("nw");
}
} // namespace clot
