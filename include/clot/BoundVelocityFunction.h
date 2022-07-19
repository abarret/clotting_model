#ifndef included_clot_BoundVelocityFunction
#define included_clot_BoundVelocityFunction

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <ibamr/AdvDiffHierarchyIntegrator.h>
#include <ibamr/INSHierarchyIntegrator.h>

#include <ibtk/CartGridFunction.h>
#include <ibtk/ibtk_utilities.h>

#include <CellVariable.h>
#include <FaceVariable.h>

/////////////////////////////// CLASS DEFINITION /////////////////////////////

namespace clot
{
/*!
 * \brief Computes the bound velocity field u_b
 */
class BoundVelocityFunction : public IBTK::CartGridFunction
{
public:
    /*!
     * \brief Constructor.
     */
    BoundVelocityFunction(const std::string& object_name, SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> input_db);

    /*!
     * \brief Destructor.
     */
    ~BoundVelocityFunction();

    /*!
     * Indicates whether the concrete CartGridFunction object is time dependent.
     */
    bool isTimeDependent() const override
    {
        return true;
    }

    void setDataOnPatchHierarchy(int data_idx,
                                 SAMRAI::tbox::Pointer<SAMRAI::hier::Variable<NDIM>> var,
                                 SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM>> hierarchy,
                                 double data_time,
                                 bool initial_time = false,
                                 int coarsest_ln = -1,
                                 int finest_ln = -1) override;

    /*!
     * Set the data on the patch interior to some initial values.
     */
    void setDataOnPatch(int data_idx,
                        SAMRAI::tbox::Pointer<SAMRAI::hier::Variable<NDIM>> var,
                        SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM>> patch,
                        double data_time,
                        bool initial_time = false,
                        SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM>> level = nullptr) override;

    /*!
     * Set the extra stress tensor data.
     *
     * \note This must be done outside of the constructor because of how the velocity function is read by the
     * integrator.
     */
    inline void setSigmaData(SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double>> sig_var,
                             SAMRAI::tbox::Pointer<IBAMR::AdvDiffHierarchyIntegrator> integrator)
    {
        d_sig_var = sig_var;
        d_sig_integrator = integrator;
        return;
    }

    /*!
     * Set the bound platelet data. Used to compute xi.
     */
    inline void setPlateletData(SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double>> phi_var,
                                SAMRAI::tbox::Pointer<IBAMR::AdvDiffHierarchyIntegrator> integrator)
    {
        d_phi_var = phi_var;
        d_phi_integrator = integrator;
        return;
    }

    /*!
     * Set the fluid velocity data.
     */
    inline void setVelocityData(SAMRAI::tbox::Pointer<SAMRAI::pdat::FaceVariable<NDIM, double>> uf_var,
                                SAMRAI::tbox::Pointer<IBAMR::INSHierarchyIntegrator> integrator)
    {
        d_uf_var = uf_var;
        d_ins_integrator = integrator;
        return;
    }

    /*!
     * Set the bond data.
     */
    inline void setBondData(SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double>> bond_var,
                            SAMRAI::tbox::Pointer<IBAMR::AdvDiffHierarchyIntegrator> integrator)
    {
        d_bond_var = bond_var;
        d_bond_integrator = integrator;
        return;
    }

private:
    BoundVelocityFunction() = delete;
    BoundVelocityFunction(const BoundVelocityFunction& from) = delete;
    BoundVelocityFunction& operator=(const BoundVelocityFunction& that) = delete;

    SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double>> d_sig_var, d_phi_var, d_bond_var;
    SAMRAI::tbox::Pointer<SAMRAI::pdat::FaceVariable<NDIM, double>> d_uf_var;

    SAMRAI::tbox::Pointer<IBAMR::AdvDiffHierarchyIntegrator> d_sig_integrator, d_phi_integrator, d_bond_integrator;
    SAMRAI::tbox::Pointer<IBAMR::INSHierarchyIntegrator> d_ins_integrator;

    double d_thresh = 1.0e-8;
    double d_vol_pl = std::numeric_limits<double>::quiet_NaN(), d_c3 = std::numeric_limits<double>::quiet_NaN(),
           d_S0 = std::numeric_limits<double>::quiet_NaN(), d_R0 = std::numeric_limits<double>::quiet_NaN();
};

} // namespace clot

#endif //#ifndef included_UFunction
