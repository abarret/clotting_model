#ifndef included_clot_DragForce
#define included_clot_DragForce

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
class DragForce : public IBTK::CartGridFunction
{
public:
    /*!
     * \brief Constructor.
     */
    DragForce(const std::string& object_name, SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> input_db);

    /*!
     * \brief Destructor.
     */
    ~DragForce();

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
     * Set the bound platelet data. Used to compute xi.
     */
    inline void setBoundPlateletData(SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double>> phi_var,
                                     SAMRAI::tbox::Pointer<IBAMR::AdvDiffHierarchyIntegrator> integrator)
    {
        d_phi_var = phi_var;
        d_phi_integrator = integrator;
        return;
    }

    /*!
     * Set the fluid velocity data.
     */
    inline void setVelocityData(SAMRAI::tbox::Pointer<SAMRAI::pdat::SideVariable<NDIM, double>> uf_var,
                                SAMRAI::tbox::Pointer<IBAMR::INSHierarchyIntegrator> integrator)
    {
        d_uf_var = uf_var;
        d_ins_integrator = integrator;
        return;
    }

    /*!
     * Set the bound velocity data.
     */
    inline void setBoundVelocityData(SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double>> ub_var,
                                     SAMRAI::tbox::Pointer<IBAMR::AdvDiffHierarchyIntegrator> integrator)
    {
        d_ub_var = ub_var;
        d_ub_integrator = integrator;
        return;
    }

    /*!
     * Set the time points used to determine the correct velocity context.
     */
    inline void setTimePoints(const double current_time, const double new_time)
    {
        d_current_time = current_time;
        d_new_time = new_time;
    }

private:
    DragForce() = delete;
    DragForce(const DragForce& from) = delete;
    DragForce& operator=(const DragForce& that) = delete;

    SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double>> d_phi_var, d_ub_var;
    SAMRAI::tbox::Pointer<SAMRAI::pdat::SideVariable<NDIM, double>> d_uf_var;

    SAMRAI::tbox::Pointer<IBAMR::AdvDiffHierarchyIntegrator> d_phi_integrator, d_ub_integrator;
    SAMRAI::tbox::Pointer<IBAMR::INSHierarchyIntegrator> d_ins_integrator;

    double d_vol_pl = std::numeric_limits<double>::quiet_NaN(), d_c3 = std::numeric_limits<double>::quiet_NaN();

    double d_current_time = std::numeric_limits<double>::quiet_NaN(),
           d_new_time = std::numeric_limits<double>::quiet_NaN();
};

} // namespace clot

#endif //#ifndef included_UFunction
