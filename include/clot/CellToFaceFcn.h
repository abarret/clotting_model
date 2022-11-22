#ifndef included_clot_CellToFaceFcn
#define included_clot_CellToFaceFcn

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
class CellToFaceFcn : public IBTK::CartGridFunction
{
public:
    /*!
     * \brief Constructor.
     */
    CellToFaceFcn(std::string object_name,
                  SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double>> Q_var,
                  SAMRAI::tbox::Pointer<IBAMR::AdvDiffHierarchyIntegrator> integrator);

    /*!
     * \brief Destructor.
     */
    ~CellToFaceFcn() = default;

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
                                 int coarsest_ln = IBTK::invalid_level_number,
                                 int finest_ln = IBTK::invalid_level_number) override;

    /*!
     * Set the data on the patch interior to some initial values.
     */
    void setDataOnPatch(int data_idx,
                        SAMRAI::tbox::Pointer<SAMRAI::hier::Variable<NDIM>> var,
                        SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM>> patch,
                        double data_time,
                        bool initial_time = false,
                        SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM>> level = nullptr) override;

private:
    CellToFaceFcn() = delete;
    CellToFaceFcn(const CellToFaceFcn& from) = delete;
    CellToFaceFcn& operator=(const CellToFaceFcn& that) = delete;

    SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double>> d_Q_var;

    SAMRAI::tbox::Pointer<IBAMR::AdvDiffHierarchyIntegrator> d_integrator;
};

} // namespace clot

#endif //#ifndef included_UFunction
