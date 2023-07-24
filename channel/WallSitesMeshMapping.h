#ifndef included_clot_WallSitesMeshMapping
#define included_clot_WallSitesMeshMapping
#include <ADS/GeneralBoundaryMeshMapping.h>

#include <ibtk/FEDataManager.h>

#include <libmesh/boundary_mesh.h>
#include <libmesh/mesh.h>

namespace clot
{
/*!
 * WallSitesMeshMapping is a class that generalizes the notion of a boundary mesh. It maintains an EquationSystems
 * object on the mesh, handles restarts of Lagrangian data, and maintains a FEDataManager for the object.
 * Implementations for this object should define how the object moves. A default implementation of no motion is
 * included.
 */
class WallSitesMeshMapping : public ADS::GeneralBoundaryMeshMapping
{
public:
    /*!
     * \brief Constructor that takes in a boundary mesh. Note that WallSitesMeshMapping assumes ownership of the
     * mesh.
     */
    WallSitesMeshMapping(std::string object_name,
                         SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> input_db,
                         std::string W_sys_name,
                         libMesh::MeshBase* vol_mesh,
                         const std::string& restart_read_dirname = "",
                         unsigned int restart_restore_number = 0);

    /*!
     * \brief Default deconstructor.
     */
    ~WallSitesMeshMapping();

    /*!
     * \brief Default constructor.
     */
    WallSitesMeshMapping() = delete;

    /*!
     * \brief Deleted copy constructor.
     */
    WallSitesMeshMapping(const WallSitesMeshMapping& from) = delete;

    /*!
     * \brief Deleted assignment operator.
     */
    WallSitesMeshMapping& operator=(const WallSitesMeshMapping& that) = delete;

    void initializeFEData(SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM>> hierarchy);

    void buildBoundaryMesh() override;

    /*!
     * \brief Spread the number of used wall sites. The provided patch index must have enough ghost cell width to
     * perform the requested spreading operation.
     */
    void spreadWallSites(int w_idx = IBTK::invalid_index);

    /*!
     * \brief Return the patch index containing the wall sites. The Eulerian description of the wall sites are updated
     * automatically when the boundary mesh moves.
     */
    inline int getWallSitesPatchIndex()
    {
        return d_w_idx;
    }

protected:
    std::string d_W_sys_name;

    // We also maintain an Eulerian description of the wall sites
    SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double>> d_w_var;
    int d_w_idx = IBTK::invalid_index;

    std::string d_kernel_fcn = "BSPLINE_6";

private:
    void commonConstructor(SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> input_db);
};

} // namespace clot
#endif
