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
                         IBTK::FEDataManager* fe_data_manager,
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

    void initializeFEData() override;

    void buildBoundaryMesh() override;

    void updateBoundaryLocation(double time, unsigned int part, bool end_of_timestep = false) override;

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
    IBTK::FEDataManager* d_base_data_manager;

    // We need a map between original mesh nodes and bdry mesh nodes
    std::map<libMesh::Node*, libMesh::Node*> d_bdry_base_node_map;

    // We also maintain an Eulerian description of the wall sites
    SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double>> d_w_var;
    int d_w_idx = IBTK::invalid_index;

private:
    void commonConstructor(SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> input_db);
};

} // namespace clot
#endif
