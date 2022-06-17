#ifndef included_clot_BoundaryMeshMapping
#define included_clot_BoundaryMeshMapping
#include <ADS/GeneralBoundaryMeshMapping.h>

#include <ibtk/FEDataManager.h>

#include <libmesh/boundary_mesh.h>
#include <libmesh/mesh.h>

namespace clot
{
/*!
 * BoundaryMeshMapping is a class that generalizes the notion of a boundary mesh. It maintains an EquationSystems
 * object on the mesh, handles restarts of Lagrangian data, and maintains a FEDataManager for the object.
 * Implementations for this object should define how the object moves. A default implementation of no motion is
 * included.
 */
class BoundaryMeshMapping : public ADS::GeneralBoundaryMeshMapping
{
public:
    /*!
     * \brief Constructor that takes in a boundary mesh. Note that BoundaryMeshMapping assumes ownership of the
     * mesh.
     */
    BoundaryMeshMapping(std::string object_name,
                        SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> input_db,
                        libMesh::MeshBase* vol_mesh,
                        IBTK::FEDataManager* fe_data_manager,
                        const std::string& restart_read_dirname = "",
                        unsigned int restart_restore_number = 0);

    /*!
     * \brief Default deconstructor.
     */
    virtual ~BoundaryMeshMapping() = default;

    /*!
     * \brief Default constructor.
     */
    BoundaryMeshMapping() = delete;

    /*!
     * \brief Deleted copy constructor.
     */
    BoundaryMeshMapping(const BoundaryMeshMapping& from) = delete;

    /*!
     * \brief Deleted assignment operator.
     */
    BoundaryMeshMapping& operator=(const BoundaryMeshMapping& that) = delete;

    void buildBoundaryMesh() override;

    void updateBoundaryLocation(double time, unsigned int part, bool end_of_timestep = false) override;

protected:
    IBTK::FEDataManager* d_base_data_manager;

    // We need a map between original mesh nodes and bdry mesh nodes
    std::map<libMesh::Node*, libMesh::Node*> d_bdry_base_node_map;
private:
    void commonConstructor(SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> input_db);
};

} // namespace clot
#endif
