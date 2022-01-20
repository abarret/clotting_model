#ifndef included_IBAMR_BoundaryMeshMapping
#define included_IBAMR_BoundaryMeshMapping
#include "ibtk/FEDataManager.h"

#include "libmesh/boundary_mesh.h"
#include "libmesh/mesh.h"

namespace IBAMR
{
/*!
 * BoundaryMeshMapping is a class that generalizes the notion of a boundary mesh. It maintains an EquationSystems
 * object on the mesh, handles restarts of Lagrangian data, and maintains a FEMeshPartitioner for the object.
 * Implementations for this object should define how the object moves. A default implementation of no motion is
 * included.
 */
class BoundaryMeshMapping
{
public:
    /*!
     * \brief Constructor that takes in a vector of boundary meshes. Note that BoundaryMeshMapping assumes
     * ownership of the meshes.
     */
    BoundaryMeshMapping(std::string object_name,
                        SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> input_db,
                        const std::vector<libMesh::MeshBase*>& vol_meshes,
                        const std::vector<IBTK::FEDataManager*>& fe_data_managers,
                        const std::string& restart_read_dirname = "",
                        unsigned int restart_restore_number = 0);

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
     * \brief Constructor that leaves object in undefined state.
     */
    BoundaryMeshMapping(std::string object_name);

    /*!
     * \brief Default constructor.
     */
    BoundaryMeshMapping() = delete;

    /*!
     * \brief Default deconstructor.
     */
    virtual ~BoundaryMeshMapping() = default;

    /*!
     * \brief Deleted copy constructor.
     */
    BoundaryMeshMapping(const BoundaryMeshMapping& from) = delete;

    /*!
     * \brief Deleted assignment operator.
     */
    BoundaryMeshMapping& operator=(const BoundaryMeshMapping& that) = delete;

    /*!
     * \name Get the FEMeshPartitioners.
     */
    //\{

    virtual inline IBTK::FEDataManager* getFEDataManager(unsigned int part = 0)
    {
        return d_bdry_data_managers[part];
    }

    virtual inline const std::vector<IBTK::FEDataManager*>& getFEDataManagers()
    {
        return d_bdry_data_managers;
    }

    virtual inline std::vector<IBTK::FEDataManager*> getFEDataManagers(std::set<unsigned int> mesh_nums)
    {
        std::vector<IBTK::FEDataManager*> data_managers;
        for (const auto& mesh_num : mesh_nums) data_managers.push_back(d_bdry_data_managers[mesh_num]);
        return data_managers;
    }

    //\}

    /*!
     * Update the location of the boundary mesh. An optional argument is provided if the location of the structure is
     * needed at the end of the timestep. By default, this function loops over parts and calls the part specific
     * function.
     *
     * \note This assumes the correct velocity field has been set up via the findVelocity function.
     */
    virtual void
    updateBoundaryLocation(double t_start, double t_end, bool end_of_timestep = false, bool initial_time = false);
    /*!
     * Update the location of the boundary mesh for a specific part. An optional argument is provided if the location of
     * the structure is needed at the end of the timestep. By default this function does nothing.
     *
     * \note This assumes the correct velocity field has been set up via the findVelocity function.
     */
    virtual void updateBoundaryLocation(double t_start, double t_end, unsigned int part, bool end_of_timestep = false, bool initial_time = false);

    /*!
     * \brief Initialize the equations systems. Note all systems should be registered with the Equation systems prior to
     * this call. This function also initialized the location of the boundary mesh.
     */
    virtual void initializeEquationSystems();

    /*!
     * \brief Set the initial conditions of the structure. The position is set to the initial position of the mesh.
     */
    virtual void setInitialConditions();

    /*!
     * \brief Find the velocity of the wall points.
     */
    virtual void findVelocity(double time, int xi_idx, int sigma_idx);

    /*!
     * \brief Find the force exerted on the fluid by the wall. The value of Xi has been set through the findVelocity
     * function.
     *
     * \note This assumes the correct velocity field has been set up via the findVelocity function.
     */
    virtual void findForce(double time, int f_idx);

    /*!
     * \brief Write data to a restart file.
     */
    virtual void writeFEDataToRestartFile(const std::string& restart_dump_dirname, unsigned int time_step_number);

protected:
    std::string d_object_name;
    SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM>> d_hierarchy;
    std::vector<libMesh::MeshBase*> d_vol_meshes;
    std::vector<IBTK::FEDataManager*> d_vol_data_managers;

    std::vector<std::shared_ptr<IBTK::FEData>> d_fe_data;
    std::vector<IBTK::FEDataManager*> d_bdry_data_managers;
    std::vector<std::unique_ptr<libMesh::BoundaryMesh>> d_bdry_meshes;
    std::vector<std::unique_ptr<libMesh::EquationSystems>> d_bdry_eq_sys_vec;

    std::string d_coords_sys_name = "COORDINATES_SYSTEM";
    std::string d_disp_sys_name = "DISPLACEMENT_SYSTEM";
    std::string d_velocity_sys_name = "VELOCITY_SYSTEM";
    std::string d_xi_sys_name = "XI_SYSTEM";
    std::string d_sigma_sys_name = "SIGMA_SYSTEM";
    std::string d_force_sys_name = "FORCING_SYSTEM";

    // Restart data
    std::string d_libmesh_restart_file_extension = "xdr";

    // Force balance data
    double d_Sw = std::numeric_limits<double>::quiet_NaN();
    double d_abs_thresh = 1.0e-8;

private:
    void commonConstructor(SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> input_db,
                           std::string restart_read_dirname,
                           unsigned int restart_restore_number);
};

} // namespace IBAMR
#endif
