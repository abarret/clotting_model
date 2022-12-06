// ---------------------------------------------------------------------
//
// Copyright (c) 2019 - 2020 by the IBAMR developers
// All rights reserved.
//
// This file is part of IBAMR.
//
// IBAMR is free software and is distributed under the 3-clause BSD
// license. The full text of the license can be found in the file
// COPYRIGHT at the top level directory of IBAMR.
//
// ---------------------------------------------------------------------

// Config files
#include <ibamr/config.h>

#include <clot/BondBoundSource.h>
#include <clot/BoundExtraStressForcing.h>
#include <clot/BoundPlateletSource.h>
#include <clot/BoundVelocitySource.h>
#include <clot/BoundaryMeshMapping.h>
#include <clot/CellToFaceFcn.h>
#include <clot/ClotParameters.h>
#include <clot/CohesionStressBoundRHS.h>
#include <clot/DragForce.h>
#include <clot/app_namespaces.h>

#include <ADS/CutCellVolumeMeshMapping.h>
#include <ADS/LSCutCellLaplaceOperator.h>
#include <ADS/LSFromMesh.h>
#include <ADS/SBAdvDiffIntegrator.h>
#include <ADS/SBIntegrator.h>

#include <ibamr/AdvDiffSemiImplicitHierarchyIntegrator.h>
#include <ibamr/CFINSForcing.h>
#include <ibamr/IBExplicitHierarchyIntegrator.h>
#include <ibamr/IIMethod.h>
#include <ibamr/INSCollocatedHierarchyIntegrator.h>
#include <ibamr/INSStaggeredHierarchyIntegrator.h>

#include <ibtk/AppInitializer.h>
#include <ibtk/CartGridFunctionSet.h>
#include <ibtk/IBTKInit.h>
#include <ibtk/IBTK_MPI.h>
#include <ibtk/muParserCartGridFunction.h>
#include <ibtk/muParserRobinBcCoefs.h>

#include <libmesh/boundary_info.h>
#include <libmesh/boundary_mesh.h>
#include <libmesh/elem.h>
#include <libmesh/equation_systems.h>
#include <libmesh/exodusII_io.h>
#include <libmesh/mesh.h>
#include <libmesh/mesh_generation.h>
#include <libmesh/mesh_modification.h>
#include <libmesh/mesh_refinement.h>
#include <libmesh/mesh_tools.h>
#include <libmesh/mesh_triangle_interface.h>

#include <petscsys.h>

#include <BergerRigoutsos.h>
#include <CartesianGridGeometry.h>
#include <LoadBalancer.h>
#include <SAMRAI_config.h>
#include <StandardTagAndInitialize.h>

#include <fstream>
#include <iostream>

// Local includes
#include "WallSitesMeshMapping.h"

// Function prototypes
void output_data(Pointer<PatchHierarchy<NDIM>> patch_hierarchy,
                 Pointer<INSHierarchyIntegrator> ins_integrator,
                 const int iteration_num,
                 const double loop_time,
                 const string& data_dump_dirname);

void
bdry_fcn(const IBTK::VectorNd& /*x*/, double& ls_val)
{
    ls_val = -5.0;
}

struct BetaFcn
{
public:
    BetaFcn(const double R0, const double beta_0, const double beta_1) : d_R0(R0), d_beta_0(beta_0), d_beta_1(beta_1)
    {
        // intentionally blank
        return;
    }
    double beta(const double eps) const
    {
        if (eps > d_R0)
        {
            return d_beta_0 * std::exp(d_beta_1 * (eps - d_R0));
            //            return d_beta_0 + d_beta_1 * (eps - d_R0);
        }
        else
        {
            return d_beta_0;
        }
    }

private:
    double d_R0 = std::numeric_limits<double>::quiet_NaN();
    double d_beta_0 = std::numeric_limits<double>::quiet_NaN();
    double d_beta_1 = std::numeric_limits<double>::quiet_NaN();
};

double
beta_wrapper(const double eps, void* ctx)
{
    auto beta_fcn = static_cast<BetaFcn*>(ctx);
    return beta_fcn->beta(eps);
}

// for now this is an amorphous blob
struct ClampVars
{
    Pointer<PatchHierarchy<NDIM>> patch_hierarchy;
    Pointer<VariableContext> context;
    // Variables to clamp go here
    // if more are needed, make this an array
    Pointer<hier::Variable<NDIM>> phi_b_var;
    Pointer<hier::Variable<NDIM>> phi_a_var;
    Pointer<hier::Variable<NDIM>> bond_var;
};

double A = 1.0;
double f = 1.0;
double t_start = 0.0;
libMesh::Point
center(const double t)
{
    // Return the displacement of the center as a function of time.
    libMesh::Point d;
    if (t >= t_start) d(0) = A * std::sin(2.0 * M_PI * (t - t_start) / f);
    return d;
}

VectorNd
center_vel(const double t)
{
    // Return the derivative of the displacement
    VectorNd d(VectorNd::Zero());
    if (t >= t_start) d[0] = 2.0 * M_PI * A / f * std::cos(2.0 * M_PI * (t - t_start) / f);
    return d;
}

namespace ModelData
{
static double kappa = 1.0e6;
static double eta = 0.0;
static double dx = -1.0;
static bool ERROR_ON_MOVE = false;

void
tether_force_function(VectorValue<double>& F,
                      const VectorValue<double>& /*n*/,
                      const VectorValue<double>& /*N*/,
                      const TensorValue<double>& /*FF*/,
                      const libMesh::Point& x,
                      const libMesh::Point& X,
                      Elem* const /*elem*/,
                      const unsigned short /*side*/,
                      const std::vector<const std::vector<double>*>& var_data,
                      const std::vector<const std::vector<VectorValue<double>>*>& /*grad_var_data*/,
                      double time,
                      void* /*ctx*/)
{
    const std::vector<double>& U = *var_data[0];
    const libMesh::Point disp = center(time) + X;
    for (unsigned int d = 0; d < NDIM; ++d) F(d) = kappa * (disp(d) - x(d)) - eta * (U[d] - center_vel(time)[d]);
    VectorValue<double> d = disp - x;
    if (ERROR_ON_MOVE && d.norm() > 0.25 * dx) TBOX_ERROR("Structure has moved too much.\n");
}

} // namespace ModelData
using namespace ModelData;

double
wall_sites_ode(const double w,
               const std::vector<double>& fl_vals,
               const std::vector<double>& /*sf_vals*/,
               double /*time*/,
               void* ctx)
{
    auto params = static_cast<BoundClotParams*>(ctx);
    return -params->Kaw * params->nw_max * w * params->nb_max * std::max(fl_vals[0], 0.0);
}

static double init_wall_sites = 1.0;
double
wall_sites_init(const VectorNd& /*X*/, const Node* const /*node*/)
{
    return init_wall_sites;
}

static std::ofstream force_stream, tether_stream;
void output_net_force(EquationSystems* eq_sys, const double loop_time);

void compute_plot_quantities(const int drag_idx,
                             const int div_sig_idx,
                             const int beta_idx,
                             Pointer<PatchHierarchy<NDIM>> hierarchy,
                             Pointer<INSHierarchyIntegrator> ins_integrator,
                             Pointer<CFINSForcing> cf_forcing,
                             Pointer<CellVariable<NDIM, double>> bond_var,
                             Pointer<AdvDiffHierarchyIntegrator> bond_integrator,
                             const int ub_idx,
                             const int phib_idx,
                             const BoundClotParams& clot_params,
                             const BetaFcn& beta_fcn);

/*******************************************************************************
 * For each run, the input filename and restart information (if needed) must   *
 * be given on the command line.  For non-restarted case, command line is:     *
 *                                                                             *
 *    executable <input file name>                                             *
 *                                                                             *
 * For restarted run, command line is:                                         *
 *                                                                             *
 *    executable <input file name> <restart directory> <restart number>        *
 *                                                                             *
 *******************************************************************************/
int
main(int argc, char* argv[])
{
    // Initialize IBAMR and libraries. Deinitialization is handled by this object as well.
    IBTKInit ibtk_init(argc, argv, MPI_COMM_WORLD);
    const LibMeshInit& init = ibtk_init.getLibMeshInit();

    { // cleanup dynamically allocated objects prior to shutdown

        // Parse command line options, set some standard options from the input
        // file, initialize the restart database (if this is a restarted run),
        // and enable file logging.
        Pointer<AppInitializer> app_initializer = new AppInitializer(argc, argv, "INS.log");
        Pointer<Database> input_db = app_initializer->getInputDatabase();

        // Get various standard options set in the input file.
        const bool dump_viz_data = app_initializer->dumpVizData();
        const bool uses_visit = dump_viz_data && app_initializer->getVisItDataWriter();
#ifdef LIBMESH_HAVE_EXODUS_API
        const bool uses_exodus = dump_viz_data && !app_initializer->getExodusIIFilename().empty();
#else
        const bool uses_exodus = false;
        if (!app_initializer->getExodusIIFilename().empty())
        {
            plog << "WARNING: libMesh was compiled without Exodus support, so no "
                 << "Exodus output will be written in this program.\n";
        }
#endif
        const std::string base_mesh_filename = app_initializer->getExodusIIFilename("mesh_");
        const std::string bdry_mesh_filename = app_initializer->getExodusIIFilename("bdry_");
        const std::string wall_mesh_filename = app_initializer->getExodusIIFilename("wall_");

        const bool dump_restart_data = app_initializer->dumpRestartData();
        const int restart_dump_interval = app_initializer->getRestartDumpInterval();
        const string restart_dump_dirname = app_initializer->getRestartDumpDirectory();
        const string restart_read_dirname = app_initializer->getRestartReadDirectory();
        const int restart_restore_num = app_initializer->getRestartRestoreNumber();

        const bool dump_postproc_data = app_initializer->dumpPostProcessingData();
        const int postproc_data_dump_interval = app_initializer->getPostProcessingDataDumpInterval();
        const string postproc_data_dump_dirname = app_initializer->getPostProcessingDataDumpDirectory();

        const bool dump_timer_data = app_initializer->dumpTimerData();
        const int timer_dump_interval = app_initializer->getTimerDumpInterval();

        // Create finite element mesh
        Mesh solid_mesh(init.comm(), NDIM);
        dx = input_db->getDouble("DX");
        ERROR_ON_MOVE = input_db->getBool("ERROR_ON_MOVE");
        // mfac is the mesh factor, how many Lagrangian nodes per Eulerian grid cell?
        const double mfac = input_db->getDouble("MFAC");
        const double ds = mfac * dx;
        std::string elem_type = input_db->getString("ELEM_TYPE");
        std::string elem_order = input_db->getString("ELEM_ORDER");
        const double R = input_db->getDouble("RADIUS");
        if (NDIM == 2 && (elem_type == "TRI3" || elem_type == "TRI6"))
        {
#ifdef LIBMESH_HAVE_TRIANGLE
            const int num_circum_nodes = ceil(2.0 * M_PI * R / ds);
            for (int k = 0; k < num_circum_nodes; ++k)
            {
                const double theta = 2.0 * M_PI * static_cast<double>(k) / static_cast<double>(num_circum_nodes);
                solid_mesh.add_point(libMesh::Point(R * cos(theta), R * sin(theta)));
            }
            TriangleInterface triangle(solid_mesh);
            triangle.triangulation_type() = TriangleInterface::GENERATE_CONVEX_HULL;
            triangle.elem_type() = Utility::string_to_enum<ElemType>(elem_type);
            triangle.desired_area() = 1.5 * sqrt(3.0) / 4.0 * ds * ds;
            triangle.insert_extra_points() = true;
            triangle.smooth_after_generating() = true;
            triangle.triangulate();
#else
            TBOX_ERROR("ERROR: libMesh appears to have been configured without support for Triangle,\n"
                       << "       but Triangle is required for TRI3 or TRI6 elements.\n");
#endif
        }
        else
        {
            // NOTE: number of segments along boundary is 4*2^r.
            const double num_circum_segments = 2.0 * M_PI * R / ds;
            const int r = log2(0.25 * num_circum_segments);
            MeshTools::Generation::build_sphere(solid_mesh, R, r, Utility::string_to_enum<ElemType>(elem_type));
        }

        // Ensure nodes on the surface are on the analytic boundary.
        MeshBase::element_iterator el_end = solid_mesh.elements_begin();
        for (auto el = solid_mesh.elements_begin(); el != el_end; ++el)
        {
            Elem* const elem = *el;
            for (unsigned int side = 0; side < elem->n_sides(); ++side)
            {
                const bool at_mesh_bdry = !elem->neighbor_ptr(side);
                if (!at_mesh_bdry) continue;
                for (unsigned int k = 0; k < elem->n_nodes(); ++k)
                {
                    if (!elem->is_node_on_side(k, side)) continue;
                    Node& n = *elem->node_ptr(k);
                    n = R * n.unit();
                }
            }
        }
        if (elem_order == "SECOND")
            solid_mesh.all_second_order(true);
        else
            solid_mesh.all_first_order();
        MeshTools::Modification::translate(solid_mesh, ds * 0.1, ds * 0.1);
        solid_mesh.prepare_for_use();
        solid_mesh.print_info();

        BoundaryMesh bdry_mesh(solid_mesh.comm(), solid_mesh.mesh_dimension() - 1);
        BoundaryInfo& bdry_info = solid_mesh.get_boundary_info();
        bdry_info.sync(bdry_mesh);
        bdry_mesh.prepare_for_use();

        MeshBase& mesh = bdry_mesh;

        // Get structure parameters
        kappa = input_db->getDouble("KAPPA");
        eta = input_db->getDouble("ETA");

        // Create major algorithm and data objects that comprise the
        // application.  These objects are configured from the input database
        // and, if this is a restarted run, from the restart database.
        Pointer<INSHierarchyIntegrator> ins_integrator = new INSStaggeredHierarchyIntegrator(
            "INSStaggeredHierarchyIntegrator",
            app_initializer->getComponentDatabase("INSStaggeredHierarchyIntegrator"));
        Pointer<SBAdvDiffIntegrator> sb_adv_diff_integrator =
            new SBAdvDiffIntegrator("SBAdvDiffIntegrator", app_initializer->getComponentDatabase("AdvDiffIntegrator"));
        Pointer<AdvDiffSemiImplicitHierarchyIntegrator> adv_diff_integrator;
        adv_diff_integrator = new AdvDiffSemiImplicitHierarchyIntegrator(
            "AdvDiffSemiImplicitHierarchyIntegrator", app_initializer->getComponentDatabase("AdvDiffIntegrator"));
        Pointer<IIMethod> ib_method_ops =
            new IIMethod("IIMethod",
                         app_initializer->getComponentDatabase("IIMethod"),
                         &mesh,
                         app_initializer->getComponentDatabase("GriddingAlgorithm")->getInteger("max_levels"),
                         true,
                         restart_read_dirname,
                         restart_restore_num);
        Pointer<IBHierarchyIntegrator> time_integrator =
            new IBExplicitHierarchyIntegrator("IBHierarchyIntegrator",
                                              app_initializer->getComponentDatabase("IBHierarchyIntegrator"),
                                              ib_method_ops,
                                              ins_integrator);
        Pointer<CartesianGridGeometry<NDIM>> grid_geometry = new CartesianGridGeometry<NDIM>(
            "CartesianGeometry", app_initializer->getComponentDatabase("CartesianGeometry"));
        Pointer<PatchHierarchy<NDIM>> patch_hierarchy = new PatchHierarchy<NDIM>("PatchHierarchy", grid_geometry);
        Pointer<StandardTagAndInitialize<NDIM>> error_detector =
            new StandardTagAndInitialize<NDIM>("StandardTagAndInitialize",
                                               time_integrator,
                                               app_initializer->getComponentDatabase("StandardTagAndInitialize"));
        Pointer<BergerRigoutsos<NDIM>> box_generator = new BergerRigoutsos<NDIM>();
        Pointer<LoadBalancer<NDIM>> load_balancer =
            new LoadBalancer<NDIM>("LoadBalancer", app_initializer->getComponentDatabase("LoadBalancer"));
        Pointer<GriddingAlgorithm<NDIM>> gridding_algorithm =
            new GriddingAlgorithm<NDIM>("GriddingAlgorithm",
                                        app_initializer->getComponentDatabase("GriddingAlgorithm"),
                                        error_detector,
                                        box_generator,
                                        load_balancer);
        ins_integrator->registerAdvDiffHierarchyIntegrator(adv_diff_integrator);
        ins_integrator->registerAdvDiffHierarchyIntegrator(sb_adv_diff_integrator);
        // Create initial condition specification objects.
        Pointer<CartGridFunction> u_init = new muParserCartGridFunction(
            "u_init", app_initializer->getComponentDatabase("VelocityInitialConditions"), grid_geometry);
        ins_integrator->registerVelocityInitialConditions(u_init);
        Pointer<CartGridFunction> p_init = new muParserCartGridFunction(
            "p_init", app_initializer->getComponentDatabase("PressureInitialConditions"), grid_geometry);
        ins_integrator->registerPressureInitialConditions(p_init);
        Pointer<CartGridFunction> phib_init =
            new muParserCartGridFunction("phibInit", app_initializer->getComponentDatabase("phibInit"), grid_geometry);
        Pointer<CartGridFunction> bond_init =
            new muParserCartGridFunction("bondInit", app_initializer->getComponentDatabase("bondInit"), grid_geometry);

        // Configure the IBFE solver
        ib_method_ops->initializeFEEquationSystems();
        std::vector<int> vars(NDIM);
        for (unsigned int d = 0; d < NDIM; ++d) vars[d] = d;
        std::vector<SystemData> sys_data = { SystemData(IIMethod::VELOCITY_SYSTEM_NAME, vars) };
        IIMethod::LagSurfaceForceFcnData surface_fcn_data(tether_force_function, sys_data);
        ib_method_ops->registerLagSurfaceForceFunction(surface_fcn_data);

        // Create boundary condition specification objects (when necessary).
        const IntVector<NDIM>& periodic_shift = grid_geometry->getPeriodicShift();
        vector<RobinBcCoefStrategy<NDIM>*> u_bc_coefs(NDIM);
        if (periodic_shift.min() > 0)
        {
            for (unsigned int d = 0; d < NDIM; ++d)
            {
                u_bc_coefs[d] = NULL;
            }
        }
        else
        {
            for (unsigned int d = 0; d < NDIM; ++d)
            {
                ostringstream bc_coefs_name_stream;
                bc_coefs_name_stream << "u_bc_coefs_" << d;
                const string bc_coefs_name = bc_coefs_name_stream.str();

                ostringstream bc_coefs_db_name_stream;
                bc_coefs_db_name_stream << "VelocityBcCoefs_" << d;
                const string bc_coefs_db_name = bc_coefs_db_name_stream.str();

                u_bc_coefs[d] = new muParserRobinBcCoefs(
                    bc_coefs_name, app_initializer->getComponentDatabase(bc_coefs_db_name), grid_geometry);
            }
            ins_integrator->registerPhysicalBoundaryConditions(u_bc_coefs);
        }

        Pointer<NodeVariable<NDIM, double>> ls_var = new NodeVariable<NDIM, double>("LS");
        sb_adv_diff_integrator->registerLevelSetVariable(ls_var);

        // Set up visualization plot file writers.
        Pointer<VisItDataWriter<NDIM>> visit_data_writer = app_initializer->getVisItDataWriter();
        if (uses_visit)
        {
            time_integrator->registerVisItDataWriter(visit_data_writer);
        }

        // Structure motion parameters
        A = input_db->getDouble("Amplitude");
        f = input_db->getDouble("Frequency");
        t_start = input_db->getDouble("T_START");

        // Create advected quantities
        Pointer<CellVariable<NDIM, double>> phi_b_var = new CellVariable<NDIM, double>("phi_b");
        Pointer<CellVariable<NDIM, double>> phi_a_var = new CellVariable<NDIM, double>("phi_a");
        Pointer<CellVariable<NDIM, double>> bond_var = new CellVariable<NDIM, double>("bond");
        Pointer<CellVariable<NDIM, double>> ub_var = new CellVariable<NDIM, double>("Bound Velocity", NDIM);

        // Pull out the database everything uses
        BoundClotParams clot_params(app_initializer->getComponentDatabase("ClotParams"));
        pout << "nw = " << clot_params.nw << "\n";
        Pointer<CellToFaceFcn> cell_to_face_fcn = new CellToFaceFcn("BoundAdvFcn", ub_var, adv_diff_integrator);

        // Set up Cohesion stress tensor
        Pointer<CFINSForcing> cohesionStressForcing =
            new CFINSForcing("CohesionStressForcing",
                             app_initializer->getComponentDatabase("CohesionStress"),
                             static_cast<Pointer<CartGridFunction>>(cell_to_face_fcn),
                             grid_geometry,
                             adv_diff_integrator,
                             visit_data_writer);
        ins_integrator->registerBodyForceFunction(cohesionStressForcing);
        Pointer<CohesionStressBoundRHS> cohesion_relax = new CohesionStressBoundRHS("CohesionRHS", clot_params);
        cohesion_relax->setBoundPlateletData(phi_b_var, sb_adv_diff_integrator);
        cohesion_relax->setActivatedPlateletData(phi_a_var, sb_adv_diff_integrator);
        cohesion_relax->setBondData(bond_var, sb_adv_diff_integrator);
        cohesion_relax->setKernel(BSPLINE_3);
        cohesionStressForcing->registerRelaxationOperator(cohesion_relax);
        Pointer<CellVariable<NDIM, double>> sig_var = cohesionStressForcing->getVariable();

        // Set up bound velocity function and register it with the sb integrator
        Pointer<FaceVariable<NDIM, double>> ub_adv_var = new FaceVariable<NDIM, double>("bound_velocity");
        sb_adv_diff_integrator->registerAdvectionVelocity(ub_adv_var);
        sb_adv_diff_integrator->setAdvectionVelocityFunction(ub_adv_var, cell_to_face_fcn);
        adv_diff_integrator->registerAdvectionVelocity(ub_adv_var);
        adv_diff_integrator->setAdvectionVelocityFunction(ub_adv_var, cell_to_face_fcn);

        // Create the bound velocity solver.
        Pointer<BoundVelocitySource> ub_src_fcn = new BoundVelocitySource("BoundVelocitySource", clot_params);
        ub_src_fcn->setBondData(bond_var, sb_adv_diff_integrator);
        ub_src_fcn->setBoundPlateletData(phi_b_var, sb_adv_diff_integrator);
        ub_src_fcn->setSigmaData(sig_var, adv_diff_integrator);
        ub_src_fcn->setVelocityData(ins_integrator->getVelocityVariable(), ins_integrator);
        ub_src_fcn->setBoundVelocityData(ub_var, adv_diff_integrator);
        adv_diff_integrator->registerTransportedQuantity(ub_var);
        adv_diff_integrator->setAdvectionVelocity(ub_var, ub_adv_var);
        std::vector<RobinBcCoefStrategy<NDIM>*> ub_bcs(NDIM, nullptr);
        if (grid_geometry->getPeriodicShift().min() == 0)
        {
            for (int d = 0; d < NDIM; ++d)
            {
                ub_bcs[d] = new muParserRobinBcCoefs(
                    "BoundVelocityBc_" + std::to_string(d),
                    app_initializer->getComponentDatabase("BoundVelocityBc_" + std::to_string(d)),
                    grid_geometry);
            }
        }
        adv_diff_integrator->setPhysicalBcCoefs(ub_var, ub_bcs);
        adv_diff_integrator->setDiffusionCoefficient(ub_var, input_db->getDouble("MU"));
        Pointer<CellVariable<NDIM, double>> ub_src_var = new CellVariable<NDIM, double>("Ub src", NDIM);
        adv_diff_integrator->registerSourceTerm(ub_src_var);
        adv_diff_integrator->setSourceTermFunction(ub_src_var, ub_src_fcn);
        adv_diff_integrator->setSourceTerm(ub_var, ub_src_var);

        // Drag force
        Pointer<DragForce> drag_force = new DragForce("DragForce", clot_params);
        drag_force->setBoundPlateletData(phi_b_var, sb_adv_diff_integrator);
        drag_force->setVelocityData(ins_integrator->getVelocityVariable(), ins_integrator);
        drag_force->setBoundVelocityData(ub_var, adv_diff_integrator);
        ins_integrator->registerBodyForceFunction(drag_force);

        // Create beta function
        double beta_0 = input_db->getDouble("BETA_0");
        double beta_1 = input_db->getDouble("BETA_1");
        double r0 = input_db->getDouble("R0");
        BetaFcn betaFcn(r0, beta_0, beta_1);
        cohesion_relax->registerBetaFcn(beta_wrapper, static_cast<void*>(&betaFcn));

        // Activated platelets
        sb_adv_diff_integrator->registerTransportedQuantity(phi_a_var);
        sb_adv_diff_integrator->restrictToLevelSet(phi_a_var, ls_var);
        Pointer<CellVariable<NDIM, double>> phi_a_src_var = new CellVariable<NDIM, double>("phi_a_src");
        sb_adv_diff_integrator->setAdvectionVelocity(phi_a_var, ins_integrator->getAdvectionVelocityVariable());
        sb_adv_diff_integrator->setDiffusionCoefficient(phi_a_var, clot_params.phi_a_diff_coef);
        Pointer<BoundPlateletSource> phi_a_src_fcn = new BoundPlateletSource("UnactivatedPlatelets", clot_params);
        phi_a_src_fcn->setSign(true);
        phi_a_src_fcn->setKernel(BSPLINE_3);
        phi_a_src_fcn->registerBetaFcn(beta_wrapper, static_cast<void*>(&betaFcn));
        phi_a_src_fcn->setActivatedPlateletData(phi_a_var, sb_adv_diff_integrator);
        phi_a_src_fcn->setBondsData(bond_var, sb_adv_diff_integrator);
        phi_a_src_fcn->setBoundPlateletData(phi_b_var, sb_adv_diff_integrator);
        phi_a_src_fcn->setStressData(sig_var, adv_diff_integrator);
        sb_adv_diff_integrator->registerSourceTerm(phi_a_src_var);
        sb_adv_diff_integrator->setSourceTermFunction(phi_a_src_var, phi_a_src_fcn);
        sb_adv_diff_integrator->setSourceTerm(phi_a_var, phi_a_src_var);
        Pointer<RobinBcCoefStrategy<NDIM>> phi_a_bcs;
        if (grid_geometry->getPeriodicShift().min() == 0)
            phi_a_bcs = new muParserRobinBcCoefs("UnactivatedPlateletsBcs",
                                                 app_initializer->getComponentDatabase("ActivatedPlateletBcs"),
                                                 grid_geometry);
        sb_adv_diff_integrator->setPhysicalBcCoef(phi_a_var, phi_a_bcs.getPointer());
        Pointer<LSCutCellLaplaceOperator> phi_a_rhs_oper = new LSCutCellLaplaceOperator(
                                              "PhiALSCutCellRHSOperator",
                                              app_initializer->getComponentDatabase("LSCutCellOperator"),
                                              false),
                                          phi_a_oper = new LSCutCellLaplaceOperator(
                                              "PhiALSCutCellOperator",
                                              app_initializer->getComponentDatabase("LSCutCellOperator"),
                                              false);
        Pointer<PETScKrylovPoissonSolver> phi_a_solver = new PETScKrylovPoissonSolver(
            "PhiAPoissonSolver", app_initializer->getComponentDatabase("PoissonSolver"), "poisson_solve_");
        phi_a_solver->setOperator(phi_a_oper);
        sb_adv_diff_integrator->setHelmholtzSolver(phi_a_var, phi_a_solver);
        sb_adv_diff_integrator->setHelmholtzRHSOperator(phi_a_var, phi_a_rhs_oper);

        // Bound platelets
        sb_adv_diff_integrator->registerTransportedQuantity(phi_b_var);
        sb_adv_diff_integrator->restrictToLevelSet(phi_b_var, ls_var);
        sb_adv_diff_integrator->setInitialConditions(phi_b_var, phib_init);
        Pointer<CellVariable<NDIM, double>> phi_b_src_var = new CellVariable<NDIM, double>("phi_b_src");
        sb_adv_diff_integrator->setAdvectionVelocity(phi_b_var, ub_adv_var);
        sb_adv_diff_integrator->setDiffusionCoefficient(phi_b_var, clot_params.phi_b_diff_coef);
        Pointer<BoundPlateletSource> phi_b_src_fcn = new BoundPlateletSource("BoundPlatelets", clot_params);
        phi_b_src_fcn->setSign(false);
        phi_b_src_fcn->setKernel(BSPLINE_3);
        phi_b_src_fcn->registerBetaFcn(beta_wrapper, static_cast<void*>(&betaFcn));
        phi_b_src_fcn->setActivatedPlateletData(phi_a_var, sb_adv_diff_integrator);
        phi_b_src_fcn->setBondsData(bond_var, sb_adv_diff_integrator);
        phi_b_src_fcn->setBoundPlateletData(phi_b_var, sb_adv_diff_integrator);
        phi_b_src_fcn->setStressData(sig_var, adv_diff_integrator);
        sb_adv_diff_integrator->registerSourceTerm(phi_b_src_var);
        sb_adv_diff_integrator->setSourceTermFunction(phi_b_src_var, phi_b_src_fcn);
        sb_adv_diff_integrator->setSourceTerm(phi_b_var, phi_b_src_var);
        Pointer<RobinBcCoefStrategy<NDIM>> phi_b_bcs;
        if (grid_geometry->getPeriodicShift().min() == 0)
            phi_b_bcs = new muParserRobinBcCoefs(
                "ActivatedPlateletsBcs", app_initializer->getComponentDatabase("BoundPlateletBcs"), grid_geometry);
        sb_adv_diff_integrator->setPhysicalBcCoef(phi_b_var, phi_b_bcs.getPointer());
        Pointer<LSCutCellLaplaceOperator> phi_b_rhs_oper = new LSCutCellLaplaceOperator(
                                              "PhiBLSCutCellRHSOperator",
                                              app_initializer->getComponentDatabase("LSCutCellOperator"),
                                              false),
                                          phi_b_oper = new LSCutCellLaplaceOperator(
                                              "PhiBLSCutCellOperator",
                                              app_initializer->getComponentDatabase("LSCutCellOperator"),
                                              false);
        Pointer<PETScKrylovPoissonSolver> phi_b_solver = new PETScKrylovPoissonSolver(
            "PhiAPoissonSolver", app_initializer->getComponentDatabase("PoissonSolver"), "poisson_solve_");
        phi_b_solver->setOperator(phi_b_oper);
        sb_adv_diff_integrator->setHelmholtzSolver(phi_b_var, phi_b_solver);
        sb_adv_diff_integrator->setHelmholtzRHSOperator(phi_b_var, phi_b_rhs_oper);

        // Bonds
        sb_adv_diff_integrator->registerTransportedQuantity(bond_var);
        sb_adv_diff_integrator->setInitialConditions(bond_var, bond_init);
        sb_adv_diff_integrator->restrictToLevelSet(bond_var, ls_var);
        Pointer<CellVariable<NDIM, double>> bond_src_var = new CellVariable<NDIM, double>("bond_src");
        sb_adv_diff_integrator->setAdvectionVelocity(bond_var, ub_adv_var);
        sb_adv_diff_integrator->setDiffusionCoefficient(bond_var, 0.0);
        Pointer<BondBoundSource> bond_src_fcn = new BondBoundSource("BondSource", clot_params);
        bond_src_fcn->registerBetaFcn(beta_wrapper, static_cast<void*>(&betaFcn));
        bond_src_fcn->setActivatedPlateletData(phi_a_var, sb_adv_diff_integrator);
        bond_src_fcn->setBondData(bond_var, sb_adv_diff_integrator);
        bond_src_fcn->setBoundPlateletData(phi_b_var, sb_adv_diff_integrator);
        bond_src_fcn->setStressData(sig_var, adv_diff_integrator);
        bond_src_fcn->setKernel(BSPLINE_3);
        sb_adv_diff_integrator->registerSourceTerm(bond_src_var);
        sb_adv_diff_integrator->setSourceTermFunction(bond_src_var, bond_src_fcn);
        sb_adv_diff_integrator->setSourceTerm(bond_var, bond_src_var);
        Pointer<RobinBcCoefStrategy<NDIM>> bond_bcs;
        if (grid_geometry->getPeriodicShift().min() == 0)
            bond_bcs = new muParserRobinBcCoefs(
                "BondVariableBcs", app_initializer->getComponentDatabase("BondVariableBcs"), grid_geometry);
        sb_adv_diff_integrator->setPhysicalBcCoef(bond_var, bond_bcs.getPointer());
        Pointer<LSCutCellLaplaceOperator> bond_rhs_oper = new LSCutCellLaplaceOperator(
                                              "BondLSCutCellRHSOperator",
                                              app_initializer->getComponentDatabase("LSCutCellOperator"),
                                              false),
                                          bond_oper = new LSCutCellLaplaceOperator(
                                              "BondLSCutCellOperator",
                                              app_initializer->getComponentDatabase("LSCutCellOperator"),
                                              false);
        Pointer<PETScKrylovPoissonSolver> bond_solver = new PETScKrylovPoissonSolver(
            "BondPoissonSolver", app_initializer->getComponentDatabase("PoissonSolver"), "poisson_solve_");
        bond_solver->setOperator(bond_oper);
        sb_adv_diff_integrator->setHelmholtzSolver(bond_var, bond_solver);
        sb_adv_diff_integrator->setHelmholtzRHSOperator(bond_var, bond_rhs_oper);

        // Wall sites
        const std::string wall_sites_name = "WallSites";
        FEDataManager* fe_data_manager = ib_method_ops->getFEDataManager();
        auto wall_mesh_mapping =
            std::make_shared<WallSitesMeshMapping>("WallSitesMesh",
                                                   app_initializer->getComponentDatabase("BoundaryMesh"),
                                                   wall_sites_name,
                                                   &mesh,
                                                   fe_data_manager,
                                                   restart_read_dirname,
                                                   restart_restore_num);
        wall_mesh_mapping->initializeEquationSystems();
        auto sb_coupling_manager =
            std::make_shared<SBSurfaceFluidCouplingManager>("CouplingManager",
                                                            app_initializer->getComponentDatabase("CouplingManager"),
                                                            wall_mesh_mapping->getMeshPartitioner());
        sb_coupling_manager->registerFluidConcentration(phi_a_var);
        sb_coupling_manager->registerSurfaceConcentration(wall_sites_name);
        sb_coupling_manager->registerSurfaceReactionFunction(
            wall_sites_name, wall_sites_ode, static_cast<void*>(&clot_params));
        sb_coupling_manager->registerFluidConcentration(phi_a_var);
        sb_coupling_manager->registerFluidSurfaceDependence(wall_sites_name, phi_a_var);
        sb_coupling_manager->registerInitialConditions(wall_sites_name, wall_sites_init);
        sb_coupling_manager->initializeFEData();
        Pointer<SBIntegrator> sb_integrator = new SBIntegrator("SBIntegrator", sb_coupling_manager);
        sb_adv_diff_integrator->registerSBIntegrator(sb_integrator, ls_var);
        sb_adv_diff_integrator->registerLevelSetSBDataManager(ls_var, sb_coupling_manager);
        auto bdry_mesh_mapping =
            std::make_shared<BoundaryMeshMapping>("BoundaryMesh",
                                                  app_initializer->getComponentDatabase("BoundaryMesh"),
                                                  &mesh,
                                                  fe_data_manager,
                                                  restart_read_dirname,
                                                  restart_restore_num);
        bdry_mesh_mapping->initializeEquationSystems();
        Pointer<CutCellMeshMapping> cut_cell_mapping =
            new CutCellVolumeMeshMapping("CutCellMapping",
                                         app_initializer->getComponentDatabase("CutCellMapping"),
                                         bdry_mesh_mapping->getMeshPartitioner(0));
        Pointer<LSFromMesh> ls_fcn = new LSFromMesh("LSFcn", patch_hierarchy, cut_cell_mapping, false);
        ls_fcn->registerBdryFcn(bdry_fcn);
        sb_adv_diff_integrator->registerLevelSetVolFunction(ls_var, ls_fcn);
        sb_adv_diff_integrator->registerGeneralBoundaryMeshMapping(bdry_mesh_mapping);
        const int w_idx = wall_mesh_mapping->getWallSitesPatchIndex();
        std::pair<WallSitesMeshMapping*, BoundaryMeshMapping*> mesh_mappings =
            std::make_pair(wall_mesh_mapping.get(), bdry_mesh_mapping.get());
        auto update_bdry_mesh = [](const double current_time,
                                   const double new_time,
                                   bool /*skip_synchronize_new_state_data*/,
                                   int /*num_cycles*/,
                                   void* ctx) {
            plog << "Updating boundary location\n";
            auto mesh_mappings = static_cast<std::pair<BoundaryMeshMapping*, GeneralBoundaryMeshMapping*>*>(ctx);
            mesh_mappings->first->updateBoundaryLocation(current_time, 0, true);
            mesh_mappings->second->updateBoundaryLocation(current_time, 0, true);
        };
        auto regrid_callback = [](SAMRAI::tbox::Pointer<SAMRAI::hier::BasePatchHierarchy<NDIM>> /*hierarchy*/,
                                  const double /*time*/,
                                  const bool /*initial_time*/,
                                  void* ctx) {
            plog << "Spreading wall sites\n";
            auto wall_mesh_mapping = static_cast<WallSitesMeshMapping*>(ctx);
            wall_mesh_mapping->spreadWallSites();
        };
        time_integrator->registerPostprocessIntegrateHierarchyCallback(update_bdry_mesh,
                                                                       static_cast<void*>(&mesh_mappings));
        time_integrator->registerRegridHierarchyCallback(regrid_callback, static_cast<void*>(wall_mesh_mapping.get()));
        visit_data_writer->registerPlotQuantity("wall_sites", "SCALAR", w_idx);

        // define variable blob, which stores the patch hierarchy, variable context, and clamped variables.
        ClampVars clamp_vars = {
            patch_hierarchy,                         // patch hierarchy
            sb_adv_diff_integrator->getNewContext(), // variable context
            phi_b_var,                               // bound platelets var
            phi_a_var,                               // activated platelets var
            bond_var                                 // bonds var
        };
        // Clamp cell values to 0 in callback
        auto clamp_cell_values =
            [](const double /*current_time*/, const double /*new_time*/, int /*num_cycles*/, void* ctx) {
                pout << "Clamping variables\n";
                auto clamp_vars = static_cast<ClampVars*>(ctx);
                // iterate levels
                for (int ln = 0; ln <= clamp_vars->patch_hierarchy->getFinestLevelNumber(); ln++)
                {
                    Pointer<PatchLevel<NDIM>> level = clamp_vars->patch_hierarchy->getPatchLevel(ln);
                    // iterate patches
                    for (PatchLevel<NDIM>::Iterator p(level); p; p++)
                    {
                        Pointer<Patch<NDIM>> patch = level->getPatch(p());
                        // get patch data
                        Pointer<CellData<NDIM, double>> phi_b_data =
                            patch->getPatchData(clamp_vars->phi_b_var, clamp_vars->context);
                        Pointer<CellData<NDIM, double>> phi_a_data =
                            patch->getPatchData(clamp_vars->phi_a_var, clamp_vars->context);
                        Pointer<CellData<NDIM, double>> bond_data =
                            patch->getPatchData(clamp_vars->bond_var, clamp_vars->context);
                        // iterate cells
                        for (CellIterator<NDIM> ci(patch->getBox()); ci; ci++)
                        {
                            const CellIndex<NDIM>& idx = ci();
                            (*phi_b_data)(idx) = std::max((*phi_b_data)(idx), 0.0);
                            (*phi_a_data)(idx) = std::max((*phi_a_data)(idx), 0.0);
                            (*bond_data)(idx) = std::max((*bond_data)(idx), 0.0);
                        }
                    }
                }
            };
        // Register clamping callback
        sb_adv_diff_integrator->registerIntegrateHierarchyCallback(clamp_cell_values, static_cast<void*>(&clamp_vars));

        EquationSystems* eq_sys = fe_data_manager->getEquationSystems();
        EquationSystems* bdry_eq_sys = bdry_mesh_mapping->getMeshPartitioner()->getEquationSystems();
        EquationSystems* wall_eq_sys = wall_mesh_mapping->getMeshPartitioner()->getEquationSystems();
        std::unique_ptr<ExodusII_IO> exodus_io(uses_exodus ? new ExodusII_IO(mesh) : nullptr);
        std::unique_ptr<ExodusII_IO> bdry_io(uses_exodus ? new ExodusII_IO(*bdry_mesh_mapping->getBoundaryMesh()) :
                                                           nullptr);
        std::unique_ptr<ExodusII_IO> wall_io(uses_exodus ? new ExodusII_IO(*wall_mesh_mapping->getBoundaryMesh()) :
                                                           nullptr);
        bool from_restart = RestartManager::getManager()->isFromRestart();
        if (exodus_io) exodus_io->append(from_restart);
        if (bdry_io) bdry_io->append(from_restart);
        if (wall_io) wall_io->append(from_restart);

        // Create extra drawing variables
        auto var_db = VariableDatabase<NDIM>::getDatabase();
        Pointer<CellVariable<NDIM, double>> drag_var = new CellVariable<NDIM, double>("DragForce", NDIM),
                                            div_sig_var = new CellVariable<NDIM, double>("DivSig", NDIM),
                                            beta_var = new CellVariable<NDIM, double>("Beta"),
                                            sig_src_var = new CellVariable<NDIM, double>("SigSrcDraw", 3);
        const int drag_idx = var_db->registerVariableAndContext(drag_var, var_db->getContext("Draw"));
        const int div_sig_idx = var_db->registerVariableAndContext(div_sig_var, var_db->getContext("Draw"));
        const int beta_idx = var_db->registerVariableAndContext(beta_var, var_db->getContext("Draw"));
        const int sig_src_idx = var_db->registerVariableAndContext(sig_src_var, var_db->getContext("Draw"));

        visit_data_writer->registerPlotQuantity("DragForce", "VECTOR", drag_idx);
        visit_data_writer->registerPlotQuantity("DivSig", "VECTOR", div_sig_idx);
        visit_data_writer->registerPlotQuantity("Beta", "SCALAR", beta_idx);
        visit_data_writer->registerPlotQuantity("sig_src_0", "SCALAR", sig_src_idx, 0);
        visit_data_writer->registerPlotQuantity("sig_src_1", "SCALAR", sig_src_idx, 1);
        visit_data_writer->registerPlotQuantity("sig_src_2", "SCALAR", sig_src_idx, 2);

        // Initialize hierarchy configuration and data on all patches.
        ib_method_ops->initializeFEData();
        bdry_mesh_mapping->initializeFEData();
        time_integrator->initializePatchHierarchy(patch_hierarchy, gridding_algorithm);
        wall_mesh_mapping->getMeshPartitioner()->setPatchHierarchy(patch_hierarchy);
        wall_mesh_mapping->initializeFEData();
        sb_coupling_manager->fillInitialConditions();

        // Make sure source terms are set correctly.
        cohesion_relax->setWIdx(w_idx);
        phi_b_src_fcn->setWIdx(w_idx);
        phi_a_src_fcn->setWIdx(w_idx);
        bond_src_fcn->setWIdx(w_idx);

        // Deallocate initialization objects.
        app_initializer.setNull();

        // Print the input database contents to the log file.
        plog << "Input database:\n";
        input_db->printClassData(plog);

        int iteration_num = time_integrator->getIntegratorStep();
        double loop_time = time_integrator->getIntegratorTime();
        // Info for visualizations
        int viz_dump_iteration_num = 1;
        double viz_dump_time_interval = input_db->getDouble("VIZ_DUMP_TIME_INTERVAL");
        double next_viz_dump_time = 0.0;
        // Account for restarts...
        while (loop_time > 0.0 &&
               (next_viz_dump_time < loop_time || MathUtilities<double>::equalEps(loop_time, next_viz_dump_time)))
        {
            next_viz_dump_time += viz_dump_time_interval;
            viz_dump_iteration_num += 1;
        }

        // Main time step loop.
        double loop_time_end = time_integrator->getEndTime();
        double dt = 0.0;
        while (!MathUtilities<double>::equalEps(loop_time, loop_time_end) && time_integrator->stepsRemaining())
        {
            if (dump_viz_data &&
                (MathUtilities<double>::equalEps(loop_time, next_viz_dump_time) || loop_time >= next_viz_dump_time))
            {
                pout << "\nWriting visualization files...\n\n";
                if (uses_visit)
                {
                    const int phib_idx =
                        var_db->mapVariableAndContextToIndex(phi_b_var, sb_adv_diff_integrator->getCurrentContext());
                    const int ub_idx =
                        var_db->mapVariableAndContextToIndex(ub_var, adv_diff_integrator->getCurrentContext());
                    adv_diff_integrator->allocatePatchData(
                        drag_idx, loop_time, 0, patch_hierarchy->getFinestLevelNumber());
                    adv_diff_integrator->allocatePatchData(
                        div_sig_idx, loop_time, 0, patch_hierarchy->getFinestLevelNumber());
                    adv_diff_integrator->allocatePatchData(
                        beta_idx, loop_time, 0, patch_hierarchy->getFinestLevelNumber());
                    adv_diff_integrator->allocatePatchData(
                        sig_src_idx, loop_time, 0, patch_hierarchy->getFinestLevelNumber());
                    cohesion_relax->setPatchDataIndex(var_db->mapVariableAndContextToIndex(
                        cohesionStressForcing->getVariable(), adv_diff_integrator->getCurrentContext()));
                    cohesion_relax->setDataOnPatchHierarchy(sig_src_idx,
                                                            sig_src_var,
                                                            patch_hierarchy,
                                                            loop_time,
                                                            false,
                                                            0,
                                                            patch_hierarchy->getFinestLevelNumber());
                    compute_plot_quantities(drag_idx,
                                            div_sig_idx,
                                            beta_idx,
                                            patch_hierarchy,
                                            ins_integrator,
                                            cohesionStressForcing,
                                            bond_var,
                                            sb_adv_diff_integrator,
                                            ub_idx,
                                            phib_idx,
                                            clot_params,
                                            betaFcn);
                    time_integrator->setupPlotData();
                    visit_data_writer->writePlotData(patch_hierarchy, iteration_num, loop_time);
                    adv_diff_integrator->deallocatePatchData(drag_idx, 0, patch_hierarchy->getFinestLevelNumber());
                    adv_diff_integrator->deallocatePatchData(div_sig_idx, 0, patch_hierarchy->getFinestLevelNumber());
                    adv_diff_integrator->deallocatePatchData(beta_idx, 0, patch_hierarchy->getFinestLevelNumber());
                    adv_diff_integrator->deallocatePatchData(sig_src_idx, 0, patch_hierarchy->getFinestLevelNumber());
                }
                if (uses_exodus)
                {
                    exodus_io->write_timestep(base_mesh_filename, *eq_sys, viz_dump_iteration_num, loop_time);
                    bdry_io->write_timestep(bdry_mesh_filename, *bdry_eq_sys, viz_dump_iteration_num, loop_time);
                    wall_io->write_timestep(wall_mesh_filename, *wall_eq_sys, viz_dump_iteration_num, loop_time);
                }

                output_net_force(eq_sys, loop_time);
                next_viz_dump_time += viz_dump_time_interval;
                viz_dump_iteration_num += 1;
            }
            iteration_num = time_integrator->getIntegratorStep();
            loop_time = time_integrator->getIntegratorTime();

            pout << "\n";
            pout << "+++++++++++++++++++++++++++++++++++++++++++++++++++\n";
            pout << "At beginning of timestep # " << iteration_num << "\n";
            pout << "Simulation time is " << loop_time << "\n";

            dt = time_integrator->getMaximumTimeStepSize();
            time_integrator->advanceHierarchy(dt);
            loop_time += dt;

            pout << "\n";
            pout << "At end       of timestep # " << iteration_num << "\n";
            pout << "Simulation time is " << loop_time << "\n";
            pout << "+++++++++++++++++++++++++++++++++++++++++++++++++++\n";
            pout << "\n";

            // At specified intervals, write visualization and restart files,
            // print out timer data, and store hierarchy data for post
            // processing.
            iteration_num += 1;
            const bool last_step = !time_integrator->stepsRemaining();
            if (dump_restart_data && (iteration_num % restart_dump_interval == 0 || last_step))
            {
                pout << "\nWriting restart files...\n\n";
                RestartManager::getManager()->writeRestartFile(restart_dump_dirname, iteration_num);
                ib_method_ops->writeFEDataToRestartFile(restart_dump_dirname, iteration_num);
                wall_mesh_mapping->writeFEDataToRestartFile(restart_dump_dirname, iteration_num);
                bdry_mesh_mapping->writeFEDataToRestartFile(restart_dump_dirname, iteration_num);
            }
            if (dump_timer_data && (iteration_num % timer_dump_interval == 0 || last_step))
            {
                pout << "\nWriting timer data...\n\n";
                TimerManager::getManager()->print(plog);
            }
            if (dump_postproc_data && (iteration_num % postproc_data_dump_interval == 0 || last_step))
            {
                output_data(patch_hierarchy, time_integrator, iteration_num, loop_time, postproc_data_dump_dirname);
            }
        }

        // Close the force stream
        if (IBTK_MPI::getRank() == 0)
        {
            force_stream.close();
            tether_stream.close();
        }

        for (unsigned int d = 0; d < NDIM; ++d)
        {
            delete ub_bcs[d];
            delete u_bc_coefs[d];
        }
    } // cleanup dynamically allocated objects prior to shutdown
} // main

void
output_data(Pointer<PatchHierarchy<NDIM>> patch_hierarchy,
            Pointer<INSHierarchyIntegrator> ins_integrator,
            const int iteration_num,
            const double loop_time,
            const string& data_dump_dirname)
{
    plog << "writing hierarchy data at iteration " << iteration_num << " to disk" << endl;
    plog << "simulation time is " << loop_time << endl;
    string file_name = data_dump_dirname + "/" + "hier_data.";
    char temp_buf[128];
    sprintf(temp_buf, "%05d.samrai.%05d", iteration_num, IBTK_MPI::getRank());
    file_name += temp_buf;
    Pointer<HDFDatabase> hier_db = new HDFDatabase("hier_db");
    hier_db->create(file_name);
    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
    ComponentSelector hier_data;
    hier_data.setFlag(var_db->mapVariableAndContextToIndex(ins_integrator->getVelocityVariable(),
                                                           ins_integrator->getCurrentContext()));
    hier_data.setFlag(var_db->mapVariableAndContextToIndex(ins_integrator->getPressureVariable(),
                                                           ins_integrator->getCurrentContext()));
    patch_hierarchy->putToDatabase(hier_db->putDatabase("PatchHierarchy"), hier_data);
    hier_db->putDouble("loop_time", loop_time);
    hier_db->putInteger("iteration_num", iteration_num);
    hier_db->close();
    return;
} // output_data

void
output_net_force(EquationSystems* equation_systems, const double loop_time)
{
    const MeshBase& mesh = equation_systems->get_mesh();
    const unsigned int dim = mesh.mesh_dimension();
    double F_integral[NDIM];
    double T_integral[NDIM];
    for (unsigned int d = 0; d < NDIM; ++d)
    {
        F_integral[d] = 0.0;
        T_integral[d] = 0.0;
    }
    System* x_system;
    System* U_system;

    x_system = &equation_systems->get_system(IIMethod::COORDS_SYSTEM_NAME);
    U_system = &equation_systems->get_system(IIMethod::VELOCITY_SYSTEM_NAME);
    NumericVector<double>* x_vec = x_system->solution.get();
    NumericVector<double>* x_ghost_vec = x_system->current_local_solution.get();
    x_vec->localize(*x_ghost_vec);
    NumericVector<double>* U_vec = U_system->solution.get();
    NumericVector<double>* U_ghost_vec = U_system->current_local_solution.get();
    U_vec->localize(*U_ghost_vec);
    const DofMap& dof_map = x_system->get_dof_map();
    std::vector<std::vector<unsigned int>> dof_indices(NDIM);

    NumericVector<double>& X_vec = x_system->get_vector("INITIAL_COORDINATES");

    std::vector<std::vector<unsigned int>> WSS_o_dof_indices(NDIM);

    std::unique_ptr<FEBase> fe(FEBase::build(dim, dof_map.variable_type(0)));
    std::unique_ptr<QBase> qrule = QBase::build(QGAUSS, dim, SEVENTH);
    fe->attach_quadrature_rule(qrule.get());
    const vector<double>& JxW = fe->get_JxW();
    const vector<vector<double>>& phi = fe->get_phi();
    const vector<vector<VectorValue<double>>>& dphi = fe->get_dphi();

    std::vector<double> U_qp_vec(NDIM);
    std::vector<const std::vector<double>*> var_data(1);
    var_data[0] = &U_qp_vec;
    std::vector<const std::vector<libMesh::VectorValue<double>>*> grad_var_data;

    TensorValue<double> FF, FF_inv_trans;
    boost::multi_array<double, 2> x_node, X_node, U_node, TAU_node;

    VectorValue<double> F, N, U, n, x, X, TAU;

    const MeshBase::const_element_iterator el_begin = mesh.active_local_elements_begin();
    const MeshBase::const_element_iterator el_end = mesh.active_local_elements_end();
    for (MeshBase::const_element_iterator el_it = el_begin; el_it != el_end; ++el_it)
    {
        Elem* const elem = *el_it;
        fe->reinit(elem);
        for (unsigned int d = 0; d < NDIM; ++d)
        {
            dof_map.dof_indices(elem, dof_indices[d], d);
        }
        get_values_for_interpolation(x_node, *x_ghost_vec, dof_indices);
        get_values_for_interpolation(U_node, *U_ghost_vec, dof_indices);
        get_values_for_interpolation(X_node, X_vec, dof_indices);

        const unsigned int n_qp = qrule->n_points();
        for (unsigned int qp = 0; qp < n_qp; ++qp)
        {
            interpolate(X, qp, X_node, phi);
            interpolate(x, qp, x_node, phi);
            jacobian(FF, qp, x_node, dphi);
            interpolate(U, qp, U_node, phi);
            for (unsigned int d = 0; d < NDIM; ++d)
            {
                U_qp_vec[d] = U(d);
            }
            tether_force_function(F, n, N, FF, x, X, elem, 0, var_data, grad_var_data, loop_time, nullptr);

            for (int d = 0; d < NDIM; ++d)
            {
                F_integral[d] += F(d) * JxW[qp];
            }
        }
    }
    SAMRAI_MPI::sumReduction(F_integral, NDIM);
    SAMRAI_MPI::sumReduction(T_integral, NDIM);
    if (SAMRAI_MPI::getRank() == 0)
    {
        force_stream << loop_time << " " << -F_integral[0] << " " << -F_integral[1] << "\n";
    }
    return;
} // postprocess_data

void
compute_plot_quantities(const int drag_idx,
                        const int div_sig_idx,
                        const int beta_idx,
                        Pointer<PatchHierarchy<NDIM>> hierarchy,
                        Pointer<INSHierarchyIntegrator> ins_integrator,
                        Pointer<CFINSForcing> cf_forcing,
                        Pointer<CellVariable<NDIM, double>> bond_var,
                        Pointer<AdvDiffHierarchyIntegrator> bond_integrator,
                        const int ub_idx,
                        const int phib_idx,
                        const BoundClotParams& clot_params,
                        const BetaFcn& beta_fcn)

{
    const int coarsest_ln = 0, finest_ln = hierarchy->getFinestLevelNumber();
    Pointer<AdvDiffHierarchyIntegrator> adv_diff_integrator = cf_forcing->getAdvDiffHierarchyIntegrator();
    auto var_db = VariableDatabase<NDIM>::getDatabase();
    const int uf_idx = var_db->mapVariableAndContextToIndex(ins_integrator->getVelocityVariable(),
                                                            ins_integrator->getCurrentContext());
    const int sig_idx = cf_forcing->getVariableIdx();

    const int bond_idx = var_db->mapVariableAndContextToIndex(bond_var, bond_integrator->getCurrentContext());
    const int bond_scr_idx = var_db->mapVariableAndContextToIndex(bond_var, bond_integrator->getScratchContext());

    // Need to fill in ghost data for divergence
    const int sig_scr_idx =
        var_db->mapVariableAndContextToIndex(cf_forcing->getVariable(), adv_diff_integrator->getScratchContext());
    const bool sig_scr_is_allocated = adv_diff_integrator->isAllocatedPatchData(sig_scr_idx, coarsest_ln, finest_ln);
    if (!sig_scr_is_allocated) adv_diff_integrator->allocatePatchData(sig_scr_idx, 0.0, coarsest_ln, finest_ln);
    const bool bond_scr_is_allocated = bond_integrator->isAllocatedPatchData(bond_scr_idx, coarsest_ln, finest_ln);
    if (!bond_scr_is_allocated) bond_integrator->allocatePatchData(bond_scr_idx, 0.0, coarsest_ln, finest_ln);
    using ITC = HierarchyGhostCellInterpolation::InterpolationTransactionComponent;
    std::vector<ITC> ghost_cell_comp(2);
    ghost_cell_comp[0] =
        ITC(sig_scr_idx, sig_idx, "CONSERVATIVE_LINEAR_REFINE", false, "NONE", "LINEAR", false, nullptr);
    ghost_cell_comp[1] =
        ITC(bond_scr_idx, bond_idx, "CONSERVATIVE_LINEAR_REFINE", false, "NONE", "LINEAR", false, nullptr);
    HierarchyGhostCellInterpolation ghost_cell_fill;
    ghost_cell_fill.initializeOperatorState(ghost_cell_comp, hierarchy, coarsest_ln, finest_ln);
    ghost_cell_fill.fillData(0.0);

    // Compute forces
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        Pointer<PatchLevel<NDIM>> level = hierarchy->getPatchLevel(ln);
        for (PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            Pointer<Patch<NDIM>> patch = level->getPatch(p());
            Pointer<CellData<NDIM, double>> ub_data = patch->getPatchData(ub_idx);
            Pointer<SideData<NDIM, double>> uf_data = patch->getPatchData(uf_idx);
            Pointer<CellData<NDIM, double>> drag_data = patch->getPatchData(drag_idx);
            Pointer<CellData<NDIM, double>> phib_data = patch->getPatchData(phib_idx);
            Pointer<CellData<NDIM, double>> bond_data = patch->getPatchData(bond_scr_idx);
            Pointer<CellData<NDIM, double>> beta_data = patch->getPatchData(beta_idx);

            Pointer<CellData<NDIM, double>> div_sig_data = patch->getPatchData(div_sig_idx);
            Pointer<CellData<NDIM, double>> sig0_data = patch->getPatchData(sig_scr_idx);

            Pointer<CellData<NDIM, double>> sig_data =
                convertToStress(*sig0_data, *bond_data, patch, clot_params.S0, clot_params.R0, true);

            Pointer<CartesianPatchGeometry<NDIM>> pgeom = patch->getPatchGeometry();
            const double* const dx = pgeom->getDx();

            for (CellIterator<NDIM> ci(patch->getBox()); ci; ci++)
            {
                const CellIndex<NDIM>& idx = ci();

                // Compute drag force
                const double th = clot_params.vol_pl * (*phib_data)(idx);
                const double xi = clot_params.drag_coef * th * th / std::pow(1.0 - th, 3.0);
                for (int d = 0; d < NDIM; ++d)
                {
                    const SideIndex<NDIM> idx_u(idx, d, 1), idx_l(idx, d, 0);
                    (*drag_data)(idx, d) = xi * (0.5 * ((*uf_data)(idx_u) + (*uf_data)(idx_l)) - (*ub_data)(idx, d));
                }

                // Compute divergence
                IntVector<NDIM> x(1, 0), y(0, 1);
                (*div_sig_data)(idx, 0) = ((*sig_data)(idx + x, 0) - (*sig_data)(idx - x, 0)) / (2.0 * dx[0]) +
                                          ((*sig_data)(idx + y, 2) - (*sig_data)(idx - y, 2)) / (2.0 * dx[1]);
                (*div_sig_data)(idx, 1) = ((*sig_data)(idx + x, 2) - (*sig_data)(idx - x, 2)) / (2.0 * dx[0]) +
                                          ((*sig_data)(idx + y, 1) - (*sig_data)(idx - y, 1)) / (2.0 * dx[1]);

                // Compute breaking rate
                const double tr = (*sig0_data)(idx, 0) + (*sig0_data)(idx, 1);
                double avg_y = 0.0;
                if (tr > 1.0e-12 && (*bond_data)(idx) > 1.0e-12)
                    avg_y = std::sqrt(tr * 2.0 / (clot_params.S0 * (*bond_data)(idx) + 1.0e-12));
                (*beta_data)(idx) = beta_fcn.beta(avg_y);
            }
        }
    }

    if (!sig_scr_is_allocated) adv_diff_integrator->deallocatePatchData(sig_scr_idx, coarsest_ln, finest_ln);
    if (!bond_scr_is_allocated) bond_integrator->deallocatePatchData(bond_scr_idx, coarsest_ln, finest_ln);
}
