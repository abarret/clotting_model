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

#include <ibamr/AdvDiffSemiImplicitHierarchyIntegrator.h>
#include <ibamr/CFINSForcing.h>
#include <ibamr/IBExplicitHierarchyIntegrator.h>
#include <ibamr/IBFEMethod.h>
#include <ibamr/INSCollocatedHierarchyIntegrator.h>
#include <ibamr/INSStaggeredHierarchyIntegrator.h>
#include <ibamr/app_namespaces.h>

#include <ibtk/AppInitializer.h>
#include <ibtk/CartGridFunctionSet.h>
#include <ibtk/IBTKInit.h>
#include <ibtk/IBTK_MPI.h>
#include <ibtk/muParserCartGridFunction.h>
#include <ibtk/muParserRobinBcCoefs.h>

#include "libmesh/mesh_modification.h"
#include "libmesh/mesh_refinement.h"
#include "libmesh/mesh_tools.h"
#include <libmesh/boundary_info.h>
#include <libmesh/boundary_mesh.h>
#include <libmesh/elem.h>
#include <libmesh/equation_systems.h>
#include <libmesh/exodusII_io.h>
#include <libmesh/mesh.h>
#include <libmesh/mesh_generation.h>
#include <libmesh/mesh_triangle_interface.h>

#include <petscsys.h>

#include <BergerRigoutsos.h>
#include <CartesianGridGeometry.h>
#include <LoadBalancer.h>
#include <SAMRAI_config.h>
#include <StandardTagAndInitialize.h>

#include <fstream>
#include <iostream>

// Local Headers
#include "BondSource.h"
#include "BoundaryMeshMapping.h"
#include "CohesionStressRHS.h"
#include "PlateletSource.h"

// Function prototypes
void output_data(Pointer<PatchHierarchy<NDIM>> patch_hierarchy,
                 Pointer<INSHierarchyIntegrator> ins_integrator,
                 const int iteration_num,
                 const double loop_time,
                 const string& data_dump_dirname);

struct BetaFcn
{
public:
    BetaFcn(const double R0, const double beta_0, const double beta_1) : d_R0(R0), d_beta_0(beta_0), d_beta_1(beta_1)
    {
        // intentionally blank
        return;
    }
    double beta(const double eps)
    {
        return d_beta_0;
        if (eps > d_R0)
        {
            return d_beta_0 * std::exp(d_beta_1 * (eps - d_R0));
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

namespace ModelData
{
static double kappa = 1.0e6;
static double eta = 0.0;
static double beta_s = 1.0e3;
static double dx = -1.0;
static bool ERROR_ON_MOVE = false;
void
tether_force_function(VectorValue<double>& F,
                      const TensorValue<double>& /*FF*/,
                      const libMesh::Point& x,
                      const libMesh::Point& X,
                      Elem* const /*elem*/,
                      const vector<const vector<double>*>& var_data,
                      const vector<const vector<VectorValue<double>>*>& /*grad_var_data*/,
                      double /*time*/,
                      void* /*ctx*/)
{
    // Tether to initial location using damped springs
    // x is current location
    // X is reference location
    // U is velocity
    const std::vector<double>& U = *var_data[0];
    for (unsigned int d = 0; d < NDIM; ++d) F(d) = kappa * (X(d) - x(d)) - eta * U[d];
    // Check to see how much structure has moved. If more than a quarter of a grid cell, quit.
    std::vector<double> d = { std::abs(x(0) - X(0)), std::abs(x(1) - X(1)) };
    if (ERROR_ON_MOVE && ((d[0] > 0.25 * dx) || (d[1] > 0.25 * dx))) TBOX_ERROR("Structure has moved too much.\n");
}

void
tether_penalty_stress_fcn(TensorValue<double>& PP,
                          const TensorValue<double>& FF,
                          const libMesh::Point& /*X*/,
                          const libMesh::Point& /*s*/,
                          Elem* const /*elem*/,
                          const vector<const vector<double>*>& /*var_data*/,
                          const vector<const vector<VectorValue<double>>*>& /*grad_var_data*/,
                          double /*time*/,
                          void* ctx)
{
    const TensorValue<double> FF_inv_trans = tensor_inverse_transpose(FF, NDIM);
    double J = FF.det();
    PP = beta_s * J * log(J) * FF_inv_trans;
    return;
}

} // namespace ModelData
using namespace ModelData;

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
        const int viz_dump_interval = app_initializer->getVizDumpInterval();
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
        const std::string exodus_filename = app_initializer->getExodusIIFilename();

        const bool dump_restart_data = app_initializer->dumpRestartData();
        const int restart_dump_interval = app_initializer->getRestartDumpInterval();
        const string restart_dump_dirname = app_initializer->getRestartDumpDirectory();

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
        solid_mesh.prepare_for_use();
        solid_mesh.print_info();

        Mesh& mesh = solid_mesh;

        // Get structure parameters
        kappa = input_db->getDouble("KAPPA");
        eta = input_db->getDouble("ETA");

        // Create major algorithm and data objects that comprise the
        // application.  These objects are configured from the input database
        // and, if this is a restarted run, from the restart database.
        Pointer<INSHierarchyIntegrator> ins_integrator = new INSStaggeredHierarchyIntegrator(
            "INSStaggeredHierarchyIntegrator",
            app_initializer->getComponentDatabase("INSStaggeredHierarchyIntegrator"));
        Pointer<AdvDiffSemiImplicitHierarchyIntegrator> adv_diff_integrator;
        adv_diff_integrator = new AdvDiffSemiImplicitHierarchyIntegrator(
            "AdvDiffSemiImplicitHierarchyIntegrator",
            app_initializer->getComponentDatabase("AdvDiffSemiImplicitHierarchyIntegrator"));
        Pointer<IBFEMethod> ib_method_ops =
            new IBFEMethod("IBFEMethod",
                           app_initializer->getComponentDatabase("IBFEMethod"),
                           &mesh,
                           app_initializer->getComponentDatabase("GriddingAlgorithm")->getInteger("max_levels"));
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
        // Create initial condition specification objects.
        Pointer<CartGridFunction> u_init = new muParserCartGridFunction(
            "u_init", app_initializer->getComponentDatabase("VelocityInitialConditions"), grid_geometry);
        ins_integrator->registerVelocityInitialConditions(u_init);
        Pointer<CartGridFunction> p_init = new muParserCartGridFunction(
            "p_init", app_initializer->getComponentDatabase("PressureInitialConditions"), grid_geometry);
        ins_integrator->registerPressureInitialConditions(p_init);

        // Configure the IBFE solver
        ib_method_ops->initializeFEEquationSystems();
        std::vector<int> vars(NDIM);
        for (unsigned int d = 0; d < NDIM; ++d) vars[d] = d;
        std::vector<SystemData> sys_data = { SystemData(IBFEMethod::VELOCITY_SYSTEM_NAME, vars) };
        IBFEMethod::LagBodyForceFcnData body_fcn_data(tether_force_function, sys_data);
        ib_method_ops->registerLagBodyForceFunction(body_fcn_data);
        beta_s = input_db->getDouble("BETA_S");
        ib_method_ops->registerPK1StressFunction(tether_penalty_stress_fcn);

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

        // Set up visualization plot file writers.
        Pointer<VisItDataWriter<NDIM>> visit_data_writer = app_initializer->getVisItDataWriter();
        if (uses_visit)
        {
            time_integrator->registerVisItDataWriter(visit_data_writer);
        }

        // Create advected quantities
        Pointer<CellVariable<NDIM, double>> phi_u_var = new CellVariable<NDIM, double>("phi_u");
        Pointer<CellVariable<NDIM, double>> phi_a_var = new CellVariable<NDIM, double>("phi_a");
        Pointer<CellVariable<NDIM, double>> bond_var = new CellVariable<NDIM, double>("bond");

        // Set up Cohesion stress tensor
        Pointer<CFINSForcing> cohesionStressForcing =
            new CFINSForcing("CohesionStressForcing",
                             app_initializer->getComponentDatabase("CohesionStress"),
                             ins_integrator,
                             grid_geometry,
                             adv_diff_integrator,
                             visit_data_writer);
        ins_integrator->registerBodyForceFunction(cohesionStressForcing);
        Pointer<CohesionStressRHS> cohesion_relax =
            new CohesionStressRHS(phi_u_var,
                                  phi_a_var,
                                  bond_var,
                                  "CohesionRHS",
                                  app_initializer->getComponentDatabase("CohesionRHS"),
                                  adv_diff_integrator);
        cohesionStressForcing->registerRelaxationOperator(cohesion_relax);

        double beta_0 = input_db->getDouble("BETA_0");
        double beta_1 = input_db->getDouble("BETA_1");
        double r_0 = input_db->getDouble("R_0");
        BetaFcn betaFcn(r_0, beta_0, beta_1);
        cohesion_relax->registerBetaFcn(beta_wrapper, static_cast<void*>(&betaFcn));
        Pointer<CellVariable<NDIM, double>> sig_var = cohesionStressForcing->getVariable();
        // Set up platelet concentrations and bond information

        // Unactivated
        adv_diff_integrator->registerTransportedQuantity(phi_u_var);
        Pointer<CellVariable<NDIM, double>> phi_u_src_var = new CellVariable<NDIM, double>("phi_u_src");
        adv_diff_integrator->setAdvectionVelocity(phi_u_var, ins_integrator->getAdvectionVelocityVariable());
        adv_diff_integrator->setDiffusionCoefficient(phi_u_var, input_db->getDouble("UNACTIVATED_DIFFUSION_COEF"));
        Pointer<PlateletSource> phi_u_src_fcn = new PlateletSource(
            phi_u_var, phi_a_var, app_initializer->getComponentDatabase("ActivatedPlatelets"), adv_diff_integrator);
        phi_u_src_fcn->setSign(false);
        phi_u_src_fcn->setKernel(BSPLINE_3);
        adv_diff_integrator->registerSourceTerm(phi_u_src_var);
        adv_diff_integrator->setSourceTermFunction(phi_u_src_var, phi_u_src_fcn);
        adv_diff_integrator->setSourceTerm(phi_u_var, phi_u_src_var);
        Pointer<RobinBcCoefStrategy<NDIM>> phi_u_bcs;
        if (grid_geometry->getPeriodicShift().min() == 0)
            phi_u_bcs = new muParserRobinBcCoefs("UnactivatedPlateletsBcs",
                                                 app_initializer->getComponentDatabase("UnactivatedPlateletBcs"),
                                                 grid_geometry);
        adv_diff_integrator->setPhysicalBcCoef(phi_u_var, phi_u_bcs.getPointer());

        // Activated
        adv_diff_integrator->registerTransportedQuantity(phi_a_var);
        Pointer<CellVariable<NDIM, double>> phi_a_src_var = new CellVariable<NDIM, double>("phi_a_src");
        adv_diff_integrator->setAdvectionVelocity(phi_a_var, ins_integrator->getAdvectionVelocityVariable());
        adv_diff_integrator->setDiffusionCoefficient(phi_a_var, input_db->getDouble("ACTIVATED_DIFFUSION_COEF"));
        Pointer<PlateletSource> phi_a_src_fcn = new PlateletSource(
            phi_u_var, phi_a_var, app_initializer->getComponentDatabase("ActivatedPlatelets"), adv_diff_integrator);
        phi_a_src_fcn->setSign(true);
        phi_a_src_fcn->setKernel(BSPLINE_3);
        adv_diff_integrator->registerSourceTerm(phi_a_src_var);
        adv_diff_integrator->setSourceTermFunction(phi_a_src_var, phi_a_src_fcn);
        adv_diff_integrator->setSourceTerm(phi_a_var, phi_a_src_var);
        Pointer<RobinBcCoefStrategy<NDIM>> phi_a_bcs;
        if (grid_geometry->getPeriodicShift().min() == 0)
            phi_a_bcs = new muParserRobinBcCoefs(
                "ActivatedPlateletsBcs", app_initializer->getComponentDatabase("ActivatedPlateletBcs"), grid_geometry);
        adv_diff_integrator->setPhysicalBcCoef(phi_a_var, phi_a_bcs.getPointer());

        // Bonds
        adv_diff_integrator->registerTransportedQuantity(bond_var);
        Pointer<CellVariable<NDIM, double>> bond_src_var = new CellVariable<NDIM, double>("bond_src");
        adv_diff_integrator->setAdvectionVelocity(bond_var, ins_integrator->getAdvectionVelocityVariable());
        adv_diff_integrator->setDiffusionCoefficient(bond_var, 0.0);
        Pointer<BondSource> bond_src_fcn = new BondSource(phi_u_var,
                                                          phi_a_var,
                                                          bond_var,
                                                          sig_var,
                                                          app_initializer->getComponentDatabase("BondVariable"),
                                                          adv_diff_integrator);
        bond_src_fcn->registerBetaFcn(beta_wrapper, static_cast<void*>(&betaFcn));
        adv_diff_integrator->registerSourceTerm(bond_src_var);
        adv_diff_integrator->setSourceTermFunction(bond_src_var, bond_src_fcn);
        adv_diff_integrator->setSourceTerm(bond_var, bond_src_var);
        Pointer<RobinBcCoefStrategy<NDIM>> bond_bcs;
        if (grid_geometry->getPeriodicShift().min() == 0)
            bond_bcs = new muParserRobinBcCoefs(
                "BondVariableBcs", app_initializer->getComponentDatabase("BondVariableBcs"), grid_geometry);
        adv_diff_integrator->setPhysicalBcCoef(bond_var, bond_bcs.getPointer());

        // Wall sites
        FEDataManager* vol_data_manager = ib_method_ops->getFEDataManager();
        BoundaryMeshMapping bdry_mesh_mapping(
            "BoundaryMesh", app_initializer->getComponentDatabase("BoundaryMesh"), &mesh, vol_data_manager);
        const int w_idx = bdry_mesh_mapping.getWallSitesPatchIndex();
        auto update_bdry_mesh = [](const double current_time,
                                   const double new_time,
                                   bool /*skip_synchronize_new_state_data*/,
                                   int /*num_cycles*/,
                                   void* ctx) {
            auto bdry_mesh_mapping = static_cast<BoundaryMeshMapping*>(ctx);
            bdry_mesh_mapping->updateBoundaryLocation(current_time, new_time, false);
        };
        auto regrid_callback = [](SAMRAI::tbox::Pointer<SAMRAI::hier::BasePatchHierarchy<NDIM>> /*hierarchy*/,
                                  const double /*time*/,
                                  const bool /*initial_time*/,
                                  void* ctx) {
            auto bdry_mesh_mapping = static_cast<BoundaryMeshMapping*>(ctx);
            bdry_mesh_mapping->spreadWallSites();
        };
        time_integrator->registerPostprocessIntegrateHierarchyCallback(update_bdry_mesh,
                                                                       static_cast<void*>(&bdry_mesh_mapping));
        time_integrator->registerRegridHierarchyCallback(regrid_callback, static_cast<void*>(&bdry_mesh_mapping));
        visit_data_writer->registerPlotQuantity("wall_sites", "SCALAR", w_idx);

        EquationSystems* eq_sys = ib_method_ops->getFEDataManager()->getEquationSystems();
        std::unique_ptr<ExodusII_IO> exodus_io(uses_exodus ? new ExodusII_IO(mesh) : nullptr);

        // Initialize hierarchy configuration and data on all patches.
        ib_method_ops->initializeFEData();
        time_integrator->initializePatchHierarchy(patch_hierarchy, gridding_algorithm);
        bdry_mesh_mapping.initializeEquationSystems();

        // Make sure source terms are set correctly.
        cohesion_relax->setOmegaIdx(w_idx);
        phi_u_src_fcn->setWIdx(w_idx);
        phi_a_src_fcn->setWIdx(w_idx);
        bond_src_fcn->setWIdx(w_idx);

        // Deallocate initialization objects.
        app_initializer.setNull();

        // Print the input database contents to the log file.
        plog << "Input database:\n";
        input_db->printClassData(plog);

        // Write out initial visualization data.
        int iteration_num = time_integrator->getIntegratorStep();
        double loop_time = time_integrator->getIntegratorTime();
        if (dump_viz_data)
        {
            pout << "\n\nWriting visualization files...\n\n";
            if (uses_visit)
            {
                time_integrator->setupPlotData();
                visit_data_writer->writePlotData(patch_hierarchy, iteration_num, loop_time);
            }
            if (uses_exodus)
            {
                exodus_io->write_timestep(exodus_filename, *eq_sys, iteration_num / viz_dump_interval + 1, loop_time);
            }
        }

        // Main time step loop.
        double loop_time_end = time_integrator->getEndTime();
        double dt = 0.0;
        while (!MathUtilities<double>::equalEps(loop_time, loop_time_end) && time_integrator->stepsRemaining())
        {
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
            if (dump_viz_data && (iteration_num % viz_dump_interval == 0 || last_step))
            {
                pout << "\nWriting visualization files...\n\n";
                if (uses_visit)
                {
                    time_integrator->setupPlotData();
                    visit_data_writer->writePlotData(patch_hierarchy, iteration_num, loop_time);
                }
                if (uses_exodus)
                {
                    exodus_io->write_timestep(
                        exodus_filename, *eq_sys, iteration_num / viz_dump_interval + 1, loop_time);
                }
            }
            if (dump_restart_data && (iteration_num % restart_dump_interval == 0 || last_step))
            {
                pout << "\nWriting restart files...\n\n";
                RestartManager::getManager()->writeRestartFile(restart_dump_dirname, iteration_num);
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
