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

#include <clot/BondSource.h>
#include <clot/BoundaryMeshMapping.h>
#include <clot/CohesionStressRHS.h>
#include <clot/PlateletSource.h>
#include <clot/app_namespaces.h>

#include <ADS/CutCellVolumeMeshMapping.h>
#include <ADS/GeneralBoundaryMeshMapping.h>
#include <ADS/LSCutCellLaplaceOperator.h>
#include <ADS/LSFromMesh.h>
#include <ADS/SBAdvDiffIntegrator.h>

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

#include <boost/math/special_functions/ellint_2.hpp>

#include <BergerRigoutsos.h>
#include <CartesianGridGeometry.h>
#include <LoadBalancer.h>
#include <SAMRAI_config.h>
#include <StandardTagAndInitialize.h>

#include <fstream>
#include <iostream>

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
    double beta(const double eps)
    {
        if (eps > d_R0)
        {
            // return d_beta_0 * std::exp(d_beta_1 * (eps - d_R0));
            return d_beta_0 * (1.0 + d_beta_1 * (eps - d_R0));
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
    Pointer<hier::Variable<NDIM>> phi_u_var;
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
                      Elem* const elem,
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
    if (ERROR_ON_MOVE && d.norm() > 0.25 * dx)
    {
        pout << "displacement: " << d << "\n";
        pout << "max allowed:  " << 0.25 * dx << "\n";
        pout << "dx:           " << dx << "\n";
        pout << "On elem:      " << *elem << "\n";
        pout << "Cur Location: " << x << "\n";
        pout << "Ref location: " << X << "\n";
        TBOX_ERROR("Structure has moved too much.\n");
    }

} // namespace ModelData
using namespace ModelData;

static std::ofstream force_stream, tether_stream;
void output_net_force(EquationSystems* eq_sys, const double loop_time);

void generateMesh(Mesh& mesh, unsigned int N, double a, double L);
std::pair<double, int> find_s_from_t(double L, double a, double t);

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
        const std::string exodus_filename = app_initializer->getExodusIIFilename();

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
        Mesh mesh(init.comm(), NDIM);
        unsigned int N = input_db->getInteger("NUM_PTS");
        double a = input_db->getDouble("DEPTH");
        double L = (input_db->getDouble("LX") - 2.0) / 2.0;
        generateMesh(mesh, N, a, L);

        // Get structure parameters
        kappa = input_db->getDouble("KAPPA");
        eta = input_db->getDouble("ETA");
        dx = input_db->getDouble("DX");
        ERROR_ON_MOVE = input_db->getBool("ERROR_ON_MOVE");

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
                                  sb_adv_diff_integrator);
        cohesionStressForcing->registerRelaxationOperator(cohesion_relax);

        double beta_0 = input_db->getDouble("BETA_0");
        double beta_1 = input_db->getDouble("BETA_1");
        double r_0 = input_db->getDouble("R_0");
        BetaFcn betaFcn(r_0, beta_0, beta_1);
        cohesion_relax->registerBetaFcn(beta_wrapper, static_cast<void*>(&betaFcn));
        Pointer<CellVariable<NDIM, double>> sig_var = cohesionStressForcing->getVariable();

        Pointer<LSCutCellLaplaceOperator> phi_u_rhs_oper = new LSCutCellLaplaceOperator(
                                              "PhiULSCutCellRHSOperator",
                                              app_initializer->getComponentDatabase("LSCutCellOperator"),
                                              false),
                                          phi_u_oper = new LSCutCellLaplaceOperator(
                                              "PhiULSCutCellOperator",
                                              app_initializer->getComponentDatabase("LSCutCellOperator"),
                                              false);
        Pointer<PETScKrylovPoissonSolver> phi_u_solver = new PETScKrylovPoissonSolver(
            "PhiUPoissonSolver", app_initializer->getComponentDatabase("PoissonSolver"), "poisson_solve_");
        phi_u_solver->setOperator(phi_u_oper);

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
        // Set up platelet concentrations and bond information

        // Unactivated
        sb_adv_diff_integrator->registerTransportedQuantity(phi_u_var);
        sb_adv_diff_integrator->restrictToLevelSet(phi_u_var, ls_var);
        Pointer<CellVariable<NDIM, double>> phi_u_src_var = new CellVariable<NDIM, double>("phi_u_src");
        sb_adv_diff_integrator->setAdvectionVelocity(phi_u_var, ins_integrator->getAdvectionVelocityVariable());
        sb_adv_diff_integrator->setDiffusionCoefficient(phi_u_var, input_db->getDouble("UNACTIVATED_DIFFUSION_COEF"));
        Pointer<PlateletSource> phi_u_src_fcn =
            new PlateletSource("UnactivatedPlatelets",
                               phi_u_var,
                               phi_a_var,
                               app_initializer->getComponentDatabase("ActivatedPlatelets"),
                               sb_adv_diff_integrator);
        phi_u_src_fcn->setSign(false);
        phi_u_src_fcn->setKernel(BSPLINE_3);
        sb_adv_diff_integrator->registerSourceTerm(phi_u_src_var);
        sb_adv_diff_integrator->setSourceTermFunction(phi_u_src_var, phi_u_src_fcn);
        sb_adv_diff_integrator->setSourceTerm(phi_u_var, phi_u_src_var);
        Pointer<RobinBcCoefStrategy<NDIM>> phi_u_bcs;
        if (grid_geometry->getPeriodicShift().min() == 0)
            phi_u_bcs = new muParserRobinBcCoefs("UnactivatedPlateletsBcs",
                                                 app_initializer->getComponentDatabase("UnactivatedPlateletBcs"),
                                                 grid_geometry);
        sb_adv_diff_integrator->setPhysicalBcCoef(phi_u_var, phi_u_bcs.getPointer());
        sb_adv_diff_integrator->setHelmholtzSolver(phi_u_var, phi_u_solver);
        sb_adv_diff_integrator->setHelmholtzRHSOperator(phi_u_var, phi_u_rhs_oper);

        // Activated
        sb_adv_diff_integrator->registerTransportedQuantity(phi_a_var);
        sb_adv_diff_integrator->restrictToLevelSet(phi_a_var, ls_var);
        Pointer<CellVariable<NDIM, double>> phi_a_src_var = new CellVariable<NDIM, double>("phi_a_src");
        sb_adv_diff_integrator->setAdvectionVelocity(phi_a_var, ins_integrator->getAdvectionVelocityVariable());
        sb_adv_diff_integrator->setDiffusionCoefficient(phi_a_var, input_db->getDouble("ACTIVATED_DIFFUSION_COEF"));
        Pointer<PlateletSource> phi_a_src_fcn =
            new PlateletSource("ActivatedPlatelets",
                               phi_u_var,
                               phi_a_var,
                               app_initializer->getComponentDatabase("ActivatedPlatelets"),
                               sb_adv_diff_integrator);
        phi_a_src_fcn->setSign(true);
        phi_a_src_fcn->setKernel(BSPLINE_3);
        sb_adv_diff_integrator->registerSourceTerm(phi_a_src_var);
        sb_adv_diff_integrator->setSourceTermFunction(phi_a_src_var, phi_a_src_fcn);
        sb_adv_diff_integrator->setSourceTerm(phi_a_var, phi_a_src_var);
        Pointer<RobinBcCoefStrategy<NDIM>> phi_a_bcs;
        if (grid_geometry->getPeriodicShift().min() == 0)
            phi_a_bcs = new muParserRobinBcCoefs(
                "ActivatedPlateletsBcs", app_initializer->getComponentDatabase("ActivatedPlateletBcs"), grid_geometry);
        sb_adv_diff_integrator->setPhysicalBcCoef(phi_a_var, phi_a_bcs.getPointer());
        sb_adv_diff_integrator->setHelmholtzSolver(phi_a_var, phi_a_solver);
        sb_adv_diff_integrator->setHelmholtzRHSOperator(phi_a_var, phi_a_rhs_oper);

        // Bonds
        sb_adv_diff_integrator->registerTransportedQuantity(bond_var);
        sb_adv_diff_integrator->restrictToLevelSet(bond_var, ls_var);
        Pointer<CellVariable<NDIM, double>> bond_src_var = new CellVariable<NDIM, double>("bond_src");
        sb_adv_diff_integrator->setAdvectionVelocity(bond_var, ins_integrator->getAdvectionVelocityVariable());
        sb_adv_diff_integrator->setDiffusionCoefficient(bond_var, 0.0);
        Pointer<BondSource> bond_src_fcn = new BondSource(phi_u_var,
                                                          phi_a_var,
                                                          bond_var,
                                                          sig_var,
                                                          app_initializer->getComponentDatabase("BondVariable"),
                                                          adv_diff_integrator,
                                                          sb_adv_diff_integrator);
        bond_src_fcn->registerBetaFcn(beta_wrapper, static_cast<void*>(&betaFcn));
        sb_adv_diff_integrator->registerSourceTerm(bond_src_var);
        sb_adv_diff_integrator->setSourceTermFunction(bond_src_var, bond_src_fcn);
        sb_adv_diff_integrator->setSourceTerm(bond_var, bond_src_var);
        Pointer<RobinBcCoefStrategy<NDIM>> bond_bcs;
        if (grid_geometry->getPeriodicShift().min() == 0)
            bond_bcs = new muParserRobinBcCoefs(
                "BondVariableBcs", app_initializer->getComponentDatabase("BondVariableBcs"), grid_geometry);
        sb_adv_diff_integrator->setPhysicalBcCoef(bond_var, bond_bcs.getPointer());
        sb_adv_diff_integrator->setHelmholtzSolver(bond_var, bond_solver);
        sb_adv_diff_integrator->setHelmholtzRHSOperator(bond_var, bond_rhs_oper);

        // Wall sites
        FEDataManager* fe_data_manager = ib_method_ops->getFEDataManager();
        auto bdry_mesh_mapping =
            std::make_shared<BoundaryMeshMapping>("BoundaryMesh",
                                                  app_initializer->getComponentDatabase("BoundaryMesh"),
                                                  &mesh,
                                                  fe_data_manager,
                                                  restart_read_dirname,
                                                  restart_restore_num);
        bdry_mesh_mapping->initializeEquationSystems();
        auto general_mesh_mapping =
            std::make_shared<GeneralBoundaryMeshMapping>("GeneralMesh",
                                                         app_initializer->getComponentDatabase("BoundaryMesh"),
                                                         &mesh,
                                                         restart_read_dirname,
                                                         restart_restore_num);
        general_mesh_mapping->initializeEquationSystems();
        Pointer<CutCellMeshMapping> cut_cell_mapping =
            new CutCellVolumeMeshMapping("CutCellMapping",
                                         app_initializer->getComponentDatabase("CutCellMapping"),
                                         general_mesh_mapping->getMeshPartitioner(0));
        Pointer<LSFromMesh> ls_fcn = new LSFromMesh("LSFcn", patch_hierarchy, cut_cell_mapping, false);
        ls_fcn->registerBdryFcn(bdry_fcn);
        sb_adv_diff_integrator->registerLevelSetVolFunction(ls_var, ls_fcn);
        sb_adv_diff_integrator->registerGeneralBoundaryMeshMapping(general_mesh_mapping);
        const int w_idx = bdry_mesh_mapping->getWallSitesPatchIndex();
        std::pair<BoundaryMeshMapping*, GeneralBoundaryMeshMapping*> mesh_mappings =
            std::make_pair(bdry_mesh_mapping.get(), general_mesh_mapping.get());
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
            auto bdry_mesh_mapping = static_cast<BoundaryMeshMapping*>(ctx);
            bdry_mesh_mapping->spreadWallSites();
        };
        time_integrator->registerPostprocessIntegrateHierarchyCallback(update_bdry_mesh,
                                                                       static_cast<void*>(&mesh_mappings));
        time_integrator->registerRegridHierarchyCallback(regrid_callback, static_cast<void*>(bdry_mesh_mapping.get()));
        visit_data_writer->registerPlotQuantity("wall_sites", "SCALAR", w_idx);

        // define variable blob, which stores the patch hierarchy, variable context, and clamped variables.
        ClampVars clamp_vars = {
            patch_hierarchy,                         // patch hierarchy
            sb_adv_diff_integrator->getNewContext(), // variable context
            phi_u_var,                               // unactivated platelets var
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
                        Pointer<CellData<NDIM, double>> phi_u_data =
                            patch->getPatchData(clamp_vars->phi_u_var, clamp_vars->context);
                        Pointer<CellData<NDIM, double>> phi_a_data =
                            patch->getPatchData(clamp_vars->phi_a_var, clamp_vars->context);
                        Pointer<CellData<NDIM, double>> bond_data =
                            patch->getPatchData(clamp_vars->bond_var, clamp_vars->context);
                        // iterate cells
                        for (CellIterator<NDIM> ci(patch->getBox()); ci; ci++)
                        {
                            const CellIndex<NDIM>& idx = ci();
                            (*phi_u_data)(idx) = std::max((*phi_u_data)(idx), 0.0);
                            (*phi_a_data)(idx) = std::max((*phi_a_data)(idx), 0.0);
                            (*bond_data)(idx) = std::max((*bond_data)(idx), 0.0);
                        }
                    }
                }
            };
        // Register clamping callback
        sb_adv_diff_integrator->registerIntegrateHierarchyCallback(clamp_cell_values, static_cast<void*>(&clamp_vars));

        EquationSystems* eq_sys = fe_data_manager->getEquationSystems();
        std::unique_ptr<ExodusII_IO> exodus_io(uses_exodus ? new ExodusII_IO(mesh) : nullptr);
        //        std::unique_ptr<ExodusII_IO> bdry_io(new ExodusII_IO(*meshes[1]));
        bool from_restart = RestartManager::getManager()->isFromRestart();
        if (exodus_io) exodus_io->append(from_restart);

        // Initialize hierarchy configuration and data on all patches.
        ib_method_ops->initializeFEData();
        general_mesh_mapping->initializeFEData();
        time_integrator->initializePatchHierarchy(patch_hierarchy, gridding_algorithm);
        bdry_mesh_mapping->getMeshPartitioner()->setPatchHierarchy(patch_hierarchy);
        bdry_mesh_mapping->initializeFEData();
        bdry_mesh_mapping->setInitialConditions();

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
                    time_integrator->setupPlotData();
                    visit_data_writer->writePlotData(patch_hierarchy, iteration_num, loop_time);
                }
                if (uses_exodus)
                {
                    exodus_io->write_timestep(exodus_filename, *eq_sys, viz_dump_iteration_num, loop_time);
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
                bdry_mesh_mapping->writeFEDataToRestartFile(restart_dump_dirname, iteration_num);
                general_mesh_mapping->writeFEDataToRestartFile(restart_dump_dirname, iteration_num);
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
generateMesh(Mesh& mesh, const unsigned int N, double a, double L)
{
    double t0 = 0.0;
    double ell_arc = 2.0 * a * boost::math::ellint_2(1.0 - 1.0 / (a * a));
    double tmax = 2.0 * L + ell_arc;
    std::vector<VectorNd> pts;
    pts.reserve(N + 1);
    for (unsigned int n = 0; n <= N; ++n)
    {
        double t = t0 + (tmax - t0) * n / static_cast<double>(N);
        if (t < L)
        {
            // Initial straight section
            pts.push_back(VectorNd(t - L - 1.0, 0.0));
        }
        else if (t < (L + ell_arc))
        {
            // On the ellipse
            std::pair<double, int> s_iter_pair = find_s_from_t(L, a, t);
            pout << "Took " << s_iter_pair.second << " iterations\n";
            double s = s_iter_pair.first;
            pts.push_back(VectorNd(std::cos(M_PI * (s - 2.0)), a * std::sin(M_PI * (s - 2.0))));
        }
        else
        {
            // On back straight section
            double s = (t + L - ell_arc) / L;
            pts.push_back(VectorNd(L * (s - 2.0) + 1.0, 0.0));
        }
    }
    pout << "pts size: " << pts.size() << "\n";

    mesh.clear();
    mesh.set_mesh_dimension(NDIM - 1);
    mesh.set_spatial_dimension(NDIM);
    mesh.reserve_nodes(N);
    mesh.reserve_elem(N - 1);
    unsigned int node_id = 0;
    for (unsigned int n = 0; n < pts.size(); ++n)
        mesh.add_point(libMesh::Point(pts[n](0), pts[n](1), 0.0), node_id++, 0);

    for (unsigned int n = 0; n < pts.size() - 1; ++n)
    {
        Elem* elem = mesh.add_elem(Elem::build_with_id(EDGE2, n));
        elem->set_node(0) = mesh.node_ptr(n);
        elem->set_node(1) = mesh.node_ptr(n + 1);
    }

    mesh.prepare_for_use(false);
    mesh.print_info(pout);
}

std::pair<double, int>
find_s_from_t(double L, double a, double t)
{
    auto fp = [a](const double x) -> double {
        return -a * M_PI * std::sqrt(1.0 - (1.0 / (a * a)) * std::sin(M_PI * x));
    };
    auto f = [t, a, L](const double x) -> double {
        return t - L -
               a * (-2.0 * boost::math::ellint_2(1.0 - 1.0 / (a * a)) +
                    boost::math::ellint_2(1.0 - 1.0 / (a * a), M_PI * x));
    };

    double s = t;
    // Do Newton's iterations
    for (int i = 0; i < 100; ++i)
    {
        double snew = s - f(s) / fp(s);
        double res = std::abs(snew - s);
        s = snew;
        if (res < 1.0e-8) return std::make_pair(s, i);
    }
    return std::make_pair(std::numeric_limits<double>::quiet_NaN(), -1);
}
