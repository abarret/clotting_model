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
#include <ADS/LSFromLevelSet.h>
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

#include <boost/math/tools/roots.hpp>

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

void compute_plot_quantities(const int drag_idx,
                             const int div_sig_idx,
                             const int beta_idx,
                             const int phib_delete_idx,
                             const int pi_idx,
                             const int grad_pi_idx,
                             Pointer<PatchHierarchy<NDIM>> hierarchy,
                             Pointer<INSHierarchyIntegrator> ins_integrator,
                             Pointer<CFINSForcing> cf_forcing,
                             Pointer<CellVariable<NDIM, double>> bond_var,
                             Pointer<AdvDiffHierarchyIntegrator> bond_integrator,
                             const int ub_idx,
                             const int phib_idx,
                             const BoundClotParams& clot_params,
                             const BetaFcn& beta_fcn);

class LSFillFcn : public CartGridFunction
{
    void setDataOnPatch(const int data_idx,
                        Pointer<hier::Variable<NDIM>> /*var*/,
                        Pointer<Patch<NDIM>> patch,
                        const double /*data_time*/,
                        const bool /*initial_time*/,
                        Pointer<PatchLevel<NDIM>> /*level*/)
    {
        Pointer<NodeData<NDIM, double>> n_data = patch->getPatchData(data_idx);
        Pointer<CellData<NDIM, double>> c_data = patch->getPatchData(data_idx);
        if (n_data) n_data->fillAll(-1.0);
        if (c_data) c_data->fillAll(-1.0);
    }

    bool isTimeDependent() const override
    {
        return false;
    }
};

class FillPhiA : public CartGridFunction
{
    void setDataOnPatch(const int data_idx,
                        Pointer<hier::Variable<NDIM>> /*var*/,
                        Pointer<Patch<NDIM>> patch,
                        const double /*data_time*/,
                        const bool /*initial_time*/,
                        Pointer<PatchLevel<NDIM>> /*level*/)
    {
        Pointer<CellData<NDIM, double>> data = patch->getPatchData(data_idx);
        data->fillAll(3.0);
    }

    bool isTimeDependent() const override
    {
        return false;
    }
};

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
        const std::string restart_dump_dirname = app_initializer->getRestartDumpDirectory();
        const std::string restart_read_dirname = app_initializer->getRestartReadDirectory();
        const int restart_restore_num = app_initializer->getRestartRestoreNumber();

        const bool dump_postproc_data = app_initializer->dumpPostProcessingData();
        const int postproc_data_dump_interval = app_initializer->getPostProcessingDataDumpInterval();
        const std::string postproc_data_dump_dirname = app_initializer->getPostProcessingDataDumpDirectory();

        const bool dump_timer_data = app_initializer->dumpTimerData();
        const int timer_dump_interval = app_initializer->getTimerDumpInterval();

        // Create finite element mesh
        Mesh solid_mesh(init.comm(), NDIM);
        Pointer<Database> mesh_db = app_initializer->getComponentDatabase("MeshDB");
        const double dx = input_db->getDouble("DX");
        // mfac is the mesh factor, how many Lagrangian nodes per Eulerian grid cell?
        const double mfac = mesh_db->getDouble("mfac");
        const double ds = mfac * dx;
        std::string elem_type = mesh_db->getString("elem_type");
        std::string elem_order = mesh_db->getString("elem_order");
        const double mesh_xlow = mesh_db->getDouble("xlow");
        const double mesh_xup = mesh_db->getDouble("xup");
        const double mesh_y = mesh_db->getDouble("y");
        MeshTools::Generation::build_line(
            solid_mesh, (mesh_xup - mesh_xlow) / ds, mesh_xlow, mesh_xup, Utility::string_to_enum<ElemType>(elem_type));
        MeshTools::Modification::translate(solid_mesh, 0.0, mesh_y + 0.5 * dx, 0.0);
        if (elem_order == "SECOND")
            solid_mesh.all_second_order(true);
        else
            solid_mesh.all_first_order();
        solid_mesh.prepare_for_use();
        solid_mesh.print_info();

        MeshBase& mesh = solid_mesh;

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
        Pointer<CartesianGridGeometry<NDIM>> grid_geometry = new CartesianGridGeometry<NDIM>(
            "CartesianGeometry", app_initializer->getComponentDatabase("CartesianGeometry"));
        Pointer<PatchHierarchy<NDIM>> patch_hierarchy = new PatchHierarchy<NDIM>("PatchHierarchy", grid_geometry);
        Pointer<StandardTagAndInitialize<NDIM>> error_detector =
            new StandardTagAndInitialize<NDIM>("StandardTagAndInitialize",
                                               ins_integrator,
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

        // Create boundary condition specification objects (when necessary).
        const IntVector<NDIM>& periodic_shift = grid_geometry->getPeriodicShift();
        std::vector<RobinBcCoefStrategy<NDIM>*> u_bc_coefs(NDIM);
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
                std::ostringstream bc_coefs_name_stream;
                bc_coefs_name_stream << "u_bc_coefs_" << d;
                const std::string bc_coefs_name = bc_coefs_name_stream.str();

                std::ostringstream bc_coefs_db_name_stream;
                bc_coefs_db_name_stream << "VelocityBcCoefs_" << d;
                const std::string bc_coefs_db_name = bc_coefs_db_name_stream.str();

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
            ins_integrator->registerVisItDataWriter(visit_data_writer);
        }

        // Create advected quantities
        Pointer<CellVariable<NDIM, double>> phi_b_var = new CellVariable<NDIM, double>("phi_b");
        Pointer<CellVariable<NDIM, double>> phi_a_var = new CellVariable<NDIM, double>("phi_a");
        Pointer<CellVariable<NDIM, double>> bond_var = new CellVariable<NDIM, double>("bond");
        Pointer<CellVariable<NDIM, double>> ub_var = new CellVariable<NDIM, double>("Bound Velocity", NDIM);

        // Pull out the database everything uses
        BoundClotParams clot_params(app_initializer->getComponentDatabase("ClotParams"));
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
        cohesion_relax->setBoundPlateletData(phi_b_var, adv_diff_integrator);
        cohesion_relax->setActivatedPlateletData(phi_a_var, adv_diff_integrator);
        cohesion_relax->setBondData(bond_var, adv_diff_integrator);
        cohesion_relax->setKernel(BSPLINE_4);
        cohesionStressForcing->registerRelaxationOperator(cohesion_relax);
        Pointer<CellVariable<NDIM, double>> sig_var = cohesionStressForcing->getVariable();

        // Set up bound velocity function and register it with the sb integrator
        Pointer<FaceVariable<NDIM, double>> ub_adv_var = new FaceVariable<NDIM, double>("bound_velocity");
        adv_diff_integrator->registerAdvectionVelocity(ub_adv_var);
        adv_diff_integrator->setAdvectionVelocityFunction(ub_adv_var, cell_to_face_fcn);

        // Create the bound velocity solver.
        Pointer<BoundVelocitySource> ub_src_fcn = new BoundVelocitySource("BoundVelocitySource", clot_params);
        ub_src_fcn->setBondData(bond_var, adv_diff_integrator);
        ub_src_fcn->setBoundPlateletData(phi_b_var, adv_diff_integrator);
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
        drag_force->setBoundPlateletData(phi_b_var, adv_diff_integrator);
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
        adv_diff_integrator->registerTransportedQuantity(phi_a_var);
        Pointer<CellVariable<NDIM, double>> phi_a_src_var = new CellVariable<NDIM, double>("phi_a_src");
        adv_diff_integrator->setAdvectionVelocity(phi_a_var, ins_integrator->getAdvectionVelocityVariable());
        adv_diff_integrator->setDiffusionCoefficient(phi_a_var, clot_params.phi_a_diff_coef);
        Pointer<BoundPlateletSource> phi_a_src_fcn = new BoundPlateletSource("UnactivatedPlatelets", clot_params);
        phi_a_src_fcn->setSign(true);
        phi_a_src_fcn->setKernel(BSPLINE_4);
        phi_a_src_fcn->registerBetaFcn(beta_wrapper, static_cast<void*>(&betaFcn));
        phi_a_src_fcn->setActivatedPlateletData(phi_a_var, adv_diff_integrator);
        phi_a_src_fcn->setBondsData(bond_var, adv_diff_integrator);
        phi_a_src_fcn->setBoundPlateletData(phi_b_var, adv_diff_integrator);
        phi_a_src_fcn->setStressData(sig_var, adv_diff_integrator);
        adv_diff_integrator->registerSourceTerm(phi_a_src_var);
        adv_diff_integrator->setSourceTermFunction(phi_a_src_var, phi_a_src_fcn);
        adv_diff_integrator->setSourceTerm(phi_a_var, phi_a_src_var);
        Pointer<RobinBcCoefStrategy<NDIM>> phi_a_bcs;
        if (grid_geometry->getPeriodicShift().min() == 0)
            phi_a_bcs = new muParserRobinBcCoefs("UnactivatedPlateletsBcs",
                                                 app_initializer->getComponentDatabase("ActivatedPlateletBcs"),
                                                 grid_geometry);
        adv_diff_integrator->setPhysicalBcCoef(phi_a_var, phi_a_bcs.getPointer());

        // Bound platelets
        adv_diff_integrator->registerTransportedQuantity(phi_b_var);
        adv_diff_integrator->setInitialConditions(phi_b_var, phib_init);
        Pointer<CellVariable<NDIM, double>> phi_b_src_var = new CellVariable<NDIM, double>("phi_b_src");
        adv_diff_integrator->setAdvectionVelocity(phi_b_var, ub_adv_var);
        adv_diff_integrator->setDiffusionCoefficient(phi_b_var, clot_params.phi_b_diff_coef);
        Pointer<BoundPlateletSource> phi_b_src_fcn = new BoundPlateletSource("BoundPlatelets", clot_params);
        phi_b_src_fcn->setSign(false);
        phi_b_src_fcn->setKernel(BSPLINE_4);
        phi_b_src_fcn->registerBetaFcn(beta_wrapper, static_cast<void*>(&betaFcn));
        phi_b_src_fcn->setActivatedPlateletData(phi_a_var, adv_diff_integrator);
        phi_b_src_fcn->setBondsData(bond_var, adv_diff_integrator);
        phi_b_src_fcn->setBoundPlateletData(phi_b_var, adv_diff_integrator);
        phi_b_src_fcn->setStressData(sig_var, adv_diff_integrator);
        adv_diff_integrator->registerSourceTerm(phi_b_src_var);
        adv_diff_integrator->setSourceTermFunction(phi_b_src_var, phi_b_src_fcn);
        adv_diff_integrator->setSourceTerm(phi_b_var, phi_b_src_var);
        Pointer<RobinBcCoefStrategy<NDIM>> phi_b_bcs;
        if (grid_geometry->getPeriodicShift().min() == 0)
            phi_b_bcs = new muParserRobinBcCoefs(
                "ActivatedPlateletsBcs", app_initializer->getComponentDatabase("BoundPlateletBcs"), grid_geometry);
        adv_diff_integrator->setPhysicalBcCoef(phi_b_var, phi_b_bcs.getPointer());

        // Bonds
        adv_diff_integrator->registerTransportedQuantity(bond_var);
        adv_diff_integrator->setInitialConditions(bond_var, bond_init);
        Pointer<CellVariable<NDIM, double>> bond_src_var = new CellVariable<NDIM, double>("bond_src");
        adv_diff_integrator->setAdvectionVelocity(bond_var, ub_adv_var);
        adv_diff_integrator->setDiffusionCoefficient(bond_var, 0.0);
        Pointer<BondBoundSource> bond_src_fcn = new BondBoundSource("BondSource", clot_params);
        bond_src_fcn->registerBetaFcn(beta_wrapper, static_cast<void*>(&betaFcn));
        bond_src_fcn->setActivatedPlateletData(phi_a_var, adv_diff_integrator);
        bond_src_fcn->setBondData(bond_var, adv_diff_integrator);
        bond_src_fcn->setBoundPlateletData(phi_b_var, adv_diff_integrator);
        bond_src_fcn->setStressData(sig_var, adv_diff_integrator);
        bond_src_fcn->setKernel(BSPLINE_4);
        adv_diff_integrator->registerSourceTerm(bond_src_var);
        adv_diff_integrator->setSourceTermFunction(bond_src_var, bond_src_fcn);
        adv_diff_integrator->setSourceTerm(bond_var, bond_src_var);
        Pointer<RobinBcCoefStrategy<NDIM>> bond_bcs;
        if (grid_geometry->getPeriodicShift().min() == 0)
            bond_bcs = new muParserRobinBcCoefs(
                "BondVariableBcs", app_initializer->getComponentDatabase("BondVariableBcs"), grid_geometry);
        adv_diff_integrator->setPhysicalBcCoef(bond_var, bond_bcs.getPointer());

        // Wall sites
        const std::string wall_sites_name = "WallSites";
        auto wall_mesh_mapping =
            std::make_shared<WallSitesMeshMapping>("WallSitesMesh",
                                                   app_initializer->getComponentDatabase("BoundaryMesh"),
                                                   wall_sites_name,
                                                   &mesh,
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
        Pointer<LSFromLevelSet> ls_fcn = new LSFromLevelSet("LSFcn", patch_hierarchy);
        Pointer<LSFillFcn> ls_fill_fcn = new LSFillFcn();
        ls_fcn->registerLSFcn(ls_fill_fcn);
        sb_adv_diff_integrator->registerLevelSetVolFunction(ls_var, ls_fcn);
        const int w_idx = wall_mesh_mapping->getWallSitesPatchIndex();
        auto update_bdry_mesh = [](const double current_time,
                                   const double new_time,
                                   bool /*skip_synchronize_new_state_data*/,
                                   int /*num_cycles*/,
                                   void* ctx)
        {
            plog << "Updating boundary location\n";
            auto mesh_mapping = static_cast<BoundaryMeshMapping*>(ctx);
            mesh_mapping->updateBoundaryLocation(current_time, 0, true);
        };
        auto regrid_callback = [](SAMRAI::tbox::Pointer<SAMRAI::hier::BasePatchHierarchy<NDIM>> /*hierarchy*/,
                                  const double /*time*/,
                                  const bool /*initial_time*/,
                                  void* ctx)
        {
            plog << "Spreading wall sites\n";
            auto wall_mesh_mapping = static_cast<WallSitesMeshMapping*>(ctx);
            wall_mesh_mapping->spreadWallSites();
        };
        ins_integrator->registerPostprocessIntegrateHierarchyCallback(update_bdry_mesh,
                                                                      static_cast<void*>(wall_mesh_mapping.get()));
        ins_integrator->registerRegridHierarchyCallback(regrid_callback, static_cast<void*>(wall_mesh_mapping.get()));
        visit_data_writer->registerPlotQuantity("wall_sites", "SCALAR", w_idx);

        // define variable blob, which stores the patch hierarchy, variable context, and clamped variables.
        ClampVars clamp_vars = {
            patch_hierarchy,                      // patch hierarchy
            adv_diff_integrator->getNewContext(), // variable context
            phi_b_var,                            // bound platelets var
            phi_a_var,                            // activated platelets var
            bond_var                              // bonds var
        };
        // Clamp cell values to 0 in callback
        auto clamp_cell_values =
            [](const double /*current_time*/, const double /*new_time*/, int /*num_cycles*/, void* ctx)
        {
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
        adv_diff_integrator->registerIntegrateHierarchyCallback(clamp_cell_values, static_cast<void*>(&clamp_vars));

        EquationSystems* wall_eq_sys = wall_mesh_mapping->getMeshPartitioner()->getEquationSystems();
        std::unique_ptr<ExodusII_IO> wall_io(uses_exodus ? new ExodusII_IO(*wall_mesh_mapping->getBoundaryMesh()) :
                                                           nullptr);
        bool from_restart = RestartManager::getManager()->isFromRestart();
        if (wall_io) wall_io->append(from_restart);

        // Create extra drawing variables
        auto var_db = VariableDatabase<NDIM>::getDatabase();
        Pointer<CellVariable<NDIM, double>> drag_var = new CellVariable<NDIM, double>("DragForce", NDIM),
                                            div_sig_var = new CellVariable<NDIM, double>("DivSig", NDIM),
                                            beta_var = new CellVariable<NDIM, double>("Beta"),
                                            sig_src_var = new CellVariable<NDIM, double>("SigSrcDraw", 3),
                                            phib_delete_var = new CellVariable<NDIM, double>("phib_delete_draw"),
                                            grad_pi_var = new CellVariable<NDIM, double>("Grad Pi", NDIM),
                                            pi_var = new CellVariable<NDIM, double>("Pi");
        const int drag_idx = var_db->registerVariableAndContext(drag_var, var_db->getContext("Draw"));
        const int div_sig_idx = var_db->registerVariableAndContext(div_sig_var, var_db->getContext("Draw"));
        const int beta_idx = var_db->registerVariableAndContext(beta_var, var_db->getContext("Draw"));
        const int sig_src_idx = var_db->registerVariableAndContext(sig_src_var, var_db->getContext("Draw"));
        const int phib_delete_idx = var_db->registerVariableAndContext(phib_delete_var, var_db->getContext("Draw"));
        const int grad_pi_idx = var_db->registerVariableAndContext(grad_pi_var, var_db->getContext("Draw"));
        const int pi_idx = var_db->registerVariableAndContext(pi_var, var_db->getContext("Pi"), IntVector<NDIM>(1));

        visit_data_writer->registerPlotQuantity("DragForce", "VECTOR", drag_idx);
        visit_data_writer->registerPlotQuantity("DivSig", "VECTOR", div_sig_idx);
        visit_data_writer->registerPlotQuantity("Beta", "SCALAR", beta_idx);
        visit_data_writer->registerPlotQuantity("sig_src_0", "SCALAR", sig_src_idx, 0);
        visit_data_writer->registerPlotQuantity("sig_src_1", "SCALAR", sig_src_idx, 1);
        visit_data_writer->registerPlotQuantity("sig_src_2", "SCALAR", sig_src_idx, 2);
        visit_data_writer->registerPlotQuantity("phib_delete", "SCALAR", phib_delete_idx);
        visit_data_writer->registerPlotQuantity("grad_pi", "VECTOR", grad_pi_idx);
        visit_data_writer->registerPlotQuantity("pi", "SCALAR", pi_idx);

        // If we are from restart, fill phi_a
        const bool fill_phi_a = true; // && from_restart && input_db->getBool("FILL_PHI_A") ;
        Pointer<FillPhiA> fill_phi_a_fcn = new FillPhiA();

        // Initialize hierarchy configuration and data on all patches.
        ins_integrator->initializePatchHierarchy(patch_hierarchy, gridding_algorithm);
        wall_mesh_mapping->getMeshPartitioner()->setPatchHierarchy(patch_hierarchy);
        wall_mesh_mapping->initializeFEData(patch_hierarchy);
        sb_coupling_manager->fillInitialConditions();
        wall_mesh_mapping->spreadWallSites(w_idx);

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

        int iteration_num = ins_integrator->getIntegratorStep();
        double loop_time = ins_integrator->getIntegratorTime();
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

        if (fill_phi_a)
        {
            auto var_db = VariableDatabase<NDIM>::getDatabase();
            const int phia_idx =
                var_db->mapVariableAndContextToIndex(phi_a_var, adv_diff_integrator->getCurrentContext());
            fill_phi_a_fcn->setDataOnPatchHierarchy(
                phia_idx, phi_a_var, patch_hierarchy, loop_time, false, 0, patch_hierarchy->getFinestLevelNumber());
        }

        // Main time step loop.
        double loop_time_end = ins_integrator->getEndTime();
        double dt = 0.0;
        while (!MathUtilities<double>::equalEps(loop_time, loop_time_end) && ins_integrator->stepsRemaining())
        {
            if (dump_viz_data &&
                (MathUtilities<double>::equalEps(loop_time, next_viz_dump_time) || loop_time >= next_viz_dump_time))
            {
                pout << "\nWriting visualization files...\n\n";
                if (uses_visit)
                {
                    const int coarsest_ln = 0;
                    const int finest_ln = patch_hierarchy->getFinestLevelNumber();
                    const int phib_idx =
                        var_db->mapVariableAndContextToIndex(phi_b_var, adv_diff_integrator->getCurrentContext());
                    const int ub_idx =
                        var_db->mapVariableAndContextToIndex(ub_var, adv_diff_integrator->getCurrentContext());
                    const int phia_idx =
                        var_db->mapVariableAndContextToIndex(phi_a_var, adv_diff_integrator->getCurrentContext());
                    adv_diff_integrator->allocatePatchData(drag_idx, loop_time, coarsest_ln, finest_ln);
                    adv_diff_integrator->allocatePatchData(div_sig_idx, loop_time, coarsest_ln, finest_ln);
                    adv_diff_integrator->allocatePatchData(beta_idx, loop_time, coarsest_ln, finest_ln);
                    adv_diff_integrator->allocatePatchData(sig_src_idx, loop_time, coarsest_ln, finest_ln);
                    adv_diff_integrator->allocatePatchData(phib_delete_idx, loop_time, coarsest_ln, finest_ln);
                    adv_diff_integrator->allocatePatchData(grad_pi_idx, loop_time, coarsest_ln, finest_ln);
                    adv_diff_integrator->allocatePatchData(pi_idx, loop_time, coarsest_ln, finest_ln);
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
                                            phib_delete_idx,
                                            pi_idx,
                                            grad_pi_idx,
                                            patch_hierarchy,
                                            ins_integrator,
                                            cohesionStressForcing,
                                            bond_var,
                                            adv_diff_integrator,
                                            ub_idx,
                                            phib_idx,
                                            clot_params,
                                            betaFcn);
                    ins_integrator->setupPlotData();
                    visit_data_writer->writePlotData(patch_hierarchy, iteration_num, loop_time);
                    adv_diff_integrator->deallocatePatchData(drag_idx, coarsest_ln, finest_ln);
                    adv_diff_integrator->deallocatePatchData(div_sig_idx, coarsest_ln, finest_ln);
                    adv_diff_integrator->deallocatePatchData(beta_idx, coarsest_ln, finest_ln);
                    adv_diff_integrator->deallocatePatchData(sig_src_idx, coarsest_ln, finest_ln);
                    adv_diff_integrator->deallocatePatchData(phib_delete_idx, coarsest_ln, finest_ln);
                    adv_diff_integrator->deallocatePatchData(pi_idx, coarsest_ln, finest_ln);
                    adv_diff_integrator->deallocatePatchData(grad_pi_idx, coarsest_ln, finest_ln);
                }
                if (uses_exodus)
                {
                    wall_io->write_timestep(wall_mesh_filename, *wall_eq_sys, viz_dump_iteration_num, loop_time);
                }

                next_viz_dump_time += viz_dump_time_interval;
                viz_dump_iteration_num += 1;
            }
            iteration_num = ins_integrator->getIntegratorStep();
            loop_time = ins_integrator->getIntegratorTime();

            pout << "\n";
            pout << "+++++++++++++++++++++++++++++++++++++++++++++++++++\n";
            pout << "At beginning of timestep # " << iteration_num << "\n";
            pout << "Simulation time is " << loop_time << "\n";

            dt = ins_integrator->getMaximumTimeStepSize();
            ins_integrator->advanceHierarchy(dt);
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
            const bool last_step = !ins_integrator->stepsRemaining();
            if (dump_restart_data && (iteration_num % restart_dump_interval == 0 || last_step))
            {
                pout << "\nWriting restart files...\n\n";
                RestartManager::getManager()->writeRestartFile(restart_dump_dirname, iteration_num);
                wall_mesh_mapping->writeFEDataToRestartFile(restart_dump_dirname, iteration_num);
            }
            if (dump_timer_data && (iteration_num % timer_dump_interval == 0 || last_step))
            {
                pout << "\nWriting timer data...\n\n";
                TimerManager::getManager()->print(plog);
            }
        }

        for (unsigned int d = 0; d < NDIM; ++d)
        {
            delete ub_bcs[d];
            delete u_bc_coefs[d];
        }
    } // cleanup dynamically allocated objects prior to shutdown
} // main

void
compute_plot_quantities(const int drag_idx,
                        const int div_sig_idx,
                        const int beta_idx,
                        const int phib_delete_idx,
                        const int pi_idx,
                        const int grad_pi_idx,
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

    // Need to fill in pi_idx first, so we can fill in ghost cells for it.
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        Pointer<PatchLevel<NDIM>> level = hierarchy->getPatchLevel(ln);
        for (PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            Pointer<Patch<NDIM>> patch = level->getPatch(p());
            Pointer<CellData<NDIM, double>> pi_data = patch->getPatchData(pi_idx);
            Pointer<CellData<NDIM, double>> phib_data = patch->getPatchData(phib_idx);
            for (CellIterator<NDIM> ci(patch->getBox()); ci; ci++)
            {
                const CellIndex<NDIM>& idx = ci();
                const double pi_fac = clot_params.pi_fac;
                const double th_b = clot_params.vol_pl * (*phib_data)(idx);
                (*pi_data)(idx) = pi_fcn(th_b, pi_fac);
            }
        }
    }

    // Need to fill in ghost data for divergence
    const int sig_scr_idx =
        var_db->mapVariableAndContextToIndex(cf_forcing->getVariable(), adv_diff_integrator->getScratchContext());
    const bool sig_scr_is_allocated = adv_diff_integrator->isAllocatedPatchData(sig_scr_idx, coarsest_ln, finest_ln);
    if (!sig_scr_is_allocated) adv_diff_integrator->allocatePatchData(sig_scr_idx, 0.0, coarsest_ln, finest_ln);
    const bool bond_scr_is_allocated = bond_integrator->isAllocatedPatchData(bond_scr_idx, coarsest_ln, finest_ln);
    if (!bond_scr_is_allocated) bond_integrator->allocatePatchData(bond_scr_idx, 0.0, coarsest_ln, finest_ln);
    using ITC = HierarchyGhostCellInterpolation::InterpolationTransactionComponent;
    std::vector<ITC> ghost_cell_comp(3);
    ghost_cell_comp[0] =
        ITC(sig_scr_idx, sig_idx, "CONSERVATIVE_LINEAR_REFINE", false, "NONE", "LINEAR", false, nullptr);
    ghost_cell_comp[1] =
        ITC(bond_scr_idx, bond_idx, "CONSERVATIVE_LINEAR_REFINE", false, "NONE", "LINEAR", false, nullptr);
    ghost_cell_comp[2] = ITC(pi_idx, "CONSERVATIVE_LINEAR_REFINE", false, "NONE", "LINEAR", false, nullptr);
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
            Pointer<CellData<NDIM, double>> phib_delete_data = patch->getPatchData(phib_delete_idx);

            Pointer<CellData<NDIM, double>> div_sig_data = patch->getPatchData(div_sig_idx);
            Pointer<CellData<NDIM, double>> sig0_data = patch->getPatchData(sig_scr_idx);

            Pointer<CellData<NDIM, double>> pi_data = patch->getPatchData(pi_idx);
            Pointer<CellData<NDIM, double>> grad_pi_data = patch->getPatchData(grad_pi_idx);

            Pointer<CellData<NDIM, double>> sig_data =
                convertToStress(*sig0_data, *bond_data, patch, clot_params.S0, clot_params.R0, true);
            // Pointer<CellData<NDIM, double>> sig_data = convertToStressTraceless(*sig0_data, patch, true);

            Pointer<CartesianPatchGeometry<NDIM>> pgeom = patch->getPatchGeometry();
            const double* const dx = pgeom->getDx();

            for (CellIterator<NDIM> ci(patch->getBox()); ci; ci++)
            {
                const CellIndex<NDIM>& idx = ci();
                const double phi_b = (*phib_data)(idx);
                const double z = (*bond_data)(idx);

                // Compute drag force
                const double th = clot_params.vol_pl * phi_b;
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
                if (tr > 1.0e-12 && z > 1.0e-12) avg_y = std::sqrt(tr * 2.0 / (clot_params.S0 * z + 1.0e-12));
                const double beta = beta_fcn.beta(avg_y);
                (*beta_data)(idx) = beta;

                double P = 1.0;
                double nb_per_pl = 1.0e4 * z / (phi_b + 1.0e-12);
                if (nb_per_pl > 2.0)
                {
                    auto lambda_fcn = [nb_per_pl](const double lambda) -> std::pair<double, double>
                    {
                        double den = 1.0 - std::exp(-lambda);
                        double fcn = lambda / den;
                        double fcn_der = 1.0 / den - std::exp(-lambda) * lambda / (den * den);
                        return std::make_pair(fcn, fcn_der);
                    };
                    boost::uintmax_t max_iters = 100; // std::numeric_limits<boost::uintmax_t>::max();
                    double lambda = boost::math::tools::newton_raphson_iterate(
                        lambda_fcn, 1.0, 1.0, std::numeric_limits<double>::max(), 10, max_iters);
                    P = lambda / (std::exp(lambda) - 1.0);
                }
                // Note we change units here.
                (*phib_delete_data)(idx) = -beta * z * P * 1.0e4;

                for (int d = 0; d < NDIM; ++d)
                {
                    IntVector<NDIM> p1 = 0;
                    p1(d) = 1;
                    (*grad_pi_data)(idx, d) = -((*pi_data)(idx + p1) - (*pi_data)(idx - p1)) / (2.0 * dx[d]);
                }
            }
        }
    }

    if (!sig_scr_is_allocated) adv_diff_integrator->deallocatePatchData(sig_scr_idx, coarsest_ln, finest_ln);
    if (!bond_scr_is_allocated) bond_integrator->deallocatePatchData(bond_scr_idx, coarsest_ln, finest_ln);
}
