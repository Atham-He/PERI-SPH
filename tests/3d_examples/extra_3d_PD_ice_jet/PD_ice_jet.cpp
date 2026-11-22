/**
 * @file PD_ice_jet.cpp
 * @brief High-velocity water jet impacting a rotating ice plate, simulated
 * by the coupled WCSPH (fluid) and NOSB-PD (structure) model with
 * Drucker-Prager plasticity and bond breaking.
 * @details Ported from Icebreaking_PD-SPH (Haotian Shi, Luwen Zhang, SJTU),
 * case PD_waterColumn_rotation_BaoyuNi_plastic. This is the full-system
 * validation of the peridynamics module: FSI coupling, plasticity,
 * fracture and fluid-solid subcycling.
 * NOTE: the plate particles are generated directly on the lattice here,
 * while the reference case relaxes and reloads them.
 */
#include "sphinxsys.h"
#include "pd_fluid_time_step.h"
using namespace SPH;

// general parameters for geometry
Real LR = 0.0036;        // liquid column radius in XZ
Real LH = 12.0 * LR;     // liquid column height in Y
Real inner_circle_radius = LR;
int resolution(20);

Real plate_Y = 0.0055;   /**< thickness of the ice plate. */
Real plate_R = 0.0575;   /**< radius of the ice plate. */

Real resolution_ref = LR / 4.0;   // particle spacing
Real BW = resolution_ref * 3.0;   // boundary width

// for material properties of the fluid
Real rho0_f = 1000.0;
Real gravity_g = 9.81;
Real U_f = 73.81;
Real c_f = SMIN(10.0 * U_f, 1580.0);
Real mu_f = 0.0;
Real k_f = 0.0;
// material properties of the ice plate
Real rho0_s = 900.0;
Real poisson = 0.33;
Real Youngs_modulus = 6.2e9;
Real sigma_t0 = 2.2e6;
Real epsilon_rate = U_f / plate_Y;
Real sigma_c0 = 10.976e6 * pow(epsilon_rate, 0.093783);
// hardening parameters for DP
Real alpha = (sigma_c0 - sigma_t0) / (sigma_c0 + sigma_t0) / sqrt(3.0);
Real flow_stress = (2.0 * sigma_c0 * sigma_t0) / (sigma_c0 + sigma_t0) / sqrt(3.0);
Real isohardening_modulus_H = 6.89e6;
Real kinhardening_modulus_H = 0.0;
// damage parameters
Real max_tension_stress = 2.2e6;
Real max_shear_stress = 6.9e6;

// define the water block shape
class WaterBlock : public ComplexShape
{
  public:
    explicit WaterBlock(const std::string &shape_name) : ComplexShape(shape_name)
    {
        Vecd translation_column(0.0, 1.5 * LH - 0.006, 0.0);
        add<TriangleMeshShapeCylinder>(Vec3d(0.0, 1.0, 0.0), inner_circle_radius,
                                       0.5 * LH, resolution, translation_column);
    }
};
// define the plate shape
class PlateShape : public ComplexShape
{
  public:
    explicit PlateShape(const std::string &shape_name) : ComplexShape(shape_name)
    {
        Vecd offset_plate(0.0, LH - 0.5 * plate_Y - 0.5 * BW, 0.0);
        add<TriangleMeshShapeCylinder>(Vec3d(0.0, 1.0, 0.0), plate_R,
                                       0.5 * plate_Y, resolution, offset_plate);
    }
};
/** rotate body initial configuration by angles (theta_X, theta_Y, theta_Z). */
class BodyInitialRotation
{
  public:
    explicit BodyInitialRotation(Real theta_X, Real theta_Y, Real theta_Z)
    {
        Eigen::AngleAxisd rot_x(theta_X, Vecd::UnitX());
        Eigen::AngleAxisd rot_y(theta_Y, Vecd::UnitY());
        Eigen::AngleAxisd rot_z(theta_Z, Vecd::UnitZ());
        Q_ = (rot_z * rot_y * rot_x).matrix();
    };
    Matd Q_;
};
class WaterInitialRotation
    : public fluid_dynamics::FluidInitialCondition, public BodyInitialRotation
{
  public:
    explicit WaterInitialRotation(SPHBody &sph_body, Real theta_X, Real theta_Y, Real theta_Z)
        : fluid_dynamics::FluidInitialCondition(sph_body),
          BodyInitialRotation(theta_X, theta_Y, theta_Z) {};

    void update(size_t index_i, Real dt)
    {
        pos_[index_i] = Q_ * pos_[index_i];
    };
};
class PlateInitialRotation
    : public solid_dynamics::ElasticDynamicsInitialCondition, public BodyInitialRotation
{
  public:
    explicit PlateInitialRotation(SPHBody &sph_body, Real theta_X, Real theta_Y, Real theta_Z)
        : solid_dynamics::ElasticDynamicsInitialCondition(sph_body),
          BodyInitialRotation(theta_X, theta_Y, theta_Z) {};

    void update(size_t index_i, Real dt)
    {
        pos_[index_i] = Q_ * pos_[index_i];
    };
};
/** initial downward jet velocity. */
class WaterInitialCondition
    : public fluid_dynamics::FluidInitialCondition
{
  public:
    explicit WaterInitialCondition(SPHBody &sph_body)
        : fluid_dynamics::FluidInitialCondition(sph_body) {};

    void update(size_t index_i, Real dt)
    {
        vel_[index_i][1] = -U_f;
    };
};

// the main program with commandline options
int main(int ac, char *av[]) try
{
    const peridynamics::PeriExecution execution(peridynamics::takePeriExecutionOption(ac,av));
    //----------------------------------------------------------------------
    //	Build up an SPHSystem.
    //----------------------------------------------------------------------
    BoundingBoxd system_domain_bounds(Vecd(-plate_R - BW, -BW, -plate_R - BW),
                                      Vecd(plate_R + BW, 2.0 * LH + BW, plate_R + BW));
    SPHSystem sph_system(system_domain_bounds, resolution_ref);
    sph_system.handleCommandlineOptions(ac, av);
    //----------------------------------------------------------------------
    //	Creating bodies with corresponding materials and particles.
    //----------------------------------------------------------------------
    FluidBody water_block(sph_system, makeShared<WaterBlock>("WaterBody"));
    water_block.defineMatterMaterial<WeaklyCompressibleFluid>(rho0_f, c_f);
    water_block.generateParticles<BaseParticles, Lattice>();

    PDBody plate(sph_system, makeShared<PlateShape>("PDBody"));
    plate.defineMatterMaterial<DruckerPragerPlasticityforPD>(
        rho0_s, Youngs_modulus, poisson, flow_stress, alpha,
        isohardening_modulus_H, kinhardening_modulus_H);
    plate.generateParticles<BaseParticles, Lattice>();
    plate.getBaseParticles().addVariableToWrite<Real>("Damage");
    plate.getBaseParticles().addVariableToWrite<Matd>("PlasticStrain");
    //----------------------------------------------------------------------
    //	Define body relation map.
    //----------------------------------------------------------------------
    InnerRelation water_block_inner(water_block);
    InnerRelation plate_inner_relation(plate);
    /** FSI contacts: the fluid kernel is always used for the neighbor search
     * (source kernel for fluid -> PD, target kernel for PD -> fluid). */
    ContactRelationSelectKernel water_plate_contact(water_block, {&plate}, true);
    ContactRelationSelectKernel plate_water_contact(plate, {&water_block}, false);
    ComplexRelation water_block_complex(water_block_inner, water_plate_contact);
    //----------------------------------------------------------------------
    //	Define the numerical methods used in the simulation.
    //----------------------------------------------------------------------
    Gravity gravity(Vec3d(0.0, -gravity_g, 0.0));
    SimpleDynamics<GravityForce<Gravity>> constant_gravity(water_block, gravity);
    SimpleDynamics<WaterInitialRotation> jet_pos_rotation(water_block, 0.0, Pi / 9.0, 0.0);
    SimpleDynamics<WaterInitialCondition> jet_vel_initialization(water_block);
    SimpleDynamics<PlateInitialRotation> plate_pos_rotation(plate, 0.0, Pi / 4.0, 0.0);

    Dynamics1Level<fluid_dynamics::Integration1stHalfWithWallRiemann> pressure_relaxation(water_block_inner, water_plate_contact);
    Dynamics1Level<fluid_dynamics::Integration2ndHalfWithWallRiemann> density_relaxation(water_block_inner, water_plate_contact);
    InteractionWithUpdate<fluid_dynamics::DensitySummationComplexFreeSurface> update_density_by_summation(water_block_inner, water_plate_contact);
    auto fluid_cfl = execution.makeFluidCfl();
    fluid_dynamics::FluidAdvectionTimeStepExecution<fluid_dynamics::AdvectionTimeStep>
        get_fluid_advection_time_step_size(water_block, fluid_cfl.get(), U_f);
    fluid_dynamics::FluidAcousticTimeStepExecution get_fluid_time_step_size(water_block, fluid_cfl.get(), 0.06);
    std::cout << "Fluid CFL backend " << (fluid_cfl ? "cuda_fp64" : "cpu_fp64") << std::endl;
    /** velocity divergence on the fluid side (after the 2nd half step).
     * Note: must be constructed before the acoustic damper,
     * which retrieves the registered VelocityDivergence variable. */
    peridynamics::FluidNeighborConfiguration fluid_configuration;
    fluid_dynamics::PDFluidAcousticWorkspace acoustic(water_block_inner, water_plate_contact,
                                                       execution.makeFluidAcoustic(), 1.0, &fluid_configuration);
    fluid_dynamics::PDFluidVelocityDivergence velocity_divergence(acoustic);
    fluid_dynamics::PDFluidAcousticDamping acoustic_damper(acoustic);
    acoustic_damper.setCoeffAcousticDamper(rho0_f, c_f, water_block.getSPHAdaptation().getKernel()->CutOffRadius());
    //----------------------------------------------------------------------
    //	Algorithms of FSI.
    //----------------------------------------------------------------------
    SimpleDynamics<NormalDirectionFromBodyShape> plate_normal_direction(plate);
    solid_dynamics::PDFluidPressureCoupling<> fluid_pressure_force_on_plate(plate_water_contact, execution, &fluid_configuration);
    std::cout<<"PD fluid pressure backend "<<(fluid_pressure_force_on_plate.onDevice()?"cuda_fp64":"cpu_fp64")<<std::endl;
    fluid_pressure_force_on_plate.setCoeffAcousticDamper(rho0_f, c_f, water_block.getSPHAdaptation().getKernel()->CutOffRadius());
    solid_dynamics::AverageVelocityAndAcceleration average_velocity_and_acceleration(plate);
    /** zero-gravity baseline for PDAccelerationPrior: the fluid pressure force
     * now accumulates onto it, so the baseline must be reset every fluid step. */
    SharedPtr<Gravity> zero_gravity_ptr = makeShared<Gravity>(Vec3d::Zero());
    SimpleDynamics<solid_dynamics::PDTimeStepInitialization> initialize_a_pd_step(plate, zero_gravity_ptr);
    //----------------------------------------------------------------------
    //	Algorithms of NOSB-PD dynamics.
    //----------------------------------------------------------------------
    auto pd_time_integration = execution.makePositionVerlet();
    solid_dynamics::NosbPDStepWithExecution<solid_dynamics::NosbPDFirstStep> nosbpd_first_step(plate, pd_time_integration.get());
    std::cout << "NOSB-PD position Verlet backend " << (pd_time_integration ? "cuda_fp64" : "cpu_fp64") << std::endl;
    solid_dynamics::NosbPDBondWorkspace pd_bonds(plate_inner_relation, execution.makeNosbPDBonds());
    solid_dynamics::NosbPDShapeWithExecution plate_shape_matrix(plate_inner_relation, &pd_bonds);
    std::cout << "NOSB-PD reference shape backend " << (plate_shape_matrix.onDevice() ? "cuda_fp64" : "cpu_fp64") << std::endl;
    const peridynamics::NosbPDBondFailureParameters failure{
        peridynamics::NosbPDBondFailureLaw::StressShear, max_tension_stress, max_shear_stress, Dimensions};
    solid_dynamics::NosbPDBondFractureWithExecution<solid_dynamics::BondBreakBySigma1andSigma3>
        check_bond_live(plate_inner_relation, &pd_bonds, failure, max_tension_stress, max_shear_stress);
    std::cout << "NOSB-PD bond failure backend " << (pd_bonds.screensBondFailure() ? "cuda_fp64" : "cpu_fp64") << std::endl;
    auto pd_plastic = execution.makeNosbPDPlastic();
    std::cout << "NOSB-PD drucker_prager stress backend " << (pd_plastic ? "cuda_fp64" : "cpu_fp64") << std::endl;
    solid_dynamics::NosbPDPlasticWithExecution nosbpd_second_step_plastic(
        plate_inner_relation, &pd_bonds, std::move(pd_plastic));
    std::cout << "NOSB-PD internal force backend " << (pd_bonds.onDevice() ? "cuda_fp64" : "cpu_fp64") << std::endl;
    solid_dynamics::NosbPDStabilizationWithExecution hourglass_control(pd_bonds);
    solid_dynamics::NosbPDDampingWithExecution numerical_damping(pd_bonds);
    std::cout << "NOSB-PD stabilization backend " << (pd_bonds.onDevice() ? "cuda_fp64" : "cpu_fp64")
              << "; damping " << (numerical_damping.onDevice() ? "cuda_fp64" : "cpu_fp64") << std::endl;
    std::cout << "NOSB-PD paired corrections " << (hourglass_control.pairsWithDamping(numerical_damping) ? "cuda_fp64" : "separate") << std::endl;
    solid_dynamics::NosbPDInternalForceWithExecution nosbpd_third_step(pd_bonds);
    solid_dynamics::NosbPDStepWithExecution<solid_dynamics::NosbPDFourthStep> nosbpd_fourth_step(plate, pd_time_integration.get());
    /** normal direction updated from the PD deformation gradient;
     * must be constructed after the NOSB-PD dynamics which register
     * the DeformationGradient variable. */
    SimpleDynamics<solid_dynamics::UpdateElasticNormalDirection> plate_update_normal(plate);
    //----------------------------------------------------------------------
    //	Output.
    //----------------------------------------------------------------------
    BodyStatesRecordingToVtp write_body_states(sph_system);
    ParticleSorting particle_sorting(water_block);
    //----------------------------------------------------------------------
    //	Prepare the simulation with cell linked list, configuration
    //	and case specified initial condition if necessary.
    //----------------------------------------------------------------------
    sph_system.initializeSystemCellLinkedLists();
    sph_system.initializeSystemConfigurations();
    constant_gravity.exec();
    jet_pos_rotation.exec();
    plate_pos_rotation.exec();
    jet_vel_initialization.exec();
    plate_normal_direction.exec();
    plate_shape_matrix.exec();
    //----------------------------------------------------------------------
    //	Setup for time-stepping control.
    //----------------------------------------------------------------------
    Real &physical_time = *sph_system.getSystemVariableDataByName<Real>("PhysicalTime");
    size_t number_of_iterations = 0;
    int screen_output_interval = 50;
    Real end_time = 0.6e-3;
    Real output_interval = end_time / 100.0;
    Real dt = 0.0;
    Real dt_s_0 = 4.0e-8;
    //----------------------------------------------------------------------
    //	Statistics for CPU time.
    //----------------------------------------------------------------------
    TickCount t1 = TickCount::now();
    TimeInterval interval;
    //----------------------------------------------------------------------
    //	First output before the main loop.
    //----------------------------------------------------------------------
    write_body_states.writeToFile(0);
    //----------------------------------------------------------------------
    //	Main loop starts here.
    //----------------------------------------------------------------------
    while (physical_time < end_time)
    {
        Real integration_time = 0.0;
        while (integration_time < output_interval)
        {
            Real Dt = get_fluid_advection_time_step_size.exec();
            update_density_by_summation.exec();
            plate_update_normal.exec();

            Real relaxation_time = 0.0;
            for (const auto configuration = fluid_configuration.reuse(); relaxation_time < Dt;)
            {
                /** fluid relaxation and force computation. */
                pressure_relaxation.exec(dt);
                acoustic_damper.exec(dt);
                initialize_a_pd_step.exec(dt);
                fluid_pressure_force_on_plate.exec(dt);
                density_relaxation.exec(dt);
                velocity_divergence.exec(dt);
                /** NOSB-PD subcycling. */
                Real dt_s_sum = 0.0;
                average_velocity_and_acceleration.initialize_displacement_.exec();
                while (dt_s_sum < dt)
                {
                    Real dt_s = dt_s_0;
                    if (dt - dt_s_sum < dt_s) dt_s = dt - dt_s_sum;

                    nosbpd_first_step.exec(dt_s);
                    check_bond_live.exec(dt_s);
                    nosbpd_second_step_plastic.exec(dt_s);
                    hourglass_control.execWithDamping(numerical_damping, dt_s);
                    nosbpd_third_step.exec(dt_s);
                    nosbpd_fourth_step.exec(dt_s);

                    dt_s_sum += dt_s;
                }
                average_velocity_and_acceleration.update_averages_.exec(dt);

                dt = get_fluid_time_step_size.exec();
                relaxation_time += dt;
                integration_time += dt;
                physical_time += dt;
            }

            if (number_of_iterations % screen_output_interval == 0)
            {
                std::cout << std::fixed << std::setprecision(9) << "N=" << number_of_iterations
                          << " Time: " << physical_time
                          << "	Dt = " << Dt << "	dt = " << dt << "\n";
            }
            number_of_iterations++;

            if (number_of_iterations % 100 == 0 && number_of_iterations != 1)
            {
                particle_sorting.exec();
            }
            water_block.updateCellLinkedList();
            plate.updateCellLinkedList();
            water_block_complex.updateConfiguration();
            plate_water_contact.updateConfiguration();
        }

        TickCount t2 = TickCount::now();
        write_body_states.writeToFile();
        TickCount t3 = TickCount::now();
        interval += t3 - t2;
    }
    TickCount t4 = TickCount::now();

    TimeInterval tt;
    tt = t4 - t1 - interval;
    std::cout << "Total wall time for computation: " << tt.seconds() << " seconds." << std::endl;

    return 0;
}
catch (const std::exception &error)
{
    std::cerr<<error.what()<<std::endl;
    return 2;
}
