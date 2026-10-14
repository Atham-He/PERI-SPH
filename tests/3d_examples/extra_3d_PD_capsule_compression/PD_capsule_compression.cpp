/*-----------------------------------------------------------------------------*
 *       SPHinXsys: 3D single capsule compression (PD shell + SPH core)        *
 *-----------------------------------------------------------------------------*
 * A liquid-filled capsule is compressed between two rigid plates.             *
 * The capsule shell is discretized by NOSB-PD particles (brittle or           *
 * visco-plastic, bond breaking), the liquid core by WCSPH particles,          *
 * and the plates are rigid bodies driven by prescribed motion.                *
 * Generic placeholder parameters are used; replace with measured values       *
 * for production runs.                                                        *
 *-----------------------------------------------------------------------------*/
#include "sphinxsys.h"
using namespace SPH;

//----------------------------------------------------------------------
//	Geometry parameters (generic: diameter 1 mm, shell ratio 2%).
//----------------------------------------------------------------------
Real capsule_radius_out = 0.5e-3;                                ///< capsule outer radius [m]
Real shell_thickness = 0.02e-3;                                  ///< shell thickness (2% of diameter) [m]
Real capsule_radius_in = capsule_radius_out - shell_thickness;   ///< capsule inner radius [m]
Real resolution_shell = shell_thickness / 3.0;                   ///< 3 particles through the thickness
Real resolution_ref = resolution_shell;                          ///< system reference spacing
Real plate_radius = 0.75e-3;                                     ///< rigid plate radius [m]
Real plate_thickness = 0.06e-3;                                  ///< rigid plate thickness [m]
Real initial_gap = 1.5 * resolution_shell;                       ///< initial plate-capsule gap, beyond the contact activation range [m]
//----------------------------------------------------------------------
//	Layout along the compression (y) axis.
//----------------------------------------------------------------------
Real lower_plate_bottom = 0.0;
Real lower_plate_top = lower_plate_bottom + plate_thickness;
Real capsule_center_y = lower_plate_top + initial_gap + capsule_radius_out;
Real upper_plate_bottom = capsule_center_y + capsule_radius_out + initial_gap;
//----------------------------------------------------------------------
//	Material parameters of the shell (generic brittle polymer).
//----------------------------------------------------------------------
Real rho0_shell = 1000.0;      ///< density [kg/m^3]
Real youngs_shell = 1.0e7;     ///< Young's modulus (10 MPa) [Pa]
Real poisson_shell = 0.4;      ///< Poisson ratio
Real max_tension_stress = 1.0e5; ///< bond-breaking critical principal stress
                                 ///< (generic low value for functional fracture test; replace with measured strength) [Pa]
//----------------------------------------------------------------------
//	Shell constitutive route.
//	false: brittle HughesWingetSolid + principal-stress bond breaking
//	       (ethyl-cellulose type shell).
//	true:  alginate-type ViscoPlasticSolidforPD (Perzyna rate-dependent J2)
//	       + plastic-strain bond breaking.
//----------------------------------------------------------------------
bool shell_viscoplastic = true;        ///< constitutive route switch
Real yield_stress_shell = 5.0e5;       ///< alginate initial yield stress (generic) [Pa]
Real hardening_shell = 5.0e6;          ///< alginate isotropic hardening modulus
                                       ///< (calibration B3: 5x, post-yield stiffness for mid-range force) [Pa]
Real viscosity_shell = 5.0e2;          ///< Perzyna viscosity (calibration B3: mid-point;
                                       ///< friction now prevents squirt-out, so less plastic softening) [Pa s]
Real critical_plastic_strain = 1.0;    ///< bond-breaking critical equivalent plastic strain
                                       ///< (calibration B1: ductile squash like hydrated alginate
                                       ///<  experiments - no premature catastrophic puncture) [-]
Real plate_friction_mu = 0.3;          ///< smooth-Coulomb plate-shell friction coefficient
                                       ///< (B2: hold the soft shell between plates, no squirt-out) [-]
Real plate_friction_veps = 5.0e-3;     ///< friction regularization velocity [m/s]
//----------------------------------------------------------------------
//	Material parameters of the liquid core.
//----------------------------------------------------------------------
Real rho0_f = 935.0;           ///< reference density of the tung-oil core [kg/m^3]
Real gravity_g = 9.8;          ///< gravity [m/s^2]
Real compression_speed = 0.1;  ///< plate speed (debug value; still low Mach) [m/s]
Real U_ref = 0.5;              ///< reference (outflow) velocity scale [m/s]
Real c_f = 25.0;               ///< sound speed -> K_f = rho0_f*c_f^2 ~ 0.58 MPa [m/s]
                               ///< reduced from 50 m/s to compensate the acoustic-CFL
                               ///< tightening caused by the fluid refinement below
                               ///< (dt_acoustic ~ 0.6*h_f/c_f, h_f halves -> c_f halves)
Real mu_f = 2.5e-1;            ///< tung-oil dynamic viscosity (~0.25 Pa s at 25 C; the
                               ///< experimental core liquid is tung oil, confirmed by user) [Pa s]
//----------------------------------------------------------------------
//	Dry-compression control.
//	When false the liquid core is omitted entirely: no fluid particles,
//	no fluid dynamics and no fluid->shell pressure load. The shell is
//	then loaded only by the plates (contact + prescribed motion), which
//	decouples the shell/contact response from the fluid coupling.
//----------------------------------------------------------------------
bool with_liquid = true;       ///< true: liquid-filled capsule (production); false: dry control
//----------------------------------------------------------------------
//	Compression schedule.
//----------------------------------------------------------------------
Real compression_distance = 0.58e-3;                        ///< total plate stroke (58% diameter: cover experimental d/D range) [m]
Real end_time = compression_distance / compression_speed;   ///< 0.0058 s
//----------------------------------------------------------------------
//	Boundary width for the system domain.
//----------------------------------------------------------------------
Real BW = resolution_ref * 4.0;

//----------------------------------------------------------------------
//	Capsule shell: hollow sphere by boolean operations.
//----------------------------------------------------------------------
class CapsuleShellShape : public ComplexShape
{
  public:
    explicit CapsuleShellShape(const std::string &shape_name) : ComplexShape(shape_name)
    {
        Vecd center(0.0, capsule_center_y, 0.0);
        add<GeometricShapeBall>(center, capsule_radius_out);
        subtract<GeometricShapeBall>(center, capsule_radius_in);
    }
};
//----------------------------------------------------------------------
//	Rigid plate: brick centered at the plate mid-plane.
//----------------------------------------------------------------------
class PlateShape : public ComplexShape
{
  public:
    PlateShape(const std::string &shape_name, Real mid_plane_y) : ComplexShape(shape_name)
    {
        Vecd halfsize(plate_radius, 0.5 * plate_thickness, plate_radius);
        add<GeometricShapeBox>(Transform(Vecd(0.0, mid_plane_y, 0.0)), halfsize);
    }
};

//----------------------------------------------------------------------
//	Rigid-plate wall material.
//
//	The two fluid-wall paths must impose the SAME mirror wall velocity
//	2*v_wall - v_i on the fluid:
//	  (a) the standard SPHinXsys wall classes (InteractionWithWall ->
//	      Integration1stHalf/2ndHalfWithWallRiemann) obtain the wall velocity
//	      through the material's VIRTUAL AverageVelocity();
//	  (b) the PD wall class VelocityDivergenceWithWallforPD reads the LITERAL
//	      state variable "AverageVelocity".
//	The base Solid material returns "Velocity" from AverageVelocity(), and
//	PositionSolidBody resets "Velocity" to zero on every step.  Registering
//	"AverageVelocity" only for the PD class therefore makes the driven piston
//	look like a stationary wall to (a) but like a moving wall to (b): the two
//	wall treatments disagree at the coupling surface and pump net energy into
//	the fluid, which destabilises the run within the first ~45 steps.
//
//	This material closes that gap by returning "AverageVelocity" from the
//	virtual AverageVelocity() as well, so both paths read the same variable
//	holding the prescribed rigid-body velocity (zero for the stationary lower
//	plate, piston speed downwards for the driven upper plate).
//----------------------------------------------------------------------
class MovingWallSolid : public Solid
{
  public:
    MovingWallSolid(Real rho0, Real youngs, const Vecd &ave_velocity)
        : Solid(rho0, youngs), ave_velocity_(ave_velocity) {}

    Vecd *AverageVelocity(BaseParticles *base_particles) override
    {
        return base_particles->registerStateVariableData<Vecd>("AverageVelocity", ave_velocity_);
    }

  protected:
    Vecd ave_velocity_;
};

// the main program with commandline options
int main(int ac, char *av[])
{
    //----------------------------------------------------------------------
    //	Build up an SPHSystem.
    //	Domain enlarged to contain the outflow jet: fluid particles escaping
    //	the capsule must keep valid neighborhoods inside the domain, instead
    //	of crowding into boundary cells (which blows up the density summation
    //	and eventually NaN-poisons the whole particle field).
    //----------------------------------------------------------------------
    BoundingBoxd system_domain_bounds(
        Vecd(-3.0e-3, -4.0e-3, -3.0e-3),
        Vecd(3.0e-3, 6.0e-3, 3.0e-3));
    SPHSystem sph_system(system_domain_bounds, resolution_ref);
    sph_system.handleCommandlineOptions(ac, av);
    //----------------------------------------------------------------------
    //	Creating bodies with corresponding materials and particles.
    //----------------------------------------------------------------------
    FluidBody liquid_core(sph_system, makeShared<GeometricShapeBall>(
                                          Vecd(0.0, capsule_center_y, 0.0), capsule_radius_in, "LiquidCore"));
    /** RESOLUTION MISMATCH FIX: the fluid spacing must NOT equal the shell
     * thickness (which would make the shell only one fluid particle thick and
     * let unresolved fluid pressure penetrate it). 1/2 refinement (relative to
     * resolution_ref = resolution_shell) gives spacing = 1.5 * resolution_shell
     * = 1e-5 m = shell_thickness / 2, i.e. several fluid particles across the
     * shell thickness. The acoustic CFL is compensated by lowering c_f. */
    liquid_core.defineAdaptationRatios(1.3, 2.0 / 3.0); ///< spacing = 1.5 * resolution_shell = 1e-5
    liquid_core.defineMatterMaterial<WeaklyCompressibleFluid>(rho0_f, c_f);
    /** NOTE: the fluid body must always carry particles (SPHinXsys raises
     * "base_particles_ is nullptr" on an empty body). In the dry control the
     * particles are generated but never integrated and never coupled to the
     * shell, so they stay frozen and exert no load. */
    liquid_core.generateParticles<BaseParticles, Lattice>();
    std::cout << "liquid generated" << std::endl;

    PDBody capsule_shell(sph_system, makeShared<CapsuleShellShape>("CapsuleShell"));
    if (shell_viscoplastic)
        capsule_shell.defineMatterMaterial<ViscoPlasticSolidforPD>(
            rho0_shell, youngs_shell, poisson_shell, yield_stress_shell, viscosity_shell, hardening_shell);
    else
        capsule_shell.defineMatterMaterial<HughesWingetSolid>(rho0_shell, youngs_shell, poisson_shell);
    capsule_shell.generateParticles<BaseParticles, Lattice>();
    std::cout << "shell generated" << std::endl;
    capsule_shell.getBaseParticles().addVariableToWrite<Real>("Damage");

    SolidBody lower_plate(sph_system, makeShared<PlateShape>("LowerPlate", 0.5 * (lower_plate_bottom + lower_plate_top)));
    lower_plate.defineAdaptationRatios(1.3, 1.0 / 3.0);
    lower_plate.defineMatterMaterial<MovingWallSolid>(rho0_shell, youngs_shell, Vecd::Zero());
    lower_plate.generateParticles<BaseParticles, Lattice>();
    std::cout << "lower generated" << std::endl;

    SolidBody upper_plate(sph_system, makeShared<PlateShape>("UpperPlate", upper_plate_bottom + 0.5 * plate_thickness));
    upper_plate.defineAdaptationRatios(1.3, 1.0 / 3.0);
    upper_plate.defineMatterMaterial<MovingWallSolid>(rho0_shell, youngs_shell, Vecd(0.0, -compression_speed, 0.0));
    upper_plate.generateParticles<BaseParticles, Lattice>();
    std::cout << "upper generated" << std::endl;
    //----------------------------------------------------------------------
    //	Define body relation map.
    //	The fluid kernel is always used for fluid-involving neighbor search.
    //----------------------------------------------------------------------
    InnerRelation liquid_inner(liquid_core);
    InnerRelation shell_inner(capsule_shell);
    ContactRelationSelectKernel liquid_contact(liquid_core, {&capsule_shell, &lower_plate, &upper_plate}, true);
    ContactRelationSelectKernel shell_liquid_contact(capsule_shell, {&liquid_core}, false);
    ContactRelation shell_plates_contact(capsule_shell, {&lower_plate, &upper_plate});
    ContactRelation lower_plate_shell_contact(lower_plate, {&capsule_shell});
    ContactRelation upper_plate_shell_contact(upper_plate, {&capsule_shell});
    ComplexRelation liquid_complex(liquid_inner, liquid_contact);
    std::cout << "relations done" << std::endl;
    //----------------------------------------------------------------------
    //	Numerical methods: fluid.
    //----------------------------------------------------------------------
    Gravity gravity(Vecd(0.0, -gravity_g, 0.0));
    SimpleDynamics<GravityForce<Gravity>> constant_gravity(liquid_core, gravity);
    SimpleDynamics<NormalDirectionFromBodyShape> shell_normal_direction(capsule_shell);
    /** wall normals for fluid-wall interaction (plates are translated, not rotated,
     * so a one-time initialization is sufficient). */
    SimpleDynamics<NormalDirectionFromBodyShape> lower_plate_normal_direction(lower_plate);
    SimpleDynamics<NormalDirectionFromBodyShape> upper_plate_normal_direction(upper_plate);

    Dynamics1Level<fluid_dynamics::Integration1stHalfWithWallRiemann> pressure_relaxation(liquid_inner, liquid_contact);
    Dynamics1Level<fluid_dynamics::Integration2ndHalfWithWallRiemann> density_relaxation(liquid_inner, liquid_contact);
    InteractionWithUpdate<fluid_dynamics::DensitySummationComplexFreeSurface> update_density_by_summation(liquid_inner, liquid_contact);
    ReduceDynamics<fluid_dynamics::AdvectionTimeStep> get_fluid_advection_time_step_size(liquid_core, U_ref);
    ReduceDynamics<fluid_dynamics::AcousticTimeStep> get_fluid_time_step_size(liquid_core);
    /** DampingPairwiseWithWall resolves the literal state variable "Velocity" on
     * every contact body of liquid_contact (shell + both plates). The shell only
     * registers "Velocity" inside its NOSB-PD dynamics, and the MovingWallSolid
     * plates register "AverageVelocity" instead, so neither is present yet at this
     * point. Register a "Velocity" state on all three here (idempotent; the later
     * dynamics reuse the same storage, as already done throughout nosbpd_dynamics)
     * to make the damping operator well-defined. */
    capsule_shell.getBaseParticles().registerStateVariableData<Vecd>("Velocity");
    lower_plate.getBaseParticles().registerStateVariableData<Vecd>("Velocity");
    upper_plate.getBaseParticles().registerStateVariableData<Vecd>("Velocity");
    /** Tung-oil viscous damping of the liquid core: this is the intended use of
     * mu_f, which so far was defined but never applied. It mirrors the library
     * fluid-shell benchmark test_2d_hydrostatic_fluid_shell, where the same
     * DampingPairwiseWithWall / FixedDampingRate pair damps the inner fluid and
     * its wall contact with the physical dynamic viscosity. Without it the core
     * has no viscous dissipation at all, so the pole-region fluid velocity
     * fluctuations feed the shell-liquid pressure coupling unchecked. */
    DampingWithRandomChoice<InteractionSplit<DampingPairwiseWithWall<Vecd, FixedDampingRate>>>
        fluid_damping(0.2, DynamicsArgs(liquid_inner, "Velocity", mu_f), DynamicsArgs(liquid_contact, "Velocity", mu_f));
    /** velocity divergence must be constructed before the acoustic damper,
     * and before PressureForceFromFluidforPD which reads it from the fluid. */
    InteractionDynamics<fluid_dynamics::VelocityDivergenceInner> velocity_divergence_inner(liquid_inner);
    /** constructed before the wall velocity divergence so that the shell's
     * AverageVelocity is registered in time. */
    SimpleDynamics<solid_dynamics::PDTimeStepInitialization> initialize_a_pd_step(capsule_shell, makeShared<Gravity>(Vecd(0.0, -gravity_g, 0.0)));
    /** Interface dissipation: the shell-liquid coupling force uses an acoustic
     * Riemann solver instead of NoRiemannSolver, so the fluid-solid interface
     * (the exact location of the pole hotspot, where the shell-plate-liquid gap
     * closes) carries a numerical dissipation term instead of relying on the
     * central pressure difference alone. */
    InteractionDynamics<solid_dynamics::PressureForceFromFluidforPD<AcousticRiemannSolver>> fluid_pressure_force_on_shell(shell_liquid_contact);
    fluid_pressure_force_on_shell.setCoeffAcousticDamper(rho0_f, c_f, liquid_core.getSPHAdaptation().getKernel()->CutOffRadius());
    /** The rigid plates carry their prescribed wall velocity in "AverageVelocity"
     * (registered by MovingWallSolid::AverageVelocity, i.e. already set up when the
     * standard wall classes above were constructed): the lower plate is stationary,
     * the upper plate moves downwards at the piston speed. VelocityDivergenceWithWallforPD
     * reads the same "AverageVelocity", so both wall paths agree. */
    InteractionDynamics<fluid_dynamics::VelocityDivergenceWithWallforPD> velocity_divergence_wall(liquid_contact);
    InteractionDynamics<fluid_dynamics::AcousticDamperInnerforPD> acoustic_damper_inner(liquid_inner);
    acoustic_damper_inner.setCoeffAcousticDamper(rho0_f, c_f, liquid_core.getSPHAdaptation().getKernel()->CutOffRadius());
    InteractionDynamics<fluid_dynamics::AcousticDamperWithWallforPD> acoustic_damper_wall(liquid_contact);
    acoustic_damper_wall.setCoeffAcousticDamper(rho0_f, c_f, liquid_core.getSPHAdaptation().getKernel()->CutOffRadius());
    /** No-penetration is provided by the standard SPHinXsys fluid-wall treatment
     * already present here: Integration1stHalf/2ndHalfWithWallRiemann plus
     * VelocityDivergenceWithWallforPD and AcousticDamperWithWallforPD, which impose
     * the mirror wall velocity 2*v_wall - v_i at the shell and thereby drive the
     * continuity/pressure response that repels fluid from the surface. The former
     * BoundingFromWall hard clamp is removed: it was a position snap-back acting on
     * a stale shell surface, and at the bottom pole it ejected fluid particles
     * (including across the lower plate) instead of preventing penetration. */
    std::cout << "fluid methods done" << std::endl;
    //----------------------------------------------------------------------
    //	Numerical methods: FSI (fluid -> shell).
    //	The PD time-step initialization and fluid pressure force are
    //	constructed above, next to the fluid wall methods (registration order).
    //----------------------------------------------------------------------
    solid_dynamics::AverageVelocityAndAcceleration average_velocity_and_acceleration(capsule_shell);
    //----------------------------------------------------------------------
    //	Numerical methods: shell-plate contact.
    //----------------------------------------------------------------------
    /** the repulsion factor must be constructed (and executed) before the
     * contact force classes reading "RepulsionFactor". */
    InteractionDynamics<solid_dynamics::ContactRepulsionFactorforPD> shell_contact_factor(shell_plates_contact);
    InteractionDynamics<solid_dynamics::ContactRepulsionFromWallforPD> plate_contact_on_shell(shell_plates_contact, plate_friction_mu, plate_friction_veps);
    InteractionDynamics<solid_dynamics::ContactRepulsionToWallforPD> shell_reaction_on_lower(lower_plate_shell_contact);
    InteractionDynamics<solid_dynamics::ContactRepulsionToWallforPD> shell_reaction_on_upper(upper_plate_shell_contact);
    std::cout << "contact methods done" << std::endl;
    //----------------------------------------------------------------------
    //	Numerical methods: NOSB-PD dynamics of the shell.
    //----------------------------------------------------------------------
    InteractionWithUpdate<solid_dynamics::NosbPDShapeMatrix> shell_shape_matrix(shell_inner);
    SimpleDynamics<solid_dynamics::NosbPDFirstStep> nosbpd_first_step(capsule_shell);
    /** constitutive update and bond-breaking criterion follow the shell route. */
    UniquePtr<BaseDynamics<void>> nosbpd_second_step, check_bondLive;
    if (shell_viscoplastic)
    {
        nosbpd_second_step = makeUnique<InteractionWithUpdate<solid_dynamics::NosbPDSecondStepPlastic>>(shell_inner);
        check_bondLive = makeUnique<InteractionDynamics<solid_dynamics::BondBreakByPlasticStrain>>(shell_inner, critical_plastic_strain);
    }
    else
    {
        nosbpd_second_step = makeUnique<InteractionWithUpdate<solid_dynamics::NosbPDSecondStep>>(shell_inner);
        /** bond breaking by the maximum principal stress criterion (brittle shell);
         * evaluated right after the stress update. */
        check_bondLive = makeUnique<InteractionDynamics<solid_dynamics::BondBreakByPrinStress>>(shell_inner, max_tension_stress);
    }
    InteractionDynamics<solid_dynamics::LittleWoodHourGlassControl> hourglass_control(shell_inner);
    InteractionDynamics<solid_dynamics::PairNumericalDampingforPD> numerical_damping(shell_inner);
    InteractionDynamics<solid_dynamics::NosbPDThirdStep> nosbpd_third_step(shell_inner);
    SimpleDynamics<solid_dynamics::NosbPDFourthStep> nosbpd_fourth_step(capsule_shell);
    /** normal direction updated from the PD deformation gradient;
     * constructed after the NOSB-PD dynamics registering DeformationGradient. */
    SimpleDynamics<solid_dynamics::UpdateElasticNormalDirection> shell_update_normal(capsule_shell);
    std::cout << "pd methods done" << std::endl;
    //----------------------------------------------------------------------
    //	Plate kinematics: upper plate driven, lower plate fixed.
    //----------------------------------------------------------------------
    Vecd upper_plate_center_end(0.0, upper_plate_bottom + 0.5 * plate_thickness - compression_distance, 0.0);
    SimpleDynamics<solid_dynamics::PositionSolidBody> drive_upper_plate(upper_plate, 0.0, end_time, upper_plate_center_end);
    //----------------------------------------------------------------------
    //	Output.
    //----------------------------------------------------------------------
    BodyStatesRecordingToVtp write_body_states(sph_system);
    ReducedQuantityRecording<QuantitySummation<Vecd>> write_upper_plate_reaction(upper_plate, "ContactRepulsionToWall");
    //----------------------------------------------------------------------
    //	Prepare the simulation.
    //----------------------------------------------------------------------
    std::cout << "particle generation done\n"
              << " liquid: " << liquid_core.getBaseParticles().TotalRealParticles()
              << " shell: " << capsule_shell.getBaseParticles().TotalRealParticles()
              << " lower: " << lower_plate.getBaseParticles().TotalRealParticles()
              << " upper: " << upper_plate.getBaseParticles().TotalRealParticles() << std::endl;
    sph_system.initializeSystemCellLinkedLists();
    std::cout << "cell linked lists done" << std::endl;
    sph_system.initializeSystemConfigurations();
    std::cout << "configurations done" << std::endl;
    if (with_liquid) constant_gravity.exec();
    shell_normal_direction.exec();
    lower_plate_normal_direction.exec();
    upper_plate_normal_direction.exec();
    std::cout << "normal direction done" << std::endl;
    shell_shape_matrix.exec();
    std::cout << "shape matrix done" << std::endl;
    //----------------------------------------------------------------------
    //	Setup for time-stepping control.
    //----------------------------------------------------------------------
    Real &physical_time = *sph_system.getSystemVariableDataByName<Real>("PhysicalTime");
    int number_of_iterations = 0;
    int screen_output_interval = 100;
    Real output_interval = end_time / 50.0;
    Real dt = 0.0;
    /** PD subcycling step. The shell acoustic CFL limit is
     *  dt_CFL = 0.6 * h_shell / c0 ~= 0.6 * 1.3 * resolution_shell / sqrt(K/rho)
     *         ~= 4.0e-8 s (c0 ~= 129 m/s, h_shell ~= 8.67e-6 m).
     *  The previous fixed 6.0e-8 exceeded this limit by ~1.5x, leaving the
     *  shell time integration marginally unstable and producing the det F
     *  compaction -> bond-break -> non-finite-velocity cascade. Use a step
     *  safely below the CFL limit. */
    Real dt_s_0 = 3.0e-8;
    /** dry control has no fluid CFL: use a fixed step matching the fluid run (~3e-7) so the
     * physical-time progression is comparable to the liquid-filled case. */
    Real dt_dry = 3.0e-7;
    //----------------------------------------------------------------------
    //	First output before the main loop.
    //----------------------------------------------------------------------
    write_body_states.writeToFile(0);
    //----------------------------------------------------------------------
    //	Main time-stepping loop.
    //----------------------------------------------------------------------
    while (physical_time < end_time)
    {
        Real integration_time = 0.0;
        while (integration_time < output_interval)
        {
            Real Dt = with_liquid ? get_fluid_advection_time_step_size.exec() : dt_dry;
            if (with_liquid) update_density_by_summation.exec();

            Real relaxation_time = 0.0;
            while (relaxation_time < Dt)
            {
                /** Refresh the shell surface (NormalDirection / SignedDistance from the
                 * current DeformationGradient) every fluid sub-step, so the fluid-wall
                 * interaction always sees the current shell instead of the surface left
                 * over from the previous output interval. This aligns the shell update
                 * with the fluid sub-step cadence. */
                shell_update_normal.exec();
                /** fluid relaxation and force computation (skipped in the dry control). */
                if (with_liquid)
                {
                    fluid_damping.exec(dt);
                    pressure_relaxation.exec(dt);
                    acoustic_damper_inner.exec(dt);
                    acoustic_damper_wall.exec(dt);
                }
                /** external loads on the shell accumulate onto PDAccelerationPrior. */
                initialize_a_pd_step.exec(dt);
                if (with_liquid) fluid_pressure_force_on_shell.exec(dt);
                shell_contact_factor.exec(dt);
                plate_contact_on_shell.exec(dt);
                if (with_liquid)
                {
                    density_relaxation.exec(dt);
                    velocity_divergence_inner.exec(dt);
                    velocity_divergence_wall.exec(dt);
                }
                /** NOSB-PD subcycling. */
                Real dt_s_sum = 0.0;
                average_velocity_and_acceleration.initialize_displacement_.exec();
                while (dt_s_sum < dt)
                {
                    Real dt_s = dt_s_0;
                    if (dt - dt_s_sum < dt_s) dt_s = dt - dt_s_sum;

                    nosbpd_first_step.exec(dt_s);
                    nosbpd_second_step->exec(dt_s);
                    hourglass_control.exec(dt_s);
                    numerical_damping.exec(dt_s);
                    nosbpd_third_step.exec(dt_s);
                    nosbpd_fourth_step.exec(dt_s);

                    dt_s_sum += dt_s;
                }
                /** bond breaking checked sparsely (eigen solve per bond is
                 * expensive); every 10 fluid steps is ample for fracture. */
                if (number_of_iterations % 10 == 0) check_bondLive->exec(dt);
                average_velocity_and_acceleration.update_averages_.exec(dt);

                dt = with_liquid ? get_fluid_time_step_size.exec() : dt_dry;
                relaxation_time += dt;
                integration_time += dt;
                physical_time += dt;

                /** plate kinematics (dt-dependent) and contact reactions. */
                drive_upper_plate.exec(dt);
                shell_reaction_on_lower.exec(dt);
                shell_reaction_on_upper.exec(dt);
            }

            number_of_iterations++;
            if (number_of_iterations % screen_output_interval == 0 || number_of_iterations <= 200)
            {
                std::cout << std::scientific << "N=" << number_of_iterations
                          << " Time: " << physical_time << "	dt = " << dt << "\n";
            }

            capsule_shell.updateCellLinkedList();
            upper_plate.updateCellLinkedList();
            if (with_liquid)
            {
                liquid_core.updateCellLinkedList();
                liquid_complex.updateConfiguration();
                shell_liquid_contact.updateConfiguration();
            }
            shell_plates_contact.updateConfiguration();
            lower_plate_shell_contact.updateConfiguration();
            upper_plate_shell_contact.updateConfiguration();
        }

        write_body_states.writeToFile(number_of_iterations);
        write_upper_plate_reaction.writeToFile(number_of_iterations);
    }
    return 0;
}
