/** Single capsule between prescribed rigid plates. SI units throughout.
 * PD shell + WCSPH core. Parameters and numerical diagnostics are saved per run.
 * Reference PD neighborhoods remain fixed; contact neighborhoods evolve.
 */
#include "sphinxsys.h"
#include "pd_case_config.h"
#include "pd_checkpoint.h"
#include <regex>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <chrono>

using namespace SPH;

Real capsule_radius_out = 0.5e-3, shell_thickness = 0.02e-3, shell_layers = 3;
Real plate_radius = 0.75e-3, plate_thickness = 0.06e-3, initial_gap = 1e-5;
Real rho0_shell = 1000, youngs_shell = 1e7, poisson_shell = 0.4;
Real max_tension_stress = 1e5, yield_stress_shell = 5e5, hardening_shell = 5e6;
Real viscosity_shell = 500, critical_plastic_strain = 1;
Real plate_friction_mu = 0.3, plate_friction_veps = 0.005;
Real rho0_f = 935, c_f = 25, mu_f = 0.25, fluid_spacing_ratio = 1.5;
Real gravity_g = 9.8, compression_speed = 0.1, compression_distance = 0.58e-3;
Real pd_horizon = 0, pd_pressure_acceleration_factor = 0;
Real U_ref = 0.5, pd_cfl = 0.1, acoustic_cfl = 0.3, acoustic_damping_factor = 0.3;
Real max_density_error = 0.1, max_fluid_mach = 0.3, hourglass_coefficient = 0.5, hourglass_damage_cutoff = 0.1;
bool shell_viscoplastic = true, with_liquid = true, reinitialize_density = false, analytic_plate_contact = true;
int output_frames = 50, diagnostic_frames = 0, bond_check_interval = 1, inject_nan_step = -1, max_pd_substeps = 6;
Real capsule_radius_in, resolution_ref, capsule_center_y, upper_plate_bottom, end_time;

class CapsuleShellShape : public ComplexShape
{
  public:
    explicit CapsuleShellShape(const std::string &name) : ComplexShape(name)
    {
        Vecd center(0, capsule_center_y, 0);
        add<GeometricShapeBall>(center, capsule_radius_out);
        subtract<GeometricShapeBall>(center, capsule_radius_in);
    }
};

class PlateShape : public ComplexShape
{
  public:
    PlateShape(const std::string &name, Real center_y) : ComplexShape(name)
    {
        add<GeometricShapeBox>(Transform(Vecd(0, center_y, 0)), Vecd(plate_radius, 0.5 * plate_thickness, plate_radius));
    }
};

static peridynamics::CaseParameters parameters()
{
    peridynamics::CaseParameters p;
#define BIND(x) p.bind(#x, x)
    BIND(capsule_radius_out); BIND(shell_thickness); BIND(shell_layers);
    BIND(plate_radius); BIND(plate_thickness); BIND(initial_gap);
    BIND(rho0_shell); BIND(youngs_shell); BIND(poisson_shell); BIND(max_tension_stress);
    BIND(yield_stress_shell); BIND(hardening_shell); BIND(viscosity_shell); BIND(critical_plastic_strain);
    BIND(plate_friction_mu); BIND(plate_friction_veps); BIND(rho0_f); BIND(c_f); BIND(mu_f);
    BIND(fluid_spacing_ratio); BIND(gravity_g); BIND(compression_speed); BIND(compression_distance);
    BIND(pd_pressure_acceleration_factor); BIND(pd_horizon); BIND(analytic_plate_contact); BIND(reinitialize_density); BIND(U_ref); BIND(pd_cfl); BIND(acoustic_cfl); BIND(acoustic_damping_factor);
    BIND(hourglass_coefficient); BIND(hourglass_damage_cutoff); BIND(max_density_error); BIND(max_fluid_mach); BIND(shell_viscoplastic); BIND(with_liquid);
    BIND(max_pd_substeps); BIND(output_frames); BIND(diagnostic_frames); BIND(bond_check_interval); BIND(inject_nan_step);
#undef BIND
    return p;
}

static void validateParameters()
{
    if (capsule_radius_out <= shell_thickness || shell_thickness <= 0 || shell_layers < 3 ||
        plate_radius <= capsule_radius_out || plate_thickness <= 0 || initial_gap < 0 ||
        compression_distance <= 0 || compression_distance >= 2 * capsule_radius_out || compression_speed <= 0 ||
        rho0_shell <= 0 || youngs_shell <= 0 || poisson_shell <= -1 || poisson_shell >= 0.5 ||
        rho0_f <= 0 || c_f <= 0 || mu_f < 0 || fluid_spacing_ratio <= 0 ||
        diagnostic_frames < 0 || output_frames < 1 || bond_check_interval < 1 || max_pd_substeps < 1 || pd_cfl <= 0 || pd_cfl > 0.2 ||
        acoustic_cfl <= 0 || acoustic_cfl > 0.6 || acoustic_damping_factor < 0 ||
        max_density_error <= 0 || max_density_error >= 1 || max_fluid_mach <= 0 ||
        plate_friction_mu < 0 || plate_friction_veps <= 0 ||
        yield_stress_shell <= 0 || hardening_shell < 0 || viscosity_shell < 0 ||
        max_tension_stress <= 0 || critical_plastic_strain <= 0 || U_ref <= 0 ||
        pd_pressure_acceleration_factor < 0 || pd_pressure_acceleration_factor > 1 ||
        hourglass_coefficient < 0 || hourglass_damage_cutoff <= 0 || hourglass_damage_cutoff > 1 || inject_nan_step < -1)
        throw std::invalid_argument("invalid capsule geometry, material, schedule or numerical controls");
    capsule_radius_in = capsule_radius_out - shell_thickness;
    resolution_ref = shell_thickness / shell_layers;
    capsule_center_y = plate_thickness + initial_gap + capsule_radius_out;
    upper_plate_bottom = capsule_center_y + capsule_radius_out + initial_gap;
    if (pd_horizon == 0) pd_horizon = 3.015 * resolution_ref;
    if (pd_horizon < 2 * resolution_ref) throw std::invalid_argument("PD horizon must span at least two particle spacings");
    end_time = compression_distance / compression_speed;
    if (diagnostic_frames == 0) diagnostic_frames = output_frames;
    if (diagnostic_frames < output_frames || diagnostic_frames % output_frames != 0)
        throw std::invalid_argument("diagnostic_frames must be a positive multiple of output_frames");
}

static Vecd sumVector(SPHBody &body, const std::string &name)
{
    auto &p = body.getBaseParticles();
    const Vecd *values = p.getVariableDataByName<Vecd>(name);
    Vecd total = Vecd::Zero();
    for (size_t i = 0; i < p.TotalRealParticles(); ++i) total += values[i];
    return total;
}

int main(int argc, char **argv)
{
    auto config = parameters();peridynamics::CheckpointOptions checkpoint_options;
    try {checkpoint_options.parse(argc,argv);config.parse(argc, argv); validateParameters();}
    catch (const std::exception &e) {std::cerr << e.what() << std::endl; return 2;}
    if(checkpoint_options.describe_only)
    {std::cout << "{\"format\":1,\"full_pd_state\":true,\"default_threads\":" << checkpoint_options.threads << "}\n";return 0;}
    Real extent = 6 * capsule_radius_out;
    peridynamics::CheckpointOutputGuard output_guard(!checkpoint_options.resume.empty());
    SPHSystem system(BoundingBoxd(Vecd(-extent, -extent, -extent), Vecd(extent, 2 * extent, extent)), resolution_ref, checkpoint_options.threads);
    output_guard.restore();
    system.handleCommandlineOptions(argc, argv);
    if (system.RestartStep() != 0)
    {
        std::cerr << "Use --pd-resume <checkpoint directory>; upstream restart_step is not a complete PD checkpoint.\n";
        return 2;
    }
    std::filesystem::create_directories("output");
    if(checkpoint_options.resume.empty()) config.save("output/effective_parameters.cfg");

    FluidBody liquid(system, makeShared<GeometricShapeBall>(Vecd(0, capsule_center_y, 0), capsule_radius_in, "LiquidCore"));
    liquid.defineAdaptationRatios(1.3, 1 / fluid_spacing_ratio);
    liquid.defineMatterMaterial<WeaklyCompressibleFluid>(rho0_f, c_f);
    liquid.addMaterialProperty<Viscosity>(mu_f);
    liquid.generateParticles<peridynamics::CheckpointParticles, Lattice>();
    PDBody shell(system, makeShared<CapsuleShellShape>("CapsuleShell"));
    shell.defineAdaptationRatios(0.5 * pd_horizon / resolution_ref);
    if (shell_viscoplastic)
        shell.defineMatterMaterial<ViscoPlasticSolidforPD>(rho0_shell, youngs_shell, poisson_shell, yield_stress_shell, viscosity_shell, hardening_shell);
    else shell.defineMatterMaterial<HughesWingetSolid>(rho0_shell, youngs_shell, poisson_shell);
    shell.generateParticles<peridynamics::CheckpointParticles, Lattice>();
    SolidBody lower(system, makeShared<PlateShape>("LowerPlate", 0.5 * plate_thickness));
    lower.defineAdaptationRatios(1.3, 1.0 / 3.0);
    lower.defineMatterMaterial<Solid>(rho0_shell, youngs_shell);
    lower.generateParticles<peridynamics::CheckpointParticles, Lattice>();
    SolidBody upper(system, makeShared<PlateShape>("UpperPlate", upper_plate_bottom + 0.5 * plate_thickness));
    upper.defineAdaptationRatios(1.3, 1.0 / 3.0);
    upper.defineMatterMaterial<Solid>(rho0_shell, youngs_shell);
    upper.generateParticles<peridynamics::CheckpointParticles, Lattice>();

    InnerRelation fluid_inner(liquid), shell_inner(shell);
    ContactRelationSelectKernel fluid_walls(liquid, {&shell, &lower, &upper}, true);
    ContactRelationSelectKernel shell_fluid(shell, {&liquid}, false);
    ContactRelationSelectKernel upper_fluid(upper, {&liquid}, false), lower_fluid(lower, {&liquid}, false);
    ContactRelation shell_plates(shell, {&lower, &upper}), lower_shell(lower, {&shell}), upper_shell(upper, {&shell});
    ComplexRelation fluid_complex(fluid_inner, fluid_walls);

    SimpleDynamics<NormalDirectionFromBodyShape> shell_normals(shell), lower_normals(lower), upper_normals(upper);
    Dynamics1Level<fluid_dynamics::Integration1stHalfWithWallforPD> pressure(fluid_inner, fluid_walls, pd_pressure_acceleration_factor);
    Dynamics1Level<fluid_dynamics::Integration2ndHalfWithWallRiemann> density(fluid_inner, fluid_walls);
    InteractionWithUpdate<fluid_dynamics::DensitySummationComplexFreeSurface> density_sum(fluid_inner, fluid_walls);
    InteractionWithUpdate<fluid_dynamics::ViscousForceWithWallforPD> fluid_viscosity(fluid_inner, fluid_walls);
    SimpleDynamics<fluid_dynamics::ContinuumVolumeUpdate> update_fluid_volume(liquid);
    SimpleDynamics<peridynamics::FluidKineticWork> fluid_work(liquid);
    ReduceDynamics<fluid_dynamics::AdvectionViscousTimeStep> advection_step(liquid, U_ref);
    ReduceDynamics<fluid_dynamics::AcousticTimeStep> acoustic_step(liquid, acoustic_cfl);
    SimpleDynamics<GravityForce<Gravity>> fluid_gravity(liquid, Gravity(Vecd(0, -gravity_g, 0)));
    SimpleDynamics<solid_dynamics::PDTimeStepInitialization> initialize_shell(shell, makeShared<Gravity>(Vecd(0, -gravity_g, 0)));
    InteractionDynamics<fluid_dynamics::VelocityDivergenceInner> divergence(fluid_inner);
    InteractionDynamics<fluid_dynamics::VelocityDivergenceWithWallforPD> wall_divergence(fluid_walls);
    InteractionDynamics<fluid_dynamics::AcousticDamperInnerforPD> inner_damping(fluid_inner);
    InteractionDynamics<fluid_dynamics::AcousticDamperWithWallforPD> wall_damping(fluid_walls);
    Real damping_c = c_f * acoustic_damping_factor / 0.3;
    Real fluid_horizon = liquid.getSPHAdaptation().getKernel()->CutOffRadius();
    inner_damping.setCoeffAcousticDamper(rho0_f, damping_c, fluid_horizon);
    wall_damping.setCoeffAcousticDamper(rho0_f, damping_c, fluid_horizon);
    InteractionDynamics<solid_dynamics::PressureForceFromFluidforPD<AcousticRiemannSolver>> shell_pressure(shell_fluid);
    shell_pressure.setCoeffAcousticDamper(rho0_f, damping_c, fluid_horizon);
    shell_pressure.setAccelerationPressureFactor(pd_pressure_acceleration_factor);
    InteractionWithUpdate<solid_dynamics::ViscousForceFromFluidforPD> shell_viscosity(shell_fluid);
    InteractionWithUpdate<solid_dynamics::PressureForceFromFluid<decltype(density)>> upper_pressure(upper_fluid), lower_pressure(lower_fluid);
    InteractionWithUpdate<solid_dynamics::ViscousForceFromFluid> upper_viscosity(upper_fluid), lower_viscosity(lower_fluid);

    InteractionDynamics<solid_dynamics::ContactRepulsionFactorforPD> contact_factor(shell_plates);
    InteractionDynamics<solid_dynamics::ContactRepulsionFromWallforPD> contact_force(shell_plates, plate_friction_mu, plate_friction_veps);
    InteractionDynamics<solid_dynamics::ContactRepulsionToWallforPD> lower_reaction(lower_shell, plate_friction_mu, plate_friction_veps);
    InteractionDynamics<solid_dynamics::ContactRepulsionToWallforPD> upper_reaction(upper_shell, plate_friction_mu, plate_friction_veps);
    SimpleDynamics<solid_dynamics::RigidBoxContactforPD> box_contact(shell, plate_friction_mu, plate_friction_veps);
    Vecd plate_half_size(plate_radius, 0.5 * plate_thickness, plate_radius);
    box_contact.applyToCurrentAcceleration(true);
    box_contact.addBox("LowerPlate", Vecd(0, 0.5 * plate_thickness, 0), plate_half_size);
    box_contact.addBox("UpperPlate", Vecd(0, upper_plate_bottom + 0.5 * plate_thickness, 0), plate_half_size);
    InteractionWithUpdate<solid_dynamics::NosbPDShapeMatrix> shape_matrix(shell_inner);
    SimpleDynamics<solid_dynamics::NosbPDFirstStep> first(shell);
    UniquePtr<BaseDynamics<void>> second, fracture;
    if (shell_viscoplastic)
    {
        second = makeUnique<InteractionWithUpdate<solid_dynamics::NosbPDSecondStepPlastic>>(shell_inner);
        fracture = makeUnique<InteractionDynamics<solid_dynamics::BondBreakByPlasticStrain>>(shell_inner, critical_plastic_strain);
    }
    else
    {
        second = makeUnique<InteractionWithUpdate<solid_dynamics::NosbPDSecondStep>>(shell_inner);
        fracture = makeUnique<InteractionDynamics<solid_dynamics::BondBreakByPrinStress>>(shell_inner, max_tension_stress);
    }
    InteractionDynamics<solid_dynamics::LittleWoodHourGlassControl> hourglass(shell_inner, hourglass_coefficient, hourglass_damage_cutoff);
    InteractionDynamics<solid_dynamics::PairNumericalDampingforPD> pd_damping(shell_inner);
    InteractionDynamics<solid_dynamics::NosbPDThirdStep> third(shell_inner);
    SimpleDynamics<solid_dynamics::NosbPDFourthStep> fourth(shell);
    ReduceDynamics<solid_dynamics::PDTimeStep> shell_step(shell, pd_cfl);
    SimpleDynamics<solid_dynamics::UpdateElasticNormalDirection> update_normals(shell);
    solid_dynamics::AverageVelocityAndAcceleration averages(shell);
    SimpleDynamics<solid_dynamics::PrescribedTranslation> move_upper(upper, 0, end_time, Vecd(0, -compression_distance, 0));

    SimpleDynamics<peridynamics::ParticleStateCheck> fluid_check(liquid), shell_check(shell);
    fluid_check.addScalar("Density", 0); fluid_check.addScalar("Pressure");
    fluid_check.addScalar("Mass", 0); shell_check.addScalar("Mass", 0);
    fluid_check.addVector("Force"); fluid_check.addVector("ForcePrior");
    shell_check.addMatrix("DeformationGradient"); shell_check.addMatrix("PDStress");
    shell_check.addVector("PDAcceleration"); shell_check.addVector("PDAccelerationPrior");
    for (const auto &name : {"Mass", "Density", "Pressure", "VolumetricMeasure"}) liquid.getBaseParticles().addVariableToWrite<Real>(name);
    for (const auto &name : {"Velocity", "AverageVelocity", "AverageAcceleration", "PDAcceleration", "PDAccelerationPrior"}) shell.getBaseParticles().addVariableToWrite<Vecd>(name);
    for (const auto &name : {"PDStress", "DeformationGradient"}) shell.getBaseParticles().addVariableToWrite<Matd>(name);
    if (shell_viscoplastic)
    {
        shell.getBaseParticles().addVariableToWrite<Matd>("PlasticStrain");
        shell.getBaseParticles().addVariableToWrite<Real>("IsotropicHardeningParam");
    }
    for (const auto &name : {"Mass", "Density", "VolumetricMeasure"}) shell.getBaseParticles().addVariableToWrite<Real>(name);
    BodyStatesRecordingToVtp states(system);
    states.setFullPrecisionOutput(true);
    Real &time = *system.getSystemVariableDataByName<Real>("PhysicalTime");
    Real pd_dt = pd_cfl * resolution_ref / DynamicCast<HughesWingetSolid>(&shell, shell.getMatterMaterial()).ReferenceSoundSpeed();
    std::ofstream diagnostics,reaction;
    std::ofstream budget,mass_record;
    Real plate_fluid_pressure_work = 0, plate_fluid_viscous_work = 0;
    Real initial_gravity_energy = 0, initial_stored_energy = 0;
    Real initial_shell_mass = 0, initial_liquid_mass = 0;
    size_t initial_shell_count = shell.getBaseParticles().TotalRealParticles();
    size_t initial_liquid_count = liquid.getBaseParticles().TotalRealParticles();
    const auto start_wall_clock = std::chrono::steady_clock::now();
    const Real plate_initial_y = upper.getBaseParticles().getVariableDataByName<Vecd>("Position")[0][1];
    StdVec<Real> shell_initial_masses, liquid_initial_masses;
    auto remember_mass = [](SPHBody &body, StdVec<Real> &values)
    {
        auto &p = body.getBaseParticles(); auto *m = p.getVariableDataByName<Real>("Mass");
        values.assign(m, m + p.TotalRealParticles());
    };
    remember_mass(shell, shell_initial_masses); remember_mass(liquid, liquid_initial_masses);
    auto sum_scalar = [](SPHBody &body, const char *name)
    {
        auto &p = body.getBaseParticles(); auto *v = p.getVariableDataByName<Real>(name); Real total = 0;
        for (size_t i = 0; i < p.TotalRealParticles(); ++i) total += v[i];
        return total;
    };
    size_t acoustic_count = 0, pd_count = 0;
    auto last_progress = start_wall_clock;
    std::ofstream progress;
    std::string phase = "initialization";
    uint64_t completed_frame=0;Real wall_offset=0;
    auto open_outputs=[&](bool append)
    {
        auto mode=std::ios::out|(append ? std::ios::app : std::ios::trunc);
        diagnostics.open("output/diagnostics.csv",mode);diagnostics << std::setprecision(17);
        reaction.open("output/UpperPlate_TotalReaction.dat",mode);reaction << std::setprecision(17);
        budget.open("output/energy_budget.csv",mode);budget << std::setprecision(17);
        mass_record.open("output/mass_budget.csv",mode);mass_record << std::setprecision(17);
        progress.open("output/progress.csv",mode);progress << std::setprecision(17);
        if(!diagnostics || !reaction || !budget || !mass_record || !progress) throw std::runtime_error("cannot open diagnostic outputs");
        if(!append)
        {
            diagnostics << std::setprecision(17) << "time_s,stroke_m,acoustic_steps,pd_steps,liquid_mass_kg,fluid_speed_max_m_s,fluid_density_error_max,damage_mean,damage_max,numerical_failures,shell_kinetic_J,liquid_kinetic_J,contact_Fy_N,fluid_Fy_N,total_Fy_N\n";
            reaction << std::setprecision(17) << "time_s stroke_m contact_Fy_N fluid_Fy_N total_Fy_N\n";
            budget << std::setprecision(17) << "time_s,plate_contact_work_J,plate_fluid_pressure_work_J,plate_fluid_viscous_work_J,kinetic_J,gravity_potential_change_J,fluid_compression_J,shell_elastic_proxy_J,contact_penalty_J,plastic_work_J,pd_damping_signed_loss_J,acoustic_inner_signed_loss_J,acoustic_wall_signed_loss_J,viscous_signed_loss_J,friction_signed_loss_J,hourglass_signed_transfer_J,mechanical_residual_J,mechanical_residual_relative,shell_kinetic_identity_residual_J,fluid_kinetic_identity_residual_J\n";
            mass_record << std::setprecision(17) << "time_s,shell_count,liquid_count,shell_mass_kg,liquid_mass_kg,total_mass_kg,shell_mass_relative_error,liquid_mass_relative_error,max_particle_mass_relative_error,liquid_volume_m3,density_volume_mass_residual_kg,liquid_mass_outside_domain_kg\n";
            progress << std::setprecision(17) << "wall_s,time_s,completion_fraction,acoustic_steps,pd_steps,last_fluid_dt_s,last_pd_dt_s\n";
        }
    };
    auto flush_outputs=[&](){diagnostics.flush();reaction.flush();budget.flush();mass_record.flush();progress.flush();
        if(!diagnostics || !reaction || !budget || !mass_record || !progress) throw std::runtime_error("cannot flush diagnostic outputs");};
    std::string checkpoint_identity=config.serialize()+"\nrecording="+std::to_string(system.StateRecording())+
        "\nexecutable="+std::to_string(peridynamics::checkpoint_detail::fileHash(std::filesystem::canonical(argv[0])));
    peridynamics::DiskCheckpoint checkpoint(checkpoint_options.directory,checkpoint_identity);
    for(SPHBody *body:{static_cast<SPHBody*>(&shell),static_cast<SPHBody*>(&liquid),static_cast<SPHBody*>(&lower),static_cast<SPHBody*>(&upper)})
    {
        // The upstream VTP writer lazily registers ParticleGroups. Register it
        // before either save or load so fresh-process schemas match exactly.
        body->getBaseParticles().getParticleGroupManager();
        checkpoint.addBody(dynamic_cast<peridynamics::CheckpointParticles&>(body->getBaseParticles()));
    }
    checkpoint.addConfiguration("ShellReference",shell_inner.inner_configuration_,shell.getBaseParticles().ParticlesBound());
    checkpoint.addConfiguration("FluidInner",fluid_inner.inner_configuration_,liquid.getBaseParticles().ParticlesBound());
    auto add_contact=[&](const std::string &name,BaseContactRelation &relation)
    {
        for(size_t k=0;k<relation.contact_configuration_.size();++k)
            checkpoint.addConfiguration(name+std::to_string(k),relation.contact_configuration_[k],relation.contact_particles_[k]->ParticlesBound());
    };
    add_contact("FluidWalls",fluid_walls);add_contact("ShellFluid",shell_fluid);add_contact("UpperFluid",upper_fluid);add_contact("LowerFluid",lower_fluid);
    add_contact("ShellPlates",shell_plates);add_contact("LowerShell",lower_shell);add_contact("UpperShell",upper_shell);
#define CHECKPOINT_SCALAR(x) checkpoint.addScalar(#x,x)
    CHECKPOINT_SCALAR(time);CHECKPOINT_SCALAR(acoustic_count);CHECKPOINT_SCALAR(pd_count);CHECKPOINT_SCALAR(completed_frame);
    CHECKPOINT_SCALAR(plate_fluid_pressure_work);CHECKPOINT_SCALAR(plate_fluid_viscous_work);
    CHECKPOINT_SCALAR(initial_gravity_energy);CHECKPOINT_SCALAR(initial_stored_energy);CHECKPOINT_SCALAR(initial_shell_mass);CHECKPOINT_SCALAR(initial_liquid_mass);
#undef CHECKPOINT_SCALAR
    checkpoint.addArray("ShellInitialMasses",shell_initial_masses);checkpoint.addArray("FluidInitialMasses",liquid_initial_masses);
    for(const char *name:{"diagnostics.csv","UpperPlate_TotalReaction.dat","energy_budget.csv","mass_budget.csv","progress.csv","effective_parameters.cfg"}) checkpoint.addOutput(name);
    auto last_checkpoint=std::chrono::steady_clock::now();size_t last_checkpoint_step=0;
    auto save_checkpoint=[&](bool force)
    {
        auto now=std::chrono::steady_clock::now();
        bool due=force || peridynamics::checkpoint_requested ||
            (checkpoint_options.seconds>0 && std::chrono::duration<Real>(now-last_checkpoint).count()>=checkpoint_options.seconds) ||
            (checkpoint_options.steps>0 && acoustic_count-last_checkpoint_step>=checkpoint_options.steps);
        if(!due) return false;
        peridynamics::checkpoint_requested=0;
        flush_outputs();
        {std::ofstream info("output/checkpoint_status.json");info << std::setprecision(17) << "{\"status\":\"writing\",\"time_s\":" << time << "}\n";}
        auto path=checkpoint.save({time,wall_offset+std::chrono::duration<Real>(now-start_wall_clock).count(),acoustic_count,pd_count,completed_frame});
        {std::ofstream info("output/checkpoint_status.json");info << std::setprecision(17) << "{\"status\":\"complete\",\"time_s\":" << time << "}\n";}
        std::cout << "CHECKPOINT committed=" << path << " time=" << std::setprecision(17) << time << " acoustic_steps=" << acoustic_count << std::endl;
        last_checkpoint=std::chrono::steady_clock::now();last_checkpoint_step=acoustic_count;
        bool stop=peridynamics::checkpoint_stop_requested || (checkpoint_options.stop_after_checkpoint && time>0);
        if(stop) std::ofstream("output/run_status.txt") << std::setprecision(17) << "checkpointed\ntime_s=" << time << "\ncheckpoint=" << path.string() << '\n';
        return stop;
    };


    auto check_fluid = [&]()
    {
        if (!with_liquid) return;
        fluid_check.exec();
        auto &p = liquid.getBaseParticles();
        auto *rho = p.getVariableDataByName<Real>("Density");
        auto *v = p.getVariableDataByName<Vecd>("Velocity");
        for (size_t i = 0; i < p.TotalRealParticles(); ++i)
        {
            if (std::abs(rho[i] / rho0_f - 1) > max_density_error)
                peridynamics::numericalFailure(liquid.Name(), i, "density deviation limit: rho/rho0=" + std::to_string(rho[i] / rho0_f));
            if (v[i].norm() > max_fluid_mach * c_f)
                peridynamics::numericalFailure(liquid.Name(), i, "Mach-number limit: speed=" + std::to_string(v[i].norm()));
        }
    };
    auto check_bounds = [&](SPHBody &body)
    {
        auto &p = body.getBaseParticles(); auto *x = p.getVariableDataByName<Vecd>("Position");
        for (size_t i = 0; i < p.TotalRealParticles(); ++i)
            if (x[i][0] < -extent || x[i][0] > extent || x[i][1] < -extent || x[i][1] > 2 * extent || x[i][2] < -extent || x[i][2] > extent)
                peridynamics::numericalFailure(body.Name(), i, "particle left computational bounds before neighbor refresh");
    };
    auto update_contacts = [&]()
    {
        check_bounds(shell); if (with_liquid) check_bounds(liquid);
        shell.updateCellLinkedList(); upper.updateCellLinkedList();
        if (with_liquid)
        {
            liquid.updateCellLinkedList(); fluid_complex.updateConfiguration(); shell_fluid.updateConfiguration();
            upper_fluid.updateConfiguration(); lower_fluid.updateConfiguration();
        }
        // RigidBoxContactforPD reads positions directly. Its analytic path never
        // consumes these three particle-contact configurations.
        if (!analytic_plate_contact)
        {
            shell_plates.updateConfiguration(); lower_shell.updateConfiguration(); upper_shell.updateConfiguration();
        }
    };
    auto apply_contact = [&](Real contact_time)
    {
        if (analytic_plate_contact)
        {
            Real displacement = std::min(contact_time / end_time, Real(1)) * compression_distance;
            Vecd velocity(0, contact_time < end_time ? -compression_speed : 0, 0);
            box_contact.setKinematics(1, Vecd(0, upper_plate_bottom + 0.5 * plate_thickness - displacement, 0), velocity);
            box_contact.exec();
        }
        else {contact_factor.exec(); contact_force.exec();}
    };
    auto measure = [&](int frame)
    {
        initialize_shell.exec();
        box_contact.applyAcceleration(false); apply_contact(time); box_contact.applyAcceleration(true);
        if (!analytic_plate_contact) {upper_reaction.exec(); lower_reaction.exec();}
        Real contact_y = (analytic_plate_contact ? sumVector(shell, "UpperPlateReaction") : sumVector(upper, "ContactRepulsionToWall"))[1], fluid_y = 0;
        if (with_liquid)
        {
            upper_pressure.exec(); upper_viscosity.exec();
            fluid_y = (sumVector(upper, "PressureForceFromFluid") + sumVector(upper, "ViscousForceFromFluid"))[1];
        }
        auto &sp = shell.getBaseParticles(); auto &fp = liquid.getBaseParticles();
        auto *d = sp.getVariableDataByName<Real>("Damage");
        auto *sv = sp.getVariableDataByName<Vecd>("Velocity"); auto *sm = sp.getVariableDataByName<Real>("Mass");
        auto *failure = sp.getVariableDataByName<int>("NumericalFailureReason");
        Real dsum = 0, dmax = 0, shell_ke = 0, fluid_ke = 0, mass = 0, vmax = 0, rho_error = 0;
        int failures = 0;
        for (size_t i = 0; i < sp.TotalRealParticles(); ++i)
        {dsum += d[i]; dmax = std::max(dmax, d[i]); shell_ke += 0.5 * sm[i] * sv[i].squaredNorm(); failures += failure[i] != 0;}
        if (with_liquid)
        {
            auto *fm = fp.getVariableDataByName<Real>("Mass"); auto *fv = fp.getVariableDataByName<Vecd>("Velocity");
            auto *frho = fp.getVariableDataByName<Real>("Density");
            for (size_t i = 0; i < fp.TotalRealParticles(); ++i)
            {mass += fm[i]; vmax = std::max(vmax, fv[i].norm()); rho_error = std::max(rho_error, std::abs(frho[i] / rho0_f - 1)); fluid_ke += 0.5 * fm[i] * fv[i].squaredNorm();}
        }
        // Independently measure the driven plate displacement from particle coordinates.
        Real stroke = plate_initial_y - upper.getBaseParticles().getVariableDataByName<Vecd>("Position")[0][1];
        Real shell_mass = 0, gravity_energy = 0, compression_energy = 0, elastic_proxy = 0, contact_energy = 0;
        Real particle_mass_error = 0, liquid_volume = 0, density_volume_mass = 0, outside_mass = 0;
        auto *sx = sp.getVariableDataByName<Vecd>("Position");
        auto *svol = sp.getVariableDataByName<Real>("VolumetricMeasure");
        auto *stress = sp.getVariableDataByName<Matd>("PDStress");
        auto *F = sp.getVariableDataByName<Matd>("DeformationGradient");
        auto *factor = sp.getVariableDataByName<Real>("RepulsionFactor");
        auto &material = DynamicCast<HughesWingetSolid>(&shell, shell.getMatterMaterial());
        if (sp.TotalRealParticles() != initial_shell_count || fp.TotalRealParticles() != initial_liquid_count)
            throw std::runtime_error("particle count changed in a closed particle inventory");
        for (size_t i = 0; i < sp.TotalRealParticles(); ++i)
        {
            shell_mass += sm[i]; gravity_energy += gravity_g * sm[i] * sx[i][1];
            particle_mass_error = std::max(particle_mass_error, std::abs(sm[i] / shell_initial_masses[i] - 1));
            Matd dev = stress[i] - stress[i].trace() / Dimensions * Matd::Identity();
            elastic_proxy += svol[i] * F[i].determinant() * (dev.squaredNorm() / (4 * material.ShearModulus()) +
                stress[i].trace() * stress[i].trace() / (2 * Dimensions * Dimensions * material.BulkModulus()));
            contact_energy += material.ContactStiffness() * svol[i] * factor[i] * factor[i];
        }
        if (with_liquid)
        {
            auto *fm = fp.getVariableDataByName<Real>("Mass"); auto *fx = fp.getVariableDataByName<Vecd>("Position");
            auto *rho = fp.getVariableDataByName<Real>("Density"); auto *vol = fp.getVariableDataByName<Real>("VolumetricMeasure");
            for (size_t i = 0; i < fp.TotalRealParticles(); ++i)
            {
                gravity_energy += gravity_g * fm[i] * fx[i][1];
                compression_energy += fm[i] * peridynamics::linearEOSSpecificEnergy(rho[i], rho0_f, c_f);
                liquid_volume += vol[i]; density_volume_mass += rho[i] * vol[i];
                particle_mass_error = std::max(particle_mass_error, std::abs(fm[i] / liquid_initial_masses[i] - 1));
                if (fx[i][0] < -extent || fx[i][0] > extent || fx[i][1] < -extent || fx[i][1] > 2 * extent || fx[i][2] < -extent || fx[i][2] > extent)
                    outside_mass += fm[i];
            }
        }
        if (analytic_plate_contact) contact_energy = sum_scalar(shell, "AnalyticContactEnergy");
        Real stored = shell_ke + fluid_ke + compression_energy + elastic_proxy + contact_energy;
        if (frame == 0)
        {
            initial_gravity_energy = gravity_energy; initial_stored_energy = stored;
            initial_shell_mass = shell_mass; initial_liquid_mass = mass;
        }
        Real contact_work = sum_scalar(shell, "ContactBoundaryWork");
        Real plastic_work = sum_scalar(shell, "PlasticWork");
        Real pd_loss = -sum_scalar(shell, "PDDampingWork");
        Real acoustic_inner_loss = with_liquid ? -sum_scalar(liquid, "AcousticDampingWorkInner") : 0;
        Real acoustic_wall_loss = with_liquid ? -sum_scalar(liquid, "AcousticDampingWorkWall") - sum_scalar(shell, "ShellAcousticWallWork") : 0;
        Real viscous_loss = with_liquid ? plate_fluid_viscous_work - sum_scalar(liquid, "FluidViscousWork") - sum_scalar(shell, "ShellViscousWork") : 0;
        Real friction_loss = sum_scalar(shell, "ContactFrictionLoss");
        Real hg_transfer = -sum_scalar(shell, "PDHourglassWork");
        Real supplied = contact_work + plate_fluid_pressure_work + plate_fluid_viscous_work;
        Real accounted = stored - initial_stored_energy + gravity_energy - initial_gravity_energy +
            plastic_work + pd_loss + acoustic_inner_loss + acoustic_wall_loss + viscous_loss + friction_loss + hg_transfer;
        Real residual = supplied - accounted;
        Real energy_scale = std::max({std::abs(supplied), std::abs(stored - initial_stored_energy), std::abs(gravity_energy - initial_gravity_energy),
            std::abs(plastic_work), std::abs(pd_loss), std::abs(acoustic_inner_loss), std::abs(acoustic_wall_loss), std::abs(viscous_loss), std::abs(friction_loss), std::abs(hg_transfer), Real(1e-30)});
        Real shell_identity = shell_ke - sum_scalar(shell, "PDInternalWork") - sum_scalar(shell, "PDExternalWork") - sum_scalar(shell, "PDHourglassWork") - sum_scalar(shell, "PDDampingWork");
        Real fluid_identity = with_liquid ? fluid_ke - sum_scalar(liquid, "FluidPressureWork") - sum_scalar(liquid, "FluidPriorWork") - sum_scalar(liquid, "AcousticDampingWorkInner") - sum_scalar(liquid, "AcousticDampingWorkWall") : 0;
        budget << time << ',' << contact_work << ',' << plate_fluid_pressure_work << ',' << plate_fluid_viscous_work << ',' << shell_ke + fluid_ke << ',' << gravity_energy - initial_gravity_energy << ',' << compression_energy << ',' << elastic_proxy << ',' << contact_energy << ',' << plastic_work << ',' << pd_loss << ',' << acoustic_inner_loss << ',' << acoustic_wall_loss << ',' << viscous_loss << ',' << friction_loss << ',' << hg_transfer << ',' << residual << ',' << residual / energy_scale << ',' << shell_identity << ',' << fluid_identity << '\n';
        mass_record << time << ',' << initial_shell_count << ',' << (with_liquid ? initial_liquid_count : 0) << ',' << shell_mass << ',' << mass << ',' << shell_mass + mass << ',' << shell_mass / initial_shell_mass - 1 << ',' << (with_liquid ? mass / initial_liquid_mass - 1 : 0) << ',' << particle_mass_error << ',' << liquid_volume << ',' << density_volume_mass - mass << ',' << outside_mass << '\n';
        budget.flush(); mass_record.flush();
        diagnostics << time << ',' << stroke << ',' << acoustic_count << ',' << pd_count << ',' << mass << ',' << vmax << ',' << rho_error << ',' << dsum / sp.TotalRealParticles() << ',' << dmax << ',' << failures << ',' << shell_ke << ',' << fluid_ke << ',' << contact_y << ',' << fluid_y << ',' << contact_y + fluid_y << '\n';
        reaction << time << ' ' << stroke << ' ' << contact_y << ' ' << fluid_y << ' ' << contact_y + fluid_y << '\n';
        diagnostics.flush(); reaction.flush();
        int snapshot_stride = diagnostic_frames / output_frames;
        if (frame % snapshot_stride == 0) states.writeToFile(frame / snapshot_stride);
        std::cout << "frame=" << frame << " time=" << time << " stroke=" << stroke << " F=" << contact_y + fluid_y << " rho_error=" << rho_error << " vmax=" << vmax
                  << " energy_residual=" << residual / energy_scale << " wall_s="
                  << std::chrono::duration<Real>(std::chrono::steady_clock::now() - start_wall_clock).count() << std::endl;
    };

    try
    {
        system.initializeSystemCellLinkedLists(); system.initializeSystemConfigurations();
        shell_normals.exec(); lower_normals.exec(); upper_normals.exec(); shape_matrix.exec(); move_upper.exec();
        if (with_liquid) {fluid_gravity.exec(); update_fluid_volume.exec();}
        if(checkpoint_options.resume.empty())
        {
            if(std::filesystem::exists(checkpoint_options.directory/"LATEST"))
                std::filesystem::rename(checkpoint_options.directory,checkpoint_options.directory.string()+".previous-"+std::to_string(std::chrono::steady_clock::now().time_since_epoch().count()));
            open_outputs(false);shell_check.exec();check_fluid();measure(0);
            if(save_checkpoint(true)) return 75;
        }
        else
        {
            auto loaded=checkpoint.restore(checkpoint_options.resume);
            if(time!=loaded.time || acoustic_count!=loaded.acoustic_steps || pd_count!=loaded.pd_steps || completed_frame!=loaded.completed_frame || completed_frame>uint64_t(diagnostic_frames))
                throw std::runtime_error("checkpoint: inconsistent integrator cursor");
            wall_offset=loaded.wall_seconds;last_checkpoint_step=acoustic_count;
            // Rebuild only spatial search structures, never the saved reference
            // bonds or lagged interaction geometry/order used by the integrator.
            shell.updateCellLinkedList();liquid.updateCellLinkedList();lower.updateCellLinkedList();upper.updateCellLinkedList();
            int last_snapshot=int(completed_frame)/(diagnostic_frames/output_frames);
            std::filesystem::path orphan("output/orphaned-after-resume-"+std::to_string(std::chrono::steady_clock::now().time_since_epoch().count()));
            std::regex pattern("(CapsuleShell|LiquidCore|LowerPlate|UpperPlate)_ite_([0-9]+)\\.vtp");
            for(auto &entry:std::filesystem::directory_iterator("output"))
            {
                std::smatch match;std::string name=entry.path().filename().string();
                if(entry.is_regular_file() && std::regex_match(name,match,pattern) && std::stoull(match[2])>uint64_t(last_snapshot))
                {std::filesystem::create_directories(orphan);std::filesystem::rename(entry.path(),orphan/name);}
            }
            open_outputs(true);std::filesystem::remove("output/run_status.txt");
            std::ofstream resume_info("output/resume_info.json");
            resume_info << std::setprecision(17) << "{\"time_s\":" << time << ",\"acoustic_steps\":" << acoustic_count << ",\"completed_frame\":" << completed_frame << ",\"generation\":\"" << checkpoint.loadedGeneration().filename().string() << "\"}\n";
            std::cout << "RESTORED time=" << std::setprecision(17) << time << " generation=" << checkpoint.loadedGeneration() << std::endl;
        }
        // Verification mode: commit restored memory before taking another step.
        if(!checkpoint_options.resume.empty() && checkpoint_options.stop_after_checkpoint && save_checkpoint(true)) return 75;
        std::cout << "particles shell=" << shell.getBaseParticles().TotalRealParticles()
                  << " liquid=" << (with_liquid ? liquid.getBaseParticles().TotalRealParticles() : 0)
                  << " PD_dt=" << pd_dt << " end_time=" << end_time << std::endl;
        for (int frame = int(completed_frame)+1; frame <= diagnostic_frames; ++frame)
        {
            Real target_time = frame == diagnostic_frames ? end_time : end_time * Real(frame) / diagnostic_frames;
            while (time < target_time)
            {
                phase = "density summation";
                if (with_liquid && reinitialize_density) {density_sum.exec(); check_fluid();}
                Real Dt = peridynamics::boundedTimeStep(with_liquid ? advection_step.exec() : pd_dt * 20, {target_time - time});
                Real advection_end = time + Dt;
                while (time < advection_end)
                {
                    Real solid_dt = shell_step.exec();
                    // Bound the explicit fluid-solid exchange as well as the fluid acoustic CFL.
                    Real candidate_dt = with_liquid ? std::min(acoustic_step.exec(), max_pd_substeps * solid_dt) : 20 * solid_dt;
                    Real dt = peridynamics::boundedTimeStep(candidate_dt, {advection_end - time, target_time - time});
                    phase = "fluid first half"; move_upper.exec(); update_normals.exec();
                    if (with_liquid)
                    {
                        // Both directions of explicit viscous traction sample the same
                        // accepted-time velocities, before pressure changes the core.
                        fluid_viscosity.exec(); shell_viscosity.exec(); upper_viscosity.exec();
                        divergence.exec(); wall_divergence.exec();
                        pressure.exec(dt); fluid_work.exec(dt);
                        upper_pressure.exec();
                        plate_fluid_pressure_work += dt * compression_speed * sumVector(upper, "PressureForceFromFluid")[1];
                        plate_fluid_viscous_work += dt * compression_speed * sumVector(upper, "ViscousForceFromFluid")[1];
                        inner_damping.exec(dt); wall_damping.exec(dt); check_fluid();
                    }
                    phase = "external shell loads";
                    initialize_shell.exec();
                    if (with_liquid)
                    {
                        shell_pressure.exec(dt);
                        auto &p = shell.getBaseParticles();
                        auto *a = p.getVariableDataByName<Vecd>("PDAccelerationPrior");
                        auto *f = p.getVariableDataByName<Vecd>("ViscousForceFromFluid");
                        auto *v = p.getVariableDataByName<Real>("VolumetricMeasure");
                        for (size_t i = 0; i < p.TotalRealParticles(); ++i) a[i] += f[i] / (rho0_shell * v[i]);
                    }
                    if (!analytic_plate_contact) apply_contact(time);
                    phase = "fluid second half";
                    if (with_liquid) {density.exec(dt); update_fluid_volume.exec(); check_fluid();}
                    phase = "PD subcycling";
                    Real elapsed = 0; averages.initialize_displacement_.exec();
                    while (elapsed < dt)
                    {
                        Real ds = peridynamics::boundedTimeStep(shell_step.exec(), {dt - elapsed});
                        first.exec(ds);
                        if (analytic_plate_contact) apply_contact(time + elapsed + 0.5 * ds);
                        second->exec(ds);
                        if (pd_count % size_t(bond_check_interval) == 0) fracture->exec(ds);
                        hourglass.exec(ds); pd_damping.exec(ds); third.exec(ds); fourth.exec(ds);
                        shell_check.exec(); ++pd_count;
                        elapsed = ds == dt - elapsed ? dt : elapsed + ds;
                    }
                    averages.update_averages_.exec(dt);
                    ++acoustic_count;
                    if (with_liquid && inject_nan_step >= 0 && acoustic_count == size_t(inject_nan_step))
                        liquid.getBaseParticles().getVariableDataByName<Vecd>("Velocity")[0][0] = std::numeric_limits<Real>::quiet_NaN();
                    phase = "accepted step checks"; check_fluid(); shell_check.exec();
                    time = dt == advection_end - time ? advection_end : time + dt;
                    move_upper.exec();
                    auto now = std::chrono::steady_clock::now();
                    if (std::chrono::duration<Real>(now - last_progress).count() >= 30)
                    {
                        progress << wall_offset + std::chrono::duration<Real>(now - start_wall_clock).count() << ',' << time << ',' << time / end_time << ',' << acoustic_count << ',' << pd_count << ',' << dt << ',' << solid_dt << '\n';
                        progress.flush(); last_progress = now;
                    }
                }
                if (advection_end == target_time) time = target_time;
                phase = "contact neighborhood refresh"; update_contacts();
                phase = "checkpoint";if(save_checkpoint(false)) return 75;
            }
            phase = "output"; measure(frame);completed_frame=frame;
        }
        if(checkpoint_options.final_checkpoint) {phase="final checkpoint";if(save_checkpoint(true)) return 75;}
        std::ofstream("output/run_status.txt") << "completed\ntime_s=" << std::setprecision(17) << time << '\n';
    }
    catch (const std::exception &e)
    {
        std::ofstream failure("output/run_status.txt");
        failure << std::setprecision(17) << "failed\ntime_s=" << time << "\nphase=" << phase
                << "\nacoustic_steps=" << acoustic_count << "\npd_steps=" << pd_count << "\nreason=" << e.what() << '\n';
        failure.close();
        try {states.writeToFile(1000000000 + acoustic_count);} catch (...) {}
        std::cerr << "FAILED at t=" << time << " phase=" << phase << ": " << e.what() << std::endl;
        return 1;
    }
    return 0;
}
