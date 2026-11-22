/** Bounded solid-sphere PD control for a future effective hydrated-gel baseline.
 * No shell cavity, dummy liquid, material fitting or fluid exudation.
 * Generic material/fracture assembly; fracture is disabled by default.
 * Observation evaluates accepted-state contact without changing integration state.
 */
#include "sphinxsys.h"
#include "pd_checkpoint.h"
#include "pd_particle_compression.h"
#include "execution_limits.h"
#include "gel_surface_contact.h"
#include "pd_contact_impulse.h"
#include "pd_execution.h"
#include "pd_time_step_execution.h"

#include <algorithm>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <regex>

using namespace SPH;
namespace fs = std::filesystem;
using gel_compression::ExecutionLimits;
using gel_compression::ResourceStop;

namespace
{
bool consumeDescriptionFlag(int &argc, char **argv)
{
    bool describe = false;
    int output = 1;
    for (int i = 1; i < argc; ++i)
        if (std::string(argv[i]) == "--describe-model")
            describe = true;
        else
            argv[output++] = argv[i];
    argc = output;
    argv[output] = nullptr;
    return describe;
}

void requireGoodStream(const std::ofstream &stream)
{
    if (!stream)
        throw std::runtime_error("gel compression: diagnostic output failed");
}

int simulate(int argc, char **argv, const peridynamics::ParticleCompressionOptions &options,
             const peridynamics::CaseParameters &inputs, const peridynamics::CheckpointOptions &restart,
             const peridynamics::PeriExecution &execution)
{
    // Fail before SPHSystem can create/rotate output, and before particle allocation.
    ExecutionLimits limits;
    if (restart.resume.empty() && (fs::exists("output") || fs::exists(restart.directory / "LATEST")))
        throw std::invalid_argument("use a new run directory; existing gel output/checkpoint is preserved");
    const auto &geometry = options.geometry;
    if (2 * geometry.radius / geometry.spacing > 26)
        throw std::invalid_argument("gel control resolution exceeds the bounded preparation budget");

    peridynamics::CheckpointOutputGuard output_guard(!restart.resume.empty());
    const Real extent = 6 * geometry.radius;
    SPHSystem system(BoundingBoxd(Vecd(-extent, -extent, -extent), Vecd(extent, 2 * extent, extent)),
                     geometry.spacing, restart.threads);
    output_guard.restore();
    system.handleCommandlineOptions(argc, argv);
    if (system.RestartStep() != 0)
        throw std::invalid_argument("upstream restart_step is not a full PD checkpoint; use --pd-resume");

    PDBody particle(system, makeShared<GeometricShapeBall>(
                                Vecd(0, options.centerHeight(), 0), geometry.radius, "GelParticle"));
    particle.defineAdaptationRatios(0.5 * options.numerics.horizon_ratio);
    const auto &material = options.material;
    peridynamics::definePDMaterial(particle, material);
    particle.generateParticles<peridynamics::CheckpointParticles, Lattice>();
    auto &particles = particle.getBaseParticles();
    const size_t count = particles.TotalRealParticles();
    if (count == 0 || count > ExecutionLimits::max_particles)
        throw std::invalid_argument("gel spatial control requires 1..8192 particles");

    InnerRelation reference(particle);
    InteractionWithUpdate<solid_dynamics::NosbPDShapeMatrix> shape(reference);
    SimpleDynamics<solid_dynamics::PDTimeStepInitialization> initialize(particle, makeShared<Gravity>(Vecd::Zero()));
    auto verlet_backend = execution.makePositionVerlet();
    solid_dynamics::NosbPDStepWithExecution<solid_dynamics::NosbPDFirstStep> first(particle, verlet_backend.get());
    auto gel_backend = material.isViscoelastic() ? execution.makeGelPoints() : nullptr;
    auto elastic_backend = material.model == "elastic" ? execution.makeNosbPDElastic() : nullptr;
    auto plastic_backend = material.isPlastic() ? execution.makeNosbPDPlastic() : nullptr;
    std::cout << "PD material backend " << (gel_backend || elastic_backend || plastic_backend ? "cuda_fp64" : "cpu_fp64") << std::endl;
    solid_dynamics::NosbPDBondWorkspace bond_workspace(reference, execution.makeNosbPDBonds());
    if (material.isViscoelastic())
        std::cout << "PD position moment backend "
                  << (gel_backend && bond_workspace.assemblesIntactPositionMoments() ? "cuda_fp64" : "cpu_fp64") << std::endl;
    auto constitutive = peridynamics::makePDConstitutiveDynamics(reference, &bond_workspace, std::move(elastic_backend),
                                                                std::move(plastic_backend), gel_backend.get());
    auto fracture = peridynamics::makePDFracture(reference, material, options.fracture, &bond_workspace);
    InteractionDynamics<solid_dynamics::LittleWoodHourGlassControl> hourglass(
        reference, options.numerics.hourglass_coefficient, 0.1);
    InteractionDynamics<solid_dynamics::PairNumericalDampingforPD> damping(reference);
    solid_dynamics::NosbPDInternalForceWithExecution third(bond_workspace);
    solid_dynamics::NosbPDStepWithExecution<solid_dynamics::NosbPDFourthStep> fourth(particle, verlet_backend.get());
    auto cfl_backend = execution.makePDCfl();
    solid_dynamics::PDTimeStepExecution timestep(particle, cfl_backend.get(), options.numerics.cfl);

    SimpleDynamics<solid_dynamics::RigidBoxContactforPD> contact(
        particle, options.contact.friction_mu, options.contact.friction_velocity);
    contact.applyToCurrentAcceleration(true);
    const Vecd half_size(geometry.plate_half_width, 0.5 * geometry.plate_thickness, geometry.plate_half_width);
    contact.addBox("LowerPlate", Vecd(0, -0.5 * geometry.plate_thickness, 0), half_size);
    const peridynamics::CompressionPlateMotion motion{
        options.initialUpperFace(), geometry.plate_thickness, options.loading.stroke, options.duration()};
    const size_t upper_box = contact.addBox("UpperPlate", motion.center(0), half_size, motion.velocity(0));
    UniquePtr<gel_compression::SurfaceContact> surface_contact;
    if (options.contact.reference_surface)
        surface_contact = makeUnique<gel_compression::SurfaceContact>(particle, geometry.radius, options.centerHeight(),
            options.contact.surface_support, options.contact.penalty_length, options.contact.surface_rings,
            options.contact.surface_azimuths, options.contact.surface_polar_levels,
            options.contact.friction_mu, options.contact.friction_velocity,
            geometry.plate_half_width);
    auto *contact_sample_position = particles.registerStateVariableData<Vecd>("GelContactSamplePosition");
    auto *contact_sample_velocity = particles.registerStateVariableData<Vecd>("GelContactSampleVelocity");
    auto *contact_base_velocity = particles.registerStateVariableData<Vecd>("GelContactBaseVelocity");
    const auto *previous_acceleration = particles.getVariableDataByName<Vecd>("PDAccelerationOld");
    contact.setForceSamplingVelocityData(contact_sample_velocity);
    auto synchronize_plate = [&](Real time)
    {
        contact.setKinematics(upper_box, motion.center(time), motion.velocity(time));
        if (surface_contact)
            surface_contact->setUpper(options.initialUpperFace() - motion.displacement(time), motion.velocity(time));
    };

    SimpleDynamics<peridynamics::ParticleStateCheck> state_check(particle);
    state_check.addScalar("Mass", 0);
    state_check.addMatrix("DeformationGradient");
    state_check.addMatrix("PDStress");
    state_check.addVector("PDAcceleration");
    state_check.addVector("PDAccelerationPrior");
    for (const char *name : {"Mass", "Density", "VolumetricMeasure", "Damage"})
        particles.addVariableToWrite<Real>(name);
    for (const char *name : {"PDStress", "DeformationGradient"})
        particles.addVariableToWrite<Matd>(name);
    for (const char *name : {"Velocity", "UpperPlateReaction", "LowerPlateReaction", "ContactRepulsionFromWall", "ContactFrictionForce", "GelContactSamplePosition", "GelContactSampleVelocity", "GelContactBaseVelocity", "PDAccelerationOld"})
        particles.addVariableToWrite<Vecd>(name);
    for (const char *name : {"ContactBoundaryWork", "ContactFrictionLoss", "PDHourglassWork", "PDDampingWork"})
        particles.addVariableToWrite<Real>(name);
    if (material.isPlastic())
    {
        particles.addVariableToWrite<Matd>("PlasticStrain");
        particles.addVariableToWrite<Real>("IsotropicHardeningParam");
    }
    if (material.viscoelastic)
    {
        for (const char *name : {"GelGreenStrain", "GelViscousStrain", "GelAcceptedStress"})
            particles.addVariableToWrite<Matd>(name);
        for (const char *name : {"GelEnergyDensity", "GelDissipationDensity", "GelHistoryTime"})
            particles.addVariableToWrite<Real>(name);
    }
    BodyStatesRecordingToVtp states(system);
    states.setFullPrecisionOutput(true);

    Real &time = *system.getSystemVariableDataByName<Real>("PhysicalTime");
    size_t steps = 0;
    uint64_t completed_frame = 0;
    Real last_contact_time = 0;
    Real wall_offset = 0;
    Real initial_stored_energy = 0;
    const Real *mass = particles.getVariableDataByName<Real>("Mass");
    const Vecd *velocity = particles.getVariableDataByName<Vecd>("Velocity");
    const Vecd *position = particles.getVariableDataByName<Vecd>("Position");
    StdVec<Real> initial_masses(mass, mass + count);
    peridynamics::ContactImpulseLedger contact_impulses;
    const Vecd *lower_reaction = particles.getVariableDataByName<Vecd>("LowerPlateReaction");
    const Vecd *upper_reaction = particles.getVariableDataByName<Vecd>("UpperPlateReaction");
    const Real *reference_volume = particles.getVariableDataByName<Real>("VolumetricMeasure");
    const Real *accepted_energy = material.viscoelastic ? particles.getVariableDataByName<Real>("GelEnergyDensity") : nullptr;
    const Real *accepted_dissipation = material.viscoelastic ? particles.getVariableDataByName<Real>("GelDissipationDensity") : nullptr;
    auto momentum = [&]()
    {
        Vecd value = Vecd::Zero();
        for (size_t i = 0; i < count; ++i) value += mass[i] * velocity[i];
        return value;
    };

    // All registered material/integrator state and reference bonds are retained.
    // This schema is not a claim that new material memories have been validated.
    particles.getParticleGroupManager();
    const std::string identity = "gel-solid-sphere-v6-pd-model-options\n" + inputs.serialize() +
        "\nrecording=" + std::to_string(system.StateRecording()) +
        "\nexecutable=" + std::to_string(peridynamics::checkpoint_detail::fileHash(fs::canonical(argv[0])));
    peridynamics::DiskCheckpoint checkpoint(restart.directory, identity);
    checkpoint.addBody(dynamic_cast<peridynamics::CheckpointParticles &>(particles));
    checkpoint.addConfiguration("GelReference", reference.inner_configuration_, particles.ParticlesBound());
    checkpoint.addScalar("time", time);
    checkpoint.addScalar("steps", steps);
    checkpoint.addScalar("completed_frame", completed_frame);
    checkpoint.addScalar("last_contact_time", last_contact_time);
    checkpoint.addScalar("initial_stored_energy", initial_stored_energy);
    checkpoint.addArray("initial_masses", initial_masses);
    if (surface_contact) surface_contact->registerCheckpoint(checkpoint);
    contact_impulses.registerCheckpoint(checkpoint);
    checkpoint.addOutput("effective_parameters.cfg");
    checkpoint.addOutput("diagnostics.csv");
    checkpoint.addOutput("step_trace.csv");

    std::ofstream diagnostics, step_trace;
    auto flush = [&]()
    {
        diagnostics.flush(); requireGoodStream(diagnostics);
        step_trace.flush(); requireGoodStream(step_trace);
    };
    size_t last_checkpoint_step = 0;
    Real last_checkpoint_wall = 0;
    auto save_checkpoint = [&](bool force)
    {
        const Real now = wall_offset + limits.elapsed();
        const bool due = force || peridynamics::checkpoint_requested ||
            (restart.steps && steps - last_checkpoint_step >= restart.steps) ||
            (restart.seconds > 0 && now - last_checkpoint_wall >= restart.seconds);
        if (!due)
            return false;
        flush();
        checkpoint.save({time, now, steps, steps, completed_frame});
        last_checkpoint_step = steps;
        last_checkpoint_wall = now;
        peridynamics::checkpoint_requested = 0;
        return peridynamics::checkpoint_stop_requested ||
            (restart.stop_after_checkpoint && (!restart.resume.empty() || time > 0));
    };

    // This observer has no initialize/force/update calls. It records current
    // accepted geometry separately from the midpoint force actually integrated.
    auto observe = [&](uint64_t frame)
    {
        const auto boxes = surface_contact ? surface_contact->observe() : contact.observe();
        Real total_mass = 0, initial_mass = 0, mass_error = 0, kinetic = 0, kick_work = 0;
        const Vecd *applied = particles.getVariableDataByName<Vecd>("UpperPlateReaction");
        Vecd applied_total = Vecd::Zero();
        const Real *internal = particles.getVariableDataByName<Real>("PDInternalWork");
        const Real *external = particles.getVariableDataByName<Real>("PDExternalWork");
        const Real *hourglass_work = particles.getVariableDataByName<Real>("PDHourglassWork");
        const Real *damping_work = particles.getVariableDataByName<Real>("PDDampingWork");
        const Real *contact_work = particles.getVariableDataByName<Real>("ContactBoundaryWork");
        const Real *friction_loss = particles.getVariableDataByName<Real>("ContactFrictionLoss");
        const Real *plastic_work = particles.getVariableDataByName<Real>("PlasticWork");
        Real applied_contact_work = 0, applied_friction_loss = 0, plastic = 0, hg_work = 0, pd_damping_work = 0;
        Real gel_energy = 0, gel_dissipation = 0;
        const Real *volume = particles.getVariableDataByName<Real>("VolumetricMeasure");
        const Real *energy = material.viscoelastic ? particles.getVariableDataByName<Real>("GelEnergyDensity") : nullptr;
        const Real *dissipation = material.viscoelastic ? particles.getVariableDataByName<Real>("GelDissipationDensity") : nullptr;
        for (size_t i = 0; i < count; ++i)
        {
            total_mass += mass[i];
            initial_mass += initial_masses[i];
            mass_error = std::max(mass_error, std::abs(mass[i] / initial_masses[i] - 1));
            kinetic += 0.5 * mass[i] * velocity[i].squaredNorm();
            kick_work += internal[i] + external[i] + hourglass_work[i] + damping_work[i];
            applied_contact_work += contact_work[i];
            applied_friction_loss += friction_loss[i];
            plastic += plastic_work[i];
            hg_work += hourglass_work[i];
            pd_damping_work += damping_work[i];
            if (material.viscoelastic)
            {
                gel_energy += volume[i] * energy[i];
                gel_dissipation += volume[i] * dissipation[i];
            }
            applied_total += applied[i];
        }
        const Real stored = kinetic + gel_energy + boxes[0].potential_energy + boxes[upper_box].potential_energy;
        if (frame == 0)
            initial_stored_energy = stored;
        const Real residual = applied_contact_work - (stored - initial_stored_energy + gel_dissipation +
            applied_friction_loss - hg_work - pd_damping_work);
        diagnostics << time << ',' << motion.displacement(time) << ',' << steps << ',' << count << ','
                    << total_mass << ',' << total_mass / initial_mass - 1 << ',' << mass_error << ','
                    << boxes[upper_box].reaction[1] << ',' << applied_total[1] << ',' << last_contact_time << ','
                    << boxes[0].reaction[1] << ',' << boxes[upper_box].max_penetration << ','
                    << boxes[upper_box].contacting_particles << ',' << kinetic << ',' << kinetic - kick_work << ','
                    << boxes[0].potential_energy + boxes[upper_box].potential_energy << ','
                    << motion.center(time)[1] - half_size[1] << ','
                    << boxes[upper_box].max_penetration / geometry.spacing << ','
                    << contact.normalPenaltyModulus() << ',' << applied_contact_work << ',' << applied_friction_loss << ','
                    << plastic << ',' << hg_work << ',' << pd_damping_work << ',' << material.viscoelastic << ','
                    << gel_energy << ',' << gel_dissipation << ',' << residual << ',' << boxes[0].max_penetration << ','
                    << (boxes[0].max_penetration + boxes[upper_box].max_penetration) / options.loading.stroke << '\n';
        flush();
        states.writeToFile(frame);
    };

    try
    {
        system.initializeSystemCellLinkedLists();
        system.initializeSystemConfigurations();
        if (material.viscoelastic && restart.resume.empty())
            solid_dynamics::canonicalizeInitialGelReference(reference, time);
        shape.exec();
        if (restart.resume.empty())
        {
            synchronize_plate(time);
            if (surface_contact) surface_contact->referenceFile(false);
            contact_impulses.initialize(momentum());
            inputs.save("output/effective_parameters.cfg");
            diagnostics.open("output/diagnostics.csv");
            step_trace.open("output/step_trace.csv");
            step_trace << std::setprecision(17)
                << "time_s,dt_s,pd_steps,last_contact_sample_time_s,material_history_time_s,pd_wave_dt_s,configured_candidate_dt_s,output_limited_dt_s,contact_limited_dt_s,"
                   "lower_last_applied_Fx_N,lower_last_applied_Fy_N,lower_last_applied_Fz_N,"
                   "upper_last_applied_Fx_N,upper_last_applied_Fy_N,upper_last_applied_Fz_N,"
                   "lower_impulse_x_Ns,lower_impulse_y_Ns,lower_impulse_z_Ns,upper_impulse_x_Ns,upper_impulse_y_Ns,upper_impulse_z_Ns,"
                   "accepted_momentum_x_kg_m_per_s,accepted_momentum_y_kg_m_per_s,accepted_momentum_z_kg_m_per_s,"
                   "momentum_balance_residual_x_kg_m_per_s,momentum_balance_residual_y_kg_m_per_s,momentum_balance_residual_z_kg_m_per_s,"
                   "accepted_kinetic_J,accepted_gel_stored_J,accepted_gel_dissipation_J,accepted_material_energy_available\n";
            diagnostics << std::setprecision(17)
                << "time_s,stroke_m,pd_steps,solid_count,solid_mass_kg,mass_relative_error,particle_mass_error,"
                   "upper_observed_Fy_N,upper_last_applied_Fy_N,last_contact_sample_time_s,lower_observed_Fy_N,"
                   "upper_max_penetration_m,upper_contact_count,kinetic_J,kick_identity_residual_J,"
                   "observed_contact_potential_J,upper_face_y_m,upper_penetration_over_spacing,"
                   "contact_penalty_modulus_Pa,applied_contact_work_J,applied_friction_loss_J,plastic_work_J,"
                   "hourglass_work_J,pd_damping_work_J,material_energy_available,gel_stored_energy_J,"
                   "gel_dissipation_J,mechanical_residual_J,lower_max_penetration_m,total_penetration_over_total_stroke\n";
            observe(0);
        }
        else
        {
            const auto loaded = checkpoint.restore(restart.resume);
            if (time != loaded.time || steps != loaded.pd_steps || steps != loaded.acoustic_steps ||
                completed_frame != loaded.completed_frame || completed_frame > uint64_t(options.loading.output_frames) ||
                time < 0 || time > options.duration() || last_contact_time < 0 || last_contact_time > time ||
                !contact_impulses.matchesCursor(time, steps))
                throw std::runtime_error("inconsistent gel checkpoint cursor");
            wall_offset = loaded.wall_seconds;
            last_checkpoint_wall = wall_offset;
            last_checkpoint_step = steps;
            particle.updateCellLinkedList(); // Do not rebuild restored reference bonds.
            synchronize_plate(time);
            if (surface_contact) surface_contact->referenceFile(true);
            const std::regex snapshot("GelParticle_ite_([0-9]+)\\.vtp");
            const fs::path orphans = "output/orphaned-after-resume-" + std::to_string(steps) + "-" +
                std::to_string(std::chrono::steady_clock::now().time_since_epoch().count());
            for (const auto &entry : fs::directory_iterator("output"))
            {
                std::smatch match;
                const std::string name = entry.path().filename().string();
                if (entry.is_regular_file() && std::regex_match(name, match, snapshot) &&
                    std::stoull(match[1]) > completed_frame)
                {
                    fs::create_directories(orphans);
                    fs::rename(entry.path(), orphans / name);
                }
            }
            diagnostics.open("output/diagnostics.csv", std::ios::app);
            diagnostics << std::setprecision(17);
            step_trace.open("output/step_trace.csv", std::ios::app);
            step_trace << std::setprecision(17);
        }
        if (save_checkpoint(true))
            return 75;

        for (uint64_t frame = completed_frame + 1; frame <= uint64_t(options.loading.output_frames); ++frame)
        {
            const Real target = frame == uint64_t(options.loading.output_frames) ? options.duration() :
                options.duration() * Real(frame) / options.loading.output_frames;
            while (time < target)
            {
                limits.check(steps, wall_offset);
                const Real wave_dt = timestep.exec();
                Real candidate = material.viscoelastic ? std::min(wave_dt, 0.1 * material.relaxation_time) : wave_dt;
                if (options.numerics.maximum_dt > 0) candidate = std::min(candidate, options.numerics.maximum_dt);
                Real dt = peridynamics::boundedTimeStep(candidate, {target - time});
                const Real output_limited_dt = dt;
                if (surface_contact)
                {
                    bool stable = false;
                    for (int iteration = 0; iteration < 16; ++iteration)
                    {
                        const Real restricted = surface_contact->stableStep(dt,
                            options.initialUpperFace() - motion.displacement(time + 0.5 * dt), options.numerics.cfl);
                        if (restricted >= dt) { stable = true; break; }
                        dt = peridynamics::boundedTimeStep(restricted, {target - time});
                    }
                    if (!stable) throw std::runtime_error("surface contact timestep restriction did not converge");
                }
                initialize.exec();
                first.exec(dt);
                last_contact_time = time + 0.5 * dt;
                synchronize_plate(last_contact_time);
                for (size_t i = 0; i < count; ++i)
                {
                    contact_sample_position[i] = position[i];
                    contact_base_velocity[i] = velocity[i];
                    contact_sample_velocity[i] = velocity[i] + 0.5 * dt * previous_acceleration[i];
                }
                if (surface_contact) surface_contact->apply(contact_sample_velocity, last_contact_time);
                else contact.exec();
                constitutive.force->exec(dt);
                if (fracture) fracture->exec(dt);
                hourglass.exec(dt);
                damping.exec(dt);
                third.exec(dt);
                fourth.exec(dt);
                if (constitutive.commit)
                    constitutive.commit->exec(dt);
                ++steps;
                state_check.exec();
                for (size_t i = 0; i < count; ++i)
                    if (position[i].cwiseAbs().maxCoeff() > extent)
                        peridynamics::numericalFailure(particle.Name(), i, "particle left bounded gel domain");
                time = dt == target - time ? target : time + dt;
                synchronize_plate(time);
                const Real history_time = material.viscoelastic ? particles.getVariableDataByName<Real>("GelHistoryTime")[0] : time;
                Vecd lower_force = Vecd::Zero(), upper_force = Vecd::Zero();
                Real accepted_kinetic = 0, accepted_gel = 0, accepted_loss = 0;
                for (size_t i = 0; i < count; ++i)
                {
                    lower_force += lower_reaction[i]; upper_force += upper_reaction[i];
                    accepted_kinetic += 0.5 * mass[i] * velocity[i].squaredNorm();
                    if (material.viscoelastic)
                    {
                        accepted_gel += reference_volume[i] * accepted_energy[i];
                        accepted_loss += reference_volume[i] * accepted_dissipation[i];
                    }
                }
                contact_impulses.advance(time, dt, last_contact_time, lower_force, upper_force);
                const Vecd accepted_momentum = momentum();
                const Vecd momentum_balance = contact_impulses.momentumResidual(accepted_momentum);
                step_trace << time << ',' << dt << ',' << steps << ',' << last_contact_time << ',' << history_time << ','
                           << wave_dt << ',' << candidate << ',' << output_limited_dt << ',' << dt;
                for (const auto &value : {lower_force, upper_force, contact_impulses.lower(), contact_impulses.upper(), accepted_momentum, momentum_balance})
                    for (int d = 0; d < Dimensions; ++d) step_trace << ',' << value[d];
                step_trace << ',' << accepted_kinetic << ',' << accepted_gel << ',' << accepted_loss << ',' << material.viscoelastic << '\n';
                if (save_checkpoint(false))
                    return 75;
            }
            observe(frame);
            completed_frame = frame;
        }
        if (restart.final_checkpoint && save_checkpoint(true))
            return 75;
        std::ofstream("output/run_status.txt") << "completed\ntime_s=" << std::setprecision(17) << time << '\n';
        return 0;
    }
    catch (const std::exception &error)
    {
        const auto *resource = dynamic_cast<const ResourceStop *>(&error);
        std::ofstream failure("output/run_status.txt");
        failure << "failed\ntime_s=" << std::setprecision(17) << time << "\npd_steps=" << steps
                << "\nreason=" << error.what() << "\nreturn_code=" << (resource ? resource->return_code : 1) << '\n';
        throw;
    }
}
} // namespace

int main(int argc, char **argv)
{
    try
    {
        const bool describe = consumeDescriptionFlag(argc, argv);
        const auto execution_mode = peridynamics::takePeriExecutionOption(argc, argv);
        peridynamics::ParticleCompressionOptions options;
        peridynamics::CheckpointOptions restart;
        restart.threads = std::max(size_t(1), std::min(restart.threads, ExecutionLimits::max_threads));
        restart.parse(argc, argv);
        auto inputs = options.bindings();
        inputs.parse(argc, argv);
        options.resolveModels(inputs);
        options.validate();
        if (restart.threads > ExecutionLimits::max_threads)
            throw std::invalid_argument("bounded gel control permits at most eight threads");
        if (describe || restart.describe_only)
        {
            std::cout << std::setprecision(17)
                      << "{\"model\":\"solid_sphere_no_separate_liquid\",\"calibrated\":false,"
                         "\"fluid_transport\":false,\"fracture\":"
                      << (options.fracture.model != "none" ? "true" : "false") << ",\"constitutive_model\":\""
                      << (options.material.viscoelastic ? "neo_hookean_green_maxwell" :
                          (options.material.model == "viscoplastic" ? "perzyna_j2_viscoplastic" :
                           (options.material.model == "j2" ? "j2_plastic" :
                            (options.material.model == "drucker_prager" ? "drucker_prager" : "objective_linear_elastic"))))
                      << "\",\"material_model\":\"" << options.material.model
                      << "\",\"hardening_model\":\"" << options.material.hardening_model
                      << "\",\"fracture_model\":\"" << options.fracture.model
                      << "\",\"contact_model\":\"" << (options.contact.reference_surface ? "reference_surface_mls" : "particle_radius") << "\","
                         "\"format\":1,\"full_pd_state\":true,\"default_threads\":" << restart.threads
                      << ",\"duration_s\":" << options.duration() << ",\"spacing_m\":" << options.geometry.spacing
                      << ",\"max_particles\":" << ExecutionLimits::max_particles
                      << ",\"execution_profile\":\"solid-particle-spatial-v1\",\"max_steps\":" << ExecutionLimits::max_steps
                      << ",\"max_threads\":" << ExecutionLimits::max_threads
                      << ",\"contact_trace_format\":1,\"restart_qualified\":false}\n";
            return 0;
        }
        return simulate(argc, argv, options, inputs, restart, peridynamics::PeriExecution(execution_mode));
    }
    catch (const ResourceStop &error)
    {
        std::cerr << "RESOURCE_STOP: " << error.what() << '\n';
        return error.return_code;
    }
    catch (const std::exception &error)
    {
        std::cerr << "GEL_CONTROL_FAILED: " << error.what() << '\n';
        return 1;
    }
}
