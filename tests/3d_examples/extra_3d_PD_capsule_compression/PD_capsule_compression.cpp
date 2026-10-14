/** Single capsule between prescribed rigid plates. SI units throughout.
 * PD shell + WCSPH core. Parameters and numerical diagnostics are saved per run.
 * Reference PD neighborhoods remain fixed; contact neighborhoods evolve.
 */
#include "sphinxsys.h"
#include "pd_case_config.h"
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>

using namespace SPH;

Real capsule_radius_out = 0.5e-3, shell_thickness = 0.02e-3, shell_layers = 3;
Real plate_radius = 0.75e-3, plate_thickness = 0.06e-3, initial_gap = 1e-5;
Real rho0_shell = 1000, youngs_shell = 1e7, poisson_shell = 0.4;
Real max_tension_stress = 1e5, yield_stress_shell = 5e5, hardening_shell = 5e6;
Real viscosity_shell = 500, critical_plastic_strain = 1;
Real plate_friction_mu = 0.3, plate_friction_veps = 0.005;
Real rho0_f = 935, c_f = 25, mu_f = 0.25, fluid_spacing_ratio = 1.5;
Real gravity_g = 9.8, compression_speed = 0.1, compression_distance = 0.58e-3;
Real U_ref = 0.5, pd_cfl = 0.1, acoustic_cfl = 0.3, acoustic_damping_factor = 0.3;
Real max_density_error = 0.1, max_fluid_mach = 0.3;
bool shell_viscoplastic = true, with_liquid = true, reinitialize_density = false;
int output_frames = 50, bond_check_interval = 1, inject_nan_step = -1, max_pd_substeps = 6;
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
    BIND(reinitialize_density); BIND(U_ref); BIND(pd_cfl); BIND(acoustic_cfl); BIND(acoustic_damping_factor);
    BIND(max_density_error); BIND(max_fluid_mach); BIND(shell_viscoplastic); BIND(with_liquid);
    BIND(max_pd_substeps); BIND(output_frames); BIND(bond_check_interval); BIND(inject_nan_step);
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
        output_frames < 1 || bond_check_interval < 1 || max_pd_substeps < 1 || pd_cfl <= 0 || pd_cfl > 0.2 ||
        acoustic_cfl <= 0 || acoustic_cfl > 0.6 || acoustic_damping_factor < 0 ||
        max_density_error <= 0 || max_density_error >= 1 || max_fluid_mach <= 0 ||
        plate_friction_mu < 0 || plate_friction_veps <= 0)
        throw std::invalid_argument("invalid capsule geometry, material, schedule or numerical controls");
    capsule_radius_in = capsule_radius_out - shell_thickness;
    resolution_ref = shell_thickness / shell_layers;
    capsule_center_y = plate_thickness + initial_gap + capsule_radius_out;
    upper_plate_bottom = capsule_center_y + capsule_radius_out + initial_gap;
    end_time = compression_distance / compression_speed;
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
    auto config = parameters();
    try {config.parse(argc, argv); validateParameters();}
    catch (const std::exception &e) {std::cerr << e.what() << std::endl; return 2;}
    Real extent = 6 * capsule_radius_out;
    SPHSystem system(BoundingBoxd(Vecd(-extent, -extent, -extent), Vecd(extent, 2 * extent, extent)), resolution_ref);
    system.handleCommandlineOptions(argc, argv);
    if (system.RestartStep() != 0)
    {
        std::cerr << "PD restart is unsupported until reference bonds and all internal variables are serialized.\n";
        return 2;
    }
    std::filesystem::create_directories("output");
    config.save("output/effective_parameters.cfg");

    FluidBody liquid(system, makeShared<GeometricShapeBall>(Vecd(0, capsule_center_y, 0), capsule_radius_in, "LiquidCore"));
    liquid.defineAdaptationRatios(1.3, 1 / fluid_spacing_ratio);
    liquid.defineMatterMaterial<WeaklyCompressibleFluid>(rho0_f, c_f);
    liquid.addMaterialProperty<Viscosity>(mu_f);
    liquid.generateParticles<BaseParticles, Lattice>();
    PDBody shell(system, makeShared<CapsuleShellShape>("CapsuleShell"));
    if (shell_viscoplastic)
        shell.defineMatterMaterial<ViscoPlasticSolidforPD>(rho0_shell, youngs_shell, poisson_shell, yield_stress_shell, viscosity_shell, hardening_shell);
    else shell.defineMatterMaterial<HughesWingetSolid>(rho0_shell, youngs_shell, poisson_shell);
    shell.generateParticles<BaseParticles, Lattice>();
    SolidBody lower(system, makeShared<PlateShape>("LowerPlate", 0.5 * plate_thickness));
    lower.defineAdaptationRatios(1.3, 1.0 / 3.0);
    lower.defineMatterMaterial<Solid>(rho0_shell, youngs_shell);
    lower.generateParticles<BaseParticles, Lattice>();
    SolidBody upper(system, makeShared<PlateShape>("UpperPlate", upper_plate_bottom + 0.5 * plate_thickness));
    upper.defineAdaptationRatios(1.3, 1.0 / 3.0);
    upper.defineMatterMaterial<Solid>(rho0_shell, youngs_shell);
    upper.generateParticles<BaseParticles, Lattice>();

    InnerRelation fluid_inner(liquid), shell_inner(shell);
    ContactRelationSelectKernel fluid_walls(liquid, {&shell, &lower, &upper}, true);
    ContactRelationSelectKernel shell_fluid(shell, {&liquid}, false);
    ContactRelationSelectKernel upper_fluid(upper, {&liquid}, false), lower_fluid(lower, {&liquid}, false);
    ContactRelation shell_plates(shell, {&lower, &upper}), lower_shell(lower, {&shell}), upper_shell(upper, {&shell});
    ComplexRelation fluid_complex(fluid_inner, fluid_walls);

    SimpleDynamics<NormalDirectionFromBodyShape> shell_normals(shell), lower_normals(lower), upper_normals(upper);
    Dynamics1Level<fluid_dynamics::Integration1stHalfWithWallRiemann> pressure(fluid_inner, fluid_walls);
    Dynamics1Level<fluid_dynamics::Integration2ndHalfWithWallRiemann> density(fluid_inner, fluid_walls);
    InteractionWithUpdate<fluid_dynamics::DensitySummationComplexFreeSurface> density_sum(fluid_inner, fluid_walls);
    InteractionWithUpdate<fluid_dynamics::ViscousForceWithWall> fluid_viscosity(fluid_inner, fluid_walls);
    SimpleDynamics<fluid_dynamics::ContinuumVolumeUpdate> update_fluid_volume(liquid);
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
    InteractionWithUpdate<solid_dynamics::ViscousForceFromFluid> shell_viscosity(shell_fluid);
    InteractionWithUpdate<solid_dynamics::PressureForceFromFluid<decltype(density)>> upper_pressure(upper_fluid), lower_pressure(lower_fluid);
    InteractionWithUpdate<solid_dynamics::ViscousForceFromFluid> upper_viscosity(upper_fluid), lower_viscosity(lower_fluid);

    InteractionDynamics<solid_dynamics::ContactRepulsionFactorforPD> contact_factor(shell_plates);
    InteractionDynamics<solid_dynamics::ContactRepulsionFromWallforPD> contact_force(shell_plates, plate_friction_mu, plate_friction_veps);
    InteractionDynamics<solid_dynamics::ContactRepulsionToWallforPD> lower_reaction(lower_shell, plate_friction_mu, plate_friction_veps);
    InteractionDynamics<solid_dynamics::ContactRepulsionToWallforPD> upper_reaction(upper_shell, plate_friction_mu, plate_friction_veps);
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
    InteractionDynamics<solid_dynamics::LittleWoodHourGlassControl> hourglass(shell_inner);
    InteractionDynamics<solid_dynamics::PairNumericalDampingforPD> pd_damping(shell_inner);
    InteractionDynamics<solid_dynamics::NosbPDThirdStep> third(shell_inner);
    SimpleDynamics<solid_dynamics::NosbPDFourthStep> fourth(shell);
    ReduceDynamics<solid_dynamics::PDTimeStep> shell_step(shell, pd_cfl);
    SimpleDynamics<solid_dynamics::UpdateElasticNormalDirection> update_normals(shell);
    solid_dynamics::AverageVelocityAndAcceleration averages(shell);
    SimpleDynamics<solid_dynamics::PrescribedTranslation> move_upper(upper, 0, end_time, Vecd(0, -compression_distance, 0));

    SimpleDynamics<peridynamics::ParticleStateCheck> fluid_check(liquid), shell_check(shell);
    fluid_check.addScalar("Density", 0); fluid_check.addScalar("Pressure");
    fluid_check.addVector("Force"); fluid_check.addVector("ForcePrior");
    shell_check.addMatrix("DeformationGradient"); shell_check.addMatrix("PDStress");
    shell_check.addVector("PDAcceleration"); shell_check.addVector("PDAccelerationPrior");
    for (const auto &name : {"Mass", "Density", "Pressure"}) liquid.getBaseParticles().addVariableToWrite<Real>(name);
    for (const auto &name : {"Velocity", "AverageVelocity", "AverageAcceleration", "PDAcceleration", "PDAccelerationPrior"}) shell.getBaseParticles().addVariableToWrite<Vecd>(name);
    for (const auto &name : {"PDStress", "DeformationGradient"}) shell.getBaseParticles().addVariableToWrite<Matd>(name);
    if (shell_viscoplastic)
    {
        shell.getBaseParticles().addVariableToWrite<Matd>("PlasticStrain");
        shell.getBaseParticles().addVariableToWrite<Real>("IsotropicHardeningParam");
    }
    BodyStatesRecordingToVtp states(system);
    Real &time = *system.getSystemVariableDataByName<Real>("PhysicalTime");
    Real pd_dt = pd_cfl * resolution_ref / DynamicCast<HughesWingetSolid>(&shell, shell.getMatterMaterial()).ReferenceSoundSpeed();
    std::ofstream diagnostics("output/diagnostics.csv"), reaction("output/UpperPlate_TotalReaction.dat");
    diagnostics << std::setprecision(17) << "time_s,stroke_m,acoustic_steps,pd_steps,liquid_mass_kg,fluid_speed_max_m_s,fluid_density_error_max,damage_mean,damage_max,numerical_failures,shell_kinetic_J,liquid_kinetic_J,contact_Fy_N,fluid_Fy_N,total_Fy_N\n";
    reaction << std::setprecision(17) << "time_s stroke_m contact_Fy_N fluid_Fy_N total_Fy_N\n";
    size_t acoustic_count = 0, pd_count = 0;
    std::string phase = "initialization";

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
    auto update_contacts = [&]()
    {
        shell.updateCellLinkedList(); upper.updateCellLinkedList();
        if (with_liquid)
        {
            liquid.updateCellLinkedList(); fluid_complex.updateConfiguration(); shell_fluid.updateConfiguration();
            upper_fluid.updateConfiguration(); lower_fluid.updateConfiguration();
        }
        shell_plates.updateConfiguration(); lower_shell.updateConfiguration(); upper_shell.updateConfiguration();
    };
    auto measure = [&](int frame)
    {
        contact_factor.exec(); initialize_shell.exec(); contact_force.exec(); upper_reaction.exec(); lower_reaction.exec();
        Real contact_y = sumVector(upper, "ContactRepulsionToWall")[1], fluid_y = 0;
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
        Real stroke = std::min(time / end_time, Real(1)) * compression_distance;
        diagnostics << time << ',' << stroke << ',' << acoustic_count << ',' << pd_count << ',' << mass << ',' << vmax << ',' << rho_error << ',' << dsum / sp.TotalRealParticles() << ',' << dmax << ',' << failures << ',' << shell_ke << ',' << fluid_ke << ',' << contact_y << ',' << fluid_y << ',' << contact_y + fluid_y << '\n';
        reaction << time << ' ' << stroke << ' ' << contact_y << ' ' << fluid_y << ' ' << contact_y + fluid_y << '\n';
        diagnostics.flush(); reaction.flush(); states.writeToFile(frame);
        std::cout << "frame=" << frame << " time=" << time << " stroke=" << stroke << " F=" << contact_y + fluid_y << " rho_error=" << rho_error << " vmax=" << vmax << std::endl;
    };

    try
    {
        system.initializeSystemCellLinkedLists(); system.initializeSystemConfigurations();
        shell_normals.exec(); lower_normals.exec(); upper_normals.exec(); shape_matrix.exec(); move_upper.exec();
        if (with_liquid) fluid_gravity.exec();
        shell_check.exec(); check_fluid(); measure(0);
        std::cout << "particles shell=" << shell.getBaseParticles().TotalRealParticles()
                  << " liquid=" << (with_liquid ? liquid.getBaseParticles().TotalRealParticles() : 0)
                  << " PD_dt=" << pd_dt << " end_time=" << end_time << std::endl;
        for (int frame = 1; frame <= output_frames; ++frame)
        {
            Real target_time = frame == output_frames ? end_time : end_time * Real(frame) / output_frames;
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
                        fluid_viscosity.exec(); divergence.exec(); wall_divergence.exec();
                        pressure.exec(dt); inner_damping.exec(dt); wall_damping.exec(dt); check_fluid();
                    }
                    phase = "external shell loads";
                    initialize_shell.exec();
                    if (with_liquid)
                    {
                        shell_pressure.exec(dt); shell_viscosity.exec();
                        auto &p = shell.getBaseParticles();
                        auto *a = p.getVariableDataByName<Vecd>("PDAccelerationPrior");
                        auto *f = p.getVariableDataByName<Vecd>("ViscousForceFromFluid");
                        auto *v = p.getVariableDataByName<Real>("VolumetricMeasure");
                        for (size_t i = 0; i < p.TotalRealParticles(); ++i) a[i] += f[i] / (rho0_shell * v[i]);
                    }
                    contact_factor.exec(); contact_force.exec();
                    phase = "fluid second half";
                    if (with_liquid) {density.exec(dt); update_fluid_volume.exec(); check_fluid();}
                    phase = "PD subcycling";
                    Real elapsed = 0; averages.initialize_displacement_.exec();
                    while (elapsed < dt)
                    {
                        Real ds = peridynamics::boundedTimeStep(shell_step.exec(), {dt - elapsed});
                        first.exec(ds); second->exec(ds);
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
                }
                if (advection_end == target_time) time = target_time;
                phase = "contact neighborhood refresh"; update_contacts();
            }
            phase = "output"; measure(frame);
        }
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
