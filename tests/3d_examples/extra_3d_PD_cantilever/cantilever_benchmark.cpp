/** Reproducible small-strain cantilever verification; see BENCHMARK.md.
 * The legacy operator is a benchmark-only transcription, not a production option.
 */
#include "sphinxsys.h"
#include "pd_case_config.h"
#include <chrono>
#include <fstream>
#include <iomanip>
#include <tbb/global_control.h>

using namespace SPH;

namespace SPH
{
class CantileverCells {};
template <> class ParticleGenerator<BaseParticles, CantileverCells> : public ParticleGenerator<BaseParticles>
{
    int n_, nx_, ns_;
    Real dp_;
  public:
    ParticleGenerator(SPHBody &body, BaseParticles &particles, int n, int nx, int ns, Real dp)
        : ParticleGenerator<BaseParticles>(body, particles), n_(n), nx_(nx), ns_(ns), dp_(dp) {}
    void prepareGeometricData() override
    {
        // Exact integer-indexed cells ensure identical topology under rescaling.
        for (int ix = -ns_; ix < nx_; ++ix)
            for (int iy = 0; iy < n_; ++iy)
                for (int iz = 0; iz < n_; ++iz)
                    addPositionAndVolumetricMeasure(Vecd((ix + 0.5) * dp_,
                        (iy + 0.5 - 0.5 * n_) * dp_, (iz + 0.5 - 0.5 * n_) * dp_), dp_ * dp_ * dp_);
    }
};
}

class LegacyHourglass : public solid_dynamics::ProjectedPairStabilization
{
  public:
    LegacyHourglass(BaseInnerRelation &relation, Real coefficient)
        : ProjectedPairStabilization(relation, coefficient) {}
    void interaction(size_t i, Real = 0)
    {
        // Icebreaking_PD-SPH elastic_dynamics.cpp:508-537, original delta^-3
        // and one-endpoint F_i. Additional diagnostics do not alter the force.
        Vecd density = Vecd::Zero();
        const auto &neighbors = inner_configuration_[i];
        Real constant = 12.0 * Kbulk_ / (Pi * std::pow(horizon_, 3));
        if (particleLive_[i] == 1)
            for (size_t n = 0; n < neighbors.current_size_; ++n)
            {
                if (!neighbors.bondLive_[n]) continue;
                size_t j = neighbors.j_[n];
                Vecd reference = -neighbors.r_ij_[n] * neighbors.e_ij_[n];
                Vecd current = pos_[j] - pos_[i];
                if (current.norm() <= SqrtEps * reference.norm())
                    throw std::runtime_error("collapsed legacy bond");
                Vecd error = pos_[i] + F_[i] * reference - pos_[j];
                density -= Chg_ * constant * error.dot(current) / reference.norm() /
                           current.norm() * Vol_[j] * current;
            }
        acc_[i] += density / rho0_;
        hourglass_force_[i] = density * Vol_[i];
    }
};

namespace
{
constexpr Real beta = 1.875104068711961;
Real rawMode(Real x)
{
    const Real sigma = (std::cosh(beta) + std::cos(beta)) / (std::sinh(beta) + std::sin(beta));
    return std::cosh(beta * x) - std::cos(beta * x) - sigma * (std::sinh(beta * x) - std::sin(beta * x));
}
Real mode(Real x) { return rawMode(x) / rawMode(1); }
Real modeDerivative(Real x)
{
    const Real sigma = (std::cosh(beta) + std::cos(beta)) / (std::sinh(beta) + std::sin(beta));
    return beta * (std::sinh(beta * x) + std::sin(beta * x) - sigma * (std::cosh(beta * x) - std::cos(beta * x))) / rawMode(1);
}
int run(int argc, char **argv)
{
    int cells = 4, length_ratio = 12, stabilization = 2, load_mode = 0, samples_per_cycle = 160, threads = 2;
    Real scale = 1, cycles = 2, cfl = 0.2, coefficient = 0.5, horizon_ratio = 3.015;
    Real amplitude_ratio = 0.001, young = 1.0e7, rho = 1000, poisson = 0.3, relaxation_ratio = 1;
    bool numerical_damping = true;
    peridynamics::CaseParameters config;
    config.bind("cells", cells); config.bind("length_ratio", length_ratio);
    config.bind("stabilization", stabilization); config.bind("load_mode", load_mode);
    config.bind("samples_per_cycle", samples_per_cycle); config.bind("threads", threads);
    config.bind("scale", scale); config.bind("cycles", cycles); config.bind("cfl", cfl);
    config.bind("coefficient", coefficient); config.bind("horizon_ratio", horizon_ratio);
    config.bind("amplitude_ratio", amplitude_ratio); config.bind("young", young);
    config.bind("rho", rho); config.bind("poisson", poisson);
    config.bind("relaxation_ratio", relaxation_ratio); config.bind("numerical_damping", numerical_damping);
    config.parse(argc, argv);
    if (cells < 3 || length_ratio < 6 || stabilization < 0 || stabilization > 2 || load_mode < 0 || load_mode > 1 ||
        scale <= 0 || cycles <= 0 || cfl <= 0 || cfl > 0.5 || coefficient < 0 || horizon_ratio < 2.01 ||
        amplitude_ratio <= 0 || amplitude_ratio > 0.01 || young <= 0 || rho <= 0 || poisson <= -1 || poisson >= 0.5 ||
        samples_per_cycle < 40 || threads < 1 || relaxation_ratio <= 0)
        throw std::invalid_argument("invalid cantilever benchmark parameters");
    config.save("effective.cfg");
    tbb::global_control thread_limit(tbb::global_control::max_allowed_parallelism, threads);
    const Real H = scale, L = length_ratio * H, dp = H / cells;
    const int nx = length_ratio * cells, ns = std::max(cells, int(std::ceil(horizon_ratio)));
    const Real support = ns * dp, horizon = horizon_ratio * dp;
    const Real area = H * H, inertia = std::pow(H, 4) / 12;
    const Real omega = beta * beta * std::sqrt(young * inertia / (rho * area)) / (L * L);
    const Real period = 2 * Pi / omega, amplitude = amplitude_ratio * H;
    const Real load_per_length = 8 * young * inertia * amplitude / std::pow(L, 4);
    const Real gravity = load_mode == 1 ? load_per_length / (rho * area) : 0;
    const Real shear_modulus = young / (2 * (1 + poisson)), kappa = 5.0 / 6.0;
    const Real damping_rate = load_mode == 1 ? 2 * relaxation_ratio * omega : 0;
    const Real margin = 2 * horizon;
    SPHSystem system(BoundingBoxd(Vecd(-support - margin, -H - margin, -H - margin),
                                 Vecd(L + margin, H + margin, H + margin)), dp);
    system.handleCommandlineOptions(argc, argv);
    auto shape = makeShared<GeometricShapeBox>(Transform(Vecd((L - support) / 2, 0, 0)),
                                               Vecd((L + support) / 2, H / 2, H / 2), "Cantilever");
    PDBody body(system, shape);
    body.defineAdaptation<PDAdaptation>(horizon_ratio / 2);
    body.defineMatterMaterial<HughesWingetSolid>(rho, young, poisson);
    body.generateParticles<BaseParticles, CantileverCells>(cells, nx, ns, dp);
    InnerRelation relation(body);
    SimpleDynamics<solid_dynamics::PDTimeStepInitialization> initialize(body, makeShared<Gravity>(Vecd::Zero()));
    InteractionWithUpdate<solid_dynamics::NosbPDShapeMatrix> shape_matrix(relation);
    ReduceDynamics<solid_dynamics::PDTimeStep> time_step(body, cfl);
    SimpleDynamics<solid_dynamics::NosbPDFirstStep> first(body);
    InteractionWithUpdate<solid_dynamics::NosbPDSecondStep> second(relation);
    InteractionDynamics<solid_dynamics::ProjectedPairStabilization> modern(relation, stabilization == 0 ? 0 : coefficient);
    InteractionDynamics<LegacyHourglass> legacy(relation, coefficient);
    InteractionDynamics<solid_dynamics::PairNumericalDampingforPD> damping(relation);
    InteractionDynamics<solid_dynamics::NosbPDThirdStep> third(relation);
    SimpleDynamics<solid_dynamics::NosbPDFourthStep> fourth(body);
    system.initializeSystemCellLinkedLists();
    system.initializeSystemConfigurations();
    shape_matrix.exec();
    auto &particles = body.getBaseParticles();
    size_t count = particles.TotalRealParticles();
    auto pos = particles.getVariableDataByName<Vecd>("Position");
    auto vel = particles.getVariableDataByName<Vecd>("Velocity");
    auto acc = particles.getVariableDataByName<Vecd>("PDAcceleration");
    auto prior = particles.getVariableDataByName<Vecd>("PDAccelerationPrior");
    auto mass = particles.getVariableDataByName<Real>("Mass");
    auto F = particles.getVariableDataByName<Matd>("DeformationGradient");
    auto hg = particles.getVariableDataByName<Vecd>("HourglassForce");
    auto internal = particles.getVariableDataByName<Vecd>("PDInternalForce");
    auto stress = particles.getVariableDataByName<Matd>("PDStress");
    auto live = particles.getVariableDataByName<int>("ParticleLive");
    auto damage = particles.getVariableDataByName<Real>("Damage");
    auto failure = particles.getVariableDataByName<int>("NumericalFailureReason");
    StdVec<Vecd> reference(pos, pos + count);
    StdVec<Vecd> force_positions(reference);
    std::vector<Real> phi(count, 0);
    Real mode_norm = 0, free_mass = 0;
    size_t bonds = 0;
    for (size_t i = 0; i < count; ++i)
    {
        bonds += relation.inner_configuration_[i].current_size_;
        if (reference[i][0] <= 0) continue;
        phi[i] = mode(reference[i][0] / L);
        mode_norm += mass[i] * phi[i] * phi[i]; free_mass += mass[i];
        if (load_mode == 0)
            vel[i] = amplitude * omega * Vecd(-reference[i][1] * modeDerivative(reference[i][0] / L) / L, phi[i], 0);
    }
    std::ofstream meta("metadata.json");
    meta << std::setprecision(17) << "{\n\"particles\":" << count << ",\n\"directed_bonds\":" << bonds
         << ",\n\"height\":" << H << ",\n\"length\":" << L << ",\n\"spacing\":" << dp
         << ",\n\"horizon\":" << horizon << ",\n\"support_length\":" << support
         << ",\n\"free_mass\":" << free_mass << ",\n\"omega_eb\":" << omega << ",\n\"period_eb\":" << period
         << ",\n\"amplitude\":" << amplitude << ",\n\"load_per_length\":" << load_per_length
         << ",\n\"static_tip_eb\":" << amplitude << ",\n\"static_tip_timoshenko\":"
         << amplitude + load_per_length * L * L / (2 * kappa * shear_modulus * area)
         << ",\n\"tip_sample_x\":" << L - dp / 2 << "\n}\n";
    meta.close();
    std::ofstream history("history.csv"), profiles("profiles.csv");
    history << std::setprecision(17);
    profiles << std::setprecision(17);
    history << "time,t_over_T,modal_over_H,tip_over_H,tip_eb_over_H,kinetic,elastic_small_strain,drag_loss,clamp_drift_over_H,max_F_minus_I,hg_net_relative,hg_torque_relative,reaction_y,internal_net_relative,live_count,damaged_count,failure_count,steps,free_acc_rms_normalized,free_vel_rms_normalized\n";
    profiles << "t_over_T,x_over_L,mean_uy_over_H,std_uy_over_H,eb_over_H,timoshenko_static_over_H\n";
    Real &time = *system.getSystemVariableDataByName<Real>("PhysicalTime");
    size_t steps = 0;
    Real drag_loss = 0, max_dt = 0, min_dt = std::numeric_limits<Real>::max(), reaction = 0;
    const auto start = std::chrono::steady_clock::now();
    auto damp_velocities = [&](Real half_dt)
    {
        if (damping_rate == 0) return;
        Real factor = std::exp(-damping_rate * half_dt);
        for (size_t i = 0; i < count; ++i)
            if (reference[i][0] > 0)
            {
                drag_loss += 0.5 * mass[i] * vel[i].squaredNorm() * (1 - factor * factor);
                vel[i] *= factor;
            }
    };
    auto record = [&]()
    {
        Real modal = 0, kinetic = 0, elastic = 0, clamp_error = 0, max_F = 0, hg_norm = 0, torque_norm = 0, int_norm = 0, acceleration_sq = 0;
        Vecd hg_net = Vecd::Zero(), torque = Vecd::Zero(), int_net = Vecd::Zero();
        int alive_count = 0, damaged_count = 0, failures = 0;
        std::vector<Real> sum(nx, 0), sum2(nx, 0);
        for (size_t i = 0; i < count; ++i)
        {
            if (!pos[i].allFinite() || !vel[i].allFinite() || !F[i].allFinite()) throw std::runtime_error("non-finite state");
            const Vecd displacement = pos[i] - reference[i];
            max_F = std::max(max_F, (F[i] - Matd::Identity()).norm());
            hg_net += hg[i]; hg_norm += hg[i].norm();
            // Use the force evaluation positions, not the post-drift positions.
            torque += force_positions[i].cross(hg[i]); torque_norm += force_positions[i].norm() * hg[i].norm();
            int_net += internal[i]; int_norm += internal[i].norm();
            alive_count += live[i] == 1; damaged_count += damage[i] > 0; failures += failure[i] != 0;
            if (reference[i][0] <= 0) { clamp_error = std::max(clamp_error, displacement.norm()); continue; }
            modal += mass[i] * phi[i] * displacement[1];
            kinetic += 0.5 * mass[i] * vel[i].squaredNorm();
            acceleration_sq += mass[i] * acc[i].squaredNorm();
            // Diagnostic only: small-strain elastic energy excludes stabilization energy.
            Matd strain = 0.5 * (F[i] + F[i].transpose()) - Matd::Identity();
            elastic += 0.5 * mass[i] / rho * (stress[i].array() * strain.array()).sum();
            int ix = int(std::floor(reference[i][0] / dp));
            sum[ix] += displacement[1]; sum2[ix] += displacement[1] * displacement[1];
        }
        if (alive_count != int(count) || damaged_count || failures || max_F > 0.05 || clamp_error > 1e-12 * H)
            throw std::runtime_error("benchmark left intact, fixed-clamp, small-strain validity envelope");
        Real tip_reference = 0;
        for (int ix = 0; ix < nx; ++ix)
        {
            Real x = (ix + 0.5) * dp, mean = sum[ix] / (cells * cells);
            Real eb_static = load_per_length * x * x * (6 * L * L - 4 * L * x + x * x) / (24 * young * inertia);
            Real eb = load_mode == 1 ? eb_static : amplitude * mode(x / L) * std::sin(omega * time);
            Real shear = load_per_length * (L * x - x * x / 2) / (kappa * shear_modulus * area);
            profiles << time / period << ',' << x / L << ',' << mean / H << ','
                     << std::sqrt(std::max(Real(0), sum2[ix] / (cells * cells) - mean * mean)) / H << ','
                     << eb / H << ',' << (eb_static + shear) / H << '\n';
            if (ix == nx - 1) tip_reference = eb;
        }
        history << time << ',' << time / period << ',' << modal / mode_norm / H << ',' << sum.back() / (cells * cells) / H
                << ',' << tip_reference / H << ',' << kinetic << ',' << elastic << ',' << drag_loss << ',' << clamp_error / H
                << ',' << max_F << ',' << (hg_norm > 0 ? hg_net.norm() / hg_norm : 0)
                << ',' << (torque_norm > 0 ? torque.norm() / torque_norm : 0) << ',' << reaction
                << ',' << (int_norm > 0 ? int_net.norm() / int_norm : 0) << ',' << alive_count << ',' << damaged_count << ',' << failures << ',' << steps
                << ',' << std::sqrt(acceleration_sq / free_mass) / (load_mode == 1 ? gravity : amplitude * omega * omega)
                << ',' << std::sqrt(2 * kinetic / free_mass) / (amplitude * omega) << '\n';
        history.flush(); profiles.flush();
    };
    record();
    const int samples = int(std::ceil(cycles * samples_per_cycle));
    std::cout << "particles=" << count << " directed_bonds=" << bonds << " T=" << period << std::endl;
    for (int sample = 1; sample <= samples; ++sample)
    {
        Real target = period * cycles * sample / samples;
        while (time < target)
        {
            Real remaining = target - time;
            Real dt = peridynamics::boundedTimeStep(time_step.exec(), {remaining});
            max_dt = std::max(max_dt, dt); min_dt = std::min(min_dt, dt);
            damp_velocities(dt / 2);
            initialize.exec(dt);
            for (size_t i = 0; i < count; ++i)
                if (reference[i][0] > 0) prior[i][1] = gravity;
            first.exec(dt); second.exec(dt);
            if (stabilization == 1) legacy.exec(dt); else modern.exec(dt);
            if (numerical_damping) damping.exec(dt);
            third.exec(dt);
            if (dt == remaining) force_positions.assign(pos, pos + count);
            reaction = 0;
            // Enforce the constraint BEFORE the second half drift and work tally.
            for (size_t i = 0; i < count; ++i)
                if (reference[i][0] <= 0) { reaction -= mass[i] * acc[i][1]; acc[i].setZero(); }
            fourth.exec(dt); damp_velocities(dt / 2);
            time = dt == remaining ? target : time + dt;
            ++steps;
            if (steps > 10000000) throw std::runtime_error("step safety limit exceeded");
        }
        record();
        if (sample % 40 == 0)
            std::cout << "t/T=" << time / period << " steps=" << steps << " wall_s="
                      << std::chrono::duration<Real>(std::chrono::steady_clock::now() - start).count() << std::endl;
    }
    std::ofstream status("completion.json");
    status << std::setprecision(17) << "{\"status\":\"completed\",\"steps\":" << steps
           << ",\"min_dt\":" << min_dt << ",\"max_dt\":" << max_dt << ",\"time\":" << time
           << ",\"wall_seconds\":" << std::chrono::duration<Real>(std::chrono::steady_clock::now() - start).count() << "}\n";
    return 0;
}
}
int main(int argc, char **argv)
{
    try { return run(argc, argv); }
    catch (const std::exception &error)
    {
        std::cerr << "CANTILEVER_BENCHMARK_FAILED: " << error.what() << std::endl;
        std::ofstream("failure.txt") << error.what() << '\n';
        return 1;
    }
}
