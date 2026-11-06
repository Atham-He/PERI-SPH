/** All-solid SPH capsule demonstrator. SI units throughout; NOT experimentally calibrated.
 * Native SPH volume core and native SPH midsurface shell retain their own inertia.
 * A case-local affine-reproducing penalty tie transfers equal/opposite interface forces.
 */
#include "sphinxsys.h"
#include <Eigen/Cholesky>
#include <Eigen/Geometry>
#include <chrono>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <memory>
#include <limits>
#include <numeric>
#include <stdexcept>

using namespace SPH;
namespace fs = std::filesystem;

struct Options
{
    bool coated = true, patch_test = false, rest_test = false, smoke_test = false;
    Real diameter = 1e-3, thickness = 1e-5, spacing = 1e-4;
    Real core_E = 1e6, shell_E = 1e7, core_nu = 0.3, shell_nu = 0.3;
    Real core_rho = 1355, shell_rho = 1000;
    Real strain = 0.01, load_time = 3e-4, hold_time = 3e-4;
    Real contact_factor = 100, tie_factor = 20, cfl = 0.2, dt_scale = 1;
    Real drag = 0, yield_stress = 2e4, hardening = 0;
    std::string material = "elastic", output = "run";
    size_t max_steps = 2000000;

    void parse(int ac, char **av)
    {
        for (int i = 1; i < ac; ++i)
        {
            std::string a(av[i]);
            if (a == "--patch-test") { patch_test = true; continue; }
            if (a == "--rest-test") { rest_test = true; strain = 0; load_time = 2e-5; hold_time = 0; continue; }
            if (a == "--smoke-test") { smoke_test = true; strain = 0.002; load_time = 1e-4; hold_time = 2e-5; continue; }
            size_t eq = a.find('=');
            if (a.substr(0, 2) != "--" || eq == std::string::npos)
                throw std::invalid_argument("use --name=value or --patch-test");
            std::string key = a.substr(2, eq - 2), val = a.substr(eq + 1);
            if (key == "case")
            {
                if (val != "bare" && val != "coated") throw std::invalid_argument("case must be bare or coated");
                coated = val == "coated";
            }
            else if (key == "material") material = val;
            else if (key == "output") output = val;
            else
            {
                size_t used = 0;
                Real x = std::stod(val, &used);
                if (used != val.size() || !std::isfinite(x)) throw std::invalid_argument("non-finite/invalid number: " + key);
                if (key == "diameter") diameter = x;
                else if (key == "thickness") thickness = x;
                else if (key == "spacing") spacing = x;
                else if (key == "core-E") core_E = x;
                else if (key == "shell-E") shell_E = x;
                else if (key == "core-nu") core_nu = x;
                else if (key == "shell-nu") shell_nu = x;
                else if (key == "core-density") core_rho = x;
                else if (key == "shell-density") shell_rho = x;
                else if (key == "strain") strain = x;
                else if (key == "load-time") load_time = x;
                else if (key == "hold-time") hold_time = x;
                else if (key == "contact-factor") contact_factor = x;
                else if (key == "tie-factor") tie_factor = x;
                else if (key == "cfl") cfl = x;
                else if (key == "dt-scale") dt_scale = x;
                else if (key == "drag") drag = x;
                else if (key == "yield-stress") yield_stress = x;
                else if (key == "hardening") hardening = x;
                else if (key == "max-steps")
                {
                    if (x < 1 || x > 1e8 || std::floor(x) != x) throw std::invalid_argument("invalid max-steps");
                    max_steps = static_cast<size_t>(x);
                }
                else throw std::invalid_argument("unknown option: " + key);
            }
        }
        if (diameter <= 0 || thickness <= 0 || thickness >= diameter / 4 || spacing < diameter / 80 || spacing > diameter / 4 ||
            core_E <= 0 || shell_E <= 0 || core_rho <= 0 || shell_rho <= 0 ||
            core_nu < 0 || core_nu >= 0.49 || shell_nu < 0 || shell_nu >= 0.49 ||
            strain < 0 || strain > 0.05 || load_time <= 0 || hold_time < 0 ||
            contact_factor <= 0 || tie_factor <= 0 || cfl <= 0 || cfl > 0.3 ||
            dt_scale <= 0 || dt_scale > 1 || drag < 0 || yield_stress <= 0 || hardening < 0)
            throw std::invalid_argument("parameters outside supported positive/small-strain ranges");
        if (material != "elastic" && material != "j2") throw std::invalid_argument("material must be elastic or j2");
    }
};

namespace SPH
{
class CapsuleCoreGenerator;
template <> class ParticleGenerator<BaseParticles, CapsuleCoreGenerator> : public ParticleGenerator<BaseParticles>
{
    Real radius_, dp_;
  public:
    ParticleGenerator(SPHBody &b, BaseParticles &p, Real r, Real dp)
        : ParticleGenerator<BaseParticles>(b, p), radius_(r), dp_(dp) {}
    void prepareGeometricData() override
    {
        int n = static_cast<int>(std::ceil(2 * radius_ / dp_));
        StdVec<Vecd> points;
        for (int i = 0; i < n; ++i)
            for (int j = 0; j < n; ++j)
                for (int k = 0; k < n; ++k)
                {
                    Vecd x((i - 0.5 * (n - 1)) * dp_, (j - 0.5 * (n - 1)) * dp_, (k - 0.5 * (n - 1)) * dp_);
                    if (x.norm() < radius_) points.push_back(x);
                }
        if (points.size() < 32) throw std::runtime_error("core resolution too low");
        // Equal reference volumes integrate the specified continuum sphere mass exactly.
        Real volume = 4 * Pi * std::pow(radius_, 3) / (3 * points.size());
        for (const auto &x : points) addPositionAndVolumetricMeasure(x, volume);
    }
};
class CapsuleShellGenerator;
template <> class ParticleGenerator<SurfaceParticles, CapsuleShellGenerator> : public ParticleGenerator<SurfaceParticles>
{
    Real radius_, dp_, thickness_;
  public:
    ParticleGenerator(SPHBody &b, SurfaceParticles &p, Real r, Real dp, Real t)
        : ParticleGenerator<SurfaceParticles>(b, p), radius_(r), dp_(dp), thickness_(t) {}
    void prepareGeometricData() override
    {
        int half = std::max(32, static_cast<int>(std::ceil(2 * Pi * radius_ * radius_ / (dp_ * dp_))));
        Real area = 4 * Pi * radius_ * radius_ / (2 * half);
        // Antipodal Fibonacci points: exact north/south pole contact, zero centroid.
        for (int i = 0; i < half; ++i)
        {
            Real z = 1.0 - Real(i) / Real(half), phi = i * Pi * (3.0 - std::sqrt(5.0));
            Real r = std::sqrt(std::max(Real(0), 1 - z * z));
            Vecd n(r * std::cos(phi), r * std::sin(phi), z);
            for (Real sign : {Real(1), Real(-1)})
            {
                addPositionAndVolumetricMeasure(sign * radius_ * n, area);
                addSurfaceProperties(sign * n, thickness_);
            }
        }
    }
};
} // namespace SPH

struct State
{
    size_t n;
    Vecd *x, *v, *prior, *force;
    Real *mass, *measure;
    Matd *F;
    StdVec<Vecd> x0;
    explicit State(SolidBody &b)
    {
        auto &p = b.getBaseParticles();
        n = p.TotalRealParticles();
        x = p.getVariableDataByName<Vecd>("Position");
        v = p.getVariableDataByName<Vecd>("Velocity");
        prior = p.registerStateVariableData<Vecd>("ForcePrior");
        force = p.registerStateVariableData<Vecd>("Force");
        mass = p.getVariableDataByName<Real>("Mass");
        measure = p.getVariableDataByName<Real>("VolumetricMeasure");
        F = p.getVariableDataByName<Matd>("DeformationGradient");
        x0.assign(x, x + n);
    }
    void clearPrior() { for (size_t i = 0; i < n; ++i) prior[i].setZero(); }
    Real totalMass() const { return std::accumulate(mass, mass + n, Real(0)); }
};

/** Fixed-reference MLS tie. Energy = sum_s k_s |x_s - sum_c W_sc x_c|^2/2.
 * Affine reproduction makes rigid rotations stress-free. The transpose scatter
 * conserves resultant force and instantaneous power; shell mass is not discarded.
 */
struct BondedInterface
{
    struct Row { std::vector<size_t> j; std::vector<Real> w; Real k; };
    State &core, &shell;
    std::vector<Row> rows;
    Real stable_dt = std::numeric_limits<Real>::infinity(), max_gap = 0, energy = 0, balance = 0, power_residual = 0;
    BondedInterface(State &c, State &s, Real dp, Real stiffness) : core(c), shell(s)
    {
        std::vector<Real> rate_c(c.n, 0);
        Real max_rate = 0;
        for (size_t i = 0; i < s.n; ++i)
        {
            Row row;
            Eigen::Matrix4d moment = Eigen::Matrix4d::Zero();
            std::vector<Eigen::Vector4d> bases;
            std::vector<Real> kernel;
            for (size_t j = 0; j < c.n; ++j)
            {
                Vecd q = (c.x0[j] - s.x0[i]) / dp;
                if (q.norm() >= 2.5) continue;
                Eigen::Vector4d b; b << 1, q[0], q[1], q[2];
                Real w = std::exp(-q.squaredNorm());
                moment += w * b * b.transpose();
                row.j.push_back(j); bases.push_back(b); kernel.push_back(w);
            }
            Eigen::LDLT<Eigen::Matrix4d> solver(moment);
            if (row.j.size() < 8 || solver.info() != Eigen::Success || solver.vectorD().minCoeff() <= 1e-10)
                throw std::runtime_error("singular core-shell mapping; refine particles or increase support");
            Eigen::Vector4d e0(1, 0, 0, 0), coeff = solver.solve(e0);
            Real sum_abs = 1, sum_w = 0;
            Vecd mapped = Vecd::Zero();
            for (size_t k = 0; k < row.j.size(); ++k)
            {
                Real w = kernel[k] * bases[k].dot(coeff);
                row.w.push_back(w); sum_w += w; sum_abs += std::abs(w);
                mapped += w * c.x0[row.j[k]];
            }
            if (std::abs(sum_w - 1) > 1e-9 || (mapped - s.x0[i]).norm() > 1e-9 * dp)
                throw std::runtime_error("affine reproduction failed");
            row.k = stiffness * s.measure[i];
            max_rate = std::max(max_rate, row.k * sum_abs / s.mass[i]);
            for (size_t k = 0; k < row.j.size(); ++k)
                rate_c[row.j[k]] += row.k * std::abs(row.w[k]) * sum_abs / c.mass[row.j[k]];
            rows.push_back(std::move(row));
        }
        for (Real r : rate_c) max_rate = std::max(max_rate, r);
        stable_dt = 0.2 / std::sqrt(max_rate);
    }
    void apply()
    {
        max_gap = energy = balance = power_residual = 0;
        Vecd sum = Vecd::Zero();
        Real force_scale = 0, power = 0;
        for (size_t i = 0; i < shell.n; ++i)
        {
            const auto &row = rows[i];
            Vecd mapped = Vecd::Zero(), mapped_v = Vecd::Zero();
            for (size_t k = 0; k < row.j.size(); ++k)
            { mapped += row.w[k] * core.x[row.j[k]]; mapped_v += row.w[k] * core.v[row.j[k]]; }
            Vecd gap = shell.x[i] - mapped, f = -row.k * gap;
            shell.prior[i] += f;
            sum += f; power += f.dot(shell.v[i]);
            for (size_t k = 0; k < row.j.size(); ++k)
            {
                Vecd reaction = -row.w[k] * f;
                core.prior[row.j[k]] += reaction;
                sum += reaction; power += reaction.dot(core.v[row.j[k]]);
            }
            power_residual += -f.dot(shell.v[i] - mapped_v); // dU/dt
            max_gap = std::max(max_gap, gap.norm());
            energy += 0.5 * row.k * gap.squaredNorm(); force_scale += f.norm();
        }
        balance = sum.norm() / std::max(force_scale, Real(1e-30));
        power_residual += power; // external interface power + dU/dt = 0
    }
    void patchTest(Real dp)
    {
        Matd R = Eigen::AngleAxisd(0.37, Vecd(1, 2, 3).normalized()).toRotationMatrix();
        Vecd shift(0.2 * dp, -0.1 * dp, 0.3 * dp);
        for (size_t i = 0; i < core.n; ++i) core.x[i] = R * core.x0[i] + shift;
        for (size_t i = 0; i < shell.n; ++i) shell.x[i] = R * shell.x0[i] + shift;
        core.clearPrior(); shell.clearPrior(); apply();
        if (max_gap / dp > 1e-10) throw std::runtime_error("rigid-motion patch failed");
        for (size_t i = 0; i < core.n; ++i)
        { core.v[i] = Vecd(std::sin(Real(i)), std::cos(Real(i)), 0.2); core.x[i] += 1e-3 * dp * core.v[i]; }
        for (size_t i = 0; i < shell.n; ++i)
        { shell.v[i] = Vecd(0.1, std::sin(Real(i)), -0.3); shell.x[i] += 1e-3 * dp * shell.v[i]; }
        core.clearPrior(); shell.clearPrior(); apply();
        Vecd torque = Vecd::Zero(); Real scale = 0;
        for (size_t i = 0; i < core.n; ++i) { torque += core.x[i].cross(core.prior[i]); scale += core.prior[i].norm(); }
        for (size_t i = 0; i < shell.n; ++i) { torque += shell.x[i].cross(shell.prior[i]); scale += shell.prior[i].norm(); }
        Real torque_residual = torque.norm() / std::max(dp * scale, Real(1e-30));
        if (balance > 1e-10 || torque_residual > 1e-10 || std::abs(power_residual) > 1e-10 * scale)
            throw std::runtime_error("force/torque/power patch failed");
        std::cout << "PATCH PASS force=" << balance << " torque=" << torque_residual << " power=" << power_residual << '\n';
        for (size_t i = 0; i < core.n; ++i) { core.x[i] = core.x0[i]; core.v[i].setZero(); }
        for (size_t i = 0; i < shell.n; ++i) { shell.x[i] = shell.x0[i]; shell.v[i].setZero(); }
        core.clearPrior(); shell.clearPrior();
    }
};

struct PlaneContact
{
    State &state;
    Real offset, stiffness, stable_dt = std::numeric_limits<Real>::infinity(), upper = 0, lower = 0, penetration = 0, energy = 0;
    PlaneContact(State &s, Real off, Real stiffness_per_area, Real depth)
        : state(s), offset(off), stiffness(stiffness_per_area)
    {
        for (size_t i = 0; i < s.n; ++i)
        {
            Real k = stiffness * s.measure[i] / depth;
            stable_dt = std::min(stable_dt, Real(0.2) * std::sqrt(s.mass[i] / k));
        }
    }
    void apply(Real bottom, Real top, Real depth)
    {
        upper = lower = penetration = energy = 0;
        for (size_t i = 0; i < state.n; ++i)
        {
            Real k = stiffness * state.measure[i] / depth;
            Real pbot = std::max(bottom - (state.x[i][2] - offset), Real(0));
            Real ptop = std::max(state.x[i][2] + offset - top, Real(0));
            state.prior[i][2] += k * (pbot - ptop);
            lower += k * pbot; upper += k * ptop;
            penetration = std::max(penetration, std::max(pbot, ptop));
            energy += 0.5 * k * (pbot * pbot + ptop * ptop);
        }
    }
};

int run(const Options &o)
{
    fs::path output_root = fs::absolute(o.output);
    if (o.patch_test || o.rest_test || o.smoke_test)
        output_root = fs::temp_directory_path() / ("sph-capsule-patch-" + std::to_string(std::chrono::steady_clock::now().time_since_epoch().count()));
    if (!o.patch_test && !o.rest_test && !o.smoke_test)
    {
        if (fs::exists(output_root) && !fs::is_empty(output_root)) throw std::runtime_error("output directory must be new or empty");
    }
    // Native IO creates/rotates ./output. Confine even those housekeeping files
    // to this new run, never the caller's existing working directory.
    fs::create_directories(output_root / "runtime");
    struct WorkingDirectory
    {
        fs::path previous = fs::current_path();
        explicit WorkingDirectory(const fs::path &p) { fs::current_path(p); }
        ~WorkingDirectory() { fs::current_path(previous); }
    } cwd(output_root / "runtime");
    Real R = o.diameter / 2, rc = o.coated ? R - o.thickness : R, rs = R - o.thickness / 2;
    SPHSystem system(BoundingBoxd(Vecd::Constant(-2 * R), Vecd::Constant(2 * R)), o.spacing, 1);
    if (!o.patch_test) IO::getEnvironment().resetOutputFolder((output_root / "vtp").string(), true);
    SolidBody core(system, makeShared<GeometricShapeBall>(Vecd::Zero(), rc, "Core"));
    core.defineAdaptationRatios(1.15, 1);
    if (o.material == "j2") core.defineMatterMaterial<HardeningPlasticSolid>(o.core_rho, o.core_E, o.core_nu, o.yield_stress, o.hardening);
    else core.defineMatterMaterial<SaintVenantKirchhoffSolid>(o.core_rho, o.core_E, o.core_nu);
    core.generateParticles<BaseParticles, CapsuleCoreGenerator>(rc, o.spacing);
    InnerRelation core_inner(core);
    InteractionWithUpdate<LinearGradientCorrectionMatrixInner> correction(core_inner);
    Dynamics1Level<solid_dynamics::Integration1stHalfPK2RightCauchy> elastic_first(core_inner);
    std::unique_ptr<Dynamics1Level<solid_dynamics::DecomposedPlasticIntegration1stHalf>> plastic_first;
    if (o.material == "j2") plastic_first = std::make_unique<Dynamics1Level<solid_dynamics::DecomposedPlasticIntegration1stHalf>>(core_inner);
    Dynamics1Level<solid_dynamics::Integration2ndHalf> core_second(core_inner);
    ReduceDynamics<solid_dynamics::AcousticTimeStep> core_dt(core, o.cfl);
    State c(core);

    std::unique_ptr<SolidBody> shell;
    std::unique_ptr<InnerRelation> shell_inner;
    std::unique_ptr<InteractionDynamics<thin_structure_dynamics::ShellCorrectConfiguration>> shell_correction;
    std::unique_ptr<Dynamics1Level<thin_structure_dynamics::ShellStressRelaxationFirstHalf>> shell_first;
    std::unique_ptr<Dynamics1Level<thin_structure_dynamics::ShellStressRelaxationSecondHalf>> shell_second;
    std::unique_ptr<ReduceDynamics<thin_structure_dynamics::ShellAcousticTimeStepSize>> shell_dt;
    std::unique_ptr<State> s;
    if (o.coated || o.patch_test)
    {
        shell = std::make_unique<SolidBody>(system, makeShared<DefaultShape>("ECShell"));
        shell->defineAdaptationRatios(1.15, 1);
        shell->defineMatterMaterial<SaintVenantKirchhoffSolid>(o.shell_rho, o.shell_E, o.shell_nu);
        shell->generateParticles<SurfaceParticles, CapsuleShellGenerator>(rs, o.spacing, o.thickness);
        shell_inner = std::make_unique<InnerRelation>(*shell);
        shell_correction = std::make_unique<InteractionDynamics<thin_structure_dynamics::ShellCorrectConfiguration>>(*shell_inner);
        shell_first = std::make_unique<Dynamics1Level<thin_structure_dynamics::ShellStressRelaxationFirstHalf>>(*shell_inner, 3, true);
        shell_second = std::make_unique<Dynamics1Level<thin_structure_dynamics::ShellStressRelaxationSecondHalf>>(*shell_inner);
        shell_dt = std::make_unique<ReduceDynamics<thin_structure_dynamics::ShellAcousticTimeStepSize>>(*shell, o.cfl);
        s = std::make_unique<State>(*shell);
    }
    system.initializeSystemCellLinkedLists();
    system.initializeSystemConfigurations();
    correction.exec();
    if (shell_correction) shell_correction->exec();
    std::unique_ptr<BondedInterface> tie;
    if (s) tie = std::make_unique<BondedInterface>(c, *s, o.spacing, o.tie_factor * o.core_E / o.spacing);
    if (o.patch_test) { tie->patchTest(o.spacing); return 0; }
    State &surface = s ? *s : c;
    Real contact_depth = s ? 1 : o.spacing;
    Real surface_offset = o.thickness / 2;
    if (!s)
    {
        Real top_node = 0;
        for (const auto &x : c.x0) top_node = std::max(top_node, x[2]);
        surface_offset = R - top_node;
    }
    PlaneContact contact(surface, surface_offset,
                         o.contact_factor * o.core_E / o.spacing, contact_depth);
    BodyStatesRecordingToVtp vtp(system);
    core.getBaseParticles().addVariableToWrite<Vecd>("Velocity");
    core.getBaseParticles().addVariableToWrite<Vecd>("ForcePrior");
    core.getBaseParticles().addVariableToWrite<Matd>("DeformationGradient");
    if (shell) shell->getBaseParticles().addVariableToWrite<Vecd>("Velocity");
    Real *plastic_alpha = nullptr;
    if (o.material == "j2")
    {
        plastic_alpha = core.getBaseParticles().getVariableDataByName<Real>("HardeningParameter");
        core.getBaseParticles().addVariableToWrite<Real>("HardeningParameter");
    }
    std::ofstream cfg(output_root / "parameters.txt");
    if (!cfg) throw std::runtime_error("cannot write parameters");
    cfg << std::setprecision(17) << "status=DEMONSTRATION_NOT_CALIBRATED\nunits=SI\ncase=" << (s ? "coated" : "bare")
        << "\nsource_sha256=" << CAPSULE_SOURCE_SHA256
        << "\ncore_material=" << o.material << "\nshell_material=" << (s ? "elastic_plane_stress" : "none")
        << "\ninterface=" << (s ? "affine_penalty_tie" : "none") << '\n'
        << "diameter=" << o.diameter << "\nthickness=" << o.thickness << "\nspacing=" << o.spacing
        << "\ncore_E=" << o.core_E << "\ncore_nu=" << o.core_nu << "\ncore_density=" << o.core_rho
        << "\nshell_E=" << o.shell_E << "\nshell_nu=" << o.shell_nu << "\nshell_density=" << o.shell_rho
        << "\nyield_stress=" << o.yield_stress << "\nhardening=" << o.hardening
        << "\nload_time=" << o.load_time << "\nhold_time=" << o.hold_time << "\nstrain=" << o.strain
        << "\nmaximum_plate_speed=" << 0.5 * Pi * o.strain * o.diameter / o.load_time
        << "\npaper_plate_speed=1e-6\ncontact_factor=" << o.contact_factor << "\ntie_factor=" << o.tie_factor
        << "\ncfl=" << o.cfl << "\ndt_scale=" << o.dt_scale << "\ndrag=" << o.drag
        << "\ncore_particles=" << c.n << "\nshell_particles=" << (s ? s->n : 0)
        << "\ncore_mass=" << c.totalMass() << "\nshell_mass=" << (s ? s->totalMass() : 0)
        << "\ncontact_surface_offset=" << surface_offset
        << "\ncontact_dt_limit=" << contact.stable_dt << "\ntie_dt_limit=" << (tie ? tie->stable_dt : 0) << '\n';
    std::ofstream csv(output_root / "history.csv");
    if (!csv) throw std::runtime_error("cannot write history");
    csv << "step,time_s,dt_s,compression_m,upper_N,lower_N,kinetic_J,core_elastic_J,tie_J,contact_J,max_tie_gap_m,max_penetration_m,min_J,interface_force_residual,interface_power_residual_W,external_work_J,drag_loss_J,max_plastic_alpha\n";
    csv << std::setprecision(17);
    Real &time = *system.getSystemVariableDataByName<Real>("PhysicalTime");
    Real end = o.load_time + o.hold_time, next_csv = 0, next_vtp = 0, work = 0, drag_loss = 0;
    Real worst_J = 1, largest_gap = 0, largest_penetration = 0, max_force = 0;
    size_t step = 0, frame = 0;
    auto compression = [&](Real t) { return 0.5 * o.strain * o.diameter * (1 - std::cos(Pi * std::min(t / o.load_time, Real(1)))); };
    auto inspect = [&](State &a, Real &kinetic, Real &U, Real &minJ, bool volume_core)
    {
        for (size_t i = 0; i < a.n; ++i)
        {
            Real J = a.F[i].determinant();
            if (!a.x[i].allFinite() || !a.v[i].allFinite() || !std::isfinite(J) || J <= 0.2)
                throw std::runtime_error("non-finite state or invalid deformation determinant");
            minJ = std::min(minJ, J);
            kinetic += 0.5 * a.mass[i] * a.v[i].squaredNorm();
            if (volume_core && o.material == "elastic")
            {
                Matd e = 0.5 * (a.F[i].transpose() * a.F[i] - Matd::Identity());
                Real mu = o.core_E / (2 * (1 + o.core_nu));
                Real lambda = o.core_E * o.core_nu / ((1 + o.core_nu) * (1 - 2 * o.core_nu));
                U += a.measure[i] * (0.5 * lambda * e.trace() * e.trace() + mu * e.squaredNorm());
            }
        }
    };
    while (true)
    {
        c.clearPrior(); if (s) s->clearPrior();
        if (tie) tie->apply();
        Real delta = compression(time);
        contact.apply(-R, R - delta, contact_depth);
        Real kinetic = 0, core_U = 0, minJ = 1;
        inspect(c, kinetic, core_U, minJ, true); if (s) inspect(*s, kinetic, core_U, minJ, false);
        largest_gap = std::max(largest_gap, tie ? tie->max_gap : 0);
        largest_penetration = std::max(largest_penetration, contact.penetration);
        worst_J = std::min(worst_J, minJ); max_force = std::max(max_force, contact.upper);
        if (time < end && end - time < 1e-12 * end) { time = end; continue; }
        Real dt = std::min(core_dt.exec(), contact.stable_dt);
        if (s) dt = std::min({dt, shell_dt->exec(), tie->stable_dt});
        if (o.drag > 0) dt = std::min(dt, Real(0.2) / o.drag);
        dt = std::min(o.dt_scale * dt, end - time);
        if (time >= next_csv || time >= end)
        {
            Real max_alpha = plastic_alpha ? *std::max_element(plastic_alpha, plastic_alpha + c.n) : 0;
            csv << step << ',' << time << ',' << dt << ',' << delta << ',' << contact.upper << ',' << contact.lower << ','
                << kinetic << ',' << core_U << ',' << (tie ? tie->energy : 0) << ',' << contact.energy << ','
                << (tie ? tie->max_gap : 0) << ',' << contact.penetration << ',' << minJ << ','
                << (tie ? tie->balance : 0) << ',' << (tie ? tie->power_residual : 0) << ',' << work << ',' << drag_loss << ',' << max_alpha << '\n';
            next_csv = time + end / 300;
        }
        if (time >= next_vtp || time >= end)
        {
            vtp.writeToFile(frame++); next_vtp = time + end / 10;
            std::cout << "step=" << step << " time=" << time << " upper=" << contact.upper << " lower=" << contact.lower
                      << " gap=" << (tie ? tie->max_gap : 0) << " minJ=" << minJ << std::endl;
        }
        if (time >= end) break;
        if (step >= o.max_steps || !std::isfinite(dt) || dt < 1e-15) throw std::runtime_error("time step/step-count guard");
        // Split the native first-half operators at the SAME midpoint for both
        // bodies. Lagging penalty forces by a half drift injects interface energy.
        auto initialize = [&](auto &method, State &a)
        { for (size_t i = 0; i < a.n; ++i) method.initialization(i, dt); };
        auto interact = [&](auto &method, State &a)
        { for (size_t i = 0; i < a.n; ++i) method.interaction(i, dt); };
        auto update = [&](auto &method, State &a)
        { for (size_t i = 0; i < a.n; ++i) method.update(i, dt); };
        if (plastic_first) initialize(*plastic_first, c); else initialize(elastic_first, c);
        if (shell_first) initialize(*shell_first, *s);
        c.clearPrior(); if (s) s->clearPrior();
        if (tie) tie->apply();
        contact.apply(-R, R - compression(time + 0.5 * dt), contact_depth);
        if (o.drag > 0)
        {
            auto damp = [&](State &a) { for (size_t i = 0; i < a.n; ++i) a.prior[i] -= o.drag * a.mass[i] * a.v[i]; };
            damp(c); if (s) damp(*s);
            drag_loss += 2 * o.drag * kinetic * dt;
        }
        work += contact.upper * (compression(time + dt) - delta);
        if (plastic_first) interact(*plastic_first, c); else interact(elastic_first, c);
        if (shell_first) interact(*shell_first, *s);
        if (plastic_first) update(*plastic_first, c); else update(elastic_first, c);
        if (shell_first) update(*shell_first, *s);
        core_second.exec(dt); if (shell_second) shell_second->exec(dt);
        time = std::min(end, time + dt); ++step;
    }
    if (o.rest_test && (max_force > 1e-8 * o.core_E * o.diameter * o.diameter ||
                        worst_J < 1 - 1e-6 || largest_gap > 1e-6 * o.spacing))
        throw std::runtime_error("zero-load equilibrium regression failed");
    if (o.smoke_test && (max_force <= 1e-10 || worst_J < 0.9 || largest_gap > 0.02 * o.spacing))
        throw std::runtime_error("loaded solid-shell smoke regression failed");
    std::ofstream summary(output_root / "summary.json");
    if (!summary) throw std::runtime_error("cannot write summary");
    summary << std::setprecision(17) << "{\n  \"status\": \"COMPLETED_DEMONSTRATION\",\n  \"steps\": " << step
            << ",\n  \"source_sha256\": \"" << CAPSULE_SOURCE_SHA256 << "\""
            << ",\n  \"final_time_s\": " << time << ",\n  \"min_J\": " << worst_J << ",\n  \"max_upper_N\": " << max_force
            << ",\n  \"max_tie_gap_m\": " << largest_gap << ",\n  \"max_penetration_m\": " << largest_penetration
            << ",\n  \"experimental_calibration\": false,\n  \"shell_energy_and_rotational_energy_audited\": false\n}\n";
    return 0;
}

int main(int ac, char **av)
{
    try { Options o; o.parse(ac, av); return run(o); }
    catch (const std::exception &e) { std::cerr << "ERROR: " << e.what() << std::endl; return 1; }
}
