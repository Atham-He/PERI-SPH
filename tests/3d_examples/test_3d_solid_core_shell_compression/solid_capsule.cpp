/** All-solid SPH capsule: fixed-input quasi-static cases and accelerated demos.
 * SI units throughout; NOT experimentally calibrated.
 * Native SPH volume core and native SPH midsurface shell retain their own inertia.
 * A case-local affine-reproducing penalty tie transfers equal/opposite interface forces.
 */
#include "sphinxsys.h"
#include "drucker_prager_cap_solid.h"
#include "solid_residual_integration.h"
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
    bool quasistatic = false;
    Real diameter = 1e-3, thickness = 1e-5, spacing = 1e-4;
    Real core_E = 1e6, shell_E = 1e7, core_nu = 0.3, shell_nu = 0.3;
    Real core_rho = 1355, shell_rho = 1000;
    Real strain = 0.01, load_time = 3e-4, hold_time = 3e-4;
    Real contact_factor = 100, tie_factor = 20, cfl = 0.2, dt_scale = 1;
    Real surface_mls_scale = 1, contact_quadrature_scale = 1, max_dt = std::numeric_limits<Real>::max();
    Real drag = 0, yield_stress = 2e4, hardening = 0;
    Real relaxation_damping = 0.002, force_tolerance = 0.05, balance_tolerance = 0.02, kinetic_tolerance = 1e-4;
    size_t load_steps = 20, ramp_steps = 1000, max_stage_steps = 200000, threads = 1;
    std::string solid_stabilization = "none";
    std::string material = "elastic", output = "run", input_file, plate_contact = "legacy", interface_surface = "legacy_mid";
    size_t max_steps = 2000000;
    std::string dpc_table;
    Real dpc_r0 = 0, dpc_solid_density = 0, constraint_E = 0, dpc_shear_dilation = 1, dpc_shear_hardening = 0;

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
            if (key == "config")
            {
                std::ifstream input(val);
                if (!input) throw std::invalid_argument("cannot read config: " + val);
                std::vector<std::string> words{"config"};
                for (std::string line; std::getline(input, line); )
                {
                    line = line.substr(0, line.find('#'));
                    size_t first = line.find_first_not_of(" \t\r"), last = line.find_last_not_of(" \t\r");
                    if (first == std::string::npos) continue;
                    line = line.substr(first, last - first + 1);
                    if (line.substr(0, 7) == "config=") throw std::invalid_argument("nested configs unsupported");
                    words.push_back("--" + line);
                }
                std::vector<char *> argv;
                for (auto &word : words) argv.push_back(word.data());
                parse(static_cast<int>(argv.size()), argv.data());
                input_file = fs::absolute(val).string();
            }
            else if (key == "case")
            {
                if (val != "bare" && val != "coated") throw std::invalid_argument("case must be bare or coated");
                coated = val == "coated";
            }
            else if (key == "solid-stabilization")
            {if(val!="none"&&val!="quadratic_residual")throw std::invalid_argument("solid-stabilization must be none or quadratic_residual");solid_stabilization=val;}
            else if (key == "material") material = val;
            else if (key == "dpc-table") dpc_table = fs::absolute(val).string();
            else if (key == "plate-contact") plate_contact = val;
            else if (key == "interface-surface") interface_surface = val;
            else if (key == "output") output = val;
            else
            {
                size_t used = 0;
                Real x = std::stod(val, &used);
                if (used != val.size() || !std::isfinite(x)) throw std::invalid_argument("non-finite/invalid number: " + key);
                if (key == "dpc-r0") dpc_r0 = x;
                else if (key == "dpc-solid-density") dpc_solid_density = x;
                else if (key == "dpc-shear-dilation") dpc_shear_dilation = x;
                else if (key == "dpc-shear-hardening") dpc_shear_hardening = x;
                else if (key == "constraint-E") constraint_E = x;
                else if (key == "diameter") diameter = x;
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
                else if (key == "surface-mls-scale") surface_mls_scale = x;
                else if (key == "contact-quadrature-scale") contact_quadrature_scale = x;
                else if (key == "max-dt") max_dt = x;
                else if (key == "drag") drag = x;
                else if (key == "yield-stress") yield_stress = x;
                else if (key == "hardening") hardening = x;
                else if (key == "quasistatic")
                { if (x != 0 && x != 1) throw std::invalid_argument("quasistatic must be 0 or 1"); quasistatic = x == 1; }
                else if (key == "relaxation-damping") relaxation_damping = x;
                else if (key == "force-tolerance") force_tolerance = x;
                else if (key == "balance-tolerance") balance_tolerance = x;
                else if (key == "kinetic-tolerance") kinetic_tolerance = x;
                else if (key == "load-steps" || key == "ramp-steps" || key == "max-stage-steps" || key == "threads")
                {
                    if (x < 1 || x > 10000000 || std::floor(x) != x) throw std::invalid_argument("invalid integer: " + key);
                    if (key == "load-steps") load_steps = static_cast<size_t>(x);
                    else if (key == "ramp-steps") ramp_steps = static_cast<size_t>(x);
                    else if (key == "max-stage-steps") max_stage_steps = static_cast<size_t>(x);
                    else threads = static_cast<size_t>(x);
                }
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
            strain < 0 || strain > 0.15 || load_time <= 0 || hold_time < 0 ||
            contact_factor <= 0 || tie_factor <= 0 || cfl <= 0 || cfl > 0.3 ||
            dt_scale <= 0 || dt_scale > 1 || drag < 0 || yield_stress <= 0 || hardening < 0)
            throw std::invalid_argument("parameters outside supported positive/pre-fracture ranges (strain <= 0.15)");
        if (material != "elastic" && material != "j2" && material != "dpc") throw std::invalid_argument("material must be elastic, j2 or dpc");
        if (constraint_E < 0) throw std::invalid_argument("constraint-E must be zero (use core E) or positive");
        if(dpc_shear_hardening<0 || (dpc_shear_hardening>0 && material!="dpc"))throw std::invalid_argument("dpc-shear-hardening requires nonnegative DPC input");
        if (material == "dpc")
        {
            auto model=dpc::Model::fromCSV(dpc_table,dpc_r0,dpc_shear_dilation,dpc_shear_hardening);auto p=model.parameters(dpc_r0);
            auto close=[](Real a,Real b){return std::abs(a-b)<=1e-9*std::max(Real(1),std::abs(b));};
            if(dpc_solid_density<=0 || !close(core_rho,dpc_solid_density*dpc_r0) || !close(core_E,p.E) || !close(core_nu,p.nu))
                throw std::invalid_argument("DPC initial E, nu and mass density must match table/r0/solid-density; no silent substitution");
        }
        if (plate_contact != "legacy" && plate_contact != "outer_surface")
            throw std::invalid_argument("plate-contact must be legacy or outer_surface");
        if (interface_surface != "legacy_mid" && interface_surface != "inner")
            throw std::invalid_argument("interface-surface must be legacy_mid or inner");
        if (surface_mls_scale < 1 || surface_mls_scale > 2 || max_dt <= 0 ||
            (surface_mls_scale != 1 && (coated || plate_contact != "outer_surface")))
            throw std::invalid_argument("surface-mls-scale in [1,2] is a bare outer_surface diagnostic; max-dt must be positive");
        if (contact_quadrature_scale < 0.25 || contact_quadrature_scale > 2 ||
            (contact_quadrature_scale != 1 && (coated || plate_contact != "outer_surface")))
            throw std::invalid_argument("contact-quadrature-scale in [0.25,2] requires bare outer_surface contact");
        if (relaxation_damping <= 0 || relaxation_damping > 0.1 || force_tolerance <= 0 || force_tolerance > .1 ||
            balance_tolerance <= 0 || balance_tolerance > .05 || kinetic_tolerance <= 0 || kinetic_tolerance > .01 || threads > 16 ||
            (quasistatic && (strain <= 0 || max_stage_steps <= ramp_steps || load_steps > 200 || drag != 0)))
            throw std::invalid_argument("invalid quasi-static solver settings");
    }
};

StdVec<Vecd> sphereDirections(Real radius, Real dp)
{
    int half = std::max(32, static_cast<int>(std::ceil(2 * Pi * radius * radius / (dp * dp))));
    StdVec<Vecd> directions;
    for (int i = 0; i < half; ++i)
    {
        Real z = 1.0 - Real(i) / Real(half), phi = i * Pi * (3.0 - std::sqrt(5.0));
        Real r = std::sqrt(std::max(Real(0), 1 - z * z));
        Vecd n(r * std::cos(phi), r * std::sin(phi), z);
        directions.push_back(n); directions.push_back(-n);
    }
    return directions;
}

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
    bool outer_sampling_;
  public:
    ParticleGenerator(SPHBody &b, SurfaceParticles &p, Real r, Real dp, Real t, bool outer_sampling = false)
        : ParticleGenerator<SurfaceParticles>(b, p), radius_(r), dp_(dp), thickness_(t), outer_sampling_(outer_sampling) {}
    void prepareGeometricData() override
    {
        auto directions = sphereDirections(radius_ + (outer_sampling_ ? thickness_ / 2 : 0), dp_);
        Real area = 4 * Pi * radius_ * radius_ / directions.size();
        for (const auto &n : directions)
        {
            addPositionAndVolumetricMeasure(radius_ * n, area);
            addSurfaceProperties(n, thickness_);
        }
    }
};
} // namespace SPH

struct State
{
    SolidBody &body;
    size_t n;
    Vecd *x, *v, *prior, *force;
    Real *mass, *measure;
    Matd *F;
    StdVec<Vecd> x0;
    explicit State(SolidBody &b) : body(b)
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

using DirectorJacobian = Eigen::Matrix<Real, 3, 2>;
DirectorJacobian directorJacobian(const Vecd &rotation)
{
    DirectorJacobian J;
    J.col(0) = thin_structure_dynamics::getVectorChangeRateAfterThinStructureRotation(Vecd(0, 0, 1), rotation, Vecd(1, 0, 0));
    J.col(1) = thin_structure_dynamics::getVectorChangeRateAfterThinStructureRotation(Vecd(0, 0, 1), rotation, Vecd(0, 1, 0));
    return J;
}

/** g=t*f/2 is conjugate to the GLOBAL UNIT director. External angular loading
 * is linear: Q=J^T*T*g, M=I*J^T*J, delta_alpha=M^-1*Q. This removes the
 * constrained normal component and preserves virtual work. Native internal
 * shell forces, constitutive stress update and gyroscopic terms are retained.
 */
class SurfaceLoadedShellFirst : public thin_structure_dynamics::ShellStressRelaxationFirstHalf
{
    Vecd *contact_director_force_, *interface_director_force_;
  public:
    explicit SurfaceLoadedShellFirst(BaseInnerRelation &inner)
        : ShellStressRelaxationFirstHalf(inner, 3, true),
          contact_director_force_(particles_->registerStateVariableData<Vecd>("ContactDirectorForce")),
          interface_director_force_(particles_->registerStateVariableData<Vecd>("InterfaceDirectorForce")) {}
    void interaction(size_t i, Real dt = 0)
    {
        ShellStressRelaxationFirstHalf::interaction(i, dt);
        Vecd director_force = contact_director_force_[i] + interface_director_force_[i];
        if (director_force.squaredNorm() == 0) return;
        DirectorJacobian J = directorJacobian(rotation_[i]);
        Eigen::Matrix2d metric = J.transpose() * J;
        if (metric.determinant() < 1e-10) throw std::runtime_error("singular shell rotation coordinates");
        Real inertia = mass_[i] * thickness_[i] * thickness_[i] / 12;
        Eigen::Vector2d Q = J.transpose() * (transformation_matrix0_[i] * director_force);
        Eigen::Vector2d extra_alpha = metric.ldlt().solve(Q) / inertia;
        dangular_vel_dt_[i][0] += extra_alpha[0]; dangular_vel_dt_[i][1] += extra_alpha[1];
        dpseudo_n_d2t_[i] += transformation_matrix0_[i].transpose() * (J * extra_alpha);
    }
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
    bool inner_surface;
    Real thickness, length_scale, reference_area_ratio, affine_error = 0;
    Real torque_residual = 0, power_relative = 0;
    Real worst_balance = 0, worst_torque_residual = 0, worst_power_relative = 0;
    Vecd *director, *initial_director, *director_force, *rotation, *omega;
    Matd *transform;
    StdVec<Vecd> reference_points, shell_points, mapped_points, forces;
    BondedInterface(State &c, State &s, Real dp, Real stiffness, bool inner = false, Real h = 0, Real area_ratio = 1)
        : core(c), shell(s), inner_surface(inner), thickness(h), length_scale(dp), reference_area_ratio(area_ratio)
    {
        auto &p = s.body.getBaseParticles();
        director = p.getVariableDataByName<Vecd>("PseudoNormal");
        initial_director = p.getVariableDataByName<Vecd>("InitialNormalDirection");
        director_force = p.getVariableDataByName<Vecd>("InterfaceDirectorForce");
        rotation = p.getVariableDataByName<Vecd>("Rotation");
        omega = p.getVariableDataByName<Vecd>("AngularVelocity");
        transform = p.getVariableDataByName<Matd>("TransformationMatrix");
        shell_points.resize(s.n); mapped_points.resize(s.n); forces.resize(s.n);
        std::vector<Real> rate_c(c.n, 0);
        Real max_rate = 0;
        for (size_t i = 0; i < s.n; ++i)
        {
            Vecd target = s.x0[i] - (inner_surface ? 0.5 * thickness : 0) * initial_director[i];
            reference_points.push_back(target);
            Row row;
            Eigen::Matrix4d moment = Eigen::Matrix4d::Zero();
            std::vector<Eigen::Vector4d> bases;
            std::vector<Real> kernel;
            for (size_t j = 0; j < c.n; ++j)
            {
                Vecd q = (c.x0[j] - target) / dp;
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
            affine_error = std::max({affine_error, std::abs(sum_w - 1), (mapped - target).norm() / dp});
            if (std::abs(sum_w - 1) > 1e-9 || (mapped - target).norm() > 1e-9 * dp)
                throw std::runtime_error("affine reproduction failed");
            row.k = stiffness * s.measure[i] * (inner_surface ? reference_area_ratio : 1);
            max_rate = std::max(max_rate, row.k * sum_abs / s.mass[i]);
            for (size_t k = 0; k < row.j.size(); ++k)
                rate_c[row.j[k]] += row.k * std::abs(row.w[k]) * sum_abs / c.mass[row.j[k]];
            rows.push_back(std::move(row));
        }
        for (Real r : rate_c) max_rate = std::max(max_rate, r);
        if (inner_surface)
        {
            // Row-block bound of sqrt(K)*B*M^-1*B^T*sqrt(K). Shell translation
            // plus unit-director inertia contributes at most 4*k/m per row;
            // the absolute shared-core terms bound cross-row coupling.
            std::vector<Real> shared_core(c.n, 0);
            for (const auto &row : rows)
                for (size_t a = 0; a < row.j.size(); ++a)
                    shared_core[row.j[a]] += std::sqrt(row.k) * std::abs(row.w[a]);
            max_rate = 0;
            for (size_t i = 0; i < s.n; ++i)
            {
                const auto &row = rows[i];
                Real rate = 4 * row.k / s.mass[i];
                for (size_t a = 0; a < row.j.size(); ++a)
                    rate += std::sqrt(row.k) * std::abs(row.w[a]) * shared_core[row.j[a]] / c.mass[row.j[a]];
                max_rate = std::max(max_rate, rate);
            }
        }
        stable_dt = 0.2 / std::sqrt(max_rate);
    }
    void apply()
    {
        max_gap = energy = balance = power_residual = torque_residual = power_relative = 0;
        Vecd sum = Vecd::Zero(), torque = Vecd::Zero();
        Real force_scale = 0, power = 0, power_scale = 0;
        for (size_t i = 0; i < shell.n; ++i) director_force[i].setZero();
        for (size_t i = 0; i < shell.n; ++i)
        {
            const auto &row = rows[i];
            Vecd mapped = Vecd::Zero(), mapped_v = Vecd::Zero();
            for (size_t k = 0; k < row.j.size(); ++k)
            { mapped += row.w[k] * core.x[row.j[k]]; mapped_v += row.w[k] * core.v[row.j[k]]; }
            Vecd inner_rate = Vecd::Zero();
            if (inner_surface)
                inner_rate = transform[i].transpose() * thin_structure_dynamics::getVectorChangeRateAfterThinStructureRotation(
                    Vecd(0, 0, 1), rotation[i], omega[i]);
            Vecd point = shell.x[i] - (inner_surface ? 0.5 * thickness : 0) * director[i];
            Vecd point_velocity = shell.v[i] - (inner_surface ? 0.5 * thickness : 0) * inner_rate;
            Vecd gap = point - mapped, f = -row.k * gap;
            shell_points[i] = point; mapped_points[i] = mapped; forces[i] = f;
            shell.prior[i] += f;
            Vecd g = -(inner_surface ? 0.5 * thickness : 0) * f;
            director_force[i] += g;
            sum += f; torque += shell.x[i].cross(f) + director[i].cross(g);
            Real ps = f.dot(shell.v[i]), pg = g.dot(inner_rate);
            power += ps + pg; power_scale += std::abs(ps) + std::abs(pg);
            for (size_t k = 0; k < row.j.size(); ++k)
            {
                Vecd reaction = -row.w[k] * f;
                core.prior[row.j[k]] += reaction;
                sum += reaction; power += reaction.dot(core.v[row.j[k]]);
                torque += core.x[row.j[k]].cross(reaction);
                power_scale += std::abs(reaction.dot(core.v[row.j[k]]));
            }
            Real spring_power = f.dot(point_velocity - mapped_v);
            power_residual += -spring_power; // dU/dt
            power_scale += std::abs(spring_power);
            max_gap = std::max(max_gap, gap.norm());
            energy += 0.5 * row.k * gap.squaredNorm(); force_scale += f.norm();
        }
        balance = sum.norm() / std::max(force_scale, Real(1e-30));
        power_residual += power; // external interface power + dU/dt = 0
        torque_residual = torque.norm() / std::max(length_scale * force_scale, Real(1e-30));
        power_relative = std::abs(power_residual) / std::max(power_scale, Real(1e-30));
        worst_balance = std::max(worst_balance, balance);
        worst_torque_residual = std::max(worst_torque_residual, torque_residual);
        worst_power_relative = std::max(worst_power_relative, power_relative);
    }
    void writeSurface(const fs::path &folder, size_t frame) const
    {
        fs::create_directories(folder);
        std::ofstream f(folder / ("interface_" + std::to_string(frame) + ".csv"));
        f << std::setprecision(17) << "x0_m,y0_m,z0_m,shell_x_m,shell_y_m,shell_z_m,core_x_m,core_y_m,core_z_m,fx_N,fy_N,fz_N\n";
        for (size_t i = 0; i < shell.n; ++i)
            f << reference_points[i][0] << ',' << reference_points[i][1] << ',' << reference_points[i][2] << ','
              << shell_points[i][0] << ',' << shell_points[i][1] << ',' << shell_points[i][2] << ','
              << mapped_points[i][0] << ',' << mapped_points[i][1] << ',' << mapped_points[i][2] << ','
              << forces[i][0] << ',' << forces[i][1] << ',' << forces[i][2] << '\n';
    }
    void patchTest(Real dp)
    {
        Matd R = Eigen::AngleAxisd(0.37, Vecd(1, 2, 3).normalized()).toRotationMatrix();
        Vecd shift(0.2 * dp, -0.1 * dp, 0.3 * dp);
        for (size_t i = 0; i < core.n; ++i) core.x[i] = R * core.x0[i] + shift;
        for (size_t i = 0; i < shell.n; ++i)
        { shell.x[i] = R * shell.x0[i] + shift; if (inner_surface) director[i] = R * initial_director[i]; }
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
        if (inner_surface)
            for (size_t i = 0; i < shell.n; ++i) torque += director[i].cross(director_force[i]);
        Real torque_residual = torque.norm() / std::max(dp * scale, Real(1e-30));
        if (balance > 1e-10 || torque_residual > 1e-10 || std::abs(power_residual) > 1e-10 * scale)
            throw std::runtime_error("force/torque/power patch failed");
        std::cout << "PATCH PASS force=" << balance << " torque=" << torque_residual << " power=" << power_residual << '\n';
        for (size_t i = 0; i < core.n; ++i) { core.x[i] = core.x0[i]; core.v[i].setZero(); }
        for (size_t i = 0; i < shell.n; ++i)
        { shell.x[i] = shell.x0[i]; shell.v[i].setZero(); if (inner_surface) director[i] = initial_director[i]; }
        core.clearPrior(); shell.clearPrior();
    }
};

struct PlaneContact
{
    State &state;
    Real offset, stiffness, stable_dt = std::numeric_limits<Real>::infinity(), upper = 0, lower = 0, penetration = 0, energy = 0;
    bool sampled, is_shell;
    Real radius, thickness, area = 0, affine_error = 0;
    Real force_residual = 0, torque_residual = 0, power_residual = 0;
    Real worst_force_residual = 0, worst_torque_residual = 0, worst_power_residual = 0;
    Vecd *director = nullptr, *director_force = nullptr, *rotation = nullptr, *omega = nullptr;
    Matd *transform = nullptr;
    StdVec<Vecd> x0, x, velocity, forces;
    std::vector<BondedInterface::Row> rows;
    PlaneContact(State &s, Real off, Real stiffness_per_area, Real depth,
                 bool sampled_surface, bool shell_surface, Real R, Real t, Real dp, Real mls_scale = 1, Real quadrature_scale = 1)
        : state(s), offset(off), stiffness(stiffness_per_area), sampled(sampled_surface),
          is_shell(shell_surface), radius(R), thickness(t)
    {
        if (sampled)
        {
            // Surface integration resolution is independent of solid particle
            // spacing and MLS support. Transpose scattering and per-area K
            // retain their definitions; quadrature area is recomputed below.
            auto directions = sphereDirections(R, quadrature_scale * dp);
            for (const auto &n : directions) x0.push_back(R * n);
            x.resize(x0.size()); velocity.resize(x0.size()); forces.resize(x0.size());
            area = 4 * Pi * R * R / x0.size();
            Real k = stiffness * area, max_rate = 0;
            if (is_shell)
            {
                if (s.n != x0.size()) throw std::runtime_error("shell/contact sampling mismatch");
                auto &p = s.body.getBaseParticles();
                director = p.getVariableDataByName<Vecd>("PseudoNormal");
                director_force = p.getVariableDataByName<Vecd>("ContactDirectorForce");
                rotation = p.getVariableDataByName<Vecd>("Rotation");
                omega = p.getVariableDataByName<Vecd>("AngularVelocity");
                transform = p.getVariableDataByName<Matd>("TransformationMatrix");
                for (size_t i = 0; i < s.n; ++i)
                {
                    affine_error = std::max(affine_error, (s.x0[i] + 0.5 * t * director[i] - x0[i]).norm() / dp);
                    // Translational mass plus director inertia m*t^2/12 gives
                    // k/m + (t/2)^2*k/(m*t^2/12) = 4*k/m.
                    max_rate = std::max(max_rate, 4 * k / s.mass[i]);
                }
            }
            else
            {
                std::vector<Real> rate(s.n, 0);
                for (const auto &target : x0)
                {
                    BondedInterface::Row row;
                    Eigen::Matrix4d moment = Eigen::Matrix4d::Zero();
                    std::vector<Eigen::Vector4d> bases;
                    std::vector<Real> kernels;
                    for (size_t j = 0; j < s.n; ++j)
                    {
                        // Controlled boundary reconstruction comparison only.
                        // Geometry, affine completion and work-conserving
                        // transpose force scatter keep the same definitions.
                        Vecd q = (s.x0[j] - target) / (mls_scale * dp);
                        if (q.norm() >= 2.5) continue;
                        Eigen::Vector4d b; b << 1, q[0], q[1], q[2];
                        Real w = std::exp(-q.squaredNorm());
                        moment += w * b * b.transpose();
                        row.j.push_back(j); bases.push_back(b); kernels.push_back(w);
                    }
                    Eigen::LDLT<Eigen::Matrix4d> solver(moment);
                    if (row.j.size() < 8 || solver.info() != Eigen::Success || solver.vectorD().minCoeff() <= 1e-10)
                        throw std::runtime_error("singular outer-surface mapping");
                    Eigen::Vector4d coeff = solver.solve(Eigen::Vector4d(1, 0, 0, 0));
                    Real sum = 0, sum_abs = 0; Vecd mapped = Vecd::Zero();
                    for (size_t a = 0; a < row.j.size(); ++a)
                    {
                        Real w = kernels[a] * bases[a].dot(coeff);
                        row.w.push_back(w); sum += w; sum_abs += std::abs(w);
                        mapped += w * s.x0[row.j[a]];
                    }
                    affine_error = std::max({affine_error, std::abs(sum - 1), (mapped - target).norm() / dp});
                    row.k = k;
                    for (size_t a = 0; a < row.j.size(); ++a)
                        rate[row.j[a]] += k * std::abs(row.w[a]) * sum_abs / s.mass[row.j[a]];
                    rows.push_back(std::move(row));
                }
                for (Real r : rate) max_rate = std::max(max_rate, r);
            }
            if (affine_error > 1e-9) throw std::runtime_error("outer surface does not reproduce initial sphere");
            stable_dt = 0.2 / std::sqrt(max_rate);
            return;
        }
        for (size_t i = 0; i < s.n; ++i)
        {
            Real k = stiffness * s.measure[i] / depth;
            stable_dt = std::min(stable_dt, Real(0.2) * std::sqrt(s.mass[i] / k));
        }
    }
    void apply(Real bottom, Real top, Real depth)
    {
        upper = lower = penetration = energy = 0;
        force_residual = torque_residual = power_residual = 0;
        if (director_force)
            for (size_t i = 0; i < state.n; ++i) director_force[i].setZero();
        if (sampled)
        {
            Vecd force_error = Vecd::Zero(), torque_error = Vecd::Zero();
            Real power_error = 0, force_scale = 0, power_scale = 0;
            for (size_t i = 0; i < x0.size(); ++i)
            {
                Vecd director_rate = Vecd::Zero();
                if (is_shell)
                {
                    director_rate = transform[i].transpose() * thin_structure_dynamics::getVectorChangeRateAfterThinStructureRotation(
                        Vecd(0, 0, 1), rotation[i], omega[i]);
                    x[i] = state.x[i] + 0.5 * thickness * director[i];
                    velocity[i] = state.v[i] + 0.5 * thickness * director_rate;
                }
                else
                {
                    x[i].setZero(); velocity[i].setZero();
                    const auto &row = rows[i];
                    for (size_t a = 0; a < row.j.size(); ++a)
                    { x[i] += row.w[a] * state.x[row.j[a]]; velocity[i] += row.w[a] * state.v[row.j[a]]; }
                }
                Real pb = std::max(bottom - x[i][2], Real(0)), pt = std::max(x[i][2] - top, Real(0));
                Real k = stiffness * area;
                Vecd f(0, 0, k * (pb - pt)); forces[i] = f;
                lower += k * pb; upper += k * pt;
                penetration = std::max({penetration, pb, pt});
                energy += 0.5 * k * (pb * pb + pt * pt);
                if (f.squaredNorm() == 0) continue;
                Vecd applied = Vecd::Zero(), torque = Vecd::Zero(); Real power = 0, absolute_body_power = 0;
                if (is_shell)
                {
                    Vecd g = 0.5 * thickness * f;
                    state.prior[i] += f; director_force[i] += g;
                    applied = f; torque = state.x[i].cross(f) + director[i].cross(g);
                    power = f.dot(state.v[i]) + g.dot(director_rate);
                    absolute_body_power = std::abs(f.dot(state.v[i])) + std::abs(g.dot(director_rate));
                }
                else
                {
                    const auto &row = rows[i];
                    for (size_t a = 0; a < row.j.size(); ++a)
                    {
                        size_t j = row.j[a]; Vecd fj = row.w[a] * f;
                        state.prior[j] += fj; applied += fj;
                        Real pj = fj.dot(state.v[j]);
                        torque += state.x[j].cross(fj); power += pj; absolute_body_power += std::abs(pj);
                    }
                }
                force_error += applied - f; torque_error += torque - x[i].cross(f);
                Real p = f.dot(velocity[i]); power_error += power - p;
                force_scale += f.norm(); power_scale += absolute_body_power + std::abs(p);
            }
            force_residual = force_error.norm() / std::max(force_scale, Real(1e-30));
            torque_residual = torque_error.norm() / std::max(radius * force_scale, Real(1e-30));
            power_residual = std::abs(power_error) / std::max(power_scale, Real(1e-30));
            worst_force_residual = std::max(worst_force_residual, force_residual);
            worst_torque_residual = std::max(worst_torque_residual, torque_residual);
            worst_power_residual = std::max(worst_power_residual, power_residual);
            return;
        }
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
    void writeSurface(const fs::path &folder, size_t frame) const
    {
        if (!sampled) return;
        fs::create_directories(folder);
        std::ofstream f(folder / ("surface_" + std::to_string(frame) + ".csv"));
        f << std::setprecision(17) << "x0_m,y0_m,z0_m,x_m,y_m,z_m,fx_N,fy_N,fz_N,reference_area_m2\n";
        for (size_t i = 0; i < x0.size(); ++i)
            f << x0[i][0] << ',' << x0[i][1] << ',' << x0[i][2] << ',' << x[i][0] << ',' << x[i][1] << ',' << x[i][2]
              << ',' << forces[i][0] << ',' << forces[i][1] << ',' << forces[i][2] << ',' << area << '\n';
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
    SPHSystem system(BoundingBoxd(Vecd::Constant(-2 * R), Vecd::Constant(2 * R)), o.spacing, o.threads);
    if (!o.patch_test) IO::getEnvironment().resetOutputFolder((output_root / "vtp").string(), true);
    SolidBody core(system, makeShared<GeometricShapeBall>(Vecd::Zero(), rc, "Core"));
    core.defineAdaptationRatios(1.15, 1);
    DruckerPragerCapSolid *dpc_material=nullptr;
    if (o.material == "dpc")
    {
        core.defineMatterMaterial<DruckerPragerCapSolid>(o.core_rho,dpc::Model::fromCSV(o.dpc_table,o.dpc_r0,o.dpc_shear_dilation,o.dpc_shear_hardening));
        dpc_material=&DynamicCast<DruckerPragerCapSolid>(&core,core.getMatterMaterial());
    }
    else if (o.material == "j2") core.defineMatterMaterial<HardeningPlasticSolid>(o.core_rho, o.core_E, o.core_nu, o.yield_stress, o.hardening);
    else core.defineMatterMaterial<SaintVenantKirchhoffSolid>(o.core_rho, o.core_E, o.core_nu);
    core.generateParticles<BaseParticles, CapsuleCoreGenerator>(rc, o.spacing);
    InnerRelation core_inner(core);
    InteractionWithUpdate<LinearGradientCorrectionMatrixInner> correction(core_inner);
    Dynamics1Level<solid_dynamics::Integration1stHalfPK2RightCauchy> elastic_first(core_inner);
    std::unique_ptr<solid_dynamics::SolidQuadraticReference> residual_reference;
    std::unique_ptr<solid_dynamics::SolidResidualIntegration1stHalf> residual_first;
    if(o.solid_stabilization=="quadratic_residual")
    {
        if(o.material!="dpc")throw std::invalid_argument("This example qualifies residual control with its DPC core only; the reusable operator remains material independent");
        residual_reference=std::make_unique<solid_dynamics::SolidQuadraticReference>(core_inner,o.core_E/(2*(1+o.core_nu)),o.core_E/(1-2*o.core_nu));
        residual_first=std::make_unique<solid_dynamics::SolidResidualIntegration1stHalf>(core_inner,*residual_reference);
        auto&p=core.getBaseParticles();p.addVariableToWrite<Real>("SolidResidualEnergyDensity");p.addVariableToWrite<Vecd>("SolidResidualForce");
    }
    std::unique_ptr<Dynamics1Level<solid_dynamics::DecomposedPlasticIntegration1stHalf>> plastic_first;
    if (o.material == "j2") plastic_first = std::make_unique<Dynamics1Level<solid_dynamics::DecomposedPlasticIntegration1stHalf>>(core_inner);
    Dynamics1Level<solid_dynamics::Integration2ndHalf> core_second(core_inner);
    ReduceDynamics<solid_dynamics::AcousticTimeStep> core_dt(core, o.cfl);
    State c(core);

    std::unique_ptr<SolidBody> shell;
    std::unique_ptr<InnerRelation> shell_inner;
    std::unique_ptr<InteractionDynamics<thin_structure_dynamics::ShellCorrectConfiguration>> shell_correction;
    std::unique_ptr<Dynamics1Level<SurfaceLoadedShellFirst>> shell_first;
    std::unique_ptr<Dynamics1Level<thin_structure_dynamics::ShellStressRelaxationSecondHalf>> shell_second;
    std::unique_ptr<ReduceDynamics<thin_structure_dynamics::ShellAcousticTimeStepSize>> shell_dt;
    std::unique_ptr<State> s;
    if (o.coated || o.patch_test)
    {
        shell = std::make_unique<SolidBody>(system, makeShared<DefaultShape>("ECShell"));
        shell->defineAdaptationRatios(1.15, 1);
        shell->defineMatterMaterial<SaintVenantKirchhoffSolid>(o.shell_rho, o.shell_E, o.shell_nu);
        shell->generateParticles<SurfaceParticles, CapsuleShellGenerator>(rs, o.spacing, o.thickness, o.plate_contact == "outer_surface");
        shell_inner = std::make_unique<InnerRelation>(*shell);
        shell_correction = std::make_unique<InteractionDynamics<thin_structure_dynamics::ShellCorrectConfiguration>>(*shell_inner);
        shell_first = std::make_unique<Dynamics1Level<SurfaceLoadedShellFirst>>(*shell_inner);
        shell_second = std::make_unique<Dynamics1Level<thin_structure_dynamics::ShellStressRelaxationSecondHalf>>(*shell_inner);
        shell_dt = std::make_unique<ReduceDynamics<thin_structure_dynamics::ShellAcousticTimeStepSize>>(*shell, o.cfl);
        s = std::make_unique<State>(*shell);
    }
    system.initializeSystemCellLinkedLists();
    system.initializeSystemConfigurations();
    correction.exec();if(residual_reference)residual_reference->build();
    if (shell_correction) shell_correction->exec();
    std::unique_ptr<BondedInterface> tie;
    if (s) tie = std::make_unique<BondedInterface>(c, *s, o.spacing, o.tie_factor * (o.constraint_E > 0 ? o.constraint_E : o.core_E) / o.spacing,
                                                 o.interface_surface == "inner", o.thickness, rc * rc / (rs * rs));
    if (o.patch_test) { tie->patchTest(o.spacing); return 0; }
    State &surface = s ? *s : c;
    Real contact_depth = s ? 1 : o.spacing;
    Real surface_offset = o.thickness / 2;
    if (!s && o.plate_contact == "legacy")
    {
        Real top_node = 0;
        for (const auto &x : c.x0) top_node = std::max(top_node, x[2]);
        surface_offset = R - top_node;
    }
    if (!s && o.plate_contact == "outer_surface") surface_offset = 0;
    PlaneContact contact(surface, surface_offset,
                         o.contact_factor * (o.constraint_E > 0 ? o.constraint_E : o.core_E) / o.spacing, contact_depth,
                         o.plate_contact == "outer_surface", bool(s), R, o.thickness, o.spacing, o.surface_mls_scale, o.contact_quadrature_scale);
    BodyStatesRecordingToVtp vtp(system);
    core.getBaseParticles().addVariableToWrite<Vecd>("Velocity");
    core.getBaseParticles().addVariableToWrite<Vecd>("ForcePrior");
    core.getBaseParticles().addVariableToWrite<Matd>("DeformationGradient");
    if (shell)
    {
        shell->getBaseParticles().addVariableToWrite<Vecd>("Velocity");
        shell->getBaseParticles().addVariableToWrite<Vecd>("PseudoNormal");
        shell->getBaseParticles().addVariableToWrite<Vecd>("ContactDirectorForce");
        shell->getBaseParticles().addVariableToWrite<Vecd>("InterfaceDirectorForce");
        shell->getBaseParticles().addVariableToWrite<Matd>("DeformationGradient");
        shell->getBaseParticles().addVariableToWrite<Matd>("BendingDeformationGradient");
        shell->getBaseParticles().addVariableToWrite<Matd>("MidSurfaceCauchyStress");
    }
    Real *plastic_alpha = nullptr;
    if (o.material == "j2")
    {
        plastic_alpha = core.getBaseParticles().getVariableDataByName<Real>("HardeningParameter");
        core.getBaseParticles().addVariableToWrite<Real>("HardeningParameter");
    }
    std::ofstream residual_history;
    if(residual_reference){residual_history.open(output_root/"solid_stabilization.csv");residual_history<<std::setprecision(17)<<"load_step,strain,residual_energy_J,core_elastic_energy_J,core_plastic_work_J,residual_fraction_of_material_energy\n";}
    std::ofstream material_history;
    Real *dpc_r=nullptr,*dpc_kappa=nullptr,*dpc_work=nullptr,*dpc_residual=nullptr,*dpc_surface=nullptr,*dpc_alpha=nullptr,*dpc_effective_d=nullptr;
    std::vector<bool> ever_cap(c.n,false);
    if(dpc_material)
    {
        auto &p=core.getBaseParticles();
        dpc_r=p.getVariableDataByName<Real>("DPCRelativeDensity");dpc_kappa=p.getVariableDataByName<Real>("DPCPlasticCompaction");
        dpc_work=p.getVariableDataByName<Real>("DPCPlasticStressWork");dpc_residual=p.getVariableDataByName<Real>("DPCYieldResidual");
        dpc_surface=p.getVariableDataByName<Real>("DPCActiveSurface");
        dpc_alpha=p.getVariableDataByName<Real>("DPCAccumulatedPlasticShear");dpc_effective_d=p.getVariableDataByName<Real>("DPCEffectiveCohesion");
        for(const std::string name:{"DPCRelativeDensity","DPCPlasticCompaction","DPCPlasticStressWork","DPCYieldResidual","DPCActiveSurface","DPCAccumulatedPlasticShear","DPCEffectiveCohesion"})p.addVariableToWrite<Real>(name);
        p.addVariableToWrite<Matd>("DPCCauchyStress");p.addVariableToWrite<Matd>("DPCInversePlasticMetric");
        fs::copy_file(o.dpc_table,output_root/"dpc_density_table.csv");
        material_history.open(output_root/"material_history.csv");
        material_history<<std::setprecision(17)<<"load_step,strain,volume_ratio,min_relative_density,max_relative_density,min_kappa,max_kappa,ever_cap_particle_fraction,plastic_stress_work_J,max_yield_residual,min_accumulated_shear,max_accumulated_shear,mean_accumulated_shear,min_effective_d_Pa,max_effective_d_Pa\n";
    }
    std::ofstream cfg(output_root / "parameters.txt");
    if (!cfg) throw std::runtime_error("cannot write parameters");
    if (!o.input_file.empty()) fs::copy_file(o.input_file, output_root / "input.cfg");
    cfg << std::setprecision(17) << "status=" << (o.quasistatic ? "CONFIGURED_ASSUMPTION_CASE_NOT_CALIBRATED" : "DEMONSTRATION_NOT_CALIBRATED")
        << "\nunits=SI\ncase=" << (s ? "coated" : "bare")
        << "\nsource_sha256=" << CAPSULE_SOURCE_SHA256
        << "\ndpc_material_sha256=" << DPC_MATERIAL_SHA256
        << "\ndpc_r0=" << o.dpc_r0 << "\ndpc_solid_density=" << o.dpc_solid_density
        << "\ndpc_flow=eta_scaled_shifted_ellipse_c_equals_b_pa\ndpc_shear_dilation=" << o.dpc_shear_dilation
        << "\ndpc_shear_hardening_Pa=" << o.dpc_shear_hardening
        << "\ndpc_shear_hardening_actual_Pa=" << (dpc_material?dpc_material->ShearHardeningModulus():0)
        << "\ndpc_reference_is_batch_calibrated=false"
        << "\nsolid_stabilization="<<o.solid_stabilization<<"\nsolid_residual_reference_dt_bound="<<(residual_reference?residual_reference->stable_dt:0)
        << "\ncore_material=" << o.material << "\nshell_material=" << (s ? "elastic_plane_stress" : "none")
        << "\ninterface=" << (s ? "affine_penalty_tie" : "none") << '\n'
        << "interface_surface=" << (s ? o.interface_surface : "none")
        << "\ninterface_affine_error=" << (tie ? tie->affine_error : 0)
        << "\ninterface_reference_area_ratio=" << (tie ? (tie->inner_surface ? tie->reference_area_ratio : 1) : 0) << '\n'
        << "plate_contact=" << o.plate_contact << "\ncontact_quadrature_points=" << contact.x0.size()
        << "\ncontact_quadrature_scale=" << o.contact_quadrature_scale
        << "\ncontact_quadrature_spacing=" << o.contact_quadrature_scale * o.spacing
        << "\nsurface_mls_scale=" << o.surface_mls_scale
        << "\nsurface_mls_support_radius=" << 2.5 * o.surface_mls_scale * o.spacing
        << "\nmaximum_dt_cap=" << o.max_dt
        << "\ncontact_reference_area_m2=" << contact.area << "\ncontact_affine_error=" << contact.affine_error << '\n'
        << "diameter=" << o.diameter << "\nthickness=" << o.thickness << "\nspacing=" << o.spacing
        << "\ncore_E=" << o.core_E << "\ncore_nu=" << o.core_nu << "\ncore_density=" << o.core_rho
        << "\nshell_E=" << o.shell_E << "\nshell_nu=" << o.shell_nu << "\nshell_density=" << o.shell_rho
        << "\nconstraint_reference_E=" << (o.constraint_E > 0 ? o.constraint_E : o.core_E)
        << (o.material == "dpc" ? "\nunused_for_dpc_yield_stress=" : "\nyield_stress=") << o.yield_stress << "\nhardening=" << o.hardening
        << (o.quasistatic ? "\nunused_dynamic_load_time=" : "\nload_time=") << o.load_time
        << (o.quasistatic ? "\nunused_dynamic_hold_time=" : "\nhold_time=") << o.hold_time << "\nstrain=" << o.strain
        << (o.quasistatic ? "\nunused_dynamic_maximum_plate_speed=" : "\nmaximum_plate_speed=")
        << 0.5 * Pi * o.strain * o.diameter / o.load_time
        << "\npaper_plate_speed=1e-6\ncontact_factor=" << o.contact_factor << "\ntie_factor=" << o.tie_factor
        << "\ncfl=" << o.cfl << "\ndt_scale=" << o.dt_scale << "\ndrag=" << o.drag
        << "\ncore_particles=" << c.n << "\nshell_particles=" << (s ? s->n : 0)
        << "\ncore_mass=" << c.totalMass() << "\nshell_mass=" << (s ? s->totalMass() : 0)
        << "\ncontact_surface_offset=" << surface_offset
        << "\nquasistatic=" << o.quasistatic << "\nload_steps=" << o.load_steps << "\nramp_steps=" << o.ramp_steps
        << "\nmax_stage_steps=" << o.max_stage_steps << "\nrelaxation_damping_per_step=" << o.relaxation_damping
        << "\nforce_tolerance=" << o.force_tolerance << "\nbalance_tolerance=" << o.balance_tolerance
        << "\nkinetic_tolerance=" << o.kinetic_tolerance << "\nthreads=" << o.threads
        << "\ntime_interpretation=" << (o.quasistatic ? "numerical_relaxation_time_not_experimental_time" : "accelerated_dynamic_time")
        << "\nplate_frame=" << (o.quasistatic ? "centered_static_equivalent_fixed_lower_z_equals_centered_z_minus_compression_over_2" : "fixed_lower_plate")
        << "\ncontact_dt_limit=" << contact.stable_dt << "\ntie_dt_limit=" << (tie ? tie->stable_dt : 0) << '\n';
    cfg.flush();
    std::ofstream csv(output_root / "history.csv");
    if (!csv) throw std::runtime_error("cannot write history");
    csv << "step," << (o.quasistatic ? "relaxation_time_s" : "time_s")
        << ",dt_s,compression_m,upper_N,lower_N,kinetic_J,core_elastic_J,tie_J,contact_J,max_tie_gap_m,max_penetration_m,min_J,interface_force_residual,interface_power_residual_W,external_work_J,drag_loss_J,max_plastic_alpha,interface_torque_residual,interface_power_relative\n";
    csv << std::setprecision(17);
    Real &time = *system.getSystemVariableDataByName<Real>("PhysicalTime");
    Real end = o.load_time + o.hold_time, next_csv = 0, next_vtp = 0, work = 0, drag_loss = 0;
    Real worst_J = 1, largest_gap = 0, largest_penetration = 0, max_force = 0;
    size_t step = 0, frame = 0;
    size_t stage = 1, stage_step = 0, consecutive_balanced = 0;
    std::ofstream equilibria;
    if (o.quasistatic)
    {
        equilibria.open(output_root / "equilibria.csv");
        if (!equilibria) throw std::runtime_error("cannot write equilibria");
        equilibria << std::setprecision(17)
                   << "load_step,relaxation_steps,compression_m,strain,experimental_displacement_time_s,upper_N,lower_N,force_residual,plate_imbalance,kinetic_ratio,max_tie_gap_m,min_J,max_plastic_alpha,interface_force_residual,interface_torque_residual,interface_power_relative\n";
    }
    auto compression = [&](Real t) { return 0.5 * o.strain * o.diameter * (1 - std::cos(Pi * std::min(t / o.load_time, Real(1)))); };
    auto step_compression = [&](Real increment)
    {
        if (!o.quasistatic) return compression(time + increment);
        Real a = std::min((Real(stage_step) + increment) / o.ramp_steps, Real(1));
        Real fraction = Real(stage - 1) + 0.5 * (1 - std::cos(Pi * a));
        return o.strain * o.diameter * fraction / o.load_steps;
    };
    // No gravity, spherical antipodal geometry, bonded isotropic materials and
    // frictionless equal planes: centered static loading is a rigid translation
    // of the fixed-lower-plate equilibrium, without a slow rigid drift mode.
    auto apply_contact = [&](Real delta)
    {
        Real center_shift = o.quasistatic ? 0.5 * delta : 0;
        contact.apply(-R + center_shift, R - delta + center_shift, contact_depth);
    };
    Real dpc_dt_bound=std::numeric_limits<Real>::infinity();
    auto inspect = [&](State &a, Real &kinetic, Real &U, Real &minJ, bool volume_core)
    {
        for (size_t i = 0; i < a.n; ++i)
        {
            Real J = a.F[i].determinant();
            if (!a.x[i].allFinite() || !a.v[i].allFinite() || !std::isfinite(J) || J <= 0.2)
                throw std::runtime_error("non-finite state or invalid deformation determinant");
            minJ = std::min(minJ, J);
            kinetic += 0.5 * a.mass[i] * a.v[i].squaredNorm();
            if(volume_core && dpc_material)
            {
                dpc_dt_bound=std::min(dpc_dt_bound,o.cfl*1.15*o.spacing/(dpc_material->acousticBound(J)+a.v[i].norm()));
                ever_cap[i]=ever_cap[i] || dpc_surface[i]==2;
            }
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
        Real delta = step_compression(0);
        apply_contact(delta);
        Real kinetic = 0, core_U = 0, minJ = 1;
        dpc_dt_bound=std::numeric_limits<Real>::infinity();
        inspect(c, kinetic, core_U, minJ, true); if (s) inspect(*s, kinetic, core_U, minJ, false);
        largest_gap = std::max(largest_gap, tie ? tie->max_gap : 0);
        largest_penetration = std::max(largest_penetration, contact.penetration);
        worst_J = std::min(worst_J, minJ); max_force = std::max(max_force, contact.upper);
        if (!o.quasistatic && time < end && end - time < 1e-12 * end) { time = end; continue; }
        Real dt = std::min({core_dt.exec(), contact.stable_dt,dpc_dt_bound,residual_reference?residual_reference->stable_dt:std::numeric_limits<Real>::infinity()});
        if (s) dt = std::min({dt, shell_dt->exec(), tie->stable_dt});
        if (o.drag > 0) dt = std::min(dt, Real(0.2) / o.drag);
        dt = std::min(dt, o.max_dt);
        dt = o.quasistatic ? o.dt_scale * dt : std::min(o.dt_scale * dt, end - time);
        if (time >= next_csv || (!o.quasistatic && time >= end))
        {
            Real max_alpha = plastic_alpha ? *std::max_element(plastic_alpha, plastic_alpha + c.n) : 0;
            csv << step << ',' << time << ',' << dt << ',' << delta << ',' << contact.upper << ',' << contact.lower << ','
                << kinetic << ',' << core_U << ',' << (tie ? tie->energy : 0) << ',' << contact.energy << ','
                << (tie ? tie->max_gap : 0) << ',' << contact.penetration << ',' << minJ << ','
                << (tie ? tie->balance : 0) << ',' << (tie ? tie->power_residual : 0) << ',' << work << ',' << drag_loss << ',' << max_alpha
                << ',' << (tie ? tie->torque_residual : 0) << ',' << (tie ? tie->power_relative : 0) << '\n';
            next_csv = time + end / 300;
        }
        if ((!o.quasistatic && (time >= next_vtp || time >= end)) || step == 0)
        {
            contact.writeSurface(output_root / "contact_surface", frame);
            if (tie) tie->writeSurface(output_root / "interface_surface", frame);
            vtp.writeToFile(frame++); next_vtp = time + end / 10;
            std::cout << "step=" << step << " time=" << time << " upper=" << contact.upper << " lower=" << contact.lower
                      << " gap=" << (tie ? tie->max_gap : 0) << " minJ=" << minJ << std::endl;
        }
        if (!o.quasistatic && time >= end) break;
        if (o.quasistatic && stage_step > o.ramp_steps && stage_step % 200 == 0)
        {
            Real sum_force2 = 0, rotational_kinetic = 0;
            for (size_t i = 0; i < c.n; ++i) sum_force2 += (c.force[i] + c.prior[i]).squaredNorm();
            if (s)
            {
                auto *omega = shell->getBaseParticles().getVariableDataByName<Vecd>("AngularVelocity");
                auto *alpha = shell->getBaseParticles().getVariableDataByName<Vecd>("AngularAcceleration");
                auto *rotation = shell->getBaseParticles().getVariableDataByName<Vecd>("Rotation");
                for (size_t i = 0; i < s->n; ++i)
                {
                    Real inertia = s->mass[i] * o.thickness * o.thickness / 12;
                    DirectorJacobian J = directorJacobian(rotation[i]);
                    Eigen::Vector2d w(omega[i][0], omega[i][1]), a(alpha[i][0], alpha[i][1]);
                    rotational_kinetic += 0.5 * inertia * (J * w).squaredNorm();
                    Eigen::Vector2d torque = inertia * (J.transpose() * J) * a;
                    sum_force2 += (s->force[i] + s->prior[i]).squaredNorm() + torque.squaredNorm() / (R * R);
                }
            }
            Real mean_force = 0.5 * (contact.upper + contact.lower);
            Real residual = std::sqrt(sum_force2) / std::max(mean_force, Real(1e-20));
            Real imbalance = std::abs(contact.upper - contact.lower) / std::max(mean_force, Real(1e-20));
            Real kinetic_ratio = (kinetic + rotational_kinetic) / std::max(mean_force * delta, Real(1e-30));
            bool balanced = mean_force > 1e-10 && residual <= o.force_tolerance &&
                            imbalance <= o.balance_tolerance && kinetic_ratio <= o.kinetic_tolerance;
            consecutive_balanced = balanced ? consecutive_balanced + 1 : 0;
            if (stage_step % 10000 == 0)
                std::cout << "load=" << stage << " relax=" << stage_step << " residual=" << residual
                          << " imbalance=" << imbalance << " kinetic=" << kinetic_ratio << std::endl;
            if (consecutive_balanced >= 3)
            {
                Real max_alpha = plastic_alpha ? *std::max_element(plastic_alpha, plastic_alpha + c.n) : 0;
                equilibria << stage << ',' << stage_step << ',' << delta << ',' << delta / o.diameter << ',' << delta / 1e-6 << ','
                           << contact.upper << ',' << contact.lower << ',' << residual << ',' << imbalance << ',' << kinetic_ratio << ','
                           << (tie ? tie->max_gap : 0) << ',' << minJ << ',' << max_alpha << ',' << (tie ? tie->balance : 0)
                           << ',' << (tie ? tie->torque_residual : 0) << ',' << (tie ? tie->power_relative : 0) << '\n';
                equilibria.flush();
                if(dpc_material)
                {
                    Real V=0,V0=0,W=0,yr=0;size_t caps=0;
                    for(size_t i=0;i<c.n;++i){V+=c.measure[i]*c.F[i].determinant();V0+=c.measure[i];W+=c.measure[i]*dpc_work[i];yr=std::max(yr,dpc_residual[i]);caps+=ever_cap[i];}
                    material_history<<stage<<','<<delta/o.diameter<<','<<V/V0<<','
                        <<*std::min_element(dpc_r,dpc_r+c.n)<<','<<*std::max_element(dpc_r,dpc_r+c.n)<<','
                        <<*std::min_element(dpc_kappa,dpc_kappa+c.n)<<','<<*std::max_element(dpc_kappa,dpc_kappa+c.n)<<','
                        <<Real(caps)/c.n<<','<<W<<','<<yr<<','<<*std::min_element(dpc_alpha,dpc_alpha+c.n)<<','
                        <<*std::max_element(dpc_alpha,dpc_alpha+c.n)<<','<<std::accumulate(dpc_alpha,dpc_alpha+c.n,Real(0))/c.n<<','
                        <<*std::min_element(dpc_effective_d,dpc_effective_d+c.n)<<','<<*std::max_element(dpc_effective_d,dpc_effective_d+c.n)<<'\n';material_history.flush();
                }
                if(residual_reference)
                {Real U=0,W=0;auto*p=core.getBaseParticles().getVariableDataByName<Real>("DPCElasticEnergyDensity");for(size_t i=0;i<c.n;++i){U+=c.measure[i]*p[i];if(dpc_work)W+=c.measure[i]*dpc_work[i];}Real S=residual_reference->totalEnergy();residual_history<<stage<<','<<delta/o.diameter<<','<<S<<','<<U<<','<<W<<','<<S/std::max(U+W,Real(1e-30))<<'\n';residual_history.flush();}
                contact.writeSurface(output_root / "contact_surface", frame);
                if (tie) tie->writeSurface(output_root / "interface_surface", frame);
                vtp.writeToFile(frame++);
                std::cout << "ACCEPTED load=" << stage << '/' << o.load_steps << " strain=" << delta / o.diameter
                          << " force=" << mean_force << " residual=" << residual << " imbalance=" << imbalance << std::endl;
                if (stage == o.load_steps) break;
                ++stage; stage_step = 0; consecutive_balanced = 0;
                continue;
            }
            if (stage_step >= o.max_stage_steps) throw std::runtime_error("load step did not converge; no complete-case result written");
        }
        if (step >= o.max_steps || !std::isfinite(dt) || dt < 1e-15) throw std::runtime_error("time step/step-count guard");
        // Split the native first-half operators at the SAME midpoint for both
        // bodies. Lagging penalty forces by a half drift injects interface energy.
        auto initialize = [&](auto &method, State &a)
        { particle_for(execution::ParallelPolicy(), a.body.LoopRange(), [&](size_t i) { method.initialization(i, dt); }); };
        auto interact = [&](auto &method, State &a)
        { particle_for(execution::ParallelPolicy(), a.body.LoopRange(), [&](size_t i) { method.interaction(i, dt); }); };
        auto update = [&](auto &method, State &a)
        { particle_for(execution::ParallelPolicy(), a.body.LoopRange(), [&](size_t i) { method.update(i, dt); }); };
        if(residual_first)initialize(*residual_first,c);else if (plastic_first) initialize(*plastic_first, c); else initialize(elastic_first, c);
        // Native shell initialization writes a shared numerical damping matrix.
        // Keep this portion serial; core and shell interactions remain parallel.
        if (shell_first)
            for (size_t i = 0; i < s->n; ++i) shell_first->initialization(i, dt);
        c.clearPrior(); if (s) s->clearPrior();
        if (tie) tie->apply();
        apply_contact(step_compression(o.quasistatic ? Real(0.5) : 0.5 * dt));
        if (o.drag > 0)
        {
            auto damp = [&](State &a) { for (size_t i = 0; i < a.n; ++i) a.prior[i] -= o.drag * a.mass[i] * a.v[i]; };
            damp(c); if (s) damp(*s);
            drag_loss += 2 * o.drag * kinetic * dt;
        }
        work += (o.quasistatic ? 0.5 * (contact.upper + contact.lower) : contact.upper) *
                (step_compression(o.quasistatic ? Real(1) : dt) - delta);
        if(residual_first)interact(*residual_first,c);else if (plastic_first) interact(*plastic_first, c); else interact(elastic_first, c);
        if (shell_first) interact(*shell_first, *s);
        if(residual_first)update(*residual_first,c);else if (plastic_first) update(*plastic_first, c); else update(elastic_first, c);
        if (shell_first) update(*shell_first, *s);
        if (o.quasistatic)
        {
            Real factor = std::exp(-o.relaxation_damping), before = 0, after = 0;
            auto damp_velocity = [&](State &a)
            {
                for (size_t i = 0; i < a.n; ++i)
                { before += 0.5 * a.mass[i] * a.v[i].squaredNorm(); a.v[i] *= factor; after += 0.5 * a.mass[i] * a.v[i].squaredNorm(); }
            };
            damp_velocity(c); if (s) damp_velocity(*s);
            if (s)
            {
                auto *omega = shell->getBaseParticles().getVariableDataByName<Vecd>("AngularVelocity");
                auto *rotation = shell->getBaseParticles().getVariableDataByName<Vecd>("Rotation");
                for (size_t i = 0; i < s->n; ++i)
                {
                    Real inertia = s->mass[i] * o.thickness * o.thickness / 12;
                    Eigen::Vector2d w(omega[i][0], omega[i][1]);
                    Real energy = 0.5 * inertia * (directorJacobian(rotation[i]) * w).squaredNorm();
                    before += energy; omega[i] *= factor; after += energy * factor * factor;
                }
            }
            drag_loss += before - after;
        }
        core_second.exec(dt); if (shell_second) shell_second->exec(dt);
        time = o.quasistatic ? time + dt : std::min(end, time + dt); ++step; ++stage_step;
    }
    if (o.rest_test && (max_force > 1e-8 * o.core_E * o.diameter * o.diameter ||
                        worst_J < 1 - 1e-6 || largest_gap > 1e-6 * o.spacing))
        throw std::runtime_error("zero-load equilibrium regression failed");
    if (o.smoke_test && (max_force <= 1e-10 || worst_J < 0.9 || largest_gap > 0.02 * o.spacing))
        throw std::runtime_error("loaded solid-shell smoke regression failed");
    std::ofstream summary(output_root / "summary.json");
    if (!summary) throw std::runtime_error("cannot write summary");
    summary << std::setprecision(17) << "{\n  \"status\": \"" << (o.quasistatic ? "COMPLETED_CONFIGURED_QUASISTATIC_CASE" : "COMPLETED_DEMONSTRATION")
            << "\",\n  \"steps\": " << step
            << ",\n  \"equilibrated_load_steps\": " << (o.quasistatic ? stage : 0)
            << ",\n  \"time_is_numerical_relaxation\": " << (o.quasistatic ? "true" : "false")
            << ",\n  \"source_sha256\": \"" << CAPSULE_SOURCE_SHA256 << "\""
            << (o.quasistatic ? ",\n  \"numerical_relaxation_time_s\": " : ",\n  \"final_time_s\": ")
            << time << ",\n  \"min_J\": " << worst_J << ",\n  \"max_upper_N\": " << max_force
            << ",\n  \"max_tie_gap_m\": " << largest_gap << ",\n  \"max_penetration_m\": " << largest_penetration
            << ",\n  \"plate_contact\": \"" << o.plate_contact << "\""
            << ",\n  \"contact_quadrature_points\": " << contact.x0.size()
            << ",\n  \"max_contact_force_mapping_residual\": " << contact.worst_force_residual
            << ",\n  \"max_contact_torque_mapping_residual\": " << contact.worst_torque_residual
            << ",\n  \"max_contact_power_mapping_residual\": " << contact.worst_power_residual
            << ",\n  \"interface_surface\": \"" << (s ? o.interface_surface : "none") << "\""
            << ",\n  \"max_interface_force_mapping_residual\": " << (tie ? tie->worst_balance : 0)
            << ",\n  \"max_interface_torque_mapping_residual\": " << (tie ? tie->worst_torque_residual : 0)
            << ",\n  \"max_interface_power_mapping_residual\": " << (tie ? tie->worst_power_relative : 0)
            << ",\n  \"experimental_calibration\": false,\n  \"shell_energy_and_rotational_energy_audited\": false\n}\n";
    return 0;
}

int main(int ac, char **av)
{
    try { Options o; o.parse(ac, av); return run(o); }
    catch (const std::exception &e) { std::cerr << "ERROR: " << e.what() << std::endl; return 1; }
}
