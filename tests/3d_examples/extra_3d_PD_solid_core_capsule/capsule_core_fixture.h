#pragma once
/** Shared core lattice, SI input parsing and state views; no SPH shell implementation. */
#include "sphinxsys.h"
#include "drucker_prager_cap_solid.h"
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

struct PDSolidCapsuleOptions
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
    std::string material = "elastic", output = "run", input_file, plate_contact = "legacy", interface_surface = "legacy_mid";
    size_t max_steps = 2000000;
    std::string dpc_table;
    Real dpc_r0 = 0, dpc_solid_density = 0, constraint_E = 0, dpc_shear_dilation = 1;

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
        if (material == "dpc")
        {
            auto model=dpc::Model::fromCSV(dpc_table,dpc_r0,dpc_shear_dilation);auto p=model.parameters(dpc_r0);
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
            (contact_quadrature_scale != 1 && plate_contact != "outer_surface"))
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
} // namespace SPH

struct State
{
    RealBody &body;
    size_t n;
    Vecd *x, *v, *prior, *force;
    Real *mass, *measure;
    Matd *F;
    StdVec<Vecd> x0;
    explicit State(RealBody &b) : body(b)
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

