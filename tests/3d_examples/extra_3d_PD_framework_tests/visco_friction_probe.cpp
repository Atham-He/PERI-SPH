/** Bounded two-coordinate research driver. No PD body is constructed. */
#include "pd_visco_friction_step.h"
#include <array>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <unistd.h>

using namespace SPH;
using namespace SPH::peridynamics;
static_assert(sizeof(Real) == 8, "S2 checkpoint requires binary64 Real");
namespace
{
const std::string header = "step,t0,t1,dt,q0,vn0,x0,qt0,vt0,q1,vn1,x1,qt1,vt1,wall0,wall1,wall_ut,jn,jt,fn,ft,material_force,U0,U1,K0,K1,Wwall,Wbody,Dmaterial,Dfriction,Dnormal,Dalgorithm,Eres,Pnres,Ptres,primal,natural,stationarity,complementarity,cone,friction_natural,slip,mode";

void writeJsonString(const std::string &value)
{
    std::cout << '"';
    for (unsigned char c : value)
        if (c == '"' || c == '\\') std::cout << '\\' << c;
        else if (c < 32) std::cout << "\\u00" << std::hex << std::setw(2) << std::setfill('0') << int(c) << std::dec << std::setfill(' ');
        else std::cout << c;
    std::cout << '"';
}

void identity()
{
    std::ifstream boot("/proc/sys/kernel/random/boot_id"), stat("/proc/self/stat"), command("/proc/self/cmdline", std::ios::binary);
    std::string boot_id, fields, word, ticks;
    std::getline(boot, boot_id); std::getline(stat, fields);
    std::istringstream tail(fields.substr(fields.rfind(')') + 2));
    for (int i = 0; i <= 19; ++i) { tail >> word; if (i == 19) ticks = word; }
    std::cout << "{\"identity\":{\"pid\":" << getpid() << ",\"ppid\":" << getppid() << ",\"start_ticks\":";
    writeJsonString(ticks); std::cout << ",\"boot_id\":"; writeJsonString(boot_id);
    std::cout << ",\"argv\":[";
    bool first = true;
    while (std::getline(command, word, '\0')) { if (!first) std::cout << ','; first = false; writeJsonString(word); }
    std::cout << "]}";
}

struct Sum
{
    Real value = 0, correction = 0;
    void add(Real x)
    {
        const Real y = x - correction, next = value + y;
        correction = (next - value) - y; value = next;
    }
};

struct History
{
    ViscoFrictionState body;
    uint64_t step = 0, cursor = 0;
    std::array<Sum, 7> sums{}; // Work, material, friction, normal, algorithm, Jn, Jt.
    Real jn = 0, jt = 0, wall = 0, wall_ut = 0;
    int32_t mode = 0;
};

template<class T> void binary(std::ostream &stream, const T &value)
{ stream.write(reinterpret_cast<const char *>(&value), sizeof(T)); }
template<class T> void binary(std::istream &stream, T &value)
{ stream.read(reinterpret_cast<char *>(&value), sizeof(T)); }

template<class Stream> void bodyFields(Stream &stream, History &s)
{
    binary(stream, s.body.normal_position); binary(stream, s.body.normal_velocity); binary(stream, s.body.maxwell_extension);
    binary(stream, s.body.tangent_position); binary(stream, s.body.tangent_velocity);
    binary(stream, s.step); binary(stream, s.cursor);
    for (auto &sum : s.sums) { binary(stream, sum.value); binary(stream, sum.correction); }
    binary(stream, s.jn); binary(stream, s.jt); binary(stream, s.wall); binary(stream, s.wall_ut); binary(stream, s.mode);
}

void checkpoint(const std::filesystem::path &path, int32_t kind, const std::array<Real, 15> &config, History history)
{
    const auto temporary = path.string() + ".tmp";
    std::ofstream output(temporary, std::ios::binary);
    const uint64_t magic = 0x5356463200000001ULL;
    binary(output, magic); binary(output, kind);
    for (Real x : config) binary(output, x);
    bodyFields(output, history);
    output.flush();
    if (!output) throw std::runtime_error("checkpoint write failed");
    output.close(); std::filesystem::rename(temporary, path);
}

History restore(const std::filesystem::path &path, int32_t kind, const std::array<Real, 15> &config)
{
    std::ifstream input(path, std::ios::binary);
    uint64_t magic{}; int32_t saved_kind{};
    binary(input, magic); binary(input, saved_kind);
    if (!input) throw std::invalid_argument("truncated checkpoint header");
    if (magic != 0x5356463200000001ULL || saved_kind != kind) throw std::invalid_argument("checkpoint format/model mismatch");
    for (Real expected : config)
    {
        Real value; binary(input, value);
        if (value != expected) throw std::invalid_argument("checkpoint immutable config mismatch");
    }
    History history; bodyFields(input, history);
    if (!input || input.peek() != std::char_traits<char>::eof()) throw std::invalid_argument("checkpoint size mismatch");
    return history;
}
}

int main(int argc, char **argv)
{
    std::cout << std::setprecision(17); identity();
    try
    {
        if (argc < 4) throw std::invalid_argument("usage: visco_friction_probe INPUT TRACE CHECKPOINT [--restore CP] [--stop STEP]");
        std::ifstream input(argv[1]);
        int32_t kind; std::array<Real, 15> config;
        if (!(input >> kind) || kind < 0 || kind > 5) throw std::invalid_argument("invalid model kind");
        for (auto &x : config) if (!(input >> x) || !std::isfinite(x)) throw std::invalid_argument("invalid config");
        input >> std::ws; if (!input.eof()) throw std::invalid_argument("unexpected config tail");
        const Real dt = config[0], total = config[1], mn = config[2], mt = config[3], displacement = config[7];
        const ScalarMaxwellParameters material{config[4], config[5], config[6]};
        material.validate();
        if (dt <= 0 || total <= 0 || total / dt > 10000 || mn <= 0 || mt <= 0 || config[8] < 0)
            throw std::invalid_argument("invalid bounded trajectory");
        const uint64_t count = std::llround(total / dt);
        if (!count || std::abs(count * dt - total) > 1e-12 * total) throw std::invalid_argument("noninteger step count");
        uint64_t stop = count;
        std::filesystem::path saved;
        for (int i = 4; i < argc; i += 2)
        {
            if (i + 1 >= argc) throw std::invalid_argument("truncated option");
            if (std::string(argv[i]) == "--restore") saved = argv[i + 1];
            else if (std::string(argv[i]) == "--stop") stop = std::stoull(argv[i + 1]);
            else throw std::invalid_argument("unknown option");
        }
        History history;
        if (!saved.empty()) history = restore(saved, kind, config);
        else
        {
            history.body.normal_position = config[11]; history.body.maxwell_extension = config[12];
            history.body.normal_velocity = config[13]; history.body.tangent_velocity = config[14] + config[9];
        }
        if (stop > count || history.step > stop) throw std::invalid_argument("invalid restart/stop cursor");
        const std::filesystem::path trace = argv[2];
        if (!saved.empty())
        {
            if (std::filesystem::file_size(trace) < history.cursor) throw std::invalid_argument("trace shorter than checkpoint cursor");
            std::filesystem::resize_file(trace, history.cursor); // Roll back only explicitly supplied new-study journal tail.
        }
        else if (std::filesystem::exists(trace)) throw std::invalid_argument("fresh trace already exists");
        std::ofstream output(trace, std::ios::binary | std::ios::app);
        if (saved.empty()) { output << header << '\n'; history.cursor = output.tellp(); }
        auto wall = [&](Real t)
        {
            if (kind == 0) return t <= 1 ? displacement * t : t <= 2 ? displacement : t <= 3 ? displacement * (3 - t) : 0;
            if (kind == 1) return displacement;
            if (kind == 4) return Real(0);
            if (kind == 5) return -std::max(Real(0), t - 8 * dt);
            const Real s = t / total;
            return displacement * s * s * s * (10 + s * (-15 + 6 * s));
        };
        for (uint64_t step = history.step + 1; step <= stop; ++step)
        {
            const auto old = history.body;
            const Real t0 = (step - 1) * dt, t1 = step * dt, middle = (t0 + t1) / 2;
            const Real w0 = wall(t0), w1 = wall(t1);
            const Real first_switch = kind == 4 ? 0.2 : 0.6, second_switch = kind == 4 ? 0.5 : 1.2;
            const Real ut = config[9] + (middle < first_switch ? 0.05 : middle < second_switch ? -0.025 : 0.1);
            ViscoFrictionResult r{};
            if (kind <= 1)
            {
                r.state = old; r.state.normal_position = w1; r.state.normal_velocity = (w1 - w0) / dt;
                r.material = advanceScalarMaxwell(old.normal_position, old.maxwell_extension, w1, dt, material);
                r.state.maxwell_extension = r.material.extension_new;
                // Prescribed material control has no inertia/contact, so actuator work is exact material work.
                r.wall_work = r.material.average_force * (w1 - w0);
                r.energy_residual = r.wall_work - (r.material.stored_new - r.material.stored_old) - r.material.physical_dissipation;
            }
            else r = advanceViscoFriction(old, mn, mt, material, dt, w0, w1, ut, config[8], config[10]);
            const Real jn = kind <= 1 ? 0 : r.normal.normal_impulse[0], jt = r.tangent.impulse;
            const std::array<Real, 43> fields = {Real(step), t0, t1, dt, old.normal_position, old.normal_velocity,
                old.maxwell_extension, old.tangent_position, old.tangent_velocity,
                r.state.normal_position, r.state.normal_velocity, r.state.maxwell_extension,
                r.state.tangent_position, r.state.tangent_velocity, w0, w1, ut,
                jn, jt, jn / dt, jt / dt, r.material.average_force, r.material.stored_old, r.material.stored_new,
                r.kinetic_old, r.kinetic_new, r.wall_work, r.body_work, r.material.physical_dissipation,
                r.friction_loss, r.normal_contact_loss, r.algorithmic_loss, r.energy_residual,
                r.momentum_normal, r.momentum_tangent, r.normal.primal_velocity_error, r.normal.natural_velocity_error,
                r.normal.stationarity_impulse_error, r.normal.complementarity_energy_error,
                r.tangent.cone_residual, r.tangent.natural_impulse_residual, r.tangent.slip, Real(r.tangent.mode)};
            output << std::setprecision(17);
            for (size_t j = 0; j < fields.size(); ++j)
            {
                if (!std::isfinite(fields[j])) throw std::runtime_error("nonfinite accepted trace field");
                if (j) output << ','; output << fields[j];
            }
            output << '\n';
            if (!output) throw std::runtime_error("accepted trace write failed");
            history.body = r.state; history.step = step; history.cursor = output.tellp();
            history.jn = jn; history.jt = jt; history.wall = w1; history.wall_ut = ut; history.mode = r.tangent.mode;
            const std::array<Real, 7> ledger = {r.wall_work + r.body_work, r.material.physical_dissipation,
                r.friction_loss, r.normal_contact_loss, r.algorithmic_loss, jn, jt};
            for (size_t j = 0; j < ledger.size(); ++j) history.sums[j].add(ledger[j]);
        }
        output.flush();
        checkpoint(argv[3], kind, config, history);
        std::cout << ",\"completed\":true,\"accepted_steps\":" << history.step << ",\"committed_bytes\":" << history.cursor
                  << ",\"prototype_only\":true}\n";
        return 0;
    }
    catch (const std::exception &error)
    {
        std::cout << ",\"error\":"; writeJsonString(error.what()); std::cout << "}\n";
        return 1;
    }
}
