/** Bounded normal-contact prototype probe. No PD body/solver is constructed. */
#include "pd_normal_projection.h"
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <unistd.h>

using namespace SPH;
using namespace SPH::peridynamics;
namespace
{
void writeJsonString(const std::string &value)
{
    std::cout << '"';
    for (unsigned char c : value)
        if (c == '"' || c == '\\') std::cout << '\\' << c;
        else if (c < 32) std::cout << "\\u00" << std::hex << std::setw(2) << std::setfill('0') << int(c) << std::dec << std::setfill(' ');
        else std::cout << c;
    std::cout << '"';
}
void array(const StdVec<Real> &values)
{
    std::cout << '[';
    for (size_t i = 0; i < values.size(); ++i) { if (i) std::cout << ','; std::cout << values[i]; }
    std::cout << ']';
}
}

int main(int argc, char **argv)
{
    std::cout << std::setprecision(17);
    // Identity is captured by the native process from its actual /proc state.
    // The wrapper also records its own OS identity and the immutable ELF SHA.
    std::ifstream boot("/proc/sys/kernel/random/boot_id"), stat("/proc/self/stat"), cmdline("/proc/self/cmdline", std::ios::binary);
    std::string boot_id, fields, token, ticks, command;
    std::getline(boot, boot_id); std::getline(stat, fields);
    std::istringstream tail(fields.substr(fields.rfind(')') + 2));
    for (int i = 0; i <= 19; ++i) { tail >> token; if (i == 19) ticks = token; }
    std::cout << "{\"identity\":{\"pid\":" << getpid() << ",\"ppid\":" << getppid() << ",\"start_ticks\":";
    writeJsonString(ticks); std::cout << ",\"boot_id\":"; writeJsonString(boot_id);
    std::cout << ",\"executable\":"; writeJsonString(std::filesystem::read_symlink("/proc/self/exe").string());
    std::cout << ",\"argv\":[";
    bool first = true;
    while (std::getline(cmdline, command, '\0')) { if (!first) std::cout << ','; writeJsonString(command); first = false; }
    std::cout << "]},\"prototype_only\":true,";
    try
    {
        if (argc != 2) throw std::invalid_argument("usage: normal_contact_probe INPUT");
        std::ifstream input(argv[1]);
        size_t count, rows; Real dt; NormalProjectionSettings settings;
        if (!(input >> count >> rows >> dt >> settings.velocity_scale >> settings.relative_tolerance >> settings.max_sweeps) || count == 0 || count > 128 || rows > 1024)
            throw std::invalid_argument("invalid bounded probe dimensions/header");
        StdVec<Real> mass(count), position(count), velocity(count);
        for (auto *values : {&mass, &position, &velocity}) for (auto &value : *values) input >> value;
        StdVec<NormalContactConstraint> constraints(rows);
        for (auto &row : constraints)
        {
            size_t width;
            if (!(input >> row.normal >> row.old_offset >> row.new_offset >> width) || width == 0 || width > count)
                throw std::invalid_argument("invalid bounded row width");
            row.targets.resize(width); row.weights.resize(width);
            for (size_t j = 0; j < width; ++j) input >> row.targets[j] >> row.weights[j];
        }
        if (!input) throw std::invalid_argument("truncated prototype input");
        input >> std::ws; if (!input.eof()) throw std::invalid_argument("unexpected prototype input tail");
        const auto r = projectNormalContact(mass, position, velocity, dt, constraints, settings);
        std::cout << "\"converged\":" << (r.converged ? "true" : "false") << ",\"sweeps\":" << r.sweeps << ",\"velocity\":"; array(r.velocity);
        std::cout << ",\"position\":"; array(r.position); std::cout << ",\"normal_impulse\":"; array(r.normal_impulse);
        std::cout << ",\"primal_velocity_error\":" << r.primal_velocity_error << ",\"natural_velocity_error\":" << r.natural_velocity_error
                  << ",\"stationarity_impulse_error\":" << r.stationarity_impulse_error << ",\"complementarity_energy_error\":" << r.complementarity_energy_error
                  << ",\"free_kinetic\":" << r.free_kinetic << ",\"accepted_kinetic\":" << r.accepted_kinetic
                  << ",\"wall_work\":" << r.wall_work << ",\"contact_step_dissipation\":" << r.contact_step_dissipation
                  << ",\"energy_residual\":" << r.energy_residual << ",\"momentum_residual\":" << r.momentum_residual << "}\n";
        return r.converged ? 0 : 2;
    }
    catch (const std::exception &error)
    {
        std::cout << "\"error\":"; writeJsonString(error.what()); std::cout << "}\n";
        return 1;
    }
}
