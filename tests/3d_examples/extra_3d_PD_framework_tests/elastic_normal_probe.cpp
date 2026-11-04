/** Bounded normal-contact prototype probe. No PD body/solver is constructed. */
#include "pd_elastic_normal_step.h"
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
        if (argc!=2) throw std::invalid_argument("usage: elastic_normal_probe INPUT");
        std::ifstream input(argv[1]);
        Real mass, stiffness, q, v, displacement, duration, dt, velocity_scale;
        if (!(input>>mass>>stiffness>>q>>v>>displacement>>duration>>dt>>velocity_scale))
            throw std::invalid_argument("truncated elastic probe input");
        input>>std::ws;
        if (!input.eof() || !std::isfinite(duration) || duration<=0 || !std::isfinite(dt) || dt<=0 ||
            !std::isfinite(displacement) || displacement<=0 || duration/dt>10000)
            throw std::invalid_argument("invalid bounded elastic trajectory");
        const size_t count=std::llround(duration/dt);
        if (!count || std::abs(count*dt-duration)>1e-12*duration)
            throw std::invalid_argument("duration must contain an integer number of steps");
        auto wall=[&](Real t) { const Real s=t/duration; return displacement*s*s*s*(10+s*(-15+6*s)); };
        NormalProjectionSettings settings; settings.velocity_scale=velocity_scale;
        std::cout<<"\"model\":\"diagonal_elastic_normal_v1\",\"steps\":[";
        for (size_t step=1; step<=count; ++step)
        {
            const Real t0=(step-1)*dt, t1=step*dt, a0=wall(t0), a1=wall(t1);
            const auto r=advanceElasticNormal({mass},{stiffness},{q},{v},dt,{{{0},{1},1,a0,a1}},settings);
            if (!r.projection.converged) throw std::runtime_error("projection did not converge; state not committed");
            if (step>1) std::cout<<',';
            std::cout<<"{\"step\":"<<step<<",\"time_old\":"<<t0<<",\"time_new\":"<<t1
                <<",\"dt\":"<<dt<<",\"q_old\":"<<q<<",\"v_old\":"<<v
                <<",\"wall_old\":"<<a0<<",\"wall_new\":"<<a1
                <<",\"q_new\":"<<r.projection.position[0]<<",\"v_new\":"<<r.projection.velocity[0]
                <<",\"impulse\":"<<r.projection.normal_impulse[0]
                <<",\"reaction_interval\":"<<r.projection.normal_impulse[0]/dt
                <<",\"elastic_force\":"<<r.elastic_force[0]
                <<",\"kinetic_old\":"<<r.kinetic_old<<",\"kinetic_new\":"<<r.kinetic_new
                <<",\"elastic_old\":"<<r.elastic_old<<",\"elastic_new\":"<<r.elastic_new
                <<",\"wall_work\":"<<r.wall_work<<",\"algorithmic_loss\":"<<r.algorithmic_loss
                <<",\"energy_residual\":"<<r.energy_residual
                <<",\"maximum_momentum_residual\":"<<r.maximum_momentum_residual
                <<",\"primal_velocity_error\":"<<r.projection.primal_velocity_error
                <<",\"natural_velocity_error\":"<<r.projection.natural_velocity_error
                <<",\"stationarity_impulse_error\":"<<r.projection.stationarity_impulse_error
                <<",\"complementarity_energy_error\":"<<r.projection.complementarity_energy_error
                <<",\"sweeps\":"<<r.projection.sweeps<<'}';
            q=r.projection.position[0]; v=r.projection.velocity[0];
        }
        std::cout<<"],\"completed\":true}\n";
        return 0;
    }
    catch (const std::exception &error)
    {
        // A partial trajectory is retained as raw failure evidence, never accepted.
        std::cerr<<error.what()<<'\n';
        return 1;
    }
}
