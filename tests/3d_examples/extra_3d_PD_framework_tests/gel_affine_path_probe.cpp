/** Bounded material-point research path. No particle body or contact solver. */
#include "pd_gel_affine_path.h"
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
void matrix(const Matd &value)
{
    std::cout << '[';
    for (size_t i = 0; i < 9; ++i) { if (i) std::cout << ','; std::cout << value(i / 3, i % 3); }
    std::cout << ']';
}
void readMatrix(std::istream &input, Matd &value)
{ for (size_t i = 0; i < 9; ++i) input >> value(i / 3, i % 3); }
}
int main(int argc, char **argv)
{
    std::cout << std::setprecision(17);
    std::ifstream boot("/proc/sys/kernel/random/boot_id"), stat("/proc/self/stat"), cmdline("/proc/self/cmdline", std::ios::binary);
    std::string boot_id, fields, token, ticks, command;
    std::getline(boot, boot_id); std::getline(stat, fields);
    std::istringstream tail(fields.substr(fields.rfind(')') + 2));
    for (int i = 0; i <= 19; ++i) { tail >> token; if (i == 19) ticks = token; }
    std::cout << "{\"type\":\"identity\",\"pid\":" << getpid() << ",\"ppid\":" << getppid() << ",\"start_ticks\":";
    writeJsonString(ticks); std::cout << ",\"boot_id\":"; writeJsonString(boot_id);
    std::cout << ",\"argv\":["; bool first = true;
    while (std::getline(cmdline, command, '\0')) { if (!first) std::cout << ','; first = false; writeJsonString(command); }
    std::cout << "]}\n";
    size_t accepted = 0;
    try
    {
        if (argc != 2) throw std::invalid_argument("usage: gel_affine_path_probe INPUT");
        std::ifstream input(argv[1]); GelParameters p; GelHistory history; Matd current;
        input >> p.bulk >> p.equilibrium_shear >> p.branch_shear >> p.relaxation_time;
        readMatrix(input, current); readMatrix(input, history.green); readMatrix(input, history.viscous);
        size_t count; input >> history.dissipation >> count;
        if (!input || !count || count > 64) throw std::invalid_argument("invalid bounded material header");
        StdVec<Real> increments(count); StdVec<Matd> endpoints(count);
        for (size_t i = 0; i < count; ++i) { input >> increments[i]; readMatrix(input, endpoints[i]); }
        if (!input) throw std::invalid_argument("truncated material path");
        input >> std::ws; if (!input.eof()) throw std::invalid_argument("unexpected material input tail");
        for (size_t step = 0; step < count; ++step)
        {
            const auto r = advanceGelAffinePath(current, endpoints[step], increments[step], p, history);
            if (std::abs(r.energy_residual) > 1e-10 * std::max(p.bulk, p.equilibrium_shear))
                throw std::runtime_error("material path work/storage/dissipation gate failed before commit");
            std::cout << "{\"type\":\"step\",\"step\":" << step + 1 << ",\"dt\":" << increments[step] << ",\"F0\":";
            matrix(current); std::cout << ",\"F1\":"; matrix(endpoints[step]);
            std::cout << ",\"green\":"; matrix(r.history.green); std::cout << ",\"viscous\":"; matrix(r.history.viscous);
            std::cout << ",\"average_pk1\":"; matrix(r.average_pk1); std::cout << ",\"tangent\":[";
            for (size_t j = 0; j < 9; ++j) { if (j) std::cout << ','; matrix(r.endpoint_tangent[j]); }
            std::cout << "],\"tangent_available\":" << (r.tangent_available ? "true" : "false")
                << ",\"stored_old\":" << r.stored_old << ",\"stored_new\":" << r.stored_new
                << ",\"physical_dissipation\":" << r.physical_dissipation << ",\"total_dissipation\":" << r.history.dissipation
                << ",\"work\":" << r.work << ",\"energy_residual\":" << r.energy_residual
                << ",\"minimum_determinant\":" << r.minimum_determinant << ",\"tangent_asymmetry\":" << r.tangent_asymmetry << "}\n";
            current = endpoints[step]; history = r.history; ++accepted;
        }
        std::cout << "{\"type\":\"terminal\",\"completed\":true,\"accepted_steps\":" << accepted << ",\"prototype_only\":true}\n";
        return 0;
    }
    catch (const std::exception &error)
    {
        std::cout << "{\"type\":\"error\",\"accepted_steps\":" << accepted << ",\"error\":";
        writeJsonString(error.what()); std::cout << "}\n"; return 1;
    }
}
