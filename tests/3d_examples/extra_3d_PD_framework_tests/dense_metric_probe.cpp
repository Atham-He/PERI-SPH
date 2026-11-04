/** Bounded full-matrix contact algebra probe; no PD body is constructed. */
#include "pd_dense_metric_contact.h"
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
    std::ifstream boot("/proc/sys/kernel/random/boot_id"), stat("/proc/self/stat"), cmdline("/proc/self/cmdline", std::ios::binary);
    std::string boot_id, fields, token, ticks, command;
    std::getline(boot, boot_id); std::getline(stat, fields);
    std::istringstream tail(fields.substr(fields.rfind(')') + 2));
    for (int i = 0; i <= 19; ++i) { tail >> token; if (i == 19) ticks = token; }
    std::cout << "{\"identity\":{\"pid\":" << getpid() << ",\"ppid\":" << getppid() << ",\"start_ticks\":";
    writeJsonString(ticks); std::cout << ",\"boot_id\":"; writeJsonString(boot_id);
    std::cout << ",\"argv\":[";
    bool first = true;
    while (std::getline(cmdline, command, '\0')) { if (!first) std::cout << ','; first = false; writeJsonString(command); }
    std::cout << "]},\"prototype_only\":true";
    try
    {
        if (argc != 2) throw std::invalid_argument("usage: dense_metric_probe INPUT");
        std::ifstream input(argv[1]);
        size_t n, rows; int mode; Real dt; NormalProjectionSettings settings;
        if (!(input >> mode >> n >> rows >> dt >> settings.velocity_scale >> settings.relative_tolerance >> settings.max_sweeps) ||
            (mode != 0 && mode != 1) || !n || n > 128 || rows > 1024)
            throw std::invalid_argument("invalid bounded input header");
        DenseContactMatrix matrix(n, StdVec<Real>(n));
        StdVec<Real> mass(n), q(n), v(n);
        for (auto &row : matrix) for (auto &value : row) input >> value;
        for (auto *vector : {&mass, &q, &v}) for (auto &value : *vector) input >> value;
        StdVec<DenseContactConstraint> constraints(rows);
        for (auto &row : constraints)
        {
            input >> row.old_offset >> row.new_offset; row.coefficients.resize(n);
            for (auto &value : row.coefficients) input >> value;
        }
        if (!input) throw std::invalid_argument("truncated input");
        input >> std::ws; if (!input.eof()) throw std::invalid_argument("unexpected input tail");
        DenseElasticContactResult elastic;
        const DenseMetricContactResult r = mode == 0 ? projectDenseMetricContact(matrix, q, v, dt, constraints, settings) :
            (elastic = advanceDenseElasticContact(mass, matrix, q, v, dt, constraints, settings)).projection;
        std::cout << ",\"converged\":" << (r.converged ? "true" : "false") << ",\"sweeps\":" << r.sweeps;
        std::cout << ",\"velocity\":"; array(r.velocity); std::cout << ",\"normal_impulse\":"; array(r.normal_impulse);
        std::cout << ",\"linear_relative_residual\":" << r.linear_relative_residual
            << ",\"primal_velocity_error\":" << r.primal_velocity_error << ",\"natural_velocity_error\":" << r.natural_velocity_error
            << ",\"stationarity_impulse_error\":" << r.stationarity_impulse_error
            << ",\"complementarity_energy_error\":" << r.complementarity_energy_error
            << ",\"physical_ledger_present\":" << (mode == 1 ? "true" : "false");
        if (mode == 1)
        {
            std::cout << ",\"position\":"; array(elastic.position); std::cout << ",\"internal_force\":"; array(elastic.internal_force);
            std::cout << ",\"momentum_residual\":"; array(elastic.momentum_residual);
            std::cout << ",\"kinetic_old\":" << elastic.kinetic_old << ",\"kinetic_new\":" << elastic.kinetic_new
                << ",\"elastic_old\":" << elastic.elastic_old << ",\"elastic_new\":" << elastic.elastic_new
                << ",\"algorithmic_loss\":" << elastic.algorithmic_loss << ",\"normal_contact_loss\":" << elastic.normal_contact_loss
                << ",\"wall_work\":" << elastic.wall_work << ",\"energy_residual\":" << elastic.energy_residual;
        }
        std::cout << "}\n"; return r.converged ? 0 : 2;
    }
    catch (const std::exception &error)
    {
        std::cout << ",\"error\":"; writeJsonString(error.what()); std::cout << "}\n"; return 1;
    }
}
