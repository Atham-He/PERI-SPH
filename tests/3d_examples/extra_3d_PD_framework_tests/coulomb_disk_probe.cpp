/** Frozen-input arithmetic observer. It never advances a physical state. */
#include "pd_coulomb_disk.h"
#include <cerrno>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <sys/stat.h>
#include <unistd.h>
using namespace SPH;
using namespace SPH::peridynamics;
namespace
{
void emitString(const std::string &s)
{
    std::cout<<'"';
    for (unsigned char c:s)
        if (c=='"' || c=='\\') std::cout<<'\\'<<c;
        else if (c<32) std::cout<<"\\u00"<<std::hex<<std::setw(2)<<std::setfill('0')<<int(c)<<std::dec<<std::setfill(' ');
        else std::cout<<c;
    std::cout<<'"';
}
void finish()
{
    std::cout<<"}\n";std::cout.flush();
    if (!std::cout) throw std::runtime_error("Coulomb observer: output failure");
    struct stat s;
    if (::fstat(STDOUT_FILENO,&s)!=0 ||
        (S_ISREG(s.st_mode) && ::fsync(STDOUT_FILENO)!=0))
        throw std::runtime_error("Coulomb observer: durable output failure");
}
void identity()
{
    std::ifstream boot("/proc/sys/kernel/random/boot_id"),stat("/proc/self/stat"),cmd("/proc/self/cmdline",std::ios::binary);
    std::string id,line,word,ticks;std::getline(boot,id);std::getline(stat,line);
    std::istringstream fields(line.substr(line.rfind(')')+2));
    for (int k=0;k<=19;++k) { fields>>word;if (k==19) ticks=word; }
    std::cout<<"{\"type\":\"identity\",\"pid\":"<<getpid()<<",\"ppid\":"<<getppid()<<",\"start_ticks\":";emitString(ticks);
    std::cout<<",\"boot_id\":";emitString(id);std::cout<<",\"argv\":[";
    bool first=true;
    while (std::getline(cmd,word,'\0')) { if (!first) std::cout<<',';first=false;emitString(word); }
    std::cout<<']';finish();
}
template<class Derived>void vector(const Eigen::MatrixBase<Derived> &v)
{
    std::cout<<'[';
    for (Eigen::Index i=0;i<v.size();++i) { if(i)std::cout<<',';std::cout<<v[i]; }
    std::cout<<']';
}
void matrix(const ContactTangentMatrix &m)
{
    std::cout<<"[["<<m(0,0)<<','<<m(0,1)<<"],["<<m(1,0)<<','<<m(1,1)<<"]]";
}
StdVec<Real> readInput(const char *path)
{
    std::ifstream input(path);StdVec<Real> result;std::string token;
    if (!input) throw std::invalid_argument("Coulomb observer: input file unavailable");
    while (input>>token)
    {
        char *end=nullptr;errno=0;
        const Real value=std::strtod(token.c_str(),&end);
        if (end!=token.c_str()+token.size() || end==token.c_str())
            throw std::invalid_argument("Coulomb observer: invalid input token");
        if (errno==ERANGE && !std::isfinite(value))
            throw std::invalid_argument("Coulomb observer: unrepresentable input token");
        result.push_back(value);
    }
    if (result.size()!=6) throw std::invalid_argument("Coulomb observer: expected six input fields");
    return result;
}
}
int main(int argc,char **argv)
{
    std::cout<<std::setprecision(17)<<std::boolalpha;std::ios::sync_with_stdio(false);
    try
    {
        if (argc==2 && std::string(argv[1])=="--capabilities")
        {
            std::cout<<"{\"type\":\"capabilities\",\"format\":1,\"kernel_only\":true,\"arithmetic_precision_bits\":256,"
                "\"envelope_binary_exponent\":-240,\"roundoff_gate\":1e-12,\"full_friction_global_solver_qualified\":false";
            finish();return 0;
        }
        identity();
        if (argc!=2) throw std::invalid_argument("Coulomb observer: usage INPUT");
        const auto c=readInput(argv[1]);
        const ContactTangentVector impulse(c[1],c[2]),velocity(c[3],c[4]);
        const auto r=coulombDiskNaturalMap(c[0],impulse,velocity,c[5]);
        std::cout<<"{\"type\":\"disk_input\",\"normal_impulse\":"<<c[0]<<",\"tangent_impulse\":";vector(impulse);
        std::cout<<",\"relative_velocity\":";vector(velocity);std::cout<<",\"mu\":"<<c[5];finish();
        const char *branch=r.branch==CoulombDiskBranch::StickInterior ? "stick_interior" :
            (r.branch==CoulombDiskBranch::SlipExterior ? "slip_exterior" : "zero_radius_origin");
        std::cout<<"{\"type\":\"disk_value\",\"branch\":";emitString(branch);
        std::cout<<",\"kernel_only\":true,\"accepted_physical_step\":false,\"representable\":"<<r.representable
            <<",\"within_budget\":"<<r.within_budget<<",\"residual\":";vector(r.residual);
        std::cout<<",\"normal_derivative\":";vector(r.normal_derivative);
        std::cout<<",\"impulse_derivative\":";matrix(r.impulse_derivative);
        std::cout<<",\"velocity_derivative\":";matrix(r.velocity_derivative);
        std::cout<<",\"residual_roundoff_bound\":";vector(r.residual_roundoff_bound);
        std::cout<<",\"maximum_derivative_roundoff_bound\":"<<r.maximum_derivative_roundoff_bound
            <<",\"friction_loss\":"<<r.friction_loss<<",\"friction_loss_roundoff_bound\":"<<r.friction_loss_roundoff_bound
            <<",\"cone_excess_upper\":"<<r.cone_excess_upper
            <<",\"maximum_dissipation_error_upper\":"<<r.maximum_dissipation_error_upper;finish();
        std::cout<<"{\"type\":\"end\",\"completed\":true,\"arithmetic_budget_satisfied\":"<<r.within_budget
            <<",\"accepted_steps\":0,\"last_accepted_time\":0,\"full_S3_complete\":false";finish();
        return r.within_budget ? 0 : 2;
    }
    catch (const std::exception &error)
    {
        std::cout<<"{\"type\":\"error\",\"accepted_physical_step\":false,\"message\":";emitString(error.what());finish();return 1;
    }
}
