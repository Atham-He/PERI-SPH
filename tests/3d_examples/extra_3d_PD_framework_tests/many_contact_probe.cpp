/** Manufactured many-row algebra, deliberately not a gel/sphere driver. */
#include "pd_many_contact.h"
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <unistd.h>
using namespace SPH;
using namespace SPH::peridynamics;
namespace
{
void jsonString(const std::string &value)
{
    std::cout << '"';
    for (unsigned char c:value)
        if(c=='"'||c=='\\')std::cout<<'\\'<<c;
        else if(c<32)std::cout<<"\\u00"<<std::hex<<std::setw(2)<<std::setfill('0')<<int(c)<<std::dec<<std::setfill(' ');
        else std::cout<<c;
    std::cout << '"';
}
void array(const StdVec<Real> &v)
{ std::cout<<'[';for(size_t i=0;i<v.size();++i){if(i)std::cout<<',';std::cout<<v[i];}std::cout<<']'; }
void matrix(const DenseContactMatrix &m)
{ std::cout<<'[';for(size_t i=0;i<m.size();++i){if(i)std::cout<<',';array(m[i]);}std::cout<<']'; }
void identity()
{
    std::ifstream boot("/proc/sys/kernel/random/boot_id"),stat("/proc/self/stat"),cmdline("/proc/self/cmdline",std::ios::binary);
    std::string id,fields,word,ticks;std::getline(boot,id);std::getline(stat,fields);
    std::istringstream tail(fields.substr(fields.rfind(')')+2));for(int i=0;i<=19;++i){tail>>word;if(i==19)ticks=word;}
    std::cout<<"{\"identity\":{\"pid\":"<<getpid()<<",\"ppid\":"<<getppid()<<",\"start_ticks\":";jsonString(ticks);
    std::cout<<",\"boot_id\":";jsonString(id);std::cout<<",\"argv\":[";bool first=true;
    while(std::getline(cmdline,word,'\0')){if(!first)std::cout<<',';first=false;jsonString(word);}std::cout<<"]},\"manufactured_operator_only\":true";
}
}
int main(int argc,char **argv)
{
    std::cout<<std::setprecision(17);identity();
    try
    {
        if(argc!=2)throw std::invalid_argument("usage: many_contact_probe INPUT");
        std::ifstream input(argv[1]);int mode=0;size_t n=0;input>>mode>>n;
        if(!input||!n||n>128)throw std::invalid_argument("invalid model dimension");
        DenseContactMatrix a(n,StdVec<Real>(n));
        if(mode==1)
        {
            StdVec<Real> rhs(n);for(auto &row:a)for(auto &v:row)input>>v;for(auto &v:rhs)input>>v;
            if(!input)throw std::invalid_argument("truncated linear input");input>>std::ws;if(!input.eof())throw std::invalid_argument("linear tail");
            const auto r=solveConsistentContactSystem(a,rhs);
            std::cout<<",\"linear_rank\":"<<r.rank<<",\"linear_residual\":"<<r.relative_residual<<",\"value\":";array(r.value);
            std::cout<<",\"status\":\"consistent_original_linear_system_solved\"}\n";return 0;
        }
        if(mode!=0)throw std::invalid_argument("unknown manufactured model");
        Real dt,mass_scale,domain_radius;ManyContactSettings settings;
        input>>dt>>mass_scale>>settings.velocity_scale>>settings.relative_tolerance>>settings.max_newton_steps>>settings.max_line_search>>domain_radius;
        StdVec<Real> q(n),beta(n),load(n);
        for(auto *values:{&q,&beta,&load})for(auto &v:*values)input>>v;
        for(auto &row:a)for(auto &v:row)input>>v;
        size_t count=0;input>>count;if(count>256)throw std::invalid_argument("too many manufactured rows");
        StdVec<DenseContactConstraint> rows(count);
        for(auto &row:rows){input>>row.old_offset>>row.new_offset;row.coefficients.resize(n);for(auto &v:row.coefficients)input>>v;}
        if(!input)throw std::invalid_argument("truncated manufactured input");input>>std::ws;if(!input.eof())throw std::invalid_argument("manufactured tail");
        if(!std::isfinite(domain_radius)||domain_radius<=0)throw std::invalid_argument("manufactured domain radius");
        for(size_t i=0;i<n;++i)
        {
            if(!std::isfinite(beta[i])||beta[i]<0||!std::isfinite(load[i]))throw std::invalid_argument("manufactured coefficients");
            for(Real v:a[i])if(!std::isfinite(v))throw std::invalid_argument("manufactured matrix");
        }
        auto equation=[&](const StdVec<Real> &v)
        {
            NonlinearContactEquation e;e.momentum.resize(n);e.jacobian=a;
            for(size_t i=0;i<n;++i){e.momentum[i]=denseContactDot(a[i],v)+beta[i]*v[i]*v[i]*v[i]-load[i];e.jacobian[i][i]+=3*beta[i]*v[i]*v[i];}
            return e;
        };
        auto in_domain=[&](const StdVec<Real> &v){for(Real x:v)if(!std::isfinite(x)||std::abs(x)>domain_radius)return false;return true;};
        const auto r=solveManyNonsymmetricContact(q,dt,mass_scale,rows,equation,in_domain,settings);
        std::cout<<",\"converged\":"<<(r.converged?"true":"false")<<",\"status\":";jsonString(r.status);
        std::cout<<",\"normalized_KKT\":"<<r.normalized_KKT<<",\"velocity\":";array(r.velocity);
        std::cout<<",\"normal_impulse\":";array(r.normal_impulse);
        if(r.converged){const auto e=equation(r.velocity);std::cout<<",\"momentum_equation\":";array(e.momentum);std::cout<<",\"jacobian\":";matrix(e.jacobian);}
        std::cout<<",\"trace\":[";
        for(size_t k=0;k<r.trace.size();++k)
        {
            if(k)std::cout<<',';const auto &t=r.trace[k];
            std::cout<<"{\"iteration\":"<<t.iteration<<",\"linear_rank\":"<<t.linear_rank<<",\"domain_rejections\":"<<t.domain_rejections
                <<",\"merit\":"<<t.merit<<",\"squared_norm\":"<<t.squared_norm<<",\"KKT\":"<<t.KKT<<",\"alpha\":"<<t.accepted_step
                <<",\"linear_residual\":"<<t.linear_residual<<",\"velocity\":";array(t.velocity);
            std::cout<<",\"normal_impulse\":";array(t.normal_impulse);std::cout<<",\"direction\":";array(t.normalized_direction);
            std::cout<<",\"normalized_residual\":";array(t.normalized_residual);std::cout<<",\"normalized_gaps\":";array(t.normalized_gaps);
            std::cout<<",\"generalized_jacobian\":";matrix(t.generalized_jacobian);std::cout<<'}';
        }
        std::cout<<"]}\n";return r.converged?0:2;
    }
    catch(const std::exception &e){std::cout<<",\"error\":";jsonString(e.what());std::cout<<"}\n";return 1;}
}
