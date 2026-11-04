/** Original simultaneous algebra controls, not an actual PD body. */
#include "coulomb_contact_probe_io.h"
using namespace SPH;using namespace SPH::peridynamics;
namespace io=bounded_probe_io;namespace fio=coupled_probe_io;
namespace
{
struct Input
{
    int mode=0,fault=0;
    size_t n=0,m=0;
    Real dt=1,mass_scale=1,domain_limit=1e6;
    ManyContactSettings settings;
    StdVec<Real> old_position,load,beta,trial_velocity,trial_impulses,linear_rhs;
    DenseContactMatrix matrix;
    StdVec<CoulombContactConstraint> contacts;
};
Input readInput(const char *path)
{
    std::ifstream stream(path);Input c;
    stream>>c.mode>>c.n>>c.m>>c.dt>>c.mass_scale>>c.settings.velocity_scale>>c.domain_limit
        >>c.settings.max_newton_steps>>c.settings.max_line_search>>c.fault;
    if(!stream || c.mode<0 || c.mode>2 || !c.n || c.n>128 || c.n%3 || c.m>256 || c.fault<0 || c.fault>3 ||
        !std::isfinite(c.domain_limit) || c.domain_limit<=0)
        throw std::invalid_argument("coupled observer: configuration");
    auto read=[&](StdVec<Real> &v,size_t n){v.resize(n);for(auto &x:v)stream>>x;};
    read(c.old_position,c.n);c.matrix.assign(c.n,StdVec<Real>(c.n));
    for(auto &row:c.matrix)for(auto &x:row)stream>>x;
    read(c.load,c.n);read(c.beta,c.n);c.contacts.resize(c.m);
    for(auto &row:c.contacts)
    {
        read(row.weights,c.n/3);for(auto &x:row.normal)stream>>x;
        for(auto &t:row.tangent)for(auto &x:t)stream>>x;
        stream>>row.old_offset>>row.new_offset>>row.mu;
        for(auto &x:row.wall_tangent_velocity)stream>>x;
    }
    read(c.trial_velocity,c.n);read(c.trial_impulses,3*c.m);read(c.linear_rhs,c.n+3*c.m);
    if(!stream)throw std::invalid_argument("coupled observer: truncated input");
    stream>>std::ws;if(!stream.eof())throw std::invalid_argument("coupled observer: input tail");
    for(const auto &row:c.matrix)for(Real x:row)if(!std::isfinite(x))throw std::invalid_argument("coupled observer: nonfinite matrix");
    for(const auto *v:{&c.load,&c.beta,&c.linear_rhs})
        for(Real x:*v)if(!std::isfinite(x))throw std::invalid_argument("coupled observer: nonfinite operator input");
    return c;
}
NonlinearContactEquation evaluate(const Input &c,const StdVec<Real> &v)
{
    if(c.fault==1 || (c.fault==2 && std::any_of(v.begin(),v.end(),[](Real x){return x!=0;})))
        throw std::runtime_error("declared frozen material evaluation fault");
    NonlinearContactEquation r;r.momentum.resize(c.n);r.jacobian=c.matrix;
    for(size_t i=0;i<c.n;++i)
    {
        r.momentum[i]=denseContactDot(c.matrix[i],v)-c.load[i]+c.beta[i]*v[i]*v[i]*v[i];
        r.jacobian[i][i]+=3*c.beta[i]*v[i]*v[i];
    }
    if(c.fault==3)r.momentum[0]=std::numeric_limits<Real>::infinity();
    return r;
}
void header(const Input &c,const CoulombContactSystem &s)
{
    io::begin("header");io::field("format",1);io::field("actual_PD_body",false);io::field("manufactured_global_control",true);
    io::field("mode",c.mode);io::field("fault",c.fault);io::field("primal_dofs",c.n);io::field("contact_count",c.m);
    io::field("unknowns",s.unknowns);io::field("old_position",c.old_position);
    io::field("dt",c.dt);io::field("mass_scale",c.mass_scale);io::field("velocity_scale",s.speed);io::field("impulse_scale",s.impulse_scale);
    io::field("domain_limit",c.domain_limit);io::field("max_newton_steps",c.settings.max_newton_steps);
    io::field("max_line_search",c.settings.max_line_search);io::field("relative_tolerance",c.settings.relative_tolerance);
    io::field("operator_matrix",c.matrix);io::field("operator_load",c.load);io::field("operator_cubic_beta",c.beta);
    io::field("operator_formula","R_i=A_i*v-load_i+beta_i*v_i^3; J=A+diag(3beta_i*v_i^2)");
    io::field("C_full_three_components",s.coefficients());std::cout<<",\"contacts\":";fio::contacts(c.contacts);
    io::field("trial_velocity",c.trial_velocity);io::field("trial_impulses",c.trial_impulses);
    io::field("linear_rhs",c.linear_rhs);io::field("accepted_physical_steps",0);io::field("S3_complete",false);io::finish();
}
}
int main(int argc,char **argv)
{
    std::cout<<std::setprecision(17);std::ios::sync_with_stdio(false);
    try
    {
        if(argc==2 && std::string(argv[1])=="--capabilities")
        {
            io::begin("capabilities");io::field("format",1);io::field("original_coupled_global_observer",true);
            io::field("actual_PD_body",false);io::field("full_S3_complete",false);io::finish();return 0;
        }
        io::identity();if(argc!=2)throw std::invalid_argument("coupled observer: usage INPUT");
        const auto c=readInput(argv[1]);const CoulombContactSystem system(c.old_position,c.dt,c.mass_scale,c.contacts,c.settings);
        header(c,system);
        const auto material=[&](const auto &v){return evaluate(c,v);};
        int code=0;
        if(c.mode==0)
        {
            const auto r=solveCoulombContact(system,material,[&](const auto &v)
                {return std::all_of(v.begin(),v.end(),[&](Real x){return std::isfinite(x) && std::abs(x)<=c.domain_limit;});},c.settings);
            io::begin("global_solve");io::field("accepted_physical_step",false);std::cout<<",\"solve\":";fio::solve(r);io::finish();
            code=r.converged ? 0 : 2;
        }
        else
        {
            const auto values=system.evaluate(c.trial_velocity,c.trial_impulses,material,true);
            io::begin("global_sample");io::field("accepted_physical_step",false);
            std::cout<<",\"original\":";fio::values(values);
            if(c.mode==2 && values.precision_budget_ok)
            {
                try
                {
                    const auto linear=solveBoundedConsistentContactSystem(values.jacobian,c.linear_rhs,1024);
                    io::field("linear_consistent",true);io::field("linear_solution",linear.value);
                    io::field("linear_rank",linear.rank);io::field("linear_backward_residual",linear.relative_residual);
                }
                catch(const std::runtime_error &error)
                { io::field("linear_consistent",false);io::field("linear_reason",error.what());code=2; }
            }
            if(!values.precision_budget_ok)code=2;
            io::finish();
        }
        io::begin("end");io::field("completed",code==0);io::field("accepted_physical_steps",0);io::field("last_accepted_time",0);
        io::field("full_S3_complete",false);io::finish();return code;
    }
    catch(const std::exception &error)
    {
        io::begin("error");io::field("message",error.what());io::field("accepted_physical_step",false);io::finish();return 1;
    }
}
