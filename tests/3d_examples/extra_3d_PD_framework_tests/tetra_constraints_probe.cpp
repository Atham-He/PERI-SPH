/** Bounded nonlinear global/reference-gradient controls, not a sphere driver. */
#include "pd_tetra_gel_constraints.h"
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
    std::cout << '['; for(size_t i=0;i<values.size();++i) { if(i)std::cout<<',';std::cout<<values[i]; } std::cout<<']';
}
void dense(const DenseContactMatrix &matrix)
{
    std::cout<<'[';for(size_t i=0;i<matrix.size();++i){if(i)std::cout<<',';array(matrix[i]);}std::cout<<']';
}
void tensor(const Matd &matrix)
{
    StdVec<Real> values(9);for(size_t i=0;i<9;++i)values[i]=matrix(i/3,i%3);array(values);
}
void vector3(const Vecd &v) { array({v[0],v[1],v[2]}); }
void solver(const PrimalDualResult &r)
{
    std::cout<<",\"converged\":"<<(r.converged?"true":"false")<<",\"status\":";writeJsonString(r.status);
    std::cout<<",\"active_mask\":"<<r.active_mask<<",\"dependent_sets_skipped\":"<<r.dependent_sets_skipped
        <<",\"normalized_KKT\":"<<r.normalized_KKT<<",\"maximum_linear_residual\":"<<r.maximum_linear_residual<<",\"trace\":[";
    for(size_t i=0;i<r.trace.size();++i)
    {
        if(i)std::cout<<',';const auto &t=r.trace[i];
        std::cout<<"{\"mask\":"<<t.active_mask<<",\"iteration\":"<<t.iteration<<",\"domain_rejections\":"<<t.domain_rejections
            <<",\"merit\":"<<t.merit<<",\"alpha\":"<<t.accepted_step<<",\"linear_residual\":"<<t.linear_residual<<'}';
    }
    std::cout<<']';
}
void evaluation(const TetraGelEvaluation &e)
{
    std::cout<<",\"momentum_equation\":";array(e.equation.momentum);
    std::cout<<",\"jacobian\":";dense(e.equation.jacobian);
    std::cout<<",\"internal_force\":";array(e.internal_force);
    std::cout<<",\"F1\":";tensor(e.deformation);
    std::cout<<",\"average_pk1\":";tensor(e.material.average_pk1);
    std::cout<<",\"green\":";tensor(e.material.history.green);std::cout<<",\"viscous\":";tensor(e.material.history.viscous);
    std::cout<<",\"total_dissipation\":"<<e.material.history.dissipation
        <<",\"material_stored_old\":"<<e.material.stored_old<<",\"material_stored_new\":"<<e.material.stored_new
        <<",\"material_loss\":"<<e.material.physical_dissipation<<",\"material_work\":"<<e.material.work
        <<",\"transpose_work_residual\":"<<e.transpose_work_residual<<",\"jacobian_asymmetry\":"<<e.jacobian_asymmetry;
}
void readTensor(std::istream &stream,Matd &v){for(size_t i=0;i<9;++i)stream>>v(i/3,i%3);}
}
int main(int argc,char **argv)
{
    std::cout<<std::setprecision(17);
    std::ifstream boot("/proc/sys/kernel/random/boot_id"),stat("/proc/self/stat"),cmdline("/proc/self/cmdline",std::ios::binary);
    std::string boot_id,fields,word,ticks;std::getline(boot,boot_id);std::getline(stat,fields);
    std::istringstream tail(fields.substr(fields.rfind(')')+2));for(int i=0;i<=19;++i){tail>>word;if(i==19)ticks=word;}
    std::cout<<"{\"identity\":{\"pid\":"<<getpid()<<",\"ppid\":"<<getppid()<<",\"start_ticks\":";writeJsonString(ticks);
    std::cout<<",\"boot_id\":";writeJsonString(boot_id);std::cout<<",\"argv\":[";bool first=true;
    while(std::getline(cmdline,word,'\0')){if(!first)std::cout<<',';first=false;writeJsonString(word);}std::cout<<"]},\"prototype_only\":true";
    try
    {
        if(argc!=2)throw std::invalid_argument("usage: tetra_constraints_probe INPUT");
        std::ifstream input(argv[1]);int mode;input>>mode;
        if(mode==1)
        {
            size_t n,count;Real dt;PrimalDualSettings settings;
            input>>n>>dt>>settings.velocity_scale>>settings.relative_tolerance>>settings.max_newton_steps>>settings.max_line_search;
            if(!n||n>8)throw std::invalid_argument("quadratic bounded dimensions");
            StdVec<Real> mass(n),q(n),v0(n);DenseContactMatrix k(n,StdVec<Real>(n));
            for(auto *values:{&mass,&q,&v0})for(auto &v:*values)input>>v;
            for(auto &row:k)for(auto &v:row)input>>v;
            input>>count;if(count>4)throw std::invalid_argument("too many constraints");
            StdVec<DenseContactConstraint> rows(count);
            for(auto &row:rows){input>>row.old_offset>>row.new_offset;row.coefficients.resize(n);for(auto &v:row.coefficients)input>>v;}
            if(!input)throw std::invalid_argument("truncated quadratic input");input>>std::ws;if(!input.eof())throw std::invalid_argument("quadratic tail");
            validateDenseStiffness(k);Real total=0;for(Real m:mass){if(!std::isfinite(m)||m<=0)throw std::invalid_argument("quadratic mass");total+=m;}
            auto equation=[&](const StdVec<Real> &v)
            {
                NonlinearContactEquation r; r.momentum.resize(n);r.jacobian.assign(n,StdVec<Real>(n));
                StdVec<Real> average(n);for(size_t i=0;i<n;++i)average[i]=q[i]+dt*v[i]/2;
                for(size_t i=0;i<n;++i){r.momentum[i]=mass[i]*(v[i]-v0[i])+dt*denseContactDot(k[i],average);
                    for(size_t j=0;j<n;++j)r.jacobian[i][j]=(i==j?mass[i]:0)+dt*dt*k[i][j]/2;}
                return r;
            };
            const auto r=solveNonsymmetricContact(q,dt,total,rows,equation,[](const StdVec<Real> &){return true;},settings);
            solver(r);if(r.converged){std::cout<<",\"velocity\":";array(r.velocity);std::cout<<",\"normal_impulse\":";array(r.normal_impulse);}
            std::cout<<",\"model\":\"quadratic_limit\"}\n";return r.converged?0:2;
        }
        if(mode!=0&&mode!=2)throw std::invalid_argument("invalid model");
        TetraGelProblem p;input>>p.dt>>p.volume>>p.material.bulk>>p.material.equilibrium_shear>>p.material.branch_shear>>p.material.relaxation_time
            >>p.settings.velocity_scale>>p.settings.relative_tolerance>>p.settings.max_newton_steps>>p.settings.max_line_search;
        for(auto *values:{&p.node_mass,&p.position,&p.velocity})for(auto &v:*values)input>>v;
        readTensor(input,p.history.green);readTensor(input,p.history.viscous);input>>p.history.dissipation;
        size_t count;input>>count;if(count>4)throw std::invalid_argument("too many constraints");p.contacts.resize(count);
        for(auto &row:p.contacts){input>>row.old_offset>>row.new_offset;row.coefficients.resize(12);for(auto &v:row.coefficients)input>>v;}
        StdVec<Real> sample(12);if(mode==2)for(auto &v:sample)input>>v;
        if(!input)throw std::invalid_argument("truncated tetra input");input>>std::ws;if(!input.eof())throw std::invalid_argument("tetra input tail");
        const auto old=p.history;std::cout<<",\"B\":";dense(tetraReferenceGradient());
        if(mode==2)
        {
            const auto e=evaluateTetraGel(p,sample);evaluation(e);
            std::cout<<",\"velocity\":";array(sample);std::cout<<",\"operator_probe_only\":true,\"converged\":false}\n";return 0;
        }
        const auto r=advanceTetraGel(p);solver(r.solve);
        std::cout<<",\"old_history_unchanged\":"<<((old.green==p.history.green&&old.viscous==p.history.viscous&&old.dissipation==p.history.dissipation)?"true":"false");
        if(r.solve.converged)
        {
            std::cout<<",\"velocity\":";array(r.velocity);std::cout<<",\"position\":";array(r.position);
            std::cout<<",\"normal_impulse\":";array(r.solve.normal_impulse);std::cout<<",\"contact_impulse\":";array(r.contact_impulse);
            std::cout<<",\"momentum_residual\":";array(r.momentum_residual);evaluation(r.evaluation);
            std::cout<<",\"kinetic_old\":"<<r.kinetic_old<<",\"kinetic_new\":"<<r.kinetic_new<<",\"algorithmic_loss\":"<<r.algorithmic_loss
                <<",\"normal_contact_loss\":"<<r.normal_contact_loss<<",\"wall_work\":"<<r.wall_work<<",\"energy_residual\":"<<r.energy_residual
                <<",\"total_momentum_residual\":";vector3(r.total_momentum_residual);
            std::cout<<",\"internal_force_sum\":";vector3(r.internal_force_sum);
            std::cout<<",\"angular_change_minus_external\":";vector3(r.angular_change_minus_external);
            std::cout<<",\"internal_torque_old_cursor\":";vector3(r.internal_torque_old_cursor);
            std::cout<<",\"angular_balance_old_cursor\":";vector3(r.angular_balance_old_cursor);
        }
        std::cout<<"}\n";return r.solve.converged?0:2;
    }
    catch(const std::exception &error){std::cout<<",\"error\":";writeJsonString(error.what());std::cout<<"}\n";return 1;}
}
