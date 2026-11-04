/* JSONL observation only; no state mutation or force/history sampling kick. */
#ifndef BOUNDED_GEL_PROBE_IO_H
#define BOUNDED_GEL_PROBE_IO_H
#include "pd_bounded_normal_scene.h"
#include <cerrno>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <sys/stat.h>
#include <unistd.h>
namespace bounded_probe_io
{
using namespace SPH;using namespace SPH::peridynamics;
inline void value(const std::string &s)
{
    std::cout<<'"';for(unsigned char c:s)
        if(c=='"'||c=='\\')std::cout<<'\\'<<c;
        else if(c<32)std::cout<<"\\u00"<<std::hex<<std::setw(2)<<std::setfill('0')<<int(c)<<std::dec<<std::setfill(' ');
        else std::cout<<c;
    std::cout<<'"';
}
inline void value(const char *s){value(std::string(s));}
template<class T,std::enable_if_t<std::is_arithmetic_v<T>,int> =0>void value(T x)
{
    if constexpr(std::is_floating_point_v<T>){if(!std::isfinite(x))throw std::runtime_error("nonfinite observer field");}
    if constexpr(std::is_same_v<T,bool>)std::cout<<(x?"true":"false");
    else if constexpr(std::is_integral_v<T> && sizeof(T)==1)std::cout<<int(x);
    else std::cout<<x;
}
template<class Derived>void value(const Eigen::MatrixBase<Derived> &m)
{
    std::cout<<'[';bool first=true;for(Eigen::Index i=0;i<m.rows();++i)for(Eigen::Index j=0;j<m.cols();++j)
    {if(!first)std::cout<<',';first=false;value(m(i,j));}std::cout<<']';
}
template<class T>void value(const StdVec<T> &v);
template<class T,size_t N>void value(const std::array<T,N> &v);
template<class T>void value(const StdVec<T> &v)
{std::cout<<'[';for(size_t i=0;i<v.size();++i){if(i)std::cout<<',';value(v[i]);}std::cout<<']';}
template<class T,size_t N>void value(const std::array<T,N> &v)
{std::cout<<'[';for(size_t i=0;i<N;++i){if(i)std::cout<<',';value(v[i]);}std::cout<<']';}
template<class T>void field(const char *key,const T &x){std::cout<<',';value(key);std::cout<<':';value(x);}
inline void begin(const char *type){std::cout<<"{\"type\":";value(type);}
inline void finish()
{
    std::cout<<"}\n";std::cout.flush();if(!std::cout)throw std::runtime_error("observer output failure");
    struct stat s;if(::fstat(STDOUT_FILENO,&s)!=0)throw std::runtime_error("observer stdout identity");
    if(S_ISREG(s.st_mode) && ::fsync(STDOUT_FILENO)!=0)throw std::runtime_error("observer durable output failure");
}
inline void identity()
{
    std::ifstream boot("/proc/sys/kernel/random/boot_id"),stat("/proc/self/stat"),cmd("/proc/self/cmdline",std::ios::binary);
    std::string id,line,word,ticks;std::getline(boot,id);std::getline(stat,line);std::istringstream fields(line.substr(line.rfind(')')+2));
    for(int k=0;k<=19;++k){fields>>word;if(k==19)ticks=word;}StdVec<std::string> argv;while(std::getline(cmd,word,'\0'))argv.push_back(word);
    begin("identity");field("pid",getpid());field("ppid",getppid());field("start_ticks",ticks);field("boot_id",id);field("argv",argv);finish();
}
inline void reference(const FixedPDReference &r)
{
    std::cout<<'{';std::cout<<"\"generation\":";value(r.generation);field("horizon",r.horizon);field("position",r.position);field("volume",r.volume);
    field("offsets",r.offsets);field("target",r.target);field("weight",r.weight);field("active",r.active);field("inverse_shape",r.inverse_shape);
    field("shape_gradient",r.shape_gradient);field("neighbor_measure",r.neighbor_measure);field("shape_eigenvalue_ratio",r.shape_eigenvalue_ratio);
    field("nonaffine_unit",r.nonaffine_unit);std::cout<<'}';
}
inline void surface(const ReferenceSurfaceTransfer &s)
{
    std::cout<<"{\"generation\":";value(s.generation);field("particle_count",s.particle_count);field("support",s.reference_support);field("reference",s.reference);
    field("areas",s.areas);field("offsets",s.offsets);field("targets",s.targets);field("weights",s.weights);std::cout<<'}';
}
inline void state(const BoundedGelState &s)
{
    std::cout<<"{\"time\":";value(s.time);field("step",s.steps);field("position",s.position);field("velocity",s.velocity);field("green",s.green);
    field("viscous",s.viscous);field("dissipation",s.dissipation);field("accepted_stress",s.accepted_stress);field("density",s.density);
    field("last_dt",s.last_dt);field("last_KKT",s.last_KKT);field("last_normal_impulse",s.last_normal_impulse);field("last_newton_iterations",s.last_newton_iterations);std::cout<<'}';
}
inline void solver(const ManyContactResult &r)
{
    std::cout<<"{\"converged\":";value(r.converged);field("status",r.status);field("velocity",r.velocity);field("normal_impulse",r.normal_impulse);
    field("normalized_KKT",r.normalized_KKT);std::cout<<",\"trace\":[";
    for(size_t k=0;k<r.trace.size();++k)
    {
        if(k)std::cout<<',';const auto &t=r.trace[k];std::cout<<"{\"iteration\":";value(t.iteration);field("linear_rank",t.linear_rank);
        field("domain_rejections",t.domain_rejections);field("merit",t.merit);field("squared_norm",t.squared_norm);field("KKT",t.KKT);
        field("alpha",t.accepted_step);field("linear_residual",t.linear_residual);field("velocity",t.velocity);field("normal_impulse",t.normal_impulse);
        field("direction",t.normalized_direction);field("normalized_residual",t.normalized_residual);field("normalized_gaps",t.normalized_gaps);
        field("accumulation_roundoff_bounds",t.accumulation_roundoff_bounds);field("generalized_jacobian",t.generalized_jacobian);std::cout<<'}';
    }
    std::cout<<"],\"precision_rejection\":{\"present\":";value(r.precision_rejection.present);
    if(r.precision_rejection.present)
    {
        const auto &p=r.precision_rejection;field("during_line_search",p.during_line_search);field("iteration",p.iteration);field("alpha",p.trial_alpha);
        field("velocity",p.velocity);field("normal_impulse",p.normal_impulse);field("supplied_momentum",p.supplied_momentum);
        field("normalized_values",p.normalized_values);field("roundoff_bounds",p.roundoff_bounds);
    }
    std::cout<<"}}";
}
inline void evaluation(const BoundedGelEvaluation &e,const BoundedGelModel &m,const BoundedGelState &old)
{
    std::cout<<"{\"momentum\":";value(e.equation.momentum);field("jacobian",e.equation.jacobian);field("material_force",e.material_force);
    field("stabilization_force",e.stabilization_force);field("internal_force",e.internal_force);field("jacobian_asymmetry",e.jacobian_asymmetry);
    field("material_stored_old",e.material_stored_old);field("material_stored_new",e.material_stored_new);field("material_work",e.material_work);
    field("physical_material_loss",e.physical_material_loss);field("stabilization_stored_old",e.stabilization_stored_old);field("stabilization_stored_new",e.stabilization_stored_new);
    field("transpose_work_residual",e.transpose_work_residual);field("stabilization_work_residual",e.stabilization_work_residual);std::cout<<",\"points\":[";
    for(size_t i=0;i<e.point.size();++i)
    {
        if(i)std::cout<<',';const auto &p=e.point[i];std::cout<<"{\"F0\":";value(m.reference.deformation(i,old.position));field("F1",e.final_deformation[i]);
        field("average_pk1",p.average_pk1);field("endpoint_tangent",p.endpoint_tangent);field("green",p.history.green);field("viscous",p.history.viscous);
        field("total_dissipation",p.history.dissipation);field("stored_old",p.stored_old);field("stored_new",p.stored_new);field("physical_dissipation",p.physical_dissipation);
        field("work",p.work);field("energy_residual",p.energy_residual);field("minimum_determinant",p.minimum_determinant);std::cout<<'}';
    }
    std::cout<<"]}";
}
} // namespace bounded_probe_io
#endif
