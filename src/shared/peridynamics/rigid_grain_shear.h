#pragma once
/** Rigid-sphere contact and a split direct-shear box. SI units.
 * Hertz normal elasticity + objective constant-stiffness Coulomb history.
 * This is NOT a Mindlin variable-stiffness implementation. Capsule contact is
 * frictionless surface-penalty contact with the exact transpose force map.
 * No strength, force scale, erosion, or prescribed capsule deformation is used.
 */
#include "pd_solid_reference_coupling.h"
#include <Eigen/Geometry>
#include <fstream>
#include <sstream>
#include <map>
#include <unordered_map>
#include <array>

namespace SPH::peridynamics
{
struct ShearBoxOptions
{
    Real length_x_m=0,length_y_m=0,initial_height_m=0,normal_pressure_Pa=0,shear_end_m=0;
    size_t shear_increments=0,normal_increments=0;
    Real grain_E_Pa=0,grain_nu=0,grain_mu=0,platen_mu=0,side_mu=0;
    Real grain_tangent_stiffness_N_m=0,capsule_grain_mu=0;
    Real normal_ramp_s=0,shear_ramp_s=0,maximum_stage_duration_s=0;
    size_t output_interval_steps=200,checkpoint_wall_seconds=300;
    void read(const std::string&path)
    {
        std::map<std::string,Real*> real{{"length_x_m",&length_x_m},{"length_y_m",&length_y_m},{"initial_height_m",&initial_height_m},
          {"normal_pressure_Pa",&normal_pressure_Pa},{"shear_end_m",&shear_end_m},{"grain_E_Pa",&grain_E_Pa},{"grain_nu",&grain_nu},
          {"grain_mu",&grain_mu},{"platen_mu",&platen_mu},{"side_mu",&side_mu},{"grain_tangent_stiffness_N_m",&grain_tangent_stiffness_N_m},
          {"capsule_grain_mu",&capsule_grain_mu},{"normal_ramp_s",&normal_ramp_s},{"shear_ramp_s",&shear_ramp_s},{"maximum_stage_duration_s",&maximum_stage_duration_s}};
        std::map<std::string,size_t*> integer{{"shear_increments",&shear_increments},{"normal_increments",&normal_increments},
          {"output_interval_steps",&output_interval_steps},{"checkpoint_wall_seconds",&checkpoint_wall_seconds}};
        std::ifstream file(path);if(!file)throw std::runtime_error("Missing shear box input");std::string line;
        while(std::getline(file,line)){line=line.substr(0,line.find('#'));if(line.empty())continue;auto e=line.find('=');if(e==std::string::npos)throw std::runtime_error("Invalid box key=value");
          auto key=line.substr(0,e),value=line.substr(e+1);if(real.count(key))*real[key]=std::stod(value);else if(integer.count(key))*integer[key]=std::stoull(value);else throw std::runtime_error("Unknown shear box input: "+key);}
        if(std::min({length_x_m,length_y_m,initial_height_m,normal_pressure_Pa,shear_end_m,grain_E_Pa,grain_tangent_stiffness_N_m,normal_ramp_s,shear_ramp_s})<=0||
           !shear_increments||!normal_increments||grain_nu<=-1||grain_nu>=.5||grain_mu<0||platen_mu<0||side_mu<0||capsule_grain_mu!=0)
          throw std::runtime_error("Invalid box input or unqualified capsule friction model");
    }
};

struct ObjectiveLinearFriction
{
    Vecd slip=Vecd::Zero(),normal=Vecd::UnitZ();bool active=false;
    Real stored=0,dissipated=0,released=0;
    void clear(){released+=stored;stored=0;slip.setZero();active=false;}
    Vecd force(const Vecd&n,const Vecd&relative_velocity,Real spin,Real fn,Real kt,Real mu,Real dt,bool commit)
    {
        if(fn<=0||mu==0){if(commit)clear();return Vecd::Zero();}
        Vecd trial=Vecd::Zero();
        if(active){Real c=std::clamp(normal.dot(n),Real(-1),Real(1));if(c<-.999999)throw std::runtime_error("Contact frame reversed in one step");
          Vecd axis=normal.cross(n);Matd skew;skew<<0,-axis[2],axis[1],axis[2],0,-axis[0],-axis[1],axis[0],0;
          trial=(Matd::Identity()+skew+skew*skew/(1+c))*slip;
          if(commit&&dt>0)trial=Eigen::AngleAxis<Real>(spin*dt,n)*trial;
          trial-=n*trial.dot(n);
          if(commit)trial+=dt*(relative_velocity-n*relative_velocity.dot(n));}
        Real magnitude=kt*trial.norm(),limit=mu*fn,loss=0;
        if(magnitude>limit){Vecd before=trial;trial*=limit/magnitude;loss=limit*(before-trial).norm();}
        if(commit){slip=trial;normal=n;active=true;stored=.5*kt*trial.squaredNorm();dissipated+=loss;}
        return -kt*trial;
    }
};

struct RigidGrain
{
    Vecd reference,x,v=Vecd::Zero(),omega=Vecd::Zero(),force=Vecd::Zero(),torque=Vecd::Zero();
    Real radius,mass,inertia;
};

/** Current-configuration spatial search. Each grain occupies one bin.
 * width=max(radius); one-cell reach covers every point/sphere overlap, and
 * two-cell reach every sphere/sphere overlap. IDs are sorted to preserve the
 * original accumulation order. No cutoff, radius or contact law is changed. */
class CurrentSphereBins
{
    using Cell=std::array<long long,3>;
    struct Hash{size_t operator()(const Cell&a)const{size_t h=0;for(auto v:a)h^=std::hash<long long>{}(v)+0x9e3779b9+(h<<6)+(h>>2);return h;}};
    std::unordered_map<Cell,StdVec<size_t>,Hash>bins,one,two;
    Real width=0;
    Cell cell(const Vecd&x)const
    {if(!x.allFinite())throw std::runtime_error("Nonfinite spatial-search point");return Cell{(long long)std::floor(x[0]/width),(long long)std::floor(x[1]/width),(long long)std::floor(x[2]/width)};}
  public:
    explicit CurrentSphereBins(const StdVec<RigidGrain>&grains)
    {
        for(const auto&g:grains)width=std::max(width,g.radius);
        if(width<=0)throw std::runtime_error("Invalid sphere-search width");
        for(size_t i=0;i<grains.size();++i)bins[cell(grains[i].x)].push_back(i);
    }
    const StdVec<size_t>&candidates(const Vecd&x,int reach)
    {
        auto key=cell(x);auto&cache=reach==1?one:two;auto found=cache.find(key);if(found!=cache.end())return found->second;
        StdVec<size_t>ids;
        for(int a=-reach;a<=reach;++a)for(int b=-reach;b<=reach;++b)for(int c=-reach;c<=reach;++c){auto q=bins.find(Cell{key[0]+a,key[1]+b,key[2]+c});if(q!=bins.end())ids.insert(ids.end(),q->second.begin(),q->second.end());}
        std::sort(ids.begin(),ids.end());return cache.emplace(key,std::move(ids)).first->second;
    }
};

class RigidGrainShearBox
{
  public:
    ShearBoxOptions options;
    ReferenceSolidSurface&capsule;
    ReferenceSolidSurface*core_surface=nullptr;
    StdVec<Real>areas,core_areas;
    StdVec<RigidGrain> grains;
    StdVec<ObjectiveLinearFriction>pairs,walls;
    StdVec<Vecd>trace_position,core_position;
    const StdVec<Vecd>*cached_core_positions=nullptr;
    Real maximum_trace_z=std::numeric_limits<Real>::infinity(),maximum_core_z=std::numeric_limits<Real>::infinity();
    bool active_pairs_cache_valid=false;
    StdVec<size_t>active_pair_indices;
    size_t pair_candidate_checks=0,surface_candidate_checks=0;
    void cacheCorePositions(const StdVec<Vecd>&x){cached_core_positions=&x;}
    const StdVec<Vecd>&corePositions()const{return cached_core_positions?*cached_core_positions:core_position;}
    void invalidateNeighborCache(){active_pairs_cache_valid=false;active_pair_indices.clear();}
    Real Ksurface,top=0,normal_load=0,elastic_energy=0,tangent_energy=0,friction_dissipation=0;
    Real max_sphere_penetration=0,max_capsule_penetration=0,max_overlap_ratio=0,stable_dt=1;
    Real force_mapping_residual=0,torque_mapping_residual=0,normal_control_residual=0;
    size_t grain_contacts=0,capsule_contacts=0,wall_contacts=0;
    Vecd lower_reaction=Vecd::Zero(),upper_reaction=Vecd::Zero(),boundary_force=Vecd::Zero();
    Vecd boundary_torque=Vecd::Zero();
    RigidGrainShearBox(ReferenceSolidSurface&surface,const StdVec<Real>&measure,Real penalty,
                      const std::string&spheres,const ShearBoxOptions&input,
                      ReferenceSolidSurface*inner=nullptr,const StdVec<Real>&inner_areas={})
      :options(input),capsule(surface),core_surface(inner),areas(measure),core_areas(inner_areas),Ksurface(penalty)
    {
        if(penalty<=0||areas.size()!=capsule.rows.size())throw std::runtime_error("Invalid sphere/trace contact input");
        std::ifstream f(spheres);if(!f)throw std::runtime_error("Missing sand geometry");std::string line;std::getline(f,line);
        while(std::getline(f,line)){if(line.empty())continue;std::replace(line.begin(),line.end(),',',' ');std::istringstream s(line);RigidGrain g;
          if(!(s>>g.x[0]>>g.x[1]>>g.x[2]>>g.radius>>g.mass>>g.inertia)||std::min({g.radius,g.mass,g.inertia})<=0)throw std::runtime_error("Invalid sand row");g.reference=g.x;grains.push_back(g);}
        if(grains.empty())throw std::runtime_error("Empty sand specimen");pairs.resize(grains.size()*grains.size());walls.resize(6*grains.size());trace_position.resize(areas.size());top=input.initial_height_m/2;
        if(core_surface&&core_surface->rows.size()!=core_areas.size())throw std::runtime_error("Invalid core contact quadrature");core_position.resize(core_areas.size());
    }
    Real hertz(Real penetration,Real radius,bool wall=false)const
    {Real effective=options.grain_E_Pa/((wall?1:2)*(1-options.grain_nu*options.grain_nu));return penetration>0?Real(4)/3*effective*std::sqrt(radius)*std::pow(penetration,1.5):0;}
    Real pressureAt(Real height)const
    {
        Real f=0;for(const auto&g:grains)f+=hertz(std::max(g.x[2]+g.radius-height,Real(0)),g.radius,true);
        if(height<maximum_trace_z)for(size_t q=0;q<areas.size();++q)f+=Ksurface*areas[q]*std::max(trace_position[q][2]-height,Real(0));
        if(height<maximum_core_z)for(size_t q=0;q<core_areas.size();++q)f+=Ksurface*core_areas[q]*std::max(corePositions()[q][2]-height,Real(0));return f;
    }
    void evaluate(Real shear,Real shear_velocity,Real target_load,Real dt,bool commit)
    {
        const size_t n=grains.size(),np=capsule.port.count,ncore=core_surface?core_surface->port.count:0,base=np+ncore;StdVec<Real>rate(base+2*n,0);
        elastic_energy=tangent_energy=friction_dissipation=0;grain_contacts=capsule_contacts=wall_contacts=0;
        max_sphere_penetration=max_capsule_penetration=max_overlap_ratio=0;
        lower_reaction.setZero();upper_reaction.setZero();boundary_force.setZero();boundary_torque.setZero();
        for(auto&g:grains){g.force.setZero();g.torque.setZero();}
        for(size_t q=0;q<areas.size();++q)trace_position[q]=capsule.position(q);
        if(!cached_core_positions)for(size_t q=0;q<core_areas.size();++q)core_position[q]=core_surface->position(q);
        if(corePositions().size()!=core_areas.size())throw std::runtime_error("Cached core trace dimensions disagree");
        maximum_trace_z=maximum_core_z=-std::numeric_limits<Real>::infinity();
        for(const auto&x:trace_position)maximum_trace_z=std::max(maximum_trace_z,x[2]);
        for(const auto&x:corePositions())maximum_core_z=std::max(maximum_core_z,x[2]);
        CurrentSphereBins spatial(grains);pair_candidate_checks=surface_candidate_checks=0;
        if(!active_pairs_cache_valid){active_pair_indices.clear();for(size_t k=0;k<pairs.size();++k)if(pairs[k].active)active_pair_indices.push_back(k);active_pairs_cache_valid=true;}
        StdVec<size_t>next_active_pairs;
        // Exact quasi-static force-controlled top platen. It carries no invented
        // mass and does not overwrite particle forces, velocities or positions.
        Real hi=-1e10;for(const auto&g:grains)hi=std::max(hi,g.x[2]+g.radius);hi=std::max({hi,maximum_trace_z,maximum_core_z});
        Real lo=-options.initial_height_m/2;
        if(target_load<0||pressureAt(lo)<target_load)throw std::runtime_error("Normal piston load cannot be bracketed");
        if(target_load==0)top=hi;else{for(int k=0;k<60;++k){Real mid=.5*(lo+hi);if(pressureAt(mid)>target_load)lo=mid;else hi=mid;}top=.5*(lo+hi);}
        normal_load=pressureAt(top);normal_control_residual=std::abs(normal_load-target_load)/std::max(target_load,Real(1e-20));
        auto grainRate=[&](size_t i,Real lever,Real stiffness,Real sum){rate[base+2*i]+=stiffness/std::sqrt(grains[i].mass)*sum;rate[base+2*i+1]+=stiffness*lever/std::sqrt(grains[i].inertia)*sum;};
        for(size_t i=0;i<n;++i)for(size_t j:spatial.candidates(grains[i].x,2))
        {
            if(j<=i)continue;++pair_candidate_checks;
            auto&a=grains[i];auto&b=grains[j];auto&history=pairs[i*n+j];Vecd d=a.x-b.x;Real length=d.norm(),pen=a.radius+b.radius-length;
            if(pen<=0){if(commit)history.clear();continue;}if(length<=0)throw std::runtime_error("Coincident grains");
            Vecd normal=d/length;Real li=(length*length+a.radius*a.radius-b.radius*b.radius)/(2*length);
            Vecd point=a.x-li*normal,ra=point-a.x,rb=point-b.x;
            Real reff=a.radius*b.radius/(a.radius+b.radius),fn=hertz(pen,reff),kn=1.5*fn/pen;
            Vecd relative=a.v+a.omega.cross(ra)-b.v-b.omega.cross(rb);
            Vecd tangential=history.force(normal,relative,.5*(a.omega+b.omega).dot(normal),fn,options.grain_tangent_stiffness_N_m,options.grain_mu,dt,commit),f=fn*normal+tangential;
            if(commit&&history.active)next_active_pairs.push_back(i*n+j);
            a.force+=f;b.force-=f;a.torque+=ra.cross(f);b.torque-=rb.cross(f);elastic_energy+=.4*fn*pen;
            ++grain_contacts;max_sphere_penetration=std::max(max_sphere_penetration,pen);max_overlap_ratio=std::max(max_overlap_ratio,pen/std::min(a.radius,b.radius));
            Real stiffness=kn+2*options.grain_tangent_stiffness_N_m;
            Real sum=1/std::sqrt(a.mass)+ra.norm()/std::sqrt(a.inertia)+1/std::sqrt(b.mass)+rb.norm()/std::sqrt(b.inertia);
            grainRate(i,ra.norm(),stiffness,sum);grainRate(j,rb.norm(),stiffness,sum);
        }
        if(commit){for(size_t k:active_pair_indices)if(!std::binary_search(next_active_pairs.begin(),next_active_pairs.end(),k))pairs[k].clear();active_pair_indices.swap(next_active_pairs);}
        auto wall=[&](size_t i,size_t id,const Vecd&normal,Real pen,const Vecd&velocity,bool lower,Real mu){
            auto&g=grains[i];auto&history=walls[6*i+id];if(pen<=0){if(commit)history.clear();return;}
            Real fn=hertz(pen,g.radius,true),kn=1.5*fn/pen;Vecd lever=-(g.radius-pen)*normal;
            Vecd tangential=history.force(normal,g.v+g.omega.cross(lever)-velocity,.5*g.omega.dot(normal),fn,options.grain_tangent_stiffness_N_m,mu,dt,commit),f=fn*normal+tangential;
            g.force+=f;g.torque+=lever.cross(f);(lower?lower_reaction:upper_reaction)-=f;boundary_force+=f;boundary_torque+=(g.x+lever).cross(f);elastic_energy+=.4*fn*pen;++wall_contacts;
            max_sphere_penetration=std::max(max_sphere_penetration,pen);max_overlap_ratio=std::max(max_overlap_ratio,pen/g.radius);
            Real stiffness=kn+2*(mu>0?options.grain_tangent_stiffness_N_m:0),sum=1/std::sqrt(g.mass)+lever.norm()/std::sqrt(g.inertia);grainRate(i,lever.norm(),stiffness,sum);
        };
        for(size_t i=0;i<n;++i)
        {
            const auto&g=grains[i];bool lower=g.x[2]<0;Real offset=lower?shear:0;Vecd moving=lower?Vecd(shear_velocity,0,0):Vecd::Zero();
            wall(i,0,Vecd::UnitX(),-options.length_x_m/2+offset-(g.x[0]-g.radius),moving,lower,options.side_mu);
            wall(i,1,-Vecd::UnitX(),g.x[0]+g.radius-(options.length_x_m/2+offset),moving,lower,options.side_mu);
            wall(i,2,Vecd::UnitY(),-options.length_y_m/2-(g.x[1]-g.radius),moving,lower,options.side_mu);
            wall(i,3,-Vecd::UnitY(),g.x[1]+g.radius-options.length_y_m/2,moving,lower,options.side_mu);
            wall(i,4,Vecd::UnitZ(),-options.initial_height_m/2-(g.x[2]-g.radius),Vecd(shear_velocity,0,0),true,options.platen_mu);
            wall(i,5,-Vecd::UnitZ(),g.x[2]+g.radius-top,Vecd::Zero(),false,options.platen_mu);
        }
        auto surfaceForces=[&](ReferenceSolidSurface&surface,const StdVec<Real>&measure,const StdVec<Vecd>&positions,size_t offset){
        auto traceRate=[&](size_t q,Real stiffness,size_t grain,bool with_grain){
            const auto&row=surface.rows[q];Real sum=0;for(size_t k=0;k<row.index.size();++k)sum+=std::abs(row.weight[k])/std::sqrt(surface.port.mass[row.index[k]]);
            if(with_grain)sum+=1/std::sqrt(grains[grain].mass);
            for(size_t k=0;k<row.index.size();++k)rate[offset+row.index[k]]+=stiffness*std::abs(row.weight[k])/std::sqrt(surface.port.mass[row.index[k]])*sum;
            if(with_grain)rate[base+2*grain]+=stiffness/std::sqrt(grains[grain].mass)*sum;
        };
        for(size_t q=0;q<measure.size();++q)
        {
            Vecd x=positions[q];Real stiffness=Ksurface*measure[q];
            for(size_t i:spatial.candidates(x,1)){++surface_candidate_checks;auto&g=grains[i];Vecd delta=x-g.x;Real r2=delta.squaredNorm();if(r2>=g.radius*g.radius)continue;Real distance=std::sqrt(r2);if(distance<=0)throw std::runtime_error("Surface at grain center");Real pen=g.radius-distance;Vecd f=stiffness*pen/distance*delta;
              surface.scatter(q,f);g.force-=f;g.torque+=(x-g.x).cross(-f);elastic_energy+=.5*stiffness*pen*pen;++capsule_contacts;max_capsule_penetration=std::max(max_capsule_penetration,pen);max_overlap_ratio=std::max(max_overlap_ratio,pen/g.radius);traceRate(q,stiffness*(1+pen/distance),i,true);}
            auto plane=[&](const Vecd&normal,Real pen,bool lower){if(pen<=0)return;Vecd f=stiffness*pen*normal;surface.scatter(q,f);(lower?lower_reaction:upper_reaction)-=f;boundary_force+=f;boundary_torque+=x.cross(f);elastic_energy+=.5*stiffness*pen*pen;max_capsule_penetration=std::max(max_capsule_penetration,pen);traceRate(q,stiffness,0,false);};
            bool lower=x[2]<0;Real offset=lower?shear:0;
            plane(Vecd::UnitX(),-options.length_x_m/2+offset-x[0],lower);plane(-Vecd::UnitX(),x[0]-(options.length_x_m/2+offset),lower);
            plane(Vecd::UnitY(),-options.length_y_m/2-x[1],lower);plane(-Vecd::UnitY(),x[1]-options.length_y_m/2,lower);
            plane(Vecd::UnitZ(),-options.initial_height_m/2-x[2],true);plane(-Vecd::UnitZ(),x[2]-top,false);
        }};
        surfaceForces(capsule,areas,trace_position,0);
        // Apply physical nonpenetration to BOTH nested material boundaries.
        // The core contact is geometric, never activated by an arbitrary damage
        // threshold; it prevents sand passing through exposed core material.
        if(core_surface)surfaceForces(*core_surface,core_areas,corePositions(),np);
        Real largest=0;for(Real value:rate)largest=std::max(largest,value);stable_dt=largest>0?.2/std::sqrt(largest):1;
        for(const auto&h:pairs){tangent_energy+=h.stored;friction_dissipation+=h.dissipated+h.released;}
        for(const auto&h:walls){tangent_energy+=h.stored;friction_dissipation+=h.dissipated+h.released;}
        // The capsule port is cleared by the caller AFTER the bonded interface
        // is collected. Thus this audit contains only external granular forces.
        Vecd force=-boundary_force,torque=-boundary_torque;Real scale=boundary_force.norm();
        for(size_t i=0;i<np;++i){Vecd f=capsule.port.storage==SolidForcePort::Storage::ForcePrior?capsule.port.external[i]:capsule.port.mass[i]*capsule.port.external[i];force+=f;torque+=capsule.port.position[i].cross(f);scale+=f.norm();}
        if(core_surface)for(size_t i=0;i<ncore;++i){auto&p=core_surface->port;Vecd f=p.storage==SolidForcePort::Storage::ForcePrior?p.external[i]:p.mass[i]*p.external[i];force+=f;torque+=p.position[i].cross(f);scale+=f.norm();}
        for(const auto&g:grains){force+=g.force;torque+=g.x.cross(g.force)+g.torque;scale+=g.force.norm();}
        force_mapping_residual=force.norm()/std::max(scale,Real(1e-30));torque_mapping_residual=torque.norm()/std::max(scale*options.length_x_m,Real(1e-30));
        if(!std::isfinite(elastic_energy)||std::max(force_mapping_residual,torque_mapping_residual)>1e-9)throw std::runtime_error("Granular coupling force/torque conservation failed");
        if(max_overlap_ratio>.05)throw std::runtime_error("Contact penetration exceeds5% grain radius; unphysical trajectory rejected");
    }
};
}
