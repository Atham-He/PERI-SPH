/* Licensed under the Apache License, Version 2.0. */
#ifndef PD_BOUNDED_NORMAL_SCENE_H
#define PD_BOUNDED_NORMAL_SCENE_H
#include "pd_bounded_gel.h"
#include "pd_surface_quadrature.h"

namespace SPH::peridynamics
{
struct BoundedNormalConfig
{
    int mode=1,rotation=0,invalid_injection=0;
    size_t steps=16,azimuths=16;
    Real dt=.01,stroke=.01,rho0=1,eta=.1,support=1.25,cap_width=1./64,relaxed_fraction=0;
    GelParameters material{5,1,.5,.3};
    Matd initial_affine=Matd::Identity(),velocity_rate=Matd::Zero();
    Vecd local_perturbation=Vecd::Zero(),local_velocity_perturbation=Vecd::Zero(),translation=Vecd::Zero();
    size_t perturb_node=16;
    void validate()const
    {
        material.validate();
        if(mode<0 || mode>2 || rotation<0 || rotation>1 || invalid_injection<0 || invalid_injection>4 || steps<1 || steps>64 ||
            azimuths<8 || azimuths>16 || !std::isfinite(dt) || dt<=0 || !std::isfinite(stroke) || stroke<=0 || stroke>=.1 ||
            !std::isfinite(rho0) || rho0<=0 || !std::isfinite(eta) || eta<0 || !std::isfinite(support) || support<=0 ||
            !std::isfinite(cap_width) || cap_width<=0 || cap_width>=.5 || !std::isfinite(relaxed_fraction) || relaxed_fraction<0 || relaxed_fraction>1 ||
            !initial_affine.allFinite() || !velocity_rate.allFinite() || !local_perturbation.allFinite() || !local_velocity_perturbation.allFinite() ||
            !translation.allFinite() || perturb_node>=33 || !std::isfinite(steps*dt))throw std::invalid_argument("bounded normal scene: configuration");
        if(mode!=0 && ((initial_affine-Matd::Identity()).norm()!=0 || velocity_rate.norm()!=0 || local_perturbation.norm()!=0 ||
            local_velocity_perturbation.norm()!=0 || relaxed_fraction!=0))throw std::invalid_argument("bounded normal scene: ignored initial-state input");
        if(mode!=2 && translation.norm()!=0)throw std::invalid_argument("bounded normal scene: ignored translation input");
    }
    Real endTime()const{return Real(steps)*dt;}
    Real upperDisplacement(Real time)const
    {
        if(!std::isfinite(time) || time<0 || time>endTime())throw std::invalid_argument("bounded normal scene: motion cursor");
        const Real s=time/endTime();return stroke*s*s*(3-2*s);
    }
};

/** Full actual sphere/contact/dense observation references are frozen at
 * independent t0. The dense map is observation only, never force/control. */
struct BoundedNormalScene
{
    BoundedNormalConfig config;
    BoundedGelModel model;
    ReferenceSurfaceTransfer surface,dense_surface;
    Matd frame=Matd::Identity();
    std::array<Vecd,2> normal{Vecd(0,1,0),Vecd(0,-1,0)};

    explicit BoundedNormalScene(const BoundedNormalConfig &c):config(c)
    {
        config.validate();if(config.rotation)frame<<0,-1,0,1,0,0,0,0,1;
        auto reference=boundedSphereReference();for(auto &x:reference)x=frame*x;
        const Real volume=4*std::acos(Real(-1))/(3*reference.size());
        model=BoundedGelModel(FixedPDReference(reference,StdVec<Real>(reference.size(),volume),1.25),config.material,config.rho0,config.eta);
        PolePanelSphereRule rule(Vecd::Zero(),1,config.azimuths,config.cap_width);
        for(auto &x:rule.position)x=frame*x;
        surface=ReferenceSurfaceTransfer(reference.data(),reference.size(),rule,config.support);
        SphereSurfaceRule dense(Vecd::Zero(),1,4,16,4);
        for(auto &x:dense.position)x=frame*x;
        dense_surface=ReferenceSurfaceTransfer(reference.data(),reference.size(),dense,config.support);
        for(auto &n:normal)n=frame*n;
        surface.validateStored(reference.data(),reference.size());dense_surface.validateStored(reference.data(),reference.size());
        if(config.invalid_injection==1)model.node_mass[0]*=2;
        if(config.invalid_injection==3)model.reference.active[0]=0;
        if(config.invalid_injection==4)
        {
            for(auto &w:surface.weights)if(w<0)w=0;
            surface.validateStored(reference.data(),reference.size());
        }
    }
    BoundedGelState initialState()const
    {
        auto state=model.initialState();
        if(config.mode==0)
        {
            for(size_t i=0;i<state.position.size();++i)
            {
                const Vecd local=frame.transpose()*model.reference.position[i];
                Vecd current=config.initial_affine*local;if(i==config.perturb_node)current+=config.local_perturbation;
                state.position[i]=frame*current;
            }
            for(size_t i=0;i<state.position.size();++i)
            {
                const Matd f=model.reference.deformation(i,state.position);state.green[i]=gelDeviatoricGreen(f);
                state.viscous[i]=config.relaxed_fraction*state.green[i];const auto response=gelResponse(f,model.material,state.history(i));
                state.accepted_stress[i]=response.cauchy;state.density[i]=config.rho0/f.determinant();
            }
        }
        if(config.mode==2)for(auto &v:state.velocity)v=frame*config.translation;
        if(config.invalid_injection==2)state.green[0](0,0)+=.01;
        return state;
    }
    StdVec<Real> sampleVelocity()const
    {
        StdVec<Vecd> v(model.reference.size());
        for(size_t i=0;i<v.size();++i)
        {
            const Vecd local=frame.transpose()*model.reference.position[i];
            Vecd current=config.velocity_rate*local;if(i==config.perturb_node)current+=config.local_velocity_perturbation;
            v[i]=frame*current;
        }
        return flattenPD(v);
    }
    StdVec<DenseContactConstraint> contacts(Real old_time,Real new_time)const
    {
        if(config.mode==2)return {};
        StdVec<DenseContactConstraint> rows;
        const Real old_u=config.upperDisplacement(old_time),new_u=config.upperDisplacement(new_time);
        for(size_t s=0;s<surface.reference.size();++s)for(size_t plate=0;plate<2;++plate)
        {
            DenseContactConstraint row;row.coefficients.assign(3*model.reference.size(),0);
            row.old_offset=-1+(plate?old_u:0);row.new_offset=-1+(plate?new_u:0);
            for(size_t p=surface.offsets[s];p<surface.offsets[s+1];++p)for(size_t a=0;a<3;++a)
                row.coefficients[3*surface.targets[p]+a]=normal[plate][a]*surface.weights[p];
            rows.push_back(std::move(row));
        }
        return rows;
    }
    StdVec<size_t> plateIDs()const
    {
        if(config.mode==2)return {};StdVec<size_t> ids;
        for(size_t s=0;s<surface.reference.size();++s){ids.push_back(0);ids.push_back(1);}return ids;
    }
    std::array<Real,2> denseGaps(const StdVec<Vecd> &q,Real time)const
    {
        std::array<Real,2> minimum{std::numeric_limits<Real>::infinity(),std::numeric_limits<Real>::infinity()};
        const Real u=config.upperDisplacement(time);
        // The union includes the exact two poles in the control map. Areas
        // are not summed over this union or used as new physical forces.
        for(const auto *map:{&surface,&dense_surface})for(size_t s=0;s<map->reference.size();++s)
        {
            const Vecd x=map->gather(s,q.data());
            for(size_t plate=0;plate<2;++plate)minimum[plate]=std::min(minimum[plate],normal[plate].dot(x)+1-(plate?u:0));
        }
        return minimum;
    }
    std::array<Real,2> denseIntervalGaps(const StdVec<Vecd> &q0,const StdVec<Vecd> &q1,Real t0,Real t1)const
    {
        if(q0.size()!=model.reference.size() || q1.size()!=q0.size() || t1<=t0)throw std::invalid_argument("bounded normal interval geometry");
        const Real theta=t0/config.endTime(),h=(t1-t0)/config.endTime();
        std::array<Real,2> minimum{std::numeric_limits<Real>::infinity(),std::numeric_limits<Real>::infinity()};
        for(const auto *map:{&surface,&dense_surface})for(size_t s=0;s<map->reference.size();++s)
        {
            const Vecd x0=map->gather(s,q0.data()),x1=map->gather(s,q1.data());
            for(size_t plate=0;plate<2;++plate)
            {
                const Real c0=normal[plate].dot(x0)+1-(plate?config.upperDisplacement(t0):0);
                const Real c1=normal[plate].dot(x1-x0)-(plate?config.stroke*(6*theta-6*theta*theta)*h:0);
                const Real c2=plate?-config.stroke*(3-6*theta)*h*h:0,c3=plate?2*config.stroke*h*h*h:0;
                auto at=[&](Real u){return c0+u*(c1+u*(c2+u*c3));};
                Real value=std::min(at(0),at(1));
                auto consider=[&](Real u){if(u>0 && u<1)value=std::min(value,at(u));};
                const Real a=3*c3,b=2*c2,c=c1;
                if(a==0){if(b!=0)consider(-c/b);}
                else
                {
                    const Real discriminant=b*b-4*a*c;
                    if(discriminant>=0)
                    {
                        const Real q=-.5*(b+std::copysign(std::sqrt(discriminant),b));
                        if(q!=0){consider(q/a);consider(c/q);}else consider(-b/(2*a));
                    }
                }
                minimum[plate]=std::min(minimum[plate],value);
            }
        }
        return minimum;
    }
};
} // namespace SPH::peridynamics
#endif
