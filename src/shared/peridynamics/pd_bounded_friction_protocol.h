/* Licensed under the Apache License, Version 2.0. */
#ifndef PD_BOUNDED_FRICTION_PROTOCOL_H
#define PD_BOUNDED_FRICTION_PROTOCOL_H

#include "pd_bounded_gel_friction.h"

namespace SPH::peridynamics
{
struct CubicPlateOffset
{
    std::array<Real,4> coefficients{-1,0,0,0};
    Real at(Real s) const
    { return coefficients[0]+s*(coefficients[1]+s*(coefficients[2]+s*coefficients[3])); }
};

inline Real minimumUnitCubic(const std::array<Real,4> &c)
{
    for(Real x:c)if(!std::isfinite(x))throw std::invalid_argument("friction protocol: nonfinite cubic");
    auto at=[&](Real s){return c[0]+s*(c[1]+s*(c[2]+s*c[3]));};
    Real value=std::min(at(0),at(1));
    auto consider=[&](Real s){if(s>0 && s<1)value=std::min(value,at(s));};
    const Real a=3*c[3],b=2*c[2],d=c[1];
    if(a==0){if(b!=0)consider(-d/b);}
    else
    {
        const Real discriminant=b*b-4*a*d;
        if(discriminant>=0)
        {
            const Real q=-.5*(b+std::copysign(std::sqrt(discriminant),b));
            if(q!=0){consider(q/a);consider(d/q);}else consider(-b/(2*a));
        }
    }
    return value;
}

struct BoundedFrictionInterval
{
    size_t phase=0,step=0;
    Real old_time=0,new_time=0,actual_dt=0;
    std::array<Real,2> old_offset{-1,-1},new_offset{-1,-1};
    std::array<CubicPlateOffset,2> normal_cubic;
    std::array<ContactTangentVector,2> old_tangent_displacement{
        ContactTangentVector::Zero(),ContactTangentVector::Zero()};
    std::array<ContactTangentVector,2> new_tangent_displacement{
        ContactTangentVector::Zero(),ContactTangentVector::Zero()};
    std::array<ContactTangentVector,2> wall_tangent_velocity{
        ContactTangentVector::Zero(),ContactTangentVector::Zero()};
};

/** Fixed16-step press, optionally8-step hold and16-step release. The static
 * scene's reference/MLS never changes. Tangent wall velocity is explicitly
 * the interval-mean prescribed displacement difference/actual dt.
 */
struct BoundedFrictionProtocol
{
    int kind=0; // 0 press16, 1 press16/hold8/release16
    Real nominal_dt=.01,normal_stroke=.01;
    ContactTangentVector tangent_stroke=ContactTangentVector::Zero();

    void validate() const
    {
        if(kind<0 || kind>1 || !std::isfinite(nominal_dt) || nominal_dt<=0 ||
            !std::isfinite(normal_stroke) || normal_stroke<=0 || normal_stroke>=.1 ||
            !tangent_stroke.allFinite() || tangent_stroke.norm()>=.1)
            throw std::invalid_argument("friction protocol: configuration");
        if(!std::isfinite(40*nominal_dt))throw std::invalid_argument("friction protocol: time range");
    }
    size_t steps() const { return kind==0 ? 16 : 40; }
    Real endTime() const { return Real(steps())*nominal_dt; }
    std::array<Real,3> phaseEnds() const { return {16*nominal_dt,24*nominal_dt,40*nominal_dt}; }

    BoundedFrictionInterval interval(size_t step,Real accepted_time) const
    {
        validate();
        if(step>=steps() || !std::isfinite(accepted_time) || accepted_time<0)
            throw std::invalid_argument("friction protocol: accepted cursor");
        const size_t phase=step<16 ? 0 : (step<24 ? 1 : 2);
        const auto ends=phaseEnds();
        const Real begin=phase==0 ? 0 : ends[phase-1],end=ends[phase];
        const bool last=step==15 || step==23 || step==39;
        const Real dt=last ? end-accepted_time : nominal_dt;
        const Real next=accepted_time+dt;
        if(accepted_time<begin || accepted_time>=end || !std::isfinite(dt) || dt<=0 || next>end)
            throw std::invalid_argument("friction protocol: phase clock mismatch");
        const Real duration=end-begin;
        const Real theta=(accepted_time-begin)/duration,delta=dt/duration;
        auto shape=[&](Real time)
        {
            if(phase==1)return Real(1);
            const Real s=(time-begin)/duration,smooth=s*s*(3-2*s);
            return phase==0 ? smooth : 1-smooth;
        };
        const Real old_shape=shape(accepted_time),new_shape=shape(next);
        BoundedFrictionInterval result;
        result.phase=phase;result.step=step;result.old_time=accepted_time;result.new_time=next;result.actual_dt=dt;
        result.old_offset[1]=-1+normal_stroke*old_shape;result.new_offset[1]=-1+normal_stroke*new_shape;
        result.old_tangent_displacement[1]=tangent_stroke*old_shape;
        result.new_tangent_displacement[1]=tangent_stroke*new_shape;
        result.wall_tangent_velocity[1]=(result.new_tangent_displacement[1]-result.old_tangent_displacement[1])/dt;
        auto &c=result.normal_cubic[1].coefficients;c[0]=result.old_offset[1];
        if(phase!=1)
        {
            const Real amplitude=phase==0 ? normal_stroke : -normal_stroke;
            c[1]=amplitude*(6*theta-6*theta*theta)*delta;
            c[2]=amplitude*(3-6*theta)*delta*delta;
            c[3]=-2*amplitude*delta*delta*delta;
        }
        return result;
    }
};

inline std::array<Real,2> boundedFrictionIntervalGaps(const BoundedNormalScene &scene,
    const StdVec<Vecd> &q0,const StdVec<Vecd> &q1,const BoundedFrictionInterval &interval,bool whole_interval=true)
{
    if(q0.size()!=scene.model.reference.size() || q1.size()!=q0.size())
        throw std::invalid_argument("friction protocol: observation dimensions");
    std::array<Real,2> minimum{std::numeric_limits<Real>::infinity(),std::numeric_limits<Real>::infinity()};
    for(const auto *map:{&scene.surface,&scene.dense_surface})for(size_t s=0;s<map->reference.size();++s)
    {
        const Vecd old=map->gather(s,q0.data()),next=map->gather(s,q1.data());
        for(size_t plate=0;plate<2;++plate)
        {
            const auto &d=interval.normal_cubic[plate].coefficients;
            const std::array<Real,4> gap{scene.normal[plate].dot(old)-d[0],
                scene.normal[plate].dot(next-old)-d[1],-d[2],-d[3]};
            minimum[plate]=std::min(minimum[plate],whole_interval ? minimumUnitCubic(gap) :
                scene.normal[plate].dot(next)-interval.new_offset[plate]);
        }
    }
    return minimum;
}
} // namespace SPH::peridynamics
#endif
