#pragma once

#include <cmath>

#if defined(__CUDACC__)
#define SPH_CONTACT_HD __host__ __device__
#else
#define SPH_CONTACT_HD
#endif

namespace SPH::peridynamics
{
template<class Real>
struct ContactVector3
{
    Real x,y,z;

    SPH_CONTACT_HD ContactVector3 operator+(ContactVector3 v)const{return {x+v.x,y+v.y,z+v.z};}
    SPH_CONTACT_HD ContactVector3 operator-(ContactVector3 v)const{return {x-v.x,y-v.y,z-v.z};}
    SPH_CONTACT_HD ContactVector3 operator*(Real a)const{return {a*x,a*y,a*z};}
    SPH_CONTACT_HD Real dot(ContactVector3 v)const{return x*v.x+(y*v.y+z*v.z);}
    SPH_CONTACT_HD ContactVector3 cross(ContactVector3 v)const
    {return {y*v.z-z*v.y,z*v.x-x*v.z,x*v.y-y*v.x};}
};

template<class Real>
struct ObjectiveFrictionState
{
    ContactVector3<Real> slip{},normal{0,0,1};
    Real stored=0,dissipated=0,released=0;
    bool active=false;
};

template<class Real>
struct ObjectiveFrictionResponse
{
    ObjectiveFrictionState<Real> state;
    ContactVector3<Real> force{};
    bool reversed=false;
};

template<class Real>
SPH_CONTACT_HD ObjectiveFrictionResponse<Real> evaluateObjectiveFriction(
    const ObjectiveFrictionState<Real> &old,ContactVector3<Real> n,ContactVector3<Real> relative,
    Real spin,Real normal_force,Real stiffness,Real mu,Real dt,bool commit)
{
    ObjectiveFrictionResponse<Real> response{old};
    if(normal_force<=0||mu==0)
    {
        if(commit)
        {
            response.state.released+=old.stored;
            response.state.stored=0;
            response.state.slip={};
            response.state.active=false;
        }
        return response;
    }
    ContactVector3<Real> trial{};
    if(old.active)
    {
        const Real cosine=old.normal.dot(n);
        const Real c=cosine<-1?-1:cosine>1?1:cosine;
        if(c<Real(-.999999))
        {
            response.reversed=true;
            return response;
        }
        // Shortest rotation transporting the stored slip to this contact plane.
        const auto axis=old.normal.cross(n);
        const auto first=axis.cross(old.slip);
        trial=old.slip+first+axis.cross(first)*(Real(1)/(1+c));
        if(commit&&dt>0)
        {
            const Real theta=spin*dt;
            const auto twist=n.cross(trial);
            trial=trial+twist*std::sin(theta)+n.cross(twist)*(1-std::cos(theta));
        }
        trial=trial-n*trial.dot(n);
        if(commit)
            trial=trial+(relative-n*relative.dot(n))*dt;
    }
    const Real magnitude=stiffness*std::sqrt(trial.dot(trial));
    const Real limit=mu*normal_force;
    Real loss=0;
    if(magnitude>limit)
    {
        const auto before=trial;
        trial=trial*(limit/magnitude);
        const auto plastic=before-trial;
        loss=limit*std::sqrt(plastic.dot(plastic));
    }
    if(commit)
    {
        response.state.slip=trial;
        response.state.normal=n;
        response.state.active=true;
        response.state.stored=Real(.5)*stiffness*trial.dot(trial);
        response.state.dissipated+=loss;
    }
    response.force=trial*(-stiffness);
    return response;
}
}

#undef SPH_CONTACT_HD
