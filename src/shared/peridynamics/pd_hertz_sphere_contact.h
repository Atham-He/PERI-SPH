#pragma once

#include "pd_objective_contact_friction.h"

#include <cmath>
#include <cstddef>
#include <memory>
#include <stdexcept>
#include <vector>

#if defined(__CUDACC__)
#define SPH_HERTZ_HD __host__ __device__
#else
#define SPH_HERTZ_HD
#endif

namespace SPH::peridynamics
{
enum class HertzPairStatus { Separated, Active, Coincident };

template <class Real>
struct HertzSphereSite
{
    Real center[3], radius;
};

struct HertzSpherePair
{
    std::size_t first, second;
};

template <class Real>
struct HertzSphereContact
{
    Real normal[3], first_arm[3], second_arm[3];
    Real penetration, force, stiffness;
    HertzPairStatus status;
};

template <class Real>
SPH_HERTZ_HD Real hertzNormalForce(Real penetration,Real effective_radius,Real effective_modulus)
{
    return penetration>0 ? Real(4)/3*effective_modulus*std::sqrt(effective_radius)*
        std::pow(penetration,1.5) : 0;
}

template <class Real>
SPH_HERTZ_HD HertzSphereContact<Real> evaluateHertzSpherePair(
    const HertzSphereSite<Real> &a,const HertzSphereSite<Real> &b,Real effective_modulus)
{
    HertzSphereContact<Real> result{};
    const Real d[3]={a.center[0]-b.center[0],a.center[1]-b.center[1],a.center[2]-b.center[2]};
    const Real distance=std::sqrt(d[0]*d[0]+(d[1]*d[1]+d[2]*d[2]));
    result.penetration=a.radius+b.radius-distance;
    if (result.penetration<=0)
        return result;
    if (distance<=0)
    {
        result.status=HertzPairStatus::Coincident;
        return result;
    }
    const Real first_distance=(distance*distance+a.radius*a.radius-b.radius*b.radius)/(2*distance);
    for (int axis=0;axis<3;++axis)
    {
        const Real normal=d[axis]/distance;
        const Real contact=a.center[axis]-first_distance*normal;
        result.normal[axis]=normal;
        result.first_arm[axis]=contact-a.center[axis];
        result.second_arm[axis]=contact-b.center[axis];
    }
    const Real effective_radius=a.radius*b.radius/(a.radius+b.radius);
    result.force=hertzNormalForce(result.penetration,effective_radius,effective_modulus);
    result.stiffness=1.5*result.force/result.penetration;
    result.status=HertzPairStatus::Active;
    return result;
}

template<class Real>
struct FrictionalSphereSite
{
    HertzSphereSite<Real> geometry;
    ContactVector3<Real> velocity,omega;
};

template<class Real>
struct FrictionalSphereParameters
{
    Real effective_modulus,stiffness,friction,dt;
    bool commit;
};

template<class Real>
struct FrictionalSphereCandidate
{
    HertzSpherePair indices;
    ObjectiveFrictionState<Real> old;
};

template<class Real>
struct FrictionalSphereResult
{
    HertzSphereContact<Real> contact;
    ObjectiveFrictionResponse<Real> friction;
};

template<class Real>
SPH_HERTZ_HD ObjectiveFrictionResponse<Real> evaluateSphereFriction(
    const FrictionalSphereSite<Real> &a,const FrictionalSphereSite<Real> &b,
    const HertzSphereContact<Real> &contact,const ObjectiveFrictionState<Real> &old,
    FrictionalSphereParameters<Real> parameters)
{
    if(contact.status==HertzPairStatus::Coincident)
        return {old};
    if(contact.status==HertzPairStatus::Separated)
        return evaluateObjectiveFriction(old,{},{},Real(0),Real(0),parameters.stiffness,
                                         parameters.friction,parameters.dt,parameters.commit);
    const ContactVector3<Real> n{contact.normal[0],contact.normal[1],contact.normal[2]};
    const ContactVector3<Real> ra{contact.first_arm[0],contact.first_arm[1],contact.first_arm[2]};
    const ContactVector3<Real> rb{contact.second_arm[0],contact.second_arm[1],contact.second_arm[2]};
    const auto relative=a.velocity+a.omega.cross(ra)-b.velocity-b.omega.cross(rb);
    const Real spin=Real(.5)*(a.omega+b.omega).dot(n);
    return evaluateObjectiveFriction(old,n,relative,spin,contact.force,parameters.stiffness,
                                     parameters.friction,parameters.dt,parameters.commit);
}

template<class Real>
SPH_HERTZ_HD FrictionalSphereResult<Real> evaluateFrictionalSpherePair(
    const FrictionalSphereSite<Real> &a,const FrictionalSphereSite<Real> &b,
    const ObjectiveFrictionState<Real> &old,FrictionalSphereParameters<Real> parameters)
{
    const auto contact=evaluateHertzSpherePair(a.geometry,b.geometry,parameters.effective_modulus);
    return {contact,evaluateSphereFriction(a,b,contact,old,parameters)};
}

class HertzSpherePairBatch
{
  public:
    virtual ~HertzSpherePairBatch() = default;
    virtual void evaluate(const std::vector<HertzSphereSite<double>> &particles,
                          const std::vector<HertzSpherePair> &candidates,double effective_modulus,
                          std::vector<HertzSphereContact<double>> &result) = 0;
    virtual bool supportsFriction() const noexcept { return false; }
    virtual void evaluateFrictional(const std::vector<FrictionalSphereSite<double>> &sites,
                                    const std::vector<FrictionalSphereCandidate<double>> &candidates,
                                    FrictionalSphereParameters<double> parameters,
                                    std::vector<FrictionalSphereResult<double>> &result)
    {
        std::vector<HertzSphereSite<double>> geometry;
        std::vector<HertzSpherePair> indices;
        geometry.reserve(sites.size());
        indices.reserve(candidates.size());
        for(const auto &site:sites) geometry.push_back(site.geometry);
        for(const auto &candidate:candidates) indices.push_back(candidate.indices);
        std::vector<HertzSphereContact<double>> contact;
        evaluate(geometry,indices,parameters.effective_modulus,contact);
        if(contact.size()!=indices.size())
            throw std::logic_error("Sphere contact backend returned an incomplete batch");
        result.resize(contact.size());
        for(std::size_t i=0;i<contact.size();++i)
        {
            const auto pair=indices[i];
            result[i]={contact[i],evaluateSphereFriction(sites.at(pair.first),sites.at(pair.second),
                                                        contact[i],candidates[i].old,parameters)};
        }
    }
};

std::unique_ptr<HertzSpherePairBatch> makeCudaHertzSpherePairs();
}

#undef SPH_HERTZ_HD
