#pragma once

#include <cmath>
#include <cstddef>
#include <memory>
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

class HertzSpherePairBatch
{
  public:
    virtual ~HertzSpherePairBatch() = default;
    virtual void evaluate(const std::vector<HertzSphereSite<double>> &particles,
                          const std::vector<HertzSpherePair> &candidates,double effective_modulus,
                          std::vector<HertzSphereContact<double>> &result) = 0;
};

std::unique_ptr<HertzSpherePairBatch> makeCudaHertzSpherePairs();
}

#undef SPH_HERTZ_HD
