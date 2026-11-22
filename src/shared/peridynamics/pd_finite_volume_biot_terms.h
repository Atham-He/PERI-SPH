#pragma once

#include <Eigen/Core>
#include <cmath>

namespace SPH::peridynamics
{
template <class Real>
struct FiniteVolumeBiotVolumetric
{
    Real energy, piola_scale;
};

template <class Real>
EIGEN_DEVICE_FUNC inline FiniteVolumeBiotVolumetric<Real> finiteVolumeBiotVolumetric(Real bulk, Real jacobian)
{
    const Real z=std::log(jacobian);
    const Real j2m1=std::expm1(2*z);
    return {Real(.25)*bulk*(j2m1-2*z),Real(.5)*bulk*j2m1};
}

template <class Real, class Stretches>
EIGEN_DEVICE_FUNC inline Real finiteVolumeBiotTangent(
    Real bulk,Real shear,const Stretches &stretch,Real J,Real damage)
{
    const Real minimum=stretch.minCoeff(),active=1-damage;
    const Real volume=bulk*J*J*stretch.cwiseInverse().squaredNorm()+
        Real(.5)*bulk*std::abs(J*J-1)/(minimum*minimum);
    const Real distortional=2*shear*(stretch.maxCoeff()+
        (stretch.array()-stretch.mean()).abs().maxCoeff())/minimum;
    return (J<1?1:active)*volume+active*distortional;
}
}
