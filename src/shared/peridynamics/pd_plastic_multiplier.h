#pragma once

#include "pd_plastic_multiplier_status.h"

#include <Eigen/Core>
#include <cmath>
#include <stdexcept>

namespace SPH::peridynamics
{
template <class Scalar, class Residual>
EIGEN_DEVICE_FUNC PlasticMultiplierStatus bracketPlasticMultiplier(
    const Residual &residual, Scalar lower, Scalar upper, Scalar scale, Scalar &multiplier)
{
    using std::abs;
    using std::isfinite;
    for (int i = 0; residual(upper) > 0.0 && i < 64; ++i)
        upper = 2.0 * upper + 1e-16;
    if (!isfinite(residual(upper)) || residual(upper) > 0.0)
        return PlasticMultiplierStatus::InvalidBracket;

    for (int i = 0; i < 100; ++i)
    {
        const Scalar mid = 0.5 * (lower + upper);
        const Scalar value = residual(mid);
        if (!isfinite(value))
            return PlasticMultiplierStatus::NonFiniteResidual;
        if (abs(value) <= 1e-11 * (scale < Scalar(1) ? Scalar(1) : scale))
        {
            multiplier = mid;
            return PlasticMultiplierStatus::Accepted;
        }
        if (value > 0.0)
            lower = mid;
        else
            upper = mid;
    }
    return PlasticMultiplierStatus::IterationLimit;
}
[[noreturn]] inline void throwPlasticMultiplierFailure(PlasticMultiplierStatus status)
{
    switch (status)
    {
    case PlasticMultiplierStatus::InvalidBracket:
        throw std::runtime_error("PD constitutive return mapping: failed to bracket plastic multiplier");
    case PlasticMultiplierStatus::NonFiniteResidual:
        throw std::runtime_error("PD constitutive return mapping: non-finite residual");
    case PlasticMultiplierStatus::IterationLimit:
        throw std::runtime_error("PD constitutive return mapping: local solve did not converge");
    case PlasticMultiplierStatus::Accepted:
        break;
    }
    throw std::logic_error("Unknown PD plastic multiplier status");
}

}
