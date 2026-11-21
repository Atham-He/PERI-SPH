#pragma once

#include "pd_plastic_multiplier.h"

namespace SPH::peridynamics
{
template <class Scalar>
struct DruckerPragerHardening
{
    Scalar linear, saturation, rate;

    EIGEN_DEVICE_FUNC Scalar operator()(Scalar alpha) const
    {
        using std::exp;
        const Scalar exponent = exp(-rate * alpha);
        return linear * alpha + saturation * (1.0 - exponent);
    }
};

template <class Scalar>
struct DruckerPragerResidual
{
    Scalar pressure_term, q, resistance, yield, alpha;
    DruckerPragerHardening<Scalar> hardening;

    EIGEN_DEVICE_FUNC Scalar operator()(Scalar gamma) const
    {
        return pressure_term + q - resistance * gamma - yield - hardening(alpha + gamma);
    }
};

/** Associated Drucker-Prager return with a smooth cone and a hydrostatic
 * apex. Stress and plastic strain are the corotated elastic trial state.
 */
template <class Matrix, class Scalar>
EIGEN_DEVICE_FUNC PlasticMultiplierStatus returnDruckerPrager(
    Matrix &stress, Matrix &plastic_strain, Scalar &alpha,
    Scalar shear, Scalar bulk, Scalar pressure_sensitivity, Scalar yield,
    const DruckerPragerHardening<Scalar> &hardening, Scalar yield_tolerance)
{
    using std::abs;
    using std::sqrt;
    const Scalar root_two = sqrt(2.0);
    const Scalar inverse_dimension = 1.0 / Matrix::RowsAtCompileTime;
    const Scalar trace = stress.trace();
    const Matrix dev = stress - trace * inverse_dimension * Matrix::Identity();
    const Scalar dev_norm = dev.norm(), pressure_term = pressure_sensitivity * trace;
    const Scalar trial_value = pressure_term + dev_norm / root_two - (hardening(alpha) + yield);
    if (!(trial_value > yield_tolerance))
        return PlasticMultiplierStatus::Accepted;

    const Scalar q_trial = dev_norm / root_two;
    const Scalar gamma_apex = q_trial / shear;
    const Scalar volume_resistance = 9.0 * bulk * pressure_sensitivity * pressure_sensitivity;
    const DruckerPragerResidual<Scalar> smooth{pressure_term, q_trial, volume_resistance + shear, yield, alpha, hardening};
    Scalar gamma = 0;
    if (smooth(gamma_apex) <= 0.0)
    {
        const Scalar scale = abs(pressure_term) + q_trial;
        const auto status = bracketPlasticMultiplier(smooth, Scalar(0), gamma_apex, scale < yield ? yield : scale, gamma);
        if (status != PlasticMultiplierStatus::Accepted)
            return status;
        const Matrix direction = dev / dev_norm;
        stress -= gamma * (3.0 * bulk * pressure_sensitivity * Matrix::Identity() + root_two * shear * direction);
        plastic_strain += gamma * (pressure_sensitivity * Matrix::Identity() + direction / root_two);
    }
    else
    {
        const DruckerPragerResidual<Scalar> apex{pressure_term, Scalar(0), volume_resistance, yield, alpha, hardening};
        const Scalar upper = gamma_apex + trial_value / volume_resistance;
        const Scalar scale = abs(pressure_term);
        const auto status = bracketPlasticMultiplier(apex, gamma_apex, upper, scale < yield ? yield : scale, gamma);
        if (status != PlasticMultiplierStatus::Accepted)
            return status;
        stress = (trace * inverse_dimension - 3.0 * bulk * pressure_sensitivity * gamma) * Matrix::Identity();
        plastic_strain += dev / (2.0 * shear) + pressure_sensitivity * gamma * Matrix::Identity();
    }

    alpha += gamma;
    return PlasticMultiplierStatus::Accepted;
}
}
