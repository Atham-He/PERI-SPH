/* Licensed under the Apache License, Version 2.0. */
#ifndef PD_HARDENING_H
#define PD_HARDENING_H

#include "base_data_type.h"
#include <cmath>
#include <initializer_list>
#include <stdexcept>

namespace SPH
{
/** Isotropic contribution q(alpha) to the initial J2 yield stress [Pa].
 * alpha is accumulated equivalent plastic strain. Perfect() disables only
 * isotropic hardening: strict perfect plasticity also requires H_kin=eta=0.
 * Parameters are immutable; changing a law requires a new material instance.
 */
class PDIsotropicHardening
{
  public:
    enum class Law { Perfect, Linear, SaturationLinear, Legacy };

    static PDIsotropicHardening Perfect() { return {Law::Perfect, 0, 0, 0}; }
    static PDIsotropicHardening Linear(Real modulus) { return {Law::Linear, modulus, 0, 0}; }
    static PDIsotropicHardening Saturation(Real modulus, Real increment, Real rate)
    { return {Law::SaturationLinear, modulus, increment, rate}; }
    /** Preserve the original constructor's constants and floating-point formula. */
    static PDIsotropicHardening Legacy(Real yield_stress, Real modulus = 0)
    {
        if (!std::isfinite(yield_stress) || yield_stress <= 0)
            throw std::invalid_argument("PD hardening: positive finite initial yield stress required");
        return {Law::Legacy, modulus, 1.589 * yield_stress - yield_stress, 16.93};
    }

    Law law() const { return law_; }
    Real linearModulus() const { return linear_; }
    Real saturationIncrement() const { return increment_; }
    Real saturationRate() const { return rate_; }

    Real value(Real alpha) const
    {
        if (!std::isfinite(alpha) || alpha < 0)
            throw std::invalid_argument("PD hardening: nonnegative finite equivalent plastic strain required");
        const Real saturation = law_ == Law::Legacy ? 1.0 - std::exp(-rate_ * alpha)
                                                     : -std::expm1(-rate_ * alpha);
        return linear_ * alpha + increment_ * saturation;
    }

  private:
    PDIsotropicHardening(Law law, Real linear, Real increment, Real rate)
        : law_(law), linear_(linear), increment_(increment), rate_(rate)
    {
        for (Real parameter : {linear, increment, rate})
            if (!std::isfinite(parameter) || parameter < 0)
                throw std::invalid_argument("PD hardening: finite nonnegative parameters required");
    }
    const Law law_;
    const Real linear_, increment_, rate_;
};
} // namespace SPH
#endif
