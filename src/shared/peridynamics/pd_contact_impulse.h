/* Licensed under the Apache License, Version 2.0. */
#ifndef PD_CONTACT_IMPULSE_H
#define PD_CONTACT_IMPULSE_H

#include "base_data_type.h"
#include <cmath>
#include <stdexcept>

namespace SPH::peridynamics
{
/** Accepted-step integral of the actual midpoint reactions on each plate.
 * This is an observation ledger, not a force filter or terminal-force estimate.
 * Compensation terms and the committed cursor are checkpointed together. */
class ContactImpulseLedger
{
  public:
    void initialize(const Vecd &momentum)
    {
        if (initialized_ || steps_ != 0 || time_ != 0 || !momentum.allFinite())
            throw std::invalid_argument("contact impulse: initial state must be finite and fresh");
        values_[4] = momentum;
        initialized_ = 1;
    }

    void advance(Real accepted_time, Real dt, Real sample_time, const Vecd &lower_reaction, const Vecd &upper_reaction)
    {
        if (initialized_ != 1 || !std::isfinite(accepted_time) || !std::isfinite(dt) || !std::isfinite(sample_time) ||
            dt <= 0 || accepted_time <= time_ || !lower_reaction.allFinite() || !upper_reaction.allFinite() ||
            std::abs(accepted_time - time_ - dt) > 1e-20 + 1e-12 * dt ||
            std::abs(sample_time - (time_ + 0.5 * dt)) > 1e-20 + 1e-12 * dt)
            throw std::invalid_argument("contact impulse: invalid accepted-step/sample cursor or reaction");
        add(0, dt * lower_reaction);
        add(1, dt * upper_reaction);
        time_ = accepted_time;
        ++steps_;
    }

    const Vecd &lower() const { return values_[0]; }
    const Vecd &upper() const { return values_[1]; }
    Vecd momentumResidual(const Vecd &accepted_momentum) const
    {
        if (initialized_ != 1 || !accepted_momentum.allFinite())
            throw std::invalid_argument("contact impulse: invalid momentum observation");
        return accepted_momentum - values_[4] + lower() + upper();
    }
    bool matchesCursor(Real time, size_t steps) const { return initialized_ == 1 && time_ == time && steps_ == steps; }

    template<class Checkpoint> void registerCheckpoint(Checkpoint &checkpoint)
    {
        checkpoint.addArray("ContactImpulseSumsCompensationInitialMomentum", values_);
        checkpoint.addScalar("ContactImpulseAcceptedTime", time_);
        checkpoint.addScalar("ContactImpulseAcceptedSteps", steps_);
        checkpoint.addScalar("ContactImpulseInitialized", initialized_);
    }

  private:
    void add(size_t side, const Vecd &increment)
    {
        const Vecd corrected = increment - values_[side + 2];
        const Vecd total = values_[side] + corrected;
        values_[side + 2] = (total - values_[side]) - corrected;
        values_[side] = total;
    }
    StdVec<Vecd> values_ = StdVec<Vecd>(5, Vecd::Zero());
    Real time_ = 0;
    size_t steps_ = 0;
    uint64_t initialized_ = 0;
};
} // namespace SPH::peridynamics
#endif
