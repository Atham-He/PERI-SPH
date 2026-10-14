#ifndef PD_NUMERICS_H
#define PD_NUMERICS_H

#include "adaptation.h"
#include "base_particles.hpp"
#include "constraint_dynamics.h"
#include "sph_system.h"
#include <algorithm>
#include <cmath>
#include <initializer_list>
#include <limits>
#include <sstream>
#include <stdexcept>

namespace SPH
{
namespace peridynamics
{
/** One accepted increment is used by every solver and the physical clock. */
inline Real boundedTimeStep(Real candidate, std::initializer_list<Real> remaining)
{
    if (!std::isfinite(candidate) || candidate <= 0.0)
        throw std::runtime_error("PERI-SPH: non-positive or non-finite proposed time step");
    for (Real limit : remaining)
    {
        if (!std::isfinite(limit) || limit <= 0.0)
            throw std::runtime_error("PERI-SPH: non-positive or non-finite remaining interval");
        candidate = std::min(candidate, limit);
    }
    return candidate;
}

[[noreturn]] inline void numericalFailure(const std::string &body, size_t particle,
                                         const std::string &quantity)
{
    std::ostringstream message;
    message << "PERI-SPH numerical failure: body=" << body << " particle=" << particle
            << " quantity=" << quantity;
    throw std::runtime_error(message.str());
}

/** Checks are explicit and configurable; no invalid particle is silently removed. */
class ParticleStateCheck : public LocalDynamics
{
  public:
    explicit ParticleStateCheck(SPHBody &body) : LocalDynamics(body)
    {
        addVector("Position");
        addVector("Velocity");
    }
    void addVector(const std::string &name)
    {
        vectors_.push_back({name, particles_->getVariableDataByName<Vecd>(name)});
    }
    void addMatrix(const std::string &name)
    {
        matrices_.push_back({name, particles_->getVariableDataByName<Matd>(name)});
    }
    void addScalar(const std::string &name, Real minimum = -std::numeric_limits<Real>::infinity())
    {
        scalars_.push_back({name, particles_->getVariableDataByName<Real>(name), minimum});
    }
    void update(size_t i, Real = 0.0)
    {
        for (const auto &v : vectors_)
            if (!v.second[i].allFinite()) numericalFailure(sph_body_->Name(), i, v.first);
        for (const auto &m : matrices_)
            if (!m.second[i].allFinite()) numericalFailure(sph_body_->Name(), i, m.first);
        for (const auto &s : scalars_)
            if (!std::isfinite(s.data[i]) || s.data[i] <= s.minimum)
                numericalFailure(sph_body_->Name(), i, s.name);
    }

  protected:
    struct Scalar { std::string name; Real *data; Real minimum; };
    StdVec<std::pair<std::string, Vecd *>> vectors_;
    StdVec<std::pair<std::string, Matd *>> matrices_;
    StdVec<Scalar> scalars_;
};
} // namespace peridynamics

namespace solid_dynamics
{
/** Conservative PD increment using wave speed, particle speed and acceleration. */
class PDTimeStep : public LocalDynamicsReduce<ReduceMax<Real>>
{
  public:
    PDTimeStep(SPHBody &body, Real cfl = 0.1)
        : LocalDynamicsReduce<ReduceMax<Real>>(body), cfl_(cfl),
          spacing_(body.getSPHAdaptation().ReferenceSpacing()),
          sound_speed_(DynamicCast<ElasticSolid>(this, body.getMatterMaterial()).ReferenceSoundSpeed()),
          velocity_(particles_->registerStateVariableData<Vecd>("Velocity")),
          acceleration_(particles_->registerStateVariableData<Vecd>("PDAcceleration"))
    {
        if (!std::isfinite(cfl) || cfl <= 0 || cfl > 0.2)
            throw std::invalid_argument("PDTimeStep requires 0 < CFL <= 0.2");
    }
    Real reduce(size_t i, Real = 0.0)
    {
        if (!velocity_[i].allFinite() || !acceleration_[i].allFinite())
            peridynamics::numericalFailure(sph_body_->Name(), i, "PD time-step state");
        return std::max((sound_speed_ + velocity_[i].norm()) / spacing_,
                        std::sqrt(acceleration_[i].norm() / spacing_));
    }
    Real outputResult(Real inverse_step) override { return cfl_ / inverse_step; }
  protected:
    Real cfl_, spacing_, sound_speed_;
    Vecd *velocity_, *acceleration_;
};

/** Exact translation from the reference position; velocity agrees with displacement.
 * Constant speed on [start,end), zero outside. Acceleration is zero within each
 * segment; callers needing a smooth start should prescribe a smooth motion law. */
class PrescribedTranslation : public MotionConstraint<SPHBody>
{
  public:
    PrescribedTranslation(SPHBody &body, Real start, Real end, const Vecd &translation)
        : MotionConstraint<SPHBody>(body), start_(start), end_(end), translation_(translation),
          time_(sph_system_->svPhysicalTime().Data()),
          average_velocity_(particles_->registerStateVariableData<Vecd>("AverageVelocity")),
          acceleration_(particles_->registerStateVariableData<Vecd>("Acceleration")),
          average_acceleration_(particles_->registerStateVariableData<Vecd>("AverageAcceleration"))
    {
        if (!std::isfinite(start) || !std::isfinite(end) || end <= start || !translation.allFinite())
            throw std::invalid_argument("PrescribedTranslation: invalid motion interval or displacement");
    }
    void update(size_t i, Real = 0.0)
    {
        Real fraction = std::clamp((*time_ - start_) / (end_ - start_), Real(0), Real(1));
        pos_[i] = pos0_[i] + fraction * translation_;
        vel_[i] = (*time_ >= start_ && *time_ < end_) ? Vecd(translation_ / (end_ - start_)) : Vecd::Zero();
        average_velocity_[i] = vel_[i];
        acceleration_[i].setZero();
        average_acceleration_[i].setZero();
    }

  protected:
    Real start_, end_;
    Vecd translation_;
    Real *time_;
    Vecd *average_velocity_, *acceleration_, *average_acceleration_;
};
} // namespace solid_dynamics
} // namespace SPH
#endif
