/* Licensed under the Apache License, Version 2.0. */
#ifndef PD_VISCO_FRICTION_STEP_H
#define PD_VISCO_FRICTION_STEP_H

#include "pd_normal_projection.h"

namespace SPH::peridynamics
{
/** Reduced SLS control only: U=K0*q^2/2+Kv*x^2/2, xdot=qdot-x/tau.
 * Not the production finite-strain gel model and not fitted parameters. */
struct ScalarMaxwellParameters
{
    Real equilibrium = 1, branch = 0.5, tau = 0.3;

    void validate() const
    {
        if (!std::isfinite(equilibrium) || equilibrium < 0 || !std::isfinite(branch) || branch < 0 ||
            !std::isfinite(tau) || tau <= 0)
            throw std::invalid_argument("scalar Maxwell: invalid parameters");
    }
};

struct ScalarMaxwellCoefficients
{
    long double decay, alpha, one_minus_alpha, square, cross, ramp_square;
};

/** Integrals over u in [0,1]: square=int exp(-2hu),
 * cross=int exp(-hu)*(1-exp(-hu))/h, ramp_square=int ((1-exp(-hu))/h)^2.
 * Series avoid subtracting nearly equal exponentials for small h. */
inline ScalarMaxwellCoefficients scalarMaxwellCoefficients(Real dt, Real tau)
{
    if (!std::isfinite(dt) || dt <= 0 || !std::isfinite(tau) || tau <= 0)
        throw std::invalid_argument("scalar Maxwell: invalid interval");
    const long double h = static_cast<long double>(dt) / tau;
    ScalarMaxwellCoefficients c{};
    c.decay = std::exp(-h);
    c.alpha = -std::expm1(-h) / h;
    c.square = -std::expm1(-2 * h) / (2 * h);
    if (h < 0.25L)
    {
        long double power = 1, factorial = 1;
        c.one_minus_alpha = 0;
        for (int j = 0; j < 32; ++j)
        {
            if (j) { power *= h; factorial *= j + 1; }
            const long double sign = j % 2 ? -1 : 1;
            // factorial=(j+1)!, valid at j=0 too.
            if (j) c.one_minus_alpha -= sign * power / factorial;
            c.cross += (std::pow(-1.L, j + 1) - std::pow(-2.L, j + 1)) * power /
                       (factorial * (j + 2));
            c.ramp_square += (-2 * std::pow(-1.L, j + 2) + std::pow(-2.L, j + 2)) * power /
                             (factorial * (j + 2) * (j + 3));
        }
    }
    else
    {
        c.one_minus_alpha = 1 - c.alpha;
        c.cross = (c.alpha - c.square) / h;
        c.ramp_square = (1 - 2 * c.alpha + c.square) / (h * h);
    }
    return c;
}

struct ScalarMaxwellInterval
{
    Real extension_new, average_force, stored_old, stored_new, physical_dissipation, energy_residual;
};

inline ScalarMaxwellInterval advanceScalarMaxwell(Real q_old, Real extension_old, Real q_new,
    Real dt, const ScalarMaxwellParameters &p)
{
    p.validate();
    if (!std::isfinite(q_old) || !std::isfinite(extension_old) || !std::isfinite(q_new))
        throw std::invalid_argument("scalar Maxwell: non-finite history");
    const auto c = scalarMaxwellCoefficients(dt, p.tau);
    const long double dq = static_cast<long double>(q_new) - q_old, x = extension_old;
    const long double next = c.decay * x + c.alpha * dq;
    const long double average = p.equilibrium * (static_cast<long double>(q_old) + q_new) / 2 +
        p.branch * (c.alpha * x + c.one_minus_alpha * p.tau * dq / dt);
    const long double integral = c.square * x * x + 2 * c.cross * x * dq + c.ramp_square * dq * dq;
    ScalarMaxwellInterval r{};
    r.extension_new = next;
    r.average_force = average;
    r.stored_old = (p.equilibrium * static_cast<long double>(q_old) * q_old + p.branch * x * x) / 2;
    r.stored_new = (p.equilibrium * static_cast<long double>(q_new) * q_new + p.branch * next * next) / 2;
    r.physical_dissipation = p.branch * static_cast<long double>(dt) / p.tau * integral;
    r.energy_residual = average * dq - (static_cast<long double>(r.stored_new) - r.stored_old) - r.physical_dissipation;
    for (Real value : {r.extension_new, r.average_force, r.stored_old, r.stored_new, r.physical_dissipation, r.energy_residual})
        if (!std::isfinite(value)) throw std::runtime_error("scalar Maxwell: non-finite result");
    if (r.physical_dissipation < 0) throw std::runtime_error("scalar Maxwell: negative physical dissipation");
    return r;
}

struct CoulombImpulseResult
{
    Real velocity, impulse, slip, physical_dissipation, cone_residual, natural_impulse_residual;
    // 0 separated, 1 sticking, 2 positive relative slip, 3 negative relative slip.
    int mode;
};

/** Fixed normal impulse, one tangent. This is the exact box-constrained
 * maximum-dissipation solve in an explicit positive scalar metric. */
inline CoulombImpulseResult solveScalarCoulomb(Real metric, Real free_velocity, Real wall_velocity,
    Real coefficient, Real normal_impulse)
{
    if (!std::isfinite(metric) || metric <= 0 || !std::isfinite(free_velocity) ||
        !std::isfinite(wall_velocity) || !std::isfinite(coefficient) || coefficient < 0 ||
        !std::isfinite(normal_impulse) || normal_impulse < 0)
        throw std::invalid_argument("scalar Coulomb: invalid input");
    const Real bound = coefficient * normal_impulse;
    if (!std::isfinite(bound)) throw std::invalid_argument("scalar Coulomb: non-finite cone");
    CoulombImpulseResult r{};
    const Real trial = metric * (wall_velocity - free_velocity);
    r.impulse = std::clamp(trial, -bound, bound);
    r.velocity = free_velocity + r.impulse / metric;
    r.slip = r.velocity - wall_velocity;
    r.physical_dissipation = -r.impulse * r.slip;
    r.cone_residual = std::max(Real(0), std::abs(r.impulse) - bound);
    r.natural_impulse_residual = std::abs(r.impulse - std::clamp(r.impulse - metric * r.slip, -bound, bound));
    r.mode = normal_impulse == 0 ? 0 : std::abs(trial) <= bound ? 1 : r.slip > 0 ? 2 : 3;
    return r;
}

struct ViscoFrictionState
{
    Real normal_position = 0, normal_velocity = 0, maxwell_extension = 0;
    Real tangent_position = 0, tangent_velocity = 0;
};

struct ViscoFrictionResult
{
    ViscoFrictionState state;
    ScalarMaxwellInterval material;
    NormalMetricProjectionResult normal;
    CoulombImpulseResult tangent;
    Real kinetic_old, kinetic_new, wall_work, body_work, algorithmic_loss;
    Real normal_contact_loss, friction_loss, energy_residual, momentum_normal, momentum_tangent;
};

inline ViscoFrictionResult advanceViscoFriction(const ViscoFrictionState &old, Real mass_normal,
    Real mass_tangent, const ScalarMaxwellParameters &material, Real dt, Real wall_old, Real wall_new,
    Real tangent_wall_velocity, Real coefficient, Real body_normal_force = 0)
{
    material.validate();
    if (!std::isfinite(mass_normal) || mass_normal <= 0 || !std::isfinite(mass_tangent) || mass_tangent <= 0 ||
        !std::isfinite(old.tangent_position) || !std::isfinite(old.tangent_velocity) || !std::isfinite(body_normal_force))
        throw std::invalid_argument("visco friction: invalid physical state");
    const auto c = scalarMaxwellCoefficients(dt, material.tau);
    const long double bias = material.equilibrium * old.normal_position + material.branch * c.alpha * old.maxwell_extension;
    const long double slope = material.equilibrium * dt / 2 + material.branch * material.tau * c.one_minus_alpha;
    const Real metric = mass_normal + dt * slope;
    const Real free_velocity = (mass_normal * old.normal_velocity + dt * (body_normal_force - bias)) / metric;
    ViscoFrictionResult r{};
    r.normal = projectNormalMetric({metric}, {old.normal_position}, {free_velocity}, dt,
        {{{0}, {1}, 1, wall_old, wall_new}});
    if (!r.normal.converged) throw std::runtime_error("visco friction: unconverged normal constraint");
    r.state = old;
    r.state.normal_position = r.normal.position[0];
    r.state.normal_velocity = r.normal.velocity[0];
    r.material = advanceScalarMaxwell(old.normal_position, old.maxwell_extension, r.state.normal_position, dt, material);
    r.state.maxwell_extension = r.material.extension_new;
    r.tangent = solveScalarCoulomb(mass_tangent, old.tangent_velocity, tangent_wall_velocity, coefficient, r.normal.normal_impulse[0]);
    r.state.tangent_velocity = r.tangent.velocity;
    r.state.tangent_position += dt * r.tangent.velocity;
    r.kinetic_old = (mass_normal * old.normal_velocity * old.normal_velocity + mass_tangent * old.tangent_velocity * old.tangent_velocity) / 2;
    r.kinetic_new = (mass_normal * r.state.normal_velocity * r.state.normal_velocity + mass_tangent * r.tangent.velocity * r.tangent.velocity) / 2;
    const Real dvn = r.state.normal_velocity - old.normal_velocity, dvt = r.tangent.velocity - old.tangent_velocity;
    r.algorithmic_loss = (mass_normal * dvn * dvn + mass_tangent * dvt * dvt) / 2;
    r.normal_contact_loss = r.normal.normal_impulse[0] * (old.normal_position - wall_old) / dt;
    r.friction_loss = r.tangent.physical_dissipation;
    r.wall_work = r.normal.normal_impulse[0] * (wall_new - wall_old) / dt + r.tangent.impulse * tangent_wall_velocity;
    r.body_work = body_normal_force * (r.state.normal_position - old.normal_position);
    r.momentum_normal = mass_normal * dvn - dt * (body_normal_force - r.material.average_force) - r.normal.normal_impulse[0];
    r.momentum_tangent = mass_tangent * dvt - r.tangent.impulse;
    r.energy_residual = r.wall_work + r.body_work - (r.kinetic_new - r.kinetic_old + r.material.stored_new - r.material.stored_old) -
        r.material.physical_dissipation - r.algorithmic_loss - r.normal_contact_loss - r.friction_loss;
    for (Real value : {r.kinetic_old, r.kinetic_new, r.algorithmic_loss, r.normal_contact_loss, r.friction_loss,
                       r.wall_work, r.body_work, r.energy_residual, r.momentum_normal, r.momentum_tangent})
        if (!std::isfinite(value)) throw std::runtime_error("visco friction: non-finite physical ledger");
    return r;
}
} // namespace SPH::peridynamics
#endif
