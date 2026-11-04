/* Licensed under the Apache License, Version 2.0. */
#ifndef PD_TETRA_GEL_CONSTRAINTS_H
#define PD_TETRA_GEL_CONSTRAINTS_H

#include "pd_gel_affine_path.h"
#include "pd_nonsymmetric_constraints.h"

namespace SPH::peridynamics
{
inline DenseContactMatrix tetraReferenceGradient()
{
    const Real gradient[4][3] = {{-1, -1, -1}, {1, 0, 0}, {0, 1, 0}, {0, 0, 1}};
    DenseContactMatrix b(9, StdVec<Real>(12, 0));
    for (size_t node = 0; node < 4; ++node)
        for (size_t i = 0; i < 3; ++i)
            for (size_t j = 0; j < 3; ++j) b[i * 3 + j][node * 3 + i] = gradient[node][j];
    return b;
}

inline Matd tetraDeformation(const StdVec<Real> &position)
{
    if (position.size() != 12) throw std::invalid_argument("tetra control: position dimensions");
    const auto b = tetraReferenceGradient(); Matd f;
    for (size_t i = 0; i < 9; ++i) f(i / 3, i % 3) = denseContactDot(b[i], position);
    return f;
}

struct TetraGelProblem
{
    StdVec<Real> position{0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1};
    StdVec<Real> velocity = StdVec<Real>(12, 0);
    StdVec<Real> node_mass{.01, .01, .01, .01};
    GelParameters material{5, 1, .5, .3};
    GelHistory history;
    Real volume = 1.0 / 6, dt = .025;
    StdVec<DenseContactConstraint> contacts;
    PrimalDualSettings settings;
};

struct TetraGelEvaluation
{
    NonlinearContactEquation equation;
    GelAffinePathResult material;
    Matd deformation;
    StdVec<Real> internal_force;
    Real transpose_work_residual = 0, jacobian_asymmetry = 0;
};

inline Real validateTetraGelProblem(const TetraGelProblem &p)
{
    p.material.validate();
    if (p.position.size() != 12 || p.velocity.size() != 12 || p.node_mass.size() != 4 ||
        !std::isfinite(p.volume) || p.volume <= 0 || !std::isfinite(p.dt) || p.dt <= 0)
        throw std::invalid_argument("tetra control: invalid physical dimensions/step");
    Real total_mass = 0;
    for (Real mass : p.node_mass) { if (!std::isfinite(mass) || mass <= 0) throw std::invalid_argument("tetra control: invalid mass"); total_mass += mass; }
    for (const auto *values : {&p.position, &p.velocity}) for (Real value : *values)
        if (!std::isfinite(value)) throw std::invalid_argument("tetra control: nonfinite state");
    return total_mass;
}

inline TetraGelEvaluation evaluateTetraGel(const TetraGelProblem &p, const StdVec<Real> &v)
{
    validateTetraGelProblem(p);
    if (v.size() != 12) throw std::invalid_argument("tetra control: trial velocity dimensions");
    for (Real value : v) if (!std::isfinite(value)) throw std::invalid_argument("tetra control: nonfinite trial velocity");
    const auto b = tetraReferenceGradient();
    const Matd initial = tetraDeformation(p.position);
    StdVec<Real> next(12); for (size_t i = 0; i < 12; ++i) next[i] = p.position[i] + p.dt * v[i];
    TetraGelEvaluation r; r.deformation = tetraDeformation(next);
    r.material = advanceGelAffinePath(initial, r.deformation, p.dt, p.material, p.history);
    r.equation.momentum.resize(12); r.equation.jacobian.assign(12, StdVec<Real>(12, 0)); r.internal_force.resize(12);
    for (size_t i = 0; i < 12; ++i)
    {
        long double force = 0;
        for (size_t a = 0; a < 9; ++a) force -= static_cast<long double>(p.volume) * b[a][i] * r.material.average_pk1(a / 3, a % 3);
        r.internal_force[i] = force;
        const Real mass = p.node_mass[i / 3];
        r.equation.momentum[i] = mass * (v[i] - p.velocity[i]) - p.dt * r.internal_force[i];
        for (size_t j = 0; j < 12; ++j)
        {
            long double derivative = i == j ? mass : 0;
            for (size_t a = 0; a < 9; ++a)
                for (size_t c = 0; c < 9; ++c)
                    derivative += static_cast<long double>(p.dt) * p.dt * p.volume * b[a][i] *
                        r.material.endpoint_tangent[c](a / 3, a % 3) * b[c][j];
            r.equation.jacobian[i][j] = derivative;
        }
    }
    long double work = 0, norm = 0, asymmetry = 0;
    for (size_t i = 0; i < 12; ++i)
    {
        work -= static_cast<long double>(p.dt) * r.internal_force[i] * v[i];
        for (size_t j = 0; j < 12; ++j)
        {
            norm += std::pow(static_cast<long double>(r.equation.jacobian[i][j]), 2);
            asymmetry += std::pow(static_cast<long double>(r.equation.jacobian[i][j]) - r.equation.jacobian[j][i], 2);
        }
    }
    r.transpose_work_residual = work - p.volume * r.material.work;
    r.jacobian_asymmetry = norm > 0 ? std::sqrt(asymmetry / norm) : 0;
    return r;
}

struct TetraGelResult
{
    PrimalDualResult solve;
    TetraGelEvaluation evaluation;
    StdVec<Real> position, velocity, contact_impulse, momentum_residual;
    GelHistory accepted_history;
    Real kinetic_old = 0, kinetic_new = 0, algorithmic_loss = 0, normal_contact_loss = 0, wall_work = 0, energy_residual = 0;
    Vecd total_momentum_residual = Vecd::Zero(), internal_force_sum = Vecd::Zero();
    Vecd angular_change_minus_external = Vecd::Zero(), internal_torque_old_cursor = Vecd::Zero(), angular_balance_old_cursor = Vecd::Zero();
};

inline TetraGelResult advanceTetraGel(const TetraGelProblem &p)
{
    const Real total_mass = validateTetraGelProblem(p);
    const Matd initial = tetraDeformation(p.position);
    auto domain = [&](const StdVec<Real> &v)
    {
        StdVec<Real> next(12);
        for (size_t i = 0; i < 12; ++i) { next[i] = p.position[i] + p.dt * v[i]; if (!std::isfinite(next[i])) return false; }
        return gelAffineMinimumDeterminant(initial, tetraDeformation(next) - initial) >= .05;
    };
    TetraGelResult r; r.position = p.position; r.velocity = p.velocity; r.accepted_history = p.history;
    r.solve = solveNonsymmetricContact(p.position, p.dt, total_mass, p.contacts,
        [&](const StdVec<Real> &v) { return evaluateTetraGel(p, v).equation; }, domain, p.settings);
    if (!r.solve.converged) return r; // No candidate history/position is committed.
    r.velocity = r.solve.velocity; r.position.resize(12); r.contact_impulse.assign(12, 0);
    r.evaluation = evaluateTetraGel(p, r.velocity); r.accepted_history = r.evaluation.material.history;
    for (size_t s = 0; s < p.contacts.size(); ++s)
    {
        const auto &row = p.contacts[s]; const Real impulse = r.solve.normal_impulse[s];
        for (size_t i = 0; i < 12; ++i) r.contact_impulse[i] += row.coefficients[i] * impulse;
        r.wall_work += impulse * (row.new_offset - row.old_offset) / p.dt;
        r.normal_contact_loss += impulse * (denseContactDot(row.coefficients, p.position) - row.old_offset) / p.dt;
    }
    r.momentum_residual.resize(12);
    for (size_t i = 0; i < 12; ++i)
    {
        const Real mass = p.node_mass[i / 3], dv = r.velocity[i] - p.velocity[i];
        r.position[i] = p.position[i] + p.dt * r.velocity[i];
        r.kinetic_old += mass * p.velocity[i] * p.velocity[i] / 2;
        r.kinetic_new += mass * r.velocity[i] * r.velocity[i] / 2; r.algorithmic_loss += mass * dv * dv / 2;
        r.momentum_residual[i] = r.evaluation.equation.momentum[i] - r.contact_impulse[i];
    }
    for (size_t a = 0; a < 4; ++a)
    {
        Vecd q0, q1, v0, v1, force, impulse;
        for (size_t i = 0; i < 3; ++i)
        {
            q0[i] = p.position[a * 3 + i]; q1[i] = r.position[a * 3 + i];
            v0[i] = p.velocity[a * 3 + i]; v1[i] = r.velocity[a * 3 + i];
            force[i] = r.evaluation.internal_force[a * 3 + i]; impulse[i] = r.contact_impulse[a * 3 + i];
        }
        r.internal_force_sum += force;
        r.total_momentum_residual += p.node_mass[a] * (v1 - v0) - impulse;
        r.internal_torque_old_cursor += q0.cross(force);
        r.angular_change_minus_external += p.node_mass[a] * (q1.cross(v1) - q0.cross(v0)) - q0.cross(impulse);
    }
    r.angular_balance_old_cursor = r.angular_change_minus_external - p.dt * r.internal_torque_old_cursor;
    r.energy_residual = r.wall_work - (r.kinetic_new - r.kinetic_old +
        p.volume * (r.evaluation.material.stored_new - r.evaluation.material.stored_old + r.evaluation.material.physical_dissipation)) -
        r.algorithmic_loss - r.normal_contact_loss;
    const Real impulse_scale = total_mass * p.settings.velocity_scale;
    const Real energy_scale = std::max(impulse_scale * p.settings.velocity_scale, p.volume * std::max(p.material.bulk, p.material.equilibrium_shear));
    Real length_scale = 1; for (Real value : p.position) length_scale = std::max(length_scale, std::abs(value));
    if (!std::isfinite(r.energy_residual) || std::abs(r.energy_residual) > 1e-10 * energy_scale ||
        std::abs(r.evaluation.transpose_work_residual) > 1e-10 * energy_scale || r.normal_contact_loss < -1e-10 * energy_scale)
        throw std::runtime_error("tetra control: precommit energy/transpose work gate");
    for (Real value : r.momentum_residual) if (!std::isfinite(value) || std::abs(value) > 1e-10 * impulse_scale)
        throw std::runtime_error("tetra control: precommit coordinate momentum gate");
    for (size_t i = 0; i < 3; ++i)
        if (std::abs(r.total_momentum_residual[i]) > 1e-10 * impulse_scale ||
            std::abs(r.internal_force_sum[i]) > 1e-10 * impulse_scale / p.dt ||
            std::abs(r.angular_balance_old_cursor[i]) > 1e-10 * impulse_scale * length_scale)
            throw std::runtime_error("tetra control: precommit global transfer/cursor ledger gate");
    return r;
}
} // namespace SPH::peridynamics
#endif
