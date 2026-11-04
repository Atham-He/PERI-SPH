/* Licensed under the Apache License, Version 2.0. */
#ifndef PD_SURFACE_QUADRATURE_H
#define PD_SURFACE_QUADRATURE_H

#include "base_data_type.h"
#include <Eigen/Eigenvalues>
#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace SPH::peridynamics
{
using SurfaceVector = Eigen::Matrix<Real, 3, 1>;

/** Massless reference-area samples, independent of the volume-particle radius.
 * Two complete hemispheres are integrated. z^2 clusters Gauss points near the
 * poles, where small-indentation plate contact first develops. */
struct SphereSurfaceRule
{
    StdVec<SurfaceVector> position;
    StdVec<Real> area;

    SphereSurfaceRule(const SurfaceVector &center, Real radius, size_t rings, size_t azimuths, size_t polar_levels = 0)
    {
        if (!center.allFinite() || !std::isfinite(radius) || radius <= 0 ||
            rings < 4 || rings > 256 || azimuths < 8 || azimuths > 256 || polar_levels > 20)
            throw std::invalid_argument("invalid sphere surface quadrature");
        const Real pi = std::acos(Real(-1));
        StdVec<Real> nodes, weights;
        for (size_t k = 0; k < rings; ++k)
        {
            Real root = std::cos(pi * (Real(k) + 0.75) / (Real(rings) + 0.5));
            Real derivative = 0;
            bool converged = false;
            for (int iteration = 0; iteration < 32; ++iteration)
            {
                Real p = 1, previous = 0;
                for (size_t n = 1; n <= rings; ++n)
                {
                    const Real older = previous;
                    previous = p;
                    p = ((2 * Real(n) - 1) * root * previous - (Real(n) - 1) * older) / Real(n);
                }
                derivative = Real(rings) * (root * p - previous) / (root * root - 1);
                const Real change = p / derivative;
                root -= change;
                if (std::abs(change) < 2e-15) { converged = true; break; }
            }
            if (!converged) throw std::runtime_error("Gauss rule did not converge");
            nodes.push_back(0.5 * (1 + root));
            weights.push_back(1 / ((1 - root * root) * derivative * derivative));
        }
        // Fixed geometric panels resolve small contact caps without adapting or
        // refreshing the reference surface during a run. levels=0 is the old rule.
        for (size_t panel = 0; panel <= polar_levels; ++panel)
        {
            const Real lower = panel == 0 ? 0 : std::ldexp(Real(1), int(panel) - int(polar_levels) - 1);
            const Real upper = std::ldexp(Real(1), int(panel) - int(polar_levels));
            for (size_t k = 0; k < rings; ++k)
            {
                const Real z = lower + (upper - lower) * nodes[k];
                const Real weight = (upper - lower) * weights[k];
                const Real mu = 1 - z * z;
                const Real radial = radius * std::sqrt((1 - mu) * (1 + mu));
                for (int side : {-1, 1})
                    for (size_t a = 0; a < azimuths; ++a)
                    {
                        const Real angle = 2 * pi * (Real(a) + 0.5) / Real(azimuths);
                        position.push_back(center + SurfaceVector(radial * std::cos(angle), side * radius * mu,
                                                                  radial * std::sin(angle)));
                        area.push_back(2 * radius * radius * z * weight * 2 * pi / Real(azimuths));
                    }
            }
        }
    }
};

/** Fixed reference MLS map x_s=sum(phi_si*x_i), also used in the transpose
 * force transfer. Signed weights are retained: clipping would destroy affine
 * reproduction, torque balance and the contact-potential gradient identity.
 * No nearest-particle fallback or reference refresh is available. */
class ReferenceSurfaceTransfer
{
  public:
    StdVec<size_t> offsets{0}, targets;
    StdVec<Real> weights, areas;
    StdVec<SurfaceVector> reference;
    size_t particle_count = 0;
    uint64_t generation = 1;

    ReferenceSurfaceTransfer(const SurfaceVector *particles, size_t count,
                             const SphereSurfaceRule &surface, Real support)
        : areas(surface.area), reference(surface.position), particle_count(count)
    {
        if (!particles || count < 4 || !std::isfinite(support) || support <= 0)
            throw std::invalid_argument("invalid reference surface support");
        for (size_t i = 0; i < count; ++i)
            if (!particles[i].allFinite()) throw std::invalid_argument("non-finite reference particle");
        using Vector4 = Eigen::Matrix<Real, 4, 1>;
        using Matrix4 = Eigen::Matrix<Real, 4, 4>;
        for (size_t s = 0; s < reference.size(); ++s)
        {
            Matrix4 moment = Matrix4::Zero();
            StdVec<Vector4> basis;
            StdVec<Real> kernel;
            for (size_t i = 0; i < count; ++i)
            {
                const SurfaceVector relative = (particles[i] - reference[s]) / support;
                const Real distance2 = relative.squaredNorm();
                if (distance2 >= 1) continue;
                Vector4 q; q << 1, relative;
                const Real value = std::pow(1 - distance2, 3);
                targets.push_back(i); basis.push_back(q); kernel.push_back(value);
                moment += value * q * q.transpose();
            }
            Eigen::SelfAdjointEigenSolver<Matrix4> spectrum(moment);
            if (spectrum.info() != Eigen::Success || spectrum.eigenvalues()[0] <=
                1e-12 * spectrum.eigenvalues()[3])
                throw std::runtime_error("rank-deficient reference surface support");
            const Vector4 coefficients = moment.ldlt().solve(Vector4(1, 0, 0, 0));
            Real partition = 0;
            SurfaceVector reproduced = SurfaceVector::Zero();
            const size_t begin = offsets.back();
            for (size_t j = 0; j < basis.size(); ++j)
            {
                const Real phi = kernel[j] * basis[j].dot(coefficients);
                weights.push_back(phi); partition += phi;
                reproduced += phi * (particles[targets[begin + j]] - reference[s]);
            }
            if (!coefficients.allFinite() || std::abs(partition - 1) > 1e-10 || reproduced.norm() / support > 1e-10)
                throw std::runtime_error("surface transfer failed affine reproduction");
            offsets.push_back(targets.size());
        }
    }

    SurfaceVector gather(size_t sample, const SurfaceVector *values) const
    {
        SurfaceVector value = SurfaceVector::Zero();
        for (size_t j = offsets.at(sample); j < offsets.at(sample + 1); ++j)
            value += weights[j] * values[targets[j]];
        return value;
    }

    void scatter(size_t sample, const SurfaceVector &value, SurfaceVector *particles) const
    {
        for (size_t j = offsets.at(sample); j < offsets.at(sample + 1); ++j)
            particles[targets[j]] += weights[j] * value;
    }

    /** Register every reference cache entry. Restore replaces stored values;
     * callers must never rebuild the map after loading an evolved checkpoint. */
    template<class Checkpoint> void registerCheckpoint(Checkpoint &checkpoint)
    {
        checkpoint.addScalar("SurfaceGeneration", generation);
        checkpoint.addScalar("SurfaceParticleCount", particle_count);
        checkpoint.addArray("SurfaceOffsets", offsets);
        checkpoint.addArray("SurfaceTargets", targets);
        checkpoint.addArray("SurfaceWeights", weights);
        checkpoint.addArray("SurfaceReference", reference);
        checkpoint.addArray("SurfaceAreas", areas);
    }
};

struct SurfacePlaneResponse
{
    SurfaceVector force = SurfaceVector::Zero(), friction = SurfaceVector::Zero();
    Real penetration = 0, energy = 0, boundary_power = 0, friction_boundary_power = 0;
};

/** Plane normal points into the permitted solid region. stiffness has units
 * Pa/m; it is a numerical penalty modulus divided by a declared physical length.
 * Reference area is unchanged by deformation (total-Lagrangian contact energy).
 * There is no finite particle contact radius and no force normalization. */
inline SurfacePlaneResponse surfacePlaneResponse(const SurfaceVector &position, const SurfaceVector &velocity,
    const SurfaceVector &normal, Real offset, const SurfaceVector &wall_velocity,
    Real area, Real stiffness, Real mu, Real veps)
{
    if (!position.allFinite() || !velocity.allFinite() || !normal.allFinite() || !wall_velocity.allFinite() ||
        !std::isfinite(offset) || !std::isfinite(area) || area <= 0 || !std::isfinite(stiffness) || stiffness <= 0 ||
        !std::isfinite(mu) || mu < 0 || !std::isfinite(veps) || veps <= 0 || std::abs(normal.norm() - 1) > 1e-12)
        throw std::invalid_argument("invalid surface-plane contact inputs");
    SurfacePlaneResponse result;
    result.penetration = std::max(Real(0), offset - normal.dot(position));
    const Real magnitude = stiffness * area * result.penetration;
    result.force = magnitude * normal;
    const SurfaceVector relative = velocity - wall_velocity;
    const SurfaceVector tangent = relative - normal * relative.dot(normal);
    const Real speed = tangent.norm();
    if (speed > 0)
        result.friction = -mu * magnitude * std::tanh(speed / veps) * tangent / speed;
    result.force += result.friction;
    result.energy = 0.5 * stiffness * area * result.penetration * result.penetration;
    result.boundary_power = result.force.dot(wall_velocity);
    result.friction_boundary_power = result.friction.dot(wall_velocity);
    return result;
}
} // namespace SPH::peridynamics
#endif
