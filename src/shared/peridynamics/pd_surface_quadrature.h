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

/** Pole-inclusive positive-area fixed spherical quadrature for bounded
 * constraint integration. Composite three-point Lobatto panels integrate
 * dmu exactly through degree3, with a predeclared narrow polar cap. Shared
 * latitudes are merged and each true pole has one positive-area sample.
 * Coarse quadrature is not a surface-convergence certificate. */
struct PolePanelSphereRule
{
    StdVec<SurfaceVector> position;
    StdVec<Real> area;
    PolePanelSphereRule(const SurfaceVector &center,Real radius,size_t azimuths=8,Real cap_width=1./64)
    {
        if(!center.allFinite() || !std::isfinite(radius) || radius<=0 || azimuths<8 || azimuths>256 ||
            !std::isfinite(cap_width) || cap_width<=0 || cap_width>=.5)throw std::invalid_argument("invalid pole-panel sphere rule");
        const Real pi=std::acos(Real(-1));
        const Real mu[5]={0,(1-cap_width)/2,1-cap_width,1-cap_width/2,1};
        const Real weights[5]={(1-cap_width)/6,2*(1-cap_width)/3,1./6,2*cap_width/3,cap_width/6};
        for(size_t k=0;k<5;++k)for(int side:{-1,1})
        {
            if(k==0 && side==1)continue;
            const Real latitude_weight=k==0?2*weights[k]:weights[k];
            if(k==4)
            {
                position.push_back(center+SurfaceVector(0,side*radius,0));area.push_back(2*pi*radius*radius*latitude_weight);continue;
            }
            const Real radial=radius*std::sqrt((1-mu[k])*(1+mu[k]));
            for(size_t a=0;a<azimuths;++a)
            {
                const Real angle=2*pi*(Real(a)+.5)/azimuths;
                position.push_back(center+SurfaceVector(radial*std::cos(angle),side*radius*mu[k],radial*std::sin(angle)));
                area.push_back(2*pi*radius*radius*latitude_weight/azimuths);
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
    Real reference_support=0;

    ReferenceSurfaceTransfer()=default;

    ReferenceSurfaceTransfer(const SurfaceVector *particles, size_t count,
                             const SphereSurfaceRule &surface, Real support)
        :ReferenceSurfaceTransfer(particles,count,surface.position,surface.area,support){}

    ReferenceSurfaceTransfer(const SurfaceVector *particles,size_t count,const PolePanelSphereRule &surface,Real support)
        :ReferenceSurfaceTransfer(particles,count,surface.position,surface.area,support){}

    ReferenceSurfaceTransfer(const SurfaceVector *particles,size_t count,const StdVec<SurfaceVector> &surface_position,
        const StdVec<Real> &surface_area,Real support)
        : areas(surface_area), reference(surface_position), particle_count(count),reference_support(support)
    {
        if (!particles || count < 4 || reference.empty() || areas.size()!=reference.size() || !std::isfinite(support) || support <= 0)
            throw std::invalid_argument("invalid reference surface support");
        for (size_t i = 0; i < count; ++i)
            if (!particles[i].allFinite()) throw std::invalid_argument("non-finite reference particle");
        using Vector4 = Eigen::Matrix<Real, 4, 1>;
        using Matrix4 = Eigen::Matrix<Real, 4, 4>;
        for (size_t s = 0; s < reference.size(); ++s)
        {
            if(!reference[s].allFinite() || !std::isfinite(areas[s]) || areas[s]<=0)throw std::invalid_argument("invalid reference surface sample/area");
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

    /** Check stored reference data without modifying/rebuilding it. The input
     * is the frozen PARTICLE REFERENCE, never accepted deformed positions. */
    void validateStored(const SurfaceVector *particles,size_t count)const
    {
        if(!particles || count!=particle_count || count<4 || generation!=1 || !std::isfinite(reference_support) || reference_support<=0 ||
            reference.empty() || areas.size()!=reference.size() || offsets.size()!=reference.size()+1 || offsets.front()!=0 || offsets.back()!=targets.size() ||
            weights.size()!=targets.size())throw std::invalid_argument("incompatible stored surface transfer");
        for(size_t i=0;i<count;++i)if(!particles[i].allFinite())throw std::invalid_argument("nonfinite stored particle reference");
        using Vector4=Eigen::Matrix<Real,4,1>;using Matrix4=Eigen::Matrix<Real,4,4>;
        for(size_t s=0;s<reference.size();++s)
        {
            if(!reference[s].allFinite() || !std::isfinite(areas[s]) || areas[s]<=0 || offsets[s]>=offsets[s+1] || offsets[s+1]>targets.size())
                throw std::invalid_argument("invalid stored surface row");
            Matrix4 moment=Matrix4::Zero();StdVec<Vector4> basis;StdVec<Real> kernel;StdVec<size_t> expected_targets;
            for(size_t i=0;i<count;++i)
            {
                const SurfaceVector relative=(particles[i]-reference[s])/reference_support;const Real distance2=relative.squaredNorm();
                if(distance2>=1)continue;Vector4 q;q<<1,relative;const Real value=std::pow(1-distance2,3);
                moment+=value*q*q.transpose();basis.push_back(q);kernel.push_back(value);expected_targets.push_back(i);
            }
            if(expected_targets.size()!=offsets[s+1]-offsets[s])throw std::invalid_argument("stored surface support membership mismatch");
            Eigen::SelfAdjointEigenSolver<Matrix4> spectrum(moment);
            if(spectrum.info()!=Eigen::Success || spectrum.eigenvalues()[0]<=1e-12*spectrum.eigenvalues()[3])
                throw std::invalid_argument("stored surface shape rank failure");
            const Vector4 coefficients=moment.ldlt().solve(Vector4(1,0,0,0));Real partition=0;SurfaceVector reproduce=SurfaceVector::Zero();
            for(size_t k=0;k<expected_targets.size();++k)
            {
                const size_t p=offsets[s]+k;const Real expected=kernel[k]*basis[k].dot(coefficients);
                if(targets[p]!=expected_targets[k] || !std::isfinite(weights[p]) || std::abs(weights[p]-expected)>1e-10*(1+std::abs(expected)))
                    throw std::invalid_argument("stored signed surface weight mismatch");
                partition+=weights[p];reproduce+=weights[p]*(particles[targets[p]]-reference[s]);
            }
            if(std::abs(partition-1)>1e-10 || reproduce.norm()/reference_support>1e-10)throw std::invalid_argument("stored surface affine mismatch");
        }
    }

    /** Register every reference cache entry. Restore replaces stored values;
     * callers must never rebuild the map after loading an evolved checkpoint. */
    template<class Checkpoint> void registerCheckpoint(Checkpoint &checkpoint)
    {
        checkpoint.addScalar("SurfaceGeneration", generation);
        checkpoint.addScalar("SurfaceParticleCount", particle_count);
        checkpoint.addScalar("SurfaceReferenceSupport",reference_support);
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
