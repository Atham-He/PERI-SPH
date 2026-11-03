/* ------------------------------------------------------------------------- *
 *                                SPHinXsys                                  *
 * ------------------------------------------------------------------------- *
 * Licensed under the Apache License, Version 2.0 (the "License");           *
 * you may not use this file except in compliance with the License.          *
 *                                                                           *
 * ------------------------------------------------------------------------- */
/**
 * @file 	pd_contact.h
 * @brief 	Solid-wall contact for peridynamic bodies, e.g. a PD capsule
 * 			compressed between rigid plates.
 * @details Penalty repulsion in the SPH-interpolation form, following
 * 			solid_dynamics::RepulsionForce (contact_repulsion.h) but acting on
 * 			PD particles. The force contribution accumulates onto
 * 			PDAccelerationPrior, whose baseline must be set beforehand by
 * 			PDTimeStepInitialization.
 */

#ifndef PD_CONTACT_H
#define PD_CONTACT_H

#include "base_body_relation.h"
#include "base_local_dynamics.h"
#include "base_particle_dynamics.h"
#include "pd_materials.h"
#include <cmath>
#include <stdexcept>

namespace SPH
{
namespace solid_dynamics
{
inline Vecd pdContactPairForce(const Vecd &normal_force, const Vecd &relative_velocity,
                               const Vecd &wall_normal, Real mu, Real veps)
{
    if (normal_force.squaredNorm() == 0.0 || mu == 0.0) return normal_force;
    Real length = wall_normal.norm();
    if (!wall_normal.allFinite() || length <= 0.0)
        throw std::runtime_error("PD friction requires initialized wall NormalDirection");
    Vecd normal = wall_normal / length;
    Real magnitude = std::abs(normal_force.dot(normal));
    Vecd tangent = relative_velocity - normal * relative_velocity.dot(normal);
    Real speed = tangent.norm();
    if (speed == 0.0) return normal_force;
    return normal_force - mu * magnitude * std::tanh(speed / veps) * tangent / speed;
}
inline Vecd pdContactPairForce(const Vecd &normal_force, const Vecd &relative_velocity, Real mu, Real veps)
{
    return pdContactPairForce(normal_force, relative_velocity, normal_force, mu, veps);
}

/**
 * @class ContactRepulsionFactorforPD
 * @brief Kernel-summation contact factor sigma in [0, ~1] for PD particles,
 * following solid_dynamics::RepulsionFactorSummation (repulsion_factor_summation.cpp).
 * The W offset at one reference spacing localizes the factor to particles
 * actually approaching the wall, so the penalty ramps in smoothly.
 */
class ContactRepulsionFactorforPD : public LocalDynamics, public DataDelegateContact
{
  public:
    explicit ContactRepulsionFactorforPD(BaseContactRelation &contact_relation);
    virtual ~ContactRepulsionFactorforPD() {};

    void interaction(size_t index_i, Real dt = 0.0);

  protected:
    StdVec<Real> offset_W_ij_;
    Real *repulsion_factor_;
    StdVec<Real *> contact_Vol_;
};

/**
 * @class ContactRepulsionFromWallforPD
 * @brief Penalty repulsion acting on PD particles from rigid wall particles.
 * The wall particles contribute with their volumetric measure only,
 * following RepulsionForce<Contact<Wall>> (contact_repulsion.cpp).
 * The contact pressure is scaled by RepulsionFactor, which must be
 * computed beforehand by ContactRepulsionFactorforPD.
 */
class ContactRepulsionFromWallforPD : public LocalDynamics, public DataDelegateContact
{
  public:
    /** friction_mu: smooth-Coulomb tangential friction coefficient (0 = frictionless,
     * backward compatible); friction_veps: regularization velocity scale [m/s]. */
    explicit ContactRepulsionFromWallforPD(BaseContactRelation &contact_relation,
                                           Real friction_mu = 0.0, Real friction_veps = 1.0e-3);
    virtual ~ContactRepulsionFromWallforPD() {};

    void interaction(size_t index_i, Real dt = 0.0);

  protected:
    HughesWingetSolid &hw_solid_;
    StdVec<Vecd *> contact_velocity_, contact_normal_;
    Real contact_stiffness_;
    Real friction_mu_, friction_veps_;
    Vecd *vel_;
    Real *Vol_;
    Real *repulsion_factor_;
    Vecd *acc_prior_;
    Vecd *repulsion_force_, *friction_force_;
    Real *boundary_power_, *friction_boundary_power_;
    StdVec<Real *> contact_Vol_;
};

/**
 * @class ContactRepulsionToWallforPD
 * @brief Reaction force exerted on rigid wall particles by the PD body.
 * Recorded per wall particle as "ContactRepulsionToWall" so that the total
 * plate reaction (force-displacement curve) can be obtained by
 * ReduceDynamics<QuantitySummation<Vecd>>.
 * The wall is assumed rigid and kinematically driven, so the force is
 * recorded only and does not enter any time integration.
 */
class ContactRepulsionToWallforPD : public LocalDynamics, public DataDelegateContact
{
  public:
    explicit ContactRepulsionToWallforPD(BaseContactRelation &contact_relation, Real friction_mu = 0.0, Real friction_veps = 1.0e-3);
    virtual ~ContactRepulsionToWallforPD() {};

    void interaction(size_t index_i, Real dt = 0.0);

  protected:
    Real *Vol_; ///< volumetric measure of wall particles
    Vecd *repulsion_force_, *velocity_, *normal_;
    Real friction_mu_, friction_veps_;
    StdVec<Vecd *> contact_velocity_;
    StdVec<Solid *> contact_solids_;
    StdVec<Real *> contact_Vol_;
    StdVec<Real *> contact_repulsion_factor_;
};

/** Signed distance and outward normal for an axis-aligned rigid box. */
inline std::pair<Real, Vecd> rigidBoxDistance(const Vecd &position, const Vecd &center, const Vecd &half_size)
{
    Vecd local = position - center;
    Vecd q = local.cwiseAbs() - half_size;
    Vecd outside = q.cwiseMax(Vecd::Zero());
    Real length = outside.norm();
    Vecd normal = Vecd::Zero();
    if (length > 0)
    {
        for (int d = 0; d < Dimensions; ++d) normal[d] = (local[d] >= 0 ? 1 : -1) * outside[d] / length;
        return {length, normal};
    }
    Eigen::Index axis; Real distance = q.maxCoeff(&axis);
    normal[axis] = local[axis] >= 0 ? 1 : -1;
    return {distance, normal};
}

/** One contact evaluation. Forces act on the deformable particle; box reaction
 * is their negative. Sampling never writes particle or integrator state. */
struct RigidBoxContactResponse
{
    Vecd force = Vecd::Zero();
    Vecd friction = Vecd::Zero();
    Real penetration = 0;
    Real potential_energy = 0;
    Real boundary_power = 0;
    Real friction_boundary_power = 0;
};

/** Accepted-state box totals, separate from the last applied integration load. */
struct RigidBoxContactObservation
{
    Vecd reaction = Vecd::Zero();
    Real potential_energy = 0;
    Real boundary_power = 0;
    Real max_penetration = 0;
    size_t contacting_particles = 0;
};

/** Conservative normal penalty against translating finite rigid boxes.
 * The particle's contact radius is dp/2, k_i=ContactStiffness*V_i/dp^2.
 * U_i=0.5*k_i*max(dp/2-signed_distance,0)^2. The force is -grad(U).
 * The source acceleration, box reaction, boundary work and regularized
 * Coulomb friction are computed from the same pair force. Rotating boxes
 * require a separate motion implementation.
 */
class RigidBoxContactforPD : public LocalDynamics
{
  public:
    RigidBoxContactforPD(SPHBody &body, Real friction_mu=0, Real friction_veps=1e-3);
    size_t addBox(const std::string &name, const Vecd &center, const Vecd &half_size, const Vecd &velocity=Vecd::Zero());
    void setKinematics(size_t box, const Vecd &center, const Vecd &velocity);
    void applyToCurrentAcceleration(bool value) {current_acceleration_ = value;}
    void applyAcceleration(bool value) {apply_acceleration_ = value;}
    RigidBoxContactResponse evaluate(size_t particle, size_t box) const;
    StdVec<RigidBoxContactObservation> observe() const;
    Real normalPenaltyModulus() const { return stiffness_; }
    void update(size_t i, Real dt=0);
  protected:
    struct Box {Vecd center, half_size, velocity; Vecd *reaction;};
    StdVec<Box> boxes_;
    bool current_acceleration_ = false, apply_acceleration_ = true;
    Real spacing_, stiffness_, rho0_, mu_, veps_;
    Vecd *acceleration_;
    Vecd *position_, *velocity_, *prior_, *force_, *friction_;
    Real *volume_, *energy_, *boundary_power_, *friction_boundary_power_;
};
} // namespace solid_dynamics
} // namespace SPH
#endif // PD_CONTACT_H
