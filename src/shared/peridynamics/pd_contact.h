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

namespace SPH
{
namespace solid_dynamics
{
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
    Real offset_W_ij_;
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
    Real contact_stiffness_;
    Real friction_mu_, friction_veps_;
    Vecd *vel_;
    Real *Vol_;
    Real *repulsion_factor_;
    Vecd *acc_prior_;
    Vecd *repulsion_force_;
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
    explicit ContactRepulsionToWallforPD(BaseContactRelation &contact_relation);
    virtual ~ContactRepulsionToWallforPD() {};

    void interaction(size_t index_i, Real dt = 0.0);

  protected:
    Real *Vol_; ///< volumetric measure of wall particles
    Vecd *repulsion_force_;
    StdVec<Solid *> contact_solids_;
    StdVec<Real *> contact_Vol_;
    StdVec<Real *> contact_repulsion_factor_;
};
} // namespace solid_dynamics
} // namespace SPH
#endif // PD_CONTACT_H
