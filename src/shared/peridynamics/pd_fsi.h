/* ------------------------------------------------------------------------- *
 *                                SPHinXsys                                  *
 * ------------------------------------------------------------------------- *
 * Licensed under the Apache License, Version 2.0 (the "License");           *
 * you may not use this file except in compliance with the License.          *
 *                                                                           *
 * ------------------------------------------------------------------------- */
/**
 * @file 	pd_fsi.h
 * @brief 	Fluid-structure interaction for peridynamic bodies.
 * @details	Ported from Icebreaking_PD-SPH (Haotian Shi, Luwen Zhang, SJTU).
 * 			- fluid side: velocity divergence and acoustic damper computed
 * 			  as standalone dynamics so that upstream fluid classes are
 * 			  not modified. Execute VelocityDivergence* after the fluid
 * 			  2nd half step, and AcousticDamper* after the 1st half step.
 * 			- solid (PD) side: pressure force from the fluid acting on
 * 			  PD particles, with an acoustic damper and interface
 * 			  pressure recording.
 * 			NOTE: for the FSI contact neighbor search the fluid kernel must
 * 			be selected, which is automatic (max smoothing length rule)
 * 			when 1.3 * df >= 1.5075 * ds, i.e. df >= 1.16 * ds.
 */

#ifndef PD_FSI_H
#define PD_FSI_H

#include "base_body_relation.h"
#include "base_local_dynamics.h"
#include "base_particle_dynamics.h"
#include "contact_body_relation.h"
#include "pd_materials.h"
#include "riemann_solver.h"
#include "adaptation.h"

namespace SPH
{
/**
 * @class NeighborBuilderContactSelectKernel
 * @brief Contact neighbor builder with an explicitly selected kernel.
 * @details The default rule picks the kernel with the larger smoothing
 * length. For PD-fluid interaction the fluid kernel must always be used
 * (as in Icebreaking_PD-SPH); otherwise, when 1.5075*ds > 1.3*df, the PD
 * kernel would be picked and the FSI force scaling would be inconsistent.
 */
class NeighborBuilderContactSelectKernel : public NeighborBuilderContact
{
  public:
    NeighborBuilderContactSelectKernel(SPHBody &body, SPHBody &contact_body, bool use_source_kernel)
        : NeighborBuilderContact(body, contact_body)
    {
        kernel_ = use_source_kernel ? body.getSPHAdaptation().getKernel()
                                    : contact_body.getSPHAdaptation().getKernel();
    };
    virtual ~NeighborBuilderContactSelectKernel() {};
};

/**
 * @class ContactRelationSelectKernel
 * @brief Contact relation using an explicitly selected kernel for the
 * neighbor search. Use use_source_kernel = true for the fluid body
 * (fluid -> PD) and false for the PD body (PD -> fluid), so that the
 * fluid kernel is always used for the FSI pair.
 */
class ContactRelationSelectKernel : public ContactRelationCrossResolution
{
  public:
    explicit ContactRelationSelectKernel(SPHBody &sph_body, RealBodyVector contact_bodies,
                                         bool use_source_kernel);
    virtual ~ContactRelationSelectKernel() {};

    virtual void updateConfiguration() override;

  protected:
    UniquePtrsKeeper<NeighborBuilderContact> neighbor_builder_contacts_keeper_;
    StdVec<NeighborBuilderContact *> get_contact_neighbors_;
};
namespace fluid_dynamics
{
/**
 * @class VelocityDivergenceInner
 * @brief Velocity divergence of fluid particles from inner interactions,
 * u_div_i = - sum_j (v_i - v_j) . e_ij dW_ij V_j.
 * Execute after the fluid 2nd half step.
 */
class VelocityDivergenceInner : public LocalDynamics, public DataDelegateInner
{
  public:
    explicit VelocityDivergenceInner(BaseInnerRelation &inner_relation);
    virtual ~VelocityDivergenceInner() {};

    void interaction(size_t index_i, Real dt = 0.0);

  protected:
    Real *Vol_;
    Vecd *vel_;
    Real *u_div_;
};

/**
 * @class VelocityDivergenceWithWallforPD
 * @brief Wall (PD body) contribution to the fluid velocity divergence,
 * vel_in_wall = 2 * vel_ave_wall - vel_i. Execute after VelocityDivergenceInner.
 */
class VelocityDivergenceWithWallforPD : public LocalDynamics, public DataDelegateContact
{
  public:
    explicit VelocityDivergenceWithWallforPD(BaseContactRelation &contact_relation);
    virtual ~VelocityDivergenceWithWallforPD() {};

    void interaction(size_t index_i, Real dt = 0.0);

  protected:
    Real *Vol_, *u_div_;
    Vecd *vel_;
    StdVec<Real *> contact_Vol_;
    StdVec<Vecd *> contact_vel_ave_;
};

/**
 * @class AcousticDamperInnerforPD
 * @brief Acoustic damper on fluid particles from inner interactions,
 * v_i += dt * coeff * sum_j (u_div_i + u_div_j) dW_ij V_j e_ij / rho_i.
 * Explicit velocity kick after the fluid 1st half step. The diagnostic force
 * is stored separately and never written into the upstream Force accumulator.
 */
class AcousticDamperInnerforPD : public LocalDynamics, public DataDelegateInner
{
  public:
    explicit AcousticDamperInnerforPD(BaseInnerRelation &inner_relation);
    virtual ~AcousticDamperInnerforPD() {};

    void setCoeffAcousticDamper(Real rho0, Real c0, Real h)
    {
        coeff_acoustic_damper_ = 0.3 * c0 * rho0 * h;
    };
    void interaction(size_t index_i, Real dt = 0.0);

  protected:
    Real coeff_acoustic_damper_ = 0.0;
    Real *Vol_, *rho_, *u_div_;
    Real *mass_;
    Vecd *force_, *vel_;
};

/**
 * @class AcousticDamperWithWallforPD
 * @brief Acoustic damper on fluid particles from wall (PD body) interactions,
 * v_i += dt * coeff * sum_j (2 * u_div_i) dW_ij V_j e_ij / rho_i.
 * Explicit velocity kick after the fluid 1st half step. The diagnostic force
 * is stored separately and never written into the upstream Force accumulator.
 */
class AcousticDamperWithWallforPD : public LocalDynamics, public DataDelegateContact
{
  public:
    explicit AcousticDamperWithWallforPD(BaseContactRelation &contact_relation);
    virtual ~AcousticDamperWithWallforPD() {};

    void setCoeffAcousticDamper(Real rho0, Real c0, Real h)
    {
        coeff_acoustic_damper_ = 0.3 * c0 * rho0 * h;
    };
    void interaction(size_t index_i, Real dt = 0.0);

  protected:
    Real coeff_acoustic_damper_ = 0.0;
    Real *Vol_, *rho_, *u_div_;
    Real *mass_;
    Vecd *force_, *vel_;
    StdVec<Real *> contact_Vol_;
};
} // namespace fluid_dynamics

namespace solid_dynamics
{
/**
 * @class PressureForceFromFluidforPD
 * @brief Pressure force from the fluid acting on PD particles.
 * @details Riemann-type wall interaction with an external-acceleration
 * corrected wall pressure, plus an acoustic damper driven by the fluid
 * velocity divergence. Force / reference mass is added to
 * "PDAccelerationPrior" after PDTimeStepInitialization; gravity and
 * other external loads are preserved.
 * The interface pressure is recorded in "PressureFromFluid" for output.
 */
template <class RiemannSolverType = NoRiemannSolver>
class PressureForceFromFluidforPD : public LocalDynamics, public DataDelegateContact
{
  public:
    explicit PressureForceFromFluidforPD(BaseContactRelation &contact_relation);
    virtual ~PressureForceFromFluidforPD() {};

    void setCoeffAcousticDamper(Real rho0, Real c0, Real h)
    {
        coeff_acoustic_damper_ = 0.3 * c0 * rho0 * h;
    };

    void interaction(size_t index_i, Real dt = 0.0);

  protected:
    HughesWingetSolid &hw_solid_;
    Real coeff_acoustic_damper_ = 0.0;
    Real *Vol_;
    Vecd *vel_ave_, *acc_prior_, *acc_ave_, *n_;
    Vecd *force_from_fluid_;
    Real *pressure_from_fluid_;
    StdVec<Fluid *> contact_fluids_;
    StdVec<Real *> contact_rho_, contact_p_, contact_u_div_, contact_mass_, contact_Vol_;
    StdVec<Vecd *> contact_vel_, contact_force_prior_;
    StdVec<RiemannSolverType> riemann_solvers_;
};
using PressureForceAccelerationFromFluidforPD = PressureForceFromFluidforPD<NoRiemannSolver>;
} // namespace solid_dynamics
} // namespace SPH
#endif // PD_FSI_H
