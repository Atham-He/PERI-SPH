/* ------------------------------------------------------------------------- *
 *                                SPHinXsys                                  *
 * ------------------------------------------------------------------------- *
 * Licensed under the Apache License, Version 2.0 (the "License");           *
 * you may not use this file except in compliance with the License.          *
 *                                                                           *
 * ------------------------------------------------------------------------- */
/**
 * @file 	nosbpd_dynamics.h
 * @brief 	Non-ordinary state-based peridynamics (NOSB-PD) dynamics.
 * @details	Ported from Icebreaking_PD-SPH (Haotian Shi, Luwen Zhang, SJTU)
 * 			to the present SPHinXsys variable-registration framework.
 * 			The algorithm follows a position-Verlet scheme:
 * 			FirstStep -> (bond breaking) -> SecondStep -> hourglass control
 * 			-> numerical damping -> ThirdStep -> FourthStep.
 * 			NOTE: the peridynamic inner configuration is built once at
 * 			initialization and must NOT be updated during the simulation,
 * 			so that the reference bond geometry and bond states persist.
 */

#ifndef NOSBPD_DYNAMICS_H
#define NOSBPD_DYNAMICS_H

#include "base_body_relation.h"
#include "base_local_dynamics.h"
#include "base_particle_dynamics.h"
#include "external_force.h"
#include "pd_materials.h"

#include <Eigen/Eigenvalues>

#include "pd_solid_time_step.h"
#include "dynamics_algorithms.h"

namespace SPH
{
namespace solid_dynamics
{
/**
 * @class PDTimeStepInitialization
 * @brief Initialize the prior acceleration of PD particles, e.g. by gravity.
 */
class PDTimeStepInitialization : public LocalDynamics
{
  public:
    explicit PDTimeStepInitialization(SPHBody &sph_body, SharedPtr<Gravity> gravity_ptr);
    virtual ~PDTimeStepInitialization() {};

    void update(size_t index_i, Real dt = 0.0);

  protected:
    SharedPtr<Gravity> gravity_ptr_;
    Vecd *pos_, *acc_prior_;
};

/**
 * @class NosbPDShapeMatrix
 * @brief Compute the shape tensor K in the reference configuration.
 * Executed only once before the time integration loop.
 */
class NosbPDShapeMatrix : public LocalDynamics, public DataDelegateInner
{
  public:
    explicit NosbPDShapeMatrix(BaseInnerRelation &inner_relation);
    virtual ~NosbPDShapeMatrix() {};

    void interaction(size_t index_i, Real dt = 0.0);
    void update(size_t index_i, Real dt = 0.0);

  protected:
    Real *Vol_;
    Matd *shape_K_, *shape_K_1_;
};

/**
 * @class NosbPDFirstStep
 * @brief Position-Verlet first half step for NOSB-PD.
 */
class NosbPDFirstStep : public LocalDynamics
{
  public:
    explicit NosbPDFirstStep(SPHBody &sph_body);
    virtual ~NosbPDFirstStep() {};

    void update(size_t index_i, Real dt = 0.0);
    void evaluateDevice(peridynamics::PositionVerletBatch &backend, Real dt);

  protected:
    HughesWingetSolid &hw_solid_;
    Real rho0_;
    Real *rho_;
    int *particleLive_;
    Vecd *pos_, *vel_, *acc_old_, *acc_, *acc_prior_;
    Matd *F_;

  private:
    void prepare(std::size_t i);
    std::vector<peridynamics::PositionVerletStart> device_input_;
    std::vector<peridynamics::ParticleVector3> device_output_;
};

/**
 * @class NosbPDSecondStep
 * @brief Deformation measure, constitutive update and force state.
 * @details interaction: approximate deformation tensors N, N_deltaU, N_half
 * by bond summation; update: deformation gradients, Hughes-Winget stress
 * update, first Piola-Kirchhoff stress and peridynamic force state T0.
 */
class NosbPDSecondStep : public LocalDynamics, public DataDelegateInner
{
  public:
    explicit NosbPDSecondStep(BaseInnerRelation &inner_relation);
    virtual ~NosbPDSecondStep() {};

    void interaction(size_t index_i, Real dt = 0.0);
    void update(size_t index_i, Real dt = 0.0);

  protected:
    HughesWingetSolid &hw_solid_;
    int *particleLive_;
    Real *Vol_;
    Vecd *pos_, *vel_;
    Matd *F_, *shape_K_1_;
    Matd *N_, *N_deltaU_, *N_half_;
    Matd *F_half_, *F_delta_, *F_1_, *F_1_half_;
    Matd *PK1_, *T0_, *stress_;
};

/**
 * @class NosbPDSecondStepPlastic
 * @brief Plastic version of NosbPDSecondStep, calling the elastoplastic
 * constitutive relation (convex cutting-plane return mapping).
 */
class NosbPDSecondStepPlastic : public LocalDynamics, public DataDelegateInner
{
  public:
    explicit NosbPDSecondStepPlastic(BaseInnerRelation &inner_relation);
    virtual ~NosbPDSecondStepPlastic() {};

    void interaction(size_t index_i, Real dt = 0.0);
    void update(size_t index_i, Real dt = 0.0);

  protected:
    PlasticSolidforPD &plastic_solid_;
    int *particleLive_;
    Real *Vol_;
    Vecd *pos_, *vel_;
    Matd *F_, *shape_K_1_;
    Matd *N_, *N_deltaU_, *N_half_;
    Matd *F_half_, *F_delta_, *F_1_, *F_1_half_;
    Matd *PK1_, *T0_, *stress_, *plastic_strain_;
    Real *plastic_work_;
};

/**
 * @class NosbPDCheckBondLive
 * @brief Base class for bond breaking criteria.
 */
class NosbPDCheckBondLive : public LocalDynamics, public DataDelegateInner
{
  public:
    explicit NosbPDCheckBondLive(BaseInnerRelation &inner_relation);
    virtual ~NosbPDCheckBondLive() {};

  protected:
    Real critical_value_;
    int *particleLive_;
    Vecd *pos_, *vel_, *acc_;
};

/**
 * @class BondBreakByPrinStress
 * @brief Bond breaking by the maximum principal stress criterion.
 */
class BondBreakByPrinStress : public NosbPDCheckBondLive
{
  public:
    explicit BondBreakByPrinStress(BaseInnerRelation &inner_relation, Real critical_value);
    virtual ~BondBreakByPrinStress() {};

    virtual bool checkBondLive(Matd &stress_eq, Real &stretch_rate);
    void interaction(size_t index_i, Real dt = 0.0);

  protected:
    int *bond_;
    Real *damage_;
    Matd *stress_;
};

/**
 * @class BondBreakBySigma1andSigma3
 * @brief Bond breaking by the maximum principal stress and maximum shear stress.
 */
class BondBreakBySigma1andSigma3 : public BondBreakByPrinStress
{
  public:
    explicit BondBreakBySigma1andSigma3(BaseInnerRelation &inner_relation,
                                        Real max_tension, Real max_shear);
    virtual ~BondBreakBySigma1andSigma3() {};

    virtual bool checkBondLive(Matd &stress_eq, Real &stretch_rate);

  protected:
    Real max_tension_, max_shear_;
};

/**
 * @class BondBreakBySigma1andNorm1
 * @brief Bond breaking by the maximum principal stress and maximum pressure.
 */
class BondBreakBySigma1andNorm1 : public BondBreakByPrinStress
{
  public:
    explicit BondBreakBySigma1andNorm1(BaseInnerRelation &inner_relation,
                                       Real max_tension, Real max_pressure);
    virtual ~BondBreakBySigma1andNorm1() {};

    virtual bool checkBondLive(Matd &stress_eq, Real &stretch_rate);

  protected:
    Real max_tension_, max_pressure_;
};

/**
 * @class BondBreakByPlasticStrain
 * @brief Bond breaking by the maximum principal plastic strain criterion.
 */
class BondBreakByPlasticStrain : public LocalDynamics, public DataDelegateInner
{
  public:
    explicit BondBreakByPlasticStrain(BaseInnerRelation &inner_relation, Real critical_value);
    virtual ~BondBreakByPlasticStrain() {};

    bool checkBondLive(Matd &plastic_strain_eq, Real &stretch_rate);
    void interaction(size_t index_i, Real dt = 0.0);

  protected:
    Real critical_value_;
    int *particleLive_, *bond_;
    Real *damage_;
    Matd *plastic_strain_;
    Vecd *pos_, *vel_, *acc_;
};

/**
 * @class NosbPDThirdStep
 * @brief First term of the NOSB-PD dynamic equation: pairwise force density.
 */
class NosbPDThirdStep : public LocalDynamics, public DataDelegateInner
{
  public:
    explicit NosbPDThirdStep(BaseInnerRelation &inner_relation);
    virtual ~NosbPDThirdStep() {};

    void interaction(size_t index_i, Real dt = 0.0);

  protected:
    HughesWingetSolid &hw_solid_;
    Real rho0_;
    int *particleLive_;
    Real *Vol_;
    Vecd *acc_, *internal_force_;
    Matd *T0_;
};

/**
 * @class NosbPDFourthStep
 * @brief Position-Verlet second half step for NOSB-PD.
 */
class NosbPDFourthStep : public LocalDynamics
{
  public:
    explicit NosbPDFourthStep(SPHBody &sph_body);
    virtual ~NosbPDFourthStep() {};

    void update(size_t index_i, Real dt = 0.0);
    void evaluateDevice(peridynamics::PositionVerletBatch &backend, Real dt);

  protected:
    int *particleLive_;
    Vecd *pos_, *vel_, *acc_;
    Real *mass_;
    Vecd *prior_, *internal_, *hg_, *damping_, *viscous_, *friction_, *acoustic_;
    Real *internal_work_, *external_work_, *hg_work_, *damping_work_, *viscous_work_, *acoustic_work_;
    Real *boundary_power_, *friction_boundary_power_, *boundary_work_, *friction_loss_;

  private:
    std::vector<peridynamics::PositionVerletFinish> device_input_;
    std::vector<peridynamics::PositionVerletResult> device_output_;
};

template <class Stage>
class NosbPDStepWithExecution : public SimpleDynamics<Stage>
{
    using Base = SimpleDynamics<Stage>;

  public:
    NosbPDStepWithExecution(SPHBody &body, peridynamics::PositionVerletBatch *backend)
        : Base(body), backend_(backend) {}

    void exec(Real dt = 0.0) override
    {
        if (!backend_)
        {
            Base::exec(dt);
            return;
        }
        this->setUpdated(this->getSPHBody());
        this->setupDynamics(dt);
        this->evaluateDevice(*backend_, dt);
    }

  private:
    peridynamics::PositionVerletBatch *backend_;
};

/**
 * @class ProjectedPairStabilization
 * @brief Scale-consistent central pair penalty for non-affine deformation.
 * Endpoint averaging and symmetric damage eligibility preserve pair momentum.
 * The legacy LittleWoodHourGlassControl name remains a source-compatible alias.
 */
class ProjectedPairStabilization : public LocalDynamics, public DataDelegateInner
{
  public:
    explicit ProjectedPairStabilization(BaseInnerRelation &inner_relation, Real coefficient = 0.5, Real damage_cutoff = 0.1);
    virtual ~ProjectedPairStabilization() {};

    void interaction(size_t index_i, Real dt = 0.0);

  protected:
    HughesWingetSolid &hw_solid_;
    Real rho0_, Chg_, Kbulk_, horizon_, damage_cutoff_;
    int *particleLive_;
    Real *damage_;
    Real *Vol_;
    Vecd *pos_, *acc_, *hourglass_force_;
    Matd *F_;
};

using LittleWoodHourGlassControl = ProjectedPairStabilization;

/**
 * @class PairNumericalDampingforPD
 * @brief Kelvin-Voigt damper with artificial viscosity for NOSB-PD.
 */
class PairNumericalDampingforPD : public LocalDynamics, public DataDelegateInner
{
  public:
    explicit PairNumericalDampingforPD(BaseInnerRelation &inner_relation);
    virtual ~PairNumericalDampingforPD() {};

    void interaction(size_t index_i, Real dt = 0.0);

  protected:
    HughesWingetSolid &hw_solid_;
    Real rho0_, c0_, numerical_dissipation_factor_, smoothing_length_, inv_W0_;
    int *particleLive_;
    Real *Vol_;
    Vecd *pos_, *vel_, *acc_, *damping_force_;
    Matd *F_;
};
} // namespace solid_dynamics
} // namespace SPH
#endif // NOSBPD_DYNAMICS_H
