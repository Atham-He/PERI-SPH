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

#include "pd_nosb_bond_batch.h"
#include "pd_nosb_elastic_batch.h"
#include "pd_nosb_j2_batch.h"
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
    Real prepareKinematics(std::size_t i, Matd &gradient);
    void publishPlasticState(std::size_t i, Real determinant, const Matd &rotated_plastic);
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

class NosbPDBondWorkspace
{
  public:
    NosbPDBondWorkspace(BaseInnerRelation &relation, std::unique_ptr<peridynamics::NosbPDBondBatch> backend)
        : relation_(relation), backend_(std::move(backend)) {}

    bool onDevice() const { return static_cast<bool>(backend_); }
    BaseInnerRelation &relation() const { return relation_; }
    void evaluateDeformation(Real dt, Matd *full, Matd *increment, Matd *midpoint);
    void evaluateInternalForce(Real reference_density, Vecd *acceleration, Vecd *force);
    void evaluateStabilization(peridynamics::NosbPDStabilizationParameters parameters, Real reference_density,
                               Vecd *acceleration, Vecd *force);
    void evaluateDamping(peridynamics::NosbPDDampingParameters parameters, Real reference_density,
                         Vecd *acceleration, Vecd *force);

  private:
    void configure(std::size_t count, std::size_t bound);
    void refreshLiveBonds(std::size_t count, std::size_t bound);
    void gatherCorrectionSites(bool damping);
    void publishCorrection(bool stabilization, Real damage_cutoff, Real reference_density, Vecd *acceleration, Vecd *force);

    BaseInnerRelation &relation_;
    std::unique_ptr<peridynamics::NosbPDBondBatch> backend_;
    std::vector<std::size_t> offsets_, ghost_sources_;
    std::vector<peridynamics::NosbPDKinematicSite> sites_;
    std::vector<peridynamics::NosbPDForceSite> force_sites_;
    std::vector<std::uint8_t> bond_live_;
    std::vector<peridynamics::NosbPDDeformationTensors> tensors_;
    std::vector<peridynamics::ParticleVector3> force_density_;
    std::vector<peridynamics::NosbPDCorrectionSite> correction_sites_;
    std::vector<peridynamics::NosbPDCorrectionForce> correction_force_;
    std::size_t bound_ = 0;
};

template <class Stage>
class NosbPDDeformationWithExecution : public InteractionWithUpdate<Stage>
{
    using Base = InteractionWithUpdate<Stage>;

  public:
    NosbPDDeformationWithExecution(BaseInnerRelation &relation, NosbPDBondWorkspace *workspace)
        : Base(relation), workspace_(workspace)
    {
        if (workspace && &workspace->relation() != &relation)
            throw std::invalid_argument("NOSB-PD deformation workspace must use the same reference relation");
    }

    bool onDevice() const { return workspace_ && workspace_->onDevice(); }

    void runMainStep(Real dt) override
    {
        if (onDevice())
            workspace_->evaluateDeformation(dt, this->N_, this->N_deltaU_, this->N_half_);
        else
            Base::runMainStep(dt);
    }

  private:
    NosbPDBondWorkspace *workspace_;
};

class NosbPDElasticWithExecution : public NosbPDDeformationWithExecution<NosbPDSecondStep>
{
  public:
    NosbPDElasticWithExecution(BaseInnerRelation &relation, NosbPDBondWorkspace *workspace,
                               std::unique_ptr<peridynamics::NosbPDElasticBatch> backend);
    void exec(Real dt = 0.0) override;

  private:
    std::unique_ptr<peridynamics::NosbPDElasticBatch> backend_;
    std::vector<peridynamics::NosbPDElasticSite> sites_;
    std::vector<peridynamics::NosbPDElasticResult> results_;
};

class NosbPDJ2WithExecution : public NosbPDDeformationWithExecution<NosbPDSecondStepPlastic>
{
  public:
    NosbPDJ2WithExecution(BaseInnerRelation &relation, NosbPDBondWorkspace *workspace,
                          std::unique_ptr<peridynamics::NosbPDJ2Batch> backend);
    void exec(Real dt = 0.0) override;

  private:
    std::unique_ptr<peridynamics::NosbPDJ2Batch> backend_;
    peridynamics::NosbPDJ2Parameters parameters_{};
    Real *alpha_ = nullptr;
    Matd *kinematic_ = nullptr;
    std::vector<Real> determinants_;
    std::vector<peridynamics::NosbPDJ2Site> sites_;
    std::vector<peridynamics::NosbPDJ2Result> results_;
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

class NosbPDInternalForceWithExecution : public InteractionDynamics<NosbPDThirdStep>
{
    using Base = InteractionDynamics<NosbPDThirdStep>;

  public:
    explicit NosbPDInternalForceWithExecution(NosbPDBondWorkspace &workspace)
        : Base(workspace.relation()), workspace_(workspace) {}

    void runMainStep(Real dt) override
    {
        if (workspace_.onDevice())
            workspace_.evaluateInternalForce(rho0_, acc_, internal_force_);
        else
            Base::runMainStep(dt);
    }

  private:
    NosbPDBondWorkspace &workspace_;
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

class NosbPDStabilizationWithExecution : public InteractionDynamics<ProjectedPairStabilization>
{
    using Base = InteractionDynamics<ProjectedPairStabilization>;

  public:
    NosbPDStabilizationWithExecution(NosbPDBondWorkspace &workspace, Real coefficient = 0.5, Real damage_cutoff = 0.1)
        : Base(workspace.relation(), coefficient, damage_cutoff), workspace_(workspace) {}
    void runMainStep(Real dt) override;

  private:
    NosbPDBondWorkspace &workspace_;
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

class NosbPDDampingWithExecution : public InteractionDynamics<PairNumericalDampingforPD>
{
    using Base = InteractionDynamics<PairNumericalDampingforPD>;

  public:
    explicit NosbPDDampingWithExecution(NosbPDBondWorkspace &workspace);
    bool onDevice() const { return supported_ && workspace_.onDevice(); }
    void runMainStep(Real dt) override;

  private:
    NosbPDBondWorkspace &workspace_;
    bool supported_;
};
} // namespace solid_dynamics
} // namespace SPH
#endif // NOSBPD_DYNAMICS_H
