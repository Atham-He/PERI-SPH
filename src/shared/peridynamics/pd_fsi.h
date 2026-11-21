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
#include "pd_execution.h"
#include "dynamics_algorithms.h"
#include "riemann_solver.h"
#include "adaptation.h"
#include "fluid_integration.hpp"
#include "viscous_dynamics.hpp"
#include "fluid_structure_interaction.hpp"

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

namespace peridynamics
{
class FluidNeighborConfiguration
{
  public:
    // Keep this scope within one fixed neighbor geometry and particle order.
    class Scope
    {
      public:
        ~Scope();
        Scope(const Scope &) = delete;
        Scope &operator=(const Scope &) = delete;

      private:
        friend class FluidNeighborConfiguration;
        explicit Scope(FluidNeighborConfiguration &configuration);
        FluidNeighborConfiguration &configuration_;
    };

    [[nodiscard]] Scope reuse() { return Scope(*this); }
    std::size_t generation() const { return depth_ ? generation_ : 0; }

  private:
    std::size_t depth_ = 0, generation_ = 0;
};

class FluidNeighborCache
{
  public:
    explicit FluidNeighborCache(const FluidNeighborConfiguration *configuration) : configuration_(configuration) {}

    bool needsConfiguration(std::size_t rows) const
    {
        const auto generation = activeGeneration();
        return generation == 0 || generation != generation_ || rows_ != rows;
    }

    void invalidate() { generation_ = 0; }
    void configured(std::size_t rows)
    {
        rows_ = rows;
        generation_ = activeGeneration();
    }

  private:
    std::size_t activeGeneration() const { return configuration_ ? configuration_->generation() : 0; }

    const FluidNeighborConfiguration *configuration_;
    std::size_t rows_ = 0, generation_ = 0;
};

class SPHViscosityWorkspace
{
  public:
    enum class Form { Inner, Wall, ReverseWall };

    void clearTargets() { targets_.clear(); }
    void addTarget(Form form, const ParticleConfiguration &configuration, const Vecd *velocity,
                   const Real *volume, std::size_t bound, Real viscosity, Real smoothing_length);
    void evaluate(FluidViscosityBatch &backend, std::size_t count, const Vecd *velocity,
                  const Real *volume, Vecd *force);

  private:
    struct Target
    {
        Form form;
        const ParticleConfiguration *configuration;
        const Vecd *velocity;
        const Real *volume;
        std::size_t bound;
        Real viscosity, smoothing_length;
    };
    std::vector<Target> targets_;
    std::vector<std::vector<std::size_t>> site_ids_;
    FluidViscousRows rows_;
    std::vector<ParticleVector3> forces_;
};

template <class LocalDynamicsType>
class ViscousInteractionWithUpdate : public InteractionWithUpdate<LocalDynamicsType>
{
    using Base = InteractionWithUpdate<LocalDynamicsType>;

  public:
    template <class... Args>
    explicit ViscousInteractionWithUpdate(FluidViscosityBatch *backend, Args &&...args)
        : Base(std::forward<Args>(args)...), backend_(backend) {}

  protected:
    void runMainStep(Real dt) override
    {
        if (backend_)
            this->evaluateDevice(*backend_, workspace_);
        else
            Base::runMainStep(dt);
    }

  private:
    FluidViscosityBatch *backend_;
    SPHViscosityWorkspace workspace_;
};
}

namespace fluid_dynamics
{
/** Wall pressure extrapolation is an explicit numerical control, shared by
 * the fluid and the reverse PD traction. factor=1 reproduces upstream;
 * factor=0 omits it on PD walls without altering wall motion. Rigid
 * non-PD walls retain the upstream factor=1.
 */
/** Explicit interface viscosity uses actual accepted-time PD velocities.
 * Pressure/density coupling keeps the original average velocity fields. */
class ViscousWallTermsforPD : public ViscousForce<Contact<Wall>, FixedViscosity, NoKernelCorrection>
{
  public:
    explicit ViscousWallTermsforPD(BaseContactRelation &walls)
        : ViscousForce<Contact<Wall>, FixedViscosity, NoKernelCorrection>(walls)
    {
        for (size_t k=0;k<contact_bodies_.size();++k)
            if (dynamic_cast<HughesWingetSolid *>(&contact_bodies_[k]->getMatterMaterial()))
                wall_vel_ave_[k]=contact_particles_[k]->registerStateVariableData<Vecd>("Velocity");
    }

    void addDeviceTargets(peridynamics::SPHViscosityWorkspace &workspace);
};
class ViscousForceWithWallforPD : public ViscousForceInner
{
  public:
    ViscousForceWithWallforPD(BaseInnerRelation &inner, BaseContactRelation &walls)
        : ViscousForceInner(inner), walls_(walls) {}
    void interaction(size_t i, Real dt=0)
    {ViscousForceInner::interaction(i,dt); walls_.interaction(i,dt);}
    void evaluateDevice(peridynamics::FluidViscosityBatch &backend,
                        peridynamics::SPHViscosityWorkspace &workspace);

  protected:
    ViscousWallTermsforPD walls_;
};
class PressureWallTermsforPD : public Integration1stHalf<Contact<Wall>, AcousticRiemannSolver, NoKernelCorrection>
{
  public:
    PressureWallTermsforPD(BaseContactRelation &relation, Real factor);
    void interaction(size_t i, Real dt=0);
  protected:
    StdVec<Real> acceleration_factors_;
};
class Integration1stHalfWithWallforPD : public Integration1stHalfInnerRiemann
{
  public:
    Integration1stHalfWithWallforPD(BaseInnerRelation &inner, BaseContactRelation &walls, Real factor=1)
        : Integration1stHalfInnerRiemann(inner), wall_terms_(walls, factor) {}
    void interaction(size_t i, Real dt=0)
    {
        Integration1stHalfInnerRiemann::interaction(i,dt); wall_terms_.interaction(i,dt);
    }
  protected:
    PressureWallTermsforPD wall_terms_;
};

class PDFluidAcousticWorkspace
{
  public:
    PDFluidAcousticWorkspace(BaseInnerRelation &inner, BaseContactRelation &wall,
                             std::unique_ptr<peridynamics::FluidAcousticBatch> backend, Real acceleration_factor,
                             const peridynamics::FluidNeighborConfiguration *configuration = nullptr);
    bool onDevice() const { return static_cast<bool>(backend_); }


  private:
    friend class PDFluidPressureHalf;
    friend class PDFluidDensityHalf;
    friend class PDFluidVelocityDivergence;
    friend class PDFluidAcousticDamping;
    void configure(std::size_t count, std::size_t bound);
    void compute(peridynamics::FluidAcousticStage stage);
    void evaluateHalf(peridynamics::FluidAcousticStage stage);

    BaseInnerRelation &inner_;
    BaseContactRelation &wall_;
    Real acceleration_factor_;
    std::unique_ptr<peridynamics::FluidAcousticBatch> backend_;
    peridynamics::FluidAcousticCoefficients coefficients_{};
    StdVec<Real> wall_acceleration_factors_;
    std::vector<peridynamics::FluidAcousticSite> fluid_sites_;
    std::vector<peridynamics::FluidAcousticWall> wall_sites_;
    std::vector<std::size_t> ghost_indices_, ghost_sources_;
    std::vector<std::vector<std::size_t>> wall_indices_;
    std::vector<std::pair<std::size_t, std::size_t>> wall_sources_;
    peridynamics::FluidNeighborCache neighbor_cache_;
    peridynamics::FluidAcousticStencil inner_stencil_, wall_stencil_;
    std::vector<peridynamics::FluidAcousticResult> results_;
};

class PDFluidPressureHalf : public Dynamics1Level<Integration1stHalfWithWallforPD>
{
    using Base = Dynamics1Level<Integration1stHalfWithWallforPD>;

  public:
    explicit PDFluidPressureHalf(PDFluidAcousticWorkspace &workspace);

  protected:
    void runMainStep(Real dt) override;

  private:
    PDFluidAcousticWorkspace &workspace_;
};

class PDFluidDensityHalf : public Dynamics1Level<Integration2ndHalfWithWallRiemann>
{
    using Base = Dynamics1Level<Integration2ndHalfWithWallRiemann>;

  public:
    explicit PDFluidDensityHalf(PDFluidAcousticWorkspace &workspace);

  protected:
    void runMainStep(Real dt) override;

  private:
    PDFluidAcousticWorkspace &workspace_;
};

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
    friend class PDFluidAcousticDamping;
    void applyDamping(std::size_t i, const Vecd &sum, Real dt);

    Real coeff_acoustic_damper_ = 0.0;
    Real *Vol_, *rho_, *u_div_;
    Real *mass_, *work_;
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
    friend class PDFluidAcousticDamping;
    void applyDamping(std::size_t i, const Vecd &sum, Real dt);

    Real coeff_acoustic_damper_ = 0.0;
    Real *Vol_, *rho_, *u_div_;
    Real *mass_, *work_;
    Vecd *force_, *vel_;
    StdVec<Real *> contact_Vol_;
};
class PDFluidVelocityDivergence
{
  public:
    explicit PDFluidVelocityDivergence(PDFluidAcousticWorkspace &workspace);
    void exec(Real dt = 0.0);

  private:
    PDFluidAcousticWorkspace &workspace_;
    InteractionDynamics<VelocityDivergenceInner> inner_;
    InteractionDynamics<VelocityDivergenceWithWallforPD> wall_;
};

class PDFluidAcousticDamping
{
  public:
    explicit PDFluidAcousticDamping(PDFluidAcousticWorkspace &workspace);
    void setCoeffAcousticDamper(Real rho0, Real c0, Real h);
    void exec(Real dt = 0.0);

  private:
    PDFluidAcousticWorkspace &workspace_;
    InteractionDynamics<AcousticDamperInnerforPD> inner_;
    InteractionDynamics<AcousticDamperWithWallforPD> wall_;
};
} // namespace fluid_dynamics

namespace solid_dynamics
{
class ViscousForceFromFluidWithExecution : public ViscousForceFromFluid
{
  public:
    using ViscousForceFromFluid::ViscousForceFromFluid;
    void evaluateDevice(peridynamics::FluidViscosityBatch &backend,
                        peridynamics::SPHViscosityWorkspace &workspace);
};

class ViscousForceFromFluidforPD : public ViscousForceFromFluidWithExecution
{
  public:
    explicit ViscousForceFromFluidforPD(BaseContactRelation &relation)
        : ViscousForceFromFluidWithExecution(relation)
    {
        vel_ave_ = particles_->registerStateVariableData<Vecd>("Velocity");
    }
};

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
    explicit PressureForceFromFluidforPD(
        BaseContactRelation &contact_relation, const peridynamics::FluidNeighborConfiguration *configuration = nullptr);
    virtual ~PressureForceFromFluidforPD() {};

    void setCoeffAcousticDamper(Real rho0, Real c0, Real h)
    {
        coeff_acoustic_damper_ = 0.3 * c0 * rho0 * h;
    };

    void setAccelerationPressureFactor(Real value)
    {
        if (!std::isfinite(value) || value < 0 || value > 1) throw std::invalid_argument("invalid pressure acceleration factor");
        acceleration_factor_=value;
    }
    void interaction(size_t index_i, Real dt = 0.0);
    void evaluateDevice(peridynamics::FluidTractionBatch &backend);

  protected:
    void publishTraction(size_t index, const Vecd &pressure_force, const Vecd &acoustic_force,
                         Real pressure_sum, Real kernel_weight);
    HughesWingetSolid &hw_solid_;
    Real coeff_acoustic_damper_ = 0.0;
    Real acceleration_factor_ = 1;
    Real *Vol_;
    Vecd *vel_ave_, *acc_prior_, *acc_ave_, *n_;
    Vecd *force_from_fluid_, *acoustic_force_from_fluid_;
    Real *pressure_from_fluid_;
    StdVec<Fluid *> contact_fluids_;
    StdVec<Real *> contact_rho_, contact_p_, contact_u_div_, contact_mass_, contact_Vol_;
    StdVec<Vecd *> contact_vel_, contact_force_prior_;
    StdVec<RiemannSolverType> riemann_solvers_;

  private:
    void configureDevice(peridynamics::FluidTractionBatch &backend, std::size_t count);
    peridynamics::FluidNeighborCache neighbor_cache_;
    peridynamics::FluidTractionBatch *configured_backend_ = nullptr;
    StdVec<std::pair<std::size_t, std::size_t>> fluid_sources_;
    StdVec<peridynamics::FluidTractionSolid> solid_rows_;
    StdVec<peridynamics::FluidTractionFluid> fluid_sites_;
    StdVec<peridynamics::FluidTractionNeighbor> traction_neighbors_;
    StdVec<peridynamics::FluidTractionResult> traction_result_;
    StdVec<size_t> row_offsets_;
    StdVec<StdVec<size_t>> fluid_indices_;
    StdVec<std::pair<Real,Real>> acoustic_coefficients_;
};

template<class RiemannSolverType = NoRiemannSolver>
class PDFluidPressureCoupling : public InteractionDynamics<PressureForceFromFluidforPD<RiemannSolverType>>
{
  public:
    using Base=InteractionDynamics<PressureForceFromFluidforPD<RiemannSolverType>>;
    PDFluidPressureCoupling(BaseContactRelation &relation, const peridynamics::PeriExecution &execution,
                            const peridynamics::FluidNeighborConfiguration *configuration = nullptr)
        : Base(relation, configuration), backend_(execution.makeFluidTraction()) {}

    bool onDevice() const noexcept { return backend_!=nullptr; }

    void runMainStep(Real dt) override
    {
        if (backend_)
            this->evaluateDevice(*backend_);
        else
            Base::runMainStep(dt);
    }

  private:
    std::unique_ptr<peridynamics::FluidTractionBatch> backend_;
};

using PressureForceAccelerationFromFluidforPD = PressureForceFromFluidforPD<NoRiemannSolver>;
} // namespace solid_dynamics
} // namespace SPH
#endif // PD_FSI_H
