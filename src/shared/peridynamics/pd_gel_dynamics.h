#ifndef PD_GEL_DYNAMICS_H
#define PD_GEL_DYNAMICS_H
#include "pd_gel_material.h"
#include "base_local_dynamics.h"
#include "base_body_relation.h"
#include "base_particle_dynamics.h"

namespace SPH::solid_dynamics
{
/** Midpoint force prediction never commits accepted history. Commit reconstructs
 * accepted geometry after the final drift, integrates the second half segment,
 * and leaves last-applied force/stress fields distinct from accepted stress. */
class GelConstitutiveUpdate : public LocalDynamics, public DataDelegateInner
{
  public:
    enum class Stage { MidpointForce, AcceptedHistory };
    GelConstitutiveUpdate(BaseInnerRelation &relation, Stage stage);
    void interaction(size_t i, Real dt = 0);

  private:
    Stage stage_;
    GelViscoelasticSolid &material_;
    Vecd *position_;
    Real *volume_, *density_, *dissipation_, *midpoint_dissipation_, *energy_, *history_time_;
    Matd *inverse_shape_, *deformation_, *stress_, *pk1_, *force_state_, *accepted_stress_;
    Matd *green_, *viscous_, *midpoint_green_, *midpoint_viscous_;
};
} // namespace SPH::solid_dynamics
#endif
