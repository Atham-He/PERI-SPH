#pragma once
#include "solid_quadratic_reference.h"

namespace SPH::solid_dynamics
{
/** Keep the native PK2 material/time update and its native SPH force.
 * Prepare the residual energy after all midpoint drifts; apply its exact
 * derivative only after every reference row has been prepared. This split
 * matches Dynamics1Level and the example's explicit three-loop midpoint path.
 */
class SolidResidualIntegration1stHalf : public Integration1stHalfPK2RightCauchy
{
    SolidQuadraticReference &reference_;
    Real *residual_energy_;
    Vecd *residual_force_;
  public:
    SolidResidualIntegration1stHalf(BaseInnerRelation &relation,SolidQuadraticReference &reference)
      :Integration1stHalfPK2RightCauchy(relation),reference_(reference),
       residual_energy_(particles_->registerStateVariableData<Real>("SolidResidualEnergyDensity")),
       residual_force_(particles_->registerStateVariableData<Vecd>("SolidResidualForce")){}
    void interaction(size_t i,Real dt=0)
    {reference_.prepare(i,pos_);residual_energy_[i]=reference_.energy[i]/reference_.volume[i];}
    void update(size_t i,Real dt=0)
    {
        Integration1stHalfPK2RightCauchy::interaction(i,dt);
        residual_force_[i]=reference_.stabilizingForce(i);force_[i]+=residual_force_[i];
        BaseIntegration1stHalf::update(i,dt);
    }
};

} // namespace SPH::solid_dynamics
