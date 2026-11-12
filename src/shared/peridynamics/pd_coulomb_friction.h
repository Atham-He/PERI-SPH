/* Licensed under the Apache License, Version 2.0. */
#ifndef PD_COULOMB_FRICTION_H
#define PD_COULOMB_FRICTION_H

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <stdexcept>

namespace SPH::peridynamics
{
/** Rate-independent elastic tangential slip with a Coulomb return map.
 * Coordinates must already be projected into a fixed contact tangent frame.
 * The caller supplies a positive compressive normal force, not a pressure.
 * Trial evaluation is const; only commit advances the slip history. */
template<class Vector>
class CoulombTangentialHistory
{
  public:
    using Scalar=typename Vector::Scalar;
    struct Result
    {
        Vector force=Vector::Zero(),anchor=Vector::Zero(),elastic_slip=Vector::Zero();
        Scalar stored=0,slip_dissipation=0,release_energy=0;
        int side=0;
        bool sliding=false;
        size_t base_revision=0;
    };
    Vector anchor=Vector::Zero(),elastic_slip=Vector::Zero();
    Scalar stored=0,slip_dissipation=0,release_energy=0;
    int side=0;
    size_t revision=0;

    Result trial(const Vector &coordinate,int contact_side,Scalar normal_force,
                 Scalar tangent_stiffness,Scalar friction) const
    {
        if(!coordinate.allFinite()||!std::isfinite(normal_force)||normal_force<0||
           !std::isfinite(tangent_stiffness)||tangent_stiffness<=0||
           !std::isfinite(friction)||friction<0||contact_side<-1||contact_side>1)
            throw std::invalid_argument("Invalid Coulomb contact state or input");
        Result r;r.base_revision=revision;r.anchor=coordinate;
        if(contact_side==0||normal_force==0||friction==0)
        {r.release_energy=stored;return r;}
        r.side=contact_side;
        if(side!=contact_side)
        {r.release_energy=stored;return r;}
        r.anchor=anchor;r.elastic_slip=coordinate-anchor;
        Scalar trial_norm=tangent_stiffness*r.elastic_slip.norm();
        Scalar limit=friction*normal_force;
        if(trial_norm>limit)
        {
            Vector trial_slip=r.elastic_slip;
            r.elastic_slip*=limit/trial_norm;
            r.anchor=coordinate-r.elastic_slip;r.sliding=true;
            r.slip_dissipation=limit*(trial_slip-r.elastic_slip).norm();
        }
        r.force=-tangent_stiffness*r.elastic_slip;
        r.stored=Scalar(.5)*tangent_stiffness*r.elastic_slip.squaredNorm();
        return r;
    }
    void commit(const Result &r)
    {
        if(r.base_revision!=revision)throw std::logic_error("Stale or duplicate Coulomb history commit");
        anchor=r.anchor;elastic_slip=r.elastic_slip;stored=r.stored;side=r.side;
        slip_dissipation+=r.slip_dissipation;release_energy+=r.release_energy;++revision;
    }
};
} // namespace SPH::peridynamics
#endif
