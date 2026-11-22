/* Licensed under the Apache License, Version 2.0. */
#ifndef PD_COULOMB_FRICTION_H
#define PD_COULOMB_FRICTION_H

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <stdexcept>

namespace SPH::peridynamics
{

#if defined(__CUDACC__)
#define SPH_COULOMB_HD __host__ __device__
#else
#define SPH_COULOMB_HD
#endif

template<class Scalar>
struct CoulombSlipUpdate
{
    Scalar stored=0, slip_dissipation=0, release_energy=0;
    int side=0;
    bool sliding=false;
};

// The fixed contact tangent frame is shared by the host history and device trial.
template<class Scalar>
SPH_COULOMB_HD CoulombSlipUpdate<Scalar> evaluateCoulombSlip(
    const Scalar *coordinate, const Scalar *old_anchor, size_t dimensions,
    int old_side, Scalar old_stored, int side, Scalar normal,
    Scalar stiffness, Scalar friction, Scalar *anchor, Scalar *slip)
{
    CoulombSlipUpdate<Scalar> result;
    for(size_t d=0;d<dimensions;++d){
        anchor[d]=coordinate[d];
        slip[d]=0;
    }
    if(side==0||normal==0||friction==0){
        result.release_energy=old_stored;
        return result;
    }
    result.side=side;
    if(old_side!=side){
        result.release_energy=old_stored;
        return result;
    }
    Scalar norm_squared=0;
    for(size_t d=0;d<dimensions;++d){
        anchor[d]=old_anchor[d];
        slip[d]=coordinate[d]-old_anchor[d];
        norm_squared+=slip[d]*slip[d];
    }
    const Scalar trial_norm=stiffness*std::sqrt(norm_squared);
    const Scalar limit=friction*normal;
    if(trial_norm>limit){
        Scalar plastic_norm_squared=0;
        norm_squared=0;
        for(size_t d=0;d<dimensions;++d){
            const Scalar trial=slip[d];
            slip[d]*=limit/trial_norm;
            anchor[d]=coordinate[d]-slip[d];
            const Scalar plastic=trial-slip[d];
            plastic_norm_squared+=plastic*plastic;
            norm_squared+=slip[d]*slip[d];
        }
        result.sliding=true;
        result.slip_dissipation=limit*std::sqrt(plastic_norm_squared);
    }
    result.stored=Scalar(.5)*stiffness*norm_squared;
    return result;
}
#undef SPH_COULOMB_HD
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
        Result r;
        r.base_revision=revision;
        const auto update=evaluateCoulombSlip(
            coordinate.data(),anchor.data(),coordinate.size(),side,stored,contact_side,normal_force,
            tangent_stiffness,friction,r.anchor.data(),r.elastic_slip.data());
        r.force=-tangent_stiffness*r.elastic_slip;
        r.stored=update.stored;
        r.slip_dissipation=update.slip_dissipation;
        r.release_energy=update.release_energy;
        r.side=update.side;
        r.sliding=update.sliding;
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
