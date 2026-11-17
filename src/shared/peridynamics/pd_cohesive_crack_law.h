#pragma once
/* Licensed under the Apache License, Version 2.0. */
#include "data_type.h"
#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace SPH::peridynamics
{
/** Irreversible bilinear cohesive traction-separation law, per reference area.
 *
 * This is the equal-mode-strength/equal-mode-toughness specialization of the
 * standard scalar-damage bilinear law of Camanho & Davila (NASA/TM-2002-211737)
 * and Turon et al. (NASA/TM-2004-213277, damage evolution Eq. 33). The absence of
 * separate mode-II/III data is an EXPLICIT mode-independent assumption: the
 * same peak effective traction ft, fracture work Gc and numerical stiffness K
 * apply to positive normal opening and tangential separation. This is not an
 * experimentally identified mixed-mode interaction law or friction model.
 *
 * jump=x_plus-x_minus, and the unit normal points from minus to plus.
 * effective_opening^2=|jump_t|^2+<jump_n>_+^2.
 * Result::traction is d(energy)/d(jump) with damage and normal FROZEN, in Pa.
 * Actual physical forces are -A_ref*traction on plus and its opposite on minus.
 * Negative normal opening has no energy/traction here: a separate unilateral
 * contact law must handle compression, without duplicating normal penalties.
 *
 * The templated storedPotential supports outer geometric AD. If the face normal
 * depends on the moving material vertices, the OUTER assembly must include its
 * derivative; this local delta derivative does not replace the geometric
 * derivative. The damage history update itself is not differentiated as a
 * recoverable potential: its release is accounted in dissipation_increment.
 *
 * This law supplies the only Gc payment for its cohesive surface. Do not pay an
 * old point-level crack-band Gc again, or add another Gc on final bond removal.
 * Topological insertion/transfer of pre-existing bulk elastic energy requires
 * its own consistent assembly and energy check; a valid local law does not
 * make a discontinuous topology switch automatically energy conserving.
 * No viscosity, residual strength, force clipping, or material rescaling.
 */
class BilinearCohesiveCrackLaw
{
  public:
    struct History
    {Real kappa=0,damage=0,dissipation=0;}; // m, dimensionless, J/m^2
    struct Result
    {
        Vecd traction=Vecd::Zero(),effective_jump=Vecd::Zero();
        Matd frozen_tangent=Matd::Zero(),algorithmic_tangent=Matd::Zero();
        Real energy=0,dissipation=0,dissipation_increment=0;
        Real normal_opening=0,effective_opening=0,kappa=0,damage=0;
        Real damage_derivative=0,frozen_tangent_bound=0,algorithmic_tangent_bound=0;
        bool loading_softening=false,fully_failed=false;
    };
    const Real strength,Gc,stiffness,delta0,delta_f;
    BilinearCohesiveCrackLaw(Real tensile_strength,Real fracture_energy,Real numerical_stiffness)
      :strength(tensile_strength),Gc(fracture_energy),stiffness(numerical_stiffness),
       delta0(tensile_strength/numerical_stiffness),delta_f(2*fracture_energy/tensile_strength)
    {
        static_assert(Dimensions==3,"Cohesive crack law uses a three-dimensional jump");
        if(!std::isfinite(strength)||!std::isfinite(Gc)||!std::isfinite(stiffness)||
           std::min({strength,Gc,stiffness})<=0||!std::isfinite(delta0)||!std::isfinite(delta_f))
            throw std::invalid_argument("Cohesive law requires positive finite ft, Gc and K");
        if(delta_f<=delta0*(1+128*std::numeric_limits<Real>::epsilon()))
            throw std::invalid_argument("Numerical cohesive K cannot resolve the prescribed ft/Gc bilinear law; do not clip material inputs");
    }
    History historyForKappa(Real kappa)const
    {
        if(!std::isfinite(kappa)||kappa<0)throw std::invalid_argument("Invalid cohesive history opening");
        History h;h.kappa=kappa;h.damage=damageAt(kappa);h.dissipation=dissipationAt(kappa);return h;
    }
    /** Existing precracks already consumed Gc before this analysis. The returned
     * historical Gc is an initial budget, not a new dissipation increment; the
     * driver must subtract its initial history from current fracture work. */
    History precrackedHistory()const{return historyForKappa(delta_f);}

    Result evaluate(const Vecd&jump,const Vecd&normal,History&history,bool advance=true)const
    {
        if(!jump.allFinite()||!normal.allFinite()||std::abs(normal.norm()-1)>1e-8)
            throw std::invalid_argument("Cohesive jump must be finite and normal must be unit length");
        validateHistory(history);Vecd n=normal.normalized();Result r;
        r.normal_opening=jump.dot(n);Matd projection=Matd::Identity();
        // At exactly zero normal opening the second derivative is not unique.
        // Choose the tensile one-sided generalized tangent; traction is C0 and
        // stored energy is C1 across this switch at fixed damage.
        if(r.normal_opening>=0)r.effective_jump=jump;
        else{projection-=n*n.transpose();r.effective_jump=jump-r.normal_opening*n;}
        r.effective_opening=r.effective_jump.stableNorm();
        if(!std::isfinite(r.effective_opening))throw std::runtime_error("Nonfinite cohesive effective opening");
        Real old_dissipation=dissipationAt(history.kappa);
        Real current_kappa=advance?std::max(history.kappa,r.effective_opening):history.kappa;
        History current=historyForKappa(current_kappa);
        r.kappa=current.kappa;r.damage=current.damage;r.dissipation=current.dissipation;
        r.dissipation_increment=advance?current.dissipation-old_dissipation:0;
        if(r.dissipation_increment<0)throw std::runtime_error("Cohesive history would heal");
        r.fully_failed=current.damage==1;
        if(!r.fully_failed)
        {
            Real retained=(1-current.damage)*stiffness;
            r.energy=.5*retained*r.effective_jump.squaredNorm();
            r.traction=retained*r.effective_jump;r.frozen_tangent=retained*projection;
            r.algorithmic_tangent=r.frozen_tangent;
            r.loading_softening=advance&&r.effective_opening>history.kappa&&current.kappa>delta0&&current.kappa<delta_f;
            if(r.loading_softening)
            {
                r.damage_derivative=delta_f*delta0/(current.kappa*current.kappa*(delta_f-delta0));
                r.algorithmic_tangent-=stiffness*r.damage_derivative/r.effective_opening*
                    (r.effective_jump*r.effective_jump.transpose());
            }
            r.frozen_tangent_bound=absoluteRowBound(r.frozen_tangent);
            r.algorithmic_tangent_bound=absoluteRowBound(r.algorithmic_tangent);
        }
        if(!std::isfinite(r.energy)||!r.traction.allFinite()||!r.frozen_tangent.allFinite()||!r.algorithmic_tangent.allFinite())
            throw std::runtime_error("Nonfinite cohesive traction/energy/tangent");
        if(advance)history=current;
        return r;
    }

    /** Recoverable surface potential with frozen history and active branch.
     * Scalar may be an outer geometry-AD type. n must be its CURRENT normalized
     * facet normal; differentiate n along with both surface traces. The caller
     * chooses normal_open from the scalar value of jump.n. This has no contact
     * compression term and no extra fracture-work term. */
    template<class Scalar>Scalar storedPotential(const std::array<Scalar,3>&jump,
        const std::array<Scalar,3>&normal,Real damage,bool normal_open)const
    {
        if(!std::isfinite(damage)||damage<0||damage>1)throw std::invalid_argument("Invalid frozen cohesive damage");
        if(damage==1)return Scalar(0);
        Scalar gn(0),equivalent2(0);for(size_t d=0;d<3;++d)gn=gn+jump[d]*normal[d];
        for(size_t d=0;d<3;++d){Scalar tangent=jump[d]-gn*normal[d];equivalent2=equivalent2+tangent*tangent;}
        if(normal_open)equivalent2=equivalent2+gn*gn;
        return Scalar(.5*(1-damage)*stiffness)*equivalent2;
    }
    Real storedPotential(const Vecd&jump,const Vecd&normal,Real damage)const
    {
        if(!jump.allFinite()||!normal.allFinite()||std::abs(normal.norm()-1)>1e-8)
            throw std::invalid_argument("Invalid frozen cohesive potential state");
        Vecd n=normal.normalized();std::array<Real,3>a{},b{};
        for(size_t d=0;d<3;++d){a[d]=jump[d];b[d]=n[d];}
        return storedPotential(a,b,damage,jump.dot(n)>=0);
    }
  private:
    Real damageAt(Real kappa)const
    {
        if(kappa<=delta0)return 0;if(kappa>=delta_f)return 1;
        return delta_f*(kappa-delta0)/(kappa*(delta_f-delta0));
    }
    /** Integral of Y dd on a loading envelope. dD/dkappa is constant for
     * bilinear softening, hence the cumulative work reaches exactly Gc. */
    Real dissipationAt(Real kappa)const
    {return Gc*(std::clamp(kappa,delta0,delta_f)-delta0)/(delta_f-delta0);}
    static Real absoluteRowBound(const Matd&matrix)
    {Real bound=0;for(size_t i=0;i<3;++i){Real row=0;for(size_t j=0;j<3;++j)row+=std::abs(matrix(i,j));bound=std::max(bound,row);}return bound;}
    void validateHistory(const History&h)const
    {
        if(!std::isfinite(h.kappa)||!std::isfinite(h.damage)||!std::isfinite(h.dissipation)||h.kappa<0||h.damage<0||h.damage>1||h.dissipation<0)
            throw std::invalid_argument("Invalid irreversible cohesive history");
        Real tolerance=256*std::numeric_limits<Real>::epsilon();
        if(std::abs(h.damage-damageAt(h.kappa))>tolerance||
           std::abs(h.dissipation-dissipationAt(h.kappa))>tolerance*std::max(Gc,Real(1)))
            throw std::invalid_argument("Cohesive derived history is inconsistent with its maximum opening");
    }
};
}
