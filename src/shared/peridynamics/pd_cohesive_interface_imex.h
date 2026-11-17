#pragma once
/* Licensed under the Apache License, Version 2.0. */
#include "data_type.h"
#include <algorithm>
#include <cmath>
#include <limits>
#include <map>
#include <sstream>
#include <stdexcept>
#include <vector>

namespace SPH::peridynamics
{
/** Coupled fixed-trace realization of the monotone/softening split of Scheme B,
 * Doyen, Ern & Piperno, Comput. Mech. 52 (2013) 401-416,
 * doi:10.1007/s00466-012-0819-2, Eqs. 42,46-48.
 *
 * This is an initially rigid EXTRINSIC equal-mode law, not the finite-initial-K
 * BilinearCohesiveCrackLaw. dc=2Gc/ft; kappa is the maximum physical effective
 * jump. P_kappa(delta)=ft*delta/max(kappa,|delta|) and
 * Theta(delta)=-ft*min(|delta|/dc,1)*delta/|delta|. At kappa=delta=0,
 * P is the closed traction ball, so finite insertion traction does not require
 * an artificial opening, an initial elastic-energy offset, or modified masses.
 * An optional softening_evaluation is an outer iterate, allowing a converged
 * fully implicit softening solve; the default remains Scheme B.
 * The caller must independently verify discrete bulk-force/energy transfer at
 * topology insertion. This solver does NOT make an inconsistent insertion valid.
 *
 * B is a FIXED signed trace Jacobian (scalar MLS weights plus optional full
 * geometric blocks); row.projector optionally
 * projects that trace (e.g. onto a closed face's tangent plane). Normals, trace
 * support, projection active-set and current geometry are frozen in this solve.
 * A changing normal/projector requires the outer objective geometric derivative
 * and a converged active-set/geometric treatment. Raw bilateral B treats all
 * directions as cohesive: do not pass compressive normal jumps while claiming
 * that compression is handled only by separate crack-face contact. A closed
 * tangent projector removes cohesive compression on a known fixed branch.
 * Contact itself is not solved here; no contact penalty or friction is added.
 *
 * qfree must be the bulk-only central-difference predictor, WITHOUT the old or
 * new cohesive correction. Theta(old jump) is added here exactly once. With
 * z=sqrt(A)*t_mono and theta=sqrt(A)*Theta(old jump), solve the convex dual
 *   min 1/2(z+theta)^T D(z+theta) - z.s + sum kappa/(2ft)|z_q|^2,
 *   |z_q|<=ft*sqrt(A_q), z_q in range(projector_q),
 * where s=sqrt(A)*B*qfree (including fixed jump offsets),
 * D=dt^2 sqrt(A) B M^{-1} B^T sqrt(A).
 * Real free-DOF inverse masses are used; zero inverse mass means prescribed DOF.
 * Shared MLS nodes/DOFs couple ALL facets. No facetwise effective-mass shortcut.
 *
 * The matrix-free row-sum Lipschitz bound is conservative, including signed MLS
 * weights, aliasing, projections, and shared DOFs. Projected-gradient/FISTA is
 * iterated to a reported projected KKT OPENING residual; nonconvergence throws.
 * Neither the converged force nor a failed solve is clipped to hide an error.
 *
 * Result provides BOTH the actual time-discrete split traction and the physical
 * end-state law traction, and their difference. If a step first crosses dc,
 * its new history is failed and end-state law traction is zero, but the lagged
 * split traction can persist for that time step; it is omitted on the next
 * call. Removing that force after computing q is inconsistent: exact release
 * timing requires an outer event-localized substep/re-solve. The split is not
 * a claim of exact finite-dt energy conservation. The RIGHT-ENDPOINT surface-work diagnostic minus Delta(psi+D) is
 * reported with its sign (not the paper's centered discrete-energy balance); the outer full energy audit must retain it alongside
 * bulk, kinetic, contact, boundary work and numerical damping. Gc is paid once.
 * A pre-existing crack's initial historical dissipation is not newly paid work.
 */
class CohesiveInterfaceIMEX
{
  public:
    struct Entry {
        size_t node=0;Real weight=0;
        Matd geometric=Matd::Zero();
        Matd block()const{return weight*Matd::Identity()+geometric;}
    };
    struct Facet
    {
        std::vector<Entry> entries;
        Real reference_area=0;
        Matd projector=Matd::Identity();
        Vecd jump_offset=Vecd::Zero(); // projected fixed affine offset, in m
        Vecd traction_halfspace_normal=Vecd::Zero(); // optional t.n >= 0, within projector range
    };
    struct History {Real kappa=0;};
    struct Options
    {
        size_t max_iterations=20000;
        Real relative_tolerance=1e-10;
        Real absolute_opening_tolerance=1e-14;
        bool accelerate=true;
    };
    struct Result
    {
        std::vector<Vecd> positions,correction,nodal_force,jump;
        std::vector<Vecd> dual,monotone_traction,softening_traction,traction,law_traction;
        std::vector<History> history;
        Real lipschitz_bound=0,projected_opening_residual=0,relative_residual=0;
        Real maximum_split_traction_error=0,maximum_dual_constraint_error=0,maximum_dual_range_error=0;
        Real stored_energy=0,previous_stored_energy=0;
        Real dissipation=0,dissipation_increment=0;
        Real interface_work=0,algorithmic_work_residual=0,dual_objective=0; // right-endpoint diagnostic, NOT an exact energy balance
        Real force_position_residual=0,fenchel_gap=0;
        size_t iterations=0,active_facets=0,indeterminate_prescribed_adhesion_facets=0;
        bool converged=false;
    };
    const Real strength,Gc,delta_f;

    CohesiveInterfaceIMEX(Real ft,Real fracture_energy,
        const std::vector<Vecd>&inverse_mass,const std::vector<Facet>&facets)
      :CohesiveInterfaceIMEX(ft,fracture_energy,inverse_mass,facets,Options{}){}
    CohesiveInterfaceIMEX(Real ft,Real fracture_energy,
        const std::vector<Vecd>&inverse_mass,const std::vector<Facet>&facets,const Options&options)
      :strength(ft),Gc(fracture_energy),delta_f(2*fracture_energy/ft),
       inverse_mass_(inverse_mass),facets_(facets),options_(options)
    {
        static_assert(Dimensions==3,"Cohesive interface solve requires 3D vector traces");
        if(!std::isfinite(ft)||!std::isfinite(Gc)||!std::isfinite(delta_f)||delta_f<=0||std::min(ft,Gc)<=0)
            throw std::invalid_argument("IMEX cohesive ft/Gc must be positive and finite");
        if(options_.max_iterations==0||!std::isfinite(options_.relative_tolerance)||
           !std::isfinite(options_.absolute_opening_tolerance)||options_.relative_tolerance<=0||
           options_.absolute_opening_tolerance<=0)
            throw std::invalid_argument("Invalid cohesive IMEX convergence options");
        for(const Vecd&m:inverse_mass_)if(!m.allFinite()||m.minCoeff()<0)
            throw std::invalid_argument("Cohesive inverse masses must be finite and nonnegative per DOF");
        for(Facet&f:facets_)
        {
            if(!std::isfinite(f.reference_area)||f.reference_area<=0||!f.projector.allFinite()||
               !f.jump_offset.allFinite()||(f.projector-f.projector.transpose()).norm()>1e-10||
               (f.projector*f.projector-f.projector).norm()>1e-10||
               (f.projector*f.jump_offset-f.jump_offset).norm()>1e-12)
                throw std::invalid_argument("Cohesive facet requires positive reference area and an orthogonal fixed projector");
            // The tolerance admits only round-off in supplied projected data.
            // Canonicalize it: a tolerated non-range offset must never become
            // a physical opening or an adhesion reaction in a zero-rank row.
            if(!f.traction_halfspace_normal.allFinite()||(f.traction_halfspace_normal.norm()!=0&&std::abs(f.traction_halfspace_normal.norm()-1)>1e-10))throw std::invalid_argument("Invalid unilateral cohesive dual normal");
            f.jump_offset=f.projector*f.jump_offset;
            std::map<size_t,Matd>merged;
            for(const Entry&e:f.entries)
            {
                if(e.node>=inverse_mass_.size()||!std::isfinite(e.weight)||!e.geometric.allFinite())throw std::invalid_argument("Invalid cohesive trace Jacobian entry");
                auto it=merged.try_emplace(e.node,Matd::Zero()).first;it->second+=e.block();
            }
            f.entries.clear();Matd translation=Matd::Zero();Real absolute_sum=0;
            for(const auto&e:merged)if(e.second.norm()!=0){f.entries.push_back({e.first,0,e.second});translation+=e.second;absolute_sum+=e.second.norm();}
            if(f.entries.empty()||translation.norm()>1e-10*std::max(absolute_sum,Real(1)))
                throw std::invalid_argument("Signed crack Jacobian must reproduce common rigid translation");
        }
    }

    std::vector<Vecd> applyB(const std::vector<Vecd>&q,bool include_offsets=false)const
    {
        requireNodes(q);std::vector<Vecd>out(facets_.size(),Vecd::Zero());
        for(size_t k=0;k<facets_.size();++k)
        {const Facet&f=facets_[k];for(const Entry&e:f.entries)out[k]+=e.block()*q[e.node];out[k]=f.projector*out[k];if(include_offsets)out[k]+=f.jump_offset;}
        return out;
    }
    std::vector<Vecd> applyBT(const std::vector<Vecd>&v)const
    {
        if(v.size()!=facets_.size())throw std::invalid_argument("Wrong cohesive transpose facet count");
        std::vector<Vecd>out(inverse_mass_.size(),Vecd::Zero());
        for(size_t k=0;k<facets_.size();++k)
        {if(!v[k].allFinite())throw std::invalid_argument("Nonfinite cohesive transpose input");Vecd value=facets_[k].projector.transpose()*v[k];for(const Entry&e:facets_[k].entries)out[e.node]+=e.block().transpose()*value;}
        return out;
    }

    Result solve(const std::vector<Vecd>&qfree,const std::vector<Vecd>&old_jump,
        const std::vector<History>&old_history,Real dt,
        const std::vector<Vecd>&warm_dual={},const std::vector<Vecd>&softening_evaluation={})const
    {
        requireNodes(qfree);const size_t count=facets_.size();
        if(old_jump.size()!=count||old_history.size()!=count||(!warm_dual.empty()&&warm_dual.size()!=count)||
           (!softening_evaluation.empty()&&softening_evaluation.size()!=count)||!std::isfinite(dt)||dt<=0||!std::isfinite(dt*dt))
            throw std::invalid_argument("Invalid cohesive solve state/time step");
        Result result;result.positions=qfree;result.jump=applyB(qfree,true);result.history=old_history;
        result.dual.assign(count,Vecd::Zero());result.softening_traction.assign(count,Vecd::Zero());
        std::vector<Vecd>s(count,Vecd::Zero()),theta(count,Vecd::Zero()),projected_old_jump(old_jump);
        std::vector<Real>sqrt_area(count),compliance(count);std::vector<bool>active(count,false);
        Real opening_scale=delta_f;
        for(size_t k=0;k<count;++k)
        {
            const Real kap=old_history[k].kappa;const Facet&f=facets_[k];sqrt_area[k]=std::sqrt(f.reference_area);
            if(!std::isfinite(kap)||kap<0||!old_jump[k].allFinite()||
               (f.projector*old_jump[k]-old_jump[k]).norm()>1e-10*std::max(delta_f,old_jump[k].norm()))
                throw std::invalid_argument("Invalid IMEX history/projected old opening");
            // Use the same projected old jump for Theta, maximum-opening
            // validation, stored energy and work; numerical range tolerance
            // must not create cohesive compression or extra fracture history.
            projected_old_jump[k]=f.projector*old_jump[k];
            Real old_r=projected_old_jump[k].stableNorm();
            if(kap<delta_f&&old_r>kap+256*std::numeric_limits<Real>::epsilon()*std::max(delta_f,kap))
                throw std::invalid_argument("Old cohesive maximum opening is smaller than its old jump");
            if(kap>=delta_f)continue;active[k]=true;++result.active_facets;
            compliance[k]=kap/strength;s[k]=sqrt_area[k]*result.jump[k];
            Vecd theta_jump=projected_old_jump[k];
            if(!softening_evaluation.empty()){
                theta_jump=softening_evaluation[k];
                if(!theta_jump.allFinite()||(f.projector*theta_jump-theta_jump).norm()>1e-10*std::max(delta_f,theta_jump.norm()))throw std::invalid_argument("Invalid outer cohesive softening evaluation");
                theta_jump=f.projector*theta_jump;
            }
            theta_jump=effectiveOpening(theta_jump,k);Real theta_radius=theta_jump.stableNorm();
            if(theta_radius>0)result.softening_traction[k]=-strength*std::min(theta_radius/delta_f,Real(1))*theta_jump/theta_radius;
            theta[k]=sqrt_area[k]*result.softening_traction[k];
            if(!warm_dual.empty()){if(!warm_dual[k].allFinite())throw std::invalid_argument("Nonfinite warm cohesive dual");result.dual[k]=project(warm_dual[k],k,sqrt_area[k]);}
            opening_scale=std::max({opening_scale,result.jump[k].stableNorm(),old_r,kap});
        }
        // A matrix-free Gershgorin upper bound, without cancellation of absolute
        // contributions. No dense DOF or facet matrix is allocated.
        std::vector<Real>metric=lipschitzBounds(active,sqrt_area,compliance,dt);result.lipschitz_bound=0;for(Real value:metric)result.lipschitz_bound=std::max(result.lipschitz_bound,value);
        auto multiplyD=[&](const std::vector<Vecd>&v)
        {
            std::vector<Vecd>weighted(count,Vecd::Zero());for(size_t k=0;k<count;++k)if(active[k])weighted[k]=sqrt_area[k]*v[k];
            std::vector<Vecd>nodes=applyBT(weighted);for(size_t i=0;i<nodes.size();++i)nodes[i]=dt*dt*inverse_mass_[i].cwiseProduct(nodes[i]);
            std::vector<Vecd>out=applyB(nodes);for(size_t k=0;k<count;++k){if(active[k])out[k]*=sqrt_area[k];else out[k].setZero();}return out;
        };
        auto gradient=[&](const std::vector<Vecd>&z)
        {std::vector<Vecd>v(count);for(size_t k=0;k<count;++k)v[k]=z[k]+theta[k];std::vector<Vecd>g=multiplyD(v);for(size_t k=0;k<count;++k)if(active[k])g[k]+=compliance[k]*z[k]-s[k];return g;};
        auto objective=[&](const std::vector<Vecd>&z)
        {std::vector<Vecd>v(count);for(size_t k=0;k<count;++k)v[k]=z[k]+theta[k];auto Dv=multiplyD(v);Real value=0;for(size_t k=0;k<count;++k)if(active[k])value+=.5*v[k].dot(Dv[k])-z[k].dot(s[k])+.5*compliance[k]*z[k].squaredNorm();return value;};
        const Real L=result.lipschitz_bound,tolerance=options_.absolute_opening_tolerance+options_.relative_tolerance*opening_scale;
        if(L==0)
        {
            // All active traces are prescribed and kappa=0. Displacements are
            // unique; zero-jump traction may be nonunique, so keep feasible warm
            // values there rather than invent an adhesion reaction.
            for(size_t k=0;k<count;++k)if(active[k])
            {Vecd direction=effectiveOpening(facets_[k].projector*s[k],k);Real norm=direction.stableNorm();
             if(norm>0)result.dual[k]=project(strength*sqrt_area[k]*(direction/norm),k,sqrt_area[k]);
             else ++result.indeterminate_prescribed_adhesion_facets;}
            result.converged=true;
        }
        else
        {
            std::vector<Vecd>y=result.dual;Real acceleration=1,current_objective=objective(result.dual);
            for(size_t iteration=1;iteration<=options_.max_iterations;++iteration)
            {
                auto makeStep=[&](const std::vector<Vecd>&point)
                {auto g=gradient(point);std::vector<Vecd>step(count,Vecd::Zero());for(size_t k=0;k<count;++k)if(active[k]){
                    if(metric[k]>0)step[k]=project(point[k]-g[k]/metric[k],k,sqrt_area[k]);
                    else{Vecd direction=effectiveOpening(facets_[k].projector*s[k],k);Real norm=direction.stableNorm();step[k]=norm>0?project(strength*sqrt_area[k]*(direction/norm),k,sqrt_area[k]):project(point[k],k,sqrt_area[k]);}}
                    return step;};
                std::vector<Vecd>candidate=makeStep(y);Real next_objective=objective(candidate);
                // Monotone restart prevents acceleration from retaining an
                // uphill step; a plain PG step uses the proven global bound.
                if(next_objective>current_objective+64*std::numeric_limits<Real>::epsilon()*std::max(std::abs(current_objective),std::abs(next_objective)))
                {y=result.dual;acceleration=1;candidate=makeStep(y);next_objective=objective(candidate);}
                auto g=gradient(candidate);Real residual=0;
                for(size_t k=0;k<count;++k)if(active[k]&&metric[k]>0)residual=std::max(residual,(metric[k]/sqrt_area[k])*(candidate[k]-project(candidate[k]-g[k]/metric[k],k,sqrt_area[k])).stableNorm());
                result.iterations=iteration;result.projected_opening_residual=residual;result.relative_residual=residual/opening_scale;
                if(residual<=tolerance){result.dual=std::move(candidate);result.converged=true;break;}
                Real next_acceleration=options_.accelerate?.5*(1+std::sqrt(1+4*acceleration*acceleration)):1;
                for(size_t k=0;k<count;++k)y[k]=candidate[k]+((acceleration-1)/next_acceleration)*(candidate[k]-result.dual[k]);
                result.dual=std::move(candidate);acceleration=next_acceleration;current_objective=next_objective;
            }
            if(!result.converged)
            {std::ostringstream error;error<<"Coupled cohesive IMEX failed to converge after "<<result.iterations<<" iterations; projected opening residual="<<result.projected_opening_residual<<" m, tolerance="<<tolerance<<" m";throw std::runtime_error(error.str());}
        }
        result.dual_objective=objective(result.dual);
        result.monotone_traction.assign(count,Vecd::Zero());result.traction.assign(count,Vecd::Zero());result.law_traction.assign(count,Vecd::Zero());
        std::vector<Vecd>surface_force(count,Vecd::Zero());
        for(size_t k=0;k<count;++k)if(active[k])
        {result.monotone_traction[k]=result.dual[k]/sqrt_area[k];result.traction[k]=result.monotone_traction[k]+result.softening_traction[k];surface_force[k]=-facets_[k].reference_area*result.traction[k];Real range_error=(result.monotone_traction[k]-facets_[k].projector*result.monotone_traction[k]).stableNorm();
         result.maximum_dual_range_error=std::max(result.maximum_dual_range_error,range_error);
         result.maximum_dual_constraint_error=std::max({result.maximum_dual_constraint_error,range_error,std::max(Real(0),result.monotone_traction[k].stableNorm()-strength),std::max(Real(0),-result.monotone_traction[k].dot(facets_[k].traction_halfspace_normal))});}
        if(result.maximum_dual_constraint_error>strength*(options_.relative_tolerance+256*std::numeric_limits<Real>::epsilon()))
            throw std::runtime_error("Cohesive IMEX dual violates its ball/range feasibility tolerance");
        result.nodal_force=applyBT(surface_force);result.correction.assign(qfree.size(),Vecd::Zero());
        for(size_t i=0;i<qfree.size();++i)
        {result.correction[i]=dt*dt*inverse_mass_[i].cwiseProduct(result.nodal_force[i]);result.positions[i]+=result.correction[i];if(!result.positions[i].allFinite())throw std::runtime_error("Nonfinite cohesive corrected position");}
        result.jump=applyB(result.positions,true);for(size_t k=0;k<count;++k)result.jump[k]=effectiveOpening(result.jump[k],k);
        for(size_t k=0;k<count;++k)
        {
            Real old_k=old_history[k].kappa,new_k=std::max(old_k,result.jump[k].stableNorm()),A=facets_[k].reference_area;
            result.history[k].kappa=new_k;
            if(active[k])
            {
                Real r=result.jump[k].stableNorm();
                Real phi=old_k>0&&r<=old_k?.5*strength*r*r/old_k:strength*(r-.5*old_k);
                result.fenchel_gap+=A*(phi+.5*old_k/strength*result.monotone_traction[k].squaredNorm()-result.monotone_traction[k].dot(result.jump[k]));
            }
            Real old_D=.5*strength*std::min(old_k,delta_f),new_D=.5*strength*std::min(new_k,delta_f);
            Real old_psi=storedEnergy(projected_old_jump[k],old_k),new_psi=storedEnergy(result.jump[k],new_k);
            result.previous_stored_energy+=A*old_psi;result.stored_energy+=A*new_psi;
            result.dissipation+=A*new_D;result.dissipation_increment+=A*(new_D-old_D);
            result.interface_work+=A*result.traction[k].dot(result.jump[k]-projected_old_jump[k]);
            if(new_k>0&&new_k<delta_f)result.law_traction[k]=strength*(1/new_k-1/delta_f)*result.jump[k];
            else if(new_k==0)result.law_traction[k]=result.traction[k]; // set-valued intact state
            result.maximum_split_traction_error=std::max(result.maximum_split_traction_error,(result.traction[k]-result.law_traction[k]).stableNorm());
        }
        result.algorithmic_work_residual=result.interface_work-(result.stored_energy-result.previous_stored_energy)-result.dissipation_increment;
        for(size_t i=0;i<qfree.size();++i)result.force_position_residual=std::max(result.force_position_residual,(result.correction[i]-dt*dt*inverse_mass_[i].cwiseProduct(result.nodal_force[i])).norm());
        if(!std::isfinite(result.stored_energy)||!std::isfinite(result.dissipation)||!std::isfinite(result.algorithmic_work_residual)||!std::isfinite(result.dual_objective))
            throw std::runtime_error("Nonfinite coupled cohesive energy budget");
        return result;
    }
  private:
    std::vector<Vecd>inverse_mass_;std::vector<Facet>facets_;Options options_;
    void requireNodes(const std::vector<Vecd>&q)const
    {if(q.size()!=inverse_mass_.size())throw std::invalid_argument("Wrong cohesive nodal vector count");for(const Vecd&v:q)if(!v.allFinite())throw std::invalid_argument("Nonfinite cohesive nodal state");}
    Vecd project(const Vecd&z,size_t k,Real sqrt_area)const
    {Vecd p=facets_[k].projector*z,n=facets_[k].projector*facets_[k].traction_halfspace_normal;Real n2=n.squaredNorm();
        // Exact proximal projection onto range(P) intersected with the tensile
        // half-ball. The default zero normal retains the original full ball.
        if(n2>1e-24&&p.dot(n)<0)p-=n*(p.dot(n)/n2);
        Real radius=strength*sqrt_area,norm=p.stableNorm();if(norm>radius)p*=radius/norm;return p;}
    Vecd effectiveOpening(Vecd jump,size_t k)const {
        Vecd n=facets_[k].projector*facets_[k].traction_halfspace_normal;Real n2=n.squaredNorm();
        if(n2>1e-24&&jump.dot(n)<0)jump-=n*(jump.dot(n)/n2);return jump;
    }
    Real storedEnergy(const Vecd&jump,Real kappa)const
    {if(kappa==0||kappa>=delta_f)return 0;return .5*strength*(1/kappa-1/delta_f)*jump.squaredNorm();}
    std::vector<Real> lipschitzBounds(const std::vector<bool>&active,const std::vector<Real>&sqrt_area,
        const std::vector<Real>&compliance,Real dt)const
    {
        std::vector<Vecd>column_abs_sum(inverse_mass_.size(),Vecd::Zero());
        for(size_t k=0;k<facets_.size();++k)if(active[k])for(const Entry&e:facets_[k].entries)
            for(size_t d=0;d<3;++d)column_abs_sum[e.node][d]+=sqrt_area[k]*(facets_[k].projector*e.block()).col(d).cwiseAbs().sum();
        // Each scalar block dominates every absolute Hessian row in that
        // facet. For the symmetric positive-semidefinite dual Hessian H,
        // diag(metric[k] I)-H is diagonally dominant and positive semidefinite.
        // This is ordinary FISTA after the fixed change zhat=sqrt(metric)*z;
        // scalar blocks preserve the exact ball/half-ball proximal projection.
        // A soft old facet's kappa/ft cannot throttle an unrelated fresh facet.
        std::vector<Real>metric(facets_.size(),0);
        for(size_t k=0;k<facets_.size();++k)if(active[k])for(size_t a=0;a<3;++a)
        {Real row=compliance[k];for(const Entry&e:facets_[k].entries)for(size_t d=0;d<3;++d)row+=dt*dt*sqrt_area[k]*std::abs((facets_[k].projector*e.block())(a,d))*inverse_mass_[e.node][d]*column_abs_sum[e.node][d];metric[k]=std::max(metric[k],row);}
        for(Real value:metric)if(!std::isfinite(value))throw std::runtime_error("Nonfinite cohesive IMEX block majorizer");return metric;
    }
};
} // namespace SPH::peridynamics
