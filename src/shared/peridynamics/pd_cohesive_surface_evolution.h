#pragma once
#include "pd_cohesive_interface_imex.h"
#include "pd_solid_reference_coupling.h"
#include <array>
#include <map>
#include <fstream>
#include <iomanip>
#include <cstdlib>
#include <sstream>
#include <type_traits>

namespace SPH::peridynamics::cohesive_geometry_detail {
template<size_t N>struct Jet {
    Real v=0;std::array<Real,N>d{};Jet(Real a=0):v(a){}
    static Jet variable(Real a,size_t i){Jet r(a);if constexpr(N>0)r.d[i]=1;return r;}
};
template<size_t N>Jet<N>operator+(const Jet<N>&a,const Jet<N>&b){Jet<N>r(a.v+b.v);for(size_t i=0;i<N;++i)r.d[i]=a.d[i]+b.d[i];return r;}
template<size_t N>Jet<N>operator-(const Jet<N>&a,const Jet<N>&b){Jet<N>r(a.v-b.v);for(size_t i=0;i<N;++i)r.d[i]=a.d[i]-b.d[i];return r;}
template<size_t N>Jet<N>operator*(const Jet<N>&a,const Jet<N>&b){Jet<N>r(a.v*b.v);for(size_t i=0;i<N;++i)r.d[i]=a.d[i]*b.v+a.v*b.d[i];return r;}
template<size_t N>Jet<N>operator/(const Jet<N>&a,const Jet<N>&b){Jet<N>r(a.v/b.v);for(size_t i=0;i<N;++i)r.d[i]=(a.d[i]-r.v*b.d[i])/b.v;return r;}
template<size_t N>Jet<N>sqrt(const Jet<N>&a){Jet<N>r(std::sqrt(a.v));for(size_t i=0;i<N;++i)r.d[i]=a.d[i]/(2*r.v);return r;}
}
namespace SPH::peridynamics {
/** Objective, unilateral finite-geometry cohesive evolution.
 * The bulk is NOSB-PD; this class only integrates extrinsic surface work.
 * delta_f=2 Gc/ft. Compression is excluded from cohesion and belongs to contact.
 * The exact moving-normal derivative is pulled through the SAME one-sided
 * traces as its potential. No frozen-normal torque or artificial birth gap.
 *
 * Each step solves q=qfree+dt^2 M^-1 f_coh(q,kappa_old) by a geometric outer
 * iteration and the coupled monotone dual solver. The softening term is also
 * converged, rather than left lagged when a surface reaches complete failure.
 * History advances only after the nonlinear residual passes. This is a fully
 * implicit cohesive drift, suitable between two explicit bulk half-kicks.
 * It is not an unmodified implementation of Doyen Scheme B, nor a claim of
 * exact finite-dt energy conservation. The driver must audit the FULL budget.
 */
class CohesiveSurfaceEvolution {
 public:
    using Field=std::vector<Vecd>;using History=CohesiveInterfaceIMEX::History;
    struct TracePoint {ReferenceSolidSurface*surface=nullptr;size_t point=0;};
    struct Quadrature {
        // plus, matching minus, then three minus-face vertices for its normal.
        std::array<TracePoint,5> trace;Real area=0;
    };
    struct Kinematics {Vecd jump=Vecd::Zero(),normal=Vecd::Zero();Real normal_gap=0;bool normal_derivative_active=true;};
    struct Geometry : Kinematics {CohesiveInterfaceIMEX::Facet linear;};
    struct Result {
        Field positions,force,jump,traction;std::vector<History>history;
        Real stored=0,dissipated=0,dissipation_increment=0;
        Real position_residual=0,traction_residual=0,torque_residual=0,force_residual=0;
        size_t outer_iterations=0,dual_iterations=0;bool converged=false;
    };
    Real ft,Gc,delta_f;std::vector<Quadrature>quadrature;
    Real absolute_tolerance=2e-14,relative_tolerance=2e-9;
    size_t max_outer_iterations=80;
    CohesiveSurfaceEvolution(Real strength,Real energy,const std::vector<Quadrature>&q):ft(strength),Gc(energy),delta_f(2*energy/strength),quadrature(q){
        if(!std::isfinite(ft)||!std::isfinite(Gc)||std::min(ft,Gc)<=0)throw std::invalid_argument("Invalid extrinsic cohesive inputs");
        for(const auto&a:q){if(!std::isfinite(a.area)||a.area<=0)throw std::invalid_argument("Invalid cohesive reference area");for(const auto&t:a.trace)if(!t.surface||t.point>=t.surface->rows.size())throw std::invalid_argument("Invalid cohesive geometry trace");}
    }
  private:
    // Position and gap checks share the force geometry, but do not need its
    // automatic derivatives or a per-node linearization.
    template<bool AssembleLinear>
    auto evaluateGeometry(size_t k,const Field&q,Real zero_opening_tolerance)const {
        using J=cohesive_geometry_detail::Jet<AssembleLinear?15:0>;
        using Evaluation=std::conditional_t<AssembleLinear,Geometry,Kinematics>;
        const auto&facet=quadrature.at(k);std::array<std::array<J,3>,5>a;
        for(size_t j=0;j<5;++j){const auto&t=facet.trace[j];t.surface->ensureCurrent();Vecd x=Vecd::Zero();const auto&row=t.surface->rows[t.point];
            for(size_t n=0;n<row.index.size();++n){if(row.index[n]>=q.size())throw std::invalid_argument("Cohesive geometry node dimensions");x+=row.weight[n]*q[row.index[n]];}
            for(size_t d=0;d<3;++d)a[j][d]=J::variable(x[d],3*j+d);
        }
        std::array<J,3>u,v,n,jump,effective;
        for(size_t d=0;d<3;++d){u[d]=a[3][d]-a[2][d];v[d]=a[4][d]-a[2][d];jump[d]=a[0][d]-a[1][d];}
        n={u[1]*v[2]-u[2]*v[1],u[2]*v[0]-u[0]*v[2],u[0]*v[1]-u[1]*v[0]};J norm2(0);
        for(const auto&z:n)norm2=norm2+z*z;
        if(norm2.v<=0||!std::isfinite(norm2.v))throw std::runtime_error("Cohesive face normal lost geometric rank");
        J norm=cohesive_geometry_detail::sqrt(norm2);for(auto&z:n)z=z/norm;J gn(0);for(size_t d=0;d<3;++d)gn=gn+jump[d]*n[d];
        Real raw_norm=0;for(size_t d=0;d<3;++d)raw_norm+=jump[d].v*jump[d].v;
        bool normal_active=gn.v>=0||(gn.v>=-zero_opening_tolerance&&std::sqrt(raw_norm)<=zero_opening_tolerance);
        for(size_t d=0;d<3;++d)effective[d]=normal_active?jump[d]:jump[d]-gn*n[d];
        Evaluation g;g.normal_gap=gn.v;g.normal_derivative_active=normal_active;
        for(size_t d=0;d<3;++d){g.jump[d]=effective[d].v;g.normal[d]=n[d].v;}
        if constexpr(AssembleLinear){
            g.linear.reference_area=facet.area;
            std::map<size_t,Matd>nodes;
            for(size_t j=0;j<5;++j){
                Matd derivative;
                for(size_t d=0;d<3;++d)for(size_t e=0;e<3;++e)derivative(d,e)=effective[d].d[3*j+e];
                const auto&t=facet.trace[j];const auto&row=t.surface->rows[t.point];
                for(size_t b=0;b<row.index.size();++b)nodes.try_emplace(row.index[b],Matd::Zero()).first->second+=row.weight[b]*derivative;
            }
            g.linear.jump_offset=g.jump;
            for(const auto&entry:nodes){g.linear.entries.push_back({entry.first,0,entry.second});g.linear.jump_offset-=entry.second*q[entry.first];}
        }
        return g;
    }
  public:
    Geometry geometry(size_t k,const Field&q,Real zero_opening_tolerance=0)const
    {return evaluateGeometry<true>(k,q,zero_opening_tolerance);}
    Kinematics kinematics(size_t k,const Field&q,Real zero_opening_tolerance=0)const
    {return evaluateGeometry<false>(k,q,zero_opening_tolerance);}
    Real stored(const Vecd&jump,Real kappa)const {return kappa>0&&kappa<delta_f?.5*ft*(1/kappa-1/delta_f)*jump.squaredNorm():0;}
    Real dissipation(Real kappa)const{return .5*ft*std::min(kappa,delta_f);}
    Result solve(const Field&qfree,const Field&qold,const std::vector<History>&old,const Field&inverse_mass,Real dt)const {
        if(qfree.size()!=qold.size()||qold.size()!=inverse_mass.size()||old.size()!=quadrature.size())throw std::invalid_argument("Cohesive evolution dimensions");
        Result result;result.positions=qfree;Field iterate=qfree,previous_jump(old.size()),warm;
        const Real tolerance=absolute_tolerance+relative_tolerance*delta_f;
        std::vector<Real> old_trace_roundoff(old.size(),0);
        for(size_t k=0;k<old.size();++k){
            previous_jump[k]=kinematics(k,qold).jump;
            const Real kappa=old[k].kappa;
            if(kappa==0&&previous_jump[k].norm()<=tolerance){previous_jump[k].setZero();continue;}
            const Real opening=previous_jump[k].stableNorm();
            if(kappa>=delta_f||opening<=kappa)continue;
            old_trace_roundoff[k]=openingTraceRoundoff(k,qold);
            if(opening-kappa>CohesiveInterfaceIMEX::historyRoundoff(delta_f,kappa)+old_trace_roundoff[k])
                throw std::runtime_error("Cohesive committed history below previous opening beyond trace roundoff");
        }
        CohesiveInterfaceIMEX::Result inner;std::ostringstream audit; audit<<std::setprecision(17);Real last_change=0,last_geometry_error=0,inner_precision=1;
        auto dumpFailure=[&](){if(const char*path=std::getenv("PERI_COHESIVE_AUDIT_FILE")){std::ofstream f(path);f<<std::setprecision(17)<<"COHESIVE_REPLAY_V1 "<<qfree.size()<<' '<<old.size()<<' '<<ft<<' '<<Gc<<' '<<dt<<'\n';for(size_t i=0;i<qfree.size();++i)f<<qfree[i].transpose()<<' '<<qold[i].transpose()<<' '<<inverse_mass[i].transpose()<<'\n';for(size_t k=0;k<old.size();++k){f<<quadrature[k].area<<' '<<old[k].kappa<<'\n';for(const auto&v:quadrature[k].trace){const auto&r=v.surface->rows[v.point];f<<r.index.size();for(size_t j=0;j<r.index.size();++j)f<<' '<<r.index[j]<<' '<<r.weight[j];f<<'\n';}}f<<audit.str();}};
        for(size_t outer=1;outer<=max_outer_iterations;++outer){std::vector<CohesiveInterfaceIMEX::Facet>linear;Field theta;
            for(size_t k=0;k<old.size();++k){Real zero=old[k].kappa==0?tolerance:0;auto g=geometry(k,iterate,zero);
                // At zero opening the correct subgradient is the tensile
                // half-ball. Treat only the existing convergence-sized zero
                // neighbourhood as that boundary; do not alternate unilateral
                // branches on round-off-sized signed gaps.
                g.linear.traction_halfspace_normal=g.normal;
                if(!g.normal_derivative_active){g.linear.projector=Matd::Identity()-g.normal*g.normal.transpose();g.linear.jump_offset=g.linear.projector*g.linear.jump_offset;g.jump=g.linear.projector*g.jump;}
                linear.push_back(g.linear);Vecd soft=g.jump;if(g.normal_gap<0)soft-=g.normal*g.normal.dot(soft);theta.push_back(soft);}
            CohesiveInterfaceIMEX::Options options;options.absolute_opening_tolerance=absolute_tolerance*.05*inner_precision;options.relative_tolerance=relative_tolerance*.05*inner_precision;
            CohesiveInterfaceIMEX solver(ft,Gc,inverse_mass,linear,options);
            Field projected_previous=previous_jump;for(size_t k=0;k<old.size();++k)projected_previous[k]=linear[k].projector*projected_previous[k];
            try{inner=solver.solve(qfree,projected_previous,old,dt,warm,theta,old_trace_roundoff);}catch(...){dumpFailure();throw;}result.dual_iterations+=inner.iterations;
            Real change=0,geometry_error=0;for(size_t i=0;i<iterate.size();++i)change=std::max(change,(inner.positions[i]-iterate[i]).norm());
            for(size_t k=0;k<old.size();++k)geometry_error=std::max(geometry_error,(kinematics(k,inner.positions).jump-inner.jump[k]).norm());
            last_change=change;last_geometry_error=geometry_error;
            if(std::getenv("PERI_COHESIVE_AUDIT_FILE")){audit<<"OUTER "<<outer<<" change "<<change<<" geometry "<<geometry_error<<" innerKKT "<<inner.projected_opening_residual<<'\n';for(size_t k=0;k<old.size();++k){auto g=kinematics(k,inner.positions);audit<<k<<" gn "<<g.normal_gap<<" jump "<<g.jump.transpose()<<" traction "<<inner.traction[k].transpose()<<" kappa "<<old[k].kappa<<'\n';}}
            iterate=inner.positions;warm=inner.dual;result.outer_iterations=outer;
            Real law_error=0;for(size_t k=0;k<old.size();++k){auto g=kinematics(k,iterate,old[k].kappa==0?tolerance:0);Real r=g.jump.norm(),kap=std::max(old[k].kappa,r);Vecd physical=inner.traction[k];
                if(kap>=delta_f)physical.setZero();else if(!(old[k].kappa==0&&r<=tolerance)&&kap>0)physical=ft*(1/kap-1/delta_f)*g.jump;
                law_error=std::max(law_error,(physical-inner.traction[k]).norm()/ft);}
            if(std::getenv("PERI_COHESIVE_AUDIT_FILE"))audit<<"LAW_ERROR "<<law_error<<" inner_precision "<<inner_precision<<'\n';
            if(change<=tolerance&&geometry_error<=tolerance){
                if(law_error<=1e-5){result.converged=true;break;}
                // Force consistency also belongs to nonlinear convergence.
                // Refine the SAME convex subproblem, not material parameters or
                // the acceptance threshold, when an emerging opening makes its
                // traction more sensitive than the generic displacement test.
                Real floor=16*std::numeric_limits<Real>::epsilon()*delta_f/(.05*(absolute_tolerance+relative_tolerance*delta_f));
                inner_precision=std::max(floor,inner_precision*std::clamp(.5e-5/law_error,Real(.1),Real(.5)));
            }
        }
        if(!result.converged){
            dumpFailure();
            std::ostringstream message;message<<std::setprecision(17)<<"Nonlinear objective cohesive drift did not converge; history unchanged; change="<<last_change<<" geometry="<<last_geometry_error<<" KKT="<<inner.projected_opening_residual;throw std::runtime_error(message.str());}
        result.positions=iterate;result.force.assign(iterate.size(),Vecd::Zero());result.history=old;
        Vecd resultant=Vecd::Zero(),torque=Vecd::Zero();Real force_scale=0,length=0;
        for(size_t k=0;k<old.size();++k){auto g=geometry(k,iterate,old[k].kappa==0?tolerance:0);Real r=g.jump.norm(),kap=std::max(old[k].kappa,r);Vecd traction;
            // The whole feasible half-ball, including its strength boundary,
            // is a valid subgradient at zero opening. Round-off at that boundary
            // must not create a new maximum opening and irreversible damage.
            bool adhered=old[k].kappa==0&&r<=tolerance;
            if(adhered){kap=0;traction=inner.traction[k];}
            else if(kap>=delta_f)traction.setZero();
            else if(kap>0)traction=ft*(1/kap-1/delta_f)*g.jump;
            else traction=inner.traction[k];
            result.traction_residual=std::max(result.traction_residual,(traction-inner.traction[k]).norm()/ft);
            result.history[k].kappa=kap;result.jump.push_back(g.jump);result.traction.push_back(traction);
            Real A=quadrature[k].area;result.stored+=A*stored(g.jump,kap);result.dissipated+=A*dissipation(kap);result.dissipation_increment+=A*(dissipation(kap)-dissipation(old[k].kappa));
            for(const auto&e:g.linear.entries)result.force[e.node]-=A*e.block().transpose()*traction;
        }
        for(size_t i=0;i<iterate.size();++i){result.position_residual=std::max(result.position_residual,(iterate[i]-qfree[i]-dt*dt*inverse_mass[i].cwiseProduct(result.force[i])).norm());resultant+=result.force[i];torque+=iterate[i].cross(result.force[i]);force_scale+=result.force[i].norm();length=std::max(length,iterate[i].norm());}
        result.force_residual=resultant.norm()/std::max(force_scale,Real(1e-30));result.torque_residual=torque.norm()/std::max(force_scale*length,Real(1e-30));
        if(result.position_residual>8*tolerance||result.traction_residual>1e-5||result.dissipation_increment< -1e-20){dumpFailure();std::ostringstream message;message<<std::setprecision(17)<<"Cohesive physical force/history residual failed; history unchanged; position="<<result.position_residual<<" traction="<<result.traction_residual<<" dissipation="<<result.dissipation_increment;throw std::runtime_error(message.str());}
        return result;
    }

  private:
    Real openingTraceRoundoff(size_t k,const Field&q)const
    {
        // Each position component is a signed dot product. Bound its rounding
        // by gamma_(n+2) times the sum of absolute products. The two face
        // positions can be much larger than their difference. The extra terms
        // cover final subtraction and comparison with the saved opening. This
        // bound is for validation only. No opening, energy or history is reset.
        const Real unit_roundoff=std::numeric_limits<Real>::epsilon()/2;
        Real bound=0;
        for(size_t side=0;side<2;++side){
            const auto&trace=quadrature[k].trace[side];
            const auto&row=trace.surface->rows[trace.point];
            Vecd magnitude=Vecd::Zero();
            for(size_t a=0;a<row.index.size();++a)
                magnitude+=std::abs(row.weight[a])*q[row.index[a]].cwiseAbs();
            const Real gamma=(Real(row.index.size())+2)*unit_roundoff;
            bound+=2*(gamma/(1-gamma))*magnitude.stableNorm();
        }
        return bound;
    }
};
}
