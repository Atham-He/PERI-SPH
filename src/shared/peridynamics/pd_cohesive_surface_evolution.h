#pragma once
#include "pd_cohesive_interface_imex.h"
#include "pd_solid_reference_coupling.h"
#include <array>
#include <map>

namespace SPH::peridynamics::cohesive_geometry_detail {
template<size_t N>struct Jet {
    Real v=0;std::array<Real,N>d{};Jet(Real a=0):v(a){}
    static Jet variable(Real a,size_t i){Jet r(a);r.d[i]=1;return r;}
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
    struct Geometry {Vecd jump=Vecd::Zero(),normal=Vecd::Zero();Real normal_gap=0;CohesiveInterfaceIMEX::Facet linear;};
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
    Geometry geometry(size_t k,const Field&q)const {
        using J=cohesive_geometry_detail::Jet<15>;
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
        for(size_t d=0;d<3;++d)effective[d]=gn.v>=0?jump[d]:jump[d]-gn*n[d];
        Geometry g;g.normal_gap=gn.v;g.linear.reference_area=facet.area;std::map<size_t,Matd>nodes;
        for(size_t d=0;d<3;++d){g.jump[d]=effective[d].v;g.normal[d]=n[d].v;}
        for(size_t j=0;j<5;++j){Matd derivative;for(size_t d=0;d<3;++d)for(size_t e=0;e<3;++e)derivative(d,e)=effective[d].d[3*j+e];const auto&t=facet.trace[j];const auto&row=t.surface->rows[t.point];
            for(size_t b=0;b<row.index.size();++b)nodes.try_emplace(row.index[b],Matd::Zero()).first->second+=row.weight[b]*derivative;
        }
        g.linear.jump_offset=g.jump;
        for(const auto&entry:nodes){g.linear.entries.push_back({entry.first,0,entry.second});g.linear.jump_offset-=entry.second*q[entry.first];}
        return g;
    }
    Real stored(const Vecd&jump,Real kappa)const {return kappa>0&&kappa<delta_f?.5*ft*(1/kappa-1/delta_f)*jump.squaredNorm():0;}
    Real dissipation(Real kappa)const{return .5*ft*std::min(kappa,delta_f);}
    Result solve(const Field&qfree,const Field&qold,const std::vector<History>&old,const Field&inverse_mass,Real dt)const {
        if(qfree.size()!=qold.size()||qold.size()!=inverse_mass.size()||old.size()!=quadrature.size())throw std::invalid_argument("Cohesive evolution dimensions");
        Result result;result.positions=qfree;Field iterate=qfree,previous_jump(old.size()),warm;
        const Real tolerance=absolute_tolerance+relative_tolerance*delta_f;
        for(size_t k=0;k<old.size();++k){previous_jump[k]=geometry(k,qold).jump;
            if(old[k].kappa==0&&previous_jump[k].norm()<=tolerance)previous_jump[k].setZero();
            if(previous_jump[k].norm()>old[k].kappa+tolerance&&old[k].kappa<delta_f)throw std::runtime_error("Cohesive committed history below previous opening");}
        CohesiveInterfaceIMEX::Result inner;
        for(size_t outer=1;outer<=max_outer_iterations;++outer){std::vector<CohesiveInterfaceIMEX::Facet>linear;Field theta;
            for(size_t k=0;k<old.size();++k){auto g=geometry(k,iterate);if(old[k].kappa==0&&g.normal_gap<0&&g.jump.norm()<=tolerance){g.linear.projector=Matd::Identity()-g.normal*g.normal.transpose();g.linear.jump_offset=g.linear.projector*g.linear.jump_offset;g.jump=g.linear.projector*g.jump;}linear.push_back(g.linear);theta.push_back(g.jump);}
            CohesiveInterfaceIMEX::Options options;options.absolute_opening_tolerance=absolute_tolerance*.05;options.relative_tolerance=relative_tolerance*.05;
            CohesiveInterfaceIMEX solver(ft,Gc,inverse_mass,linear,options);
            inner=solver.solve(qfree,previous_jump,old,dt,warm,theta);result.dual_iterations+=inner.iterations;
            Real change=0,geometry_error=0;for(size_t i=0;i<iterate.size();++i)change=std::max(change,(inner.positions[i]-iterate[i]).norm());
            for(size_t k=0;k<old.size();++k)geometry_error=std::max(geometry_error,(geometry(k,inner.positions).jump-inner.jump[k]).norm());
            iterate=inner.positions;warm=inner.dual;result.outer_iterations=outer;
            if(change<=tolerance&&geometry_error<=tolerance){result.converged=true;break;}
        }
        if(!result.converged)throw std::runtime_error("Nonlinear objective cohesive drift did not converge; history unchanged");
        result.positions=iterate;result.force.assign(iterate.size(),Vecd::Zero());result.history=old;
        Vecd resultant=Vecd::Zero(),torque=Vecd::Zero();Real force_scale=0,length=0;
        for(size_t k=0;k<old.size();++k){auto g=geometry(k,iterate);Real r=g.jump.norm(),kap=std::max(old[k].kappa,r);Vecd traction;
            // Interior dual + zero opening is the set-valued intact state, not
            // irreversible damage caused by linear-solver round-off.
            bool adhered=old[k].kappa==0&&r<=tolerance&&inner.monotone_traction[k].norm()<ft*(1-1e-7);
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
        if(result.position_residual>8*tolerance||result.traction_residual>1e-5||result.dissipation_increment< -1e-20)throw std::runtime_error("Cohesive physical force/history residual failed; history unchanged");
        return result;
    }
};
}
