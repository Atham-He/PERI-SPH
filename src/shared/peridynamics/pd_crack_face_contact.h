#pragma once
/* Licensed under the Apache License, Version 2.0. */
#include "pd_solid_reference_coupling.h"
#include <array>

namespace SPH::peridynamics::crack_contact_detail
{
/** Exact local derivatives in the slave-point and three master-vertex
 * coordinates. Used only for the conservative contact tangent bound. */
struct Jet12
{
    Real value=0;
    std::array<Real,12> gradient{};
    std::array<Real,144> hessian{};
    Jet12(Real x=0):value(x){}
    static Jet12 variable(Real x,size_t i){Jet12 r(x);r.gradient[i]=1;return r;}
};
inline Jet12 operator+(const Jet12&a,const Jet12&b)
{Jet12 r(a.value+b.value);for(size_t i=0;i<12;++i)r.gradient[i]=a.gradient[i]+b.gradient[i];for(size_t i=0;i<144;++i)r.hessian[i]=a.hessian[i]+b.hessian[i];return r;}
inline Jet12 operator-(const Jet12&a,const Jet12&b)
{Jet12 r(a.value-b.value);for(size_t i=0;i<12;++i)r.gradient[i]=a.gradient[i]-b.gradient[i];for(size_t i=0;i<144;++i)r.hessian[i]=a.hessian[i]-b.hessian[i];return r;}
inline Jet12 operator*(const Jet12&a,const Jet12&b)
{Jet12 r(a.value*b.value);for(size_t i=0;i<12;++i)r.gradient[i]=a.gradient[i]*b.value+a.value*b.gradient[i];for(size_t i=0;i<12;++i)for(size_t j=0;j<12;++j)r.hessian[12*i+j]=a.hessian[12*i+j]*b.value+a.gradient[i]*b.gradient[j]+a.gradient[j]*b.gradient[i]+a.value*b.hessian[12*i+j];return r;}
inline Jet12 inverse(const Jet12&a)
{Jet12 r(1/a.value);Real first=-1/(a.value*a.value),second=2/(a.value*a.value*a.value);for(size_t i=0;i<12;++i)r.gradient[i]=first*a.gradient[i];for(size_t i=0;i<12;++i)for(size_t j=0;j<12;++j)r.hessian[12*i+j]=first*a.hessian[12*i+j]+second*a.gradient[i]*a.gradient[j];return r;}
inline Jet12 operator/(const Jet12&a,const Jet12&b){return a*inverse(b);}
using JVec=std::array<Jet12,3>;
inline JVec subtract(const JVec&a,const JVec&b){return {a[0]-b[0],a[1]-b[1],a[2]-b[2]};}
inline JVec cross(const JVec&a,const JVec&b){return {a[1]*b[2]-a[2]*b[1],a[2]*b[0]-a[0]*b[2],a[0]*b[1]-a[1]*b[0]};}
inline Jet12 dot(const JVec&a,const JVec&b){return a[0]*b[0]+a[1]*b[1]+a[2]*b[2];}
inline Jet12 energy(const std::array<Vecd,4>&x,Real stiffness)
{
    std::array<JVec,4>a;
    for(size_t i=0;i<4;++i)for(size_t d=0;d<3;++d)a[i][d]=Jet12::variable(x[i][d],3*i+d);
    JVec normal=cross(subtract(a[2],a[1]),subtract(a[3],a[1]));
    Jet12 numerator=dot(subtract(a[0],a[1]),normal);
    return Jet12(.5*stiffness)*numerator*numerator/dot(normal,normal);
}
}

namespace SPH::peridynamics
{
/** Frictionless signed-gap contact for KNOWN opposite crack-face patches.
 *
 * The caller supplies a slave surface quadrature, master triangular geometry,
 * and one-sided material traces. Both surfaces may belong to the SAME body.
 * No global fragment-id test, material deletion, fracture work, or tangential
 * spring is introduced. Master triangles must point towards the slave side.
 *
 * For an interior normal projection p=sum beta_k x_k and gap g=n.(x_s-a),
 * differentiating the plane identity gives exactly
 *   dg=n.(dx_s-sum beta_k dx_k).
 * The normal-rotation terms have therefore NOT been omitted from the first
 * derivative. U=1/2 K A_ref min(g,0)^2 gives the forces below. The local second
 * derivative (including all normal-rotation terms) is differentiated exactly,
 * pulled back through both material traces and bounded in physical masses.
 *
 * This is NOT a general closed-body collision/edge-contact search. Projections
 * outside the supplied facets are counted and receive no fictitious radial
 * force. Supply sufficient master margin for the validated sliding range.
 * Candidate lists can be rebuilt by the caller for efficient moving geometry;
 * an empty outer list means scan all supplied master triangles.
 */
class FrictionlessCrackFaceContact
{
  public:
    using Triangle=std::array<size_t,3>;
    struct ContactRecord
    {size_t slave_point,master_triangle;Real gap;Vecd normal,barycentric,slave_force;};
    ReferenceSolidSurface &slave,&master;
    StdVec<Triangle> triangles;
    StdVec<Real> reference_areas;
    Real penalty;
    Real energy=0,max_penetration=0,normal_force_sum=0;
    Real force_residual=0,torque_residual=0,power_residual=0;
    Real force_gradient_relative_error=0,tangential_force_sum=0;
    Real energy_rate=0,mechanical_power=0,stable_dt=std::numeric_limits<Real>::infinity();
    Real mass_normalized_tangent_bound=0;
    size_t active_contacts=0,outside_projections=0,boundary_projections=0;
    StdVec<ContactRecord> contacts;

    FrictionlessCrackFaceContact(ReferenceSolidSurface&s,ReferenceSolidSurface&m,
        const StdVec<Triangle>&faces,const StdVec<Real>&areas,Real stiffness_per_area,
        const StdVec<StdVec<size_t>>&candidate_triangles={})
      :slave(s),master(m),triangles(faces),reference_areas(areas),penalty(stiffness_per_area),candidates(candidate_triangles)
    {
        static_assert(Dimensions==3,"Crack-face contact requires 3D geometry");
        if(triangles.empty()||areas.size()!=slave.rows.size()||penalty<=0||!std::isfinite(penalty))
            throw std::invalid_argument("Crack contact requires faces, slave areas and positive finite penalty");
        for(Real a:areas)if(a<=0||!std::isfinite(a))throw std::invalid_argument("Invalid crack quadrature area");
        for(const auto&t:triangles)
            if(*std::max_element(t.begin(),t.end())>=master.rows.size()||t[0]==t[1]||t[0]==t[2]||t[1]==t[2])
                throw std::invalid_argument("Invalid master crack-face triangle");
        if(slave.points.size()!=slave.rows.size()||master.points.size()!=master.rows.size())throw std::invalid_argument("Crack trace point/row dimensions disagree");
        validateCandidates();
        rates_slave.resize(slave.port.count);rates_master.resize(master.port.count);
        master_x.resize(master.rows.size());master_v.resize(master.rows.size());
        alias=slave.port.position==master.port.position;
        if(alias&&(slave.port.count!=master.port.count||slave.port.mass!=master.port.mass))
            throw std::invalid_argument("Inconsistent aliases of one physical crack body");
        reference_length=0;
        for(const auto&t:triangles)for(size_t k=0;k<3;++k)
            reference_length=std::max(reference_length,(master.points[t[k]]-master.points[t[(k+1)%3]]).norm());
        if(reference_length<=0)throw std::invalid_argument("Degenerate reference crack faces");
    }
    void setCandidates(const StdVec<StdVec<size_t>>&list){candidates=list;validateCandidates();}
    void apply()
    {
        if(slave.rows.size()!=reference_areas.size()||master.rows.size()!=master_x.size())throw std::runtime_error("Crack trace dimensions changed after contact construction");
        energy=max_penetration=normal_force_sum=force_residual=torque_residual=power_residual=0;
        force_gradient_relative_error=tangential_force_sum=energy_rate=mechanical_power=0;
        active_contacts=outside_projections=boundary_projections=0;contacts.clear();
        for(auto&r:rates_slave)r.setZero();for(auto&r:rates_master)r.setZero();
        Vecd resultant=Vecd::Zero(),moment=Vecd::Zero();Real force_scale=0,power_scale=0;
        for(size_t q=0;q<master.rows.size();++q){validateRow(master,q);master.kinematics(q,master_x[q],master_v[q]);if(!master_x[q].allFinite()||!master_v[q].allFinite())throw std::runtime_error("Nonfinite master crack trace");}
        for(size_t q=0;q<slave.rows.size();++q)
        {
            validateRow(slave,q);Vecd x,v;slave.kinematics(q,x,v);if(!x.allFinite()||!v.allFinite())throw std::runtime_error("Nonfinite slave crack trace");
            size_t chosen=triangles.size();Real best=std::numeric_limits<Real>::infinity(),gap=0;
            Vecd bary=Vecd::Zero(),normal=Vecd::Zero();
            auto consider=[&](size_t index)
            {
                const auto&t=triangles[index];const Vecd&a=master_x[t[0]];Vecd u=master_x[t[1]]-a,w=master_x[t[2]]-a;
                Vecd cross=u.cross(w);Real u2=u.squaredNorm(),w2=w.squaredNorm(),uw=u.dot(w),den=cross.squaredNorm();
                if(!std::isfinite(den)||den<=64*std::numeric_limits<Real>::epsilon()*u2*w2)
                    throw std::runtime_error("Crack-face contact triangle lost geometric rank");
                Vecd n=cross/std::sqrt(den),d=x-a;Real b=(w2*d.dot(u)-uw*d.dot(w))/den,c=(u2*d.dot(w)-uw*d.dot(u))/den;
                Vecd beta(1-b-c,b,c);const Real tol=64*std::numeric_limits<Real>::epsilon();
                if(beta.minCoeff() < -tol)return;
                Real g=n.dot(d),distance=std::abs(g);
                if(distance<best){best=distance;chosen=index;gap=g;bary=beta;normal=n;}
            };
            if(candidates.empty())for(size_t f=0;f<triangles.size();++f)consider(f);
            else for(size_t f:candidates[q])consider(f);
            if(chosen==triangles.size()){++outside_projections;continue;}
            if(bary.minCoeff()<=64*std::numeric_limits<Real>::epsilon())++boundary_projections;
            const auto&t=triangles[chosen];Real k=penalty*reference_areas[q];
            std::array<Vecd,4> coordinates{x,master_x[t[0]],master_x[t[1]],master_x[t[2]]};
            // Keep a prospective closed-contact bound even while open, so the
            // first contact does not start from an unconstrained/infinite dt.
            if(gap>=0)coordinates[0]-=gap*normal;
            auto local=crack_contact_detail::energy(coordinates,k);
            std::array<const ReferenceSolidSurface::Row*,4> trace{&slave.rows[q],&master.rows[t[0]],&master.rows[t[1]],&master.rows[t[2]]};
            accumulateBound(local,trace);
            if(gap>=0)continue;
            Vecd f=-k*gap*normal;
            std::array<Vecd,4> forces{f,-bary[0]*f,-bary[1]*f,-bary[2]*f};
            std::array<Vecd,4> velocities{v,master_v[t[0]],master_v[t[1]],master_v[t[2]]};
            Real error2=0;
            for(size_t a=0;a<4;++a)for(size_t d=0;d<3;++d)error2+=std::pow(forces[a][d]+local.gradient[3*a+d],2);
            force_gradient_relative_error=std::max(force_gradient_relative_error,std::sqrt(error2)/
                (k*std::max(-gap,std::numeric_limits<Real>::epsilon()*reference_length)));
            slave.scatter(q,f);for(size_t a=0;a<3;++a)master.scatter(t[a],forces[a+1]);
            for(size_t a=0;a<4;++a)
            {resultant+=forces[a];moment+=coordinates[a].cross(forces[a]);Real work=forces[a].dot(velocities[a]);mechanical_power+=work;power_scale+=std::abs(work);force_scale+=forces[a].norm();}
            Vecd relative_velocity=v-bary[0]*master_v[t[0]]-bary[1]*master_v[t[1]]-bary[2]*master_v[t[2]];
            energy_rate+=k*gap*normal.dot(relative_velocity);energy+=local.value;
            max_penetration=std::max(max_penetration,-gap);normal_force_sum+=f.dot(normal);
            tangential_force_sum+=(f-normal*f.dot(normal)).norm();++active_contacts;
            contacts.push_back({q,chosen,gap,normal,bary,f});
        }
        mass_normalized_tangent_bound=0;
        for(size_t i=0;i<rates_slave.size();++i)
            for(size_t d=0;d<3;++d)mass_normalized_tangent_bound=std::max(mass_normalized_tangent_bound,rates_slave[i][d]+(alias?rates_master[i][d]:0));
        if(!alias)for(const auto&r:rates_master)for(size_t d=0;d<3;++d)mass_normalized_tangent_bound=std::max(mass_normalized_tangent_bound,r[d]);
        stable_dt=mass_normalized_tangent_bound>0?.2/std::sqrt(mass_normalized_tangent_bound):std::numeric_limits<Real>::infinity();
        force_residual=resultant.norm()/std::max(force_scale,Real(1e-30));
        torque_residual=moment.norm()/std::max(reference_length*force_scale,Real(1e-30));
        power_residual=std::abs(mechanical_power+energy_rate)/std::max({power_scale,std::abs(energy_rate),Real(1e-30)});
        if(!std::isfinite(energy)||!std::isfinite(max_penetration)||!std::isfinite(mass_normalized_tangent_bound))
            throw std::runtime_error("Nonfinite crack-face contact response");
    }
  private:
    StdVec<StdVec<size_t>> candidates;
    StdVec<Vecd> rates_slave,rates_master,master_x,master_v;
    bool alias=false;Real reference_length=0;
    void validateCandidates()const
    {
        if(!candidates.empty()&&candidates.size()!=slave.rows.size())throw std::invalid_argument("Crack contact candidate dimensions disagree");
        for(const auto&row:candidates)for(size_t f:row)if(f>=triangles.size())throw std::invalid_argument("Invalid crack contact candidate");
    }
    static void validateRow(const ReferenceSolidSurface&s,size_t q)
    {
        const auto&r=s.rows[q];if(r.index.empty()||r.index.size()!=r.weight.size())throw std::runtime_error("Unresolved one-sided crack surface trace");
        Real sum=0;for(size_t a=0;a<r.index.size();++a)
        {if(r.index[a]>=s.port.count||!std::isfinite(r.weight[a]))throw std::runtime_error("Invalid crack surface trace");sum+=r.weight[a];}
        if(std::abs(sum-1)>1e-8)throw std::runtime_error("Crack surface trace lost partition of unity");
    }
    void accumulateBound(const crack_contact_detail::Jet12&local,
                         const std::array<const ReferenceSolidSurface::Row*,4>&trace)
    {
        for(Real value:local.hessian)if(!std::isfinite(value))throw std::runtime_error("Nonfinite crack-contact Hessian");
        std::array<Real,4> norm{};
        for(size_t a=0;a<4;++a){const auto&port=a==0?slave.port:master.port;const auto&r=*trace[a];
            for(size_t j=0;j<r.index.size();++j)norm[a]+=std::abs(r.weight[j])/std::sqrt(port.mass[r.index[j]]);}
        for(size_t a=0;a<4;++a)for(size_t d=0;d<3;++d)
        {
            Real row=0;size_t i=3*a+d;
            for(size_t j=0;j<12;++j)row+=std::abs(local.hessian[12*i+j])*norm[j/3];
            const auto&port=a==0?slave.port:master.port;const auto&r=*trace[a];auto&rate=a==0?rates_slave:rates_master;
            for(size_t j=0;j<r.index.size();++j)rate[r.index[j]][d]+=std::abs(r.weight[j])*row/std::sqrt(port.mass[r.index[j]]);
        }
    }
};
}
