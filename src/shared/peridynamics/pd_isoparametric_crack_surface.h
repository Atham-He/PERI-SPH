#pragma once
#include "pd_fracture_visible_surface.h"
#include "pd_crack_projection.h"
#include <array>
#include <map>
namespace SPH::peridynamics {
struct BilinearCrackGeometry {
    using Param=Eigen::Vector2d;
    struct Value {std::array<Real,4>N,Nu,Nv,Nuv;Vecd x=Vecd::Zero(),u=Vecd::Zero(),v=Vecd::Zero(),uv=Vecd::Zero();};
    static Value evaluate(const std::array<Vecd,4>&x,const Param&p){Real s=p[0],t=p[1];Value a;
        a.N={(1-s)*(1-t),s*(1-t),s*t,(1-s)*t};a.Nu={-(1-t),1-t,t,-t};a.Nv={-(1-s),-s,s,1-s};a.Nuv={1,-1,1,-1};
        Vecd e=x[1]-x[0],f=x[3]-x[0],warp=(x[2]-x[3])-e;a.x=x[0]+s*e+t*f+s*t*warp;a.u=e+t*warp;a.v=f+s*warp;a.uv=warp;return a;
    }
    static Param closest(const std::array<Vecd,4>&global,const Vecd&point,Param p=Param(.5,.5));
};
/** One surface finite element whose corner DOFs are one-sided bulk traces.
 * All cohesive quadrature and contact geometry use the SAME isoparametric map.
 * This avoids inventing a gap between coincident curved surfaces by measuring
 * one high-order trace against an unrelated flat chord plane.
 */
class BilinearCrackSurface:public ReferenceSolidSurface {
    VisibleCrackSurface&corner;StdVec<std::array<Real,4>>coefficient;size_t generation=0,owner_version=0;
 public:
    using Param=BilinearCrackGeometry::Param;StdVec<Param>natural;
    explicit BilinearCrackSurface(VisibleCrackSurface&source,const StdVec<Vecd>&quadrature)
      :ReferenceSolidSurface(source.port,{},1,DeferredBuildTag{}),corner(source){
        if(source.points.size()!=4)throw std::invalid_argument("Crack surface needs four material corner traces");
        std::array<Vecd,4>x;for(size_t i=0;i<4;++i){x[i]=source.points[i];std::array<Real,4>n{};n[i]=1;coefficient.push_back(n);points.push_back(x[i]);}
        StdVec<std::array<Real,4>>normal_probes;
        for(const auto&q:quadrature){auto p=BilinearCrackGeometry::closest(x,q);auto a=BilinearCrackGeometry::evaluate(x,p);Real scale=std::max((x[1]-x[0]).norm(),(x[3]-x[0]).norm());
            if(p.minCoeff()< -1e-8||p.maxCoeff()>1+1e-8||(a.x-q).norm()>1e-9*scale)throw std::runtime_error("Crack quadrature is outside its reference isoparametric patch");
            natural.push_back(p);coefficient.push_back(a.N);points.push_back(a.x);
            auto u=a.N,v=a.N;for(size_t i=0;i<4;++i){u[i]+=a.Nu[i];v[i]+=a.Nv[i];}normal_probes.push_back(u);normal_probes.push_back(v);
        }
        for(const auto&c:normal_probes){coefficient.push_back(c);Vecd q=Vecd::Zero();for(size_t i=0;i<4;++i)q+=c[i]*x[i];points.push_back(q);}rebuild();
    }
    void rebuild(){corner.ensureCurrent();rows.clear();for(const auto&c:coefficient){std::map<size_t,Real>sum;for(size_t i=0;i<4;++i)for(size_t a=0;a<corner.rows[i].index.size();++a)sum[corner.rows[i].index[a]]+=c[i]*corner.rows[i].weight[a];Row row;for(auto e:sum)if(e.second!=0){row.index.push_back(e.first);row.weight.push_back(e.second);}rows.push_back(std::move(row));}owner_version=corner.mappingVersion();++generation;}
    void rebuildChanged(){corner.rebuildChanged();rebuild();}
    void captureTransferValues(){corner.captureTransferValues();}
    const StdVec<Vecd>&transferValues()const{return corner.transferValues();}
    const StdVec<Vecd>&transferVelocities()const{return corner.transferVelocities();}
    void restoreTransferValues(const StdVec<Vecd>&v,const StdVec<Vecd>&velocity={}){corner.restoreTransferValues(v,velocity);}
    void ensureCurrent()const override{corner.ensureCurrent();if(owner_version!=corner.mappingVersion())throw std::runtime_error("Isoparametric crack map is stale");}
    size_t mappingVersion()const override{return generation;}
    size_t normalU(size_t q)const{return 4+natural.size()+2*q;}
    size_t normalV(size_t q)const{return normalU(q)+1;}
};

/** Symmetric two-pass frictionless contact on the ACTUAL crack surface bank.
 * Each pass has half the reference area, preserving the input penalty K.
 * A current spatial hash and outward-normal test permit sliding onto adjacent
 * crack patches without relying on global fragment IDs (a crack may still be
 * connected around its tip). Closest projection is onto a bilinear master.
 * The exact condensed distance Hessian includes the parameter/normal motion.
 * No node deletion, friction spring, force clipping or residual cohesion.
 */
class CrackSurfaceContactNetwork {
    CrackProjectionBatch *projection_batch_;
    std::vector<CrackProjectionMaster> projection_masters_;
    std::vector<CrackProjectionPoint> projection_points_;
    std::vector<Vecd> projection_velocities_;
    std::vector<size_t> projection_point_offsets_;
    std::vector<CrackProjectionPair> projection_pairs_;
    std::vector<CrackProjectionResult> projection_results_;

    struct ContactPoint
    {
        size_t slave, quadrature, master;
        Vecd position, velocity, normal;
        BilinearCrackGeometry::Value geometry;
        Real gap, stiffness;
    };
    std::vector<ContactPoint> selected_contacts_;
    std::vector<CrackLocalContactInput> contact_inputs_;
    std::vector<CrackLocalContactResult> contact_results_;
    std::uint64_t trace_generation_ = 0;
    std::vector<size_t> trace_side_offsets_;
    std::vector<std::pair<const BilinearCrackSurface *, size_t>> trace_versions_;
    std::vector<CrackTraceAddress> contact_addresses_;
    std::vector<CrackTraceParticleState> trace_particle_input_;
    std::vector<CrackTraceParticleResult> trace_particle_output_;
    void updateDeviceTraceMap();
    void applyDeviceTracePullback(SolidForcePort &port);

 public:
    struct Side {BilinearCrackSurface*surface=nullptr;const StdVec<Real>*areas=nullptr;Real outward_sign=1;size_t reference_patch=~size_t(0);};
    StdVec<Side>sides;Real penalty,spacing;
    Real energy=0,max_penetration=0,normal_force_sum=0,stable_dt=std::numeric_limits<Real>::infinity();
    Real force_residual=0,torque_residual=0,power_residual=0,tangential_force_sum=0;
    size_t active_contacts=0,outside_projections=0,neighbor_contacts=0;StdVec<Vecd>rates;
    size_t device_contact_evaluations=0,host_contact_fallbacks=0;
    bool device_trace_pullback=false;
    explicit CrackSurfaceContactNetwork(const StdVec<Side>&s,Real K,Real dp,CrackProjectionBatch *projection_batch=nullptr):projection_batch_(projection_batch),sides(s),penalty(K),spacing(dp){if(K<=0||dp<=0||s.empty())throw std::invalid_argument("Invalid crack-contact network");for(const auto&a:s)if(!a.surface||!a.areas||a.areas->size()!=a.surface->natural.size()||std::abs(a.outward_sign)!=1||a.surface->port.position!=s[0].surface->port.position||a.surface->port.mass!=s[0].surface->port.mass)throw std::invalid_argument("Crack-contact sides must belong to one physical material port");rates.resize(s[0].surface->port.count);}
    void apply();
};
}
