#pragma once
#include "pd_fracture_visible_surface.h"
#include <array>
#include <unordered_map>
#include <map>
namespace SPH::peridynamics {
struct BilinearCrackGeometry {
    using Param=Eigen::Vector2d;
    struct Value {std::array<Real,4>N,Nu,Nv,Nuv;Vecd x=Vecd::Zero(),u=Vecd::Zero(),v=Vecd::Zero(),uv=Vecd::Zero();};
    static Value evaluate(const std::array<Vecd,4>&x,const Param&p){Real s=p[0],t=p[1];Value a;
        a.N={(1-s)*(1-t),s*(1-t),s*t,(1-s)*t};a.Nu={-(1-t),1-t,t,-t};a.Nv={-(1-s),-s,s,1-s};a.Nuv={1,-1,1,-1};
        for(size_t i=0;i<4;++i){a.x+=a.N[i]*x[i];a.u+=a.Nu[i]*x[i];a.v+=a.Nv[i]*x[i];a.uv+=a.Nuv[i]*x[i];}return a;
    }
    static Param closest(const std::array<Vecd,4>&x,const Vecd&q,Param p=Param(.5,.5)){
        for(size_t iteration=0;iteration<40;++iteration){auto a=evaluate(x,p);Vecd r=q-a.x;Eigen::Matrix2d H;H<<a.u.squaredNorm(),a.u.dot(a.v)-r.dot(a.uv),a.u.dot(a.v)-r.dot(a.uv),a.v.squaredNorm();
            Real scale=std::max(a.u.squaredNorm(),a.v.squaredNorm());if(!std::isfinite(scale)||scale<=0)throw std::runtime_error("Crack contact surface has no geometric rank");
            Param gradient(a.u.dot(r),a.v.dot(r));Eigen::SelfAdjointEigenSolver<Eigen::Matrix2d> e(H);
            if(e.info()!=Eigen::Success)throw std::runtime_error("Crack projection factorization");
            if(e.eigenvalues().minCoeff()<=1e-13*scale){H<<a.u.squaredNorm(),a.u.dot(a.v),a.u.dot(a.v),a.v.squaredNorm();}
            if(H.determinant()<=1e-14*scale*scale)throw std::runtime_error("Crack contact patch is geometrically unresolved");
            Param step=H.ldlt().solve(gradient);if(!step.allFinite())throw std::runtime_error("Nonfinite closest crack projection");
            if(step.norm()<2e-12)return p;
            Real old=r.squaredNorm(),alpha=1;bool accepted=false;
            for(size_t line=0;line<25;++line){Param trial=p+alpha*step;Real value=(q-evaluate(x,trial).x).squaredNorm();if(value<=old*(1+1e-14)){p=trial;accepted=true;break;}alpha*=.5;}
            if(!accepted){if(gradient.norm()<1e-12*scale)return p;throw std::runtime_error("Closest crack projection failed its descent check");}
        }
        throw std::runtime_error("Closest crack-surface projection did not converge");
    }
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
 public:
    struct Side {BilinearCrackSurface*surface=nullptr;const StdVec<Real>*areas=nullptr;Real outward_sign=1;size_t reference_patch=~size_t(0);};
    StdVec<Side>sides;Real penalty,spacing;
    Real energy=0,max_penetration=0,normal_force_sum=0,stable_dt=std::numeric_limits<Real>::infinity();
    Real force_residual=0,torque_residual=0,power_residual=0,tangential_force_sum=0;
    size_t active_contacts=0,outside_projections=0,neighbor_contacts=0;StdVec<Vecd>rates;
    explicit CrackSurfaceContactNetwork(const StdVec<Side>&s,Real K,Real dp):sides(s),penalty(K),spacing(dp){if(K<=0||dp<=0||s.empty())throw std::invalid_argument("Invalid crack-contact network");for(const auto&a:s)if(!a.surface||!a.areas||a.areas->size()!=a.surface->natural.size()||std::abs(a.outward_sign)!=1||a.surface->port.position!=s[0].surface->port.position||a.surface->port.mass!=s[0].surface->port.mass)throw std::invalid_argument("Crack-contact sides must belong to one physical material port");rates.resize(s[0].surface->port.count);}
    void apply(){
        energy=max_penetration=normal_force_sum=force_residual=torque_residual=power_residual=tangential_force_sum=0;active_contacts=outside_projections=neighbor_contacts=0;for(auto&r:rates)r.setZero();
        struct Master {std::array<Vecd,4>x,v;};StdVec<Master>masters(sides.size());
        using Cell=std::array<long long,3>;struct Hash{size_t operator()(const Cell&a)const{size_t h=0;for(auto v:a)h^=std::hash<long long>{}(v)+0x9e3779b9+(h<<6)+(h>>2);return h;}};
        std::unordered_map<Cell,StdVec<size_t>,Hash>grid;Real cell_size=2*spacing;
        auto cell=[&](const Vecd&x){return Cell{(long long)std::floor(x[0]/cell_size),(long long)std::floor(x[1]/cell_size),(long long)std::floor(x[2]/cell_size)};};
        for(size_t m=0;m<sides.size();++m){auto&s=*sides[m].surface;s.ensureCurrent();Vecd lo=Vecd::Constant(std::numeric_limits<Real>::infinity()),hi=-lo;
            for(size_t a=0;a<4;++a){s.kinematics(a,masters[m].x[a],masters[m].v[a]);lo=lo.cwiseMin(masters[m].x[a]);hi=hi.cwiseMax(masters[m].x[a]);}Real margin=std::max(2*spacing,(hi-lo).norm());Cell a=cell(lo-Vecd::Constant(margin)),b=cell(hi+Vecd::Constant(margin));
            for(long long i=a[0];i<=b[0];++i)for(long long j=a[1];j<=b[1];++j)for(long long k=a[2];k<=b[2];++k)grid[{i,j,k}].push_back(m);
        }
        Vecd balance=Vecd::Zero(),moment=Vecd::Zero();Real magnitude=0,length=spacing,power=0,power_scale=0,energy_rate=0;
        for(size_t si=0;si<sides.size();++si){auto&slave=*sides[si].surface;for(size_t q=0;q<slave.natural.size();++q){Vecd x,v;slave.kinematics(4+q,x,v);Vecd tangentU=slave.position(slave.normalU(q))-x,tangentV=slave.position(slave.normalV(q))-x;Vecd outward=sides[si].outward_sign*tangentU.cross(tangentV).normalized();
                if(!outward.allFinite())throw std::runtime_error("Slave crack normal is unresolved");auto it=grid.find(cell(x));if(it==grid.end())continue;
                size_t chosen=sides.size();Real distance=std::numeric_limits<Real>::infinity(),gap=0;BilinearCrackGeometry::Value selected;Vecd normal;
                for(size_t m:it->second){if(m==si)continue;auto p=BilinearCrackGeometry::closest(masters[m].x,x);if(p.minCoeff()< -1e-9||p.maxCoeff()>1+1e-9){++outside_projections;continue;}auto a=BilinearCrackGeometry::evaluate(masters[m].x,p);Vecd n=sides[m].outward_sign*a.u.cross(a.v).normalized();if(!n.allFinite())throw std::runtime_error("Master crack normal is unresolved");if(n.dot(outward)>=-1e-8)continue;Real g=n.dot(x-a.x);if(std::abs(g)<distance){distance=std::abs(g);chosen=m;gap=g;selected=a;normal=n;}}
                if(chosen==sides.size())continue;auto&master=*sides[chosen].surface;const auto&a=selected;Real stiffness=.5*penalty*sides[si].areas->at(q);Vecd residual=Vecd::Zero();if(gap<0)residual=gap*normal;
                Eigen::Matrix2d Hp;Hp<<a.u.squaredNorm(),a.u.dot(a.v)-residual.dot(a.uv),a.u.dot(a.v)-residual.dot(a.uv),a.v.squaredNorm();Eigen::SelfAdjointEigenSolver<Eigen::Matrix2d>e(Hp);if(e.info()!=Eigen::Success||e.eigenvalues().minCoeff()<=1e-13*e.eigenvalues().maxCoeff())throw std::runtime_error("Closest crack contact tangent is not a local minimum");
                Eigen::Matrix<Real,15,2>coupling;std::array<Real,5>c{1,-a.N[0],-a.N[1],-a.N[2],-a.N[3]};
                for(size_t i=0;i<5;++i){Vecd u=-c[i]*a.u,vv=-c[i]*a.v;if(i){u-=a.Nu[i-1]*residual;vv-=a.Nv[i-1]*residual;}coupling.template block<3,1>(3*i,0)=u;coupling.template block<3,1>(3*i,1)=vv;}
                Eigen::Matrix<Real,15,15>H=-coupling*Hp.inverse()*coupling.transpose();for(size_t i=0;i<5;++i)for(size_t j=0;j<5;++j)H.template block<3,3>(3*i,3*j)+=c[i]*c[j]*Matd::Identity();H*=stiffness;
                std::array<const ReferenceSolidSurface::Row*,5>trace{&slave.rows[4+q],&master.rows[0],&master.rows[1],&master.rows[2],&master.rows[3]};std::array<Real,5>weighted{};
                for(size_t i=0;i<5;++i)for(size_t j=0;j<trace[i]->index.size();++j)weighted[i]+=std::abs(trace[i]->weight[j])/std::sqrt(slave.port.mass[trace[i]->index[j]]);
                for(size_t i=0;i<5;++i)for(size_t d=0;d<3;++d){Real row=0;for(size_t j=0;j<15;++j)row+=std::abs(H(3*i+d,j))*weighted[j/3];for(size_t j=0;j<trace[i]->index.size();++j)rates[trace[i]->index[j]][d]+=std::abs(trace[i]->weight[j])*row/std::sqrt(slave.port.mass[trace[i]->index[j]]);}
                if(gap>=0)continue;if(sides[si].reference_patch!=sides[chosen].reference_patch)++neighbor_contacts;Vecd f=-stiffness*gap*normal;slave.scatter(4+q,f);Vecd relative=v;
                balance+=f;moment+=x.cross(f);magnitude+=f.norm();length=std::max(length,x.norm());power+=f.dot(v);power_scale+=std::abs(f.dot(v));
                for(size_t j=0;j<4;++j){Vecd reaction=-a.N[j]*f;master.scatter(j,reaction);balance+=reaction;moment+=masters[chosen].x[j].cross(reaction);magnitude+=reaction.norm();power+=reaction.dot(masters[chosen].v[j]);power_scale+=std::abs(reaction.dot(masters[chosen].v[j]));relative-=a.N[j]*masters[chosen].v[j];}
                energy_rate+=stiffness*gap*normal.dot(relative);energy+=.5*stiffness*gap*gap;max_penetration=std::max(max_penetration,-gap);normal_force_sum+=f.norm();tangential_force_sum+=(f-normal*f.dot(normal)).norm();++active_contacts;
            }}
        Real bound=0;for(const auto&r:rates)bound=std::max(bound,r.maxCoeff());stable_dt=bound>0?.2/std::sqrt(bound):std::numeric_limits<Real>::infinity();force_residual=balance.norm()/std::max(magnitude,Real(1e-30));torque_residual=moment.norm()/std::max(magnitude*length,Real(1e-30));power_residual=std::abs(power+energy_rate)/std::max(power_scale,Real(1e-30));
    }
};
}
