/* Licensed under the Apache License, Version 2.0. */
#ifndef PD_SOLID_REFERENCE_COUPLING_H
#define PD_SOLID_REFERENCE_COUPLING_H

#include "base_particles.hpp"
#include "base_body.h"
#include <Eigen/Cholesky>
#include <Eigen/Eigenvalues>
#include <array>
#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <unordered_map>

namespace SPH::peridynamics
{
/** A force port, not a constitutive adapter. PD stores external acceleration;
 * solid SPH stores external force. Both ports receive physical forces in N.
 * The caller owns clearing and the common integration time layer. */
class SolidForcePort
{
  public:
    enum class Storage { ForcePrior, PDAccelerationPrior };
    RealBody &body;
    size_t count;
    Vecd *position, *velocity, *external;
    Real *mass, *volume;
    StdVec<Vecd> reference;
    Storage storage;
    SolidForcePort(RealBody &b, Storage s) : body(b), storage(s)
    {
        auto &p=b.getBaseParticles();count=p.TotalRealParticles();
        position=p.getVariableDataByName<Vecd>("Position");
        velocity=p.getVariableDataByName<Vecd>("Velocity");
        mass=p.getVariableDataByName<Real>("Mass");
        volume=p.getVariableDataByName<Real>("VolumetricMeasure");
        external=p.registerStateVariableData<Vecd>(s==Storage::ForcePrior ? "ForcePrior" : "PDAccelerationPrior");
        reference.assign(position,position+count);
        for(size_t i=0;i<count;++i)
            if(!std::isfinite(mass[i]) || mass[i]<=0 || !std::isfinite(volume[i]) || volume[i]<=0)
                throw std::invalid_argument("Solid coupling requires positive physical masses and volumes");
    }
    void clear() {for(size_t i=0;i<count;++i) external[i].setZero();}
    void addForce(size_t i,const Vecd &force)
    {external[i] += storage==Storage::ForcePrior ? force : force/mass[i];}
};

/** Fixed-reference polynomial-reproducing MLS surface evaluation. It supports
 * nonuniform particle volumes, including explicitly resolved thickness layers.
 * It supplies no bending law and does not turn a volume body into a shell. */
class ReferenceSolidSurface
{
  public:
    struct Row {StdVec<size_t> index;StdVec<Real> weight;};
    SolidForcePort &port;
    StdVec<Vecd> points;
    StdVec<Row> rows;
    Real affine_error=0, polynomial_error=0, largest_condition_number=0;
    const unsigned reproduction_order;
    ReferenceSolidSurface(SolidForcePort &p,const StdVec<Vecd> &targets,Real spacing,Real support_ratio=2.5,unsigned order=1)
        : port(p),points(targets),reproduction_order(order)
    {
        static_assert(Dimensions==3,"ReferenceSolidSurface currently implements three-dimensional affine reconstruction");
        if(spacing<=0 || support_ratio<=1 || !std::isfinite(spacing*support_ratio))
            throw std::invalid_argument("Invalid reference-surface support");
        if(order!=1 && order!=2)throw std::invalid_argument("Reference surface supports polynomial order 1 or 2");
        const Real support=spacing*support_ratio;
        using Cell=std::array<long long,3>;
        struct Hash {size_t operator()(const Cell &a) const
            {size_t h=0;for(auto x:a) h^=std::hash<long long>{}(x)+0x9e3779b9+(h<<6)+(h>>2);return h;}};
        std::unordered_map<Cell,StdVec<size_t>,Hash> cells;
        auto cell=[&](const Vecd &x){return Cell{(long long)std::floor(x[0]/support),
            (long long)std::floor(x[1]/support),(long long)std::floor(x[2]/support)};};
        for(size_t j=0;j<port.count;++j) cells[cell(port.reference[j])].push_back(j);
        rows.reserve(points.size());
        for(const auto &target:points)
        {
            Row row;StdVec<Eigen::Vector4d> basis;StdVec<Real> kernel;
            Eigen::Matrix4d moment=Eigen::Matrix4d::Zero();auto center=cell(target);
            for(int a=-1;a<=1;++a) for(int b=-1;b<=1;++b) for(int c=-1;c<=1;++c)
            {
                auto found=cells.find(Cell{center[0]+a,center[1]+b,center[2]+c});
                if(found==cells.end()) continue;
                for(size_t j:found->second)
                {
                    Vecd q=(port.reference[j]-target)/spacing;
                    if(q.squaredNorm()>=support_ratio*support_ratio) continue;
                    Eigen::Vector4d phi;phi<<1,q[0],q[1],q[2];
                    Real w=std::exp(-q.squaredNorm())*port.volume[j]/std::pow(spacing,3);
                    moment.noalias()+=w*phi*phi.transpose();
                    row.index.push_back(j);basis.push_back(phi);kernel.push_back(w);
                }
            }
            if(order==2)
            {
                // Degree two requires at least three resolved normal levels on
                // a planar volume support. Reject rank loss rather than adding
                // a diagonal regularizer that silently changes reproduction.
                using V10=Eigen::Matrix<Real,10,1>;using M10=Eigen::Matrix<Real,10,10>;
                Real total=0;Matd covariance=Matd::Zero();
                for(size_t a=0;a<row.index.size();++a)
                {Vecd q=(port.reference[row.index[a]]-target)/spacing;
                 total+=kernel[a];covariance+=kernel[a]*q*q.transpose();}
                if(row.index.size()<10 || total<=0)throw std::runtime_error("Insufficient quadratic surface support");
                covariance/=total;Eigen::SelfAdjointEigenSolver<Matd> ceig(covariance);
                if(ceig.info()!=Eigen::Success || ceig.eigenvalues().minCoeff()<=0)
                    throw std::runtime_error("Invalid quadratic surface covariance");
                Matd L=ceig.eigenvectors()*ceig.eigenvalues().cwiseSqrt().cwiseInverse().asDiagonal()*ceig.eigenvectors().transpose();
                auto basis2=[&](size_t j)->V10 {Vecd z=L*((port.reference[j]-target)/spacing);
                    V10 phi;phi<<1,z[0],z[1],z[2],z[0]*z[0],z[1]*z[1],z[2]*z[2],z[0]*z[1],z[0]*z[2],z[1]*z[2];return phi;};
                M10 m=M10::Zero();
                for(size_t a=0;a<row.index.size();++a){V10 phi=basis2(row.index[a]);m+=kernel[a]*phi*phi.transpose();}
                Eigen::SelfAdjointEigenSolver<M10> eig2(m);
                if(eig2.info()!=Eigen::Success || eig2.eigenvalues().minCoeff()<=1e-12*eig2.eigenvalues().maxCoeff())
                    throw std::runtime_error("Rank-deficient quadratic volume-to-surface mapping; increase resolved normal support");
                largest_condition_number=std::max(largest_condition_number,eig2.eigenvalues().maxCoeff()/eig2.eigenvalues().minCoeff());
                V10 e0=V10::Zero();e0[0]=1;V10 coeff=m.ldlt().solve(e0),reproduced_basis=V10::Zero();
                if(!coeff.allFinite())throw std::runtime_error("Non-finite quadratic reference surface coefficients");
                Real sum=0;Vecd reproduced=Vecd::Zero();
                for(size_t a=0;a<row.index.size();++a)
                {V10 phi=basis2(row.index[a]);Real w=kernel[a]*phi.dot(coeff);row.weight.push_back(w);
                 sum+=w;reproduced+=w*port.reference[row.index[a]];reproduced_basis+=w*phi;}
                affine_error=std::max({affine_error,std::abs(sum-1),(reproduced-target).norm()/spacing});
                polynomial_error=std::max(polynomial_error,(reproduced_basis-e0).cwiseAbs().maxCoeff());
                if(!std::isfinite(polynomial_error) || polynomial_error>1e-8 || !std::isfinite(affine_error) || affine_error>1e-8)
                    throw std::runtime_error("Quadratic reference surface reproduction failed");
                rows.push_back(std::move(row));continue;
            }
            // Preserve the legacy affine arithmetic and default behavior.
            Eigen::SelfAdjointEigenSolver<Eigen::Matrix4d> eig(moment);
            if(row.index.size()<8 || eig.info()!=Eigen::Success ||
                eig.eigenvalues().minCoeff()<=1e-12*eig.eigenvalues().maxCoeff())
                throw std::runtime_error("Rank-deficient volume-to-surface mapping; resolve physical thickness/support");
            largest_condition_number=std::max(largest_condition_number,
                eig.eigenvalues().maxCoeff()/eig.eigenvalues().minCoeff());
            Eigen::Vector4d e0(1,0,0,0),coeff=moment.ldlt().solve(e0);
            Real sum=0;Vecd reproduced=Vecd::Zero();
            for(size_t a=0;a<row.index.size();++a)
            {Real w=kernel[a]*basis[a].dot(coeff);row.weight.push_back(w);sum+=w;reproduced+=w*port.reference[row.index[a]];}
            affine_error=std::max({affine_error,std::abs(sum-1),(reproduced-target).norm()/spacing});
            if(!std::isfinite(affine_error) || affine_error>1e-8)
                throw std::runtime_error("Reference surface does not reproduce affine motion");
            rows.push_back(std::move(row));
        }
    }
    Vecd position(size_t i) const {Vecd x=Vecd::Zero();const auto &r=rows[i];
        for(size_t a=0;a<r.index.size();++a)x+=r.weight[a]*port.position[r.index[a]];return x;}
    Vecd velocity(size_t i) const {Vecd v=Vecd::Zero();const auto &r=rows[i];
        for(size_t a=0;a<r.index.size();++a)v+=r.weight[a]*port.velocity[r.index[a]];return v;}
    void scatter(size_t i,const Vecd &f) {const auto &r=rows[i];
        for(size_t a=0;a<r.index.size();++a)port.addForce(r.index[a],r.weight[a]*f);}

    /** Absolute row-sum bound of sqrt(K) B M^-1 B^T sqrt(K). */
    static Real stableBound(const ReferenceSolidSurface &a,const ReferenceSolidSurface *b,
                            const StdVec<Real> &stiffness,Real safety=.2)
    {
        if(stiffness.size()!=a.rows.size() || (b && b->rows.size()!=a.rows.size()))
            throw std::invalid_argument("Surface constraint dimensions disagree");
        StdVec<Real> sa(a.port.count,0),sb(b ? b->port.count:0,0);
        for(size_t i=0;i<stiffness.size();++i)
        {
            if(!std::isfinite(stiffness[i]) || stiffness[i]<=0)throw std::invalid_argument("Positive constraint stiffness required");
            Real k=std::sqrt(stiffness[i]);const auto &r=a.rows[i];
            for(size_t j=0;j<r.index.size();++j)sa[r.index[j]]+=k*std::abs(r.weight[j]);
            if(b){const auto &s=b->rows[i];for(size_t j=0;j<s.index.size();++j)sb[s.index[j]]+=k*std::abs(s.weight[j]);}
        }
        Real largest=0;
        for(size_t i=0;i<stiffness.size();++i)
        {
            Real rate=0,k=std::sqrt(stiffness[i]);const auto &r=a.rows[i];
            for(size_t j=0;j<r.index.size();++j)rate+=k*std::abs(r.weight[j])*sa[r.index[j]]/a.port.mass[r.index[j]];
            if(b){const auto &s=b->rows[i];for(size_t j=0;j<s.index.size();++j)rate+=k*std::abs(s.weight[j])*sb[s.index[j]]/b->port.mass[s.index[j]];}
            largest=std::max(largest,rate);
        }
        return safety/std::sqrt(largest);
    }
};

/** One reference interface energy:U=sum(k_i|W_a*x_a-W_b*x_b|^2/2).
 * Both bodies receive the transpose pullback of the SAME physical force.
 * No fluid pressure, velocity overwrite or duplicate repulsion is added. */
class BondedSolidInterface
{
  public:
    ReferenceSolidSurface &a,&b;
    StdVec<Real> stiffness;
    StdVec<Vecd> x_a,x_b,force_a;
    Real stable_dt,energy=0,max_gap=0,force_residual=0,torque_residual=0,power_residual=0;
    Real worst_force=0,worst_torque=0,worst_power=0,length_scale;
    BondedSolidInterface(ReferenceSolidSurface &first,ReferenceSolidSurface &second,
                         const StdVec<Real> &areas,Real stiffness_per_area,Real length)
        : a(first),b(second),length_scale(length)
    {
        if(a.rows.size()!=b.rows.size() || a.points.size()!=areas.size() || length<=0)
            throw std::invalid_argument("Invalid bonded-solid interface dimensions");
        for(size_t i=0;i<areas.size();++i)
        {if((a.points[i]-b.points[i]).norm()>1e-9*length || areas[i]<=0)
            throw std::invalid_argument("Bonded-solid reference surfaces do not coincide");
         stiffness.push_back(areas[i]*stiffness_per_area);}
        x_a.resize(areas.size());x_b.resize(areas.size());force_a.resize(areas.size());
        stable_dt=ReferenceSolidSurface::stableBound(a,&b,stiffness);
    }
    void apply()
    {
        energy=max_gap=0;Vecd sum=Vecd::Zero(),torque=Vecd::Zero();Real power=0,fs=0,ps=0;
        for(size_t i=0;i<stiffness.size();++i)
        {
            Vecd xa=a.position(i),xb=b.position(i),va=a.velocity(i),vb=b.velocity(i);
            Vecd gap=xa-xb,f=-stiffness[i]*gap;
            x_a[i]=xa;x_b[i]=xb;force_a[i]=f;
            a.scatter(i,f);b.scatter(i,-f);
            energy+=.5*stiffness[i]*gap.squaredNorm();max_gap=std::max(max_gap,gap.norm());
            for(int side=0;side<2;++side)
            {
                const auto &map=side==0?a:b;const auto &row=map.rows[i];Real sign=side==0?1:-1;
                for(size_t j=0;j<row.index.size();++j)
                {
                    Vecd fj=sign*row.weight[j]*f;size_t id=row.index[j];
                    sum+=fj;torque+=map.port.position[id].cross(fj);
                    Real p=fj.dot(map.port.velocity[id]);power+=p;fs+=fj.norm();ps+=std::abs(p);
                }
            }
            Real energy_rate=-f.dot(va-vb);power+=energy_rate;ps+=std::abs(energy_rate);
        }
        force_residual=sum.norm()/std::max(fs,Real(1e-30));
        torque_residual=torque.norm()/std::max(length_scale*fs,Real(1e-30));
        power_residual=std::abs(power)/std::max(ps,Real(1e-30));
        worst_force=std::max(worst_force,force_residual);worst_torque=std::max(worst_torque,torque_residual);worst_power=std::max(worst_power,power_residual);
    }
};

/** Frictionless analytic planes acting on a reconstructed physical surface. */
class SolidSurfacePlaneContact
{
  public:
    ReferenceSolidSurface &surface;
    StdVec<Real> areas,stiffness;
    StdVec<Vecd> position,force;
    Real stable_dt,upper=0,lower=0,energy=0,max_penetration=0;
    Real worst_force=0,worst_torque=0,worst_power=0,length_scale;
    SolidSurfacePlaneContact(ReferenceSolidSurface &s,const StdVec<Real> &area,Real stiffness_per_area,Real length)
        : surface(s),areas(area),length_scale(length)
    {
        if(areas.size()!=s.rows.size())throw std::invalid_argument("Contact area count differs");
        for(Real a:areas){if(a<=0)throw std::invalid_argument("Contact areas must be positive");stiffness.push_back(a*stiffness_per_area);}
        position.resize(area.size());force.resize(area.size());
        stable_dt=ReferenceSolidSurface::stableBound(surface,nullptr,stiffness);
    }
    void apply(Real bottom,Real top)
    {
        upper=lower=energy=max_penetration=0;Vecd ferr=Vecd::Zero(),terr=Vecd::Zero();Real perr=0,fs=0,ps=0;
        for(size_t i=0;i<stiffness.size();++i)
        {
            Vecd x=surface.position(i),v=surface.velocity(i);
            Real pb=std::max(bottom-x[2],Real(0)),pt=std::max(x[2]-top,Real(0));
            Vecd f(0,0,stiffness[i]*(pb-pt));position[i]=x;force[i]=f;
            upper+=stiffness[i]*pt;lower+=stiffness[i]*pb;
            energy+=.5*stiffness[i]*(pb*pb+pt*pt);max_penetration=std::max({max_penetration,pb,pt});
            surface.scatter(i,f);ferr-=f;terr-=x.cross(f);perr-=f.dot(v);ps+=std::abs(f.dot(v));
            const auto &r=surface.rows[i];
            for(size_t j=0;j<r.index.size();++j)
            {size_t id=r.index[j];Vecd fj=r.weight[j]*f;ferr+=fj;terr+=surface.port.position[id].cross(fj);
             Real p=fj.dot(surface.port.velocity[id]);perr+=p;fs+=fj.norm();ps+=std::abs(p);}
        }
        worst_force=std::max(worst_force,ferr.norm()/std::max(fs,Real(1e-30)));
        worst_torque=std::max(worst_torque,terr.norm()/std::max(length_scale*fs,Real(1e-30)));
        worst_power=std::max(worst_power,std::abs(perr)/std::max(ps,Real(1e-30)));
    }
};
} // namespace SPH::peridynamics
#endif
