/* Licensed under the Apache License, Version 2.0. */
#ifndef PD_FIXED_REFERENCE_H
#define PD_FIXED_REFERENCE_H

#include "data_type.h"
#include <Eigen/Eigenvalues>
#include <algorithm>
#include <cstdint>
#include <stdexcept>

namespace SPH::peridynamics
{
/** Actual bounded total-Lagrangian correspondence reference, not an affine
 * tetra surrogate. Every directed neighbor packet, shape inverse, gradient
 * and nonaffine quadratic cache is persistent restart state. Construction is
 * fresh-reference only; restoration must replace stored arrays, not rebuild.
 */
struct FixedPDReference
{
    uint64_t generation=1;
    Real horizon=0;
    StdVec<Vecd> position,shape_gradient;
    StdVec<Matd> inverse_shape;
    StdVec<Real> volume,weight,neighbor_measure,nonaffine_unit;
    StdVec<size_t> offsets{0},target;
    StdVec<uint8_t> active;
    StdVec<Real> shape_eigenvalue_ratio;

    FixedPDReference()=default; // Empty archive receiver; no reference refresh.
    FixedPDReference(StdVec<Vecd> reference,StdVec<Real> volumes,Real delta)
        :horizon(delta),position(std::move(reference)),volume(std::move(volumes))
    {
        const size_t n=position.size();
        if(n<4 || n>128 || volume.size()!=n || !std::isfinite(delta) || delta<=0)
            throw std::invalid_argument("fixed PD reference: bounded dimensions/horizon");
        for(size_t i=0;i<n;++i)
            if(!position[i].allFinite() || !std::isfinite(volume[i]) || volume[i]<=0)
                throw std::invalid_argument("fixed PD reference: invalid position/volume");
        shape_gradient.assign(n*n,Vecd::Zero());inverse_shape.resize(n);neighbor_measure.resize(n);
        shape_eigenvalue_ratio.resize(n);nonaffine_unit.assign(n*n,0);
        for(size_t i=0;i<n;++i)
        {
            Matd shape=Matd::Zero();Real measure=0;
            for(size_t j=0;j<n;++j)if(i!=j)
            {
                const Vecd xi=position[j]-position[i];const Real radius2=xi.squaredNorm();
                if(radius2==0)throw std::invalid_argument("fixed PD reference: duplicate particle");
                if(radius2>=delta*delta)continue;
                const Real w=std::pow(1-radius2/(delta*delta),3);
                target.push_back(j);weight.push_back(w);active.push_back(1);
                shape+=w*volume[j]*xi*xi.transpose();measure+=w*volume[j];
            }
            offsets.push_back(target.size());
            Eigen::SelfAdjointEigenSolver<Matd> spectrum(shape);
            if(spectrum.info()!=Eigen::Success || !(spectrum.eigenvalues()[2]>0) ||
                spectrum.eigenvalues()[0]<=1e-12*spectrum.eigenvalues()[2] || !std::isfinite(measure) || measure<=0)
                throw std::invalid_argument("fixed PD reference: rank deficient shape");
            inverse_shape[i]=shape.ldlt().solve(Matd::Identity());neighbor_measure[i]=measure;
            shape_eigenvalue_ratio[i]=spectrum.eigenvalues()[0]/spectrum.eigenvalues()[2];
            for(size_t p=offsets[i];p<offsets[i+1];++p)
            {
                const size_t j=target[p];
                shape_gradient[i*n+j]=weight[p]*volume[j]*inverse_shape[i]*(position[j]-position[i]);
                shape_gradient[i*n+i]-=shape_gradient[i*n+j];
            }
        }
        // H_ij q=(q_j-q_i)-F_i xi_ij; this scalar row is shared by all3
        // physical components. The fixed scalar cache L=sum a_ij H^T H,
        // multiplied later by eta*mu, is positive semidefinite. Division by
        // neighbor measure gives the stabilization coefficient correct units.
        for(size_t i=0;i<n;++i)for(size_t p=offsets[i];p<offsets[i+1];++p)
        {
            const size_t j=target[p];const Vecd xi=position[j]-position[i];
            const Real coefficient=volume[i]*weight[p]*volume[j]/(delta*delta*neighbor_measure[i]);
            StdVec<Real> h(n);
            for(size_t k=0;k<n;++k)h[k]=Real(k==j)-Real(k==i)-shape_gradient[i*n+k].dot(xi);
            for(size_t a=0;a<n;++a)for(size_t b=0;b<n;++b)nonaffine_unit[a*n+b]+=coefficient*h[a]*h[b];
        }
        validate();
    }

    size_t size() const{return position.size();}
    Matd deformation(size_t i,const StdVec<Vecd> &current) const
    {
        const size_t n=size();if(i>=n || current.size()!=n)throw std::invalid_argument("fixed PD gradient: dimensions");
        Matd f;for(size_t a=0;a<3;++a)for(size_t b=0;b<3;++b)
        {
            long double sum=0;for(size_t k=0;k<n;++k)sum+=static_cast<long double>(current[k][a])*shape_gradient[i*n+k][b];
            f(a,b)=sum;
        }
        if(!f.allFinite())throw std::runtime_error("fixed PD gradient: nonfinite deformation");return f;
    }
    /** No cache repair: structurally and kinematically validate stored data.
     * Independent qualification also reconstructs every cached number. */
    void validate() const
    {
        const size_t n=size();
        if(generation!=1 || n<4 || n>128 || volume.size()!=n || inverse_shape.size()!=n ||
            shape_gradient.size()!=n*n || nonaffine_unit.size()!=n*n || neighbor_measure.size()!=n || shape_eigenvalue_ratio.size()!=n ||
            offsets.size()!=n+1 || offsets.front()!=0 || offsets.back()!=target.size() || weight.size()!=target.size() || active.size()!=target.size() ||
            !std::isfinite(horizon) || horizon<=0)throw std::invalid_argument("fixed PD reference: incompatible stored dimensions/generation");
        StdVec<uint8_t> reached(n,0);reached[0]=1;
        for(size_t sweep=0;sweep<n;++sweep)for(size_t i=0;i<n;++i)
        {
            if(offsets[i]>offsets[i+1] || offsets[i+1]>target.size())throw std::invalid_argument("fixed PD reference: offsets");
            for(size_t p=offsets[i];p<offsets[i+1];++p)
            {
                const size_t j=target[p];
                if(j>=n || j==i || (p>offsets[i] && j<=target[p-1]) || active[p]!=1 || !std::isfinite(weight[p]) || weight[p]<=0)
                    throw std::invalid_argument("fixed PD reference: invalid directed packet");
                if(reached[i])reached[j]=1;
            }
        }
        if(std::find(reached.begin(),reached.end(),uint8_t(0))!=reached.end())throw std::invalid_argument("fixed PD reference: disconnected graph");
        for(size_t i=0;i<n;++i)
        {
            if(!position[i].allFinite() || !inverse_shape[i].allFinite() || !std::isfinite(volume[i]) || volume[i]<=0 ||
                !std::isfinite(neighbor_measure[i]) || neighbor_measure[i]<=0 || !std::isfinite(shape_eigenvalue_ratio[i]) || shape_eigenvalue_ratio[i]<=1e-12)
                throw std::invalid_argument("fixed PD reference: invalid stored geometry/shape");
            Vecd partition=Vecd::Zero();
            for(size_t j=0;j<n;++j)
            {
                if(!shape_gradient[i*n+j].allFinite() || !std::isfinite(nonaffine_unit[i*n+j]))throw std::invalid_argument("fixed PD reference: nonfinite cache");
                partition+=shape_gradient[i*n+j];
                if(std::abs(nonaffine_unit[i*n+j]-nonaffine_unit[j*n+i])>1e-12*(1+std::abs(nonaffine_unit[i*n+j])))
                    throw std::invalid_argument("fixed PD reference: asymmetric quadratic cache");
            }
            if(partition.norm()*horizon>1e-10 || (deformation(i,position)-Matd::Identity()).norm()>1e-10)
                throw std::invalid_argument("fixed PD reference: affine reproduction");
            for(size_t p=offsets[i];p<offsets[i+1];++p)
            {
                const size_t j=target[p];const auto begin=target.begin()+offsets[j],end=target.begin()+offsets[j+1];
                if(!std::binary_search(begin,end,i))throw std::invalid_argument("fixed PD reference: missing reciprocal packet");
            }
            for(size_t j=0;j<n;++j)if(i!=j)
            {
                const bool required=(position[j]-position[i]).squaredNorm()<horizon*horizon;
                const bool present=std::binary_search(target.begin()+offsets[i],target.begin()+offsets[i+1],j);
                if(required!=present)throw std::invalid_argument("fixed PD reference: incomplete intact neighborhood");
            }
        }
        // Read-only validation of ALL stored numerical caches. Reconstructing
        // expected values in scratch is permitted; replacing a restored cache
        // or refreshing an accepted reference is deliberately unavailable.
        StdVec<Real> expected_l(n*n,0);
        for(size_t i=0;i<n;++i)
        {
            Matd shape=Matd::Zero();Real measure=0;StdVec<Vecd> expected_g(n,Vecd::Zero());
            for(size_t p=offsets[i];p<offsets[i+1];++p)
            {
                const size_t j=target[p];const Vecd xi=position[j]-position[i];
                const Real radius2=xi.squaredNorm(),expected_weight=std::pow(1-radius2/(horizon*horizon),3);
                if(!(radius2>0 && radius2<horizon*horizon) || std::abs(weight[p]-expected_weight)>1e-12*(1+std::abs(expected_weight)))
                    throw std::invalid_argument("fixed PD reference: stored kernel mismatch");
                shape+=weight[p]*volume[j]*xi*xi.transpose();measure+=weight[p]*volume[j];
                expected_g[j]=weight[p]*volume[j]*inverse_shape[i]*xi;expected_g[i]-=expected_g[j];
            }
            Eigen::SelfAdjointEigenSolver<Matd> spectrum(shape);
            if(spectrum.info()!=Eigen::Success || spectrum.eigenvalues()[0]<=1e-12*spectrum.eigenvalues()[2] ||
                std::abs(spectrum.eigenvalues()[0]/spectrum.eigenvalues()[2]-shape_eigenvalue_ratio[i])>1e-10 ||
                (shape*inverse_shape[i]-Matd::Identity()).norm()>1e-10 || std::abs(measure-neighbor_measure[i])>1e-12*measure)
                throw std::invalid_argument("fixed PD reference: stored inverse/measure mismatch");
            for(size_t j=0;j<n;++j)if((expected_g[j]-shape_gradient[i*n+j]).norm()*horizon>1e-10)
                throw std::invalid_argument("fixed PD reference: stored gradient mismatch");
            for(size_t p=offsets[i];p<offsets[i+1];++p)
            {
                const size_t j=target[p];const Vecd xi=position[j]-position[i];
                const Real coefficient=volume[i]*weight[p]*volume[j]/(horizon*horizon*neighbor_measure[i]);StdVec<Real> h(n);
                for(size_t k=0;k<n;++k)h[k]=Real(k==j)-Real(k==i)-shape_gradient[i*n+k].dot(xi);
                for(size_t a=0;a<n;++a)for(size_t b=0;b<n;++b)expected_l[a*n+b]+=coefficient*h[a]*h[b];
            }
        }
        for(size_t i=0;i<n*n;++i)if(std::abs(nonaffine_unit[i]-expected_l[i])>1e-12*(1+std::abs(expected_l[i])))
            throw std::invalid_argument("fixed PD reference: stored nonaffine cache mismatch");
    }

    template<class Checkpoint>void registerCheckpoint(Checkpoint &cp)
    {
        cp.addScalar("PDReferenceGeneration",generation);cp.addScalar("PDReferenceHorizon",horizon);
        cp.addArray("PDReferencePosition",position);cp.addArray("PDReferenceVolume",volume);
        cp.addArray("PDReferenceOffsets",offsets);cp.addArray("PDReferenceTarget",target);cp.addArray("PDReferenceWeight",weight);
        cp.addArray("PDReferenceActive",active);cp.addArray("PDShapeInverse",inverse_shape);cp.addArray("PDShapeGradient",shape_gradient);
        cp.addArray("PDNeighborMeasure",neighbor_measure);cp.addArray("PDShapeEigenvalueRatio",shape_eigenvalue_ratio);
        cp.addArray("PDNonaffineUnitCache",nonaffine_unit);
    }
};

inline Real nonaffineEnergy(const FixedPDReference &ref,const StdVec<Vecd> &q,Real eta_mu)
{
    if(q.size()!=ref.size() || !std::isfinite(eta_mu) || eta_mu<0)throw std::invalid_argument("PD stabilization: dimensions/coefficient");
    for(const auto &x:q)if(!x.allFinite())throw std::invalid_argument("PD stabilization: nonfinite state");
    // Evaluate positive residual squares directly, not a cancellation-prone
    // q^T L q in a nearly affine configuration. The force still uses Lq.
    long double energy=0;
    for(size_t i=0;i<ref.size();++i)
    {
        const Matd f=ref.deformation(i,q);
        for(size_t p=ref.offsets[i];p<ref.offsets[i+1];++p)
        {
            const size_t j=ref.target[p];const Vecd r=q[j]-q[i]-f*(ref.position[j]-ref.position[i]);
            energy+=static_cast<long double>(eta_mu)*ref.volume[i]*ref.weight[p]*ref.volume[j]*r.squaredNorm()/
                (2*ref.horizon*ref.horizon*ref.neighbor_measure[i]);
        }
    }
    if(!std::isfinite(energy) || !std::isfinite(Real(energy)))throw std::runtime_error("PD stabilization: nonfinite energy");return energy;
}

inline StdVec<Vecd> nonaffineForce(const FixedPDReference &ref,const StdVec<Vecd> &q,Real eta_mu)
{
    if(q.size()!=ref.size() || !std::isfinite(eta_mu) || eta_mu<0)throw std::invalid_argument("PD stabilization: dimensions/coefficient");
    for(const auto &x:q)if(!x.allFinite())throw std::invalid_argument("PD stabilization: nonfinite state");
    StdVec<Vecd> result(q.size(),Vecd::Zero());
    for(size_t i=0;i<q.size();++i)for(size_t a=0;a<3;++a)
    {
        long double sum=0;for(size_t j=0;j<q.size();++j)sum-=static_cast<long double>(eta_mu)*ref.nonaffine_unit[i*q.size()+j]*q[j][a];
        result[i][a]=sum;
        if(!std::isfinite(result[i][a]))throw std::runtime_error("PD stabilization: nonfinite force");
    }
    return result;
}

/** Canonical initial33point full sphere, explicit normalized numerical
 * geometry. This is not an experimental diameter or a scalable mesher. */
inline StdVec<Vecd> boundedSphereReference()
{
    StdVec<Vecd> points;
    for(int i=-2;i<=2;++i)for(int j=-2;j<=2;++j)for(int k=-2;k<=2;++k)
        if(i*i+j*j+k*k<=4)points.emplace_back(.5*i,.5*j,.5*k);
    return points;
}
} // namespace SPH::peridynamics
#endif
