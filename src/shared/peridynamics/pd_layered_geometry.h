#pragma once
#include "base_particles.hpp"
#include <Eigen/Eigenvalues>
#include <algorithm>
#include <limits>
#include <stdexcept>

namespace SPH::peridynamics
{
/** Explicit reference thickness columns of a volume quadrature. Directors are
 * reference geometric data, not current-position constraints. No particles,
 * masses, material constants or physical thicknesses are modified. */
struct LayeredParticleGeometry
{
    StdVec<StdVec<size_t>> columns;
    StdVec<Vecd> director,reference;
    StdVec<Real> xi,thickness;
    StdVec<StdVec<Real>> derivative;
    StdVec<size_t> column_of,layer_of;
    LayeredParticleGeometry(BaseParticles&p,const StdVec<StdVec<size_t>>&c,const StdVec<Vecd>&n,const StdVec<Real>&nodes)
        :columns(c),director(n),xi(nodes)
    {
        const size_t count=p.TotalRealParticles(),L=xi.size();
        if(L<3||columns.empty()||director.size()!=columns.size())throw std::invalid_argument("Layered geometry needs >=3 levels and explicit columns/directors");
        for(size_t l=0;l<L;++l)if(!std::isfinite(xi[l])||(l&&xi[l]<=xi[l-1]))throw std::invalid_argument("Invalid normal quadrature coordinates");
        auto*x=p.getVariableDataByName<Vecd>("Position");reference.assign(x,x+count);
        column_of.assign(count,std::numeric_limits<size_t>::max());layer_of=column_of;thickness.resize(columns.size());
        for(size_t k=0;k<columns.size();++k)
        {
            if(columns[k].size()!=L||!director[k].allFinite()||std::abs(director[k].norm()-1)>1e-10)throw std::invalid_argument("Invalid thickness column/director");
            for(size_t l=0;l<L;++l){size_t i=columns[k][l];if(i>=count||column_of[i]!=std::numeric_limits<size_t>::max())throw std::invalid_argument("Duplicate/missing thickness-column particle");column_of[i]=k;layer_of[i]=l;}
            Vecd d=(x[columns[k].back()]-x[columns[k].front()])/(xi.back()-xi.front());
            Real half=d.dot(director[k]);if(half<=0||(d-half*director[k]).norm()>1e-9*half)throw std::invalid_argument("Thickness columns must be straight and director-aligned");
            thickness[k]=2*half;
            for(size_t l=0;l<L;++l)if((x[columns[k][l]]-x[columns[k][0]]-(xi[l]-xi[0])*d).norm()>1e-9*half)throw std::invalid_argument("Inconsistent normal quadrature coordinates");
        }
        for(size_t i:column_of)if(i==std::numeric_limits<size_t>::max())throw std::invalid_argument("Every volume particle must belong to one thickness column");
        StdVec<Real>bary(L,1);for(size_t l=0;l<L;++l)for(size_t m=0;m<L;++m)if(m!=l)bary[l]/=xi[l]-xi[m];
        derivative.assign(L,StdVec<Real>(L,0));
        for(size_t l=0;l<L;++l)for(size_t m=0;m<L;++m)if(m!=l){derivative[l][m]=bary[m]/(bary[l]*(xi[l]-xi[m]));derivative[l][l]-=derivative[l][m];}
    }
    StdVec<Real> trace(Real z)const
    {
        if(!std::isfinite(z)||z< -1||z>1)throw std::invalid_argument("Physical normal trace coordinate must lie in [-1,1]");
        StdVec<Real>w(xi.size(),1);for(size_t l=0;l<xi.size();++l)for(size_t m=0;m<xi.size();++m)if(m!=l)w[l]*=(z-xi[m])/(xi[l]-xi[m]);return w;
    }
};

namespace layered_detail
{
template<int N>struct MomentInverse
{
    using Matrix=Eigen::Matrix<Real,N,N>;
    Matrix inverse;Real condition;int rank=0;
    explicit MomentInverse(const Matrix&M)
    {
        Eigen::SelfAdjointEigenSolver<Matrix> e(M);
        if(e.info()!=Eigen::Success||e.eigenvalues().maxCoeff()<=0||e.eigenvalues().minCoeff()< -1e-10*e.eigenvalues().maxCoeff())throw std::runtime_error("Invalid surface polynomial moment");
        auto values=e.eigenvalues();Real largest=values.maxCoeff(),smallest=largest;
        for(int k=0;k<N;++k)if(values[k]>1e-12*largest){smallest=std::min(smallest,values[k]);values[k]=1/values[k];++rank;}else values[k]=0;
        inverse=e.eigenvectors()*values.asDiagonal()*e.eigenvectors().transpose();condition=largest/smallest;
    }
};
inline Matd whitening(const Matd&covariance,Real spacing)
{
    Eigen::SelfAdjointEigenSolver<Matd> e(covariance);
    if(e.info()!=Eigen::Success||e.eigenvalues().maxCoeff()<=0)throw std::runtime_error("Invalid surface covariance");
    Vecd scale=e.eigenvalues();Real largest=scale.maxCoeff();
    for(int k=0;k<3;++k)scale[k]=scale[k]>1e-12*largest?1/std::sqrt(scale[k]):1/spacing;
    return e.eigenvectors()*scale.asDiagonal()*e.eigenvectors().transpose();
}
inline Eigen::Matrix<Real,9,1> basis9(const Vecd&z)
{Eigen::Matrix<Real,9,1>b;b<<z[0],z[1],z[2],z[0]*z[0],z[1]*z[1],z[2]*z[2],z[0]*z[1],z[0]*z[2],z[1]*z[2];return b;}
inline Eigen::Matrix<Real,10,1> basis10(const Vecd&z)
{Eigen::Matrix<Real,10,1>b;b<<1,basis9(z);return b;}
// The nullspace here is analytic surface-polynomial redundancy (e.g. the
// sphere equation), not lost physical bonds. Every requested functional must
// lie in the represented range; no diagonal regularization or silent fallback.
}
}
