/* Optional 3-D linear reproducing-kernel basis; no constitutive parameters. */
#pragma once
#include "base_data_type.h"
#include <Eigen/Cholesky>
#include <algorithm>
#include <vector>
#include <stdexcept>
namespace SPH::solid_dynamics
{
class SnsniRKCloud
{
 public:
    struct Entry {size_t j;Real w;Vecd gradient;};
    StdVec<Vecd> X;StdVec<Real> measure;
    Real dp,support;Vecd origin;Eigen::Vector3i dimensions;
    std::vector<std::vector<size_t>> bins;
    SnsniRKCloud(const StdVec<Vecd>&points,Real spacing,const StdVec<Real>&volumes={},Real radius=2.5)
        :X(points),measure(volumes),dp(spacing),support(radius)
    {
        static_assert(Dimensions==3,"S-NSNI implementation is three-dimensional");
        if(X.empty()||dp<=0||support<=1)throw std::invalid_argument("Invalid RK cloud");
        origin=X[0];Vecd hi=X[0];for(auto&x:X){origin=origin.cwiseMin(x);hi=hi.cwiseMax(x);}
        dimensions=((hi-origin)/dp).array().floor().cast<int>().matrix()+Eigen::Vector3i::Constant(2);
        bins.resize(size_t(dimensions[0])*dimensions[1]*dimensions[2]);
        if(measure.empty())measure.assign(X.size(),1);if(measure.size()!=X.size())throw std::invalid_argument("RK volume dimensions");
        for(size_t j=0;j<X.size();++j){if(measure[j]<=0)throw std::invalid_argument("Positive RK sampling weights required");auto c=cell(X[j]);bins[key(c)].push_back(j);}
    }
    Eigen::Vector3i cell(const Vecd&x)const{return ((x-origin)/dp).array().floor().cast<int>();}
    size_t key(const Eigen::Vector3i&c)const{return (size_t(c[0])*dimensions[1]+c[1])*dimensions[2]+c[2];}
    static Real kernel(Real r){r=std::abs(r);return r<.5?2./3.-4*r*r+4*r*r*r:(r<1?4./3.*std::pow(1-r,3):0);}
    std::vector<Entry> shape(const Vecd&x)const
    {
        using V4=Eigen::Vector4d;using M4=Eigen::Matrix4d;
        Eigen::Vector3i low=cell(x-Vecd::Constant(support*dp)).cwiseMax(0);
        Eigen::Vector3i high=cell(x+Vecd::Constant(support*dp)).cwiseMin(dimensions-Eigen::Vector3i::Ones());
        std::vector<size_t> candidates;for(int a=low[0];a<=high[0];++a)for(int b=low[1];b<=high[1];++b)for(int c=low[2];c<=high[2];++c){auto&v=bins[key(Eigen::Vector3i(a,b,c))];candidates.insert(candidates.end(),v.begin(),v.end());}
        std::sort(candidates.begin(),candidates.end());
        std::vector<size_t> ids;std::vector<V4> basis;std::vector<Real>w;M4 M=M4::Zero();
        for(auto j:candidates){Vecd z=(X[j]-x)/(support*dp);Real q=measure[j]*kernel(z[0])*kernel(z[1])*kernel(z[2]);if(q<=0)continue;V4 p;p<<1,z;M+=q*p*p.transpose();ids.push_back(j);basis.push_back(p);w.push_back(q);}
        Eigen::LDLT<M4> inverse(M);if(inverse.info()!=Eigen::Success||inverse.vectorD().minCoeff()<=0)throw std::runtime_error("Rank-deficient RK support");
        M4 rhs=M4::Identity();rhs.bottomRightCorner<3,3>()/=(support*dp);M4 coeff=inverse.solve(rhs);
        std::vector<Entry> out;out.reserve(ids.size());for(size_t j=0;j<ids.size();++j){V4 p=w[j]*coeff.transpose()*basis[j];out.push_back({ids[j],p[0],p.tail<3>()});}return out;
    }
};
}
