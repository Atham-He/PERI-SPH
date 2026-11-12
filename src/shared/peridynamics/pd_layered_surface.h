#pragma once
#include "pd_solid_reference_coupling.h"
#include "pd_layered_geometry.h"

namespace SPH::peridynamics
{
/** Exact normal polynomial trace followed by a quadratic surface MLS value.
 * The flattened rows still use the existing force/velocity transpose ports.
 * No extra surface particles, penalty correction, or force scale is added. */
class LayeredSolidSurface : public ReferenceSolidSurface
{
  public:
    LayeredSolidSurface(SolidForcePort&p,const StdVec<Vecd>&targets,const LayeredParticleGeometry&g,Real normal_coordinate,Real spacing,Real support_ratio=2.5)
        :ReferenceSolidSurface(p,targets,2,DeferredBuildTag{})
    {
        using namespace layered_detail;using V10=Eigen::Matrix<Real,10,1>;using M10=Eigen::Matrix<Real,10,10>;
        if(spacing<=0||support_ratio<=1)throw std::invalid_argument("Invalid layered surface support");
        auto trace=g.trace(normal_coordinate);StdVec<Vecd>surface(g.columns.size(),Vecd::Zero());StdVec<Real>measure(g.columns.size(),0);
        for(size_t c=0;c<g.columns.size();++c)for(size_t l=0;l<trace.size();++l){size_t j=g.columns[c][l];surface[c]+=trace[l]*p.reference[j];measure[c]+=p.volume[j];}
        Real support=spacing*support_ratio;using Cell=std::array<long long,3>;
        struct Hash{size_t operator()(const Cell&a)const{size_t h=0;for(auto x:a)h^=std::hash<long long>{}(x)+0x9e3779b9+(h<<6)+(h>>2);return h;}};
        auto cell=[&](const Vecd&x){return Cell{(long long)std::floor(x[0]/support),(long long)std::floor(x[1]/support),(long long)std::floor(x[2]/support)};};
        std::unordered_map<Cell,StdVec<size_t>,Hash>cells;for(size_t c=0;c<surface.size();++c)cells[cell(surface[c])].push_back(c);
        rows.reserve(targets.size());
        for(const auto&target:targets)
        {
            StdVec<size_t>columns;StdVec<Real>kernel;Matd covariance=Matd::Zero();Real total=0;auto center=cell(target);
            for(int a=-1;a<=1;++a)for(int b=-1;b<=1;++b)for(int c=-1;c<=1;++c)
            {auto found=cells.find(Cell{center[0]+a,center[1]+b,center[2]+c});if(found==cells.end())continue;
             for(size_t j:found->second){Vecd q=(surface[j]-target)/spacing;if(q.squaredNorm()>=support_ratio*support_ratio)continue;Real w=std::exp(-q.squaredNorm())*measure[j]/std::pow(spacing,3);columns.push_back(j);kernel.push_back(w);total+=w;covariance+=w*q*q.transpose();}}
            if(columns.size()<6||total<=0)throw std::runtime_error("Insufficient layered surface support");
            Matd L=whitening(covariance/total,1);M10 M=M10::Zero();
            for(size_t a=0;a<columns.size();++a){V10 b=basis10(L*((surface[columns[a]]-target)/spacing));M+=kernel[a]/total*b*b.transpose();}
            MomentInverse<10> inverse(M);largest_condition_number=std::max(largest_condition_number,inverse.condition);V10 e0=V10::Zero();e0[0]=1;
            if(inverse.rank<6||(M*inverse.inverse*e0-e0).norm()>1e-9)throw std::runtime_error("Surface value functional is not represented by layered support");
            V10 coefficient=inverse.inverse*e0;Row row;Real sum=0;Vecd reproduced=Vecd::Zero();V10 poly=V10::Zero();
            for(size_t a=0;a<columns.size();++a)
            {size_t c=columns[a];V10 b=basis10(L*((surface[c]-target)/spacing));Real angular=kernel[a]/total*b.dot(coefficient);
             for(size_t l=0;l<trace.size();++l){size_t j=g.columns[c][l];Real w=angular*trace[l];row.index.push_back(j);row.weight.push_back(w);sum+=w;reproduced+=w*p.reference[j];poly+=w*basis10((p.reference[j]-target)/spacing);}}
            affine_error=std::max({affine_error,std::abs(sum-1),(reproduced-target).norm()/spacing});polynomial_error=std::max(polynomial_error,(poly-e0).cwiseAbs().maxCoeff());
            if(affine_error>1e-8||polynomial_error>1e-8)throw std::runtime_error("Layered surface quadratic reproduction failed");rows.push_back(std::move(row));
        }
    }
};
}
