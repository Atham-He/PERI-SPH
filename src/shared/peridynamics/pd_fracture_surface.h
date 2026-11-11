#pragma once
#include "pd_solid_reference_coupling.h"
#include "pd_fracture_correspondence.h"
#include <queue>

namespace SPH::peridynamics
{
/** A surface quadrature point keeps its original material anchor. Its support
 * is restricted to the anchor's live-bond component WITHIN the local support.
 * Lost/unsupported owners produce an eroded surface entry, never a spring
 * across a crack. Quadratic -> full-rank affine fallback is explicit. */
class FractureAwareSurface : public ReferenceSolidSurface
{
    solid_dynamics::PDFractureReference &reference;
    Real dp;StdVec<Row> original;StdVec<size_t> anchor;StdVec<StdVec<size_t>> affected;
    StdVec<size_t> membership,visited;size_t token=0;
  public:
    StdVec<int> row_order;
    size_t materialAnchor(size_t q)const{return anchor.at(q);}
    size_t inactive_points=0,affine_points=0;
    FractureAwareSurface(SolidForcePort &p,const StdVec<Vecd> &targets,Real spacing,Real support,unsigned order,
                         solid_dynamics::PDFractureReference &r)
      :ReferenceSolidSurface(p,targets,spacing,support,order),reference(r),dp(spacing),original(rows),
       affected(p.count),membership(p.count,0),visited(p.count,0)
    {
        anchor.resize(points.size());row_order.assign(points.size(),order);
        for(size_t q=0;q<points.size();++q)
        {Real best=std::numeric_limits<Real>::infinity();for(size_t j:original[q].index)
         {affected[j].push_back(q);Real d=(port.reference[j]-points[q]).squaredNorm();if(d<best){best=d;anchor[q]=j;}}}
    }
    void rebuildChanged()
    {
        StdVec<bool> update(points.size(),false);
        for(size_t j:reference.changed_particles)for(size_t q:affected[j])update[q]=true;
        for(size_t q=0;q<points.size();++q)if(update[q])rebuildRow(q);
        inactive_points=affine_points=0;
        for(int order:row_order){inactive_points+=order==0;affine_points+=order==1;}
    }
  private:
    void rebuildRow(size_t q)
    {
        ++token;if(token==0)throw std::runtime_error("Surface visibility stamp overflow");
        for(size_t j:original[q].index)membership[j]=token;
        Row result;StdVec<size_t> queue;
        if(reference.degree[anchor[q]]>0){queue.push_back(anchor[q]);visited[anchor[q]]=token;}
        for(size_t k=0;k<queue.size();++k)
        {
            size_t i=queue[k];const auto &nb=reference.neighbors(i);
            for(size_t n=0;n<nb.current_size_;++n)
            {size_t j=nb.j_[n];if(nb.bondLive_[n]&&membership[j]==token&&visited[j]!=token){visited[j]=token;queue.push_back(j);}}
        }
        // Keep the original order to preserve deterministic quadrature sums.
        for(size_t j:original[q].index)if(visited[j]==token)result.index.push_back(j);
        if(reproduction_order==2&&fit(q,result,2)){rows[q]=std::move(result);row_order[q]=2;return;}
        result.weight.clear();if(fit(q,result,1)){rows[q]=std::move(result);row_order[q]=1;return;}
        rows[q]=Row{};row_order[q]=0;
    }
    bool fit(size_t q,Row &r,unsigned order)
    {
        size_t dim=order==2?10:4;if(r.index.size()<dim)return false;
        StdVec<Vecd> z;StdVec<Real> w;Real sum=0;Matd cov=Matd::Zero();
        for(size_t j:r.index){Vecd a=(port.reference[j]-points[q])/dp;Real k=std::exp(-a.squaredNorm())*port.volume[j]/std::pow(dp,3);
            z.push_back(a);w.push_back(k);sum+=k;cov+=k*a*a.transpose();}
        cov/=sum;Eigen::SelfAdjointEigenSolver<Matd> eig(cov);
        if(eig.info()!=Eigen::Success||eig.eigenvalues().minCoeff()<=1e-6*eig.eigenvalues().maxCoeff())return false;
        Matd L=eig.eigenvectors()*eig.eigenvalues().cwiseSqrt().cwiseInverse().asDiagonal()*eig.eigenvectors().transpose();
        Eigen::MatrixXd phi(r.index.size(),dim);Eigen::MatrixXd moment=Eigen::MatrixXd::Zero(dim,dim);
        for(size_t k=0;k<z.size();++k)
        {Vecd a=L*z[k];phi(k,0)=1;for(int d=0;d<3;++d)phi(k,d+1)=a[d];
         if(order==2){phi(k,4)=a[0]*a[0];phi(k,5)=a[1]*a[1];phi(k,6)=a[2]*a[2];phi(k,7)=a[0]*a[1];phi(k,8)=a[0]*a[2];phi(k,9)=a[1]*a[2];}
         moment+=w[k]*phi.row(k).transpose()*phi.row(k);}
        Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> me(moment);
        if(me.info()!=Eigen::Success||me.eigenvalues().minCoeff()<=(order==2?1e-8:1e-6)*me.eigenvalues().maxCoeff())return false;
        Eigen::VectorXd e0=Eigen::VectorXd::Zero(dim);e0[0]=1;Eigen::VectorXd c=moment.ldlt().solve(e0),reproduced=Eigen::VectorXd::Zero(dim);
        if(!c.allFinite())return false;r.weight.clear();
        for(size_t k=0;k<z.size();++k){Real a=w[k]*phi.row(k).dot(c);r.weight.push_back(a);reproduced+=a*phi.row(k).transpose();}
        return (reproduced-e0).cwiseAbs().maxCoeff()<1e-8;
    }
};

/** Plate contact at core locations exposed by erosion of the shell map. */
class ExposedCorePlaneContact
{
    ReferenceSolidSurface &core,&shell;
    StdVec<Real> stiffness;
  public:
    Real stable_dt,upper=0,lower=0,energy=0;
    ExposedCorePlaneContact(ReferenceSolidSurface &c,ReferenceSolidSurface &s,const StdVec<Real>&areas,Real K)
      :core(c),shell(s)
    {if(c.rows.size()!=s.rows.size()||areas.size()!=c.rows.size())throw std::invalid_argument("Exposed contact dimensions");
     for(Real a:areas)stiffness.push_back(a*K);stable_dt=ReferenceSolidSurface::stableBound(core,nullptr,stiffness);}
    void apply(Real bottom,Real top)
    {
        upper=lower=energy=0;
        for(size_t q=0;q<stiffness.size();++q)if(shell.rows[q].index.empty())
        {Vecd x=core.position(q);Real lo=std::max(bottom-x[2],Real(0)),hi=std::max(x[2]-top,Real(0));
         core.scatter(q,Vecd(0,0,stiffness[q]*(lo-hi)));upper+=stiffness[q]*hi;lower+=stiffness[q]*lo;
         energy+=.5*stiffness[q]*(lo*lo+hi*hi);}
    }
};
}
