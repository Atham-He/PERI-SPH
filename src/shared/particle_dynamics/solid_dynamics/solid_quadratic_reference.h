/* Optional quadratic-consistent total-Lagrangian solid-SPH reference operator.
 * The reference kernel and neighborhood are retained. Forces are the exact
 * transpose of the reference gradient plus the derivative of a nonnegative
 * weighted polynomial-residual energy. Constitutive updates are independent.
 * No pressure/force offset, clipping, particle shifting or fitted multiplier.
 */
#pragma once
#include "elastic_dynamics.h"
#include "base_particles.hpp"
#include <Eigen/Cholesky>
#include <Eigen/Eigenvalues>
#include <algorithm>
#include <vector>

namespace SPH::solid_dynamics
{
class SolidQuadraticReference : public LocalDynamics, public DataDelegateInner
{
  public:
    using V9=Eigen::Matrix<Real,9,1>;
    using M9=Eigen::Matrix<Real,9,9>;
    using A39=Eigen::Matrix<Real,3,9>;
    struct Edge {size_t j,reverse;Vecd gradient;V9 basis,coefficient;Real weight;Vecd residual,stabilizing_force;};
    std::vector<std::vector<Edge>> rows;
    StdVec<Real> length2,energy,scale;
    Real max_condition=0,max_reproduction_error=0,stable_dt=0;
    Real *volume,*mass;Vecd *position;
    Real shear,elastic_norm;
    SolidQuadraticReference(BaseInnerRelation &relation,Real shear_modulus,Real elastic_operator_norm)
      :LocalDynamics(relation.getSPHBody()),DataDelegateInner(relation),
       volume(particles_->getVariableDataByName<Real>("VolumetricMeasure")),
       mass(particles_->getVariableDataByName<Real>("Mass")),
       position(particles_->getVariableDataByName<Vecd>("Position")),
       shear(shear_modulus),elastic_norm(elastic_operator_norm)
    {if(!(shear>0&&elastic_norm>0))throw std::invalid_argument("Positive solid elastic bounds required");}
    void build()
    {
        static_assert(Dimensions==3,"Quadratic solid SPH currently supports 3-D");
        size_t count=particles_->TotalRealParticles();rows.resize(count);length2.resize(count);energy.resize(count);scale.resize(count);
        for(size_t i=0;i<count;++i)
        {
            const auto &nb=inner_configuration_[i];Real sum=0;Matd covariance=Matd::Zero();
            for(size_t k=0;k<nb.current_size_;++k)
            {Vecd xi=position[nb.j_[k]]-position[i];Real w=-nb.dW_ij_[k]*volume[nb.j_[k]]/nb.r_ij_[k];
             if(!(w>0))throw std::runtime_error("Positive reference radial-kernel weights required");sum+=w;covariance+=w*xi*xi.transpose();}
            if(!(sum>0))throw std::runtime_error("Empty quadratic solid support");covariance/=sum;
            Eigen::SelfAdjointEigenSolver<Matd> eig(covariance);
            if(eig.info()!=Eigen::Success||eig.eigenvalues().minCoeff()<=0)throw std::runtime_error("Invalid quadratic solid covariance");
            Matd L=eig.eigenvectors()*eig.eigenvalues().cwiseSqrt().cwiseInverse().asDiagonal()*eig.eigenvectors().transpose();
            auto basis=[&](const Vecd&xi)->V9{Vecd z=L*xi;V9 b;b<<z[0],z[1],z[2],z[0]*z[0],z[1]*z[1],z[2]*z[2],z[0]*z[1],z[0]*z[2],z[1]*z[2];return b;};
            M9 moment=M9::Zero();
            for(size_t k=0;k<nb.current_size_;++k)
            {V9 b=basis(position[nb.j_[k]]-position[i]);Real w=-nb.dW_ij_[k]*volume[nb.j_[k]]/(nb.r_ij_[k]*sum);moment+=w*b*b.transpose();}
            Eigen::SelfAdjointEigenSolver<M9> meig(moment);
            if(meig.info()!=Eigen::Success||meig.eigenvalues().minCoeff()<=1e-12*meig.eigenvalues().maxCoeff())throw std::runtime_error("Rank-deficient quadratic solid reference support");
            max_condition=std::max(max_condition,meig.eigenvalues().maxCoeff()/meig.eigenvalues().minCoeff());
            Eigen::LDLT<M9> solve(moment);Eigen::Matrix<Real,9,3> reproduced=Eigen::Matrix<Real,9,3>::Zero(),expected=Eigen::Matrix<Real,9,3>::Zero();expected.topRows<3>()=L;
            rows[i].reserve(nb.current_size_);length2[i]=covariance.trace();scale[i]=shear/length2[i];
            for(size_t k=0;k<nb.current_size_;++k)
            {size_t j=nb.j_[k];V9 b=basis(position[j]-position[i]);Real w=-nb.dW_ij_[k]*volume[j]/(nb.r_ij_[k]*sum);V9 c=w*solve.solve(b);Vecd g=L.transpose()*c.head<3>();
             rows[i].push_back({j,0,g,b,c,w,Vecd::Zero(),Vecd::Zero()});reproduced+=b*g.transpose();}
            max_reproduction_error=std::max(max_reproduction_error,(reproduced-expected).norm()/expected.norm());
        }
        if(max_reproduction_error>1e-9)throw std::runtime_error("Quadratic solid reproduction failed");
        StdVec<Real> bound(count,0),hg_bound(count,0);
        for(size_t i=0;i<count;++i)
        {
            Vecd self=Vecd::Zero();Real norms=0;
            for(auto &e:rows[i])
            {auto &other=rows[e.j];auto it=std::find_if(other.begin(),other.end(),[&](const Edge&r){return r.j==i;});
             if(it==other.end())throw std::runtime_error("Asymmetric total-Lagrangian solid support");e.reverse=size_t(it-other.begin());self-=e.gradient;norms+=e.gradient.norm();
             Real k=volume[i]*scale[i]*e.weight;hg_bound[i]+=2*k;hg_bound[e.j]+=2*k;}
            norms+=self.norm();Real factor=volume[i]*elastic_norm*norms;bound[i]+=factor*self.norm();for(const auto&e:rows[i])bound[e.j]+=factor*e.gradient.norm();
        }
        Real maximum=0;for(size_t i=0;i<count;++i)maximum=std::max(maximum,(bound[i]+hg_bound[i])/mass[i]);stable_dt=.2/std::sqrt(maximum);
    }
    Matd gradient(size_t i,const Vecd*field)const
    {Matd g=Matd::Zero();for(const auto&e:rows[i])g+=(field[e.j]-field[i])*e.gradient.transpose();return g;}
    void prepare(size_t i,const Vecd*field)
    {
        A39 fit=A39::Zero();for(const auto&e:rows[i])fit+=(field[e.j]-field[i])*e.coefficient.transpose();
        A39 orthogonality=A39::Zero();energy[i]=0;
        for(auto&e:rows[i]){e.residual=field[e.j]-field[i]-fit*e.basis;orthogonality+=e.weight*e.residual*e.basis.transpose();energy[i]+=.5*volume[i]*scale[i]*e.weight*e.residual.squaredNorm();}
        // Retain the full residual-projector transpose even in floating point.
        for(auto&e:rows[i])e.stabilizing_force=volume[i]*scale[i]*(e.weight*e.residual-orthogonality*e.coefficient);
    }
    Vecd force(size_t i,const Matd*P)const
    {Vecd f=Vecd::Zero();for(const auto&e:rows[i]){const auto&r=rows[e.j][e.reverse];f+=volume[i]*P[i]*e.gradient-volume[e.j]*P[e.j]*r.gradient+e.stabilizing_force-r.stabilizing_force;}return f;}
    Real totalEnergy()const {Real sum=0;for(Real e:energy)sum+=e;return sum;}
};
} // namespace SPH::solid_dynamics
