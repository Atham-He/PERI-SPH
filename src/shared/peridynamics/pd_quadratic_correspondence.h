#pragma once
/** Optional quadratic-reproducing NOSB correspondence for intact 3-D solids.
 * Reference bond differences, full degree-two reproduction and the EXACT
 * transpose force operator are used together. No material constants or
 * thicknesses are rescaled. The legacy first-order operator is unchanged.
 * Bond failure requires rebuilding the reference fit and is rejected here.
 */
#include "nosbpd_dynamics.h"
#include <Eigen/Cholesky>
#include <Eigen/Eigenvalues>
#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <vector>

namespace SPH { namespace solid_dynamics {
class PDQuadraticReference : public LocalDynamics, public DataDelegateInner
{
  public:
    struct Edge { size_t j, reverse; Vecd gradient; };
    std::vector<std::vector<Edge>> rows;
    Real max_moment_condition=0, max_reproduction_error=0;
    Real elastic_dt_bound=0;
    explicit PDQuadraticReference(BaseInnerRelation &relation)
      : LocalDynamics(relation.getSPHBody()), DataDelegateInner(relation),
        volume(particles_->getVariableDataByName<Real>("VolumetricMeasure")),
        position(particles_->getVariableDataByName<Vecd>("Position")) {}

    void build()
    {
        static_assert(Dimensions==3,"Quadratic correspondence currently supports 3-D");
        using V9=Eigen::Matrix<Real,9,1>;using M9=Eigen::Matrix<Real,9,9>;
        size_t count=particles_->TotalRealParticles();rows.resize(count);
        for(size_t i=0;i<count;++i)
        {
            const auto &nb=inner_configuration_[i];
            Real sum=0;Matd covariance=Matd::Zero();
            for(size_t n=0;n<nb.current_size_;++n)
            {
                size_t j=nb.j_[n];Vecd xi=position[j]-position[i];
                Real a=nb.W_ij_[n]*volume[j];sum+=a;covariance+=a*xi*xi.transpose();
            }
            covariance/=sum;
            Eigen::SelfAdjointEigenSolver<Matd> eig(covariance);
            if(eig.info()!=Eigen::Success || eig.eigenvalues().minCoeff()<=0)
                throw std::runtime_error("Quadratic PD: invalid reference covariance");
            // Reference whitening resolves the physical thickness/tangent scale
            // disparity without changing the Euclidean horizon or quadrature.
            Matd L=eig.eigenvectors()*eig.eigenvalues().cwiseSqrt().cwiseInverse().asDiagonal()*eig.eigenvectors().transpose();
            auto basis=[&](const Vecd &xi)->V9 {
                Vecd z=L*xi;V9 b;b<<z[0],z[1],z[2],z[0]*z[0],z[1]*z[1],z[2]*z[2],z[0]*z[1],z[0]*z[2],z[1]*z[2];return b;
            };
            M9 moment=M9::Zero();
            for(size_t n=0;n<nb.current_size_;++n)
            {
                size_t j=nb.j_[n];V9 b=basis(position[j]-position[i]);
                moment+=(nb.W_ij_[n]*volume[j]/sum)*b*b.transpose();
            }
            Eigen::SelfAdjointEigenSolver<M9> meig(moment);
            if(meig.info()!=Eigen::Success || meig.eigenvalues().minCoeff()<=1e-12*meig.eigenvalues().maxCoeff())
                throw std::runtime_error("Quadratic PD: rank-deficient degree-two reference support");
            max_moment_condition=std::max(max_moment_condition,meig.eigenvalues().maxCoeff()/meig.eigenvalues().minCoeff());
            Eigen::LDLT<M9> solve(moment);
            Eigen::Matrix<Real,9,3> expected=Eigen::Matrix<Real,9,3>::Zero();expected.topRows<3>()=L;
            Eigen::Matrix<Real,9,3> reproduced=Eigen::Matrix<Real,9,3>::Zero();
            rows[i].reserve(nb.current_size_);
            for(size_t n=0;n<nb.current_size_;++n)
            {
                size_t j=nb.j_[n];V9 b=basis(position[j]-position[i]);V9 coeff=solve.solve(b);
                Vecd g=(nb.W_ij_[n]*volume[j]/sum)*L.transpose()*coeff.head<3>();
                if(!g.allFinite())throw std::runtime_error("Quadratic PD: non-finite reference gradient");
                rows[i].push_back({j,0,g});reproduced+=b*g.transpose();
            }
            max_reproduction_error=std::max(max_reproduction_error,(reproduced-expected).norm()/expected.norm());
        }
        if(max_reproduction_error>1e-9)throw std::runtime_error("Quadratic PD: polynomial reproduction failed");
        std::vector<Real> stiffness_row(count,0);
        auto &mat=DynamicCast<HughesWingetSolid>(this,sph_body_->getMatterMaterial());
        // Conservative block row-sum bound for the INITIAL elastic Hessian.
        // Other finite-strain, acoustic, stabilization and constraint bounds
        // remain enforced separately by the caller.
        Real rho=mat.ReferenceDensity(),C=3*mat.getBulkModulusforHGC(3);
        for(size_t i=0;i<count;++i)
        {
            Vecd self=Vecd::Zero();Real sum_norm=0;
            for(auto &e:rows[i])
            {
                auto &other=rows[e.j];auto found=std::find_if(other.begin(),other.end(),[&](const Edge &r){return r.j==i;});
                if(found==other.end())throw std::runtime_error("Quadratic PD: asymmetric reference neighborhood");
                e.reverse=size_t(found-other.begin());self-=e.gradient;sum_norm+=e.gradient.norm();
            }
            sum_norm+=self.norm();Real factor=volume[i]*C*sum_norm;
            stiffness_row[i]+=factor*self.norm();
            for(const auto &e:rows[i])stiffness_row[e.j]+=factor*e.gradient.norm();
        }
        Real bound=0;for(size_t i=0;i<count;++i)bound=std::max(bound,stiffness_row[i]/(rho*volume[i]));
        elastic_dt_bound=.2/std::sqrt(bound);
    }

    void checkIntact(size_t i) const
    {
        const auto &nb=inner_configuration_[i];
        if(nb.current_size_!=rows[i].size())throw std::runtime_error("Quadratic PD: reference neighbors changed");
        for(size_t n=0;n<nb.current_size_;++n)
            if(!nb.bondLive_[n] || nb.j_[n]!=rows[i][n].j)throw std::runtime_error("Quadratic PD supports intact reference bonds only");
    }
    Real *volume;Vecd *position;
};

class NosbPDQuadraticSecondStep : public NosbPDSecondStep
{
    PDQuadraticReference &reference_;
  public:
    NosbPDQuadraticSecondStep(BaseInnerRelation &relation,PDQuadraticReference &reference)
      :NosbPDSecondStep(relation),reference_(reference){}
    void interaction(size_t i,Real dt=0)
    {
        reference_.checkIntact(i);
        Matd f=Matd::Zero(),increment=Matd::Zero();
        for(const auto &e:reference_.rows[i])
        {
            f+=(pos_[e.j]-pos_[i])*e.gradient.transpose();
            increment+=dt*(vel_[e.j]-vel_[i])*e.gradient.transpose();
        }
        // Feed the existing, unchanged objective material update its expected
        // tensors; its update() then obtains exactly the corrected gradients.
        Matd K=shape_K_1_[i].inverse();
        N_[i]=f*K;N_deltaU_[i]=increment*K;N_half_[i]=(f-.5*increment)*K;
    }
};

class NosbPDQuadraticThirdStep : public LocalDynamics
{
    PDQuadraticReference &reference_;Real *volume_,*mass_;Vecd *acc_,*force_;Matd *P_;
  public:
    NosbPDQuadraticThirdStep(BaseInnerRelation &relation,PDQuadraticReference &reference)
      :LocalDynamics(relation.getSPHBody()),reference_(reference),
       volume_(particles_->getVariableDataByName<Real>("VolumetricMeasure")),
       mass_(particles_->getVariableDataByName<Real>("Mass")),
       acc_(particles_->getVariableDataByName<Vecd>("PDAcceleration")),
       force_(particles_->getVariableDataByName<Vecd>("PDInternalForce")),
       P_(particles_->getVariableDataByName<Matd>("PDPK1")){}
    void interaction(size_t i,Real dt=0)
    {
        Vecd force=Vecd::Zero();
        for(const auto &e:reference_.rows[i])
        {
            const Vecd &reverse=reference_.rows[e.j][e.reverse].gradient;
            force+=volume_[i]*P_[i]*e.gradient-volume_[e.j]*P_[e.j]*reverse;
        }
        force_[i]=force;acc_[i]+=force/mass_[i];
    }
};
}} // namespace SPH::solid_dynamics
