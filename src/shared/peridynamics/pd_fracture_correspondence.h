#pragma once
/** Live-bond polynomial correspondence for strength-controlled fracture.
 * The reference coordinates remain fixed. Only the live support changes.
 * All force rows retain the original reciprocal indexing, with zero gradient
 * on a broken bond. The existing EXACT transpose force operator is reused.
 * A rank-complete affine fallback is explicit. Rank-deficient live supports
 * stop the calculation; numerical rank loss must not cause material failure.
 * This is not a calibrated fracture-energy regularization.
 */
#include "pd_quadratic_correspondence.h"

namespace SPH::solid_dynamics
{
class PDFractureReference : public PDQuadraticReference
{
  public:
    StdVec<Vecd> X;
    StdVec<int> degree;
    StdVec<size_t> changed_particles;
    size_t topology_version=0,affine_rows=0,ballistic_particles=0,resolution_removed_bonds=0;
    Real maximum_affine_error=0;
    PDFractureReference(BaseInnerRelation &relation,bool erode=false)
      :PDQuadraticReference(relation),
       particle_live(particles_->getVariableDataByName<int>("ParticleLive")),
       damage(particles_->getVariableDataByName<Real>("Damage")),
       bonds(particles_->getVariableDataByName<int>("BondNumber")),
       stress(particles_->getVariableDataByName<Matd>("PDStress")),
       reproduction(particles_->registerStateVariableData<int>("PDReproductionOrder")),
       eroded(particles_->registerStateVariableData<int>("ResolutionEroded"))
    {if(erode)throw std::invalid_argument("Numerical resolution erosion is excluded from physical fracture calibration");}
    void initialize()
    {
        X.assign(position,position+particles_->TotalRealParticles());PDQuadraticReference::build();
        degree.assign(rows.size(),2);mask.resize(rows.size());
        for(size_t i=0;i<rows.size();++i){mask[i].assign(rows[i].size(),true);reproduction[i]=2;eroded[i]=0;}
    }
    const Neighborhood &neighbors(size_t i)const{return inner_configuration_[i];}
    bool live(size_t i,size_t n)const{return inner_configuration_[i].bondLive_[n];}
    void updateAfterFailure()
    {
        using V9=Eigen::Matrix<Real,9,1>;using M9=Eigen::Matrix<Real,9,9>;
        if(X.size()!=rows.size())throw std::logic_error("Fracture reference not initialized");
        bool repeat=true;size_t passes=0;
        while(repeat)
        {
            if(++passes>rows.size()+1)throw std::runtime_error("Fragment resolution update did not terminate");
            repeat=false;affine_rows=ballistic_particles=0;maximum_affine_error=0;
            StdVec<size_t> unsupported;
            for(size_t i=0;i<rows.size();++i)
            {
                const auto &nb=inner_configuration_[i];size_t live_count=0;Real sum=0;Matd covariance=Matd::Zero();
                for(size_t n=0;n<rows[i].size();++n)
                {
                    rows[i][n].gradient.setZero();if(!nb.bondLive_[n])continue;
                    ++live_count;size_t j=rows[i][n].j;Vecd xi=X[j]-X[i];Real a=nb.W_ij_[n]*volume[j];sum+=a;covariance+=a*xi*xi.transpose();
                }
                if(live_count==0)
                {degree[i]=reproduction[i]=0;particle_live[i]=0;stress[i].setZero();++ballistic_particles;continue;}
                covariance/=sum;Eigen::SelfAdjointEigenSolver<Matd> ce(covariance);
                if(ce.info()!=Eigen::Success||ce.eigenvalues().minCoeff()<=1e-6*ce.eigenvalues().maxCoeff())
                {unsupported.push_back(i);continue;}
                Matd L=ce.eigenvectors()*ce.eigenvalues().cwiseSqrt().cwiseInverse().asDiagonal()*ce.eigenvectors().transpose();
                auto basis=[&](size_t j)->V9 {Vecd z=L*(X[j]-X[i]);V9 b;b<<z[0],z[1],z[2],z[0]*z[0],z[1]*z[1],z[2]*z[2],z[0]*z[1],z[0]*z[2],z[1]*z[2];return b;};
                M9 moment=M9::Zero();
                for(size_t n=0;n<rows[i].size();++n)if(nb.bondLive_[n])
                {V9 b=basis(rows[i][n].j);moment+=(nb.W_ij_[n]*volume[rows[i][n].j]/sum)*b*b.transpose();}
                Eigen::SelfAdjointEigenSolver<M9> me(moment);
                // Keep the intact bending fit away from cracks. A damaged row
                // uses one-sided affine correspondence: imposing smooth
                // quadratic continuation next to a displacement discontinuity
                // creates amplified stencils and needlessly small time steps.
                bool quadratic=live_count==rows[i].size()&&live_count>=9&&me.info()==Eigen::Success&&me.eigenvalues().minCoeff()>1e-8*me.eigenvalues().maxCoeff();
                Eigen::LDLT<M9> solve(moment);Matd inv_cov=covariance.inverse(),reproduced=Matd::Zero();
                // A nominally invertible moment can still lose numerical
                // reproduction. Test the fit itself before accepting its order.
                if(quadratic)
                {
                    Eigen::Matrix<Real,9,3> check=Eigen::Matrix<Real,9,3>::Zero();
                    for(size_t n=0;n<rows[i].size();++n)if(nb.bondLive_[n])
                    {auto &e=rows[i][n];Real a=nb.W_ij_[n]*volume[e.j]/sum;V9 b=basis(e.j),c=solve.solve(b);
                     e.gradient=a*L.transpose()*c.head<3>();check+=a*b*c.head<3>().transpose();}
                    Eigen::Matrix<Real,9,3> exact=Eigen::Matrix<Real,9,3>::Zero();exact.topRows<3>().setIdentity();
                    quadratic=check.allFinite()&&(check-exact).norm()<1e-8;
                }
                if(!quadratic)for(size_t n=0;n<rows[i].size();++n)if(nb.bondLive_[n])
                {auto &e=rows[i][n];e.gradient=(nb.W_ij_[n]*volume[e.j]/sum)*inv_cov*(X[e.j]-X[i]);}
                for(auto &e:rows[i])reproduced+=(X[e.j]-X[i])*e.gradient.transpose();
                Real error=(reproduced-Matd::Identity()).norm();
                if(!reproduced.allFinite()||error>1e-7){unsupported.push_back(i);continue;}
                degree[i]=reproduction[i]=quadratic?2:1;affine_rows+=!quadratic;
                maximum_affine_error=std::max(maximum_affine_error,error);
            }
            if(!unsupported.empty())
            {
                throw std::runtime_error("Unresolved live PD support: refine or use a separately verified discretization; no bonds or interface constraints were removed to bypass rank loss");
            }
        }
        changed_particles.clear();
        for(size_t i=0;i<rows.size();++i)
        {
            auto &nb=inner_configuration_[i];size_t count=0;bool changed=false;
            for(size_t n=0;n<rows[i].size();++n){changed|=bool(mask[i][n])!=bool(nb.bondLive_[n]);mask[i][n]=nb.bondLive_[n];count+=nb.bondLive_[n];
                if(nb.bondLive_[n]!=inner_configuration_[rows[i][n].j].bondLive_[rows[i][n].reverse])throw std::runtime_error("Nonreciprocal fracture topology");}
            if(changed)changed_particles.push_back(i);
            bonds[i]=int(count);damage[i]=rows[i].empty()?0:1-Real(count)/rows[i].size();particle_live[i]=count?1:0;
        }
        updateBound();++topology_version;
    }
    size_t topologyVersion()const override{return topology_version;}
    void checkIntact(size_t i)const override
    {
        const auto &nb=inner_configuration_[i];if(nb.current_size_!=rows[i].size())throw std::runtime_error("Fracture reference neighbors changed");
        for(size_t n=0;n<rows[i].size();++n)if(nb.j_[n]!=rows[i][n].j||bool(nb.bondLive_[n])!=bool(mask[i][n]))
            throw std::runtime_error("Fracture topology changed without rebuilding the gradient");
    }
  private:
    StdVec<StdVec<bool>> mask;int *particle_live,*bonds;Real *damage;Matd *stress;int *reproduction,*eroded;
    void updateBound()
    {
        StdVec<Real> stiffness(rows.size(),0);auto &m=DynamicCast<HughesWingetSolid>(this,sph_body_->getMatterMaterial());
        Real C=3*m.getBulkModulusforHGC(3),rho=m.ReferenceDensity();
        for(size_t i=0;i<rows.size();++i)
        {Vecd self=Vecd::Zero();Real norm=0;for(auto &e:rows[i]){self-=e.gradient;norm+=e.gradient.norm();}
         norm+=self.norm();Real a=volume[i]*C*norm;stiffness[i]+=a*self.norm();for(auto &e:rows[i])stiffness[e.j]+=a*e.gradient.norm();}
        Real largest=0;for(size_t i=0;i<rows.size();++i)largest=std::max(largest,stiffness[i]/(rho*volume[i]));
        elastic_dt_bound=largest>0?.2/std::sqrt(largest):std::numeric_limits<Real>::infinity();
    }
};
}
