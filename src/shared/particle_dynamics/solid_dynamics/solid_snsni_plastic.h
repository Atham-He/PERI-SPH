/* Optional finite-strain S-NSNI-based RKPM operator.
 * Wang et al. 10.1007/s00466-024-02494-0, Eqs.57,59-60,72-76.
 * Adaptations: native multiplicative Hencky DPC, independent centroid material
 * sites, complete pushed-forward cell moments. No added hourglass penalty,
 * force scaling, yield clipping or artificial stress is used.
 */
#pragma once
#include "solid_snsni_reference.h"
#include "dpc_snsni_adapter.h"
namespace SPH::solid_dynamics
{
class SnsniPlasticSolid
{
 public:
    using Field=StdVec<Vecd>;
    using Triple=std::array<Matd,3>;
    SnsniReference&reference;
    dpc::Model model;
    std::vector<dpc::State> history;
    std::vector<dpc::SnsniMaterialStep> pending;
    std::vector<std::vector<SnsniRKCloud::Entry>> sampling_gradient;
    std::vector<Matd> F,inverseF,correction,force_gradient;
    std::vector<Triple> displacement_hessian,force_hessian;
    StdVec<Real> accumulated_plastic_work;
    Field internal;
    size_t evaluations=0,commits=0;
    bool pending_commit=false;
    Real linear_frequency_bound_squared=0,acoustic_dt=0;
    SnsniPlasticSolid(SnsniReference&r,const dpc::Model&m):reference(r),model(m)
    {
        const size_t n=r.cells.size();history.resize(n);pending.resize(n);F.resize(n);inverseF.resize(n);correction.resize(n);force_gradient.resize(n);displacement_hessian.resize(n);force_hessian.resize(n);accumulated_plastic_work.assign(n,0);sampling_gradient.resize(n);
        StdVec<Vecd> sites;StdVec<Real> volumes;
        for(auto&c:r.cells){sites.push_back(c.center);volumes.push_back(c.volume/std::pow(r.dp,3));}
        SnsniRKCloud samples(sites,r.dp,volumes);
        tbb::parallel_for(size_t(0),n,[&](size_t l){sampling_gradient[l]=samples.shape(sites[l]);});
        // Block row-sum bound of M^-1/2 K_elastic M^-1/2. The Frobenius
        // strain bound uses the largest Kelvin elastic eigenvalue.
        auto a=model.parameters(model.r0());Real C=std::max(3*a.K,2*a.G);
        StdVec<Real> rate(r.X.size(),0);
        for(const auto&cell:r.cells){Eigen::SelfAdjointEigenSolver<Matd> eig(cell.moment);
            Real sum=0;Vecd higher=Vecd::Zero();
            for(const auto&e:cell.rows){Real im=1/std::sqrt(r.mass[e.j]);sum+=im*e.gradient.norm();for(int k=0;k<3;++k)higher[k]+=im*(e.hessian*eig.eigenvectors().col(k)).norm();}
            for(const auto&e:cell.rows){Real im=1/std::sqrt(r.mass[e.j]);rate[e.j]+=C*cell.volume*im*e.gradient.norm()*sum;
                for(int k=0;k<3;++k)rate[e.j]+=C*eig.eigenvalues()[k]*im*(e.hessian*eig.eigenvectors().col(k)).norm()*higher[k];}}
        linear_frequency_bound_squared=*std::max_element(rate.begin(),rate.end());
        if(!(linear_frequency_bound_squared>0))throw std::runtime_error("Invalid S-NSNI frequency bound");
    }
    void evaluate(const Field&d,Real rho0,Real cfl)
    {
        if(pending_commit)throw std::logic_error("Commit or explicitly discard pending DPC history before the next evaluation");
        const size_t n=history.size();
        tbb::parallel_for(size_t(0),n,[&](size_t l){auto&h=displacement_hessian[l];for(auto&a:h)a.setZero();F[l]=Matd::Identity();
            for(const auto&e:reference.cells[l].rows){F[l]+=d[e.j]*e.gradient.transpose();for(int a=0;a<3;++a)h[a]+=d[e.j][a]*e.hessian;}
            if(!F[l].allFinite()||F[l].determinant()<=0)throw std::runtime_error("Invalid S-NSNI cell Jacobian, no deletion or clamping");inverseF[l]=F[l].inverse();
            Triple derivativeF,derivativeCp;for(auto&a:derivativeCp)a.setZero();Vecd derivativeKappa=Vecd::Zero();
            // Difference form enforces the derivative of a constant exactly;
            // algebraically it is the partition-of-unity RK derivative.
            for(const auto&e:sampling_gradient[l]){Matd delta=history[e.j].inverse_Cp-history[l].inverse_Cp;Real dk=history[e.j].kappa-history[l].kappa;
                for(int k=0;k<3;++k)derivativeCp[k]+=e.gradient[k]*delta;derivativeKappa+=dk*e.gradient;}
            Triple cp_spatial;Vecd kappa_spatial=inverseF[l].transpose()*derivativeKappa;
            for(int k=0;k<3;++k){derivativeF[k].setZero();cp_spatial[k].setZero();
                for(int j=0;j<3;++j){cp_spatial[k]+=inverseF[l](j,k)*derivativeCp[j];for(int a=0;a<3;++a)derivativeF[k].row(a)+=inverseF[l](j,k)*h[a].col(j).transpose();}}
            pending[l]=dpc::snsniMaterialStep(model,F[l],history[l],derivativeF,cp_spatial,kappa_spatial);
            correction[l]=pending[l].correction;
        });
        tbb::parallel_for(size_t(0),n,[&](size_t l){const auto&cell=reference.cells[l];const auto&inv=inverseF[l];
            Triple derivative;for(auto&a:derivative)a.setZero();
            for(const auto&e:sampling_gradient[l]){Matd delta=correction[e.j]-correction[l];for(int k=0;k<3;++k)derivative[k]+=e.gradient[k]*delta;}
            Triple sigma_gradient;
            for(int k=0;k<3;++k){sigma_gradient[k]=pending[l].trial.gradient[k];for(int j=0;j<3;++j)sigma_gradient[k]-=inverseF[l](j,k)*derivative[j];}
            const Real J=F[l].determinant();Matd moment=J*F[l]*cell.moment*F[l].transpose(),geometric=Matd::Zero();
            for(int a=0;a<3;++a){auto&T=force_hessian[l][a];T.setZero();
                for(int k=0;k<3;++k)T+=(inv*sigma_gradient[k].row(a).transpose())*(inv*moment.row(k).transpose()).transpose();
                for(int b=0;b<3;++b)geometric(a,b)=(T.array()*displacement_hessian[l][b].array()).sum();}
            // Eq.60 is fully retained: H_x=F^-T(H_X-grad_x(N).hess(u))F^-1.
            // Contract first, then gather, to avoid recomputing it per DOF.
            force_gradient[l]=(J*cell.volume*pending[l].material.sigma-geometric)*inv.transpose();
        });
        internal.resize(reference.X.size());
        tbb::parallel_for(size_t(0),internal.size(),[&](size_t j){Vecd f=Vecd::Zero();for(auto[l,k]:reference.transpose[j]){const auto&e=reference.cells[l].rows[k];f+=force_gradient[l]*e.gradient;
            for(int a=0;a<3;++a)f[a]+=(force_hessian[l][a].array()*e.hessian.array()).sum();}internal[j]=f;});
        acoustic_dt=std::numeric_limits<Real>::infinity();
        for(size_t l=0;l<n;++l){auto a=model.parameters(pending[l].material.relative_density,pending[l].material.state.alpha);
            Eigen::SelfAdjointEigenSolver<Matd> eig(inverseF[l].transpose()*inverseF[l]);
            Real inverse_stretch=std::sqrt(eig.eigenvalues().maxCoeff());
            Real sound=std::sqrt((a.K+4*a.G/3)*F[l].determinant()/rho0);
            acoustic_dt=std::min(acoustic_dt,cfl*reference.dp/(inverse_stretch*sound));}
        ++evaluations;pending_commit=true;
    }
    void commit()
    {
        if(!pending_commit)throw std::logic_error("S-NSNI material history committed twice");
        for(size_t l=0;l<history.size();++l){history[l]=pending[l].material.state;accumulated_plastic_work[l]+=pending[l].material.plastic_stress_work;}
        ++commits;pending_commit=false;
    }
    void discard(){if(!pending_commit)throw std::logic_error("No pending S-NSNI state");pending_commit=false;}
};

/** Explicit Newmark beta=0,gamma=1/2 (kick-drift-kick), constant lumped M.
 * Optional symmetric exact damping splits retain the pre-existing numerical
 * relaxation rate; the undamped scheme is the central-difference method.
 */
class SnsniCentralDifference
{
 public:
    using Field=StdVec<Vecd>;Field displacement,velocity,acceleration;
    const StdVec<Real>&mass;Real damping_loss=0;size_t steps=0;
    explicit SnsniCentralDifference(const StdVec<Real>&m):displacement(m.size(),Vecd::Zero()),velocity(displacement),acceleration(displacement),mass(m){}
    Real kinetic()const{Real v=0;for(size_t j=0;j<mass.size();++j)v+=.5*mass[j]*velocity[j].squaredNorm();return v;}
    template<class Evaluate>void step(Real dt,Real damping,Evaluate&&evaluate)
    {
        Real factor=std::exp(-.5*damping*dt);
        tbb::parallel_for(size_t(0),mass.size(),[&](size_t j){velocity[j]*=factor;velocity[j]+=.5*dt*acceleration[j];displacement[j]+=dt*velocity[j];});
        Field force;evaluate(displacement,velocity,force);
        tbb::parallel_for(size_t(0),mass.size(),[&](size_t j){acceleration[j]=force[j]/mass[j];velocity[j]+=.5*dt*acceleration[j];velocity[j]*=factor;});++steps;
    }
};
}
