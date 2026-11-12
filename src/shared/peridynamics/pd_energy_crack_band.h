#pragma once
/** Energy-regularized continuum damage for quadratic NOSB correspondence.
 * This is a diffuse crack-band representation, NOT an extrinsic cohesive facet
 * or an energy-based BB-PD critical-stretch formula. No bonds/particles are
 * numerically deleted. Volume compression is retained after tensile failure.
 * The support-defined band width is 2*horizon. Its work-of-separation accuracy
 * must be verified on the actual quadrature before quantitative calibration.
 */
#include "pd_quadratic_correspondence.h"
#include <Eigen/SVD>

namespace SPH::solid_dynamics
{
struct CorotatedBiotEnergy
{
    Real E,nu,K,mu;
    struct Result {Matd positive_P,negative_P;Real positive,negative,max_effective_tension,J;Vecd normal;};
    CorotatedBiotEnergy(Real young,Real poisson):E(young),nu(poisson),K(E/(3*(1-2*nu))),mu(E/(2*(1+nu)))
    {if(E<=0||nu<=-1||nu>=.5)throw std::invalid_argument("Invalid Biot elastic constants");}
    Result evaluate(const Matd&F)const
    {
        Real J=F.determinant();if(!F.allFinite()||!std::isfinite(J)||J<=0)throw std::runtime_error("Invalid energy-material deformation gradient");
        Eigen::JacobiSVD<Matd> svd(F,Eigen::ComputeFullU|Eigen::ComputeFullV);
        if(svd.info()!=Eigen::Success||svd.singularValues().minCoeff()<=0)throw std::runtime_error("Biot polar decomposition failed");
        Matd R=svd.matrixU()*svd.matrixV().transpose();Vecd e=svd.singularValues()-Vecd::Ones();Real tr=e.sum();
        Matd strain=svd.matrixV()*e.asDiagonal()*svd.matrixV().transpose(),dev=strain-tr/3*Matd::Identity();
        Result a;a.positive=.5*K*std::pow(std::max(tr,Real(0)),2)+mu*dev.squaredNorm();a.negative=.5*K*std::pow(std::min(tr,Real(0)),2);
        a.positive_P=R*(K*std::max(tr,Real(0))*Matd::Identity()+2*mu*dev);a.negative_P=R*(K*std::min(tr,Real(0))*Matd::Identity());
        Eigen::Index direction;e.maxCoeff(&direction);a.normal=svd.matrixV().col(direction);
        a.max_effective_tension=K*tr+2*mu*(e[direction]-tr/3);a.J=J;return a;
    }
};


/** Finite-volume correction of the Biot distortional solid. The compressive
 * part MUST depend on J, not tr(U-I): the latter decreases under isochoric
 * shear in a compressed, fully damaged material and supplies negative shear
 * stiffness. The volumetric energy below is the standard convex-in-J form
 * K/4*(J^2-1-2*log(J)), with its logarithmic non-interpenetration barrier.
 * It has the same infinitesimal K and changes no elastic input constants.
 * Distortional Biot energy remains objective; its exact derivative is coaxial
 * with U. No strain/force clipping or residual tensile stiffness is used.
 */
struct FiniteVolumeBiotEnergy
{
    using Result=CorotatedBiotEnergy::Result;
    Real E,nu,K,mu;
    FiniteVolumeBiotEnergy(Real young,Real poisson):E(young),nu(poisson),K(E/(3*(1-2*nu))),mu(E/(2*(1+nu)))
    {if(E<=0||nu<=-1||nu>=.5)throw std::invalid_argument("Invalid finite-volume Biot elastic constants");}
    Result evaluate(const Matd&F)const
    {
        Real J=F.determinant();if(!F.allFinite()||!std::isfinite(J)||J<=0)throw std::runtime_error("Invalid finite-volume deformation gradient");
        Eigen::JacobiSVD<Matd> svd(F,Eigen::ComputeFullU|Eigen::ComputeFullV);
        Vecd l=svd.singularValues();if(svd.info()!=Eigen::Success||l.minCoeff()<=0)throw std::runtime_error("Finite-volume polar decomposition failed");
        Matd R=svd.matrixU()*svd.matrixV().transpose();Vecd e=l-Vecd::Constant(l.mean());
        Matd dev=svd.matrixV()*e.asDiagonal()*svd.matrixV().transpose();
        Real z=std::log(J),j2m1=std::expm1(2*z),vol=.25*K*(j2m1-2*z);
        Matd pv=.5*K*j2m1*F.inverse().transpose();
        Result a;a.positive=mu*dev.squaredNorm()+(J>=1?vol:0);a.negative=J<1?vol:0;
        a.positive_P=R*(2*mu*dev);a.negative_P=Matd::Zero();
        if(J>=1)a.positive_P+=pv;else a.negative_P=pv;
        Vecd principal_cauchy=(.5*K*j2m1*Vecd::Ones()+2*mu*e.cwiseProduct(l))/J;
        Eigen::Index direction;a.max_effective_tension=principal_cauchy.maxCoeff(&direction);a.normal=svd.matrixV().col(direction);a.J=J;return a;
    }
    Real tangentBound(const Matd&F,Real damage)const
    {
        Eigen::JacobiSVD<Matd> svd(F);Vecd l=svd.singularValues();Real a=l.minCoeff(),J=F.determinant(),g=1-damage;
        if(a<=0||J<=0)throw std::runtime_error("Invalid tangent-bound deformation");
        Real vol=K*J*J*l.cwiseInverse().squaredNorm()+.5*K*std::abs(J*J-1)/(a*a);
        Real dev=2*mu*(l.maxCoeff()+(l.array()-l.mean()).abs().maxCoeff())/a;
        return (J<1?1:g)*vol+g*dev;
    }
    /** Stress-free lateral contraction of an undamaged uniaxial specimen.
     * Solved from this energy, not a small-strain Poisson extrapolation. */
    Matd uniaxial(Real axial)const
    {
        if(axial<1)throw std::invalid_argument("Tensile calibration requires axial stretch >= 1");
        auto lateral_stress=[&](Real b){Real J=axial*b*b;return .5*K*(J*J-1)/b-2*mu*(axial-b)/3;};
        Real lo=.05,hi=std::max(Real(2),axial);
        if(lateral_stress(lo)>=0||lateral_stress(hi)<=0)throw std::runtime_error("Uniaxial lateral root is not bracketed");
        for(int k=0;k<80;++k){Real b=.5*(lo+hi);if(lateral_stress(b)>0)hi=b;else lo=b;}
        Matd F=Matd::Identity();F(0,0)=axial;F(1,1)=F(2,2)=.5*(lo+hi);return F;
    }
    Real initiationEnergy(Real nominal_strength)const
    {
        if(nominal_strength<=0)throw std::invalid_argument("Positive nominal tensile strength required");
        Real lo=1,hi=1+std::max(Real(.01),2*nominal_strength/E);
        auto traction=[&](Real a){auto r=evaluate(uniaxial(a));return r.positive_P(0,0)+r.negative_P(0,0);};
        for(int k=0;k<40&&traction(hi)<nominal_strength;++k)hi=1+2*(hi-1);
        if(traction(hi)<nominal_strength)throw std::runtime_error("Nominal tensile strength could not be bracketed");
        for(int k=0;k<80;++k){Real m=.5*(lo+hi);if(traction(m)>nominal_strength)hi=m;else lo=m;}
        auto r=evaluate(uniaxial(.5*(lo+hi)));if(r.negative>1e-12*r.positive)throw std::runtime_error("Uniaxial calibration contains a retained compressive energy");return r.positive;
    }
};

struct LinearEnergyCrackBand
{
    Real E,ft,Gc,width,k0,kf,dissipation_slope;
    struct History {Real kappa=0,damage=0,dissipation=0;};
    LinearEnergyCrackBand(Real young,Real strength,Real energy,Real band_width,Real initiation_energy=-1)
      :E(young),ft(strength),Gc(energy),width(band_width)
    {
        if(!std::isfinite(E*ft*Gc*width)||std::min({E,ft,Gc,width})<=0)throw std::invalid_argument("Positive finite E, tensile strength, Gc and band width required");
        if(initiation_energy!=-1 && (!std::isfinite(initiation_energy)||initiation_energy<=0))throw std::invalid_argument("Invalid initiation energy");
        k0=initiation_energy>0?std::sqrt(2*initiation_energy/E):ft/E;kf=2*Gc/(E*k0*width);
        if(kf<=k0*(1+1e-10))throw std::invalid_argument("Fracture band cannot resolve this (ft,Gc): terminal energy-equivalent strain must exceed initiation strain; refine the discretization, do not clip material parameters");
        dissipation_slope=.5*E*k0*kf/(kf-k0);
    }
    void advance(History &s,Real positive_energy,bool tensile)const
    {
        if(!std::isfinite(positive_energy)||positive_energy<0)throw std::runtime_error("Invalid crack-driving energy");
        if(!tensile)return;
        Real k=std::max(s.kappa,std::sqrt(2*positive_energy/E));
        Real old=std::clamp(s.kappa,k0,kf),now=std::clamp(k,k0,kf);
        s.dissipation+=dissipation_slope*(now-old);s.kappa=k;
        s.damage=k<=k0?0:k>=kf?1:kf*(k-k0)/(k*(kf-k0));
        if(s.damage<0||s.damage>1||s.dissipation<0||!std::isfinite(s.dissipation))throw std::runtime_error("Invalid irreversible damage state");
    }
};

/** Full degree-two least-squares residual. The residual is orthogonal to the
 * polynomial space, so its exact derivative needs no additional fit derivative.
 * This supplies positive energy to zero-gradient oscillations and vanishes for
 * rigid, affine AND quadratic displacement fields. Its energy participates in
 * the same damage driving force/budget; it is never an undegraded crack bridge.
 */
class PDEnergyCrackBand : public LocalDynamics,public DataDelegateInner
{
  public:
    using V9=Eigen::Matrix<Real,9,1>;using M9=Eigen::Matrix<Real,9,9>;using A39=Eigen::Matrix<Real,3,9>;
    struct Entry {V9 basis,coefficient;Real weight;Vecd residual,stabilizing_force;};
    StdVec<StdVec<Entry>> projection;
    StdVec<LinearEnergyCrackBand::History> history;
    PDQuadraticReference &reference;FiniteVolumeBiotEnergy elastic;LinearEnergyCrackBand fracture;
    Real stabilizer,horizon,stabilizer_dt=0;
    Real *damage,*dissipation_density,*stored_density,*hg_density,*kappa_field,*tangent_bound;int *material_failed;
    Vecd *crack_normal;
    PDEnergyCrackBand(BaseInnerRelation&r,PDQuadraticReference&ref,Real E,Real nu,Real ft,Real Gc,Real hg,Real delta)
      :LocalDynamics(r.getSPHBody()),DataDelegateInner(r),reference(ref),elastic(E,nu),fracture(E,ft,Gc,2*delta,elastic.initiationEnergy(ft)),horizon(delta),
       x(particles_->getVariableDataByName<Vecd>("Position")),volume(particles_->getVariableDataByName<Real>("VolumetricMeasure")),
       mass(particles_->getVariableDataByName<Real>("Mass")),F(particles_->getVariableDataByName<Matd>("DeformationGradient")),
       P(particles_->getVariableDataByName<Matd>("PDPK1")),stress(particles_->getVariableDataByName<Matd>("PDStress")),
       acc(particles_->getVariableDataByName<Vecd>("PDAcceleration")),internal(particles_->getVariableDataByName<Vecd>("PDInternalForce")),
       hourglass(particles_->getVariableDataByName<Vecd>("HourglassForce")),damping(particles_->getVariableDataByName<Vecd>("PDNumericalDampingForce"))
    {
        if(hg<0||delta<=0)throw std::invalid_argument("Invalid energy stabilization");stabilizer=hg*elastic.K/(delta*delta);
        damage=particles_->getVariableDataByName<Real>("Damage");
        tangent_bound=particles_->registerStateVariableData<Real>("ShellTangentModulusBound");
        dissipation_density=particles_->registerStateVariableData<Real>("FractureDissipationDensity");stored_density=particles_->registerStateVariableData<Real>("ShellStoredEnergyDensity");
        hg_density=particles_->registerStateVariableData<Real>("ShellStabilizationEnergyDensity");kappa_field=particles_->registerStateVariableData<Real>("FractureHistoryKappa");
        material_failed=particles_->registerStateVariableData<int>("MaterialTensileFailure");crack_normal=particles_->registerStateVariableData<Vecd>("CrackNormalReference");
    }
    void initialize()
    {
        size_t count=particles_->TotalRealParticles();history.resize(count);projection.resize(count);Real max_rate=0;
        for(size_t i=0;i<count;++i)
        {
            const auto&nb=inner_configuration_[i];Real sum=0;Matd cov=Matd::Zero();
            for(size_t n=0;n<nb.current_size_;++n){Real a=nb.W_ij_[n]*volume[nb.j_[n]];Vecd z=x[nb.j_[n]]-x[i];sum+=a;cov+=a*z*z.transpose();}
            cov/=sum;Eigen::SelfAdjointEigenSolver<Matd> ce(cov);if(ce.info()!=Eigen::Success||ce.eigenvalues().minCoeff()<=0)throw std::runtime_error("Invalid energetic projection covariance");
            Matd L=ce.eigenvectors()*ce.eigenvalues().cwiseSqrt().cwiseInverse().asDiagonal()*ce.eigenvectors().transpose();
            M9 moment=M9::Zero();projection[i].resize(nb.current_size_);
            for(size_t n=0;n<nb.current_size_;++n){auto&e=projection[i][n];Vecd z=L*(x[nb.j_[n]]-x[i]);e.basis<<z[0],z[1],z[2],z[0]*z[0],z[1]*z[1],z[2]*z[2],z[0]*z[1],z[0]*z[2],z[1]*z[2];
                e.weight=nb.W_ij_[n]*volume[nb.j_[n]]/sum;moment+=e.weight*e.basis*e.basis.transpose();}
            Eigen::LDLT<M9> solve(moment);Eigen::Matrix<Real,9,9> check=M9::Zero();
            for(auto&e:projection[i]){e.coefficient=e.weight*solve.solve(e.basis);check+=e.basis*e.coefficient.transpose();}
            if((check-M9::Identity()).norm()>1e-8)throw std::runtime_error("Energetic quadratic projection failed");
        }
        // The least-squares residual Hessian is bounded above by the unprojected
        // weighted graph Laplacian, giving a conservative positive spectral bound.
        StdVec<Real> rate(count,0);
        for(size_t i=0;i<count;++i)for(size_t n=0;n<reference.rows[i].size();++n)
        {auto&e=reference.rows[i][n];Real k=volume[i]*stabilizer*projection[i][n].weight;rate[i]+=2*k/mass[i];rate[e.j]+=2*k/mass[e.j];}
        for(Real a:rate)max_rate=std::max(max_rate,a);stabilizer_dt=max_rate>0?.2/std::sqrt(max_rate):std::numeric_limits<Real>::infinity();
    }
    void evaluate(size_t i,bool advance=true)
    {
        reference.checkIntact(i);A39 A=A39::Zero();Matd f=Matd::Zero();
        for(size_t n=0;n<reference.rows[i].size();++n){auto&e=reference.rows[i][n];Vecd y=x[e.j]-x[i];A+=y*projection[i][n].coefficient.transpose();f+=y*e.gradient.transpose();}
        F[i]=f;auto trial=elastic.evaluate(f);Real hg=0;
        for(size_t n=0;n<reference.rows[i].size();++n){auto&e=projection[i][n];e.residual=x[reference.rows[i][n].j]-x[i]-A*e.basis;hg+=.5*stabilizer*e.weight*e.residual.squaredNorm();}
        Real old=history[i].damage;
        if(advance)fracture.advance(history[i],trial.positive+hg,trial.max_effective_tension>1e-12*elastic.E);
        if(old==0&&history[i].damage>0)crack_normal[i]=trial.normal;
        tangent_bound[i]=elastic.tangentBound(f,history[i].damage);
        Real g=1-history[i].damage;P[i]=g*trial.positive_P+trial.negative_P;stress[i]=P[i]*f.transpose()/trial.J;
        damage[i]=history[i].damage;dissipation_density[i]=history[i].dissipation;kappa_field[i]=history[i].kappa;material_failed[i]=damage[i]==1;
        stored_density[i]=g*(trial.positive+hg)+trial.negative;hg_density[i]=g*hg;
        for(auto&e:projection[i])e.stabilizing_force=volume[i]*g*stabilizer*e.weight*e.residual;
    }
    void forces(size_t i)
    {
        Vecd f=Vecd::Zero(),h=Vecd::Zero();for(size_t n=0;n<reference.rows[i].size();++n)
        {auto&e=reference.rows[i][n];f+=volume[i]*P[i]*e.gradient-volume[e.j]*P[e.j]*reference.rows[e.j][e.reverse].gradient;
         h+=projection[i][n].stabilizing_force-projection[e.j][e.reverse].stabilizing_force;}
        internal[i]=f;hourglass[i]=h;damping[i].setZero();acc[i]+=(f+h)/mass[i];
    }
    Real finiteStrainTimeStep(Real initial_elastic_dt)const
    {
        Real bound=3*elastic.K;for(size_t i=0;i<history.size();++i)bound=std::max(bound,tangent_bound[i]);
        return initial_elastic_dt*std::sqrt(3*elastic.K/bound);
    }
    struct Statistics {Real stored=0,stabilization=0,dissipated=0,max_damage=0;size_t initiated=0,failed=0;};
    Statistics statistics()const
    {Statistics s;for(size_t i=0;i<history.size();++i){s.stored+=volume[i]*stored_density[i];s.stabilization+=volume[i]*hg_density[i];s.dissipated+=volume[i]*dissipation_density[i];
     s.max_damage=std::max(s.max_damage,damage[i]);s.initiated+=damage[i]>0;s.failed+=damage[i]==1;}return s;}
  private:
    Vecd *x,*acc,*internal,*hourglass,*damping;Real *volume,*mass;Matd *F,*P,*stress;
};
}
