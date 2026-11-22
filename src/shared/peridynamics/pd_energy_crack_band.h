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
#include "pd_directional_crack_energy.h"
#include "pd_layered_geometry.h"
#include "pd_finite_volume_biot_terms.h"

namespace SPH::solid_dynamics
{
struct CorotatedBiotEnergy
{
    Real E,nu,K,mu;
    struct Result {Matd positive_P,negative_P;Real positive,negative,max_effective_tension,J;Vecd normal,principal_stretches;};
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
        a.max_effective_tension=K*tr+2*mu*(e[direction]-tr/3);a.J=J;a.principal_stretches=svd.singularValues();return a;
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
        const auto volumetric=peridynamics::finiteVolumeBiotVolumetric(K,J);
        Matd pv=volumetric.piola_scale*F.inverse().transpose();
        Result a;a.positive=mu*dev.squaredNorm()+(J>=1?volumetric.energy:0);a.negative=J<1?volumetric.energy:0;
        a.positive_P=R*(2*mu*dev);a.negative_P=Matd::Zero();
        if(J>=1)a.positive_P+=pv;else a.negative_P=pv;
        Vecd principal_cauchy=(volumetric.piola_scale*Vecd::Ones()+2*mu*e.cwiseProduct(l))/J;
        Eigen::Index direction;a.max_effective_tension=principal_cauchy.maxCoeff(&direction);a.normal=svd.matrixV().col(direction);a.J=J;a.principal_stretches=svd.singularValues();return a;
    }
    Real tangentBound(const Matd&F,Real damage)const
    {
        Eigen::JacobiSVD<Matd> svd(F);return tangentBound(svd.singularValues(),F.determinant(),damage);
    }
    /** Reuse singular values already obtained by the energy evaluation. */
    Real tangentBound(const Vecd&l,Real J,Real damage)const
    {
        if(l.minCoeff()<=0||J<=0)throw std::runtime_error("Invalid tangent-bound deformation");
        return peridynamics::finiteVolumeBiotTangent(K,mu,l,J,damage);
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
    struct Entry {V9 basis,coefficient;Real weight,crack_mode=0;Vecd reference_difference,residual,stabilizing_force;};
    StdVec<StdVec<Entry>> projection;
    StdVec<LinearEnergyCrackBand::History> history;
    PDQuadraticReference &reference;FiniteVolumeBiotEnergy elastic;DirectionalCrackEnergy<FiniteVolumeBiotEnergy> crack_energy;LinearEnergyCrackBand fracture;
    bool directional_contact,live_fracture_support;StdVec<Vecd>reference_positions;StdVec<int>projection_degree;StdVec<Real>crack_projection_norm;StdVec<int>crack_projection_ready;
    Real stabilizer,horizon,stabilizer_dt=0;size_t projection_topology_version=~size_t(0);
    Real *damage,*dissipation_density,*stored_density,*hg_density,*kappa_field,*tangent_bound;int *material_failed;
    Vecd *crack_normal;
    PDEnergyCrackBand(BaseInnerRelation&r,PDQuadraticReference&ref,Real E,Real nu,Real ft,Real Gc,Real hg,Real delta,bool directional=false,bool live_support=false)
      :LocalDynamics(r.getSPHBody()),DataDelegateInner(r),reference(ref),elastic(E,nu),crack_energy(elastic),fracture(E,ft,Gc,2*delta,directional?crack_energy.initiationEnergy(ft):elastic.initiationEnergy(ft)),directional_contact(directional),live_fracture_support(live_support),horizon(delta),
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
        size_t count=particles_->TotalRealParticles();
        if(reference_positions.empty())reference_positions.assign(x,x+count);
        if(reference_positions.size()!=count)throw std::runtime_error("Energetic reference particle count changed");
        history.resize(count);projection.resize(count);projection_degree.resize(count);crack_projection_norm.assign(count,0);crack_projection_ready.assign(count,0);Real max_rate=0;
        for(size_t i=0;i<count;++i)
        {
            reference.checkIntact(i); // Reject new masks with stale gradients.
            const auto&nb=inner_configuration_[i];Real sum=0;Matd cov=Matd::Zero();
            for(size_t n=0;n<nb.current_size_;++n){if(!nb.bondLive_[n]){if(!live_fracture_support)throw std::runtime_error("Continuous damage requires intact reference bonds; enable the qualified live-support path");continue;}Real a=nb.W_ij_[n]*volume[nb.j_[n]];Vecd z=reference_positions[nb.j_[n]]-reference_positions[i];sum+=a;cov+=a*z*z.transpose();}
            if(sum<=0)throw std::runtime_error("No energetic live support; no particle was erased");cov/=sum;Eigen::SelfAdjointEigenSolver<Matd> ce(cov);if(ce.info()!=Eigen::Success||ce.eigenvalues().minCoeff()<=0)throw std::runtime_error("Invalid energetic projection covariance");
            Matd L=ce.eigenvectors()*ce.eigenvalues().cwiseSqrt().cwiseInverse().asDiagonal()*ce.eigenvectors().transpose();
            M9 moment=M9::Zero();projection[i].resize(nb.current_size_);
            for(size_t n=0;n<nb.current_size_;++n){auto&e=projection[i][n];e.reference_difference=reference_positions[nb.j_[n]]-reference_positions[i];e.basis.setZero();e.coefficient.setZero();e.weight=0;e.crack_mode=0;e.stabilizing_force.setZero();if(!nb.bondLive_[n])continue;Vecd z=L*e.reference_difference;e.basis<<z[0],z[1],z[2],z[0]*z[0],z[1]*z[1],z[2]*z[2],z[0]*z[1],z[0]*z[2],z[1]*z[2];
                e.weight=nb.W_ij_[n]*volume[nb.j_[n]]/sum;moment+=e.weight*e.basis*e.basis.transpose();}
            peridynamics::layered_detail::MomentInverse<9> inverse(moment);bool quadratic=inverse.rank==9;projection_degree[i]=quadratic?2:1;
            if(!quadratic&&!live_fracture_support)throw std::runtime_error("Intact energetic quadratic projection lost rank");
            Eigen::LDLT<M9> solve(moment);Matd affine_moment=moment.topLeftCorner<3,3>();Eigen::LDLT<Matd> affine(affine_moment);M9 check=M9::Zero(),expected=M9::Zero();if(quadratic)expected.setIdentity();else expected.topLeftCorner<3,3>().setIdentity();
            for(size_t n=0;n<projection[i].size();++n){auto&e=projection[i][n];if(!nb.bondLive_[n])continue;if(quadratic)e.coefficient=e.weight*solve.solve(e.basis);else{e.basis.tail<6>().setZero();e.coefficient.head<3>()=e.weight*affine.solve(e.basis.head<3>());}check+=e.basis*e.coefficient.transpose();}
            if((check-expected).norm()>1e-8)throw std::runtime_error("Energetic live polynomial projection failed");
        }
        // The least-squares residual Hessian is bounded above by the unprojected
        // weighted graph Laplacian, giving a conservative positive spectral bound.
        StdVec<Real> rate(count,0);
        for(size_t i=0;i<count;++i)for(size_t n=0;n<reference.rows[i].size();++n)
        {auto&e=reference.rows[i][n];Real k=volume[i]*stabilizer*projection[i][n].weight;rate[i]+=2*k/mass[i];rate[e.j]+=2*k/mass[e.j];}
        for(Real a:rate)max_rate=std::max(max_rate,a);stabilizer_dt=max_rate>0?.2/std::sqrt(max_rate):std::numeric_limits<Real>::infinity();projection_topology_version=reference.topologyVersion();
    }
    /** Condense a centered Heaviside jump from the existing polynomial
     * least-squares residual. This preserves stabilization of non-crack modes
     * without retaining a spring across a representable displacement jump.
     * q is orthogonal to every existing polynomial basis; no new fit is used. */
    void buildCrackProjection(size_t i,const Vecd&N)
    {
        auto sign=[&](const Entry&e){Real d=N.dot(e.reference_difference),tol=32*std::numeric_limits<Real>::epsilon()*e.reference_difference.norm();return std::abs(d)<=tol?Real(0):(d>0?Real(.5):Real(-.5));};
        V9 coefficient=V9::Zero();for(const auto&e:projection[i])coefficient+=sign(e)*e.coefficient;
        Real norm=0;for(auto&e:projection[i]){e.crack_mode=sign(e)-e.basis.dot(coefficient);norm+=e.weight*e.crack_mode*e.crack_mode;}
        if(norm<=64*std::numeric_limits<Real>::epsilon()){norm=0;for(auto&e:projection[i])e.crack_mode=0;}
        crack_projection_norm[i]=norm;crack_projection_ready[i]=1;
    }
    void validateProjectionTopology()const{if(projection_topology_version!=reference.topologyVersion())throw std::runtime_error("Energetic projection is stale after a fracture event; rebuild while preserving material history");}
    void evaluateDirectional(size_t i,bool advance)
    {
        validateProjectionTopology();
        reference.checkIntact(i);A39 A=A39::Zero();Matd f=Matd::Zero();
        for(size_t n=0;n<reference.rows[i].size();++n){auto&e=reference.rows[i][n];Vecd y=x[e.j]-x[i];A+=y*projection[i][n].coefficient.transpose();f+=y*e.gradient.transpose();}
        F[i]=f;auto trial=elastic.evaluate(f);Real intact=trial.positive+trial.negative,hg=0,old=history[i].damage;
        Vecd N=old>0?crack_normal[i]:trial.normal;
        for(size_t n=0;n<reference.rows[i].size();++n){auto&e=projection[i][n];e.residual=x[reference.rows[i][n].j]-x[i]-A*e.basis;hg+=.5*stabilizer*e.weight*e.residual.squaredNorm();}
        if(old==0||!crack_projection_ready[i])buildCrackProjection(i,N);
        Vecd jump=Vecd::Zero();for(const auto&e:projection[i])jump+=e.weight*e.crack_mode*e.residual;
        if(crack_projection_norm[i]>0)jump/=crack_projection_norm[i];Real cracked_hg=0;for(const auto&e:projection[i])cracked_hg+=.5*stabilizer*e.weight*(e.residual-jump*e.crack_mode).squaredNorm();
        Real retained=crack_energy.energy(f,N),release=intact-retained;
        Real tol=128*std::numeric_limits<Real>::epsilon()*elastic.E*(1+f.squaredNorm());
        if(release < -tol)throw std::runtime_error("Directional effective energy exceeds intact elastic energy");
        release=std::max(release,Real(0));Real hg_release=.5*stabilizer*crack_projection_norm[i]*jump.squaredNorm();
        if(advance)fracture.advance(history[i],release+hg_release,trial.max_effective_tension>1e-12*elastic.E);
        Real d=history[i].damage,g=1-d;if(old==0&&d>0)crack_normal[i]=N;
        Matd full=trial.positive_P+trial.negative_P;
        if(d>0){auto cracked=crack_energy.evaluate(f,N);P[i]=g*full+d*cracked.P;tangent_bound[i]=g*elastic.tangentBound(trial.principal_stretches,trial.J,0)+d*cracked.tangent_bound;}
        else{P[i]=full;tangent_bound[i]=elastic.tangentBound(trial.principal_stretches,trial.J,0);}
        stress[i]=P[i]*f.transpose()/trial.J;damage[i]=d;dissipation_density[i]=history[i].dissipation;kappa_field[i]=history[i].kappa;material_failed[i]=d==1;
        hg_density[i]=g*hg+d*cracked_hg;stored_density[i]=g*intact+d*retained+hg_density[i];
        for(auto&e:projection[i])e.stabilizing_force=volume[i]*stabilizer*e.weight*(e.residual-d*jump*e.crack_mode);
    }
    void evaluate(size_t i,bool advance=true)
    {
        validateProjectionTopology();
        if(directional_contact){evaluateDirectional(i,advance);return;}
        reference.checkIntact(i);A39 A=A39::Zero();Matd f=Matd::Zero();
        for(size_t n=0;n<reference.rows[i].size();++n){auto&e=reference.rows[i][n];Vecd y=x[e.j]-x[i];A+=y*projection[i][n].coefficient.transpose();f+=y*e.gradient.transpose();}
        F[i]=f;auto trial=elastic.evaluate(f);Real hg=0;
        for(size_t n=0;n<reference.rows[i].size();++n){auto&e=projection[i][n];e.residual=x[reference.rows[i][n].j]-x[i]-A*e.basis;hg+=.5*stabilizer*e.weight*e.residual.squaredNorm();}
        Real old=history[i].damage;
        if(advance)fracture.advance(history[i],trial.positive+hg,trial.max_effective_tension>1e-12*elastic.E);
        if(old==0&&history[i].damage>0)crack_normal[i]=trial.normal;
        tangent_bound[i]=elastic.tangentBound(trial.principal_stretches,trial.J,history[i].damage);
        Real g=1-history[i].damage;P[i]=g*trial.positive_P+trial.negative_P;stress[i]=P[i]*f.transpose()/trial.J;
        damage[i]=history[i].damage;dissipation_density[i]=history[i].dissipation;kappa_field[i]=history[i].kappa;material_failed[i]=damage[i]==1;
        stored_density[i]=g*(trial.positive+hg)+trial.negative;hg_density[i]=g*hg;
        for(auto&e:projection[i])e.stabilizing_force=volume[i]*g*stabilizer*e.weight*e.residual;
    }
    void forces(size_t i)
    {
        validateProjectionTopology();reference.checkIntact(i);
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
