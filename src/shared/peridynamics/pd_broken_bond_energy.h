#pragma once
/** Full existing elastic energy on the live, one-sided NOSB-PD reference.
 * An already cut surface carries no material or stabilization bridge.
 * Fracture work and cohesive/contact energies belong to their own surface
 * laws; this class neither duplicates Gc nor degrades a second bulk energy.
 */
#include "pd_energy_crack_band.h"
#include "pd_layered_fracture_reference.h"
namespace SPH::solid_dynamics {
class PDBrokenBondEnergy : public LocalDynamics,public DataDelegateInner {
  public:
    using V9=Eigen::Matrix<Real,9,1>;using M9=Eigen::Matrix<Real,9,9>;using A39=Eigen::Matrix<Real,3,9>;
    struct Entry {V9 basis=V9::Zero(),coefficient=V9::Zero();Real weight=0;Vecd residual=Vecd::Zero(),force=Vecd::Zero();};
    PDLayeredFractureReference &reference;FiniteVolumeBiotEnergy elastic;
    StdVec<StdVec<Entry>> projection;StdVec<int>residual_order;
    Real stabilizer,stabilizer_dt=0,maximum_projection_error=0;size_t version=~size_t(0);
    Real *stored_density,*hg_density,*tangent_bound;Matd *F,*P,*stress;Vecd *internal,*hourglass;
    PDBrokenBondEnergy(BaseInnerRelation&r,PDLayeredFractureReference&ref,Real E,Real nu,Real hg,Real horizon)
      :LocalDynamics(r.getSPHBody()),DataDelegateInner(r),reference(ref),elastic(E,nu),
       x(particles_->getVariableDataByName<Vecd>("Position")),volume(particles_->getVariableDataByName<Real>("VolumetricMeasure")),
       mass(particles_->getVariableDataByName<Real>("Mass")),acc(particles_->getVariableDataByName<Vecd>("PDAcceleration")),
       F(particles_->getVariableDataByName<Matd>("DeformationGradient")),P(particles_->getVariableDataByName<Matd>("PDPK1")),
       stress(particles_->getVariableDataByName<Matd>("PDStress")),internal(particles_->getVariableDataByName<Vecd>("PDInternalForce")),
       hourglass(particles_->getVariableDataByName<Vecd>("HourglassForce")) {
        if(hg<0||horizon<=0)throw std::invalid_argument("Invalid live-bond stabilizer");stabilizer=hg*elastic.K/(horizon*horizon);
        stored_density=particles_->registerStateVariableData<Real>("ShellStoredEnergyDensity");hg_density=particles_->registerStateVariableData<Real>("ShellStabilizationEnergyDensity");tangent_bound=particles_->registerStateVariableData<Real>("ShellTangentModulusBound");
    }
    void rebuild() {
        if(reference.topologyPending())throw std::runtime_error("Rebuild the layered gradient before its live residual");
        using namespace peridynamics::layered_detail;size_t count=particles_->TotalRealParticles();projection.resize(count);residual_order.resize(count);maximum_projection_error=0;StdVec<Real>rate(count,0);
        for(size_t i=0;i<count;++i){
            const auto&nb=reference.neighbors(i);projection[i].assign(nb.current_size_,Entry{});Real total=0;Matd covariance=Matd::Zero();
            for(size_t n=0;n<nb.current_size_;++n)if(reference.live(i,n)){size_t j=nb.j_[n];Vecd z=reference.X[j]-reference.X[i];Real a=nb.W_ij_[n]*volume[j];total+=a;covariance+=a*z*z.transpose();}
            if(total<=0)throw std::runtime_error("No live bulk support; unsupported fragments require a separate verified kinematic model");
            covariance/=total;Eigen::SelfAdjointEigenSolver<Matd>ce(covariance);if(ce.info()!=Eigen::Success||ce.eigenvalues().minCoeff()<=0)throw std::runtime_error("Live residual covariance lacks 3D affine support");
            Matd L=ce.eigenvectors()*ce.eigenvalues().cwiseSqrt().cwiseInverse().asDiagonal()*ce.eigenvectors().transpose();M9 moment=M9::Zero();
            for(size_t n=0;n<nb.current_size_;++n)if(reference.live(i,n)){auto&e=projection[i][n];Vecd z=L*(reference.X[nb.j_[n]]-reference.X[i]);e.basis<<z[0],z[1],z[2],z[0]*z[0],z[1]*z[1],z[2]*z[2],z[0]*z[1],z[0]*z[2],z[1]*z[2];e.weight=nb.W_ij_[n]*volume[nb.j_[n]]/total;moment+=e.weight*e.basis*e.basis.transpose();}
            MomentInverse<9> inverse(moment);bool quadratic=inverse.rank==9;residual_order[i]=quadratic?2:1;M9 expected=M9::Zero();if(quadratic)expected.setIdentity();else expected.topLeftCorner<3,3>().setIdentity();
            Matd affine=moment.topLeftCorner<3,3>();Eigen::LDLT<Matd>asolve(affine);M9 check=M9::Zero();
            for(size_t n=0;n<nb.current_size_;++n)if(reference.live(i,n)){auto&e=projection[i][n];if(quadratic)e.coefficient=e.weight*inverse.inverse*e.basis;else{e.basis.tail<6>().setZero();e.coefficient.head<3>()=e.weight*asolve.solve(e.basis.head<3>());}check+=e.basis*e.coefficient.transpose();Real k=volume[i]*stabilizer*e.weight;rate[i]+=2*k/mass[i];rate[nb.j_[n]]+=2*k/mass[nb.j_[n]];}
            Real error=(check-expected).norm();maximum_projection_error=std::max(maximum_projection_error,error);if(error>1e-8)throw std::runtime_error("Live residual polynomial reproduction failed");
        }
        Real largest=0;for(Real a:rate)largest=std::max(largest,a);stabilizer_dt=largest>0?.2/std::sqrt(largest):std::numeric_limits<Real>::infinity();version=reference.topology_version;
    }
    void evaluate(size_t i) {
        if(version!=reference.topology_version)throw std::runtime_error("Live elastic residual not rebuilt after bond topology change");reference.checkIntact(i);A39 A=A39::Zero();Matd f=Matd::Zero();
        for(size_t n=0;n<reference.rows[i].size();++n){auto&e=reference.rows[i][n];Vecd y=x[e.j]-x[i];A+=y*projection[i][n].coefficient.transpose();f+=y*e.gradient.transpose();}
        F[i]=f;auto state=elastic.evaluate(f);P[i]=state.positive_P+state.negative_P;stress[i]=P[i]*f.transpose()/state.J;Real hg=0;
        for(size_t n=0;n<projection[i].size();++n){auto&e=projection[i][n];e.residual=x[reference.rows[i][n].j]-x[i]-A*e.basis;hg+=.5*stabilizer*e.weight*e.residual.squaredNorm();e.force=volume[i]*stabilizer*e.weight*e.residual;}
        hg_density[i]=hg;stored_density[i]=state.positive+state.negative+hg;tangent_bound[i]=elastic.tangentBound(state.principal_stretches,state.J,0);
    }
    void forces(size_t i) {
        if(version!=reference.topology_version||reference.topologyPending())throw std::runtime_error("Live force assembly is stale after a fracture event");reference.checkIntact(i);
        Vecd f=Vecd::Zero(),h=Vecd::Zero();for(size_t n=0;n<reference.rows[i].size();++n){auto&e=reference.rows[i][n];f+=volume[i]*P[i]*e.gradient-volume[e.j]*P[e.j]*reference.rows[e.j][e.reverse].gradient;h+=projection[i][n].force-projection[e.j][e.reverse].force;}
        internal[i]=f;hourglass[i]=h;acc[i]+=(f+h)/mass[i];
    }
    Real energy()const{Real result=0;for(size_t i=0;i<projection.size();++i)result+=volume[i]*stored_density[i];return result;}
    Real stabilizationEnergy()const{Real result=0;for(size_t i=0;i<projection.size();++i)result+=volume[i]*hg_density[i];return result;}
    Real finiteStrainTimeStep()const{Real bound=3*elastic.K;for(size_t i=0;i<projection.size();++i)bound=std::max(bound,tangent_bound[i]);return reference.elastic_dt_bound*std::sqrt(3*elastic.K/bound);}
  private:Vecd*x;Real*volume,*mass;Vecd*acc;
};
}
