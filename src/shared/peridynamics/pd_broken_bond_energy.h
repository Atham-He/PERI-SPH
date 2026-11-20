#pragma once
/** Full existing elastic energy on the live, one-sided NOSB-PD reference.
 * An already cut surface carries no material or stabilization bridge.
 * Fracture work and cohesive/contact energies belong to their own surface
 * laws; this class neither duplicates Gc nor degrades a second bulk energy.
 */
#include "pd_energy_crack_band.h"
#include "pd_layered_fracture_reference.h"
#include "pd_spatial_derivatives.h"
namespace SPH::solid_dynamics {
class PDBrokenBondEnergy : public LocalDynamics,public DataDelegateInner {
  public:
    using V9=Eigen::Matrix<Real,9,1>;using M9=Eigen::Matrix<Real,9,9>;using A39=Eigen::Matrix<Real,3,9>;
    struct Entry {V9 basis=V9::Zero(),coefficient=V9::Zero();Real weight=0;Vecd residual=Vecd::Zero(),force=Vecd::Zero();};
    PDLayeredFractureReference &reference;FiniteVolumeBiotEnergy elastic;
    StdVec<StdVec<Entry>> projection;StdVec<int>residual_order;
    Real stabilizer,stabilizer_dt=0,maximum_projection_error=0;size_t version=~size_t(0);
    Real *stored_density,*hg_density,*tangent_bound;Matd *F,*P,*stress;Vecd *internal,*hourglass;
    PDBrokenBondEnergy(BaseInnerRelation&r,PDLayeredFractureReference&ref,Real E,Real nu,Real hg,Real horizon,
                       std::unique_ptr<peridynamics::SpatialDerivativeBatch> derivatives = nullptr)
      :LocalDynamics(r.getSPHBody()),DataDelegateInner(r),reference(ref),elastic(E,nu),
       x(particles_->getVariableDataByName<Vecd>("Position")),volume(particles_->getVariableDataByName<Real>("VolumetricMeasure")),
       mass(particles_->getVariableDataByName<Real>("Mass")),acc(particles_->getVariableDataByName<Vecd>("PDAcceleration")),
       F(particles_->getVariableDataByName<Matd>("DeformationGradient")),P(particles_->getVariableDataByName<Matd>("PDPK1")),
       stress(particles_->getVariableDataByName<Matd>("PDStress")),internal(particles_->getVariableDataByName<Vecd>("PDInternalForce")),
       hourglass(particles_->getVariableDataByName<Vecd>("HourglassForce")), derivatives_(std::move(derivatives)) {
        if(hg<0||horizon<=0)throw std::invalid_argument("Invalid live-bond stabilizer");stabilizer=hg*elastic.K/(horizon*horizon);
        stored_density=particles_->registerStateVariableData<Real>("ShellStoredEnergyDensity");hg_density=particles_->registerStateVariableData<Real>("ShellStabilizationEnergyDensity");tangent_bound=particles_->registerStateVariableData<Real>("ShellTangentModulusBound");
    }
    void rebuild() {
        derivative_stencil_current_ = false;
        device_residual_current_ = false;
        if(reference.topologyPending())throw std::runtime_error("Rebuild the layered gradient before its live residual");
        using namespace peridynamics::layered_detail;size_t count=particles_->TotalRealParticles();projection.resize(count);residual_order.resize(count);maximum_projection_error=0;StdVec<Real>rate(count,0);
        for(size_t i=0;i<count;++i){
            const auto&nb=reference.neighbors(i);projection[i].assign(nb.current_size_,Entry{});Real total=0,original_total=0;for(size_t n=0;n<nb.current_size_;++n)original_total+=nb.W_ij_[n]*volume[nb.j_[n]];Matd covariance=Matd::Zero();
            for(size_t n=0;n<nb.current_size_;++n)if(reference.live(i,n)){size_t j=nb.j_[n];Vecd z=reference.X[j]-reference.X[i];Real a=nb.W_ij_[n]*volume[j];total+=a;covariance+=a*z*z.transpose();}
            if(total<=0)throw std::runtime_error("No live bulk support; unsupported fragments require a separate verified kinematic model");
            covariance/=total;Eigen::SelfAdjointEigenSolver<Matd>ce(covariance);if(ce.info()!=Eigen::Success||ce.eigenvalues().minCoeff()<=0)throw std::runtime_error("Live residual covariance lacks 3D affine support");
            Matd L=ce.eigenvectors()*ce.eigenvalues().cwiseSqrt().cwiseInverse().asDiagonal()*ce.eigenvectors().transpose();M9 moment=M9::Zero();
            for(size_t n=0;n<nb.current_size_;++n)if(reference.live(i,n)){auto&e=projection[i][n];Vecd z=L*(reference.X[nb.j_[n]]-reference.X[i]);e.basis<<z[0],z[1],z[2],z[0]*z[0],z[1]*z[1],z[2]*z[2],z[0]*z[1],z[0]*z[2],z[1]*z[2];e.weight=nb.W_ij_[n]*volume[nb.j_[n]]/original_total;moment+=e.weight*e.basis*e.basis.transpose();}
            // The LS fitted values are unique even when the quadratic
            // coefficients are not. Retain their represented span; do not
            // introduce a stiffer affine residual on damaged supports.
            MomentInverse<9> inverse(moment);residual_order[i]=2;
            for(size_t n=0;n<nb.current_size_;++n)if(reference.live(i,n)){auto&e=projection[i][n];e.coefficient=e.weight*inverse.inverse*e.basis;Real k=volume[i]*stabilizer*e.weight;rate[i]+=2*k/mass[i];rate[nb.j_[n]]+=2*k/mass[nb.j_[n]];}
            Real error=(moment*inverse.inverse*moment-moment).norm()/moment.norm();maximum_projection_error=std::max(maximum_projection_error,error);if(error>1e-8)throw std::runtime_error("Live residual least-squares projection failed");
        }
        Real largest=0;for(Real a:rate)largest=std::max(largest,a);stabilizer_dt=largest>0?.2/std::sqrt(largest):std::numeric_limits<Real>::infinity();version=reference.topology_version;
    }
    bool hasDeviceDerivatives() const noexcept { return derivatives_ != nullptr; }
    // Batched evaluation is required when a device backend owns the force field.
    void evaluate(size_t i)
    {
        if (derivatives_)
            throw std::logic_error("Use evaluateAll for the device live residual");
        evaluatePoint(i, nullptr);
    }

    void evaluateAll()
    {
        if (version != reference.topology_version || reference.topologyPending())
            throw std::runtime_error("Live elastic residual not rebuilt after bond topology change");
        const size_t count = projection.size();
        if (derivatives_)
            evaluateDeviceDerivatives(count);
        tbb::parallel_for(size_t(0), count, [&](size_t i)
        {
            evaluatePoint(i, derivatives_ ? &device_derivatives_[i] : nullptr);
        });
        device_residual_current_ = derivatives_ != nullptr;
    }

  private:
    void evaluatePoint(size_t i, const peridynamics::SpatialDerivatives3 *derivatives)
    {
        if (version != reference.topology_version)
            throw std::runtime_error("Live elastic residual not rebuilt after bond topology change");
        reference.checkIntact(i);
        Matd f = Matd::Zero();
        Real hg = 0;
        if (derivatives)
        {
            for (int a = 0; a < 3; ++a)
                for (int b = 0; b < 3; ++b)
                    f(a, b) = derivatives->gradient[a][b];
            hg = device_residual_.energy_density[i];
        }
        else
        {
            A39 A = A39::Zero();
            for (size_t n = 0; n < reference.rows[i].size(); ++n)
            {
                const auto &e = reference.rows[i][n];
                const Vecd y = x[e.j] - x[i];
                A += y * projection[i][n].coefficient.transpose();
                f += y * e.gradient.transpose();
            }
            for (size_t n = 0; n < projection[i].size(); ++n)
            {
                auto &e = projection[i][n];
                e.residual = x[reference.rows[i][n].j] - x[i] - A * e.basis;
                hg += .5 * stabilizer * e.weight * e.residual.squaredNorm();
                e.force = volume[i] * stabilizer * e.weight * e.residual;
            }
        }
        F[i] = f;
        const auto state = elastic.evaluate(f);
        P[i] = state.positive_P + state.negative_P;
        stress[i] = P[i] * f.transpose() / state.J;
        hg_density[i] = hg;
        stored_density[i] = state.positive + state.negative + hg;
        tangent_bound[i] = elastic.tangentBound(state.principal_stretches, state.J, 0);
    }
  public:
    void forces(size_t i)
    {
        if (derivatives_)
            throw std::logic_error("Use recomputeForces for device live forces");
        assembleParticleForce(i, nullptr);
    }

    void recomputeForces()
    {
        if (version != reference.topology_version || reference.topologyPending())
            throw std::runtime_error("Live force assembly is stale after a fracture event");
        const auto count = projection.size();
        if (derivatives_)
        {
            if (!device_residual_current_)
                throw std::logic_error("Evaluate the device live residual before force assembly");
            assembleDeviceStressForce(count);
        }
        tbb::parallel_for(size_t(0), count, [&](size_t i)
        {
            acc[i].setZero();
            assembleParticleForce(i, derivatives_ ? &device_stress_force_[i] : nullptr);
        });
    }
    Real energy()const{Real result=0;for(size_t i=0;i<projection.size();++i)result+=volume[i]*stored_density[i];return result;}
    Real stabilizationEnergy()const{Real result=0;for(size_t i=0;i<projection.size();++i)result+=volume[i]*hg_density[i];return result;}
    Real finiteStrainTimeStep()const{Real bound=3*elastic.K;for(size_t i=0;i<projection.size();++i)bound=std::max(bound,tangent_bound[i]);return reference.elastic_dt_bound*std::sqrt(3*elastic.K/bound);}
  private:
    Vecd*x;Real*volume,*mass;Vecd*acc;
    std::unique_ptr<peridynamics::SpatialDerivativeBatch> derivatives_;
    bool derivative_stencil_current_ = false;
    bool device_residual_current_ = false;
    std::vector<peridynamics::ParticleVector3> device_field_;
    std::vector<peridynamics::SpatialDerivatives3> device_derivatives_, device_force_coefficients_;
    std::vector<peridynamics::ParticleVector3> device_stress_force_;
    peridynamics::WeightedDifferenceResult device_residual_;

    void assembleParticleForce(size_t i, const peridynamics::ParticleVector3 *adjoint)
    {
        if(version!=reference.topology_version||reference.topologyPending())throw std::runtime_error("Live force assembly is stale after a fracture event");
        reference.checkIntact(i);
        Vecd f = Vecd::Zero(), h = Vecd::Zero();
        if (adjoint)
        {
            f = -Vecd(adjoint->component[0], adjoint->component[1], adjoint->component[2]);
            const auto &residual_force = device_residual_.forces[i].component;
            h = Vecd(residual_force[0], residual_force[1], residual_force[2]);
        }
        else
        {
            for (size_t n = 0; n < reference.rows[i].size(); ++n)
            {
                const auto &e = reference.rows[i][n];
                f += volume[i] * P[i] * e.gradient - volume[e.j] * P[e.j] * reference.rows[e.j][e.reverse].gradient;
                h += projection[i][n].force - projection[e.j][e.reverse].force;
            }
        }
        internal[i] = f;
        hourglass[i] = h;
        acc[i] += (f + h) / mass[i];
    }

    void prepareDerivativeStencil(size_t count)
    {
        if (!derivative_stencil_current_)
        {
            peridynamics::SpatialDerivativeMap map{};
            map.particle_count = count;
            map.offsets.push_back(0);
            for (size_t i = 0; i < count; ++i)
            {
                map.origins.push_back(i);
                for (size_t n = 0; n < reference.rows[i].size(); ++n)
                {
                    if (!reference.live(i, n))
                        continue;
                    const auto &e = reference.rows[i][n];
                    const auto &r = projection[i][n];
                    map.neighbors.push_back(e.j);
                    map.gradients.push_back({{e.gradient[0], e.gradient[1], e.gradient[2]}});
                    peridynamics::WeightedDifferenceEntry weights{};
                    for (int term = 0; term < 9; ++term)
                    {
                        weights.basis[term] = r.basis[term];
                        weights.coefficient[term] = r.coefficient[term];
                    }
                    weights.energy_weight = .5 * stabilizer * r.weight;
                    weights.force_weight = volume[i] * stabilizer * r.weight;
                    map.weighted_differences.push_back(weights);
                }
                map.offsets.push_back(map.neighbors.size());
            }
            derivatives_->configure(map);
            derivative_stencil_current_ = true;
        }
    }

    void assembleDeviceStressForce(size_t count)
    {
        prepareDerivativeStencil(count);
        device_force_coefficients_.resize(count);
        for (size_t i = 0; i < count; ++i)
            for (int a = 0; a < 3; ++a)
                for (int b = 0; b < 3; ++b)
                    device_force_coefficients_[i].gradient[a][b] = volume[i] * P[i](a, b);
        derivatives_->assembleAdjoint(device_force_coefficients_, device_stress_force_);
        if (device_stress_force_.size() != count)
            throw std::logic_error("Spatial adjoint backend returned an incomplete live PD force field");
    }

    void evaluateDeviceDerivatives(size_t count)
    {
        prepareDerivativeStencil(count);
        device_field_.resize(count);
        for (size_t i = 0; i < count; ++i)
            device_field_[i] = {{x[i][0], x[i][1], x[i][2]}};
        device_residual_current_ = false;
        derivatives_->evaluate(device_field_, device_derivatives_, &device_residual_);
        if (device_derivatives_.size() != count || device_residual_.energy_density.size() != count ||
            device_residual_.forces.size() != count)
            throw std::logic_error("Spatial derivative backend returned an incomplete live PD field");
    }
};
}
