/* Optional finite-strain S-NSNI-based RKPM operator.
 * Wang et al. 10.1007/s00466-024-02494-0, Eqs.57,59-60,72-76.
 * Adaptations: native multiplicative Hencky DPC, independent centroid material
 * sites, complete pushed-forward cell moments. No added hourglass penalty,
 * force scaling, yield clipping or artificial stress is used.
 */
#pragma once
#include "solid_snsni_reference.h"
#include "dpc_snsni_adapter.h"
#include "pd_execution.h"
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
    SnsniPlasticSolid(SnsniReference &r, const dpc::Model &m, const peridynamics::PeriExecution &execution)
        : SnsniPlasticSolid(r, m, execution.makeSpatialDerivatives(), execution.makeDpcPoints(m),
                            execution.makeSpatialDerivatives()) {}

    SnsniPlasticSolid(SnsniReference&r,const dpc::Model&m,
                      std::unique_ptr<peridynamics::SpatialDerivativeBatch> derivatives = nullptr,
                      std::unique_ptr<dpc::PointBatch> material_points = nullptr,
                      std::unique_ptr<peridynamics::SpatialDerivativeBatch> sampling_derivatives = nullptr)
      : reference(r), model(m), derivatives_(std::move(derivatives)), material_points_(std::move(material_points)),
        sampling_derivatives_(std::move(sampling_derivatives))
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
        if (derivatives_)
            configureDerivatives();
        if (sampling_derivatives_)
            configureSamplingDerivatives();
    }
    bool hasDeviceSamplingGradients() const noexcept { return sampling_derivatives_ != nullptr; }
    bool hasDeviceDerivatives() const noexcept { return derivatives_ != nullptr; }
    bool hasDeviceMaterial() const noexcept { return material_points_ != nullptr; }
    void evaluate(const Field&d,Real rho0,Real cfl)
    {
        if(pending_commit)throw std::logic_error("Commit or explicitly discard pending DPC history before the next evaluation");
        const size_t n=history.size();
        if (derivatives_)
            evaluateDeviceDerivatives(d);
        if (material_points_)
            point_input_.resize(n);
        if (sampling_derivatives_)
            evaluateSamplingField(n, history_channels, [&](size_t l, double *values)
            {
                packSamplingTensor(history[l].inverse_Cp, values);
                values[tensor_channels]=history[l].kappa;
            });
        tbb::parallel_for(size_t(0),n,[&](size_t l){auto&h=displacement_hessian[l];
            if (derivatives_)
            {
                const auto &cell = device_derivatives_[l];
                for (int i = 0; i < 3; ++i)
                    for (int j = 0; j < 3; ++j)
                    {
                        F[l](i, j) = cell.gradient[i][j];
                        for (int a = 0; a < 3; ++a)
                            h[a](i, j) = cell.hessian[a][i][j];
                    }
            }
            else
            {
                for(auto&a:h)a.setZero();
                F[l]=Matd::Identity();
                for(const auto&e:reference.cells[l].rows){F[l]+=d[e.j]*e.gradient.transpose();for(int a=0;a<3;++a)h[a]+=d[e.j][a]*e.hessian;}
            }
            if(!F[l].allFinite()||F[l].determinant()<=0)throw std::runtime_error("Invalid S-NSNI cell Jacobian, no deletion or clamping");inverseF[l]=F[l].inverse();
            Triple derivativeF,derivativeCp;
            Vecd derivativeKappa;
            if (sampling_derivatives_)
            {
                derivativeCp=sampledTensorGradient(l, history_channels);
                const auto &kappa=sampled_gradients_[l*history_channels+tensor_channels].component;
                derivativeKappa=Vecd(kappa[0],kappa[1],kappa[2]);
            }
            else
            {
                for(auto &a:derivativeCp) a.setZero();
                derivativeKappa.setZero();
                // Difference form enforces the derivative of a constant exactly.
                for(const auto &e:sampling_gradient[l])
                {
                    const Matd delta=history[e.j].inverse_Cp-history[l].inverse_Cp;
                    const Real dk=history[e.j].kappa-history[l].kappa;
                    for(int k=0;k<3;++k) derivativeCp[k]+=e.gradient[k]*delta;
                    derivativeKappa+=dk*e.gradient;
                }
            }
            Triple cp_spatial;Vecd kappa_spatial=inverseF[l].transpose()*derivativeKappa;
            for(int k=0;k<3;++k){derivativeF[k].setZero();cp_spatial[k].setZero();
                for(int j=0;j<3;++j){cp_spatial[k]+=inverseF[l](j,k)*derivativeCp[j];for(int a=0;a<3;++a)derivativeF[k].row(a)+=inverseF[l](j,k)*h[a].col(j).transpose();}}
            if (material_points_)
                packMaterialInput(l,derivativeF,cp_spatial,kappa_spatial);
            else
            {
                pending[l]=dpc::snsniMaterialStep(model,F[l],history[l],derivativeF,cp_spatial,kappa_spatial);
                correction[l]=pending[l].correction;
            }
        });
        if (material_points_)
            evaluateDeviceMaterial(n,rho0,cfl*reference.dp);
        if (sampling_derivatives_)
            evaluateSamplingField(n, tensor_channels, [&](size_t l, double *values)
            {
                packSamplingTensor(correction[l], values);
            });
        tbb::parallel_for(size_t(0),n,[&](size_t l){const auto&cell=reference.cells[l];const auto&inv=inverseF[l];
            Triple derivative;
            if (sampling_derivatives_)
                derivative=sampledTensorGradient(l,tensor_channels);
            else
            {
                for(auto &a:derivative) a.setZero();
                for(const auto &e:sampling_gradient[l])
                {
                    const Matd delta=correction[e.j]-correction[l];
                    for(int k=0;k<3;++k) derivative[k]+=e.gradient[k]*delta;
                }
            }
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
        if (derivatives_)
            assembleDeviceInternalForce(n);
        else
            tbb::parallel_for(size_t(0),internal.size(),[&](size_t j){Vecd f=Vecd::Zero();for(auto[l,k]:reference.transpose[j]){const auto&e=reference.cells[l].rows[k];f+=force_gradient[l]*e.gradient;
                for(int a=0;a<3;++a)f[a]+=(force_hessian[l][a].array()*e.hessian.array()).sum();}internal[j]=f;});
        acoustic_dt=std::numeric_limits<Real>::infinity();
        if (material_points_)
        {
            for (const auto &point : point_output_)
                acoustic_dt=std::min(acoustic_dt,point.acoustic_dt);
        }
        else
        {
            for (size_t l=0;l<n;++l)
            {
                const auto dpc_parameters=model.parameters(pending[l].material.relative_density,pending[l].material.state.alpha);
                const Real local_step=dpc::acousticStep(F[l],inverseF[l],dpc_parameters,rho0,cfl*reference.dp);
                acoustic_dt=std::min(acoustic_dt,local_step);
            }
        }
        ++evaluations;pending_commit=true;
    }
    void commit()
    {
        if(!pending_commit)throw std::logic_error("S-NSNI material history committed twice");
        for(size_t l=0;l<history.size();++l){history[l]=pending[l].material.state;accumulated_plastic_work[l]+=pending[l].material.plastic_stress_work;}
        ++commits;pending_commit=false;
    }
    void discard(){if(!pending_commit)throw std::logic_error("No pending S-NSNI state");pending_commit=false;}

 private:
    std::unique_ptr<peridynamics::SpatialDerivativeBatch> derivatives_;
    std::vector<peridynamics::ParticleVector3> device_field_;
    std::vector<peridynamics::SpatialDerivatives3> device_derivatives_, device_force_coefficients_;
    std::vector<peridynamics::ParticleVector3> device_internal_force_;
    std::unique_ptr<dpc::PointBatch> material_points_;
    std::vector<dpc::SnsniPointInput> point_input_;
    std::vector<dpc::SnsniPointOutput> point_output_;
    std::unique_ptr<peridynamics::SpatialDerivativeBatch> sampling_derivatives_;
    std::vector<double> sampling_values_;
    std::vector<peridynamics::ParticleVector3> sampled_gradients_;
    static constexpr size_t tensor_channels=9, history_channels=tensor_channels+1;

    void configureSamplingDerivatives()
    {
        peridynamics::SpatialDerivativeMap map{};
        map.particle_count=sampling_gradient.size();
        map.build_adjoint=false;
        map.offsets.push_back(0);
        for (size_t l=0;l<sampling_gradient.size();++l)
        {
            map.origins.push_back(l);
            for (const auto &entry:sampling_gradient[l])
            {
                map.neighbors.push_back(entry.j);
                map.gradients.push_back({{entry.gradient[0],entry.gradient[1],entry.gradient[2]}});
            }
            map.offsets.push_back(map.neighbors.size());
        }
        sampling_derivatives_->configure(map);
    }

    template<class Fill>
    void evaluateSamplingField(size_t cells, size_t channels, const Fill &fill)
    {
        sampling_values_.resize(cells*channels);
        tbb::parallel_for(size_t(0),cells,[&](size_t l)
        {
            fill(l,sampling_values_.data()+l*channels);
        });
        sampling_derivatives_->evaluateChannelGradients(sampling_values_,channels,sampled_gradients_);
        if (sampled_gradients_.size()!=cells*channels)
            throw std::logic_error("Spatial backend returned an incomplete S-NSNI sampling field");
    }

    static void packSamplingTensor(const Matd &matrix, double *values)
    {
        for (int row=0;row<3;++row)
            for (int col=0;col<3;++col)
                values[3*row+col]=matrix(row,col);
    }

    Triple sampledTensorGradient(size_t cell, size_t channels) const
    {
        Triple gradient;
        for (int row=0;row<3;++row)
            for (int col=0;col<3;++col)
            {
                const auto &value=sampled_gradients_[cell*channels+3*row+col].component;
                for (int direction=0;direction<3;++direction)
                    gradient[direction](row,col)=value[direction];
            }
        return gradient;
    }

    void packMaterialInput(size_t cell,const Triple &grad_F,const Triple &grad_Cp,const Vecd &grad_kappa)
    {
        auto &input = point_input_[cell];
        const auto &old = history[cell];
        for (int row = 0; row < 3; ++row)
            for (int col = 0; col < 3; ++col)
            {
                input.material.deformation[row][col] = F[cell](row,col);
                input.material.history.inverse_Cp[row][col] = old.inverse_Cp(row,col);
                for (int direction=0;direction<3;++direction)
                {
                    input.deformation_gradient[direction][row][col]=grad_F[direction](row,col);
                    input.inverse_Cp_gradient[direction][row][col]=grad_Cp[direction](row,col);
                }
            }
        for (int direction=0;direction<3;++direction)
            input.kappa_gradient[direction]=grad_kappa[direction];
        input.material.history.kappa = old.kappa;
        input.material.history.alpha = old.alpha;
    }

    void evaluateDeviceMaterial(size_t cells,Real rho0,Real cfl_length)
    {
        material_points_->evaluateSnsni(point_input_,rho0,cfl_length,point_output_);
        if (point_output_.size() != cells)
            throw std::logic_error("DPC backend returned an incomplete material field");
        // No accepted history is mutated. Diagnose a failed cell with the
        // same deformation and original history as the CPU model.
        for (size_t l = 0; l < cells; ++l)
        {
            const auto &point = point_output_[l].material;
            if (point.failure != dpc::PointFailure::None)
            {
                dpc::Result diagnostic;
                if (point.failure == dpc::PointFailure::BadReturn)
                    diagnostic = unpackMaterialResult(point);
                model.throwPointError(point.failure, F[l], history[l], diagnostic);
            }
        }
        tbb::parallel_for(size_t(0), cells, [&](size_t l)
        {
            auto &step = pending[l];
            const auto &point=point_output_[l];
            step.material = unpackMaterialResult(point.material);
            for (int row=0;row<3;++row)
                for (int col=0;col<3;++col)
                {
                    step.trial.sigma(row,col)=point.trial_sigma[row][col];
                    for (int direction=0;direction<3;++direction)
                        step.trial.gradient[direction](row,col)=point.trial_sigma_gradient[direction][row][col];
                }
            step.correction = step.trial.sigma - step.material.sigma;
            correction[l] = step.correction;
        });
    }

    static dpc::Result unpackMaterialResult(const dpc::PointOutput &point)
    {
        dpc::Result result;
        for (int row = 0; row < 3; ++row)
            for (int col = 0; col < 3; ++col)
            {
                result.sigma(row,col) = point.sigma[row][col];
                result.state.inverse_Cp(row,col) = point.history.inverse_Cp[row][col];
            }
        result.state.kappa = point.history.kappa;
        result.state.alpha = point.history.alpha;
        result.p = point.p;
        result.q = point.q;
        result.relative_density = point.relative_density;
        result.plastic_shear = point.plastic_shear;
        result.plastic_stress_work = point.plastic_stress_work;
        result.elastic_energy = point.elastic_energy;
        result.yield_residual = point.yield_residual;
        result.effective_d = point.effective_d;
        result.surface = point.surface;
        return result;
    }

    void assembleDeviceInternalForce(size_t cells)
    {
        device_force_coefficients_.resize(cells);
        for (size_t l = 0; l < cells; ++l)
        {
            auto &coefficient = device_force_coefficients_[l];
            for (int a = 0; a < 3; ++a)
                for (int i = 0; i < 3; ++i)
                {
                    coefficient.gradient[a][i] = force_gradient[l](a, i);
                    for (int j = 0; j < 3; ++j)
                        coefficient.hessian[a][i][j] = force_hessian[l][a](i, j);
                }
        }
        derivatives_->assembleAdjoint(device_force_coefficients_, device_internal_force_);
        if (device_internal_force_.size() != internal.size())
            throw std::logic_error("Spatial adjoint backend returned an incomplete S-NSNI particle field");
        for (size_t i = 0; i < internal.size(); ++i)
        {
            const auto &force = device_internal_force_[i].component;
            internal[i] = Vecd(force[0], force[1], force[2]);
        }
    }


    void configureDerivatives()
    {
        peridynamics::SpatialDerivativeMap map{};
        map.particle_count = reference.X.size();
        map.offsets.push_back(0);
        map.add_identity = true;
        for (const auto &cell : reference.cells)
        {
            for (const auto &e : cell.rows)
            {
                map.neighbors.push_back(e.j);
                map.gradients.push_back({{e.gradient[0], e.gradient[1], e.gradient[2]}});
                peridynamics::ShapeHessian3 hessian{};
                for (int i = 0; i < 3; ++i)
                    for (int j = 0; j < 3; ++j)
                        hessian.derivative[i][j] = e.hessian(i, j);
                map.hessians.push_back(hessian);
            }
            map.offsets.push_back(map.neighbors.size());
        }
        derivatives_->configure(map);
    }

    void evaluateDeviceDerivatives(const Field &displacement)
    {
        device_field_.resize(displacement.size());
        for (size_t i = 0; i < displacement.size(); ++i)
            device_field_[i] = {{displacement[i][0], displacement[i][1], displacement[i][2]}};
        derivatives_->evaluate(device_field_, device_derivatives_);
        if (device_derivatives_.size() != history.size())
            throw std::logic_error("Spatial derivative backend returned an incomplete S-NSNI cell field");
    }
};

/** Explicit Newmark beta=0,gamma=1/2 (kick-drift-kick), constant lumped M.
 * Optional symmetric exact damping splits retain the pre-existing numerical
 * relaxation rate; the undamped scheme is the central-difference method.
 */
class SnsniCentralDifference
{
  public:
    using Field=StdVec<Vecd>;
    Field displacement,velocity,acceleration;
    const StdVec<Real> &mass;
    Real damping_loss=0;
    size_t steps=0;

    explicit SnsniCentralDifference(const StdVec<Real> &m)
        : displacement(m.size(),Vecd::Zero()), velocity(displacement), acceleration(displacement), mass(m) {}

    SnsniCentralDifference(const StdVec<Real> &m, const peridynamics::PeriExecution &execution)
        : SnsniCentralDifference(m)
    {
        device_time_step_=execution.makeSolidTimeStep(mass);
    }

    bool hasDeviceTimeStep() const noexcept { return device_time_step_!=nullptr; }

    Real kinetic() const
    {
        Real energy=0;
        for(size_t j=0;j<mass.size();++j)
            energy+=.5*mass[j]*velocity[j].squaredNorm();
        return energy;
    }

    template<class Evaluate>
    void step(Real dt, Real damping, Evaluate &&evaluate)
    {
        const Real factor=std::exp(-.5*damping*dt);
        if (device_time_step_)
            predictOnDevice(dt,factor,acceleration,peridynamics::SolidKickInput::Acceleration,displacement,velocity);
        else
            tbb::parallel_for(size_t(0),mass.size(),[&](size_t j)
            {
                velocity[j]*=factor;
                velocity[j]+=.5*dt*acceleration[j];
                displacement[j]+=dt*velocity[j];
            });

        Field force;
        evaluate(displacement,velocity,force);
        if (device_time_step_)
            finishOnDevice(dt,factor,velocity,force,Field{},mass.size(),peridynamics::SolidKickInput::Acceleration);
        else
            tbb::parallel_for(size_t(0),mass.size(),[&](size_t j)
            {
                acceleration[j]=force[j]/mass[j];
                velocity[j]+=.5*dt*acceleration[j];
                velocity[j]*=factor;
            });
        ++steps;
    }

    // Split form for intermediate contact operators and trial step rejection.
    Real predictForceStep(Real dt, Real factor, const Field &force, Field &next_displacement, Field &half_velocity)
    {
        if (device_time_step_)
            return predictOnDevice(dt,factor,force,peridynamics::SolidKickInput::Force,next_displacement,half_velocity);

        Real damped_kinetic=0;
        half_velocity=velocity;
        next_displacement=displacement;
        for (size_t j=0;j<mass.size();++j)
        {
            half_velocity[j]*=factor;
            damped_kinetic+=.5*mass[j]*half_velocity[j].squaredNorm();
            half_velocity[j]+=.5*dt*force[j]/mass[j];
            next_displacement[j]+=dt*half_velocity[j];
        }
        return damped_kinetic;
    }

    Real finishForceStep(Real dt, Real factor, const Field &half_velocity, const Field &force,
                         const Field &additional_force, size_t additional_start)
    {
        if (additional_start>mass.size() || additional_force.size()>mass.size()-additional_start)
            throw std::invalid_argument("Additional solid force lies outside the particle range");
        if (device_time_step_)
            return finishOnDevice(dt,factor,half_velocity,force,additional_force,additional_start,peridynamics::SolidKickInput::Force);

        Real undamped_kinetic=0;
        for (size_t j=0;j<mass.size();++j)
        {
            velocity[j]=half_velocity[j]+.5*dt*force[j]/mass[j];
            undamped_kinetic+=.5*mass[j]*velocity[j].squaredNorm();
            velocity[j]*=factor;
            const Vecd additional=additionalForce(additional_force,additional_start,j);
            acceleration[j]=(force[j]+additional)/mass[j];
        }
        return undamped_kinetic;
    }

  private:
    using Input=peridynamics::SolidKickInput;
    std::unique_ptr<peridynamics::SolidTimeStepBatch> device_time_step_;
    std::vector<peridynamics::SolidDriftInput> drift_input_;
    std::vector<peridynamics::SolidDriftResult> drift_output_;
    std::vector<peridynamics::SolidFinishInput> finish_input_;
    std::vector<peridynamics::SolidFinishResult> finish_output_;

    static peridynamics::ParticleVector3 pack(const Vecd &v) { return {{v[0],v[1],v[2]}}; }
    static Vecd unpack(const peridynamics::ParticleVector3 &v)
    {
        return Vecd(v.component[0],v.component[1],v.component[2]);
    }

    static Vecd additionalForce(const Field &force, size_t start, size_t particle)
    {
        if (particle>=start && particle-start<force.size())
            return force[particle-start];
        return Vecd::Zero();
    }

    Real predictOnDevice(Real dt, Real factor, const Field &drive, Input mode,
                         Field &next_displacement, Field &half_velocity)
    {
        const size_t count=mass.size();
        if (drive.size()!=count)
            throw std::invalid_argument("Solid drift drive does not match the particle count");
        drift_input_.resize(count);
        for (size_t j=0;j<count;++j)
            drift_input_[j]={pack(displacement[j]),pack(velocity[j]),pack(drive[j])};
        device_time_step_->predict(drift_input_,dt,factor,mode,drift_output_);
        if (drift_output_.size()!=count)
            throw std::logic_error("Solid time stepping backend returned an incomplete drift");
        next_displacement.resize(count);
        half_velocity.resize(count);
        Real damped_kinetic=0;
        for (size_t j=0;j<count;++j)
        {
            next_displacement[j]=unpack(drift_output_[j].displacement);
            half_velocity[j]=unpack(drift_output_[j].half_velocity);
            damped_kinetic+=drift_output_[j].damped_kinetic;
        }
        return damped_kinetic;
    }

    Real finishOnDevice(Real dt, Real factor, const Field &half_velocity, const Field &force,
                        const Field &additional_force, size_t additional_start, Input mode)
    {
        const size_t count=mass.size();
        if (half_velocity.size()!=count || force.size()!=count)
            throw std::invalid_argument("Solid final kick fields do not match the particle count");
        finish_input_.resize(count);
        for (size_t j=0;j<count;++j)
            finish_input_[j]={pack(half_velocity[j]),pack(force[j]),pack(additionalForce(additional_force,additional_start,j))};
        device_time_step_->finish(finish_input_,dt,factor,mode,finish_output_);
        if (finish_output_.size()!=count)
            throw std::logic_error("Solid time stepping backend returned an incomplete final kick");
        Real undamped_kinetic=0;
        for (size_t j=0;j<count;++j)
        {
            velocity[j]=unpack(finish_output_[j].velocity);
            acceleration[j]=unpack(finish_output_[j].acceleration);
            undamped_kinetic+=finish_output_[j].undamped_kinetic;
        }
        return undamped_kinetic;
    }
};
}
