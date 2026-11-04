/* Licensed under the Apache License, Version 2.0. */
#ifndef PD_BOUNDED_GEL_H
#define PD_BOUNDED_GEL_H

#include "pd_fixed_reference.h"
#include "pd_gel_affine_path.h"
#include "pd_many_contact.h"

namespace SPH::peridynamics
{
using PDGradientMatrix=Eigen::Matrix<Real,9,Eigen::Dynamic>;
inline PDGradientMatrix pdGradientMatrix(const FixedPDReference &ref,size_t point)
{
    const size_t n=ref.size();if(point>=n)throw std::invalid_argument("bounded PD gradient: point");
    PDGradientMatrix b=PDGradientMatrix::Zero(9,3*n);
    for(size_t k=0;k<n;++k)for(size_t a=0;a<3;++a)for(size_t c=0;c<3;++c)b(3*a+c,3*k+a)=ref.shape_gradient[point*n+k][c];
    return b;
}
inline StdVec<Real> flattenPD(const StdVec<Vecd> &values)
{
    StdVec<Real> result(3*values.size());for(size_t i=0;i<values.size();++i)for(size_t a=0;a<3;++a)result[3*i+a]=values[i][a];return result;
}
inline StdVec<Vecd> unflattenPD(const StdVec<Real> &values)
{
    if(values.empty() || values.size()%3)throw std::invalid_argument("bounded PD vector: dimensions");
    StdVec<Vecd> result(values.size()/3);for(size_t i=0;i<result.size();++i)for(size_t a=0;a<3;++a)
    {if(!std::isfinite(values[3*i+a]))throw std::invalid_argument("bounded PD vector: nonfinite");result[i][a]=values[3*i+a];}return result;
}

struct BoundedGelState
{
    StdVec<Vecd> position,velocity;
    StdVec<Matd> green,viscous,accepted_stress;
    StdVec<Real> dissipation,density,last_normal_impulse;
    Real time=0,last_dt=0,last_KKT=0;
    uint64_t steps=0,last_newton_iterations=0;
    GelHistory history(size_t i)const{return {green.at(i),viscous.at(i),dissipation.at(i)};}

    template<class Checkpoint>void registerCheckpoint(Checkpoint &cp)
    {
        cp.addArray("GelPosition",position);cp.addArray("GelVelocity",velocity);cp.addArray("GelGreen",green);cp.addArray("GelViscous",viscous);
        cp.addArray("GelPhysicalDissipationDensity",dissipation);cp.addArray("GelAcceptedCauchyStress",accepted_stress);cp.addArray("GelAcceptedDensity",density);
        cp.addArray("LastNormalImpulse",last_normal_impulse);cp.addScalar("AcceptedTime",time);cp.addScalar("AcceptedStep",steps);
        cp.addScalar("LastActualDt",last_dt);cp.addScalar("LastKKT",last_KKT);cp.addScalar("LastNewtonIterations",last_newton_iterations);
    }
};

struct BoundedGelModel
{
    FixedPDReference reference;
    GelParameters material{5,1,.5,.3};
    Real reference_density=1,nonaffine_eta=.1;
    StdVec<Real> node_mass;
    BoundedGelModel()=default;
    BoundedGelModel(FixedPDReference ref,GelParameters parameters,Real rho0,Real eta)
        :reference(std::move(ref)),material(parameters),reference_density(rho0),nonaffine_eta(eta)
    {
        for(Real volume:reference.volume)node_mass.push_back(rho0*volume);validate();
    }
    void validate()const
    {
        reference.validate();material.validate();
        if(3*reference.size()>128 || node_mass.size()!=reference.size() || !std::isfinite(reference_density) || reference_density<=0 ||
            !std::isfinite(nonaffine_eta) || nonaffine_eta<0)throw std::invalid_argument("bounded gel model: dimensions/parameters");
        for(size_t i=0;i<node_mass.size();++i)if(!std::isfinite(node_mass[i]) || node_mass[i]<=0 ||
            std::abs(node_mass[i]-reference_density*reference.volume[i])>1e-12*node_mass[i])
            throw std::invalid_argument("bounded gel model: physical mass mismatch");
    }
    BoundedGelState initialState()const
    {
        validate();BoundedGelState state;state.position=reference.position;state.velocity.assign(reference.size(),Vecd::Zero());
        for(size_t i=0;i<reference.size();++i)
        {
            const Matd f=reference.deformation(i,state.position);GelHistory h;h.green=gelDeviatoricGreen(f);
            state.green.push_back(h.green);state.viscous.push_back(h.viscous);state.dissipation.push_back(0);
            state.accepted_stress.push_back(gelResponse(f,material,h).cauchy);state.density.push_back(reference_density/f.determinant());
        }
        return state;
    }
    template<class Checkpoint>void registerCheckpoint(Checkpoint &cp)
    {
        reference.registerCheckpoint(cp);cp.addArray("PhysicalNodeMass",node_mass);cp.addScalar("ReferenceDensity",reference_density);
        cp.addScalar("GelBulk",material.bulk);cp.addScalar("GelEquilibriumShear",material.equilibrium_shear);
        cp.addScalar("GelBranchShear",material.branch_shear);cp.addScalar("GelRelaxationTime",material.relaxation_time);
        cp.addScalar("NonaffineEta",nonaffine_eta);
    }
};

struct BoundedGelEvaluation
{
    NonlinearContactEquation equation;
    StdVec<GelAffinePathResult> point;
    StdVec<Matd> final_deformation;
    StdVec<Vecd> material_force,stabilization_force,internal_force;
    Real material_stored_old=0,material_stored_new=0,material_work=0,physical_material_loss=0;
    Real stabilization_stored_old=0,stabilization_stored_new=0,transpose_work_residual=0,stabilization_work_residual=0,jacobian_asymmetry=0;
};

/** Immutable accepted state and fixed reference throughout nonlinear trials. */
struct BoundedGelStepContext
{
    const BoundedGelModel &model;const BoundedGelState &old;Real dt,total_mass=0;
    StdVec<Matd> initial_deformation;
    BoundedGelStepContext(const BoundedGelModel &m,const BoundedGelState &s,Real step):model(m),old(s),dt(step)
    {
        model.validate();const size_t n=model.reference.size();
        if(!std::isfinite(dt) || dt<=0 || !std::isfinite(s.time) || s.time<0 || s.steps>=98304 || s.position.size()!=n || s.velocity.size()!=n ||
            s.green.size()!=n || s.viscous.size()!=n || s.dissipation.size()!=n || s.density.size()!=n || s.accepted_stress.size()!=n)
            throw std::invalid_argument("bounded gel state: dimensions/time");
        for(size_t i=0;i<n;++i)
        {
            if(!s.position[i].allFinite() || !s.velocity[i].allFinite())throw std::invalid_argument("bounded gel state: nonfinite kinematics");
            const Matd f=model.reference.deformation(i,s.position);const auto response=gelResponse(f,model.material,s.history(i));
            if((s.green[i]-gelDeviatoricGreen(f)).norm()>1e-12*(1+s.green[i].norm()) ||
                !std::isfinite(s.density[i]) || std::abs(s.density[i]-model.reference_density/f.determinant())>1e-10*model.reference_density ||
                !s.accepted_stress[i].allFinite() || (s.accepted_stress[i]-response.cauchy).norm()>1e-10*(1+response.cauchy.norm()))
                throw std::invalid_argument("bounded gel state: inconsistent accepted material cursor");
            initial_deformation.push_back(f);total_mass+=model.node_mass[i];
        }
    }
    bool inDomain(const StdVec<Real> &velocity)const
    {
        if(velocity.size()!=3*old.position.size())return false;
        StdVec<Vecd> q=old.position;for(size_t i=0;i<q.size();++i)for(size_t a=0;a<3;++a)
        {q[i][a]+=dt*velocity[3*i+a];if(!std::isfinite(q[i][a]))return false;}
        for(size_t i=0;i<q.size();++i)if(gelAffineMinimumDeterminant(initial_deformation[i],model.reference.deformation(i,q)-initial_deformation[i])<.05)return false;
        return true;
    }
    BoundedGelEvaluation evaluate(const StdVec<Real> &velocity)const
    {
        const auto v=unflattenPD(velocity);const auto &ref=model.reference;const size_t n=ref.size(),dofs=3*n;
        if(v.size()!=n)throw std::invalid_argument("bounded gel trial: dimensions");
        StdVec<Vecd> next=old.position,average=old.position;
        for(size_t i=0;i<n;++i){next[i]+=dt*v[i];average[i]+=.5*dt*v[i];}
        BoundedGelEvaluation out;out.material_force.assign(n,Vecd::Zero());
        Eigen::MatrixXd j=Eigen::MatrixXd::Zero(dofs,dofs);Eigen::VectorXd force=Eigen::VectorXd::Zero(dofs);
        for(size_t i=0;i<n;++i)
        {
            const Matd final=ref.deformation(i,next);out.final_deformation.push_back(final);
            out.point.push_back(advanceGelAffinePath(initial_deformation[i],final,dt,model.material,old.history(i)));
            const auto &p=out.point.back();const auto b=pdGradientMatrix(ref,i);
            Eigen::Matrix<Real,9,1> stress;Eigen::Matrix<Real,9,9> tangent;
            for(size_t a=0;a<9;++a)
            {
                stress[a]=p.average_pk1(a/3,a%3);
                for(size_t c=0;c<9;++c)tangent(a,c)=p.endpoint_tangent[c](a/3,a%3);
            }
            force.noalias()-=ref.volume[i]*b.transpose()*stress;
            j.noalias()+=dt*dt*ref.volume[i]*b.transpose()*tangent*b;
            out.material_stored_old+=ref.volume[i]*p.stored_old;out.material_stored_new+=ref.volume[i]*p.stored_new;
            out.material_work+=ref.volume[i]*p.work;out.physical_material_loss+=ref.volume[i]*p.physical_dissipation;
        }
        const Real eta_mu=model.nonaffine_eta*model.material.equilibrium_shear;
        out.stabilization_force=nonaffineForce(ref,average,eta_mu);out.internal_force.resize(n);
        out.stabilization_stored_old=nonaffineEnergy(ref,old.position,eta_mu);out.stabilization_stored_new=nonaffineEnergy(ref,next,eta_mu);
        out.equation.momentum.resize(dofs);out.equation.jacobian.assign(dofs,StdVec<Real>(dofs));
        long double material_work=0,stabilization_work=0;
        for(size_t i=0;i<n;++i)for(size_t a=0;a<3;++a)
        {
            const size_t row=3*i+a;out.material_force[i][a]=force[row];out.internal_force[i][a]=force[row]+out.stabilization_force[i][a];
            out.equation.momentum[row]=model.node_mass[i]*(v[i][a]-old.velocity[i][a])-dt*out.internal_force[i][a];
            material_work-=static_cast<long double>(dt)*force[row]*v[i][a];stabilization_work-=static_cast<long double>(dt)*out.stabilization_force[i][a]*v[i][a];
            j(row,row)+=model.node_mass[i];
            for(size_t k=0;k<n;++k)j(row,3*k+a)+=.5*dt*dt*eta_mu*ref.nonaffine_unit[i*n+k];
        }
        for(size_t i=0;i<dofs;++i)for(size_t k=0;k<dofs;++k)out.equation.jacobian[i][k]=j(i,k);
        out.transpose_work_residual=material_work-out.material_work;
        out.stabilization_work_residual=stabilization_work-(out.stabilization_stored_new-out.stabilization_stored_old);
        out.jacobian_asymmetry=j.norm()>0?(j-j.transpose()).norm()/j.norm():0;
        if(!j.allFinite() || !std::isfinite(out.transpose_work_residual) || !std::isfinite(out.stabilization_work_residual))
            throw std::runtime_error("bounded gel trial: nonfinite assembled operator");
        return out;
    }
};

struct BoundedGelStepResult
{
    bool accepted=false;
    ManyContactResult solve;
    BoundedGelState next;
    BoundedGelEvaluation evaluation;
    StdVec<Vecd> contact_impulse,momentum_residual;
    std::array<Vecd,2> particle_plate_impulse{Vecd::Zero(),Vecd::Zero()};
    std::array<Real,2> plate_work{0,0};
    Real kinetic_old=0,kinetic_new=0,algorithmic_loss=0,normal_contact_loss=0,energy_residual=0;
    Vecd total_momentum_residual=Vecd::Zero(),internal_force_sum=Vecd::Zero(),angular_drift_minus_external=Vecd::Zero(),internal_old_torque=Vecd::Zero();
};

inline BoundedGelStepResult advanceBoundedGel(const BoundedGelModel &model,const BoundedGelState &old,Real dt,
    const StdVec<DenseContactConstraint> &rows,const StdVec<size_t> &plate,const ManyContactSettings &settings={})
{
    const BoundedGelStepContext context(model,old,dt);const size_t n=model.reference.size();
    if(plate.size()!=rows.size())throw std::invalid_argument("bounded gel contacts: plate IDs");
    for(size_t id:plate)if(id>1)throw std::invalid_argument("bounded gel contacts: unknown plate");
    BoundedGelStepResult r;r.next=old;
    r.solve=solveManyNonsymmetricContact(flattenPD(old.position),dt,context.total_mass,rows,
        [&](const auto &v){return context.evaluate(v).equation;},[&](const auto &v){return context.inDomain(v);},settings);
    if(!r.solve.converged)return r;
    const Real impulse_scale=context.total_mass*settings.velocity_scale;
    r.evaluation=context.evaluate(r.solve.velocity);r.next.velocity=unflattenPD(r.solve.velocity);
    r.contact_impulse.assign(n,Vecd::Zero());r.momentum_residual.assign(n,Vecd::Zero());
    for(size_t s=0;s<rows.size();++s)
    {
        const Real impulse=r.solve.normal_impulse[s];const auto &row=rows[s];
        r.plate_work[plate[s]]+=impulse*(row.new_offset-row.old_offset)/dt;
        r.normal_contact_loss+=impulse*(denseContactDot(row.coefficients,flattenPD(old.position))-row.old_offset)/dt;
    }
    for(size_t id=0;id<2;++id)for(size_t a=0;a<3;++a)
    {
        StdVec<Real> coefficients,impulses;
        for(size_t s=0;s<rows.size();++s)if(plate[s]==id)for(size_t i=0;i<n;++i)
        {coefficients.push_back(rows[s].coefficients[3*i+a]);impulses.push_back(r.solve.normal_impulse[s]);}
        const auto balance=exactContactBalance(0,coefficients,impulses,impulse_scale);
        if(!balance.within_budget)throw std::runtime_error("bounded gel: plate impulse precision budget");
        r.particle_plate_impulse[id][a]=-balance.normalized_value*impulse_scale;
    }
    for(size_t i=0;i<n;++i)
    {
        r.next.position[i]=old.position[i]+dt*r.next.velocity[i];const Real mass=model.node_mass[i];
        const Vecd dv=r.next.velocity[i]-old.velocity[i];
        r.kinetic_old+=.5*mass*old.velocity[i].squaredNorm();r.kinetic_new+=.5*mass*r.next.velocity[i].squaredNorm();r.algorithmic_loss+=.5*mass*dv.squaredNorm();
        for(size_t a=0;a<3;++a)
        {
            StdVec<Real> coefficients(rows.size());for(size_t s=0;s<rows.size();++s)coefficients[s]=rows[s].coefficients[3*i+a];
            const auto scatter=exactContactBalance(0,coefficients,r.solve.normal_impulse,impulse_scale);
            const auto residual=exactContactBalance(r.evaluation.equation.momentum[3*i+a],coefficients,r.solve.normal_impulse,impulse_scale);
            if(!scatter.within_budget || !residual.within_budget)throw std::runtime_error("bounded gel: precommit impulse precision budget");
            r.contact_impulse[i][a]=-scatter.normalized_value*impulse_scale;r.momentum_residual[i][a]=residual.normalized_value;
        }
        const auto &p=r.evaluation.point[i];r.next.green[i]=p.history.green;r.next.viscous[i]=p.history.viscous;r.next.dissipation[i]=p.history.dissipation;
        const auto response=gelResponse(r.evaluation.final_deformation[i],model.material,p.history);
        r.next.accepted_stress[i]=response.cauchy;r.next.density[i]=model.reference_density/r.evaluation.final_deformation[i].determinant();
        r.total_momentum_residual+=mass*dv-r.contact_impulse[i];r.internal_force_sum+=r.evaluation.internal_force[i];
        r.angular_drift_minus_external+=mass*(r.next.position[i].cross(r.next.velocity[i])-old.position[i].cross(old.velocity[i]))-old.position[i].cross(r.contact_impulse[i]);
        r.internal_old_torque+=old.position[i].cross(r.evaluation.internal_force[i]);
    }
    const auto &e=r.evaluation;
    r.energy_residual=r.plate_work[0]+r.plate_work[1]-(r.kinetic_new-r.kinetic_old+e.material_stored_new-e.material_stored_old+
        e.stabilization_stored_new-e.stabilization_stored_old+e.physical_material_loss+r.algorithmic_loss+r.normal_contact_loss);
    Real total_volume=0,length_scale=0;for(size_t i=0;i<n;++i){total_volume+=model.reference.volume[i];length_scale=std::max(length_scale,old.position[i].norm());}
    const Real energy_scale=std::max(impulse_scale*settings.velocity_scale,total_volume*std::max(model.material.bulk,model.material.equilibrium_shear));
    if(!std::isfinite(r.energy_residual) || std::abs(r.energy_residual)>1e-10*energy_scale ||
        std::abs(e.transpose_work_residual)>1e-10*energy_scale || std::abs(e.stabilization_work_residual)>1e-10*energy_scale ||
        r.normal_contact_loss < -1e-10*energy_scale)throw std::runtime_error("bounded gel: precommit energy/transpose work gate");
    for(size_t i=0;i<n;++i)for(size_t a=0;a<3;++a)if(!std::isfinite(r.momentum_residual[i][a]) || std::abs(r.momentum_residual[i][a])>1e-10)
        throw std::runtime_error("bounded gel: precommit coordinate momentum gate");
    if(r.total_momentum_residual.norm()>1e-10*impulse_scale || r.internal_force_sum.norm()>1e-10*impulse_scale/dt ||
        (r.angular_drift_minus_external-dt*r.internal_old_torque).norm()>1e-10*impulse_scale*std::max(length_scale,model.reference.horizon))
        throw std::runtime_error("bounded gel: precommit global/cursor ledger gate");
    r.next.time=old.time+dt;r.next.steps=old.steps+1;r.next.last_dt=dt;r.next.last_KKT=r.solve.normalized_KKT;
    r.next.last_normal_impulse=r.solve.normal_impulse;r.next.last_newton_iterations=r.solve.trace.empty()?0:r.solve.trace.back().iteration;
    r.accepted=true;return r;
}
} // namespace SPH::peridynamics
#endif
