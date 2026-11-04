/* Licensed under the Apache License, Version 2.0. */
#ifndef PD_BOUNDED_GEL_FRICTION_H
#define PD_BOUNDED_GEL_FRICTION_H

#include "pd_bounded_normal_scene.h"
#include "pd_coulomb_contact.h"

namespace SPH::peridynamics
{
enum class ContactLedger : size_t
{
    LowerImpulseX=0,LowerImpulseY,LowerImpulseZ,UpperImpulseX,UpperImpulseY,UpperImpulseZ,
    LowerWork,UpperWork,MaxwellLoss,NormalLoss,FrictionLoss,AlgorithmicLoss,Count
};

inline void addCompensatedContactValue(Real x,Real &sum,Real &correction)
{
    const Real next=sum+x;
    correction+=std::abs(sum)>=std::abs(x) ? (sum-next)+x : (x-next)+sum;
    sum=next;
    if (!std::isfinite(sum) || !std::isfinite(correction))
        throw std::runtime_error("bounded gel friction: nonfinite cumulative ledger");
}

/** All persistent diagnostics are explicit state, including compensation.
 * These registration hooks alone are NOT a complete restart qualification.
 */
struct BoundedGelFrictionState
{
    BoundedGelState gel;
    StdVec<ContactTangentVector> last_tangent_impulse,tangent_slip_sum,tangent_slip_correction;
    StdVec<Real> ledger_sum=StdVec<Real>(size_t(ContactLedger::Count),0);
    StdVec<Real> ledger_correction=StdVec<Real>(size_t(ContactLedger::Count),0);

    Real ledger(ContactLedger index) const
    { const size_t i=size_t(index);return ledger_sum.at(i)+ledger_correction.at(i); }

    template<class Checkpoint>void registerCheckpoint(Checkpoint &cp)
    {
        gel.registerCheckpoint(cp);
        cp.addArray("LastTwoTangentImpulse",last_tangent_impulse);
        cp.addArray("TangentSlipSum",tangent_slip_sum);cp.addArray("TangentSlipCorrection",tangent_slip_correction);
        cp.addArray("ContactCumulativeLedgerSum",ledger_sum);cp.addArray("ContactCumulativeLedgerCorrection",ledger_correction);
    }
};

inline BoundedGelFrictionState initialBoundedGelFrictionState(const BoundedNormalScene &scene)
{
    BoundedGelFrictionState state;state.gel=scene.initialState();
    const size_t count=scene.plateIDs().size();
    state.gel.last_normal_impulse.assign(count,0);
    state.last_tangent_impulse.assign(count,ContactTangentVector::Zero());
    state.tangent_slip_sum=state.last_tangent_impulse;state.tangent_slip_correction=state.last_tangent_impulse;
    return state;
}

inline StdVec<CoulombContactConstraint> boundedGelFrictionContacts(const BoundedNormalScene &scene,
    Real old_time,Real new_time,const std::array<Real,2> &mu,
    const std::array<ContactTangentVector,2> &wall_tangent_velocity)
{
    const auto normal=scene.contacts(old_time,new_time);
    if (scene.config.mode==2) return {};
    StdVec<CoulombContactConstraint> rows;rows.reserve(normal.size());
    for(size_t s=0;s<scene.surface.reference.size();++s)for(size_t plate=0;plate<2;++plate)
    {
        CoulombContactConstraint row;row.weights.assign(scene.model.reference.size(),0);
        for(size_t p=scene.surface.offsets[s];p<scene.surface.offsets[s+1];++p)
            row.weights[scene.surface.targets[p]]=scene.surface.weights[p];
        row.normal=scene.normal[plate];row.tangent[0]=scene.frame*ContactSpatialVector::UnitX();
        row.tangent[1]=scene.frame*ContactSpatialVector::UnitZ();
        row.old_offset=normal.at(2*s+plate).old_offset;row.new_offset=normal.at(2*s+plate).new_offset;
        row.mu=mu[plate];row.wall_tangent_velocity=wall_tangent_velocity[plate];
        rows.push_back(std::move(row));
    }
    return rows;
}

struct BoundedGelFrictionStepResult
{
    bool accepted=false,candidate_ready=false;
    std::string status;
    BoundedGelFrictionState next;
    CoulombContactResult solve;
    BoundedGelEvaluation evaluation;
    StdVec<Vecd> contact_impulse,normalized_momentum_residual;
    std::array<Vecd,2> particle_plate_impulse{Vecd::Zero(),Vecd::Zero()};
    std::array<Vecd,2> normal_plate_impulse{Vecd::Zero(),Vecd::Zero()};
    std::array<Vecd,2> tangent_plate_impulse{Vecd::Zero(),Vecd::Zero()};
    std::array<Real,2> normal_plate_work{0,0},tangent_plate_work{0,0};
    Real kinetic_old=0,kinetic_new=0,algorithmic_loss=0,normal_contact_loss=0,friction_loss=0,energy_residual=0;
    Vecd total_momentum_residual=Vecd::Zero(),internal_force_sum=Vecd::Zero();
    Vecd angular_drift_minus_external=Vecd::Zero(),internal_old_torque=Vecd::Zero();
};

/** True local-PD material context plus simultaneous normal/two-tangent solve.
 * The old state remains untouched on EVERY failure. The caller additionally
 * applies whole-interval geometry/output gates before committing this result.
 */
inline BoundedGelFrictionStepResult advanceBoundedGelFriction(const BoundedGelModel &model,
    const BoundedGelFrictionState &old,Real dt,const StdVec<CoulombContactConstraint> &rows,
    const StdVec<size_t> &plate,const ManyContactSettings &settings={})
{
    using namespace coulomb_contact_detail;
    const BoundedGelStepContext context(model,old.gel,dt);const size_t n=model.reference.size(),m=rows.size();
    if (plate.size()!=m || old.gel.last_normal_impulse.size()!=m || old.last_tangent_impulse.size()!=m ||
        old.tangent_slip_sum.size()!=m || old.tangent_slip_correction.size()!=m ||
        old.ledger_sum.size()!=size_t(ContactLedger::Count) || old.ledger_correction.size()!=old.ledger_sum.size())
        throw std::invalid_argument("bounded gel friction: history dimensions");
    for(size_t s=0;s<m;++s)
        if (plate[s]>1 || !std::isfinite(old.gel.last_normal_impulse[s]) || !old.last_tangent_impulse[s].allFinite() ||
            !old.tangent_slip_sum[s].allFinite() || !old.tangent_slip_correction[s].allFinite())
            throw std::invalid_argument("bounded gel friction: contact history");
    for(size_t i=0;i<old.ledger_sum.size();++i)
        if (!std::isfinite(old.ledger_sum[i]) || !std::isfinite(old.ledger_correction[i]))
            throw std::invalid_argument("bounded gel friction: cumulative history");
    const CoulombContactSystem system(flattenPD(old.gel.position),dt,context.total_mass,rows,settings);
    BoundedGelFrictionStepResult r;r.next=old;
    r.solve=solveCoulombContact(system,[&](const auto &v){return context.evaluate(v).equation;},
        [&](const auto &v){return context.inDomain(v);},settings);
    if (!r.solve.converged) { r.status=r.solve.status;return r; }
    r.evaluation=context.evaluate(r.solve.velocity);r.next.gel.velocity=unflattenPD(r.solve.velocity);
    const auto &g=system.coefficients();const auto &impulse=r.solve.impulses;
    r.next.gel.time=old.gel.time+dt;r.next.gel.steps=old.gel.steps+1;
    r.next.gel.last_dt=dt;r.next.gel.last_KKT=r.solve.normalized_KKT;
    r.next.gel.last_newton_iterations=r.solve.trace.back().iteration;
    for(size_t s=0;s<m;++s)
    {
        r.next.gel.last_normal_impulse[s]=impulse[3*s];
        r.next.last_tangent_impulse[s]={impulse[3*s+1],impulse[3*s+2]};
    }
    for(size_t i=0;i<n;++i)
    {
        r.next.gel.position[i]=old.gel.position[i]+dt*r.next.gel.velocity[i];
        const auto &p=r.evaluation.point[i];
        r.next.gel.green[i]=p.history.green;r.next.gel.viscous[i]=p.history.viscous;r.next.gel.dissipation[i]=p.history.dissipation;
        const auto response=gelResponse(r.evaluation.final_deformation[i],model.material,p.history);
        r.next.gel.accepted_stress[i]=response.cauchy;r.next.gel.density[i]=model.reference_density/r.evaluation.final_deformation[i].determinant();
    }
    r.candidate_ready=true;
    Real volume=0,length=0;
    for(size_t i=0;i<n;++i){volume+=model.reference.volume[i];length=std::max(length,old.gel.position[i].norm());}
    const Real E=std::max(system.impulse_scale*system.speed,volume*std::max(model.material.bulk,model.material.equilibrium_shear));
    const Real L=std::max(length,model.reference.horizon);
    auto roundQuantity=[&](const Rational &value,Real scale)
    {
        const auto rounded=exact_contact_detail::roundRatio(numerator(value),denominator(value),0);
        if (exactReal(rounded.absolute_error_bound)/exactReal(scale)>exactReal(1e-12))
            throw std::runtime_error("bounded gel friction: ledger conversion precision budget");
        return rounded.value;
    };
    try
    {
        std::array<Rational,2> work_n{},work_t{};
        Rational normal_loss=0,friction_loss=0,k0=0,k1=0,algorithmic=0;
        const auto q0=flattenPD(old.gel.position);
        for(size_t s=0;s<m;++s)
        {
            const Rational lambda=exactReal(impulse[3*s]),actual_dt=exactReal(dt);
            work_n[plate[s]]+=lambda*(exactReal(rows[s].new_offset)-exactReal(rows[s].old_offset))/actual_dt;
            normal_loss+=lambda*(exactDot(g[3*s],q0)-exactReal(rows[s].old_offset))/actual_dt;
            for(size_t a=0;a<2;++a)
            {
                const Rational tau=exactReal(impulse[3*s+1+a]),wall=exactReal(rows[s].wall_tangent_velocity[a]);
                const Rational relative=exactDot(g[3*s+1+a],r.solve.velocity)-wall;
                work_t[plate[s]]+=tau*wall;friction_loss-=tau*relative;
                const Real increment=roundQuantity(actual_dt*relative,L);
                addCompensatedContactValue(increment,r.next.tangent_slip_sum[s][a],r.next.tangent_slip_correction[s][a]);
            }
        }
        r.contact_impulse.assign(n,Vecd::Zero());r.normalized_momentum_residual.assign(n,Vecd::Zero());
        for(size_t i=0;i<n;++i)
        {
            const Real mass=model.node_mass[i];const Vecd delta=r.next.gel.velocity[i]-old.gel.velocity[i];
            for(size_t a=0;a<3;++a)
            {
                k0+=exactReal(mass)*exactReal(old.gel.velocity[i][a])*exactReal(old.gel.velocity[i][a])/2;
                k1+=exactReal(mass)*exactReal(r.next.gel.velocity[i][a])*exactReal(r.next.gel.velocity[i][a])/2;
                // dv is the difference of original endpoint binary64 values.
                const Rational dv=exactReal(r.next.gel.velocity[i][a])-exactReal(old.gel.velocity[i][a]);
                algorithmic+=exactReal(mass)*dv*dv/2;
                StdVec<Real> c(3*m);for(size_t s=0;s<3*m;++s)c[s]=g[s][3*i+a];
                const auto scatter=exactContactBalance(0,c,impulse,system.impulse_scale);
                const auto residual=exactContactBalance(r.evaluation.equation.momentum[3*i+a],c,impulse,system.impulse_scale);
                if (!scatter.within_budget || !residual.within_budget)
                    throw std::runtime_error("bounded gel friction: impulse precision budget");
                r.contact_impulse[i][a]=roundQuantity(exactDot(c,impulse),system.impulse_scale);
                r.normalized_momentum_residual[i][a]=residual.normalized_value;
            }
            r.total_momentum_residual+=mass*delta-r.contact_impulse[i];
            r.internal_force_sum+=r.evaluation.internal_force[i];
            r.angular_drift_minus_external+=mass*(r.next.gel.position[i].cross(r.next.gel.velocity[i])-
                old.gel.position[i].cross(old.gel.velocity[i]))-old.gel.position[i].cross(r.contact_impulse[i]);
            r.internal_old_torque+=old.gel.position[i].cross(r.evaluation.internal_force[i]);
        }
        for(size_t id=0;id<2;++id)
        {
            r.normal_plate_work[id]=roundQuantity(work_n[id],E);r.tangent_plate_work[id]=roundQuantity(work_t[id],E);
            for(size_t a=0;a<3;++a)
            {
                Rational normal=0,tangent=0;
                for(size_t s=0;s<m;++s)if(plate[s]==id)for(size_t i=0;i<n;++i)
                {
                    normal+=exactReal(g[3*s][3*i+a])*exactReal(impulse[3*s]);
                    for(size_t b=1;b<3;++b)tangent+=exactReal(g[3*s+b][3*i+a])*exactReal(impulse[3*s+b]);
                }
                r.normal_plate_impulse[id][a]=roundQuantity(normal,system.impulse_scale);
                r.tangent_plate_impulse[id][a]=roundQuantity(tangent,system.impulse_scale);
                r.particle_plate_impulse[id][a]=roundQuantity(normal+tangent,system.impulse_scale);
            }
        }
        r.kinetic_old=roundQuantity(k0,E);r.kinetic_new=roundQuantity(k1,E);r.algorithmic_loss=roundQuantity(algorithmic,E);
        r.normal_contact_loss=roundQuantity(normal_loss,E);r.friction_loss=roundQuantity(friction_loss,E);
        const auto &e=r.evaluation;
        r.energy_residual=r.normal_plate_work[0]+r.normal_plate_work[1]+r.tangent_plate_work[0]+r.tangent_plate_work[1]-
            (r.kinetic_new-r.kinetic_old+e.material_stored_new-e.material_stored_old+
             e.stabilization_stored_new-e.stabilization_stored_old+e.physical_material_loss+
             r.algorithmic_loss+r.normal_contact_loss+r.friction_loss);
        if (!std::isfinite(r.energy_residual) || std::abs(r.energy_residual)>1e-10*E ||
            std::abs(e.transpose_work_residual)>1e-10*E || std::abs(e.stabilization_work_residual)>1e-10*E ||
            r.normal_contact_loss<-1e-10*E || r.friction_loss<-1e-10*E)
        { r.status="original_physical_energy_or_dissipation_gate";return r; }
        for(const auto &v:r.normalized_momentum_residual)for(size_t a=0;a<3;++a)
            if (!std::isfinite(v[a]) || std::abs(v[a])>1e-10)
            { r.status="original_coordinate_momentum_gate";return r; }
        if (r.total_momentum_residual.norm()>1e-10*system.impulse_scale ||
            r.internal_force_sum.norm()>1e-10*system.impulse_scale/dt ||
            (r.angular_drift_minus_external-dt*r.internal_old_torque).norm()>
                1e-10*system.impulse_scale*std::max(length,model.reference.horizon))
        { r.status="original_global_or_angular_cursor_ledger_gate";return r; }
        auto add=[&](size_t index,Real x)
        { addCompensatedContactValue(x,r.next.ledger_sum[index],r.next.ledger_correction[index]); };
        for(size_t id=0;id<2;++id)
        {
            for(size_t a=0;a<3;++a)add(3*id+a,r.particle_plate_impulse[id][a]);
            add(size_t(ContactLedger::LowerWork)+id,r.normal_plate_work[id]+r.tangent_plate_work[id]);
        }
        add(size_t(ContactLedger::MaxwellLoss),e.physical_material_loss);
        add(size_t(ContactLedger::NormalLoss),r.normal_contact_loss);
        add(size_t(ContactLedger::FrictionLoss),r.friction_loss);
        add(size_t(ContactLedger::AlgorithmicLoss),r.algorithmic_loss);
    }
    catch(const std::overflow_error &error) { r.status=error.what();return r; }
    catch(const std::runtime_error &error) { r.status=error.what();return r; }
    r.accepted=true;r.status="original_PD_normal_and_two_tangent_step_gates_passed";
    return r;
}
} // namespace SPH::peridynamics
#endif
