/* Full original friction observer, reusable without constructing a fresh scene. */
#ifndef BOUNDED_GEL_FRICTION_OBSERVER_H
#define BOUNDED_GEL_FRICTION_OBSERVER_H
#include "coulomb_contact_probe_io.h"
namespace bounded_friction_observer
{
using namespace SPH;using namespace SPH::peridynamics;
namespace io=bounded_probe_io;namespace fio=coupled_probe_io;
struct Input
{
    BoundedNormalConfig scene;
    BoundedFrictionProtocol motion;
    std::array<Real,2> mu{0,0};
};
inline Input readInput(const char *path)
{
    std::ifstream input(path);Input c;c.scene.azimuths=8;
    input>>c.motion.kind>>c.scene.rotation>>c.scene.azimuths>>c.scene.dt>>c.scene.stroke
        >>c.motion.tangent_stroke[0]>>c.motion.tangent_stroke[1]>>c.mu[0]>>c.mu[1]>>c.scene.rho0
        >>c.scene.material.bulk>>c.scene.material.equilibrium_shear>>c.scene.material.branch_shear
        >>c.scene.material.relaxation_time>>c.scene.eta>>c.scene.support>>c.scene.cap_width;
    if(!input)throw std::invalid_argument("bounded friction observer: expected17input fields");
    input>>std::ws;if(!input.eof())throw std::invalid_argument("bounded friction observer: input tail");
    c.scene.steps=16;c.scene.mode=1;c.scene.validate();
    c.motion.nominal_dt=c.scene.dt;c.motion.normal_stroke=c.scene.stroke;c.motion.validate();
    for(Real mu:c.mu)if(!std::isfinite(mu) || mu<0)throw std::invalid_argument("bounded friction observer: coefficient");
    return c;
}
template<class Scene> inline void header(const Input &c,const Scene &scene,const BoundedGelFrictionState &initial,
    const StdVec<size_t> &plate_ids)
{
    io::begin("header");io::field("format",1);io::field("actual_PD_body",true);io::field("actual_volume_nodes",33);
    io::field("normalized_benchmark",true);io::field("normal_only",false);io::field("full_S3_complete",false);
    io::field("mode",1);io::field("steps_requested",c.motion.steps());io::field("nominal_dt",c.motion.nominal_dt);
    io::field("end_time",c.motion.endTime());io::field("stroke",c.scene.stroke);io::field("frame",scene.frame);
    io::field("parameters",StdVec<Real>{c.scene.material.bulk,c.scene.material.equilibrium_shear,c.scene.material.branch_shear,c.scene.material.relaxation_time});
    io::field("rho0",c.scene.rho0);io::field("eta",c.scene.eta);io::field("physical_mass",scene.model.node_mass);
    io::field("solver_velocity_scale",1);io::field("plate_normal",scene.normal);
    io::field("plate_tangent",std::array<Vecd,2>{scene.frame*Vecd::UnitX(),scene.frame*Vecd::UnitZ()});
    io::field("mu",c.mu);io::field("tangent_stroke",c.motion.tangent_stroke);
    io::field("protocol_kind",c.motion.kind);io::field("phase_end_times",c.motion.phaseEnds());
    io::field("phase_step_counts",StdVec<size_t>{16,8,16});
    if(c.motion.kind==2)
    {
        io::field("tangent_cycle_end_times",EarlyTangentCycle{c.motion.nominal_dt}.ends());
        io::field("tangent_cycle_factors",StdVec<Real>{0,1,-1,0});
        io::field("tangent_cycle_cursor","independent_C1_early_cycle_then_zero");
    }
    io::field("wall_velocity_cursor","actual_interval_mean_prescribed_displacement_difference_over_actual_dt");
    io::field("dense_penetration_stroke_gate",.01);
    io::field("geometry_cursor","minimum_over_each_linear_particle_and_true_prescribed_cubic_wall_interval");
    io::field("area_union_is_integral",false);
    std::cout<<",\"reference\":";io::reference(scene.model.reference);
    std::cout<<",\"surface\":";io::surface(scene.surface);std::cout<<",\"dense_surface\":";io::surface(scene.dense_surface);
    const auto interval=c.motion.interval(0,0);
    const auto rows=boundedGelFrictionContactsAtOffsets(scene,interval.old_offset,interval.new_offset,c.mu,interval.wall_tangent_velocity);
    Real mass=0;for(Real x:scene.model.node_mass)mass+=x;
    const CoulombContactSystem system(flattenPD(initial.gel.position),interval.actual_dt,mass,rows);
    io::field("normalization_mass_scale",mass);io::field("normalization_impulse_scale",system.impulse_scale);
    io::field("primal_dofs",system.dofs);io::field("contact_count",system.contacts_count);io::field("unknowns",system.unknowns);
    DenseContactMatrix normal;for(size_t s=0;s<rows.size();++s)normal.push_back(system.coefficients()[3*s]);
    io::field("C",normal);io::field("C_full_three_components",system.coefficients());io::field("plate_ids",plate_ids);
    std::cout<<",\"initial_contacts\":";fio::contacts(rows);
    std::cout<<",\"initial_state\":";io::state(initial.gel);std::cout<<",\"initial_friction_history\":";fio::history(initial);
    io::finish();
}
template<class Scene> inline void step(const Scene &scene,const BoundedGelFrictionState &old,const BoundedGelFrictionStepResult &r,
    const BoundedFrictionInterval &interval,const StdVec<CoulombContactConstraint> &rows,
    bool committed,const std::array<Real,2> &whole,const std::array<Real,2> &endpoint,const char *reason)
{
    io::begin(committed ? "step" : "rejected_step");io::field("committed",committed);io::field("local_step_gates_passed",r.accepted);
    io::field("candidate_ready",r.candidate_ready);io::field("reason",reason);
    std::cout<<",\"interval\":";fio::interval(interval);std::cout<<",\"contacts\":";fio::contacts(rows);
    std::cout<<",\"solve\":";fio::solve(r.solve);
    if(r.candidate_ready)
    {
        std::cout<<",\"state\":";io::state(r.next.gel);std::cout<<",\"friction_history\":";fio::history(r.next);
        std::cout<<",\"evaluation\":";io::evaluation(r.evaluation,scene.model,old.gel);
        io::field("contact_impulse",r.contact_impulse);io::field("normalized_momentum_residual",r.normalized_momentum_residual);
        io::field("particle_plate_impulse",r.particle_plate_impulse);io::field("normal_plate_impulse",r.normal_plate_impulse);
        io::field("tangent_plate_impulse",r.tangent_plate_impulse);
        const std::array<Vecd,2> reaction{-r.particle_plate_impulse[0]/interval.actual_dt,-r.particle_plate_impulse[1]/interval.actual_dt};
        io::field("plate_reaction_interval",reaction);io::field("force_cursor","impulse_over_actual_interval_dt_not_instantaneous_terminal_force");
        io::field("normal_plate_work",r.normal_plate_work);io::field("tangent_plate_work",r.tangent_plate_work);
        io::field("kinetic_old",r.kinetic_old);io::field("kinetic_new",r.kinetic_new);io::field("algorithmic_loss",r.algorithmic_loss);
        io::field("normal_contact_loss",r.normal_contact_loss);io::field("friction_loss",r.friction_loss);io::field("energy_residual",r.energy_residual);
        io::field("total_momentum_residual",r.total_momentum_residual);io::field("internal_force_sum",r.internal_force_sum);
        io::field("angular_drift_minus_external",r.angular_drift_minus_external);io::field("internal_old_torque",r.internal_old_torque);
        io::field("geometry_gate_applicable",r.accepted);
        if(r.accepted){io::field("dense_minimum_plate_gap",whole);io::field("dense_endpoint_minimum_plate_gap",endpoint);}
    }
    io::finish();
}
}
#endif
