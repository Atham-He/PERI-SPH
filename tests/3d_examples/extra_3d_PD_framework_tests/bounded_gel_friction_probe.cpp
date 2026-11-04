/** Full actual-PD friction observations; no calibration or restart release. */
#include "coulomb_contact_probe_io.h"
using namespace SPH;using namespace SPH::peridynamics;
namespace io=bounded_probe_io;
namespace fio=coupled_probe_io;
namespace
{
struct Input
{
    BoundedNormalConfig scene;
    BoundedFrictionProtocol motion;
    std::array<Real,2> mu{0,0};
};
Input readInput(const char *path)
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
void header(const Input &c,const BoundedNormalScene &scene,const BoundedGelFrictionState &initial)
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
    io::field("C",normal);io::field("C_full_three_components",system.coefficients());io::field("plate_ids",scene.plateIDs());
    std::cout<<",\"initial_contacts\":";fio::contacts(rows);
    std::cout<<",\"initial_state\":";io::state(initial.gel);std::cout<<",\"initial_friction_history\":";fio::history(initial);
    io::finish();
}
void step(const BoundedNormalScene &scene,const BoundedGelFrictionState &old,const BoundedGelFrictionStepResult &r,
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
void errorStderr(const std::exception &error,Real time,uint64_t steps)
{
    std::cerr<<std::setprecision(17)<<"{\"type\":\"native_error\",\"last_durable_accepted_time\":"<<time
        <<",\"last_durable_accepted_steps\":"<<steps<<",\"message\":\"";
    for(unsigned char c:std::string(error.what()))
        if(c=='"' || c=='\\')std::cerr<<'\\'<<c;
        else if(c<32)std::cerr<<"\\u00"<<std::hex<<std::setw(2)<<std::setfill('0')<<int(c)<<std::dec<<std::setfill(' ');
        else std::cerr<<c;
    std::cerr<<"\"}\n";std::cerr.flush();
}
}
int main(int argc,char **argv)
{
    std::cout<<std::setprecision(17);std::ios::sync_with_stdio(false);
    Real durable_time=0;uint64_t durable_steps=0;
    try
    {
        if(argc==2 && std::string(argv[1])=="--capabilities")
        {
            io::begin("capabilities");io::field("format",1);io::field("actual_PD_friction_observer",true);
            io::field("piecewise_wall_protocol",true);io::field("full_restart_qualified",false);io::field("S3_complete",false);
            io::finish();return 0;
        }
        io::identity();if(argc!=2)throw std::invalid_argument("bounded friction observer: usage INPUT");
        const auto c=readInput(argv[1]);const BoundedNormalScene scene(c.scene);
        auto state=initialBoundedGelFrictionState(scene);header(c,scene,state);
        for(size_t k=0;k<c.motion.steps();++k)
        {
            const auto interval=c.motion.interval(k,state.gel.time);
            const auto rows=boundedGelFrictionContactsAtOffsets(scene,interval.old_offset,interval.new_offset,c.mu,interval.wall_tangent_velocity);
            auto r=advanceBoundedGelFriction(scene.model,state,interval.actual_dt,rows,scene.plateIDs());
            std::array<Real,2> whole{0,0},endpoint{0,0};
            if(r.accepted)
            {
                whole=boundedFrictionIntervalGaps(scene,state.gel.position,r.next.gel.position,interval);
                endpoint=boundedFrictionIntervalGaps(scene,state.gel.position,r.next.gel.position,interval,false);
            }
            const bool geometry=r.accepted && std::min(whole[0],whole[1])>=-.01*c.scene.stroke;
            const std::string reason=!r.accepted ? r.status : (geometry ? "all_applicable_step_gates_passed" : "whole_interval_dense_penetration");
            step(scene,state,r,interval,rows,geometry,whole,endpoint,reason.c_str());
            if(!geometry)
            {
                io::begin("end");io::field("completed",false);io::field("accepted_steps",durable_steps);
                io::field("last_accepted_time",durable_time);io::field("reason",reason);io::finish();return r.accepted ? 3 : 2;
            }
            // The durable validated state record is the acceptance point.
            // Scalar cursor updates cannot allocate; moving vectors applies it.
            durable_time=r.next.gel.time;durable_steps=r.next.gel.steps;
            state=std::move(r.next);
        }
        io::begin("end");io::field("completed",true);io::field("accepted_steps",durable_steps);
        io::field("last_accepted_time",durable_time);io::finish();return 0;
    }
    catch(const std::exception &error)
    {
        errorStderr(error,durable_time,durable_steps);
        try
        {
            io::begin("error");io::field("message",error.what());io::field("last_durable_accepted_time",durable_time);
            io::field("last_durable_accepted_steps",durable_steps);io::field("accepted_physical_step",false);io::finish();
        }
        catch(...) {}
        return 1;
    }
}
