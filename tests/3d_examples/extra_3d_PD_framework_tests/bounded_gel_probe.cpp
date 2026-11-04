/** Full accepted-step actualPD normal sphere observation, not calibration. */
#include "bounded_gel_probe_io.h"
using namespace SPH;using namespace SPH::peridynamics;
namespace io=bounded_probe_io;
namespace
{
BoundedNormalConfig readInput(const char *filename)
{
    std::ifstream input(filename);BoundedNormalConfig c;
    input>>c.mode>>c.steps>>c.dt>>c.stroke>>c.rotation>>c.azimuths>>c.rho0>>c.material.bulk>>c.material.equilibrium_shear
        >>c.material.branch_shear>>c.material.relaxation_time>>c.eta>>c.support>>c.cap_width>>c.relaxed_fraction>>c.perturb_node>>c.invalid_injection;
    for(size_t i=0;i<9;++i)input>>c.initial_affine(i/3,i%3);
    for(size_t i=0;i<9;++i)input>>c.velocity_rate(i/3,i%3);
    for(auto *v:{&c.local_perturbation,&c.local_velocity_perturbation,&c.translation})for(size_t a=0;a<3;++a)input>>(*v)[a];
    if(!input)throw std::invalid_argument("truncated bounded sphere input");input>>std::ws;if(!input.eof())throw std::invalid_argument("bounded sphere input tail");
    c.validate();return c;
}
void header(const BoundedNormalScene &scene,const BoundedGelState &initial)
{
    const auto &c=scene.config;io::begin("header");io::field("format",1);io::field("actual_volume_nodes",scene.model.reference.size());
    io::field("normalized_benchmark",true);io::field("time_units","normalized_model_time");io::field("normal_only",true);io::field("full_S3_complete",false);
    io::field("mode",c.mode);io::field("steps_requested",c.steps);io::field("nominal_dt",c.dt);io::field("end_time",c.endTime());io::field("stroke",c.stroke);
    io::field("frame",scene.frame);io::field("parameters",StdVec<Real>{c.material.bulk,c.material.equilibrium_shear,c.material.branch_shear,c.material.relaxation_time});
    io::field("rho0",c.rho0);io::field("eta",c.eta);io::field("physical_mass",scene.model.node_mass);io::field("solver_velocity_scale",1);
    io::field("plate_normal",scene.normal);io::field("dense_penetration_stroke_gate",.01);
    io::field("geometry_cursor","minimum_over_entire_linear_particle_path_and_true_cubic_wall_interval");
    io::field("area_union_is_integral",false);std::cout<<",\"reference\":";io::reference(scene.model.reference);
    std::cout<<",\"surface\":";io::surface(scene.surface);std::cout<<",\"dense_surface\":";io::surface(scene.dense_surface);
    const auto rows=scene.contacts(0,0);DenseContactMatrix coefficients;for(const auto &row:rows)coefficients.push_back(row.coefficients);
    io::field("C",coefficients);io::field("plate_ids",scene.plateIDs());std::cout<<",\"initial_state\":";io::state(initial);io::finish();
}
void stepRecord(const BoundedNormalScene &scene,const BoundedGelState &old,const BoundedGelStepResult &r,
    const StdVec<DenseContactConstraint> &rows,Real dt,bool commit,const std::array<Real,2> &gaps,const std::array<Real,2> &endpoint_gaps,const char *reason)
{
    io::begin(commit?"step":"rejected_step");io::field("committed",commit);io::field("local_root_accepted",r.accepted);io::field("reason",reason);
    io::field("old_time",old.time);io::field("actual_dt",dt);io::field("candidate_time",old.time+dt);
    io::field("upper_old",scene.config.upperDisplacement(old.time));io::field("upper_new",scene.config.upperDisplacement(old.time+dt));
    StdVec<Real> old_offsets,new_offsets;for(const auto &row:rows){old_offsets.push_back(row.old_offset);new_offsets.push_back(row.new_offset);}
    io::field("old_offsets",old_offsets);io::field("new_offsets",new_offsets);std::cout<<",\"solve\":";io::solver(r.solve);
    if(r.accepted)
    {
        std::cout<<",\"state\":";io::state(r.next);std::cout<<",\"evaluation\":";io::evaluation(r.evaluation,scene.model,old);
        io::field("contact_impulse",r.contact_impulse);io::field("normalized_momentum_residual",r.momentum_residual);io::field("particle_plate_impulse",r.particle_plate_impulse);
        std::array<Vecd,2> reaction{-r.particle_plate_impulse[0]/dt,-r.particle_plate_impulse[1]/dt};io::field("plate_reaction_interval",reaction);
        io::field("force_cursor","impulse_over_actual_interval_dt_not_instantaneous_terminal_force");io::field("plate_work",r.plate_work);
        io::field("kinetic_old",r.kinetic_old);io::field("kinetic_new",r.kinetic_new);io::field("algorithmic_loss",r.algorithmic_loss);
        io::field("normal_contact_loss",r.normal_contact_loss);io::field("friction_implemented",false);io::field("energy_residual",r.energy_residual);
        io::field("total_momentum_residual",r.total_momentum_residual);io::field("internal_force_sum",r.internal_force_sum);
        io::field("angular_drift_minus_external",r.angular_drift_minus_external);io::field("internal_old_torque",r.internal_old_torque);
        io::field("geometry_gate_applicable",scene.config.mode==1);
        if(scene.config.mode==1){io::field("dense_minimum_plate_gap",gaps);io::field("dense_endpoint_minimum_plate_gap",endpoint_gaps);}
        else std::cout<<",\"dense_minimum_plate_gap\":null,\"dense_endpoint_minimum_plate_gap\":null";
    }
    io::finish();
}
}
int main(int argc,char **argv)
{
    std::cout<<std::setprecision(17);std::ios::sync_with_stdio(false);
    try
    {
        if(argc==2 && std::string(argv[1])=="--capabilities")
        {
            io::begin("capabilities");io::field("format",1);io::field("actual_PD_normal_observer",true);io::field("full_restart_qualified",false);io::field("full_S3_complete",false);io::finish();return 0;
        }
        io::identity();if(argc!=2)throw std::invalid_argument("usage: bounded_gel_probe INPUT");
        const auto c=readInput(argv[1]);const BoundedNormalScene scene(c);auto state=scene.initialState();header(scene,state);
        if(c.mode==0)
        {
            const BoundedGelStepContext context(scene.model,state,c.dt);const auto v=scene.sampleVelocity();
            const auto e=context.evaluate(v);io::begin("operator_sample");io::field("accepted_physical_step",false);io::field("dt",c.dt);io::field("velocity",v);
            long double kinetic=0;for(size_t i=0;i<scene.model.node_mass.size();++i)for(size_t a=0;a<3;++a)
                kinetic+=.5L*scene.model.node_mass[i]*v[3*i+a]*v[3*i+a];
            io::field("trial_kinetic_actual_mass",Real(kinetic));
            std::cout<<",\"evaluation\":";io::evaluation(e,scene.model,state);io::finish();
            io::begin("end");io::field("completed",true);io::field("accepted_steps",0);io::field("last_accepted_time",state.time);io::finish();return 0;
        }
        for(size_t k=0;k<c.steps;++k)
        {
            const Real dt=k+1==c.steps?c.endTime()-state.time:c.dt;
            const auto rows=scene.contacts(state.time,state.time+dt);const auto r=advanceBoundedGel(scene.model,state,dt,rows,scene.plateIDs());
            if(!r.accepted)
            {
                stepRecord(scene,state,r,rows,dt,false,{0,0},{0,0},r.solve.status.c_str());io::begin("end");io::field("completed",false);
                io::field("accepted_steps",state.steps);io::field("last_accepted_time",state.time);io::field("reason",r.solve.status);io::finish();return 2;
            }
            const auto endpoint_gaps=c.mode==1?scene.denseGaps(r.next.position,r.next.time):std::array<Real,2>{0,0};
            const auto gaps=c.mode==1?scene.denseIntervalGaps(state.position,r.next.position,state.time,r.next.time):std::array<Real,2>{0,0};
            const bool geometry=c.mode!=1 || std::min(gaps[0],gaps[1])>=-.01*c.stroke;
            stepRecord(scene,state,r,rows,dt,geometry,gaps,endpoint_gaps,geometry?"all_applicable_native_gates_passed":"dense_surface_penetration");
            if(!geometry)
            {
                io::begin("end");io::field("completed",false);io::field("accepted_steps",state.steps);io::field("last_accepted_time",state.time);
                io::field("reason","dense_surface_penetration");io::finish();return 3;
            }
            state=r.next;
        }
        io::begin("end");io::field("completed",true);io::field("accepted_steps",state.steps);io::field("last_accepted_time",state.time);io::finish();return 0;
    }
    catch(const std::exception &error)
    {
        try{io::begin("error");io::field("message",error.what());io::field("physical_release",false);io::finish();}catch(...){}
        return 1;
    }
}
