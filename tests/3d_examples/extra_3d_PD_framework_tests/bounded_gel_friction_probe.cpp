/** Full actual-PD friction observations; no calibration or restart release. */
#include "bounded_gel_friction_observer.h"
using namespace SPH;using namespace SPH::peridynamics;
namespace io=bounded_probe_io;
namespace fio=coupled_probe_io;
namespace
{
using namespace bounded_friction_observer;
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
        auto state=initialBoundedGelFrictionState(scene);header(c,scene,state,scene.plateIDs());
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
