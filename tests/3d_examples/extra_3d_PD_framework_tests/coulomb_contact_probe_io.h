/* Full original observations; no material or physical state mutation. */
#ifndef COULOMB_CONTACT_PROBE_IO_H
#define COULOMB_CONTACT_PROBE_IO_H
#include "bounded_gel_probe_io.h"
#include "pd_bounded_friction_protocol.h"
namespace coupled_probe_io
{
using namespace SPH;using namespace SPH::peridynamics;
namespace io=bounded_probe_io;
inline DenseContactMatrix matrix2(const ContactTangentMatrix &m)
{ return {{m(0,0),m(0,1)},{m(1,0),m(1,1)}}; }
inline void disk(const CoulombDiskValue &r)
{
    std::cout<<"{\"branch\":";
    io::value(r.branch==CoulombDiskBranch::StickInterior ? "stick_interior" :
        (r.branch==CoulombDiskBranch::SlipExterior ? "slip_exterior" : "zero_radius_origin"));
    io::field("representable",r.representable);io::field("within_budget",r.within_budget);
    io::field("residual",r.residual);io::field("normal_derivative",r.normal_derivative);
    io::field("impulse_derivative",matrix2(r.impulse_derivative));io::field("velocity_derivative",matrix2(r.velocity_derivative));
    io::field("residual_roundoff_bound",r.residual_roundoff_bound);
    io::field("maximum_derivative_roundoff_bound",r.maximum_derivative_roundoff_bound);
    io::field("friction_loss",r.friction_loss);io::field("friction_loss_roundoff_bound",r.friction_loss_roundoff_bound);
    io::field("cone_excess_upper",r.cone_excess_upper);io::field("maximum_dissipation_error_upper",r.maximum_dissipation_error_upper);
    std::cout<<'}';
}
inline void contact(const CoulombContactConstraint &c)
{
    std::cout<<"{\"weights\":";io::value(c.weights);io::field("normal",c.normal);io::field("tangent",c.tangent);
    io::field("old_offset",c.old_offset);io::field("new_offset",c.new_offset);io::field("mu",c.mu);
    io::field("wall_tangent_velocity",c.wall_tangent_velocity);std::cout<<'}';
}
inline void contacts(const StdVec<CoulombContactConstraint> &rows)
{
    std::cout<<'[';for(size_t i=0;i<rows.size();++i){if(i)std::cout<<',';contact(rows[i]);}std::cout<<']';
}
inline void values(const CoulombContactValues &r)
{
    std::cout<<"{\"residual\":";io::value(r.residual);io::field("residual_roundoff_bound",r.residual_roundoff_bound);
    io::field("supplied_momentum",r.supplied_momentum);io::field("supplied_jacobian",r.supplied_jacobian);
    io::field("normalized_gaps",r.normalized_gaps);io::field("gap_roundoff_bound",r.gap_roundoff_bound);
    io::field("normalized_relative_velocity",r.normalized_relative_velocity);
    io::field("relative_velocity_roundoff_bound",r.relative_velocity_roundoff_bound);
    io::field("jacobian",r.jacobian);io::field("merit",r.merit);io::field("squared_norm",r.squared_norm);io::field("KKT",r.KKT);
    io::field("maximum_jacobian_roundoff_bound",r.maximum_jacobian_roundoff_bound);
    io::field("precision_budget_ok",r.precision_budget_ok);
    std::cout<<",\"disk\":[";for(size_t s=0;s<r.disk.size();++s){if(s)std::cout<<',';disk(r.disk[s]);}std::cout<<"]}";
}
inline void solve(const CoulombContactResult &r)
{
    std::cout<<"{\"converged\":";io::value(r.converged);io::field("status",r.status);
    io::field("velocity",r.velocity);io::field("impulses",r.impulses);io::field("normalized_KKT",r.normalized_KKT);
    std::cout<<",\"trace\":[";
    for(size_t i=0;i<r.trace.size();++i)
    {
        if(i)std::cout<<',';const auto &t=r.trace[i];
        std::cout<<"{\"iteration\":";io::value(t.iteration);io::field("linear_rank",t.linear_rank);
        io::field("domain_rejections",t.domain_rejections);io::field("alpha",t.alpha);io::field("linear_residual",t.linear_residual);
        io::field("velocity",t.velocity);io::field("impulses",t.impulses);io::field("normalized_direction",t.normalized_direction);
        std::cout<<",\"original\":";values(t.original);std::cout<<'}';
    }
    std::cout<<"],\"precision_rejection\":{\"present\":";io::value(r.precision_rejection.present);
    if(r.precision_rejection.present)
    {
        const auto &p=r.precision_rejection;io::field("during_line_search",p.during_line_search);io::field("iteration",p.iteration);
        io::field("alpha",p.alpha);io::field("velocity",p.velocity);io::field("impulses",p.impulses);
        std::cout<<",\"original\":";values(p.original);
    }
    std::cout<<"},\"evaluation_failure\":{\"present\":";io::value(r.evaluation_failure.present);
    if(r.evaluation_failure.present)
    {
        const auto &p=r.evaluation_failure;io::field("during_line_search",p.during_line_search);io::field("iteration",p.iteration);
        io::field("alpha",p.alpha);io::field("velocity",p.velocity);io::field("impulses",p.impulses);io::field("reason",p.reason);
    }
    std::cout<<"}}";
}
inline void history(const BoundedGelFrictionState &s)
{
    std::cout<<"{\"last_tangent_impulse\":";io::value(s.last_tangent_impulse);
    io::field("relative_tangent_travel_sum",s.tangent_slip_sum);io::field("relative_tangent_travel_correction",s.tangent_slip_correction);
    io::field("relative_travel_includes_open_rows",true);
    io::field("ledger_sum",s.ledger_sum);io::field("ledger_correction",s.ledger_correction);
    io::field("ledger_order","lower_Ixyz,upper_Ixyz,lower_work,upper_work,Maxwell,normal,friction,algorithmic");std::cout<<'}';
}
inline void interval(const BoundedFrictionInterval &i)
{
    std::cout<<"{\"phase\":";io::value(i.phase);io::field("step",i.step);
    io::field("old_time",i.old_time);io::field("new_time",i.new_time);io::field("actual_dt",i.actual_dt);
    io::field("old_offset",i.old_offset);io::field("new_offset",i.new_offset);
    io::field("normal_cubic",std::array<std::array<Real,4>,2>{i.normal_cubic[0].coefficients,i.normal_cubic[1].coefficients});
    io::field("old_tangent_displacement",i.old_tangent_displacement);io::field("new_tangent_displacement",i.new_tangent_displacement);
    io::field("wall_tangent_velocity",i.wall_tangent_velocity);
    io::field("wall_velocity_cursor","actual_interval_mean_prescribed_displacement_difference_over_actual_dt");std::cout<<'}';
}
} // namespace coupled_probe_io
#endif
