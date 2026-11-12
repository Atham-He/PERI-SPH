#include <Eigen/Core>
#include "pd_coulomb_friction.h"
#include <cmath>
#include <iostream>
#include <stdexcept>

using History=SPH::peridynamics::CoulombTangentialHistory<Eigen::Vector3d>;
static void check(bool condition,const char *message){if(!condition)throw std::runtime_error(message);}
static void near(double actual,double expected,const char *message)
{check(std::abs(actual-expected)<1e-12*std::max(1.,std::abs(expected)),message);}
int main()
{
    constexpr double Kt=1000,mu=.5,Fn=10;
    auto x=[](double a){return Eigen::Vector3d(a,0,0);};History h;
    h.commit(h.trial(x(0),1,Fn,Kt,mu));
    auto sticking=h.trial(x(.002),1,Fn,Kt,mu);
    near(sticking.force[0],-2,"Elastic tangential force");near(sticking.stored,.002,"Elastic energy");
    check(h.revision==1&&h.stored==0,"Trial must not change history");
    h.commit(sticking);
    auto sliding=h.trial(x(.010),1,Fn,Kt,mu);
    near(sliding.force[0],-5,"Coulomb cap");near(sliding.stored,.0125,"Sliding elastic energy");
    near(sliding.slip_dissipation,.025,"Plastic slip dissipation");check(sliding.sliding,"Sliding flag");
    h.commit(sliding);double D=h.slip_dissipation;
    auto repeated=h.trial(x(.010),1,Fn,Kt,mu);near(repeated.slip_dissipation,0,"Repeated state slip");
    h.commit(repeated);near(h.slip_dissipation,D,"Idempotent repeated state");
    bool rejected=false;try{h.commit(repeated);}catch(const std::logic_error&){rejected=true;}
    check(rejected,"Duplicate stale commit must be rejected");
    auto unloading=h.trial(x(.009),1,Fn,Kt,mu);near(unloading.force[0],-4,"Elastic unloading");
    near(unloading.slip_dissipation,0,"No slip on elastic unloading");h.commit(unloading);
    auto reversal=h.trial(x(-.005),1,Fn,Kt,mu);near(reversal.force[0],5,"Reverse sliding");
    check(reversal.slip_dissipation>0,"Reverse slip dissipates");h.commit(reversal);
    auto normal_unload=h.trial(x(-.005),1,2,Kt,mu);near(normal_unload.force[0],1,"Changing normal cone");
    check(normal_unload.slip_dissipation>=0,"Normal unload cannot create negative slip dissipation");h.commit(normal_unload);
    auto detachment=h.trial(x(-.005),0,0,Kt,mu);near(detachment.release_energy,.0005,"Released spring energy is explicit");
    h.commit(detachment);near(h.stored,0,"Detached spring cleared");near(h.release_energy,.0005,"Separate release ledger");
    auto again=h.trial(x(10),-1,Fn,Kt,mu);near(again.force.norm(),0,"Recontact has no inherited spring offset");h.commit(again);
    History frictionless;frictionless.commit(frictionless.trial(x(0),1,Fn,Kt,0));
    near(frictionless.trial(x(1),1,Fn,Kt,0).force.norm(),0,"Zero friction");
    // A closed stick/slip cycle: external work is the exact piecewise-linear
    // tangential law integral, rather than an implementation-mirroring check.
    History cycle;cycle.commit(cycle.trial(x(0),1,Fn,Kt,mu));
    for(double point:{.010,-.010,.010})cycle.commit(cycle.trial(x(point),1,Fn,Kt,mu));
    near(cycle.slip_dissipation,.125,"Known cyclic plastic work");
    near(cycle.slip_dissipation+cycle.stored,.1375,"Known cyclic external work");
    std::cout<<"{\"status\":\"COULOMB_STICK_SLIP_HISTORY_CHECKS_PASSED\",\"cycle_dissipation_J\":"
             <<cycle.slip_dissipation<<",\"sliding_cap_N\":5,\"separate_release_energy_J\":"<<h.release_energy<<"}\n";
}
