/** Prescribed-motion arithmetic only; zero physical accepted steps. */
#include "../src/shared/peridynamics/pd_early_tangent_cycle.h"
#include <iomanip>
#include <iostream>
#include <limits>
int main()
{
    using SPH::peridynamics::EarlyTangentCycle;
    EarlyTangentCycle cycle;
    int refused=0;
    for(double dt:{0.,-1.,std::numeric_limits<double>::infinity(),std::numeric_limits<double>::max()})
        try{EarlyTangentCycle{dt}.validate();return 2;}catch(const std::invalid_argument&){++refused;}
    for(double t:{-1.,std::numeric_limits<double>::quiet_NaN(),std::numeric_limits<double>::infinity()})
        try{cycle.at(t);return 3;}catch(const std::invalid_argument&){++refused;}
    std::cout<<std::setprecision(17)<<"{\"physical_steps\":0,\"refused\":"<<refused<<",\"samples\":[";
    bool first=true;
    for(double dt:{.01,.003125})
        for(int i=0;i<=1000;++i)
        {
            const double time=double(i)/1000*(40*dt);
            if(!first)std::cout<<',';first=false;
            std::cout<<'['<<dt<<','<<time<<','<<EarlyTangentCycle{dt}.at(time)<<']';
        }
    std::cout<<"]}\n";
}
