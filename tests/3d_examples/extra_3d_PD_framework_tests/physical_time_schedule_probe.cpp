#include "pd_physical_time_schedule.h"
#include <charconv>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <iterator>
#include <string>

using namespace SPH::peridynamics;
uint64_t number(const char *s)
{
    uint64_t x=0;const auto end=s+std::strlen(s);auto p=std::from_chars(s,end,x);
    if(p.ec!=std::errc{} || p.ptr!=end)throw std::invalid_argument("probe: unsigned integer");
    return x;
}
void sample(const PhysicalPlateSample &s)
{
    std::cout<<'['<<s.normal_offset<<','<<s.normal_velocity;
    for(auto x:s.tangent_displacement)std::cout<<','<<x;
    for(auto x:s.tangent_velocity)std::cout<<','<<x;
    std::cout<<']';
}
std::vector<uint8_t> read(const char *path)
{
    std::ifstream f(path,std::ios::binary|std::ios::ate);if(!f)throw std::invalid_argument("probe: packet open");
    const auto size=f.tellg();
    if(size<0 || size>1024)throw std::invalid_argument("probe: packet size bound");
    f.seekg(0);
    std::vector<uint8_t> b(std::istreambuf_iterator<char>(f),{});
    return b;
}
int main(int argc,char **argv)
{
    std::cout<<std::setprecision(17);const PhysicalTimeSchedule schedule;
    try
    {
        if(argc<2)throw std::invalid_argument("probe: missing mode");
        const std::string mode=argv[1];
        if(mode=="sample" && argc==3)
        {
            const auto t=number(argv[2]);std::cout<<"{\"tick\":"<<t<<",\"time\":"<<schedule.time(t)<<",\"wall\":";
            sample(schedule.sample(t));std::cout<<"}\n";return 0;
        }
        if(mode=="restore" && argc==4)
        {
            const auto step=number(argv[2]);uint64_t live=12345;const auto before=schedule.packet(step,live);
            try
            {
                const auto staged=PhysicalTimeSchedule::restoreCursor(read(argv[3]),schedule,step);
                live=staged;
                std::cout<<"{\"restored\":true,\"tick\":"<<live<<",\"new_physical_steps\":0}\n";return 0;
            }
            catch(const std::exception &e)
            {
                if(schedule.packet(step,live)!=before)throw std::runtime_error("probe: rejected restore changed live packet");
                std::cerr<<e.what()<<'\n';
                std::cout<<"{\"restored\":false,\"live_unchanged\":true,\"new_physical_steps\":0}\n";return 1;
            }
        }
        if(mode=="grid" && argc==4)
        {
            const auto step=number(argv[2]);auto tick=number(argv[3]);uint64_t count=0;
            do
            {
                if(++count>100000)throw std::invalid_argument("probe: interval output bound");
                const auto v=schedule.interval(tick,step);
                std::cout<<"{\"old_tick\":"<<v.old_tick<<",\"new_tick\":"<<v.new_tick
                    <<",\"phase\":"<<v.normal_phase<<",\"old_time\":"<<v.old_time<<",\"new_time\":"<<v.new_time
                    <<",\"actual_dt\":"<<v.actual_dt<<",\"old_wall\":";sample(v.old_wall);
                std::cout<<",\"new_wall\":";sample(v.new_wall);std::cout<<",\"normal_cubic\":[";
                for(size_t k=0;k<4;++k){if(k)std::cout<<',';std::cout<<v.normal_cubic[k];}
                std::cout<<"],\"tangent_interval_velocity\":["<<v.tangent_interval_velocity[0]<<','<<v.tangent_interval_velocity[1]<<"]}\n";
                tick=v.new_tick;
            }while(tick<schedule.endTick());
            return 0;
        }
        throw std::invalid_argument("probe: unsupported mode/argument count");
    }
    catch(const std::exception &e){std::cerr<<e.what()<<'\n';return 1;}
}
