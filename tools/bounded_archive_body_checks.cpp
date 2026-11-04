/** Existing nonzero-history data codec only: no scene construction or solve. */
#include "../src/shared/peridynamics/pd_bounded_archive_body.h"
#include <fcntl.h>
#include <fstream>
#include <functional>
#include <iomanip>
#include <iostream>
#include <sstream>
using namespace SPH;
using namespace SPH::peridynamics;
using namespace SPH::peridynamics::canonical_archive;
namespace
{
void require(bool ok,const char *why){if(!ok)throw std::runtime_error(why);}
std::string quote(const std::string &s)
{
    std::ostringstream out;out<<'"';
    for(unsigned char c:s)
        if(c=='"'||c=='\\')out<<'\\'<<c;
        else if(c<32)out<<"\\u00"<<std::hex<<std::setw(2)<<std::setfill('0')<<unsigned(c)<<std::dec;
        else out<<c;
    out<<'"';return out.str();
}
Bytes load(const std::string &path)
{
    std::ifstream stream(path,std::ios::binary|std::ios::ate);require(bool(stream),"archive input open");
    const auto length=stream.tellg();require(length>=0 && uint64_t(length)<=maximum_bytes,"archive input bound");
    Bytes bytes(size_t(length),0);stream.seekg(0);stream.read(reinterpret_cast<char*>(bytes.data()),length);
    require(bool(stream),"archive input read");return bytes;
}
void store(const std::string &path,const Bytes &bytes)
{
    Descriptor target(::open(path.c_str(),O_WRONLY|O_CREAT|O_EXCL|O_CLOEXEC,0600));size_t offset=0;
    while(offset<bytes.size())
    {
        const auto n=::write(target.get(),bytes.data()+offset,bytes.size()-offset);
        if(n<0 && errno==EINTR)continue;
        require(n>0,"archive evidence write");offset+=size_t(n);
    }
    require(::fsync(target.get())==0,"archive evidence fsync");
}
}
int main(int argc,char **argv)
{
    try
    {
        require(argc==3,"usage INPUT OUTPUT_DIR");const std::string directory=argv[2];
        std::ifstream stat_file("/proc/self/stat"),boot_file("/proc/sys/kernel/random/boot_id");
        std::string stat_line,boot;std::getline(stat_file,stat_line);std::getline(boot_file,boot);
        const auto closing=stat_line.rfind(')');require(closing!=std::string::npos,"native identity stat");
        std::istringstream stat_fields(stat_line.substr(closing+2));std::string token,start_ticks;
        for(size_t i=0;i<=19;++i){require(bool(stat_fields>>token),"native identity fields");if(i==19)start_ticks=token;}
        require(!boot.empty(),"native identity boot");const auto pid=::getpid(),ppid=::getppid();
        BoundedArchiveBody empty;
        require(empty.model.reference.position.empty() && empty.model.reference.shape_gradient.empty() &&
            empty.surface.reference.empty() && empty.dense_surface.reference.empty() && empty.accepted.gel.position.empty(),"receiver not empty");
        const auto input=load(argv[1]);const uint64_t generation=37;
        const auto baseline=decodeBoundedArchiveBody(input,generation);
        require(baseline.accepted.gel.steps==40 && baseline.accepted.gel.time==.4,"nonzero-history fixture cursor");
        const auto encoded=encodeBoundedArchiveBody(baseline,generation);require(encoded==input,"full body bytes differ");
        store(directory+"/roundtrip.bin",encoded);
        std::vector<std::pair<std::string,std::string>> refusals;
        auto refuse=[&](const std::string &label,const std::string &reason,const std::function<void(BoundedArchiveBody&)> &damage)
        {
            auto candidate=baseline;damage(candidate);const auto bytes=encodeBoundedArchiveBody(candidate,generation);
            store(directory+"/"+label+".bin",bytes);
            try{auto restored=decodeBoundedArchiveBody(bytes,generation);(void)restored;}
            catch(const std::exception &error)
            {
                require(std::string(error.what()).find(reason)!=std::string::npos,"wrong semantic refusal reason");
                refusals.emplace_back(label,error.what());
                require(encodeBoundedArchiveBody(baseline,generation)==input,"baseline changed on failed decode");return;
            }
            throw std::runtime_error("invalid semantic candidate accepted: "+label);
        };
        refuse("active_mask","invalid directed packet",[](auto &b){b.model.reference.active[0]=0;});
        refuse("kernel_weight","stored kernel mismatch",[](auto &b){b.model.reference.weight[0]*=1.1;});
        refuse("shape_inverse","stored inverse/measure mismatch",[](auto &b){b.model.reference.inverse_shape[0](0,0)+=.01;});
        refuse("shape_gradient_null_affine","stored gradient mismatch",[](auto &b)
        {
            const std::array<Real,3> coordinate{-.5,0,.5},change{.001,-.002,.001};
            for(size_t k=0;k<3;++k)
            {
                const Vecd expected=b.frame*Vecd(coordinate[k],0,0);size_t found=0;
                for(;found<b.model.reference.size();++found)if((b.model.reference.position[found]-expected).norm()==0)break;
                require(found<b.model.reference.size(),"fixture collinear reference missing");
                b.model.reference.shape_gradient[found][0]+=change[k];
            }
        });
        refuse("nonaffine_cache","stored nonaffine cache mismatch",[](auto &b){b.model.reference.nonaffine_unit[0]+=.01;});
        refuse("signed_MLS","stored signed surface weight mismatch",[](auto &b)
        {auto it=std::min_element(b.surface.weights.begin(),b.surface.weights.end());require(*it< -1e-8,"signed fixture absent");*it=0;});
        refuse("physical_mass","physical mass mismatch",[](auto &b){b.model.node_mass[0]*=2;});
        refuse("material_green","accepted material cursor",[](auto &b){b.accepted.gel.green[0](0,0)+=.01;b.accepted.gel.green[0](1,1)-=.01;});
        refuse("material_viscous","invalid gel material history",[](auto &b){b.accepted.gel.viscous[0](0,1)+=.01;});
        refuse("material_dissipation","invalid gel dissipation history",[](auto &b){b.accepted.gel.dissipation[0]=-.01;});
        refuse("accepted_stress","accepted material cursor",[](auto &b){b.accepted.gel.accepted_stress[0](0,0)+=.01;});
        refuse("accepted_density","accepted material cursor",[](auto &b){b.accepted.gel.density[0]*=1.01;});
        refuse("normal_dual","stored dual cone",[](auto &b){b.accepted.gel.last_normal_impulse[0]=-.01;});
        refuse("tangent_dual","stored dual cone",[](auto &b){b.accepted.last_tangent_impulse[0][0]+=.01;});
        refuse("Maxwell_compensation","Maxwell cumulative history",[](auto &b){b.accepted.ledger_correction[size_t(ContactLedger::MaxwellLoss)]+=.01;});
        refuse("contact_order","contact order",[](auto &b){b.plate_ids[0]=1;});
        refuse("contact_frame","stored contact frame",[](auto &b){b.tangent[0][0]+=.01;});
        Builder extra(generation);FieldWriter fields(extra);auto copy=baseline;copy.registerCheckpoint(fields);extra.integers("unknown_field",{},{1});
        const auto extra_bytes=extra.finish();store(directory+"/unknown_schema.bin",extra_bytes);
        try{auto candidate=decodeBoundedArchiveBody(extra_bytes,generation);(void)candidate;throw std::runtime_error("extra schema accepted");}
        catch(const std::invalid_argument &e){require(std::string(e.what()).find("exact schema names")!=std::string::npos,"wrong schema refusal");refusals.emplace_back("unknown_schema",e.what());}
        require(encodeBoundedArchiveBody(baseline,generation)==input,"baseline bytes changed");
        size_t nonzero_viscous=0,nonzero_tangent=0,nonzero_travel_correction=0;
        for(const auto &x:baseline.accepted.gel.viscous)nonzero_viscous+=x.squaredNorm()!=0;
        for(const auto &x:baseline.accepted.last_tangent_impulse)nonzero_tangent+=x.squaredNorm()!=0;
        for(const auto &x:baseline.accepted.tangent_slip_correction)nonzero_travel_correction+=x.squaredNorm()!=0;
        require(nonzero_viscous && nonzero_tangent && nonzero_travel_correction,"nonzero-history fixture incomplete");
        std::cout<<"{\"status\":\"nonzero_body_component_codec_checks_passed_not_full_session_CP\",\"physical_advances\":0,\"reference_rebuilt\":false,\"generation\":"<<generation
            <<",\"bytes\":"<<input.size()<<",\"nonzero_viscous_nodes\":"<<nonzero_viscous<<",\"nonzero_tangent_rows\":"<<nonzero_tangent
            <<",\"nonzero_travel_corrections\":"<<nonzero_travel_correction<<",\"refusals\":[";
        for(size_t i=0;i<refusals.size();++i){if(i)std::cout<<',';std::cout<<'['<<quote(refusals[i].first)<<','<<quote(refusals[i].second)<<']';}
        std::cout<<"],\"native_identity\":{\"pid\":"<<pid<<",\"ppid\":"<<ppid
            <<",\"start_ticks\":"<<quote(start_ticks)<<",\"boot_id\":"<<quote(boot)<<",\"argv\":[";
        for(int i=0;i<argc;++i){if(i)std::cout<<',';std::cout<<quote(argv[i]);}
        std::cout<<"]},\"full_restart_qualified\":false}\n";return 0;
    }
    catch(const std::exception &e){std::cerr<<quote(e.what())<<'\n';return 1;}
}
