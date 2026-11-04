/** Detached session contract only, no physical advance or live commit. */
#include "../src/shared/peridynamics/pd_bounded_session_archive.h"
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
void require(bool ok,const std::string &why){if(!ok)throw std::runtime_error(why);}
std::string quote(const std::string &text)
{
    std::ostringstream out;out<<'"';for(unsigned char c:text)
        if(c=='"'||c=='\\')out<<'\\'<<c;
        else if(c<32)out<<"\\u00"<<std::hex<<std::setw(2)<<std::setfill('0')<<unsigned(c)<<std::dec;
        else out<<c;
    out<<'"';return out.str();
}
Bytes load(const std::string &path)
{
    std::ifstream stream(path,std::ios::binary|std::ios::ate);require(bool(stream),"input open");
    const auto length=stream.tellg();require(length>=0 && uint64_t(length)<=maximum_bytes,"input bound");
    Bytes bytes(size_t(length),0);stream.seekg(0);stream.read(reinterpret_cast<char*>(bytes.data()),length);
    require(bool(stream),"input read");return bytes;
}
void store(const std::string &path,const Bytes &bytes)
{
    Descriptor target(::open(path.c_str(),O_WRONLY|O_CREAT|O_EXCL|O_CLOEXEC,0600));size_t offset=0;
    while(offset<bytes.size())
    {
        const auto n=::write(target.get(),bytes.data()+offset,bytes.size()-offset);
        if(n<0 && errno==EINTR)continue;
        require(n>0,"exclusive evidence write");offset+=size_t(n);
    }
    require(::fsync(target.get())==0,"exclusive evidence fsync");
}
Digest unhex(const std::string &s)
{
    require(s.size()==64,"expected independent trusted fingerprint");Digest digest{};
    auto digit=[](char c){if(c>='0'&&c<='9')return c-'0';if(c>='a'&&c<='f')return c-'a'+10;throw std::invalid_argument("hex");};
    for(size_t i=0;i<32;++i)digest[i]=uint8_t(16*digit(s[2*i])+digit(s[2*i+1]));return digest;
}
}
int main(int argc,char **argv)
{
    try
    {
        require(argc==4,"usage INPUT OUTPUT_DIR EXTERNAL_TRUSTED_CONTRACT_SHA256");
        const auto original=load(argv[1]);const std::string directory=argv[2];const auto trusted=unhex(argv[3]);
        const auto baseline=decodeBoundedSessionCandidate(original,37,trusted);
        require(baseline.body.accepted.gel.steps==40 && baseline.next_phase==3,"terminal nonzero fixture");
        require(encodeBoundedSessionCandidate(baseline,37)==original,"whole session candidate bytes");
        store(directory+"/roundtrip.bin",original);
        std::vector<std::pair<std::string,std::string>> refused;
        auto expect=[&](const std::string &label,const Bytes &bytes,const char *reason,uint64_t generation,const Digest &fingerprint)
        {
            try{auto candidate=decodeBoundedSessionCandidate(bytes,generation,fingerprint);(void)candidate;}
            catch(const std::exception &e)
            {
                require(std::string(e.what()).find(reason)!=std::string::npos,"wrong refusal: "+label+": "+e.what());
                refused.emplace_back(label,e.what());
                require(encodeBoundedSessionCandidate(baseline,37)==original,"baseline changed after refusal");return;
            }
            throw std::runtime_error("invalid session candidate accepted: "+label);
        };
        auto damage=[&](const std::string &label,const char *reason,const std::function<void(BoundedSessionCandidate&)> &change)
        {
            auto candidate=baseline;change(candidate);const auto bytes=encodeBoundedSessionCandidate(candidate,37);
            store(directory+"/"+label+".bin",bytes);expect(label,bytes,reason,37,trusted);
        };
        damage("source_identity","external trusted contract mismatch",[](auto &s){s.contract.source_manifest_sha256[0]^=1;});
        damage("ELF_identity","external trusted contract mismatch",[](auto &s){s.contract.ELF_sha256[0]^=1;});
        damage("build_runtime_identity","external trusted contract mismatch",[](auto &s){s.contract.build_runtime_manifest_sha256[0]^=1;});
        damage("input_identity","external trusted contract mismatch",[](auto &s){s.contract.input_sha256[0]^=1;});
        damage("protocol_tangent","external trusted contract mismatch",[](auto &s){s.contract.protocol.tangent_stroke[0]*=1.01;});
        damage("cap_width","external trusted contract mismatch",[](auto &s){s.contract.cap_width*=1.01;});
        damage("solver_gate","incompatible frozen solver or geometry gate",[](auto &s){s.contract.solver.relative_tolerance*=2;});
        damage("geometry_gate","incompatible frozen solver or geometry gate",[](auto &s){s.contract.penetration_stroke_gate*=2;});
        damage("algorithm_version","incompatible algorithm version",[](auto &s){++s.contract.equation_version;});
        damage("static_friction","static body fingerprint mismatch",[](auto &s){s.body.mu[0]*=2;});
        damage("static_surface_area","static body fingerprint mismatch",[](auto &s){s.body.surface.areas[0]*=1.01;});
        damage("accepted_time_ulp","exact accepted protocol clock",[](auto &s){s.body.accepted.gel.time=std::nextafter(s.body.accepted.gel.time,1);});
        damage("actual_dt_ulp","exact accepted protocol clock",[](auto &s){s.body.accepted.gel.last_dt=std::nextafter(s.body.accepted.gel.last_dt,1);});
        damage("accepted_boundary","accepted-only boundary",[](auto &s){s.accepted_boundary=0;});
        damage("next_phase","exact accepted protocol clock",[](auto &s){s.next_phase=2;});
        damage("output_records","structural output cursor",[](auto &s){++s.output.records;});
        damage("output_step","structural output cursor",[](auto &s){--s.output.last_step;});
        damage("output_time_ulp","structural output cursor",[](auto &s){s.output.last_time=std::nextafter(s.output.last_time,1);});
        damage("output_bytes","structural output cursor",[](auto &s){s.output.bytes=maximum_output_prefix+1;});
        Builder extra(37);FieldWriter fields(extra);auto copy=baseline;copy.registerCheckpoint(fields);extra.integers("unregistered",{},{1});
        const auto schema=extra.finish();store(directory+"/unknown_schema.bin",schema);
        expect("unknown_schema",schema,"exact schema names",37,trusted);
        auto wrong=trusted;wrong[0]^=1;
        expect("wrong_external_trust",original,"external trusted contract mismatch",37,wrong);
        expect("wrong_explicit_generation",original,"schema or explicit generation mismatch",38,trusted);
        const auto rounding=std::fegetround();require(std::fesetround(FE_DOWNWARD)==0,"set engineering round mode");
        try{expect("different_rounding_mode",original,"unsupported current floating point profile",37,trusted);}
        catch(...){std::fesetround(rounding);throw;}
        require(std::fesetround(rounding)==0,"restore engineering round mode");
#if defined(__x86_64__)
        const auto csr=_mm_getcsr();_mm_setcsr(csr|0x8000);
        try{expect("flush_subnormal_mode",original,"unsupported current floating point profile",37,trusted);}
        catch(...){_mm_setcsr(csr);throw;}
        _mm_setcsr(csr);
#endif
        require(refused.size()==24,"contract control count");
        // Legal protocol clocks must be checked at every accepted cursor,
        // including exact phase-end truncation, not just the final40-step time.
        size_t clock_checks=0;std::vector<std::pair<uint64_t,uint64_t>> clock_bits;
        for(int kind=0;kind<=2;++kind)
        {
            BoundedFrictionProtocol p=baseline.contract.protocol;p.kind=kind;Real time=0;
            for(size_t step=0;step<p.steps();++step)
            {
                const auto interval=p.interval(step,time);time=interval.new_time;
                require(interval.actual_dt>0 && time<=p.endTime(),"phase clock construction");
                clock_bits.emplace_back(doubleBits(time),doubleBits(interval.actual_dt));++clock_checks;
            }
            require(time==p.endTime(),"phase endpoint mismatch");
        }
        std::ifstream stat_file("/proc/self/stat"),boot_file("/proc/sys/kernel/random/boot_id");
        std::string line,boot,token,ticks;std::getline(stat_file,line);std::getline(boot_file,boot);
        std::istringstream stat(line.substr(line.rfind(')')+2));
        for(size_t i=0;i<=19;++i){require(bool(stat>>token),"native identity");if(i==19)ticks=token;}
        std::cout<<"{\"status\":\"detached_session_contract_component_checks_passed_NOT_output_or_restart_qualification\",\"physical_advances\":0,"
            <<"\"output_prefix_verified\":false,\"live_commits\":0,\"full_restart_qualified\":false,\"archive_bytes\":"<<original.size()
            <<",\"clock_checks\":"<<clock_checks<<",\"clock_bits\":[";
        for(size_t i=0;i<clock_bits.size();++i){if(i)std::cout<<',';std::cout<<'['<<clock_bits[i].first<<','<<clock_bits[i].second<<']';}
        std::cout<<"],\"refusals\":[";
        for(size_t i=0;i<refused.size();++i){if(i)std::cout<<',';std::cout<<'['<<quote(refused[i].first)<<','<<quote(refused[i].second)<<']';}
        std::cout<<"],\"native_identity\":{\"pid\":"<<::getpid()<<",\"ppid\":"<<::getppid()<<",\"start_ticks\":"<<quote(ticks)
            <<",\"boot_id\":"<<quote(boot)<<",\"argv\":[";
        for(int i=0;i<argc;++i){if(i)std::cout<<',';std::cout<<quote(argv[i]);}
        std::cout<<"]}}\n";return 0;
    }
    catch(const std::exception &e){std::cerr<<quote(e.what())<<'\n';return 1;}
}
