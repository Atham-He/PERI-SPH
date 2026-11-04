/** Readonly existing-data mechanical stream qualification; zero Newton calls. */
#include "../src/shared/peridynamics/pd_bounded_mechanical_stream.h"
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
using namespace SPH;
using namespace SPH::peridynamics;
using namespace SPH::peridynamics::canonical_archive;
namespace
{
void require(bool ok,const std::string &why){if(!ok)throw std::runtime_error(why);}
std::string quote(const std::string &s)
{
    std::ostringstream out;out<<'"';for(unsigned char c:s)
        if(c=='"'||c=='\\')out<<'\\'<<c;
        else if(c<32)out<<"\\u00"<<std::hex<<std::setw(2)<<std::setfill('0')<<unsigned(c)<<std::dec;
        else out<<c;
    out<<'"';return out.str();
}
Bytes load(const std::string &path)
{
    std::ifstream file(path,std::ios::binary|std::ios::ate);require(bool(file),"candidate input open");
    const auto size=file.tellg();require(size>=0 && uint64_t(size)<=maximum_bytes,"candidate input bound");
    Bytes bytes(size_t(size),0);file.seekg(0);file.read(reinterpret_cast<char*>(bytes.data()),size);
    require(bool(file),"candidate input read");return bytes;
}
Digest unhex(const std::string &s)
{
    require(s.size()==64,"trusted SHA256 length");Digest result{};
    auto digit=[](char c){if(c>='0'&&c<='9')return c-'0';if(c>='a'&&c<='f')return c-'a'+10;throw std::invalid_argument("hex");};
    for(size_t i=0;i<32;++i)result[i]=uint8_t(16*digit(s[2*i])+digit(s[2*i+1]));return result;
}
void leaf(const std::string &s)
{require(!s.empty() && s.size()<128 && s.find_first_not_of("abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789_.-")==std::string::npos,"unsafe case leaf");}
}
int main(int argc,char **argv)
{
    try
    {
        require(argc==3,"usage ROOT EXTERNAL_TRUSTED_CONTRACT_SHA256");const std::string root=argv[1];const auto trusted=unhex(argv[2]);
        const auto input=load(root+"/session_candidate.bin");const auto candidate=decodeBoundedSessionCandidate(input,37,trusted);
        VerifiedBoundedMechanicalPrefix baseline(root+"/physical.phy",candidate,37,trusted,true);
        Descriptor stage(::open((root+"/staged").c_str(),O_RDONLY|O_DIRECTORY|O_CLOEXEC|O_NOFOLLOW));
        baseline.stageExclusiveCopy(candidate,37,stage.get(),"physical.phy");
        VerifiedBoundedMechanicalPrefix later(root+"/later_tail.phy",candidate,37,trusted,false);
        later.stageExclusiveCopy(candidate,37,stage.get(),"from_later_prefix.phy");
        std::vector<std::pair<std::string,std::string>> refusals;
        auto expect=[&](const std::string &label,const std::string &reason,const auto &action)
        {
            try{action();}
            catch(const std::exception &e)
            {
                require(std::string(e.what()).find(reason)!=std::string::npos,"wrong refusal: "+label+": "+e.what());
                refusals.emplace_back(label,e.what());
                require(encodeBoundedSessionCandidate(candidate,37)==input,"baseline candidate changed");return;
            }
            throw std::runtime_error("invalid physical stream accepted: "+label);
        };
        std::ifstream cases(root+"/cases.tsv");require(bool(cases),"cases open");std::string line;
        size_t file_cases=0;
        while(std::getline(cases,line))
        {
            std::istringstream fields(line);std::string label,archive,stream,reason;
            require(bool(std::getline(fields,label,'\t')) && bool(std::getline(fields,archive,'\t')) &&
                bool(std::getline(fields,stream,'\t')) && bool(std::getline(fields,reason)),"case columns");
            leaf(label);leaf(archive);leaf(stream);
            expect(label,reason,[&]
            {
                const auto bad=decodeBoundedSessionCandidate(load(root+"/evidence/"+archive),37,trusted);
                VerifiedBoundedMechanicalPrefix proof(root+"/evidence/"+stream,bad,37,trusted,true);(void)proof;
            });++file_cases;
        }
        require(file_cases==22,"physical stream file-control count");
        expect("exact_snapshot_later_tail","explicit byte cursor or snapshot size",[&]
        {VerifiedBoundedMechanicalPrefix proof(root+"/later_tail.phy",candidate,37,trusted,true);(void)proof;});
        auto changed=candidate;changed.body.accepted.ledger_correction[size_t(ContactLedger::AlgorithmicLoss)]+=.001;
        expect("changed_candidate_after_verification","candidate changed after verification",[&]
        {baseline.stageExclusiveCopy(changed,37,stage.get(),"must_not_exist.phy");});
        VerifiedOutputPrefix bytes(root+"/physical.phy",candidate.output.bytes,bounded_session_detail::digestArray(candidate.output.prefix_sha256),true);
        expect("semantic_read_bound","semantic read bound",[&]
        {uint8_t value=0;bytes.readAt(candidate.output.bytes,&value,1);});
        require(refusals.size()==25,"all physical stream controls");
        std::ifstream stat_file("/proc/self/stat"),boot_file("/proc/sys/kernel/random/boot_id");
        std::string stat,boot,word,ticks;std::getline(stat_file,stat);std::getline(boot_file,boot);
        std::istringstream fields(stat.substr(stat.rfind(')')+2));
        for(size_t i=0;i<=19;++i){require(bool(fields>>word),"native identity");if(i==19)ticks=word;}
        std::cout<<"{\"status\":\"converted_existing_physical_prefix_semantics_and_scratch_proposals_passed_NOT_physical_restart\","
            <<"\"records\":"<<candidate.output.records<<",\"accepted_steps\":"<<candidate.body.accepted.gel.steps
            <<",\"bytes\":"<<candidate.output.bytes<<",\"positive_prefixes\":2,\"nonlinear_solver_calls\":0,\"physical_advances\":0,"
            <<"\"live_commits\":0,\"filesystem_transaction_qualified\":false,\"full_restart_qualified\":false,\"refusals\":[";
        for(size_t i=0;i<refusals.size();++i){if(i)std::cout<<',';std::cout<<'['<<quote(refusals[i].first)<<','<<quote(refusals[i].second)<<']';}
        std::cout<<"],\"native_identity\":{\"pid\":"<<::getpid()<<",\"ppid\":"<<::getppid()<<",\"start_ticks\":"<<quote(ticks)
            <<",\"boot_id\":"<<quote(boot)<<",\"argv\":[";
        for(int i=0;i<argc;++i){if(i)std::cout<<',';std::cout<<quote(argv[i]);}
        std::cout<<"]}}\n";return 0;
    }
    catch(const std::exception &e){std::cerr<<quote(e.what())<<'\n';return 1;}
}
