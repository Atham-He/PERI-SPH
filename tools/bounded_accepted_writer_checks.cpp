/** Replay saved accepted records through the REAL durable writer, no Newton. */
#include "../src/shared/peridynamics/pd_bounded_accepted_writer.h"
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
using namespace SPH;
using namespace SPH::peridynamics;
using namespace SPH::peridynamics::canonical_archive;
namespace gen=bounded_generation_detail;
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
Digest unhex(const std::string &s)
{
    require(s.size()==64,"external contract length");Digest d{};
    auto digit=[](char c){if(c>='0'&&c<='9')return c-'0';if(c>='a'&&c<='f')return c-'a'+10;throw std::invalid_argument("hex");};
    for(size_t i=0;i<32;++i)d[i]=uint8_t(16*digit(s[2*i])+digit(s[2*i+1]));return d;
}
void mkdir(const std::string &path){if(::mkdir(path.c_str(),0700)!=0)gen::fail("engineering replay mkdir");}
bool exists(const std::string &path)
{struct stat s{};if(::lstat(path.c_str(),&s)==0)return true;require(errno==ENOENT,"engineering stat");return false;}
void readExact(int fd,uint64_t offset,uint8_t *data,size_t count,uint64_t bound)
{
    require(offset<=bound && count<=bound-offset,"saved record read bound");size_t read=0;
    while(read<count){const auto n=::pread(fd,data+read,count-read,off_t(offset+read));if(n<0&&errno==EINTR)continue;require(n>0,"saved record read");read+=size_t(n);}
}
BoundedMechanicalRecord recordAt(const std::string &root,uint64_t step,const Digest &trusted)
{
    Descriptor metadata(gen::openDirectory(root));Descriptor file(::openat(metadata.get(),"reference_stream.phy",O_RDONLY|O_NOFOLLOW|O_CLOEXEC));
    const auto before=gen::regular(file.get(),maximum_output_prefix);const auto header=boundedPhysicalStreamHeader(trusted);
    Bytes actual(header.size());readExact(file.get(),0,actual.data(),actual.size(),before.st_size);require(actual==header,"saved record trusted header");
    uint64_t offset=48;
    for(uint64_t i=0;i<=step;++i)
    {
        std::array<uint8_t,8> length_bytes{};readExact(file.get(),offset,length_bytes.data(),8,before.st_size);offset+=8;
        uint64_t length=0;for(size_t k=0;k<8;++k)length|=uint64_t(length_bytes[k])<<(8*k);
        require(length>=72 && length<=maximum_bytes && offset<=uint64_t(before.st_size) && length<=uint64_t(before.st_size)-offset,"saved record length");
        if(i==step)
        {Bytes bytes(size_t(length),0);readExact(file.get(),offset,bytes.data(),bytes.size(),before.st_size);gen::unchanged(file.get(),before,maximum_output_prefix);return decodeBoundedMechanicalRecord(bytes,step);}
        offset+=length;
    }
    throw std::runtime_error("saved record absent");
}
BoundedSessionCandidate initial(const std::string &root,const Digest &trusted)
{Descriptor d(gen::openDirectory(root));return decodeBoundedSessionCandidate(gen::readAt(d.get(),"initial_candidate.bin",maximum_bytes).value(),1,trusted);}
BoundedRunningContract bindRunning(const std::string &root,const BoundedSessionContract &contract)
{return BoundedRunningContract(contract,root+"/source_manifest.json",root+"/build_runtime_manifest.json",root+"/original.input");}
void start(BoundedGenerationStore &store,const std::string &root,const BoundedSessionCandidate &candidate,const Digest &trusted)
{
    VerifiedBoundedMechanicalPrefix proof(root+"/initial.phy",candidate,1,trusted,true);
    const auto result=store.publish(candidate,1,proof,"g_1_initial",std::nullopt);
    require(result.status==BoundedPublishStatus::PublishedDurable,"initial writer checkpoint publication");
}
void appendRange(BoundedGenerationStore &store,const std::string &root,const BoundedRunningContract &running,const Digest &trusted,uint64_t first,uint64_t last)
{
    for(uint64_t step=first;step<=last;++step)
    {
        const auto record=recordAt(root,step,trusted);const auto result=BoundedAcceptedWriter::append(store,record,running);
        require(result.status==BoundedAppendStatus::AcceptedDurable && result.witness.committed &&
            result.witness.written_bytes==result.witness.frame_bytes && store.live()->candidate().body.accepted.gel.steps==step,
            "saved accepted-record commit failed");
    }
}
Bytes seal(BoundedGenerationStore &store,const Digest &trusted,uint64_t generation,const std::string &branch)
{
    const auto &candidate=store.live()->candidate();VerifiedBoundedMechanicalPrefix proof(store.live()->workingSourcePath(),candidate,generation,trusted,true);
    const auto result=store.publish(candidate,generation,proof,branch,store.currentPointerBytes());
    require(result.status==BoundedPublishStatus::PublishedDurable,"replay checkpoint seal");return store.currentPointerBytes().value();
}
void exported(const std::string &root,const std::string &name,const Bytes &bytes)
{Descriptor d(gen::openDirectory(root));gen::writeExclusive(d.get(),name,bytes);}
void nativeIdentity(int argc,char **argv)
{
    std::ifstream stat_file("/proc/self/stat"),boot_file("/proc/sys/kernel/random/boot_id");std::string line,boot,word,ticks;
    std::getline(stat_file,line);std::getline(boot_file,boot);std::istringstream stat(line.substr(line.rfind(')')+2));
    for(size_t i=0;i<=19;++i){require(bool(stat>>word),"native identity");if(i==19)ticks=word;}
    std::cout<<"\"native_identity\":{\"pid\":"<<::getpid()<<",\"ppid\":"<<::getppid()<<",\"start_ticks\":"<<quote(ticks)
        <<",\"boot_id\":"<<quote(boot)<<",\"argv\":[";
    for(int i=0;i<argc;++i){if(i)std::cout<<',';std::cout<<quote(argv[i]);}std::cout<<"]}";
}
void replayMode(const std::string &mode,const std::string &root,const Digest &trusted,int argc,char **argv)
{
    const auto begin=initial(root,trusted);const auto running=bindRunning(root,begin.contract);
    const bool continuous=mode=="--continuous",first=mode=="--split-first",rollback=mode=="--rollback";
    const auto store_path=root+(continuous?"/continuous":"/split");BoundedGenerationStore store(store_path,trusted);
    uint64_t from=1,to=first?16:40;
    if(continuous || first)start(store,root,begin,trusted);
    else
    {
        Descriptor d(gen::openDirectory(root));
        const auto selected_bytes=gen::readAt(d.get(),"split_selection_16.bin",2048).value();
        const auto expected=rollback?gen::readAt(d.get(),"split_selection_100.bin",2048).value():selected_bytes;
        require(store.currentPointerBytes()==std::optional<Bytes>(expected),"explicit replay resume expected current");
        const auto loaded=store.openExplicit(decodeBoundedSelection(selected_bytes,16));
        const auto published=store.publish(loaded.candidate,16,*loaded.prefix,rollback?"g_16_rollback":"g_16_resumed",expected);
        require(published.status==BoundedPublishStatus::PublishedDurable,"explicit replay resume publish");
        exported(root,rollback?"rollback_noadvance.bin":"split_noadvance.bin",encodeBoundedSessionCandidate(store.live()->candidate(),16));
        from=17;
    }
    appendRange(store,root,running,trusted,from,to);
    const auto selected=seal(store,trusted,first?16:100,first?"g_16_mid":(rollback?"g_100_rollback_final":"g_100_complete"));
    const auto prefix=continuous?"continuous":(first?"split_first":(rollback?"rollback":"split_resumed"));
    exported(root,std::string(prefix)+"_final.bin",encodeBoundedSessionCandidate(store.live()->candidate(),first?16:100));
    if(first)exported(root,"split_selection_16.bin",selected);
    else if(!continuous && !rollback)exported(root,"split_selection_100.bin",selected);
    std::cout<<"{\"status\":\"saved_record_durable_replay_phase_passed_NOT_new_physical_trajectory\",\"mode\":"<<quote(mode)
        <<",\"accepted_record_replays\":"<<(to-from+1)<<",\"terminal_saved_step\":"<<to
        <<",\"running_replay_ELF_identity_bound\":true,\"nonlinear_solver_calls\":0,\"new_physical_steps\":0,\"full_restart_qualified\":false,";
    nativeIdentity(argc,argv);std::cout<<"}\n";
}
std::string witnessJSON(const BoundedAppendWitness &w)
{
    std::ostringstream out;out<<std::setprecision(17)<<"{\"stage\":"<<unsigned(w.stage)<<",\"failure_stage\":"<<unsigned(w.failure_stage)
        <<",\"old_step\":"<<w.old_step<<",\"proposed_step\":"<<w.proposed_step<<",\"old_time\":"<<w.old_time<<",\"proposed_time\":"<<w.proposed_time
        <<",\"old_committed_bytes\":"<<w.old_committed_bytes<<",\"frame_bytes\":"<<w.frame_bytes<<",\"written_bytes\":"<<w.written_bytes
        <<",\"error_number\":"<<w.error_number<<",\"injected\":"<<(w.injected?"true":"false")<<",\"committed\":"<<(w.committed?"true":"false")<<'}';return out.str();
}
void faults(const std::string &root,const Digest &trusted,int argc,char **argv)
{
    const auto begin=initial(root,trusted);const auto running=bindRunning(root,begin.contract);const auto valid=recordAt(root,1,trusted);
    const auto original=encodeBoundedSessionCandidate(begin,1);std::vector<std::string> reports;size_t blocked=0;
    const std::vector<std::string> labels{"wrong_ELF","wrong_source","wrong_build","wrong_input","wrong_bound_contract","wrong_FP",
        "invalid_reaction","wrong_old_hash","extra_old_tail","before_append","partial_write","full_write_no_sync","sync_error"};
    for(const auto &label:labels)
    {
        const auto path=root+"/faults/"+label;mkdir(path);BoundedGenerationStore store(path,trusted);start(store,root,begin,trusted);
        const auto head=store.currentPointerBytes();const auto *owner=store.live();std::string reason;bool binding_case=false,returned=false;
        int status=-1;auto record=valid;auto fault=BoundedAppendFault::None;
        if(label=="invalid_reaction")record.observation.reaction[0][0]+=.01;
        if(label=="before_append")fault=BoundedAppendFault::BeforeAppend;
        if(label=="partial_write")fault=BoundedAppendFault::AfterPartialWrite;
        if(label=="full_write_no_sync")fault=BoundedAppendFault::AfterFullWriteBeforeSync;
        if(label=="sync_error")fault=BoundedAppendFault::FsyncError;
        if(label=="wrong_old_hash" || label=="extra_old_tail")
        {
            Descriptor file(::open(store.live()->workingSourcePath().c_str(),O_RDWR|O_NOFOLLOW|O_CLOEXEC));
            if(label=="wrong_old_hash"){const uint8_t byte='Q';require(::pwrite(file.get(),&byte,1,0)==1,"engineering prefix mutation");}
            else{const std::string tail="UNSELECTED\n";require(::pwrite(file.get(),tail.data(),tail.size(),off_t(begin.output.bytes))==ssize_t(tail.size()),"engineering extra tail");}
            gen::sync(file.get(),"engineering prefix mutation fsync");
        }
        try
        {
            if(label=="wrong_ELF" || label=="wrong_source" || label=="wrong_build" || label=="wrong_input")
            {
                binding_case=true;auto contract=begin.contract;
                if(label=="wrong_ELF")contract.ELF_sha256[0]^=1;
                if(label=="wrong_source")contract.source_manifest_sha256[0]^=1;
                if(label=="wrong_build")contract.build_runtime_manifest_sha256[0]^=1;
                if(label=="wrong_input")contract.input_sha256[0]^=1;
                const auto bad=bindRunning(root,contract);(void)bad;
            }
            else if(label=="wrong_bound_contract")
            {
                auto other=begin.contract;other.protocol.tangent_stroke[0]*=1.01;const auto bad=bindRunning(root,other);
                const auto result=BoundedAcceptedWriter::append(store,record,bad);(void)result;
            }
            else if(label=="wrong_FP")
            {
                const auto mode=std::fegetround();require(std::fesetround(FE_DOWNWARD)==0,"engineering rounding set");
                try{const auto result=BoundedAcceptedWriter::append(store,record,running);(void)result;}
                catch(...){std::fesetround(mode);throw;}require(std::fesetround(mode)==0,"engineering rounding restore");
            }
            else
            {
                const auto result=BoundedAcceptedWriter::append(store,record,running,fault);status=int(result.status);returned=true;
                require(result.status!=BoundedAppendStatus::AcceptedDurable && !result.witness.committed,"fault accepted saved candidate");
                reason="injected_append_failure";
            }
        }
        catch(const std::exception &e){reason=e.what();}
        const std::map<std::string,std::string> expected{{"wrong_ELF","actual executable digest mismatch"},{"wrong_source","source manifest digest mismatch"},
            {"wrong_build","build runtime manifest digest mismatch"},{"wrong_input","original input digest mismatch"},{"wrong_bound_contract","bound contract mismatch"},
            {"wrong_FP","unsupported current floating point profile"},{"invalid_reaction","physical observation byte mismatch"},
            {"wrong_old_hash","old committed prefix digest mismatch"},{"extra_old_tail","exact committed working length"},{"before_append","injected before append"},
            {"partial_write","injected_append_failure"},{"full_write_no_sync","injected_append_failure"},{"sync_error","injected_append_failure"}};
        require(reason.find(expected.at(label))!=std::string::npos,"wrong accepted-writer refusal: "+label+": "+reason);
        require(store.currentPointerBytes()==head && store.live()==owner && encodeBoundedSessionCandidate(owner->candidate(),1)==original,"fault changed old accepted owner/HEAD");
        require(owner->candidate().body.accepted.gel.steps==0 && !owner->lastAppend().committed,"fault accepted record");
        const auto witness=witnessJSON(owner->lastAppend());bool blocked_append=false,blocked_publish=false;
        if(!binding_case)
        {
            require(owner->writerStopped(),"failed accepted writer was not stopped");
            try{const auto result=BoundedAcceptedWriter::append(store,valid,running);(void)result;}
            catch(const std::exception &e){blocked_append=std::string(e.what()).find("stopped; explicit recovery required")!=std::string::npos;}
            VerifiedBoundedMechanicalPrefix proof(path+"/g_1_initial/prefix.phy",begin,1,trusted,true);
            try{const auto result=store.publish(begin,1,proof,"g_forbidden",head);(void)result;}
            catch(const std::exception &e){blocked_publish=std::string(e.what()).find("accepted writer stopped; explicit recovery required")!=std::string::npos;}
            require(blocked_append && blocked_publish && !exists(path+"/g_forbidden"),"failed writer continuation/publish was allowed");
            require(witnessJSON(owner->lastAppend())==witness,"first append failure witness changed");blocked+=2;
        }
        reports.push_back("{\"label\":"+quote(label)+",\"reason\":"+quote(reason)+",\"returned_failure\":"+(returned?"true":"false")+
            ",\"status\":"+std::to_string(status)+",\"writer_stopped\":"+(owner->writerStopped()?"true":"false")+",\"witness\":"+witness+"}");
    }
    require(reports.size()==13 && blocked==18,"writer control count");
    std::cout<<"{\"status\":\"accepted_writer_identity_prepare_and_partial_tail_controls_passed_NOT_physical_restart\",\"primary_controls\":13,"
        <<"\"blocked_continuation_controls\":18,\"nonlinear_solver_calls\":0,\"new_physical_steps\":0,\"full_restart_qualified\":false,\"controls\":[";
    for(size_t i=0;i<reports.size();++i){if(i)std::cout<<',';std::cout<<reports[i];}std::cout<<"],";nativeIdentity(argc,argv);std::cout<<"}\n";
}
}
int main(int argc,char **argv)
{
    try
    {
        require(argc==4,"usage MODE ROOT TRUSTED_SHA256");const std::string mode=argv[1],root=argv[2];const auto trusted=unhex(argv[3]);
        if(mode=="--faults")faults(root,trusted,argc,argv);
        else if(mode=="--continuous"||mode=="--split-first"||mode=="--resume"||mode=="--rollback")replayMode(mode,root,trusted,argc,argv);
        else throw std::invalid_argument("saved-record replay mode");return 0;
    }
    catch(const std::exception &e){std::cerr<<quote(e.what())<<'\n';return 1;}
}
