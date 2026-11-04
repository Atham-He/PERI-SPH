/** Real F-filesystem generation transactions, no physical advance/solver. */
#include "../src/shared/peridynamics/pd_bounded_generation_store.h"
#include <fstream>
#include <functional>
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
    require(s.size()==64,"external fingerprint length");Digest out{};
    auto digit=[](char c){if(c>='0'&&c<='9')return c-'0';if(c>='a'&&c<='f')return c-'a'+10;throw std::invalid_argument("hex");};
    for(size_t i=0;i<32;++i)out[i]=uint8_t(16*digit(s[2*i])+digit(s[2*i+1]));return out;
}
bool exists(const std::string &path)
{struct stat value{};if(::lstat(path.c_str(),&value)==0)return true;require(errno==ENOENT,"unexpected path stat error");return false;}
void mkdir(const std::string &path){if(::mkdir(path.c_str(),0700)!=0)gen::fail("engineering mkdir");}
void identity(int argc,char **argv)
{
    std::ifstream file("/proc/self/stat"),boot_file("/proc/sys/kernel/random/boot_id");std::string line,boot,word,ticks;
    std::getline(file,line);std::getline(boot_file,boot);std::istringstream stat(line.substr(line.rfind(')')+2));
    for(size_t i=0;i<=19;++i){require(bool(stat>>word),"native identity fields");if(i==19)ticks=word;}
    std::cout<<"\"native_identity\":{\"pid\":"<<::getpid()<<",\"ppid\":"<<::getppid()<<",\"start_ticks\":"<<quote(ticks)
        <<",\"boot_id\":"<<quote(boot)<<",\"argv\":[";
    for(int i=0;i<argc;++i){if(i)std::cout<<',';std::cout<<quote(argv[i]);}std::cout<<"]}";
}
void publishPhase(const std::string &root,const Digest &trusted,int argc,char **argv)
{
    Descriptor metadata(gen::openDirectory(root));const auto input=gen::readAt(metadata.get(),"session_candidate.bin",maximum_bytes).value();
    const auto candidate=decodeBoundedSessionCandidate(input,37,trusted);
    VerifiedBoundedMechanicalPrefix proof(root+"/physical.phy",candidate,37,trusted,true);
    BoundedGenerationStore store(root+"/store",trusted);require(!store.currentPointerBytes() && !store.live(),"new store not empty");
    const auto first=store.publish(candidate,37,proof,"g_37_initial",std::nullopt);
    require(first.status==BoundedPublishStatus::PublishedDurable && first.error_number==0 && !first.injected,"actual initial durable publish");
    auto baseline=store.currentPointerBytes();require(baseline.has_value(),"initial HEAD missing");
    const auto selection=decodeBoundedSelection(*baseline,37);require(selection.branch=="g_37_initial","initial branch selection");
    gen::writeExclusive(metadata.get(),"selection_initial.bin",*baseline);
    require(encodeBoundedSessionCandidate(store.live()->candidate(),37)==input,"initial owned bytes");
    std::vector<std::pair<std::string,std::string>> refused;
    auto expect=[&](const std::string &label,const std::string &reason,const auto &action,const std::string &absent_branch=std::string())
    {
        const auto *old=store.live();const auto owned=encodeBoundedSessionCandidate(old->candidate(),old->generation());
        const auto head=store.currentPointerBytes();
        try{action();}
        catch(const std::exception &e)
        {
            require(std::string(e.what()).find(reason)!=std::string::npos,"wrong refusal: "+label+": "+e.what());
            require(store.currentPointerBytes()==head && store.live()==old &&
                encodeBoundedSessionCandidate(store.live()->candidate(),store.live()->generation())==owned,"pre-publication failure changed live");
            if(!absent_branch.empty())require(!exists(root+"/store/"+absent_branch),"invalid input created private staging");
            refused.emplace_back(label,e.what());return;
        }
        throw std::runtime_error("invalid generation action accepted: "+label);
    };
    expect("duplicate_owner_lock","exclusive owner lock",[&]{BoundedGenerationStore other(root+"/store",trusted);});
    expect("unsafe_branch","unsafe generation branch",[&]{const auto r=store.publish(candidate,37,proof,"../escape",baseline);(void)r;});
    auto stale=*baseline;stale[0]^=1;
    expect("stale_expected","expected-current pointer mismatch",[&]{const auto r=store.publish(candidate,37,proof,"g_stale",stale);(void)r;},"g_stale");
    expect("wrong_explicit_generation","schema or explicit generation mismatch",[&]{const auto bad=decodeBoundedSelection(*baseline,38);(void)bad;});
    auto wrong_archive=selection;wrong_archive.archive_sha256[0]^=1;
    expect("wrong_selected_archive_digest","selected archive digest",[&]{const auto bad=store.openExplicit(wrong_archive);(void)bad;});
    auto wrong_contract=selection;wrong_contract.contract_sha256[0]^=1;
    expect("wrong_contract_selection","external contract selection mismatch",[&]{const auto bad=store.openExplicit(wrong_contract);(void)bad;});
    const Bytes protected_bytes{'k','e','e','p','_','p','r','o','t','e','c','t','e','d'};
    gen::writeExclusive(metadata.get(),"protected.bin",protected_bytes);
    require(::symlink((root+"/protected.bin").c_str(),(root+"/store/g_symlink").c_str())==0,"engineering branch symlink create");
    expect("symlink_branch","exclusive branch mkdir",[&]{const auto r=store.publish(candidate,37,proof,"g_symlink",baseline);(void)r;});
    require(gen::readAt(metadata.get(),"protected.bin",128).value()==protected_bytes,"protected symlink target changed");
    mkdir(root+"/store/g_exists");
    expect("existing_directory","exclusive branch mkdir",[&]{const auto r=store.publish(candidate,37,proof,"g_exists",baseline);(void)r;});
    auto invalid=candidate;invalid.body.accepted.ledger_correction[size_t(ContactLedger::MaxwellLoss)]+=.01;
    expect("semantic_bad_candidate","Maxwell cumulative history",[&]{const auto r=store.publish(invalid,37,proof,"g_bad_semantic",baseline);(void)r;},"g_bad_semantic");
    auto changed=candidate;changed.body.accepted.ledger_correction[size_t(ContactLedger::AlgorithmicLoss)]+=.001;
    expect("changed_after_proof","candidate changed after verification",[&]{const auto r=store.publish(changed,37,proof,"g_bad_binding",baseline);(void)r;},"g_bad_binding");
    expect("injected_before_stage","injected before stage",[&]{const auto r=store.publish(candidate,37,proof,"g_fault_stage",baseline,BoundedPublishFault::BeforeStage);(void)r;},"g_fault_stage");
    expect("injected_after_files","injected after private files",[&]{const auto r=store.publish(candidate,37,proof,"g_fault_files",baseline,BoundedPublishFault::AfterFiles);(void)r;});
    expect("injected_before_rename","injected before pointer rename",[&]{const auto r=store.publish(candidate,37,proof,"g_fault_rename",baseline,BoundedPublishFault::BeforeRename);(void)r;});
    Descriptor directory(gen::openDirectory(root+"/store"));gen::writeExclusive(directory.get(),"pending_g_pending.bin",protected_bytes);
    expect("existing_pending_pointer","exclusive file create",[&]{const auto r=store.publish(candidate,37,proof,"g_pending",baseline);(void)r;});
    mkdir(root+"/store/g_badarchive");
    require(::symlink((root+"/store/g_37_initial/state.bin").c_str(),(root+"/store/g_badarchive/state.bin").c_str())==0,"engineering archive symlink create");
    auto symbolic=selection;symbolic.branch="g_badarchive";
    expect("symlink_selected_archive","selected file open",[&]{const auto bad=store.openExplicit(symbolic);(void)bad;});
    mkdir(root+"/symlink_head_store");
    require(::symlink((root+"/selection_initial.bin").c_str(),(root+"/symlink_head_store/HEAD").c_str())==0,"engineering HEAD symlink create");
    expect("symlink_HEAD","selected file open",[&]
    {BoundedGenerationStore symbolic_store(root+"/symlink_head_store",trusted);const auto r=symbolic_store.publish(candidate,37,proof,"g_no_follow",std::nullopt);(void)r;});
    require(!exists(root+"/symlink_head_store/g_no_follow"),"symlink HEAD caused staging");
    auto wrong_prefix=selection;wrong_prefix.physical_prefix_sha256[0]^=1;
    expect("selected_prefix_identity","selected physical prefix identity",[&]{const auto bad=store.openExplicit(wrong_prefix);(void)bad;});
    Builder badshape(37);badshape.integers("SelectionSchema",{},{1});badshape.bytes("Branch",{2,3},{'g','_','a','b','c','d'});
    badshape.bytes("ArchiveSHA256",{32},Bytes(selection.archive_sha256.begin(),selection.archive_sha256.end()));
    badshape.bytes("PhysicalPrefixSHA256",{32},Bytes(selection.physical_prefix_sha256.begin(),selection.physical_prefix_sha256.end()));
    badshape.bytes("ContractSHA256",{32},Bytes(selection.contract_sha256.begin(),selection.contract_sha256.end()));
    const auto shaped=badshape.finish();gen::writeExclusive(metadata.get(),"selection_badshape.bin",shaped);
    expect("branch_path_shape","selected branch dimensions",[&]{const auto bad=decodeBoundedSelection(shaped,37);(void)bad;});
    const auto next=store.publish(candidate,37,proof,"g_37_next",baseline);
    require(next.status==BoundedPublishStatus::PublishedDurable && !store.live()->publicationUncertain(),"second durable publish");
    baseline=store.currentPointerBytes();gen::writeExclusive(metadata.get(),"selection_next.bin",*baseline);
    require(store.live()->branch()=="g_37_next" && encodeBoundedSessionCandidate(store.live()->candidate(),37)==input,"second owner bytes");
    // A synthetic abandoned-working-file tail, NOT a later accepted PD step.
    const std::string tail="synthetic_uncommitted_later_tail\n";
    Descriptor abandoned(::open((root+"/store/g_37_initial/active.phy").c_str(),O_WRONLY|O_APPEND|O_NOFOLLOW|O_CLOEXEC));
    require(::write(abandoned.get(),tail.data(),tail.size())==ssize_t(tail.size()),"synthetic tail append");gen::sync(abandoned.get(),"synthetic tail fsync");
    VerifiedBoundedMechanicalPrefix proof38(root+"/physical.phy",candidate,38,trusted,true);
    const auto uncertain=store.publish(candidate,38,proof38,"g_38_uncertain",baseline,BoundedPublishFault::AfterRenameDirectorySync);
    require(uncertain.status==BoundedPublishStatus::PublishedDirectoryDurabilityUncertain && uncertain.error_number==EIO && uncertain.injected,
        "controlled post-publish uncertainty result");
    require(store.live()->branch()=="g_38_uncertain" && store.live()->generation()==38 && store.live()->publicationUncertain(),
        "post-publish uncertainty did not select new owner");
    const auto selected_uncertain=store.currentPointerBytes();require(selected_uncertain!=baseline,"post-publish HEAD not selected");
    require(gen::readAt(directory.get(),"previous_g_38_uncertain.bin",2048).value()==*baseline,"previous pointer not preserved");
    gen::writeExclusive(metadata.get(),"selection_uncertain.bin",*selected_uncertain);
    expect("publication_uncertain_block","publication uncertain; explicit recovery required",[&]
    {const auto r=store.publish(candidate,37,proof,"g_forbidden_after_uncertainty",selected_uncertain);(void)r;},"g_forbidden_after_uncertainty");
    require(refused.size()==19,"generation negative-control count");
    std::cout<<"{\"status\":\"generation_publish_fault_controls_passed_NOT_physical_restart\",\"durable_publishes\":2,"
        <<"\"injected_post_publish_uncertainty\":true,\"uncertain_owner_known\":true,\"further_publish_refused\":true,"
        <<"\"source_generation\":37,\"terminal_selected_generation\":38,\"nonlinear_solver_calls\":0,\"physical_advances\":0,"
        <<"\"full_restart_qualified\":false,\"refusals\":[";
    for(size_t i=0;i<refused.size();++i){if(i)std::cout<<',';std::cout<<'['<<quote(refused[i].first)<<','<<quote(refused[i].second)<<']';}
    std::cout<<"],";identity(argc,argv);std::cout<<"}\n";
}
void restorePhase(const std::string &root,const Digest &trusted,int argc,char **argv)
{
    Descriptor metadata(gen::openDirectory(root));
    const auto initial=gen::readAt(metadata.get(),"selection_initial.bin",2048).value();
    const auto expected=gen::readAt(metadata.get(),"selection_uncertain.bin",2048).value();
    const auto selected=decodeBoundedSelection(initial,37);
    BoundedGenerationStore store(root+"/store",trusted);
    require(store.currentPointerBytes()==std::optional<Bytes>(expected) && !store.live(),"explicit recovery expected current");
    const auto loaded=store.openExplicit(selected);
    require(loaded.generation==37,"older generation was not explicitly selected");
    const auto result=store.publish(loaded.candidate,37,*loaded.prefix,"g_37_restored",expected);
    require(result.status==BoundedPublishStatus::PublishedDurable && result.error_number==0 && !store.live()->publicationUncertain(),"explicit restored publication");
    require(store.live()->branch()=="g_37_restored" && store.live()->generation()==37,"restored owner identity");
    const auto original=gen::readAt(metadata.get(),"session_candidate.bin",maximum_bytes).value();
    require(loaded.archive_bytes==original && encodeBoundedSessionCandidate(store.live()->candidate(),37)==original,"separate-process noadvance archive bytes");
    gen::writeExclusive(metadata.get(),"selection_restored.bin",store.currentPointerBytes().value());
    std::cout<<"{\"status\":\"explicit_older_generation_separate_process_restore_passed_NOT_physical_restart\","
        <<"\"explicit_generation\":37,\"ignored_current_generation\":38,\"archive_bytes\":"<<original.size()
        <<",\"physical_prefix_bytes\":"<<loaded.candidate.output.bytes<<",\"durable_publishes\":1,\"nonlinear_solver_calls\":0,"
        <<"\"physical_advances\":0,\"full_restart_qualified\":false,";identity(argc,argv);std::cout<<"}\n";
}
}
int main(int argc,char **argv)
{
    try
    {
        require(argc==4,"usage --publish|--restore ROOT TRUSTED_SHA256");const std::string mode=argv[1],root=argv[2];const auto trusted=unhex(argv[3]);
        if(mode=="--publish")publishPhase(root,trusted,argc,argv);
        else if(mode=="--restore")restorePhase(root,trusted,argc,argv);
        else throw std::invalid_argument("engineering mode");return 0;
    }
    catch(const std::exception &e){std::cerr<<quote(e.what())<<'\n';return 1;}
}
