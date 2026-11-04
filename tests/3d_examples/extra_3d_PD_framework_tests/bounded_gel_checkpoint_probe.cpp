/** Actual bounded PD checkpoint driver. Qualification is supplied by external
 * frozen campaigns, never by this executable's capability or terminal flags. */
#include "bounded_gel_friction_observer.h"
#include "pd_bounded_accepted_writer.h"
#include <charconv>
using namespace SPH;using namespace SPH::peridynamics;
using namespace SPH::peridynamics::canonical_archive;
namespace gen=bounded_generation_detail;
namespace obs=bounded_friction_observer;
namespace io=bounded_probe_io;
namespace
{
void require(bool ok,const char *why){if(!ok)throw std::invalid_argument(why);}
uint64_t integer(const char *value)
{
    uint64_t result=0;const char *last=value+std::strlen(value);
    const auto parsed=std::from_chars(value,last,result);
    require(parsed.ec==std::errc{} && parsed.ptr==last,"checkpoint driver: unsigned integer");return result;
}
Digest unhex(const std::string &s)
{
    require(s.size()==64,"checkpoint driver: trusted digest length");Digest d{};
    auto digit=[](char c){if(c>='0'&&c<='9')return c-'0';if(c>='a'&&c<='f')return c-'a'+10;throw std::invalid_argument("checkpoint driver: digest hex");};
    for(size_t i=0;i<32;++i)d[i]=uint8_t(16*digit(s[2*i])+digit(s[2*i+1]));return d;
}
std::string hex(const Digest &digest)
{
    static const char alphabet[]="0123456789abcdef";std::string out;
    for(uint8_t byte:digest){out.push_back(alphabet[byte>>4]);out.push_back(alphabet[byte&15]);}return out;
}
StdVec<uint8_t> hashPath(const std::string &path,bool executable=false)
{
    Descriptor file(::open(path.c_str(),O_RDONLY|O_CLOEXEC|(executable?0:O_NOFOLLOW)));
    const auto before=gen::regular(file.get(),executable?256*1024*1024:maximum_bytes);
    const auto digest=gen::hashFile(file.get(),uint64_t(before.st_size));return {digest.begin(),digest.end()};
}
void writeAt(const std::string &root,const std::string &leaf,const Bytes &bytes)
{Descriptor directory(gen::openDirectory(root));gen::writeExclusive(directory.get(),leaf,bytes);gen::sync(directory.get(),"checkpoint driver: output parent fsync");}
Bytes readAt(const std::string &root,const std::string &leaf,uint64_t bound=maximum_bytes)
{Descriptor directory(gen::openDirectory(root));return gen::readAt(directory.get(),leaf,bound).value();}
BoundedRunningContract running(const std::string &root,const BoundedSessionContract &contract)
{return BoundedRunningContract(contract,root+"/source_manifest.json",root+"/build_runtime_manifest.json",root+"/original.input");}

/** Initialization is a separate, zero-advance preflight. Its entire output and
 * external contract are frozen BEFORE any actual physical command is allowed.
 * Neither resume nor noadvance calls this function or constructs a scene. */
void initialize(const std::string &root)
{
    const auto input=obs::readInput((root+"/original.input").c_str());const BoundedNormalScene scene(input.scene);
    BoundedSessionCandidate candidate;auto &body=candidate.body;
    body.model=scene.model;body.surface=scene.surface;body.dense_surface=scene.dense_surface;
    body.frame=scene.frame;body.normal=scene.normal;body.tangent={scene.frame*Vecd::UnitX(),scene.frame*Vecd::UnitZ()};
    body.mu=input.mu;body.plate_ids=scene.plateIDs();body.accepted=initialBoundedGelFrictionState(scene);
    auto &contract=candidate.contract;contract.rotation=input.scene.rotation;contract.azimuths=input.scene.azimuths;
    contract.cap_width=input.scene.cap_width;contract.protocol=input.motion;
    contract.source_manifest_sha256=hashPath(root+"/source_manifest.json");contract.ELF_sha256=hashPath("/proc/self/exe",true);
    contract.build_runtime_manifest_sha256=hashPath(root+"/build_runtime_manifest.json");contract.input_sha256=hashPath(root+"/original.input");
    const auto static_digest=bounded_session_detail::staticBodyDigest(body);contract.static_body_sha256.assign(static_digest.begin(),static_digest.end());
    const auto bound=running(root,contract);bound.requireFor(contract);const auto trusted=boundedContractDigest(contract);
    BoundedMechanicalRecord initial;initial.body=body;bounded_mechanical_detail::validateInitial(initial);
    const auto header=boundedPhysicalStreamHeader(trusted);auto physical=header;const auto record=encodeBoundedMechanicalRecord(initial);
    appendU64(physical,record.size());physical.insert(physical.end(),record.begin(),record.end());
    const auto header_sha=sha256(header,header.size()),prefix_sha=sha256(physical,physical.size());
    candidate.output.records=1;candidate.output.bytes=physical.size();
    candidate.output.header_sha256.assign(header_sha.begin(),header_sha.end());candidate.output.prefix_sha256.assign(prefix_sha.begin(),prefix_sha.end());
    candidate.validateCandidate(trusted);
    writeAt(root,"initial_candidate.bin",encodeBoundedSessionCandidate(candidate,1));writeAt(root,"initial.phy",physical);
    writeAt(root,"contract.bin",boundedContractBytes(contract));
    // Read all emitted records back through the actual detached semantic loader.
    const auto decoded=decodeBoundedSessionCandidate(readAt(root,"initial_candidate.bin"),1,trusted);
    VerifiedBoundedMechanicalPrefix proof(root+"/initial.phy",decoded,1,trusted,true);
    io::begin("initialization");io::field("contract_sha256",hex(trusted));io::field("physical_bytes",physical.size());
    io::field("new_physical_steps",0);io::field("full_restart_qualified",false);io::finish();
}

struct DurableDriverCursor
{
    uint64_t durable_step=0,durable_bytes=0,debug_step=0,solver_calls=0;
    Real durable_time=0,debug_time=0;
    BoundedAppendWitness last_append{};
    int publication_error_number=0;
    unsigned publication_stage=0; // 1 initial,2 explicit resume,3 later seal
    void accepted(const BoundedOwnedCheckpoint *owner)noexcept
    {
        if(!owner)return;const auto &c=owner->candidate();
        durable_step=c.output.last_step;durable_time=c.output.last_time;durable_bytes=c.output.bytes;
        if(owner->lastAppend().stage!=BoundedAppendStage::Idle)last_append=owner->lastAppend();
    }
};
// Allocation-free first failure witness, emitted while the owner still exists.
// This remains useful if full Newton JSON formatting failed after fsync/commit.
void diagnostic(const std::exception &error,const DurableDriverCursor &cursor,const BoundedOwnedCheckpoint *owner)noexcept
{
    std::fprintf(stderr,"{\"type\":\"checkpoint_driver_failure\",\"mechanical_durable_step\":%llu,\"mechanical_durable_time\":%.17g,\"mechanical_durable_bytes\":%llu,\"last_full_debug_step\":%llu,\"last_full_debug_time\":%.17g,\"actual_solver_calls\":%llu,\"publication_error_number\":%d,\"publication_stage\":%u,\"message\":\"",
        (unsigned long long)cursor.durable_step,cursor.durable_time,(unsigned long long)cursor.durable_bytes,
        (unsigned long long)cursor.debug_step,cursor.debug_time,(unsigned long long)cursor.solver_calls,cursor.publication_error_number,cursor.publication_stage);
    for(const unsigned char *p=reinterpret_cast<const unsigned char*>(error.what());*p;++p)
    {if(*p=='"'||*p=='\\')std::fputc('\\',stderr);if(*p<32)std::fprintf(stderr,"\\u%04x",unsigned(*p));else std::fputc(*p,stderr);}
    std::fputs("\",\"append_witness\":",stderr);
    if(!owner)std::fputs("null",stderr);
    else
    {
        const auto &w=owner->lastAppend().stage==BoundedAppendStage::Idle ? cursor.last_append : owner->lastAppend();
        std::fprintf(stderr,"{\"stage\":%u,\"failure_stage\":%u,\"old_step\":%llu,\"proposed_step\":%llu,\"old_time\":%.17g,\"proposed_time\":%.17g,\"old_committed_bytes\":%llu,\"frame_bytes\":%llu,\"written_bytes\":%llu,\"error_number\":%d,\"injected\":%s,\"committed\":%s,\"writer_stopped\":%s,\"publication_uncertain\":%s}",
            unsigned(w.stage),unsigned(w.failure_stage),(unsigned long long)w.old_step,(unsigned long long)w.proposed_step,w.old_time,w.proposed_time,
            (unsigned long long)w.old_committed_bytes,(unsigned long long)w.frame_bytes,(unsigned long long)w.written_bytes,w.error_number,
            w.injected?"true":"false",w.committed?"true":"false",owner->writerStopped()?"true":"false",owner->publicationUncertain()?"true":"false");
    }
    std::fputs("}\n",stderr);std::fflush(stderr);
}
void seal(BoundedGenerationStore &store,const std::string &store_path,const Digest &trusted,const std::string &label,DurableDriverCursor &cursor)
{
    const auto &candidate=store.live()->candidate();const uint64_t step=candidate.body.accepted.gel.steps,generation=step?step:1;
    VerifiedBoundedMechanicalPrefix proof(store.live()->workingSourcePath(),candidate,generation,trusted,true);
    const auto branch="g_"+label+"_s"+std::to_string(step);
    cursor.publication_stage=3;
    const auto result=store.publish(candidate,generation,proof,branch,store.currentPointerBytes());
    cursor.publication_error_number=result.error_number;cursor.accepted(store.live());
    require(result.status==BoundedPublishStatus::PublishedDurable,"checkpoint driver: checkpoint publication durability uncertain");
    writeAt(store_path,"selected_"+label+"_s"+std::to_string(step)+".bin",store.currentPointerBytes().value());
    io::begin("checkpoint");io::field("generation",generation);io::field("accepted_step",step);
    io::field("branch",branch);io::field("full_restart_qualified",false);io::finish();
}

int execute(int argc,char **argv)
{
    const std::string mode=argv[1];const bool fresh=mode=="--fresh",validate=mode=="--validate";
    require((fresh && argc==7)||(!fresh && !validate && mode=="--resume" && argc==10)||(validate && argc==8),"checkpoint driver: execution argument count");
    const std::string root=argv[2],store_path=argv[3];const Digest trusted=unhex(argv[4]);
    uint64_t target=0,generation=0;std::string label;std::optional<Bytes> expected;Bytes selected;
    if(fresh){target=integer(argv[5]);label=argv[6];}
    else
    {
        // Selection filename is an explicit safe leaf under this store; never
        // search HEAD/branches for a replacement on any validation failure.
        selected=readAt(store_path,argv[5],2048);generation=integer(argv[6]);
        if(!validate){target=integer(argv[8]);label=argv[9];}
    }
    if(!validate){gen::branchLeaf("g_"+label+"_s40");require(target>0 && target<=40,"checkpoint driver: prospective terminal step bound");}
    BoundedGenerationStore store(store_path,trusted);DurableDriverCursor cursor;
    try
    {
        std::unique_ptr<BoundedRunningContract> bound;
        if(fresh)
        {
            require(!store.currentPointerBytes(),"checkpoint driver: fresh store already selected");
            const auto initial=decodeBoundedSessionCandidate(readAt(root,"initial_candidate.bin"),1,trusted);
            require(initial.body.accepted.gel.steps==0,"checkpoint driver: fresh initial step");
            bound.reset(new BoundedRunningContract(running(root,initial.contract)));
            VerifiedBoundedMechanicalPrefix proof(root+"/initial.phy",initial,1,trusted,true);
            cursor.publication_stage=1;
            const auto result=store.publish(initial,1,proof,"g_"+label+"_initial",std::nullopt);
            cursor.publication_error_number=result.error_number;cursor.accepted(store.live());
            require(result.status==BoundedPublishStatus::PublishedDurable,"checkpoint driver: initial publication durability uncertain");
            cursor.accepted(store.live());
            // Reusing the original full header preserves independent reference
            // compatibility. Only this fresh branch emits a t0 header.
            const auto input=obs::readInput((root+"/original.input").c_str());
            // Header needs stored plateIDs; no fresh reference construction.
            const auto &body=store.live()->candidate().body;const BoundedStoredSceneView view(body);
            obs::header(input,view,body.accepted,body.plate_ids);
        }
        else
        {
            expected=store.currentPointerBytes();require(bool(expected),"checkpoint driver: missing expected HEAD");
            require(sha256(*expected,expected->size())==unhex(argv[7]),"checkpoint driver: expected HEAD digest mismatch");
            const auto selection=decodeBoundedSelection(selected,generation);const auto loaded=store.openExplicit(selection);
            bound.reset(new BoundedRunningContract(running(root,loaded.candidate.contract)));
            require(validate || loaded.candidate.body.accepted.gel.steps<target,"checkpoint driver: resumed target must advance");
            if(validate)
            {
                // No private branch, live owner or physical output mutation.
                require(store.currentPointerBytes()==expected,"checkpoint driver: HEAD changed during noadvance validation");
                io::begin("validated_noadvance");io::field("accepted_step",loaded.candidate.body.accepted.gel.steps);
                io::field("archive_sha256",hex(sha256(loaded.archive_bytes,loaded.archive_bytes.size())));
                io::field("new_physical_steps",0);io::field("full_restart_qualified",false);io::finish();return 0;
            }
            cursor.publication_stage=2;
            const auto result=store.publish(loaded.candidate,generation,*loaded.prefix,"g_"+label+"_resume",expected);
            cursor.publication_error_number=result.error_number;cursor.accepted(store.live());
            require(result.status==BoundedPublishStatus::PublishedDurable,"checkpoint driver: resume publication durability uncertain");
            cursor.accepted(store.live());
            io::begin("resume_header");io::field("selected_generation",generation);io::field("accepted_step",cursor.durable_step);
            io::field("accepted_time",cursor.durable_time);io::field("new_initial_state",false);io::field("contract_sha256",hex(trusted));io::finish();
        }
        require(target<=store.live()->candidate().contract.protocol.steps(),"checkpoint driver: terminal protocol limit");
        while(store.live()->candidate().body.accepted.gel.steps<target)
        {
            // OWN the old packet. append/publish swaps and destroys the previous
            // owner; no borrowed scene/state may dangle in subsequent debug IO.
            const auto old=store.live()->candidate().body;const auto contract=store.live()->candidate().contract;
            bound->requireFor(contract);const BoundedStoredSceneView scene(old);
            const auto interval=contract.protocol.interval(old.accepted.gel.steps,old.accepted.gel.time);
            const auto rows=boundedGelFrictionContactsAtOffsets(scene,interval.old_offset,interval.new_offset,old.mu,interval.wall_tangent_velocity);
            ++cursor.solver_calls;const auto result=advanceBoundedGelFriction(old.model,old.accepted,interval.actual_dt,rows,old.plate_ids,contract.solver);
            std::array<Real,2> whole{0,0},endpoint{0,0};
            if(result.accepted)
            {
                whole=boundedFrictionIntervalGaps(scene,old.accepted.gel.position,result.next.gel.position,interval);
                endpoint=boundedFrictionIntervalGaps(scene,old.accepted.gel.position,result.next.gel.position,interval,false);
            }
            const bool geometry=result.accepted && std::min(whole[0],whole[1])>=-contract.penetration_stroke_gate*contract.protocol.normal_stroke;
            if(!geometry)
            {
                obs::step(scene,old.accepted,result,interval,rows,false,whole,endpoint,result.accepted?"whole_interval_dense_penetration":result.status.c_str());
                throw std::runtime_error(result.accepted?"checkpoint driver: whole interval geometry gate":result.status);
            }
            BoundedMechanicalRecord record;record.body=old;record.body.accepted=result.next;record.kind=1;
            record.interval_phase=interval.phase;record.old_time=interval.old_time;
            record.observation=boundedMechanicalObservation(result,interval.actual_dt,whole,endpoint);
            const auto append=BoundedAcceptedWriter::append(store,record,*bound);
            cursor.accepted(store.live()); // primitive cursor FIRST, before JSON/allocation
            require(append.status==BoundedAppendStatus::AcceptedDurable,"checkpoint driver: accepted physical write failed");
            obs::step(scene,old.accepted,result,interval,rows,true,whole,endpoint,"all_applicable_step_gates_passed");
            cursor.debug_step=cursor.durable_step;cursor.debug_time=cursor.durable_time;
            const auto k=cursor.durable_step;if(k==8 || k==20 || k==32 || k==40 || k==target)seal(store,store_path,trusted,label,cursor);
        }
        io::begin("end");io::field("completed",true);io::field("accepted_steps",cursor.durable_step);
        io::field("last_accepted_time",cursor.durable_time);io::field("actual_solver_calls",cursor.solver_calls);
        io::field("full_restart_qualified",false);io::field("S3_complete",false);io::finish();return 0;
    }
    catch(const std::exception &error)
    {
        cursor.accepted(store.live());diagnostic(error,cursor,store.live());return 1;
    }
}
}
int main(int argc,char **argv)
{
    std::ios::sync_with_stdio(false);std::cout<<std::setprecision(17);
    std::cout.exceptions(std::ios::badbit|std::ios::failbit);
    try
    {
        require(argc>=2,"checkpoint driver: missing mode");const std::string mode=argv[1];
        if(mode=="--capabilities")
        {
            require(argc==2,"checkpoint driver: capability arguments");io::begin("capabilities");
            io::field("actual_PD_checkpoint_driver",true);io::field("initialization_only_preflight_available",true);
            io::field("explicit_selected_generation",true);io::field("full_restart_qualified",false);io::field("S3_complete",false);io::finish();return 0;
        }
        if(mode=="--diagnostic-fixture")
        {
            require(argc==2,"checkpoint driver: diagnostic fixture arguments");
            DurableDriverCursor fixture;fixture.durable_step=7;fixture.durable_time=.07;fixture.durable_bytes=123456;
            fixture.debug_step=6;fixture.debug_time=.06;fixture.solver_calls=7;
            fixture.publication_error_number=EIO;fixture.publication_stage=3;
            diagnostic(std::runtime_error("controlled formatter fixture: quote\" backslash\\ newline\n"),fixture,nullptr);
            io::begin("diagnostic_fixture");io::field("controlled_formatter_only",true);io::field("new_physical_steps",0);io::finish();return 0;
        }
        io::identity();
        if(mode=="--initialize"){require(argc==3,"checkpoint driver: initialize ROOT");initialize(argv[2]);return 0;}
        return execute(argc,argv);
    }
    catch(const std::exception &error){diagnostic(error,DurableDriverCursor{},nullptr);return 1;}
}
