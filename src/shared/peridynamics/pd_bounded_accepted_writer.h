/* Licensed under the Apache License, Version 2.0. */
#ifndef PD_BOUNDED_ACCEPTED_WRITER_H
#define PD_BOUNDED_ACCEPTED_WRITER_H
#include "pd_bounded_generation_store.h"
#include <type_traits>

namespace SPH::peridynamics
{
/** Bind the ACTUAL running executable and frozen producer manifests. A tool
 * which can inspect a different producer's archive does not thereby gain
 * permission to advance or append accepted records under that identity.
 * The driver also verifies all actual dependencies named by its manifests.
 */
class BoundedRunningContract
{
    canonical_archive::Digest contract_digest_;
    std::unique_ptr<canonical_archive::Descriptor> executable_,source_,build_,input_;
    struct stat executable_stamp_{},source_stamp_{},build_stamp_{},input_stamp_{};
    static int openManifest(const std::string &path)
    {
        bounded_session_detail::require(!path.empty() && path[0]=='/',"running contract: absolute manifest path required");
        const int fd=::open(path.c_str(),O_RDONLY|O_NOFOLLOW|O_NONBLOCK|O_CLOEXEC);
        if(fd<0)bounded_generation_detail::fail("running manifest open");return fd;
    }
    static struct stat verify(int fd,uint64_t bound,const canonical_archive::Digest &expected,const char *reason)
    {
        const auto before=bounded_generation_detail::regular(fd,bound);
        bounded_session_detail::require(bounded_generation_detail::hashFile(fd,uint64_t(before.st_size))==expected,reason);return before;
    }
  public:
    BoundedRunningContract(const BoundedSessionContract &contract,const std::string &source_manifest,
        const std::string &build_runtime_manifest,const std::string &original_input)
        :contract_digest_(boundedContractDigest(contract)),
         executable_(new canonical_archive::Descriptor(::open("/proc/self/exe",O_RDONLY|O_CLOEXEC))),
         source_(new canonical_archive::Descriptor(openManifest(source_manifest))),
         build_(new canonical_archive::Descriptor(openManifest(build_runtime_manifest))),
         input_(new canonical_archive::Descriptor(openManifest(original_input)))
    {
        using namespace bounded_session_detail;contract.validateProfile();
        executable_stamp_=verify(executable_->get(),256*1024*1024,digestArray(contract.ELF_sha256),"running contract: actual executable digest mismatch");
        source_stamp_=verify(source_->get(),canonical_archive::maximum_bytes,digestArray(contract.source_manifest_sha256),"running contract: source manifest digest mismatch");
        build_stamp_=verify(build_->get(),canonical_archive::maximum_bytes,digestArray(contract.build_runtime_manifest_sha256),"running contract: build runtime manifest digest mismatch");
        input_stamp_=verify(input_->get(),canonical_archive::maximum_bytes,digestArray(contract.input_sha256),"running contract: original input digest mismatch");
    }
    void requireFor(const BoundedSessionContract &contract)const
    {
        contract.validateProfile();bounded_session_detail::require(boundedContractDigest(contract)==contract_digest_,"running contract: bound contract mismatch");
        bounded_generation_detail::unchanged(executable_->get(),executable_stamp_,256*1024*1024);
        bounded_generation_detail::unchanged(source_->get(),source_stamp_,canonical_archive::maximum_bytes);
        bounded_generation_detail::unchanged(build_->get(),build_stamp_,canonical_archive::maximum_bytes);
        bounded_generation_detail::unchanged(input_->get(),input_stamp_,canonical_archive::maximum_bytes);
    }
};

enum class BoundedAppendFault{None,BeforeAppend,AfterPartialWrite,AfterFullWriteBeforeSync,FsyncError};
enum class BoundedAppendStatus{AcceptedDurable,WriteFailed,ReadbackFailed,SyncFailed};
struct BoundedAppendResult{BoundedAppendStatus status;BoundedAppendWitness witness;};

/** One serial cooperating accepted-output writer. No Newton is called here.
 * All proposal verification/allocation/encoding/hash preparation happens
 * before any append. Successful fsync is followed only by noexcept state/
 * cursor ownership swap and primitive witness updates.
 */
class BoundedAcceptedWriter
{
    using Bytes=canonical_archive::Bytes;using Digest=canonical_archive::Digest;
    struct Prepared
    {
        std::unique_ptr<BoundedSessionCandidate> next;
        Bytes frame;
        struct stat before{};
    };
    static Prepared prepare(BoundedOwnedCheckpoint &owner,const BoundedMechanicalRecord &record,
        const BoundedRunningContract &running,const Digest &trusted)
    {
        using namespace canonical_archive;using namespace bounded_session_detail;
        const auto &old=*owner.candidate_;running.requireFor(old.contract);old.validateCandidate(trusted);
        record.body.validateStored();bounded_mechanical_detail::validateStatic(record.body,old.contract);
        BoundedMechanicalRecord previous;previous.body=old.body;
        bounded_mechanical_detail::validateAdjacent(previous,record,old.contract);
        const auto encoded=encodeBoundedMechanicalRecord(record);Prepared prepared;
        appendU64(prepared.frame,encoded.size());prepared.frame.insert(prepared.frame.end(),encoded.begin(),encoded.end());
        require(old.output.bytes<=maximum_output_prefix && prepared.frame.size()<=maximum_output_prefix-old.output.bytes,
            "accepted writer: output byte budget");
        prepared.next.reset(new BoundedSessionCandidate(old));auto &next=*prepared.next;
        next.body=record.body;next.output.records=old.output.records+1;next.output.bytes=old.output.bytes+prepared.frame.size();
        next.output.last_step=record.body.accepted.gel.steps;next.output.last_time=record.body.accepted.gel.time;
        const auto step=next.body.accepted.gel.steps;
        next.next_phase=step==next.contract.protocol.steps()?3:(step<16?0:(step<24?1:2));
        prepared.before=bounded_generation_detail::regular(owner.working_->get(),maximum_output_prefix);
        require(uint64_t(prepared.before.st_size)==old.output.bytes,"accepted writer: exact committed working length");
        Sha256Stream original_hash,combined_hash;std::array<uint8_t,65536> buffer{};uint64_t offset=0;
        while(offset<old.output.bytes)
        {
            const size_t wanted=size_t(std::min(uint64_t(buffer.size()),old.output.bytes-offset));
            const auto n=::pread(owner.working_->get(),buffer.data(),wanted,off_t(offset));
            if(n<0 && errno==EINTR)continue;if(n<=0)bounded_generation_detail::fail("old committed prefix read");
            original_hash.update(buffer.data(),size_t(n));combined_hash.update(buffer.data(),size_t(n));offset+=uint64_t(n);
        }
        require(original_hash.finish()==digestArray(old.output.prefix_sha256),"accepted writer: old committed prefix digest mismatch");
        combined_hash.update(prepared.frame.data(),prepared.frame.size());const auto digest=combined_hash.finish();
        next.output.prefix_sha256.assign(digest.begin(),digest.end());
        next.validateCandidate(trusted);bounded_generation_detail::unchanged(owner.working_->get(),prepared.before,maximum_output_prefix);
        return prepared;
    }
  public:
    [[nodiscard]] static BoundedAppendResult append(BoundedGenerationStore &store,const BoundedMechanicalRecord &record,
        const BoundedRunningContract &running,BoundedAppendFault fault=BoundedAppendFault::None)
    {
        using namespace bounded_session_detail;using namespace canonical_archive;
        require(bool(store.live_),"accepted writer: no owned checkpoint");auto &owner=*store.live_;
        require(!owner.publication_uncertain_,"accepted writer: publication uncertain");
        require(!owner.writer_stopped_,"accepted writer: stopped; explicit recovery required");
        const auto &old=*owner.candidate_;owner.append_witness_={};auto &w=owner.append_witness_;
        w.stage=BoundedAppendStage::Preparing;w.old_committed_bytes=old.output.bytes;w.old_step=old.body.accepted.gel.steps;
        w.old_time=old.body.accepted.gel.time;w.proposed_step=record.body.accepted.gel.steps;w.proposed_time=record.body.accepted.gel.time;
        Prepared prepared;
        try
        {
            prepared=prepare(owner,record,running,store.trusted_);w.frame_bytes=prepared.frame.size();
            if(fault==BoundedAppendFault::BeforeAppend){w.injected=true;throw std::runtime_error("accepted writer: injected before append");}
            bounded_generation_detail::unchanged(owner.working_->get(),prepared.before,maximum_output_prefix);
        }
        catch(...){owner.writer_stopped_=true;w.failure_stage=w.stage;w.stage=BoundedAppendStage::Stopped;throw;}
        // No user-space allocation/open/encoding/hash follows this point.
        owner.writer_stopped_=true;w.stage=BoundedAppendStage::Writing;
        const size_t target=fault==BoundedAppendFault::AfterPartialWrite?std::max(size_t(1),prepared.frame.size()/2):prepared.frame.size();
        auto failed=[&](BoundedAppendStatus status,int error,bool injected)noexcept
        {w.failure_stage=w.stage;w.stage=BoundedAppendStage::Stopped;w.error_number=error;w.injected=injected;return BoundedAppendResult{status,w};};
        while(w.written_bytes<target)
        {
            const auto n=::write(owner.working_->get(),prepared.frame.data()+w.written_bytes,target-size_t(w.written_bytes));
            if(n<0 && errno==EINTR)continue;if(n<=0)return failed(BoundedAppendStatus::WriteFailed,n<0?errno:EIO,false);
            w.written_bytes+=uint64_t(n);
        }
        if(fault==BoundedAppendFault::AfterPartialWrite || fault==BoundedAppendFault::AfterFullWriteBeforeSync)
            return failed(BoundedAppendStatus::WriteFailed,EIO,true);
        w.stage=BoundedAppendStage::Readback;struct stat after{};
        if(::fstat(owner.working_->get(),&after)!=0)return failed(BoundedAppendStatus::ReadbackFailed,errno,false);
        if(after.st_dev!=prepared.before.st_dev || after.st_ino!=prepared.before.st_ino || after.st_size<0 ||
            uint64_t(after.st_size)!=w.old_committed_bytes+w.frame_bytes)
            return failed(BoundedAppendStatus::ReadbackFailed,EIO,false);
        std::array<uint8_t,65536> buffer{};uint64_t offset=0;
        while(offset<w.frame_bytes)
        {
            const size_t wanted=size_t(std::min(uint64_t(buffer.size()),w.frame_bytes-offset));
            const auto n=::pread(owner.working_->get(),buffer.data(),wanted,off_t(w.old_committed_bytes+offset));
            if(n<0 && errno==EINTR)continue;if(n<=0)return failed(BoundedAppendStatus::ReadbackFailed,n<0?errno:EIO,false);
            if(std::memcmp(buffer.data(),prepared.frame.data()+offset,size_t(n))!=0)return failed(BoundedAppendStatus::ReadbackFailed,EIO,false);
            offset+=uint64_t(n);
        }
        w.stage=BoundedAppendStage::Syncing;
        if(fault==BoundedAppendFault::FsyncError)return failed(BoundedAppendStatus::SyncFailed,EIO,true);
        if(::fsync(owner.working_->get())!=0)return failed(BoundedAppendStatus::SyncFailed,errno,false);
        static_assert(std::is_nothrow_destructible_v<BoundedSessionCandidate>,"accepted state cleanup must not throw");
        owner.candidate_.swap(prepared.next);owner.writer_stopped_=false;w.stage=BoundedAppendStage::Accepted;w.committed=true;
        return {BoundedAppendStatus::AcceptedDurable,w};
    }
};
}
#endif
