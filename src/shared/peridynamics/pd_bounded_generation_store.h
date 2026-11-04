/* Licensed under the Apache License, Version 2.0. */
#ifndef PD_BOUNDED_GENERATION_STORE_H
#define PD_BOUNDED_GENERATION_STORE_H
#include "pd_bounded_mechanical_stream.h"
#include <memory>
#include <sys/file.h>

namespace SPH::peridynamics
{
namespace bounded_generation_detail
{
using namespace canonical_archive;
inline void fail(const char *operation)
{const int error=errno;throw std::runtime_error(std::string("generation store: ")+operation+" errno="+std::to_string(error));}
inline void leaf(const std::string &name)
{
    bounded_session_detail::require(!name.empty() && name.size()<=128 && name!="." && name!="..",
        "generation store: unsafe file leaf");
    for(unsigned char c:name)bounded_session_detail::require(
        (c>='a'&&c<='z')||(c>='A'&&c<='Z')||(c>='0'&&c<='9')||c=='_'||c=='-'||c=='.',
        "generation store: unsafe file leaf");
}
inline void branchLeaf(const std::string &name)
{
    bounded_session_detail::require(name.size()>=3 && name.size()<=96 && name.substr(0,2)=="g_",
        "generation store: unsafe generation branch");
    for(unsigned char c:name)bounded_session_detail::require(
        (c>='a'&&c<='z')||(c>='A'&&c<='Z')||(c>='0'&&c<='9')||c=='_'||c=='-',
        "generation store: unsafe generation branch");
}
inline struct stat regular(int fd,uint64_t bound)
{
    struct stat s{};
    if(::fstat(fd,&s)!=0)fail("file stat");
    bounded_session_detail::require(S_ISREG(s.st_mode) && s.st_nlink==1 && s.st_size>=0 && uint64_t(s.st_size)<=bound,
        "generation store: regular file profile");return s;
}
inline void unchanged(int fd,const struct stat &before,uint64_t bound)
{
    const auto after=regular(fd,bound);
    bounded_session_detail::require(before.st_dev==after.st_dev && before.st_ino==after.st_ino &&
        before.st_size==after.st_size && before.st_mtim.tv_sec==after.st_mtim.tv_sec && before.st_mtim.tv_nsec==after.st_mtim.tv_nsec &&
        before.st_ctim.tv_sec==after.st_ctim.tv_sec && before.st_ctim.tv_nsec==after.st_ctim.tv_nsec,
        "generation store: file changed while reading");
}
inline void sync(int fd,const char *why){if(::fsync(fd)!=0)fail(why);}
inline int openDirectory(const std::string &path)
{
    bounded_session_detail::require(!path.empty() && path[0]=='/',"generation store: absolute root required");
    const int fd=::open(path.c_str(),O_RDONLY|O_DIRECTORY|O_NOFOLLOW|O_CLOEXEC);
    if(fd<0)fail("root directory open");return fd;
}
inline int openBranch(int root,const std::string &name)
{
    branchLeaf(name);const int fd=::openat(root,name.c_str(),O_RDONLY|O_DIRECTORY|O_NOFOLLOW|O_CLOEXEC);
    if(fd<0)fail("branch directory open");return fd;
}
inline std::optional<Bytes> readAt(int directory,const std::string &name,uint64_t bound,bool allow_missing=false)
{
    leaf(name);const int raw=::openat(directory,name.c_str(),O_RDONLY|O_NOFOLLOW|O_CLOEXEC|O_NONBLOCK);
    if(raw<0)
    {
        if(allow_missing && errno==ENOENT)return std::nullopt;
        fail("selected file open");
    }
    Descriptor file(raw);const auto before=regular(file.get(),bound);Bytes bytes(size_t(before.st_size),0);size_t offset=0;
    while(offset<bytes.size())
    {
        const auto n=::pread(file.get(),bytes.data()+offset,bytes.size()-offset,off_t(offset));
        if(n<0 && errno==EINTR)continue;if(n<=0)fail("selected file read");offset+=size_t(n);
    }
    unchanged(file.get(),before,bound);return bytes;
}
inline void writeExclusive(int directory,const std::string &name,const Bytes &bytes)
{
    leaf(name);const int raw=::openat(directory,name.c_str(),O_WRONLY|O_CREAT|O_EXCL|O_NOFOLLOW|O_CLOEXEC,0600);
    if(raw<0)fail("exclusive file create");Descriptor file(raw);size_t offset=0;
    while(offset<bytes.size())
    {
        const auto n=::write(file.get(),bytes.data()+offset,bytes.size()-offset);
        if(n<0 && errno==EINTR)continue;if(n<=0)fail("exclusive file write");offset+=size_t(n);
    }
    bounded_session_detail::require(uint64_t(regular(file.get(),maximum_output_prefix).st_size)==bytes.size(),
        "generation store: exclusive file length");sync(file.get(),"exclusive file fsync");
}
inline Digest hashFile(int fd,uint64_t expected_bytes)
{
    const auto before=regular(fd,maximum_output_prefix);
    bounded_session_detail::require(uint64_t(before.st_size)==expected_bytes,"generation store: working copy length");
    Sha256Stream hash;std::array<uint8_t,65536> buffer{};uint64_t offset=0;
    while(offset<expected_bytes)
    {
        const auto wanted=size_t(std::min(uint64_t(buffer.size()),expected_bytes-offset));
        const auto n=::pread(fd,buffer.data(),wanted,off_t(offset));
        if(n<0 && errno==EINTR)continue;if(n<=0)fail("working copy read");hash.update(buffer.data(),size_t(n));offset+=uint64_t(n);
    }
    unchanged(fd,before,maximum_output_prefix);return hash.finish();
}
}

struct BoundedGenerationSelection
{
    uint64_t generation=0;
    std::string branch;
    canonical_archive::Digest archive_sha256{},physical_prefix_sha256{},contract_sha256{};
};
inline canonical_archive::Bytes encodeBoundedSelection(const BoundedGenerationSelection &selection)
{
    using namespace canonical_archive;bounded_generation_detail::branchLeaf(selection.branch);Builder builder(selection.generation);
    builder.integers("SelectionSchema",{},{1});builder.bytes("Branch",{selection.branch.size()},Bytes(selection.branch.begin(),selection.branch.end()));
    builder.bytes("ArchiveSHA256",{32},Bytes(selection.archive_sha256.begin(),selection.archive_sha256.end()));
    builder.bytes("PhysicalPrefixSHA256",{32},Bytes(selection.physical_prefix_sha256.begin(),selection.physical_prefix_sha256.end()));
    builder.bytes("ContractSHA256",{32},Bytes(selection.contract_sha256.begin(),selection.contract_sha256.end()));return builder.finish();
}
inline BoundedGenerationSelection decodeBoundedSelection(const canonical_archive::Bytes &bytes,uint64_t explicit_generation)
{
    using namespace canonical_archive;Candidate candidate(bytes,explicit_generation);
    candidate.requireExactNames({"SelectionSchema","Branch","ArchiveSHA256","PhysicalPrefixSHA256","ContractSHA256"});
    bounded_session_detail::require(candidate.integers("SelectionSchema",{}).at(0)==1,"generation store: selection schema");
    BoundedGenerationSelection selection;selection.generation=explicit_generation;
    // Bound the explicit byte-array shape before reading; never infer a path
    // from a latest pointer or search for another generation on failure.
    const auto &shape=candidate.dimensions("Branch");
    bounded_session_detail::require(shape.size()==1 && shape[0]>=3 && shape[0]<=96,"generation store: selected branch dimensions");
    const auto branch=candidate.bytes("Branch",shape);selection.branch=std::string(branch.begin(),branch.end());
    bounded_generation_detail::branchLeaf(selection.branch);
    auto digest=[&](const char *name){const auto data=candidate.bytes(name,{32});Digest out;std::copy(data.begin(),data.end(),out.begin());return out;};
    selection.archive_sha256=digest("ArchiveSHA256");selection.physical_prefix_sha256=digest("PhysicalPrefixSHA256");
    selection.contract_sha256=digest("ContractSHA256");return selection;
}

class BoundedGenerationStore;
/** Owned checkpoint state/output handles. No physical advance API yet: the
 * eventual native driver must separately bind its running ELF/contract.
 */
class BoundedOwnedCheckpoint
{
    friend class BoundedGenerationStore;
    BoundedSessionCandidate candidate_;
    uint64_t generation_;
    std::string branch_;
    std::shared_ptr<canonical_archive::Descriptor> root_,lock_;
    std::unique_ptr<canonical_archive::Descriptor> directory_,working_;
    bool publication_uncertain_=false;
    BoundedOwnedCheckpoint(const BoundedSessionCandidate &candidate,uint64_t generation,std::string branch,
        std::shared_ptr<canonical_archive::Descriptor> root,std::shared_ptr<canonical_archive::Descriptor> lock,int directory)
        :candidate_(candidate),generation_(generation),branch_(std::move(branch)),root_(std::move(root)),lock_(std::move(lock)),
         directory_(new canonical_archive::Descriptor(::fcntl(directory,F_DUPFD_CLOEXEC,0))),
         working_(new canonical_archive::Descriptor(::openat(directory_->get(),"active.phy",O_RDWR|O_APPEND|O_NOFOLLOW|O_CLOEXEC)))
    {
        bounded_session_detail::require(bounded_generation_detail::hashFile(working_->get(),candidate.output.bytes)==
            bounded_session_detail::digestArray(candidate.output.prefix_sha256),"generation store: working copy digest");
    }
  public:
    const BoundedSessionCandidate &candidate()const noexcept{return candidate_;}
    uint64_t generation()const noexcept{return generation_;}
    const std::string &branch()const noexcept{return branch_;}
    bool publicationUncertain()const noexcept{return publication_uncertain_;}
};

struct BoundedSelectedCheckpoint
{
    BoundedSessionCandidate candidate;
    canonical_archive::Bytes archive_bytes;
    uint64_t generation;
    std::unique_ptr<VerifiedBoundedMechanicalPrefix> prefix;
};
enum class BoundedPublishFault{None,BeforeStage,AfterFiles,BeforeRename,AfterRenameDirectorySync};
enum class BoundedPublishStatus{PublishedDurable,PublishedDirectoryDurabilityUncertain};
struct BoundedPublishResult
{
    BoundedPublishStatus status;
    int error_number=0;
    bool injected=false;
};

/** Cooperating single-owner store, explicit selection only. All branch and
 * file operations are directory-FD relative. Inputs/root ancestors and the
 * external trusted contract are chosen by the frozen driver, not by archives.
 */
class BoundedGenerationStore
{
    using Descriptor=canonical_archive::Descriptor;
    std::shared_ptr<Descriptor> root_,lock_;
    canonical_archive::Digest trusted_;
    std::unique_ptr<BoundedOwnedCheckpoint> live_;
  public:
    BoundedGenerationStore(const BoundedGenerationStore &)=delete;
    BoundedGenerationStore &operator=(const BoundedGenerationStore &)=delete;
    BoundedGenerationStore(BoundedGenerationStore &&)=delete;
    BoundedGenerationStore &operator=(BoundedGenerationStore &&)=delete;
    BoundedGenerationStore(const std::string &absolute_existing_root,const canonical_archive::Digest &external_trusted_contract)
        :root_(new Descriptor(bounded_generation_detail::openDirectory(absolute_existing_root))),
         lock_(new Descriptor(::openat(root_->get(),"SESSION.lock",O_RDWR|O_CREAT|O_NOFOLLOW|O_CLOEXEC|O_NONBLOCK,0600))),trusted_(external_trusted_contract)
    {
        bounded_session_detail::require(bounded_generation_detail::regular(lock_->get(),0).st_size==0,"generation store: lock file profile");
        if(::flock(lock_->get(),LOCK_EX|LOCK_NB)!=0)bounded_generation_detail::fail("exclusive owner lock");
        bounded_generation_detail::sync(lock_->get(),"lock file fsync");
    }
    std::optional<canonical_archive::Bytes> currentPointerBytes()const
    {return bounded_generation_detail::readAt(root_->get(),"HEAD",2048,true);}
    const BoundedOwnedCheckpoint *live()const noexcept{return live_.get();}

    BoundedSelectedCheckpoint openExplicit(const BoundedGenerationSelection &selection)const
    {
        using namespace canonical_archive;using namespace bounded_generation_detail;
        bounded_session_detail::require(selection.contract_sha256==trusted_,"generation store: external contract selection mismatch");
        Descriptor directory(openBranch(root_->get(),selection.branch));const auto archive=readAt(directory.get(),"state.bin",maximum_bytes).value();
        bounded_session_detail::require(sha256(archive,archive.size())==selection.archive_sha256,"generation store: selected archive digest");
        auto candidate=decodeBoundedSessionCandidate(archive,selection.generation,trusted_);
        bounded_session_detail::require(bounded_session_detail::digestArray(candidate.output.prefix_sha256)==selection.physical_prefix_sha256,
            "generation store: selected physical prefix identity");
        const auto path=std::string("/proc/self/fd/")+std::to_string(directory.get())+"/prefix.phy";
        auto proof=std::unique_ptr<VerifiedBoundedMechanicalPrefix>(new VerifiedBoundedMechanicalPrefix(path,candidate,selection.generation,trusted_,true));
        return {std::move(candidate),archive,selection.generation,std::move(proof)};
    }

    [[nodiscard]] BoundedPublishResult publish(const BoundedSessionCandidate &candidate,uint64_t generation,
        const VerifiedBoundedMechanicalPrefix &proof,const std::string &new_branch,
        const std::optional<canonical_archive::Bytes> &expected_current,
        BoundedPublishFault fault=BoundedPublishFault::None)
    {
        using namespace canonical_archive;using namespace bounded_generation_detail;
        bounded_session_detail::require(!live_ || !live_->publication_uncertain_,"generation store: publication uncertain; explicit recovery required");
        branchLeaf(new_branch);candidate.validateCandidate(trusted_);proof.requireBoundCandidate(candidate,generation);
        bounded_session_detail::require(currentPointerBytes()==expected_current,"generation store: expected-current pointer mismatch");
        if(fault==BoundedPublishFault::BeforeStage)throw std::runtime_error("generation store: injected before stage");
        const auto archive=encodeBoundedSessionCandidate(candidate,generation);
        BoundedGenerationSelection selection{generation,new_branch,sha256(archive,archive.size()),
            bounded_session_detail::digestArray(candidate.output.prefix_sha256),trusted_};
        const auto pointer=encodeBoundedSelection(selection);
        const auto pending=std::string("pending_")+new_branch+".bin",previous=std::string("previous_")+new_branch+".bin";
        if(::mkdirat(root_->get(),new_branch.c_str(),0700)!=0)fail("exclusive branch mkdir");
        Descriptor directory(openBranch(root_->get(),new_branch));
        writeExclusive(directory.get(),"state.bin",archive);
        bounded_session_detail::require(readAt(directory.get(),"state.bin",maximum_bytes).value()==archive,"generation store: stored archive bytes");
        proof.stageExclusiveCopy(candidate,generation,directory.get(),"prefix.phy");
        proof.stageExclusiveCopy(candidate,generation,directory.get(),"active.phy");
        auto prepared=std::unique_ptr<BoundedOwnedCheckpoint>(new BoundedOwnedCheckpoint(candidate,generation,new_branch,root_,lock_,directory.get()));
        sync(directory.get(),"branch directory fsync");sync(root_->get(),"new branch parent fsync");
        if(fault==BoundedPublishFault::AfterFiles)throw std::runtime_error("generation store: injected after private files");
        if(expected_current)writeExclusive(root_->get(),previous,*expected_current);
        writeExclusive(root_->get(),pending,pointer);sync(root_->get(),"pending pointer parent fsync");
        bounded_session_detail::require(currentPointerBytes()==expected_current,"generation store: expected-current changed before publish");
        if(fault==BoundedPublishFault::BeforeRename)throw std::runtime_error("generation store: injected before pointer rename");
        if(::renameat(root_->get(),pending.c_str(),root_->get(),"HEAD")!=0)fail("atomic pointer rename");
        // Everything which may allocate/open/validate is BEFORE rename. After
        // a successful rename, even a durability error selects the new owner.
        // Return an explicit uncertainty result, never claim live unchanged.
        const bool injected=fault==BoundedPublishFault::AfterRenameDirectorySync;
        const int code=injected ? -1 : ::fsync(root_->get());const int error=code==0?0:(injected?EIO:errno);
        prepared->publication_uncertain_=code!=0;live_.swap(prepared); // unique_ptr swap is noexcept
        return {code==0?BoundedPublishStatus::PublishedDurable:BoundedPublishStatus::PublishedDirectoryDurabilityUncertain,error,injected};
    }
};
}
#endif
