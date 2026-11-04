/* Licensed under the Apache License, Version 2.0. */
#ifndef PD_OUTPUT_PREFIX_H
#define PD_OUTPUT_PREFIX_H
#include "pd_sha256.h"
#include <fcntl.h>
#include <string>
#include <sys/stat.h>

namespace SPH::peridynamics::canonical_archive
{
inline constexpr uint64_t maximum_output_prefix=uint64_t(4)*1024*1024*1024;
static_assert(sizeof(off_t)>=8,"Explicit output byte cursors require 64-bit file offsets");

/** Read-only byte proof, not a mechanical/output-cursor semantic validator.
 * A live file may contain a later tail. The explicit selected byte cursor is
 * hashed, never selected by a LATEST pointer or inferred from that tail.
 * A separate immutable snapshot must have exactly the selected byte count.
 */
class VerifiedOutputPrefix
{
    Descriptor source_;
    struct stat before_{};
    uint64_t count_;
    Digest digest_;
    static struct stat regular(int fd)
    {
        struct stat result{};
        if(::fstat(fd,&result)!=0 || !S_ISREG(result.st_mode) || result.st_size<0)
            throw std::invalid_argument("output prefix: not a readable regular file");
        return result;
    }
    void unchanged()const
    {
        const auto after=regular(source_.get());
        if(before_.st_dev!=after.st_dev || before_.st_ino!=after.st_ino || before_.st_size!=after.st_size ||
            before_.st_mtim.tv_sec!=after.st_mtim.tv_sec || before_.st_mtim.tv_nsec!=after.st_mtim.tv_nsec ||
            before_.st_ctim.tv_sec!=after.st_ctim.tv_sec || before_.st_ctim.tv_nsec!=after.st_ctim.tv_nsec)
            throw std::runtime_error("output prefix: source changed during validation or staging");
    }
    Digest transfer(int destination)const
    {
        std::array<uint8_t,65536> buffer{};uint64_t offset=0;Sha256Stream hash;
        while(offset<count_)
        {
            const size_t wanted=size_t(std::min(uint64_t(buffer.size()),count_-offset));
            ssize_t received;
            do{received=::pread(source_.get(),buffer.data(),wanted,off_t(offset));}while(received<0 && errno==EINTR);
            if(received<=0)throw std::runtime_error("output prefix: premature read end");
            hash.update(buffer.data(),size_t(received));
            if(destination>=0)
            {
                size_t sent=0;
                while(sent<size_t(received))
                {
                    const ssize_t n=::write(destination,buffer.data()+sent,size_t(received)-sent);
                    if(n<0 && errno==EINTR)continue;
                    if(n<=0)throw std::runtime_error("output prefix: private staging write");
                    sent+=size_t(n);
                }
            }
            offset+=uint64_t(received);
        }
        unchanged();return hash.finish();
    }
  public:
    VerifiedOutputPrefix(const std::string &path,uint64_t explicit_bytes,const Digest &expected,bool exact_snapshot)
        :source_(::open(path.c_str(),O_RDONLY|O_CLOEXEC|O_NOFOLLOW|O_NONBLOCK)),count_(explicit_bytes),digest_(expected)
    {
        before_=regular(source_.get());
        if(count_>maximum_output_prefix || count_>uint64_t(before_.st_size) ||
            (exact_snapshot && count_!=uint64_t(before_.st_size)))
            throw std::invalid_argument("output prefix: explicit byte cursor or snapshot size");
        if(transfer(-1)!=digest_)throw std::invalid_argument("output prefix: selected digest mismatch");
    }
    uint64_t bytes()const{return count_;}
    const Digest &digest()const{return digest_;}
    void assertUnchanged()const{unchanged();}

    // Read semantic records through the SAME verified descriptor. Reopening a
    // pathname after hash validation would introduce a path-replacement race.
    void readAt(uint64_t offset,uint8_t *destination,size_t size)const
    {
        if(offset>count_ || uint64_t(size)>count_-offset || (size && !destination))
            throw std::invalid_argument("output prefix: semantic read bound");
        unchanged();size_t received=0;
        while(received<size)
        {
            const auto n=::pread(source_.get(),destination+received,size-received,off_t(offset+received));
            if(n<0 && errno==EINTR)continue;
            if(n<=0)throw std::runtime_error("output prefix: semantic read end");
            received+=size_t(n);
        }
        unchanged();
    }

    // The transaction owner supplies its private staging directory. Creation
    // is exclusively a NEW leaf: this API cannot overwrite a live descriptor.
    // Failed staging is preserved. Session/pointer commit is a separate layer.
    void stageExclusiveCopy(int private_directory,const std::string &leaf)const
    {
        struct stat directory{};
        if(private_directory<0 || ::fstat(private_directory,&directory)!=0 || !S_ISDIR(directory.st_mode) ||
            leaf.empty() || leaf.size()>128 || leaf=="." || leaf=="..")
            throw std::invalid_argument("output prefix: private directory or leaf");
        for(unsigned char c:leaf)
            if(!((c>='A'&&c<='Z')||(c>='a'&&c<='z')||(c>='0'&&c<='9')||c=='_'||c=='-'||c=='.'))
                throw std::invalid_argument("output prefix: unsafe staging leaf");
        unchanged();
        Descriptor target(::openat(private_directory,leaf.c_str(),O_WRONLY|O_CREAT|O_EXCL|O_NOFOLLOW|O_CLOEXEC,0600));
        if(transfer(target.get())!=digest_)throw std::runtime_error("output prefix: changed staged digest");
        const auto after=regular(target.get());
        if(uint64_t(after.st_size)!=count_)throw std::runtime_error("output prefix: staged length");
        if(::fsync(target.get())!=0)throw std::runtime_error("output prefix: staging fsync");
    }
};
}
#endif
