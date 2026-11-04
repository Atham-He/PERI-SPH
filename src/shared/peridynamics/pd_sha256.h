/* Licensed under the Apache License, Version 2.0. */
#ifndef PD_SHA256_H
#define PD_SHA256_H
#include <algorithm>
#include <array>
#include <cerrno>
#include <cstdint>
#include <cstring>
#include <limits>
#include <stdexcept>
#include <string>
#if !defined(__linux__)
#error "This digest profile requires the frozen Linux AF_ALG provider"
#endif
#include <linux/if_alg.h>
#include <sys/socket.h>
#include <unistd.h>

namespace SPH::peridynamics::canonical_archive
{
using Digest=std::array<uint8_t,32>;
class Descriptor
{
    int value_;
  public:
    explicit Descriptor(int fd):value_(fd)
    {
        if(fd<0)
        {
            const int error=errno;
            throw std::runtime_error("archive: descriptor creation errno="+std::to_string(error)+": "+std::strerror(error));
        }
    }
    ~Descriptor(){::close(value_);}
    Descriptor(const Descriptor&)=delete;
    Descriptor &operator=(const Descriptor&)=delete;
    int get()const{return value_;}
};

/** The same installed SHA256 provider for archive bytes and large output
 * prefixes. Bounded buffers; no replacement provider or algorithm fallback.
 * A finalized stream cannot be reused, even after an I/O exception.
 */
class Sha256Stream
{
    static int algorithm()
    {
        const int fd=::socket(AF_ALG,SOCK_SEQPACKET|SOCK_CLOEXEC,0);
        if(fd<0)throw std::runtime_error("archive: SHA256 socket");
        sockaddr_alg address{};address.salg_family=AF_ALG;
        std::memcpy(address.salg_type,"hash",5);std::memcpy(address.salg_name,"sha256",7);
        if(::bind(fd,reinterpret_cast<sockaddr*>(&address),sizeof(address))!=0)
        {::close(fd);throw std::runtime_error("archive: SHA256 provider unavailable");}
        return fd;
    }
    Descriptor algorithm_,operation_;
    bool finalized_=false;
    uint64_t bytes_=0;
  public:
    Sha256Stream():algorithm_(algorithm()),
        operation_(::accept4(algorithm_.get(),nullptr,nullptr,SOCK_CLOEXEC)){}
    uint64_t bytes()const{return bytes_;}
    void update(const uint8_t *data,size_t count)
    {
        if(finalized_ || (count && !data) || count>std::numeric_limits<uint64_t>::max()-bytes_)
            throw std::invalid_argument("archive: SHA256 update state or size");
        size_t offset=0;
        while(offset<count)
        {
            const size_t chunk=std::min(size_t(16384),count-offset);
            const ssize_t sent=::send(operation_.get(),data+offset,chunk,MSG_MORE);
            if(sent<0 && errno==EINTR)continue;
            if(sent<=0){finalized_=true;throw std::runtime_error("archive: SHA256 update");}
            offset+=size_t(sent);bytes_+=uint64_t(sent);
        }
    }
    Digest finish()
    {
        if(finalized_)throw std::invalid_argument("archive: SHA256 already finalized");
        finalized_=true;
        ssize_t sent;
        do{sent=::send(operation_.get(),"",0,0);}while(sent<0 && errno==EINTR);
        if(sent!=0)throw std::runtime_error("archive: SHA256 finalization");
        Digest result{};ssize_t received;
        do{received=::recv(operation_.get(),result.data(),result.size(),0);}while(received<0 && errno==EINTR);
        if(received!=ssize_t(result.size()))throw std::runtime_error("archive: SHA256 output");
        return result;
    }
};
}
#endif
