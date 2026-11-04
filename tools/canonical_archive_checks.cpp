/** Byte-envelope engineering only, not a mechanical checkpoint qualification. */
#include "pd_canonical_archive.h"
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <sys/resource.h>
using namespace SPH::peridynamics::canonical_archive;
namespace
{
void require(bool value,const char *reason)
{if(!value)throw std::runtime_error(reason);}
template<class Action>void rejects(Action action,const std::string &reason)
{
    try{action();}catch(const std::invalid_argument &error)
    {require(std::string(error.what()).find(reason)!=std::string::npos,"unexpected refusal");return;}
    throw std::runtime_error("expected refusal missing");
}
std::string hex(const Digest &digest)
{
    const char *digits="0123456789abcdef";std::string out;
    for(uint8_t x:digest){out+=digits[x>>4];out+=digits[x&15];}return out;
}
void integerAt(Bytes &bytes,size_t offset,uint64_t value)
{for(size_t k=0;k<8;++k)bytes.at(offset+k)=uint8_t(value>>(8*k));}
void resign(Bytes &bytes)
{
    const auto digest=sha256(bytes,bytes.size()-40);
    std::copy(digest.begin(),digest.end(),bytes.end()-40);
}
size_t dataOffset(const Bytes &bytes,const std::string &wanted)
{
    Cursor c(bytes,bytes.size()-40);c.take(8);c.integer();c.integer();const auto fields=c.integer();
    for(uint64_t i=0;i<fields;++i)
    {
        const auto n=c.integer();const auto name=c.take(n);c.byte();const auto rank=c.byte();
        for(size_t k=0;k<rank;++k)c.integer();
        const auto size=c.integer();const auto offset=c.offset();
        if(std::string(name.begin(),name.end())==wanted)return offset;
        c.take(size);
    }
    throw std::runtime_error("test field missing");
}
}
int main(int argc,char **argv)
{
    try
    {
        require(argc==2,"usage OUTPUT");
        std::ifstream boot("/proc/sys/kernel/random/boot_id"),stat("/proc/self/stat");
        std::string id,line,word,ticks;std::getline(boot,id);std::getline(stat,line);
        std::istringstream fields(line.substr(line.rfind(')')+2));
        for(size_t k=0;k<=19;++k){fields>>word;if(k==19)ticks=word;}
        std::cout<<"{\"type\":\"identity\",\"pid\":"<<getpid()<<",\"ppid\":"<<getppid()
            <<",\"boot_id\":"<<std::quoted(id)<<",\"start_ticks\":"<<std::quoted(ticks)
            <<",\"argv\":["<<std::quoted(argv[0])<<","<<std::quoted(argv[1])<<"]}\n"<<std::flush;
        const Bytes abc{'a','b','c'},empty{};
        require(hex(sha256(abc,abc.size()))=="ba7816bf8f01cfea414140de5dae2223b00361a396177a9cb410ff61f20015ad","known SHA256");
        require(hex(sha256(empty,0))=="e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855","empty SHA256");
        Bytes large(1024*1024);for(size_t i=0;i<large.size();++i)large[i]=uint8_t(17*i+3);
        const auto large_digest=sha256(large,large.size());
        const std::vector<double> real{0.,-0.,.1,std::numeric_limits<double>::denorm_min(),
            std::numeric_limits<double>::max(),std::numeric_limits<double>::min()};
        Builder a(7);a.integers("accepted/step",{}, {19});a.reals("accepted/values",{2,3},real);
        a.bytes("reference/mask",{4},{1,0,1,1});const auto bytes=a.finish();
        Builder b(7);b.bytes("reference/mask",{4},{1,0,1,1});b.reals("accepted/values",{2,3},real);
        b.integers("accepted/step",{}, {19});require(bytes==b.finish(),"canonical registration order");
        Candidate c(bytes,7);c.requireExactNames({"accepted/step","accepted/values","reference/mask"});
        require(c.canonicalBytes()==bytes,"noadvance bytes");
        require(c.integers("accepted/step",{})==std::vector<uint64_t>{19},"integer decode");
        require(c.bytes("reference/mask",{4})==Bytes({1,0,1,1}),"byte decode");
        const auto got=c.reals("accepted/values",{2,3});
        for(size_t k=0;k<real.size();++k)require(doubleBits(got[k])==doubleBits(real[k]),"floating bits");
        rejects([&]{Candidate bad(bytes,8);},"generation");
        rejects([&]{c.requireExactNames({"accepted/step"});},"schema");
        rejects([&]{c.reals("accepted/values",{3,2});},"dimensions");
        rejects([&]{c.integers("accepted/values",{2,3});},"type");
        rejects([&]{Builder bad(0);},"generation");
        rejects([&]{a.integers("accepted/step",{}, {1});},"duplicate");
        rejects([&]{Builder bad(1);bad.reals("x",{2},{1});},"dimensions");
        rejects([&]{Builder bad(1);bad.reals("x",{}, {std::numeric_limits<double>::quiet_NaN()});},"nonfinite");
        rejects([&]{elements({maximum_elements,2});},"bound");
        rejects([&]{Builder bad(1);bad.bytes("bad name",{}, {0});},"alphabet");
        auto changed=bytes;changed[56]^=1;
        rejects([&]{Candidate bad(changed,7);},"SHA256");
        auto short_=bytes;short_.resize(short_.size()-8);
        rejects([&]{Candidate bad(short_,7);},"footer");
        changed=bytes;integerAt(changed,24,maximum_fields+1);resign(changed);
        rejects([&]{Candidate bad(changed,7);},"field count");
        changed=bytes;integerAt(changed,8,2);resign(changed);
        rejects([&]{Candidate bad(changed,7);},"schema");
        changed=bytes;const auto offset=dataOffset(bytes,"accepted/values");
        integerAt(changed,offset,UINT64_C(0x7ff8000000000000));resign(changed);
        rejects([&]{Candidate bad(changed,7);},"nonfinite");
        changed=bytes;integerAt(changed,55,16);resign(changed);
        rejects([&]{Candidate bad(changed,7);},"field size");
        // Detached getter output cannot mutate the candidate's original bytes.
        auto detached=got;detached[0]=1;require(c.canonicalBytes()==bytes,"detached decode mutated archive");
        std::ofstream output(argv[1],std::ios::binary|std::ios::out);
        output.write(reinterpret_cast<const char*>(bytes.data()),bytes.size());output.close();
        require(bool(output),"fixture write");
        rusage usage{};require(getrusage(RUSAGE_SELF,&usage)==0,"resource witness");
        std::cout<<"{\"status\":\"byte_envelope_engineering_passed\","
            "\"fixture_bytes\":"<<bytes.size()<<",\"large_SHA256\":\""<<hex(large_digest)<<"\","
            "\"peak_RSS_KiB\":"<<usage.ru_maxrss<<",\"physical_steps\":0,\"full_restart_qualified\":false}\n";return 0;
    }
    catch(const std::exception &error)
    {std::cerr<<error.what()<<'\n';return 1;}
}
