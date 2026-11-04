/** Standalone byte engineering; no PD call, material history or live commit. */
#include "../src/shared/peridynamics/pd_output_prefix.h"
#include "../src/shared/peridynamics/pd_canonical_archive.h"
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <functional>
using namespace SPH::peridynamics::canonical_archive;
namespace
{
void require(bool ok,const char *message){if(!ok)throw std::runtime_error(message);}
std::string hex(const Digest &d)
{std::ostringstream out;for(auto x:d)out<<std::hex<<std::setw(2)<<std::setfill('0')<<unsigned(x);return out.str();}
Digest unhex(const std::string &s)
{
    require(s.size()==64,"expected64digestchars");Digest d{};
    auto digit=[](char c)->unsigned{if(c>='0'&&c<='9')return unsigned(c-'0');if(c>='a'&&c<='f')return unsigned(c-'a'+10);throw std::runtime_error("hex alphabet");};
    for(size_t k=0;k<32;++k)d[k]=uint8_t(16*digit(s[2*k])+digit(s[2*k+1]));
    return d;
}
std::string quote(const std::string &s)
{
    std::ostringstream out;out<<'"';
    for(unsigned char c:s)
        if(c=='"'||c=='\\')out<<'\\'<<c;
        else if(c<32)out<<"\\u00"<<std::hex<<std::setw(2)<<std::setfill('0')<<unsigned(c)<<std::dec;
        else out<<c;
    out<<'"';return out.str();
}
Digest chunks(const std::string &path)
{
    std::ifstream in(path,std::ios::binary);require(bool(in),"streamsource");
    Sha256Stream h;std::array<uint8_t,1031> buffer{};
    while(in.read(reinterpret_cast<char*>(buffer.data()),buffer.size()) || in.gcount())
        h.update(buffer.data(),size_t(in.gcount()));
    require(in.eof(),"streamread");return h.finish();
}
void mutate(const std::string &path,bool append)
{
    Descriptor f(::open(path.c_str(),O_RDWR|O_CLOEXEC|O_NOFOLLOW));
    if(append)require(::lseek(f.get(),0,SEEK_END)>=0,"seekappend");
    const uint8_t value=123;require(::write(f.get(),&value,1)==1,"mutatefixture");
    require(::fsync(f.get())==0,"fsyncfixture");
}
}
int main(int argc,char **argv)
{
    try
    {
        require(argc==4,"usage ROOT BYTES DIGEST");const std::string root=argv[1];
        size_t used=0;const uint64_t count=std::stoull(argv[2],&used);require(used==std::string(argv[2]).size(),"bytetail");
        const auto expected=unhex(argv[3]);Descriptor stage(::open((root+"/private").c_str(),O_RDONLY|O_DIRECTORY|O_CLOEXEC|O_NOFOLLOW));
        std::vector<std::pair<std::string,std::string>> rejected;
        const std::map<std::string,std::string> reasons{
            {"wrong_digest","selected digest mismatch"},
            {"cursor_past_eof","explicit byte cursor or snapshot size"},
            {"cursor_past_profile_bound","explicit byte cursor or snapshot size"},
            {"later_tail_not_an_exact_snapshot","explicit byte cursor or snapshot size"},
            {"changed_selected_byte","selected digest mismatch"},
            {"truncated_snapshot","explicit byte cursor or snapshot size"},
            {"symlink_source","descriptor creation errno=40"},
            {"directory_source","not a readable regular file"},
            {"existing_destination","descriptor creation errno=17"},
            {"symlink_destination","descriptor creation errno=17"},
            {"parent_traversal","unsafe staging leaf"},
            {"directory_dot","private directory or leaf"},
            {"bad_destination_descriptor","private directory or leaf"},
            {"source_changed_before_staging","source changed during validation or staging"},
            {"source_tail_changed_after_proof","source changed during validation or staging"},
            {"reuse_finalized_stream","SHA256 update state or size"},
            {"repeat_finalization","SHA256 already finalized"},
            {"null_nonempty_update","SHA256 update state or size"}};
        auto refuse=[&](const std::string &name,const std::function<void()> &operation)
        {
            try{operation();}catch(const std::exception &e)
            {
                if(std::string(e.what()).find(reasons.at(name))==std::string::npos)
                    throw std::runtime_error("wrong refusal reason for "+name+": "+e.what());
                rejected.emplace_back(name,e.what());return;
            }
            throw std::runtime_error("unexpected acceptance: "+name);
        };
        const VerifiedOutputPrefix prefix(root+"/live_with_tail.bin",count,expected,false);
        prefix.stageExclusiveCopy(stage.get(),"selected_prefix.bin");
        require(chunks(root+"/private/selected_prefix.bin")==expected,"staged prefix SHA");
        const VerifiedOutputPrefix snapshot(root+"/snapshot.bin",count,expected,true);
        snapshot.stageExclusiveCopy(stage.get(),"snapshot_copy.bin");
        auto wrong=expected;wrong[7]^=1;
        refuse("wrong_digest",[&]{VerifiedOutputPrefix p(root+"/snapshot.bin",count,wrong,true);});
        refuse("cursor_past_eof",[&]{VerifiedOutputPrefix p(root+"/snapshot.bin",count+1,expected,false);});
        refuse("cursor_past_profile_bound",[&]{VerifiedOutputPrefix p(root+"/snapshot.bin",maximum_output_prefix+1,expected,false);});
        refuse("later_tail_not_an_exact_snapshot",[&]{VerifiedOutputPrefix p(root+"/live_with_tail.bin",count,expected,true);});
        refuse("changed_selected_byte",[&]{VerifiedOutputPrefix p(root+"/corrupt.bin",count,expected,true);});
        refuse("truncated_snapshot",[&]{VerifiedOutputPrefix p(root+"/truncated.bin",count,expected,true);});
        refuse("symlink_source",[&]{VerifiedOutputPrefix p(root+"/source_link",count,expected,true);});
        refuse("directory_source",[&]{VerifiedOutputPrefix p(root+"/private",0,expected,false);});
        refuse("existing_destination",[&]{prefix.stageExclusiveCopy(stage.get(),"protected.bin");});
        refuse("symlink_destination",[&]{prefix.stageExclusiveCopy(stage.get(),"protected_link");});
        refuse("parent_traversal",[&]{prefix.stageExclusiveCopy(stage.get(),"../protected.bin");});
        refuse("directory_dot",[&]{prefix.stageExclusiveCopy(stage.get(),".");});
        refuse("bad_destination_descriptor",[&]{prefix.stageExclusiveCopy(-1,"new.bin");});
        {
            VerifiedOutputPrefix p(root+"/mutable.bin",count,expected,true);mutate(root+"/mutable.bin",false);
            refuse("source_changed_before_staging",[&]{p.stageExclusiveCopy(stage.get(),"must_not_exist_mutation.bin");});
        }
        {
            VerifiedOutputPrefix p(root+"/appendable.bin",count,expected,true);mutate(root+"/appendable.bin",true);
            refuse("source_tail_changed_after_proof",[&]{p.stageExclusiveCopy(stage.get(),"must_not_exist_append.bin");});
        }
        Sha256Stream empty;empty.update(nullptr,0);const auto zero=empty.finish();
        refuse("reuse_finalized_stream",[&]{empty.update(nullptr,0);});
        refuse("repeat_finalization",[&]{empty.finish();});
        Sha256Stream bad;refuse("null_nonempty_update",[&]{bad.update(nullptr,1);});
        const VerifiedOutputPrefix empty_prefix(root+"/empty.bin",0,zero,true);
        empty_prefix.stageExclusiveCopy(stage.get(),"empty_copy.bin");
        Builder archive(29);archive.reals("test/floats",{4},{-0.,.1,std::numeric_limits<double>::denorm_min(),std::numeric_limits<double>::max()});
        const auto bytes=archive.finish();Candidate decoded(bytes,29);require(decoded.reals("test/floats",{4}).size()==4,"canonical adapter");
        require(rejected.size()==18,"refusal count");
        std::ifstream boot("/proc/sys/kernel/random/boot_id");std::string boot_id;std::getline(boot,boot_id);
        std::ifstream stat("/proc/self/stat");std::string line;std::getline(stat,line);std::istringstream fields(line.substr(line.rfind(')')+2));
        std::string token,ticks;for(size_t k=0;k<=19;++k){fields>>token;if(k==19)ticks=token;}
        std::cout<<"{\"status\":\"byte_prefix_checks_passed_not_mechanical_CP\",\"physical_steps\":0,\"pid\":"<<::getpid()
            <<",\"ppid\":"<<::getppid()<<",\"boot_id\":"<<quote(boot_id)<<",\"start_ticks\":"<<quote(ticks)<<",\"argv\":[";
        for(int i=0;i<argc;++i){if(i)std::cout<<',';std::cout<<quote(argv[i]);}
        std::cout<<"],\"prefix_digest\":"<<quote(hex(prefix.digest()))<<",\"empty_digest\":"<<quote(hex(zero))
            <<",\"whole_file_digest\":"<<quote(hex(chunks(root+"/live_with_tail.bin")))<<",\"archive_hex\":\"";
        for(auto b:bytes)std::cout<<std::hex<<std::setw(2)<<std::setfill('0')<<unsigned(b);
        std::cout<<"\",\"refusals\":[";
        for(size_t i=0;i<rejected.size();++i){if(i)std::cout<<',';std::cout<<"["<<quote(rejected[i].first)<<','<<quote(rejected[i].second)<<']';}
        std::cout<<"]}\n";return 0;
    }
    catch(const std::exception &e){std::cerr<<quote(e.what())<<'\n';return 1;}
}
