#ifndef PD_CHECKPOINT_H
#define PD_CHECKPOINT_H

#include "base_particles.hpp"
#include "base_body.h"
#include "particle_generator_lattice.h"
#include "neighborhood.h"
#include <array>
#include <chrono>
#include <csignal>
#include <cstdint>
#include <limits>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <functional>
#include <iomanip>
#include <sstream>
#include <tuple>
#include <thread>
#include <type_traits>
#ifndef _WIN32
#include <fcntl.h>
#include <unistd.h>
#endif

#if defined(__GNUC__) || defined(__clang__)
#define PERI_CHECKPOINT_NOINLINE __attribute__((noinline))
#else
#define PERI_CHECKPOINT_NOINLINE
#endif

namespace SPH::peridynamics
{
/** Opt-in particle container: enumerate every registered state, including IDs
 * and constitutive scratch/history. No upstream particle API is changed. */
class CheckpointParticles : public BaseParticles
{
  public:
    using BaseParticles::BaseParticles;
    DiscreteVariables &checkpointVariables() { return all_discrete_variables_; }
    SingleVariables &checkpointSingles() { return all_singular_variables_; }
};
} // namespace SPH::peridynamics
namespace SPH
{
template<> class ParticleGenerator<peridynamics::CheckpointParticles,Lattice> : public ParticleGenerator<BaseParticles,Lattice>
{
  public:
    using ParticleGenerator<BaseParticles,Lattice>::ParticleGenerator;
};
}
namespace SPH::peridynamics
{
namespace checkpoint_detail
{
namespace fs = std::filesystem;
constexpr uint64_t offset = 14695981039346656037ULL, prime = 1099511628211ULL;
constexpr uint64_t footer = 0x50455249444f4e45ULL;
inline void hashBytes(uint64_t &hash, const char *data, size_t n)
{for (size_t i=0;i<n;++i) {hash ^= static_cast<unsigned char>(data[i]); hash *= prime;}}
inline uint64_t fileHash(const fs::path &path, uint64_t limit=UINT64_MAX)
{
    std::ifstream in(path,std::ios::binary); if(!in) throw std::runtime_error("checkpoint: cannot read "+path.string());
    std::array<char,65536> buffer; uint64_t hash=offset;
    while(limit && in) {size_t n=std::min<uint64_t>(buffer.size(),limit);in.read(buffer.data(),n);size_t got=in.gcount();hashBytes(hash,buffer.data(),got);limit-=got;}
    return hash;
}
inline void durable(const fs::path &path, bool directory=false)
{
#ifndef _WIN32
    int fd=::open(path.c_str(),directory ? O_RDONLY|O_DIRECTORY : O_RDONLY);
    if(fd<0) throw std::runtime_error("checkpoint: open for fsync failed: "+path.string());
    int result=::fsync(fd);::close(fd);
    if(result) throw std::runtime_error("checkpoint: fsync failed: "+path.string());
#else
    throw std::runtime_error("durable PD checkpoints currently require a POSIX filesystem");
#endif
}
inline void verifyFile(const fs::path &path)
{
    uint64_t size=fs::file_size(path),hash=0,magic=0;
    if(size<32) throw std::runtime_error("checkpoint: truncated state");
    std::ifstream in(path,std::ios::binary);in.seekg(size-16);in.read(reinterpret_cast<char*>(&hash),8);in.read(reinterpret_cast<char*>(&magic),8);
    if(!in || magic!=footer || fileHash(path,size-16)!=hash) throw std::runtime_error("checkpoint: checksum/commit footer mismatch");
}
class Archive
{
  public:
    enum Mode {Write,Validate,Restore};
    Archive(const fs::path &path,Mode mode): mode_(mode)
    {
        if(mode==Write) out_.open(path,std::ios::binary|std::ios::trunc);
        else in_.open(path,std::ios::binary);
        if(!(mode==Write ? bool(out_) : bool(in_))) throw std::runtime_error("checkpoint: cannot open archive");
    }
    Mode mode() const {return mode_;}
    void physical(bool enabled) {physical_enabled_=enabled;}
    uint64_t physicalHash() const {return physical_hash_;}
    void bytes(char *data,size_t n)
    {
        if(mode_==Write) {out_.write(data,n);if(!out_) throw std::runtime_error("checkpoint: write failed");}
        else {in_.read(data,n);if(!in_) throw std::runtime_error("checkpoint: truncated archive");}
        hashBytes(hash_,data,n);if(physical_enabled_) hashBytes(physical_hash_,data,n);
    }
    template<class T> void scalar(T &value)
    {
        static_assert(std::is_arithmetic_v<T>);
        T temporary=value;bytes(reinterpret_cast<char*>(&temporary),sizeof(T));
        if constexpr(std::is_floating_point_v<T>) if(!std::isfinite(temporary)) throw std::runtime_error("checkpoint: non-finite scalar");
        if(mode_==Restore) value=temporary;
    }
    uint64_t number(uint64_t value) {bytes(reinterpret_cast<char*>(&value),sizeof(value));return value;}
    void expect(uint64_t value)
    {
        auto stored=number(value);
        if(stored!=value) throw std::runtime_error("checkpoint: incompatible schema/count after "+context_+" (expected "+std::to_string(value)+", stored "+std::to_string(stored)+", offset "+std::to_string(static_cast<uint64_t>(in_.tellg()))+")");
    }
    void label(const std::string &name)
    {
        expect(name.size());std::string read=name; if(!read.empty()) bytes(read.data(),read.size());
        if(read!=name) throw std::runtime_error("checkpoint: incompatible variable/configuration: "+name);
        context_=name;
    }
    template<class T> PERI_CHECKPOINT_NOINLINE void array(T *data,size_t count,size_t max_index=SIZE_MAX)
    {
        if constexpr(std::is_arithmetic_v<T>)
        {
            std::array<T,4096> buffer;
            for(size_t at=0;at<count;at+=buffer.size())
            {
                size_t n=std::min(buffer.size(),count-at);
                if(mode_==Write) std::copy_n(data+at,n,buffer.data());
                bytes(reinterpret_cast<char*>(buffer.data()),n*sizeof(T));
                for(size_t k=0;k<n;++k)
                {
                    if constexpr(std::is_floating_point_v<T>) {if(!std::isfinite(buffer[k])) throw std::runtime_error("checkpoint: non-finite particle state");}
                    if constexpr(std::is_unsigned_v<T>) {if(max_index!=SIZE_MAX && buffer[k]>=max_index) throw std::runtime_error("checkpoint: invalid neighbor index");}
                }
                if(mode_==Restore) std::copy_n(buffer.data(),n,data+at);
            }
        }
        else
        {
            constexpr size_t width=T::SizeAtCompileTime;
            static_assert(width>0 && width<1024,"checkpoint requires fixed-size Eigen values");
            std::array<Real,width*128> buffer;
            for(size_t at=0;at<count;at+=128)
            {
                size_t n=std::min(size_t(128),count-at);
                if(mode_==Write) for(size_t k=0;k<n;++k) std::copy_n(data[at+k].data(),width,buffer.data()+k*width);
                bytes(reinterpret_cast<char*>(buffer.data()),n*width*sizeof(Real));
                for(size_t k=0;k<n*width;++k) if(!std::isfinite(buffer[k])) throw std::runtime_error("checkpoint: non-finite matrix/vector");
                if(mode_==Restore) for(size_t k=0;k<n;++k) std::copy_n(buffer.data()+k*width,width,data[at+k].data());
            }
        }
    }
    template<class T> PERI_CHECKPOINT_NOINLINE void discrete(StdVec<DiscreteVariable<T>*> &variables)
    {
        expect(variables.size());
        for(auto *v:variables) {label(v->Name());expect(v->getSize());expect(v->getWidth());array(v->Data(),v->getTotalSize());}
    }
    template<class T> PERI_CHECKPOINT_NOINLINE void singles(StdVec<SingleVariable<T>*> &variables)
    {expect(variables.size());for(auto *v:variables) {label(v->Name());array(v->Data(),1);}}
    PERI_CHECKPOINT_NOINLINE void configuration(ParticleConfiguration &values,size_t target)
    {
        expect(values.size());
        for(auto &n:values)
        {
            uint64_t current=number(n.current_size_),allocated=number(n.allocated_size_);
            if(current>allocated || allocated>target) throw std::runtime_error("checkpoint: invalid neighborhood size");
            if(mode_==Restore)
            {
                n.j_.resize(allocated);n.W_ij_.resize(allocated);n.dW_ij_.resize(allocated);
                n.r_ij_.resize(allocated);n.e_ij_.resize(allocated);n.bondLive_.resize(allocated);
                n.current_size_=current;n.allocated_size_=allocated;
            }
            if(mode_==Write && (n.j_.size()!=allocated || n.bondLive_.size()!=allocated)) throw std::runtime_error("checkpoint: inconsistent neighborhood storage");
            array(n.j_.data(),allocated,target);array(n.W_ij_.data(),allocated);array(n.dW_ij_.data(),allocated);
            array(n.r_ij_.data(),allocated);array(n.e_ij_.data(),allocated);
            for(size_t j=0;j<allocated;++j)
            {
                uint8_t live=mode_==Write ? uint8_t(n.bondLive_[j]) : 0;bytes(reinterpret_cast<char*>(&live),1);
                if(live>1) throw std::runtime_error("checkpoint: invalid bond flag");
                if(mode_==Restore) n.bondLive_[j]=live;
            }
        }
    }
    void attachment(const fs::path &output,const std::string &name)
    {
        label(name);uint64_t size=mode_==Write ? fs::file_size(output/name) : 0;size=number(size);
        if(size>uint64_t(1)<<32) throw std::runtime_error("checkpoint: output prefix too large");
        std::ifstream source;std::ofstream destination;
        if(mode_==Write) source.open(output/name,std::ios::binary);
        if(mode_==Restore) destination.open(output/(name+".resume.tmp"),std::ios::binary|std::ios::trunc);
        std::array<char,65536> buffer;
        while(size)
        {
            size_t n=std::min<uint64_t>(buffer.size(),size);
            if(mode_==Write) {source.read(buffer.data(),n);if(!source) throw std::runtime_error("checkpoint: output prefix changed");}
            bytes(buffer.data(),n);
            if(mode_==Restore) {destination.write(buffer.data(),n);if(!destination) throw std::runtime_error("checkpoint: output restore failed");}
            size-=n;
        }
    }
    uint64_t finish()
    {
        if(mode_==Write) {out_.write(reinterpret_cast<char*>(&hash_),8);uint64_t magic=footer;out_.write(reinterpret_cast<char*>(&magic),8);out_.flush();if(!out_) throw std::runtime_error("checkpoint: flush failed");out_.close();}
        return hash_;
    }
  private:
    std::string context_="archive header";
    Mode mode_;std::ifstream in_;std::ofstream out_;uint64_t hash_=offset,physical_hash_=offset;bool physical_enabled_=false;
};
} // namespace checkpoint_detail

struct CheckpointProgress {Real time=0,wall_seconds=0;uint64_t acoustic_steps=0,pd_steps=0,completed_frame=0;};
/** Complete CPU state and neighborhood ordering. Publication is generation-atomic;
 * validation and checksums finish before any live state or output is restored. */
class DiskCheckpoint
{
  public:
    DiskCheckpoint(std::filesystem::path root,std::string compatibility,std::filesystem::path output="output")
        :root_(std::filesystem::absolute(root)),compatibility_(std::move(compatibility)),output_(std::move(output)) {}
    void addBody(CheckpointParticles &p) {bodies_.push_back(&p);}
    void addConfiguration(std::string name,ParticleConfiguration &configuration,size_t target)
    {configurations_.push_back({std::move(name),&configuration,target});}
    template<class T> void addScalar(std::string name,T &value)
    {globals_.push_back([name,value_ptr=&value](checkpoint_detail::Archive &a){a.label(name);a.scalar(*value_ptr);});}
    template<class T> void addArray(std::string name,StdVec<T> &value)
    {globals_.push_back([name,value_ptr=&value](checkpoint_detail::Archive &a){a.label(name);a.expect(value_ptr->size());a.array(value_ptr->data(),value_ptr->size());});}
    void addOutput(std::string filename)
    {if(std::filesystem::path(filename).filename()!=filename) throw std::runtime_error("checkpoint: unsafe output filename");outputs_.push_back(std::move(filename));}
    const std::filesystem::path &loadedGeneration() const {return loaded_;}
    std::filesystem::path save(CheckpointProgress progress,bool fail_before_commit=false)
    {
        namespace fs=std::filesystem;using namespace checkpoint_detail;
        fs::create_directories(root_);
        std::ostringstream id;id<<"checkpoint-"<<std::setw(20)<<std::setfill('0')<<progress.acoustic_steps<<"-"<<std::chrono::steady_clock::now().time_since_epoch().count();
        fs::path temporary=root_/(".tmp-"+id.str()),committed=root_/id.str();fs::create_directory(temporary);
        Archive a(temporary/"state.bin",Archive::Write);transfer(a,progress);uint64_t hash=a.finish();durable(temporary/"state.bin");
        std::ofstream manifest(temporary/"manifest.json");
        manifest<<std::setprecision(17)<<"{\n  \"format\": 1, \"status\": \"complete\",\n  \"time_s\": "<<progress.time<<", \"acoustic_steps\": "<<progress.acoustic_steps<<", \"pd_steps\": "<<progress.pd_steps<<",\n  \"completed_frame\": "<<progress.completed_frame<<", \"wall_seconds\": "<<progress.wall_seconds<<",\n  \"state_bytes\": "<<fs::file_size(temporary/"state.bin")<<", \"fnv1a64\": \""<<std::hex<<hash<<"\", \"physical_state_fnv1a64\": \""<<physical_hash_<<"\"\n}\n";
        manifest.close();if(!manifest) throw std::runtime_error("checkpoint: manifest failed");durable(temporary/"manifest.json");durable(temporary,true);
        if(fail_before_commit) throw std::runtime_error("checkpoint: injected interruption before commit");
        fs::rename(temporary,committed);durable(root_,true);
        publishLatest(root_,id.str());
        auto all=generations(root_);for(size_t i=2;i<all.size();++i)
            if(all[i]!=committed && all[i].parent_path()==root_ && !fs::is_symlink(all[i])) fs::remove_all(all[i]);
        return committed;
    }
    CheckpointProgress restore(const std::filesystem::path &requested)
    {
        namespace fs=std::filesystem;using namespace checkpoint_detail;
        fs::path selected=fs::absolute(requested);StdVec<fs::path> candidates;
        if(fs::is_regular_file(selected/"state.bin")) candidates.push_back(selected);
        else
        {
            candidates=generations(selected);std::ifstream latest(selected/"LATEST");std::string name;std::getline(latest,name);
            if(!name.empty() && fs::path(name).filename()==name)
            {auto first=selected/name;auto found=std::find(candidates.begin(),candidates.end(),first);if(found!=candidates.end()) std::rotate(candidates.begin(),found,found+1);}
        }
        std::string failures;StdVec<fs::path> rejected;
        for(const auto &path:candidates)
        {
            try {verifyFile(path/"state.bin");CheckpointProgress check;Archive validation(path/"state.bin",Archive::Validate);transfer(validation,check);}
            catch(const std::exception &error) {failures+=path.filename().string()+": "+error.what()+"; ";rejected.push_back(path);continue;}
            CheckpointProgress progress;Archive apply(path/"state.bin",Archive::Restore);transfer(apply,progress);
            for(const auto &name:outputs_) {durable(output_/(name+".resume.tmp"));fs::rename(output_/(name+".resume.tmp"),output_/name);}
            durable(output_,true);loaded_=path;
            // Quarantine rejected generations only after successful recovery.
            // Otherwise retention could evict the last valid state on the next save.
            if(!fs::is_regular_file(selected/"state.bin"))
            {
                for(const auto &bad:rejected) fs::rename(bad,selected/("rejected-"+bad.filename().string()));
                publishLatest(selected,path.filename().string());
            }
            return progress;
        }
        throw std::runtime_error("checkpoint: no compatible committed generation: "+failures);
    }
  private:
    struct Configuration {std::string name;ParticleConfiguration *data;size_t target;};
    std::filesystem::path root_,output_,loaded_;std::string compatibility_;
    uint64_t physical_hash_=0;
    StdVec<CheckpointParticles*> bodies_;StdVec<Configuration> configurations_;
    StdVec<std::function<void(checkpoint_detail::Archive&)>> globals_;StdVec<std::string> outputs_;
    static void publishLatest(const std::filesystem::path &root,const std::string &name)
    {
        namespace fs=std::filesystem;using namespace checkpoint_detail;
        {std::ofstream latest(root/"LATEST.tmp");latest<<name<<'\n';latest.close();if(!latest) throw std::runtime_error("checkpoint: LATEST failed");}
        durable(root/"LATEST.tmp");fs::rename(root/"LATEST.tmp",root/"LATEST");durable(root,true);
    }
    static StdVec<std::filesystem::path> generations(const std::filesystem::path &root)
    {
        namespace fs=std::filesystem;StdVec<fs::path> paths;
        if(!fs::is_directory(root)) return paths;
        for(auto &entry:fs::directory_iterator(root)) if(entry.is_directory() && !entry.is_symlink() && entry.path().filename().string().rfind("checkpoint-",0)==0 && fs::is_regular_file(entry.path()/"manifest.json") && fs::is_regular_file(entry.path()/"state.bin")) paths.push_back(entry.path());
        std::sort(paths.begin(),paths.end(),std::greater<fs::path>());return paths;
    }
    PERI_CHECKPOINT_NOINLINE void transfer(checkpoint_detail::Archive &a,CheckpointProgress &progress)
    {
        a.label("PERI-SPH-CPU-CHECKPOINT-1");a.expect(0x0102030405060708ULL);a.expect(sizeof(Real));a.expect(sizeof(size_t));a.expect(Dimensions);a.label(compatibility_);
        a.scalar(progress.time);a.scalar(progress.wall_seconds);a.scalar(progress.acoustic_steps);a.scalar(progress.pd_steps);a.scalar(progress.completed_frame);
        a.physical(true);a.expect(bodies_.size());
        for(auto *p:bodies_)
        {
            a.label(p->getBodyName());a.expect(p->ParticlesBound());a.expect(p->TotalRealParticles());
            std::apply([&](auto &...list){(a.discrete(list),...);},p->checkpointVariables());
            std::apply([&](auto &...list){(a.singles(list),...);},p->checkpointSingles());
            uint64_t dirty=a.number(p->getSPHBody().checkNewlyUpdated());
            if(dirty>1) throw std::runtime_error("checkpoint: invalid body update flag");
            if(a.mode()==checkpoint_detail::Archive::Restore) {if(dirty) p->getSPHBody().setNewlyUpdated();else p->getSPHBody().setNotNewlyUpdated();}
        }
        a.expect(configurations_.size());for(auto &configuration:configurations_) {a.label(configuration.name);a.configuration(*configuration.data,configuration.target);}
        a.expect(globals_.size());for(auto &global:globals_) global(a);
        physical_hash_=a.physicalHash();a.physical(false);
        a.expect(outputs_.size());for(const auto &name:outputs_) a.attachment(output_,name);
    }
};

/** Preserve outputs while upstream IO::initEnvironment initializes its default
 * output/backup folders. A durable rename journal also recovers an interrupted
 * initialization. No existing result is deleted by this opt-in guard. */
class CheckpointOutputGuard
{
  public:
    explicit CheckpointOutputGuard(bool enabled): enabled_(enabled)
    {
        if(!enabled_) return;
        recover();
        std::string suffix=std::to_string(std::chrono::steady_clock::now().time_since_epoch().count());
        for(const char *name:{"output","output_backup"})
            if(std::filesystem::exists(name)) moves_.push_back({name,std::string(name)+".pd-preserved-"+suffix});
        {std::ofstream journal(marker_);journal<<"PERI-OUTPUT-GUARD-1\n";for(auto &move:moves_) journal<<move.first<<' '<<move.second<<'\n';journal.close();if(!journal) throw std::runtime_error("checkpoint: output journal failed");}
        checkpoint_detail::durable(marker_);checkpoint_detail::durable(".",true);
        for(auto &move:moves_) std::filesystem::rename(move.first,move.second);
        checkpoint_detail::durable(".",true);
    }
    void restore()
    {
        if(!enabled_) return;
        recover();enabled_=false;
    }
    ~CheckpointOutputGuard(){if(enabled_) {try{recover();}catch(...) {}}}
  private:
    bool enabled_;const std::string marker_=".peri-output-restore-journal";
    StdVec<std::pair<std::string,std::string>> moves_;
    void recover()
    {
        namespace fs=std::filesystem;
        if(!fs::exists(marker_)) return;
        std::ifstream journal(marker_);std::string magic;std::getline(journal,magic);
        if(magic!="PERI-OUTPUT-GUARD-1") throw std::runtime_error("checkpoint: invalid output recovery journal");
        std::string original,saved;
        while(journal>>original>>saved)
        {
            if((original!="output" && original!="output_backup") || fs::path(saved).filename()!=saved || saved.rfind(original+".pd-preserved-",0)!=0)
                throw std::runtime_error("checkpoint: unsafe output recovery path");
            if(!fs::exists(saved)) continue;
            if(fs::is_symlink(saved) || !fs::is_directory(saved)) throw std::runtime_error("checkpoint: invalid preserved output");
            if(fs::exists(original))
            {
                auto artifacts=original+".init-artifacts-"+std::to_string(std::chrono::steady_clock::now().time_since_epoch().count());
                fs::rename(original,artifacts);
            }
            fs::rename(saved,original);
        }
        journal.close();fs::remove(marker_);checkpoint_detail::durable(".",true);
    }
};

inline volatile std::sig_atomic_t checkpoint_requested=0,checkpoint_stop_requested=0;
inline void checkpointSignal(int value)
{
    checkpoint_requested=1;
#ifndef _WIN32
    if(value!=SIGUSR1) checkpoint_stop_requested=1;
#else
    checkpoint_stop_requested=1;
#endif
}
struct CheckpointOptions
{
    std::filesystem::path directory="checkpoints",resume;
    Real seconds=600;uint64_t steps=0;size_t threads=std::thread::hardware_concurrency();bool stop_after_checkpoint=false,final_checkpoint=true,describe_only=false;
    void parse(int &argc,char **argv)
    {
        int out=1;
        for(int i=1;i<argc;++i)
        {
            std::string option(argv[i]);size_t equal=option.find('=');std::string key=option.substr(0,equal);
            auto value=[&](){if(equal!=std::string::npos) return option.substr(equal+1);if(++i>=argc) throw std::invalid_argument(key+" requires a value");return std::string(argv[i]);};
            if(key=="--pd-checkpoint-capabilities") describe_only=true;
            else if(key=="--pd-resume") resume=value();
            else if(key=="--pd-checkpoint-dir") directory=value();
            else if(key=="--pd-checkpoint-seconds")
            {auto token=value();size_t used=0;seconds=std::stod(token,&used);if(used!=token.size()) throw std::invalid_argument("invalid checkpoint seconds");}
            else if(key=="--pd-checkpoint-steps" || key=="--pd-threads")
            {
                auto token=value();if(token.empty() || token.find_first_not_of("0123456789")!=std::string::npos) throw std::invalid_argument("invalid checkpoint integer");
                auto count=std::stoull(token);
                if(key=="--pd-checkpoint-steps") steps=count;
                else {if(count>std::numeric_limits<size_t>::max()) throw std::invalid_argument("thread count overflow");threads=count;}
            }
            else if(key=="--pd-stop-after-checkpoint") stop_after_checkpoint=true;
            else if(key=="--pd-no-final-checkpoint") final_checkpoint=false;
            else argv[out++]=argv[i];
        }
        argc=out;argv[out]=nullptr;
        if(!std::isfinite(seconds) || seconds<0 || threads<1) throw std::invalid_argument("invalid checkpoint interval");
#ifndef _WIN32
        std::signal(SIGUSR1,checkpointSignal);std::signal(SIGUSR2,checkpointSignal);std::signal(SIGTERM,checkpointSignal);std::signal(SIGINT,checkpointSignal);
#endif
    }
};
} // namespace SPH::peridynamics
#undef PERI_CHECKPOINT_NOINLINE
#endif
