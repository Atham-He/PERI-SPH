#pragma once
// Persistence only. No constitutive, force, loading or timestep changes.
#include <chrono>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <type_traits>
#include <fcntl.h>
#include <unistd.h>
namespace recovery {
namespace fs=std::filesystem;
class Archive {
    std::fstream file;
    uint64_t hash=1469598103934665603ULL;
 public:
    const bool reading;
    Archive(const fs::path&p,bool read):reading(read) {
        file.open(p,std::ios::binary|(read?std::ios::in:(std::ios::out|std::ios::trunc)));
        if(!file)throw std::runtime_error("Cannot open checkpoint: "+p.string());
    }
    void bytes(void*p,size_t n) {
        if(reading)file.read(static_cast<char*>(p),n);else file.write(static_cast<char*>(p),n);
        if(!file)throw std::runtime_error("Incomplete checkpoint I/O");
        for(size_t i=0;i<n;++i){hash^=static_cast<unsigned char*>(p)[i];hash*=1099511628211ULL;}
    }
    template<class T,std::enable_if_t<std::is_arithmetic_v<T>,int> =0>void value(T&v){bytes(&v,sizeof(v));}
    template<class S,int R,int C,int O,int MR,int MC>void value(Eigen::Matrix<S,R,C,O,MR,MC>&v){bytes(v.data(),sizeof(S)*v.size());}
    template<class T,class A>void value(std::vector<T,A>&v) {
        uint64_t n=v.size();value(n);if(n>8000000)throw std::runtime_error("Invalid checkpoint vector length");
        if(reading)v.resize(n);for(auto&x:v)value(x);
    }
    template<class T,size_t N>void value(std::array<T,N>&v){for(auto&x:v)value(x);}
    void value(std::string&v){uint64_t n=v.size();value(n);if(n>1048576)throw std::runtime_error("Invalid checkpoint string");if(reading)v.resize(n);bytes(v.data(),n);}
    void value(SPH::dpc::State&v){value(v.inverse_Cp);value(v.kappa);value(v.alpha);}
    void value(SPH::dpc::Result&v){value(v.state);value(v.sigma);value(v.p);value(v.q);value(v.relative_density);value(v.plastic_shear);value(v.plastic_stress_work);value(v.elastic_energy);value(v.yield_residual);value(v.effective_d);value(v.surface);}
    void value(SPH::dpc::TrialField&v){value(v.sigma);value(v.gradient);}
    void value(SPH::dpc::SnsniMaterialStep&v){value(v.material);value(v.trial);value(v.correction);}
    template<class T>auto value(T&v)->decltype(v.kappa,v.damage,v.dissipation,void()){fields(v.kappa,v.damage,v.dissipation);}
    template<class...T>void fields(T&...v){(value(v),...);}
    template<class T>void span(T*p,size_t n){for(size_t i=0;i<n;++i)value(p[i]);}
    void finish() {
        uint64_t expected=hash;
        if(reading){file.read(reinterpret_cast<char*>(&expected),sizeof(expected));if(!file||expected!=hash||file.peek()!=std::char_traits<char>::eof())throw std::runtime_error("Checkpoint checksum/trailer mismatch");}
        else{file.write(reinterpret_cast<char*>(&expected),sizeof(expected));file.flush();if(!file)throw std::runtime_error("Checkpoint flush failed");}
        file.close();
    }
};
template<class B>void solid(Archive&a,B&b) {
    a.fields(b.history,b.pending,b.F,b.inverseF,b.correction,b.force_gradient,
        b.displacement_hessian,b.force_hessian,b.accumulated_plastic_work,b.internal,
        b.evaluations,b.commits,b.pending_commit,b.acoustic_dt);
    if(b.pending_commit)throw std::runtime_error("Only committed material states may be checkpointed");
}
template<class T>void integrator(Archive&a,T&t){a.fields(t.displacement,t.velocity,t.acceleration,t.damping_loss,t.steps);}
template<class C>void contact(Archive&a,C&c){
    if(c.friction!=0)throw std::runtime_error("Checkpoint scope requires the frozen frictionless contact");
    a.fields(c.upper,c.lower,c.energy,c.max_penetration,c.tangent_energy,c.slip_dissipation,c.release_energy,
        c.max_cone_ratio,c.tangential_force_sum,c.worst_force,c.worst_torque,c.worst_power);
}
class Control {
    std::chrono::steady_clock::time_point last=std::chrono::steady_clock::now();
    uint64_t serial=0;
    static void durable(const fs::path&p){int fd=::open(p.c_str(),O_RDONLY);if(fd<0)throw std::runtime_error("Checkpoint fsync open failed");int code=::fsync(fd);::close(fd);if(code)throw std::runtime_error("Checkpoint fsync failed");}
 public:
    std::string resume,fingerprint;
    uint64_t stop_step=0,schema_version=2;bool stop_after_equilibrium=false;
    bool take(const std::string&s){
        if(s=="--checkpoint-stop-after-equilibrium")stop_after_equilibrium=true;
        else if(s.rfind("--resume-state=",0)==0)resume=s.substr(15);
        else if(s.rfind("--checkpoint-fingerprint=",0)==0)fingerprint=s.substr(25);
        else if(s.rfind("--checkpoint-stop-step=",0)==0)stop_step=std::stoull(s.substr(23));
        else return false;return true;
    }
    bool enabled()const{return !fingerprint.empty();}
    bool due(size_t step,bool accepted=false)const{return enabled()&&(accepted||(stop_step&&step>=stop_step)||std::chrono::steady_clock::now()-last>=std::chrono::minutes(10));}
    template<class Transfer>void save(const fs::path&out,const std::vector<std::string>&csv,size_t step,Transfer&&transfer){
        if(!enabled())return;
        auto slot=serial++%2;auto target=out/("checkpoint_"+std::to_string(slot)+".bin"),temp=fs::path(target.string()+".tmp");
        Archive a(temp,false);uint64_t version=schema_version;std::string magic="PERI_SNSNI_RESTART",identity=fingerprint,folder=fs::absolute(out).string();
        a.fields(magic,version,identity,folder);auto names=csv;a.value(names);
        std::vector<uint64_t>sizes;for(auto&n:names)sizes.push_back(fs::file_size(out/n));a.value(sizes);
        transfer(a);a.finish();durable(temp);fs::rename(temp,target);
        auto manifest=out/"checkpoint.latest.json",mtmp=fs::path(manifest.string()+".tmp");
        {std::ofstream m(mtmp);m<<"{\"version\":"<<schema_version<<",\"file\":"<<std::quoted(target.string())<<",\"step\":"<<step<<",\"fingerprint\":"<<std::quoted(fingerprint)<<"}\n";m.flush();if(!m)throw std::runtime_error("Checkpoint manifest write failed");}
        durable(mtmp);fs::rename(mtmp,manifest);last=std::chrono::steady_clock::now();
        std::cout<<"CHECKPOINT step="<<step<<" file="<<target<<std::endl;
    }
    template<class Transfer>void load(const fs::path&out,const std::vector<std::string>&csv,Transfer&&transfer){
        if(resume.empty())return;
        if(!enabled())throw std::runtime_error("Resume requires a pinned input/source fingerprint");
        Archive a(resume,true);std::string magic,identity,folder;uint64_t version=0;
        a.fields(magic,version,identity,folder);
        if(magic!="PERI_SNSNI_RESTART"||version!=schema_version||identity!=fingerprint)throw std::runtime_error("Checkpoint identity/version mismatch");
        std::vector<std::string>names;std::vector<uint64_t>sizes;a.fields(names,sizes);
        if(names!=csv||sizes.size()!=names.size()||fs::absolute(folder)==fs::absolute(out))throw std::runtime_error("Resume must use a fresh output segment and matching CSV schema");
        transfer(a);a.finish();
        for(size_t i=0;i<names.size();++i){
            fs::path source=fs::path(folder)/names[i];if(fs::file_size(source)<sizes[i])throw std::runtime_error("Checkpoint CSV prefix missing");
            std::ifstream src(source,std::ios::binary);std::ofstream dst(out/names[i],std::ios::binary|std::ios::trunc);std::array<char,65536>buffer;uint64_t left=sizes[i];
            while(left){size_t n=std::min<uint64_t>(left,buffer.size());src.read(buffer.data(),n);dst.write(buffer.data(),n);if(!src||!dst)throw std::runtime_error("Checkpoint CSV prefix copy failed");left-=n;}
        }
        last=std::chrono::steady_clock::now();std::cout<<"RESUMED checkpoint="<<resume<<std::endl;
    }
};
inline Control checkpoint;
}
