// Optional reference-cell RKPM/S-NSNI elastic applicability operator.
// Wang et al., DOI 10.1007/s00466-024-02494-0, Eqs.18-19,43,48-55.
// Centroid quadrature and complete geometric moment tensors are explicit
// adaptations. This is NOT a qualified native SPH or finite-plastic integrator.
#pragma once
#include "solid_snsni_rk.h"
#include <tbb/parallel_for.h>
#include <fstream>
#include <numeric>
#include <map>
#include <cstdint>
namespace SPH::solid_dynamics {
struct SnsniReference
{
    using RKEntry=SnsniRKCloud::Entry;
    struct Entry {size_t j;Vecd gradient=Vecd::Zero();Matd hessian=Matd::Zero();};
    struct Q {Vecd x;Real w;};
    struct Face {Vecd normal;bool boundary;std::vector<Q> q;std::vector<RKEntry> integral;};
    struct Cell {Vecd center;Real volume;Matd moment;std::vector<int64_t> faces;std::vector<Entry> rows;std::vector<RKEntry> values;};
    using Field=StdVec<Vecd>;
    const Real dp,support=2.5;
    Real mu,lambda,patch_error=0,energy=0,minimum_mass=0;
    std::vector<Cell> cells;
    std::vector<Face> faces;
    SnsniRKCloud cloud;
    Field X,diagonal;
    StdVec<Real> mass,pressure;
    std::vector<std::vector<std::pair<size_t,size_t>>> transpose;
    std::vector<Matd> stress;
    std::vector<std::array<Matd,3>> stress_gradient;
    std::vector<Real> cell_energy;
    static Real kernel(Real r)
    {r=std::abs(r);return r<.5?2./3.-4*r*r+4*r*r*r:(r<1?4./3.*std::pow(1-r,3):0);}
    SnsniReference(const Field&positions,Real spacing,Real shear,Real lame,Real rho,const std::string&path)
        :dp(spacing),mu(shear),lambda(lame),cloud(positions,spacing),X(positions)
    {
        std::ifstream file(path,std::ios::binary);if(!file)throw std::runtime_error("Missing conforming reference cells");
        auto take=[&](auto &value){file.read(reinterpret_cast<char*>(&value),sizeof(value));if(!file)throw std::runtime_error("Truncated reference cells");};
        auto vector=[&](Vecd&v){for(int k=0;k<3;++k)take(v[k]);};
        uint64_t nc,nf;take(nc);take(nf);cells.resize(nc);faces.resize(nf);
        for(auto&c:cells){vector(c.center);take(c.volume);for(int a=0;a<3;++a)for(int b=0;b<3;++b)take(c.moment(a,b));uint32_t n;take(n);c.faces.resize(n);for(auto&id:c.faces)take(id);}
        for(auto&f:faces){vector(f.normal);uint32_t boundary,n;take(boundary);take(n);f.boundary=boundary;f.q.resize(n);for(auto&q:f.q){vector(q.x);take(q.w);}}
        tbb::parallel_for(size_t(0),faces.size(),[&](size_t f){std::map<size_t,RKEntry> sum;
            for(const auto&q:faces[f].q)for(const auto&e:shape(q.x)){auto it=sum.find(e.j);if(it==sum.end())it=sum.emplace(e.j,RKEntry{e.j,0,Vecd::Zero()}).first;it->second.w+=q.w*e.w;it->second.gradient+=q.w*e.gradient;}
            for(auto&[j,e]:sum)faces[f].integral.push_back(e);});
        tbb::parallel_for(size_t(0),cells.size(),[&](size_t l){auto&c=cells[l];std::map<size_t,Entry> sum;
            for(auto id:c.faces){const auto&f=faces[std::abs(id)-1];Vecd n=(id>0?1.:-1.)*f.normal;
                for(const auto&e:f.integral){auto[it,inserted]=sum.try_emplace(e.j);auto&a=it->second;a.j=e.j;a.gradient+=e.w*n/c.volume;
                    a.hessian+=.5*(e.gradient*n.transpose()+n*e.gradient.transpose())/c.volume;}}
            for(auto&[j,e]:sum)c.rows.push_back(e);c.values=shape(c.center);});
        faces.shrink_to_fit();transpose.resize(X.size());mass.assign(X.size(),0);diagonal.assign(X.size(),Vecd::Zero());pressure.resize(cells.size());stress.resize(cells.size());stress_gradient.resize(cells.size());cell_energy.resize(cells.size());
        for(size_t l=0;l<cells.size();++l){auto&c=cells[l];Vecd sum=Vecd::Zero();Matd affine=Matd::Zero();
            std::array<Matd,3> h_affine={Matd::Zero(),Matd::Zero(),Matd::Zero()};
            for(size_t k=0;k<c.rows.size();++k){const auto&e=c.rows[k];transpose[e.j].push_back({l,k});sum+=e.gradient;affine+=X[e.j]*e.gradient.transpose();
                for(int a=0;a<3;++a)h_affine[a]+=X[e.j][a]*e.hessian;
                Vecd diag=c.volume*(mu*e.gradient.squaredNorm()*Vecd::Ones()+(lambda+mu)*e.gradient.cwiseProduct(e.gradient));
                for(int a=0;a<3;++a)for(int b=0;b<3;++b)diag+=c.moment(a,b)*(mu*e.hessian.col(a).dot(e.hessian.col(b))*Vecd::Ones()+(lambda+mu)*e.hessian.col(a).cwiseProduct(e.hessian.col(b)));
                diagonal[e.j]+=diag;}
            patch_error=std::max({patch_error,dp*sum.norm(),(affine-Matd::Identity()).norm()});
            for(auto&h:h_affine)patch_error=std::max(patch_error,dp*h.norm());
            for(const auto&e:c.values)mass[e.j]+=rho*c.volume*e.w;}
        minimum_mass=*std::min_element(mass.begin(),mass.end());
        if(minimum_mass<=0||patch_error>1e-8)throw std::runtime_error("RK mass positivity or affine/rigid patch failed");
        std::cout<<"SNSNI cells="<<cells.size()<<" faces="<<faces.size()<<" patch="<<patch_error<<" min_mass="<<minimum_mass<<std::endl;
    }
    std::vector<RKEntry> shape(const Vecd&x)const{return cloud.shape(x);}
    Matd elastic(const Matd&gradient)const
    {Matd strain=.5*(gradient+gradient.transpose());return lambda*strain.trace()*Matd::Identity()+2*mu*strain;}
    Vecd value(size_t l,const Field&u)const
    {Vecd v=Vecd::Zero();for(const auto&e:cells[l].values)v+=e.w*u[e.j];return v;}
    void apply(const Field&u,Field&force)
    {
        tbb::parallel_for(size_t(0),cells.size(),[&](size_t l){const auto&c=cells[l];Matd grad=Matd::Zero();std::array<Matd,3> h={Matd::Zero(),Matd::Zero(),Matd::Zero()};
            for(const auto&e:c.rows){grad+=u[e.j]*e.gradient.transpose();for(int k=0;k<3;++k)h[k]+=u[e.j]*e.hessian.col(k).transpose();}
            stress[l]=elastic(grad);pressure[l]=-stress[l].trace()/3;cell_energy[l]=.5*c.volume*(stress[l].array()*grad.array()).sum();
            for(int k=0;k<3;++k)stress_gradient[l][k]=elastic(h[k]);
            for(int k=0;k<3;++k)for(int m=0;m<3;++m)cell_energy[l]+=.5*c.moment(k,m)*(stress_gradient[l][k].array()*h[m].array()).sum();});
        force.resize(X.size());tbb::parallel_for(size_t(0),X.size(),[&](size_t j){Vecd f=Vecd::Zero();for(auto[l,k]:transpose[j]){const auto&c=cells[l];const auto&e=c.rows[k];f+=c.volume*stress[l]*e.gradient;
            for(int a=0;a<3;++a)for(int b=0;b<3;++b)f+=c.moment(a,b)*stress_gradient[l][a]*e.hessian.col(b);}force[j]=f;});
        energy=std::accumulate(cell_energy.begin(),cell_energy.end(),Real(0));
    }
};

} // namespace SPH::solid_dynamics
