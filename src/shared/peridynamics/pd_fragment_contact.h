#pragma once
#include "pd_fracture_surface.h"
#include <unordered_map>
#include <array>
#include <fstream>

namespace SPH::peridynamics
{
/** Numerical contact proxies for cracked PD material and mass-retaining debris.
 * Existing live bonds are excluded from self-contact. A broken original pair
 * has an exclusion distance no larger than its original separation, so failure
 * cannot introduce initial overlap solely by changing the contact model.
 * Debris-core contact uses the closest point on a moving triangulated SPH
 * material surface, with barycentric AND MLS transpose force pullback.
 * Strength-controlled fracture and these proxies require separate physical
 * validation; no friction or fracture energy is silently supplied here.
 */
class PDFragmentContact
{
    using Cell=std::array<long long,3>;
    struct Hash {size_t operator()(const Cell&a)const{size_t h=0;for(auto x:a)h^=std::hash<long long>{}(x)+0x9e3779b9+(h<<6)+(h>>2);return h;}};
    SolidForcePort &pd;solid_dynamics::PDFractureReference &reference;
    ReferenceSolidSurface &core;const StdVec<Vecd> &vertices;
    StdVec<Matd> metric;StdVec<Real> radius;StdVec<std::array<size_t,3>> triangles;
    StdVec<StdVec<std::pair<size_t,size_t>>> lookup;
    Real K,cell_size;StdVec<Real> rate_pd,rate_core;
    Vecd sum=Vecd::Zero(),moment=Vecd::Zero();Real power=0,power_scale=0,force_scale=0,energy_rate=0,boundary_power=0;
  public:
    Real upper=0,lower=0,energy=0,stable_dt=std::numeric_limits<Real>::infinity();
    Real force_residual=0,torque_residual=0,power_residual=0;
    size_t self_contacts=0,core_contacts=0,plane_contacts=0;
    PDFragmentContact(SolidForcePort &p,solid_dynamics::PDFractureReference &r,ReferenceSolidSurface &c,
        const StdVec<Vecd>&v,const StdVec<Matd>&cell_metric,Real penalty,const std::string &mesh)
      :pd(p),reference(r),core(c),vertices(v),metric(cell_metric),K(penalty),rate_pd(p.count),rate_core(c.port.count)
    {
        if(metric.size()!=p.count||K<=0)throw std::invalid_argument("Fragment contact requires metrics and positive stiffness");
        radius.resize(p.count);Real largest=0;
        for(size_t i=0;i<p.count;++i){Eigen::SelfAdjointEigenSolver<Matd> e(metric[i]);
            if(e.info()!=Eigen::Success||e.eigenvalues().minCoeff()<=0)throw std::invalid_argument("Invalid fragment proxy metric");
            radius[i]=std::sqrt(e.eigenvalues().minCoeff());largest=std::max(largest,std::sqrt(e.eigenvalues().maxCoeff()));}
        cell_size=2*largest;lookup.resize(p.count);
        for(size_t i=0;i<p.count;++i){const auto &nb=reference.neighbors(i);for(size_t n=0;n<nb.current_size_;++n)lookup[i].push_back({nb.j_[n],n});std::sort(lookup[i].begin(),lookup[i].end());}
        std::ifstream f(mesh);size_t a,b,d;if(!f)throw std::invalid_argument("Missing core surface triangle file");
        while(f>>a>>b>>d){if(std::max({a,b,d})>=c.rows.size()||a==b||a==d||b==d)throw std::invalid_argument("Invalid core triangle");triangles.push_back({a,b,d});}
        if(triangles.empty())throw std::invalid_argument("Core contact triangulation is empty");
    }
    void apply(Real bottom,Real top,Real bottom_velocity=0,Real top_velocity=0)
    {
        upper=lower=energy=0;self_contacts=core_contacts=plane_contacts=0;
        sum.setZero();moment.setZero();power=power_scale=force_scale=energy_rate=boundary_power=0;
        std::fill(rate_pd.begin(),rate_pd.end(),0);std::fill(rate_core.begin(),rate_core.end(),0);
        std::unordered_map<Cell,StdVec<size_t>,Hash> cells;
        for(size_t i=0;i<pd.count;++i)cells[cell(pd.position[i])].push_back(i);
        for(size_t i=0;i<pd.count;++i)
        {
            Cell c=cell(pd.position[i]);
            for(int a=-1;a<=1;++a)for(int b=-1;b<=1;++b)for(int d=-1;d<=1;++d)
            {auto it=cells.find(Cell{c[0]+a,c[1]+b,c[2]+d});if(it==cells.end())continue;
             for(size_t j:it->second)if(j>i&&!hasLiveBond(i,j))pair(i,j);}
            if(reference.degree[i]!=0)continue;
            Real k=K*pd.volume[i]/(2*radius[i]);
            Real lo=std::max(bottom+radius[i]-pd.position[i][2],Real(0)),hi=std::max(pd.position[i][2]+radius[i]-top,Real(0));
            if(lo>0||hi>0)
            {Vecd f(0,0,k*(lo-hi));addPD(i,f);sum-=f;moment-=pd.position[i].cross(f); // analytical plate reaction
             Real boundary=f[2]*(lo>0?bottom_velocity:top_velocity);boundary_power+=boundary;
             energy_rate+=k*(lo*(bottom_velocity-pd.velocity[i][2])+hi*(pd.velocity[i][2]-top_velocity));
             energy+=.5*k*(lo*lo+hi*hi);upper+=k*hi;lower+=k*lo;rate_pd[i]+=k/pd.mass[i];++plane_contacts;}
        }
        if(reference.ballistic_particles)debrisCore();
        Real largest=0;for(Real a:rate_pd)largest=std::max(largest,a);for(Real a:rate_core)largest=std::max(largest,a);
        stable_dt=largest>0?.2/std::sqrt(largest):std::numeric_limits<Real>::infinity();
        Real L=cell_size;force_residual=sum.norm()/std::max(force_scale,Real(1e-30));
        torque_residual=moment.norm()/std::max(L*force_scale,Real(1e-30));
        power_residual=std::abs(power+energy_rate-boundary_power)/std::max({power_scale,std::abs(energy_rate),std::abs(boundary_power),Real(1e-30)});
    }
  private:
    Cell cell(const Vecd&x)const{return {(long long)std::floor(x[0]/cell_size),(long long)std::floor(x[1]/cell_size),(long long)std::floor(x[2]/cell_size)};}
    bool hasLiveBond(size_t i,size_t j)const
    {auto &v=lookup[i];auto p=std::lower_bound(v.begin(),v.end(),std::pair<size_t,size_t>{j,0});return p!=v.end()&&p->first==j&&reference.live(i,p->second);}
    void record(const Vecd&x,const Vecd&v,const Vecd&f){sum+=f;moment+=x.cross(f);power+=f.dot(v);power_scale+=std::abs(f.dot(v));force_scale+=f.norm();}
    void addPD(size_t i,const Vecd&f){pd.addForce(i,f);record(pd.position[i],pd.velocity[i],f);}
    void addCore(size_t q,const Vecd&f)
    {const auto&r=core.rows[q];for(size_t k=0;k<r.index.size();++k){size_t j=r.index[k];Vecd a=r.weight[k]*f;core.port.addForce(j,a);record(core.port.position[j],core.port.velocity[j],a);}}
    void pair(size_t i,size_t j)
    {
        Vecd xi=pd.reference[j]-pd.reference[i];Real initial=xi.norm();Vecd direction=xi/initial;
        Real D=std::min(initial,std::sqrt(direction.dot(metric[i]*direction))+std::sqrt(direction.dot(metric[j]*direction)));
        Vecd dx=pd.position[j]-pd.position[i];Real distance=dx.norm();if(distance>=D)return;
        if(distance<=1e-10*cell_size)throw std::runtime_error("Collapsed fragment contact pair");
        Real penetration=D-distance,V=2*pd.volume[i]*pd.volume[j]/(pd.volume[i]+pd.volume[j]),k=K*V/D;
        Vecd f=-k*penetration*dx/distance;addPD(i,f);addPD(j,-f);
        energy+=.5*k*penetration*penetration;energy_rate-=f.dot(pd.velocity[i]-pd.velocity[j]);
        Real a=2*k*(1+penetration/distance);rate_pd[i]+=a/pd.mass[i];rate_pd[j]+=a/pd.mass[j];++self_contacts;
    }
    static Vecd closest(const Vecd&p,const Vecd&a,const Vecd&b,const Vecd&c,Vecd&weights)
    {
        Vecd ab=b-a,ac=c-a,ap=p-a;Real d1=ab.dot(ap),d2=ac.dot(ap);
        if(d1<=0&&d2<=0){weights=Vecd(1,0,0);return a;}
        Vecd bp=p-b;Real d3=ab.dot(bp),d4=ac.dot(bp);if(d3>=0&&d4<=d3){weights=Vecd(0,1,0);return b;}
        Real vc=d1*d4-d3*d2;if(vc<=0&&d1>=0&&d3<=0){Real v=d1/(d1-d3);weights=Vecd(1-v,v,0);return a+v*ab;}
        Vecd cp=p-c;Real d5=ab.dot(cp),d6=ac.dot(cp);if(d6>=0&&d5<=d6){weights=Vecd(0,0,1);return c;}
        Real vb=d5*d2-d1*d6;if(vb<=0&&d2>=0&&d6<=0){Real w=d2/(d2-d6);weights=Vecd(1-w,0,w);return a+w*ac;}
        Real va=d3*d6-d5*d4;if(va<=0&&(d4-d3)>=0&&(d5-d6)>=0){Real w=(d4-d3)/((d4-d3)+(d5-d6));weights=Vecd(0,1-w,w);return b+w*(c-b);}
        Real den=va+vb+vc;if(std::abs(den)<=1e-40)throw std::runtime_error("Degenerate core contact triangle");
        Real v=vb/den,w=vc/den;weights=Vecd(1-v-w,v,w);return a+v*ab+w*ac;
    }
    void debrisCore()
    {
        std::unordered_map<Cell,StdVec<size_t>,Hash> grid;
        for(size_t q=0;q<triangles.size();++q)
        {
            auto t=triangles[q];Vecd lo=vertices[t[0]].cwiseMin(vertices[t[1]]).cwiseMin(vertices[t[2]]),hi=vertices[t[0]].cwiseMax(vertices[t[1]]).cwiseMax(vertices[t[2]]);
            Cell a=cell(lo),b=cell(hi);for(long long x=a[0];x<=b[0];++x)for(long long y=a[1];y<=b[1];++y)for(long long z=a[2];z<=b[2];++z)grid[{x,y,z}].push_back(q);
        }
        StdVec<size_t> marked(triangles.size(),pd.count);
        for(size_t i=0;i<pd.count;++i)if(reference.degree[i]==0)
        {
            Cell c=cell(pd.position[i]);Real best=std::numeric_limits<Real>::infinity();size_t best_triangle=triangles.size();Vecd bary=Vecd::Zero(),point=Vecd::Zero();
            for(int a=-2;a<=2;++a)for(int b=-2;b<=2;++b)for(int d=-2;d<=2;++d)
            {auto it=grid.find({c[0]+a,c[1]+b,c[2]+d});if(it==grid.end())continue;
             for(size_t q:it->second)if(marked[q]!=i){marked[q]=i;auto t=triangles[q];Vecd w;Vecd p=closest(pd.position[i],vertices[t[0]],vertices[t[1]],vertices[t[2]],w);
                Real dist=(p-pd.position[i]).squaredNorm();if(dist<best){best=dist;best_triangle=q;bary=w;point=p;}}}
            if(best_triangle==triangles.size())continue;
            auto t=triangles[best_triangle];Vecd normal=(vertices[t[1]]-vertices[t[0]]).cross(vertices[t[2]]-vertices[t[0]]);
            Vecd dx=pd.position[i]-point;Real distance=std::sqrt(best);if(normal.squaredNorm()<1e-40)throw std::runtime_error("Collapsed core surface face");
            Real sign=normal.dot(dx)>=0?1:-1,penetration=radius[i]-sign*distance;if(penetration<=0)continue;
            if(distance<=1e-10*cell_size)throw std::runtime_error("Debris lies exactly on a core contact face");
            Real k=K*pd.volume[i]/(2*radius[i]);Vecd f=k*penetration*sign*dx/distance;addPD(i,f);Vecd vq=Vecd::Zero();Real norm=1;
            for(size_t a=0;a<3;++a){addCore(t[a],-bary[a]*f);vq+=bary[a]*core.velocity(t[a]);for(Real w:core.rows[t[a]].weight)norm+=bary[a]*std::abs(w);}
            energy+=.5*k*penetration*penetration;energy_rate-=f.dot(pd.velocity[i]-vq);
            Real factor=k*(1+penetration/distance);rate_pd[i]+=factor*norm/pd.mass[i];
            for(size_t a=0;a<3;++a){const auto&r=core.rows[t[a]];for(size_t b=0;b<r.index.size();++b)rate_core[r.index[b]]+=factor*norm*bary[a]*std::abs(r.weight[b])/core.port.mass[r.index[b]];}
            ++core_contacts;
        }
    }
};
}
