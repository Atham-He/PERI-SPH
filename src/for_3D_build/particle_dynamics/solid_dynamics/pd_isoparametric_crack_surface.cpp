#include "pd_isoparametric_crack_surface.h"

#include <iomanip>
#include <sstream>
#include <unordered_map>

namespace SPH::peridynamics
{
BilinearCrackGeometry::Param BilinearCrackGeometry::closest(const std::array<Vecd,4>&global,const Vecd&point,Param p){
    // Translate once before the Newton solve. Thin rotated faces must not
    // lose stationarity through subtraction of large global coordinates.
    std::array<Vecd,4>x;for(size_t i=0;i<4;++i)x[i]=global[i]-global[0];Vecd q=point-global[0];
    Real coordinate_scale=std::max({q.norm(),x[1].norm(),x[2].norm(),x[3].norm(),Real(1e-30)});
    for(size_t iteration=0;iteration<40;++iteration){auto a=evaluate(x,p);Vecd r=q-a.x;Eigen::Matrix2d H;H<<a.u.squaredNorm(),a.u.dot(a.v)-r.dot(a.uv),a.u.dot(a.v)-r.dot(a.uv),a.v.squaredNorm();
        Real scale=std::max(a.u.squaredNorm(),a.v.squaredNorm());if(!std::isfinite(scale)||scale<=0)throw std::runtime_error("Crack contact surface has no geometric rank");
        Param gradient(a.u.dot(r),a.v.dot(r));
        Real roundoff=64*std::numeric_limits<Real>::epsilon()*coordinate_scale;
        if(std::abs(gradient[0])<=roundoff*a.u.norm()&&std::abs(gradient[1])<=roundoff*a.v.norm())return p;
        Eigen::SelfAdjointEigenSolver<Eigen::Matrix2d> e(H);
        if(e.info()!=Eigen::Success)throw std::runtime_error("Crack projection factorization");
        if(e.eigenvalues().minCoeff()<=1e-13*scale){H<<a.u.squaredNorm(),a.u.dot(a.v),a.u.dot(a.v),a.v.squaredNorm();}
        if(H.determinant()<=1e-14*scale*scale)throw std::runtime_error("Crack contact patch is geometrically unresolved");
        Param step=H.ldlt().solve(gradient);if(!step.allFinite())throw std::runtime_error("Nonfinite closest crack projection");
        if(step.norm()<2e-12)return p;
        Real alpha=1;bool accepted=false;
        for(size_t line=0;line<25;++line){Param d=alpha*step;
            // Exact difference of squared distances in the local bilinear
            // increment. Comparing two nearly equal |r|^2 values stalls
            // thin rotated faces before their tangential Newton correction.
            Vecd dx=d[0]*a.u+d[1]*a.v+d[0]*d[1]*a.uv;
            Real change=dx.squaredNorm()-2*r.dot(dx);
            Real rounding=64*std::numeric_limits<Real>::epsilon()*(dx.squaredNorm()+2*r.cwiseAbs().dot(dx.cwiseAbs()));
            if(change<=rounding){p+=d;accepted=true;break;}alpha*=.5;}
        if(!accepted){if(gradient.norm()<1e-12*scale)return p;throw std::runtime_error("Closest crack projection failed its descent check");}
    }
    std::ostringstream message;message<<std::setprecision(17)<<"Closest crack-surface projection did not converge; local_q "<<q.transpose()<<" param "<<p.transpose()<<" corners ";for(const auto&v:x)message<<v.transpose()<<" / ";throw std::runtime_error(message.str());
}

void CrackSurfaceContactNetwork::apply(){
    energy=max_penetration=normal_force_sum=force_residual=torque_residual=power_residual=tangential_force_sum=0;active_contacts=outside_projections=neighbor_contacts=0;for(auto&r:rates)r.setZero();
    struct Master {std::array<Vecd,4>x,v;};StdVec<Master>masters(sides.size());
    using Cell=std::array<long long,3>;struct Hash{size_t operator()(const Cell&a)const{size_t h=0;for(auto v:a)h^=std::hash<long long>{}(v)+0x9e3779b9+(h<<6)+(h>>2);return h;}};
    std::unordered_map<Cell,StdVec<size_t>,Hash>grid;Real cell_size=2*spacing;
    auto cell=[&](const Vecd&x){return Cell{(long long)std::floor(x[0]/cell_size),(long long)std::floor(x[1]/cell_size),(long long)std::floor(x[2]/cell_size)};};
    for(size_t m=0;m<sides.size();++m){auto&s=*sides[m].surface;s.ensureCurrent();Vecd lo=Vecd::Constant(std::numeric_limits<Real>::infinity()),hi=-lo;
        for(size_t a=0;a<4;++a){s.kinematics(a,masters[m].x[a],masters[m].v[a]);lo=lo.cwiseMin(masters[m].x[a]);hi=hi.cwiseMax(masters[m].x[a]);}Real margin=std::max(2*spacing,(hi-lo).norm());Cell a=cell(lo-Vecd::Constant(margin)),b=cell(hi+Vecd::Constant(margin));
        for(long long i=a[0];i<=b[0];++i)for(long long j=a[1];j<=b[1];++j)for(long long k=a[2];k<=b[2];++k)grid[{i,j,k}].push_back(m);
    }
    // Use the grid and the existing master order on both execution paths.
    // Projection is independent of contact force assembly and trace scattering.
    size_t next_projection=0;
    if(projection_batch_){
        projection_inputs_.clear();
        for(size_t si=0;si<sides.size();++si){
            const auto&slave=*sides[si].surface;
            for(size_t q=0;q<slave.natural.size();++q){
                Vecd x,v;
                slave.kinematics(4+q,x,v);
                const auto it=grid.find(cell(x));
                if(it==grid.end())continue;
                for(size_t m:it->second){
                    if(m==si)continue;
                    CrackProjectionInput input;
                    for(size_t d=0;d<3;++d){
                        for(size_t a=0;a<4;++a)input.coordinates[a][d]=masters[m].x[a][d];
                        input.coordinates[4][d]=x[d];
                    }
                    projection_inputs_.push_back(input);
                }
            }
        }
        projection_batch_->project(projection_inputs_,projection_results_);
        if(projection_results_.size()!=projection_inputs_.size())
            throw std::logic_error("Crack projection backend returned an incomplete batch");
    }
    auto closest=[&](size_t master,const Vecd&point){
        if(projection_batch_){
            const auto&p=projection_results_.at(next_projection++);
            if(p.status==0&&std::isfinite(p.s)&&std::isfinite(p.t))
                return BilinearCrackGeometry::Param(p.s,p.t);
            // The scalar solver handles exceptional geometry and preserves its diagnostics.
        }
        return BilinearCrackGeometry::closest(masters[master].x,point);
    };
    Vecd balance=Vecd::Zero(),moment=Vecd::Zero();Real magnitude=0,length=spacing,power=0,power_scale=0,energy_rate=0;
    for(size_t si=0;si<sides.size();++si){auto&slave=*sides[si].surface;for(size_t q=0;q<slave.natural.size();++q){Vecd x,v;slave.kinematics(4+q,x,v);Vecd tangentU=slave.position(slave.normalU(q))-x,tangentV=slave.position(slave.normalV(q))-x;Vecd outward=sides[si].outward_sign*tangentU.cross(tangentV).normalized();
            if(!outward.allFinite())throw std::runtime_error("Slave crack normal is unresolved");auto it=grid.find(cell(x));if(it==grid.end())continue;
            size_t chosen=sides.size();Real distance=std::numeric_limits<Real>::infinity(),gap=0;BilinearCrackGeometry::Value selected;Vecd normal;
            for(size_t m:it->second){if(m==si)continue;auto p=closest(m,x);if(p.minCoeff()< -1e-9||p.maxCoeff()>1+1e-9){++outside_projections;continue;}auto a=BilinearCrackGeometry::evaluate(masters[m].x,p);Vecd n=sides[m].outward_sign*a.u.cross(a.v).normalized();if(!n.allFinite())throw std::runtime_error("Master crack normal is unresolved");if(n.dot(outward)>=-1e-8)continue;Real g=n.dot(x-a.x);if(std::abs(g)<distance){distance=std::abs(g);chosen=m;gap=g;selected=a;normal=n;}}
            if(chosen==sides.size())continue;auto&master=*sides[chosen].surface;const auto&a=selected;Real stiffness=.5*penalty*sides[si].areas->at(q);Vecd residual=Vecd::Zero();if(gap<0)residual=gap*normal;
            Eigen::Matrix2d Hp;Hp<<a.u.squaredNorm(),a.u.dot(a.v)-residual.dot(a.uv),a.u.dot(a.v)-residual.dot(a.uv),a.v.squaredNorm();Eigen::SelfAdjointEigenSolver<Eigen::Matrix2d>e(Hp);if(e.info()!=Eigen::Success||e.eigenvalues().minCoeff()<=1e-13*e.eigenvalues().maxCoeff())throw std::runtime_error("Closest crack contact tangent is not a local minimum");
            Eigen::Matrix<Real,15,2>coupling;std::array<Real,5>c{1,-a.N[0],-a.N[1],-a.N[2],-a.N[3]};
            for(size_t i=0;i<5;++i){Vecd u=-c[i]*a.u,vv=-c[i]*a.v;if(i){u-=a.Nu[i-1]*residual;vv-=a.Nv[i-1]*residual;}coupling.template block<3,1>(3*i,0)=u;coupling.template block<3,1>(3*i,1)=vv;}
            Eigen::Matrix<Real,15,15>H=-coupling*Hp.inverse()*coupling.transpose();for(size_t i=0;i<5;++i)for(size_t j=0;j<5;++j)H.template block<3,3>(3*i,3*j)+=c[i]*c[j]*Matd::Identity();H*=stiffness;
            std::array<const ReferenceSolidSurface::Row*,5>trace{&slave.rows[4+q],&master.rows[0],&master.rows[1],&master.rows[2],&master.rows[3]};std::array<Real,5>weighted{};
            for(size_t i=0;i<5;++i)for(size_t j=0;j<trace[i]->index.size();++j)weighted[i]+=std::abs(trace[i]->weight[j])/std::sqrt(slave.port.mass[trace[i]->index[j]]);
            for(size_t i=0;i<5;++i)for(size_t d=0;d<3;++d){Real row=0;for(size_t j=0;j<15;++j)row+=std::abs(H(3*i+d,j))*weighted[j/3];for(size_t j=0;j<trace[i]->index.size();++j)rates[trace[i]->index[j]][d]+=std::abs(trace[i]->weight[j])*row/std::sqrt(slave.port.mass[trace[i]->index[j]]);}
            if(gap>=0)continue;if(sides[si].reference_patch!=sides[chosen].reference_patch)++neighbor_contacts;Vecd f=-stiffness*gap*normal;slave.scatter(4+q,f);Vecd relative=v;
            balance+=f;moment+=x.cross(f);magnitude+=f.norm();length=std::max(length,x.norm());power+=f.dot(v);power_scale+=std::abs(f.dot(v));
            for(size_t j=0;j<4;++j){Vecd reaction=-a.N[j]*f;master.scatter(j,reaction);balance+=reaction;moment+=masters[chosen].x[j].cross(reaction);magnitude+=reaction.norm();power+=reaction.dot(masters[chosen].v[j]);power_scale+=std::abs(reaction.dot(masters[chosen].v[j]));relative-=a.N[j]*masters[chosen].v[j];}
            energy_rate+=stiffness*gap*normal.dot(relative);energy+=.5*stiffness*gap*gap;max_penetration=std::max(max_penetration,-gap);normal_force_sum+=f.norm();tangential_force_sum+=(f-normal*f.dot(normal)).norm();++active_contacts;
        }}
    Real bound=0;for(const auto&r:rates)bound=std::max(bound,r.maxCoeff());stable_dt=bound>0?.2/std::sqrt(bound):std::numeric_limits<Real>::infinity();force_residual=balance.norm()/std::max(magnitude,Real(1e-30));torque_residual=moment.norm()/std::max(magnitude*length,Real(1e-30));power_residual=std::abs(power+energy_rate)/std::max(power_scale,Real(1e-30));
}
}
