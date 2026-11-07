/* General 3-D density-dependent Drucker-Prager/Cap material-point kernel.
 * SI units. Pressure is positive in compression; yield stresses are Cauchy.
 * Finite multiplicative plasticity with isotropic elastic Hencky strain.
 * K=1, transition alpha=0. Associated cap; nonassociated shear potential
 * g_s=b^2(p-pa)^2+(q+c)^2, c=b*pa (Frenning, Powder Technol. 423, 118503).
 * The shift is derived, not a fitted dilation parameter. Density-dependent
 * elasticity is an effective extension; plastic stress work below is NOT a
 * proof of a complete thermodynamic balance with evolving elastic moduli.
 */
#ifndef SPH_DRUCKER_PRAGER_CAP_KERNEL_H
#define SPH_DRUCKER_PRAGER_CAP_KERNEL_H
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#include <algorithm>
#include <cmath>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <vector>

namespace SPH { namespace dpc {
using Tensor = Eigen::Matrix3d;
struct Node { double r, E, nu, d, beta, R, pb; };
struct Parameters : Node
{
    double K, G, b, pa, c;
};
struct State
{
    Tensor inverse_Cp = Tensor::Identity();
    double kappa = 0;
};
struct Result
{
    State state;
    Tensor sigma = Tensor::Zero();
    double p = 0, q = 0, relative_density = 0, plastic_shear = 0;
    double plastic_stress_work = 0, elastic_energy = 0, yield_residual = 0;
    int surface = 0; // 0 elastic, 1 shear, 2 cap, 3 tensile apex
};

class Model
{
    std::vector<Node> nodes_;
    double r0_;
    static bool finite(double x) { return std::isfinite(x); }
  public:
    explicit Model(std::vector<Node> nodes, double r0) : nodes_(std::move(nodes)), r0_(r0)
    {
        if (nodes_.size() < 2 || !finite(r0_)) throw std::invalid_argument("DPC requires a density table and initial relative density");
        for (size_t i=0; i<nodes_.size(); ++i)
        {
            const auto &a=nodes_[i];
            if (!(finite(a.r)&&finite(a.E)&&finite(a.nu)&&finite(a.d)&&finite(a.beta)&&finite(a.R)&&finite(a.pb)) ||
                a.r<=0 || a.r>=1 || a.E<=0 || a.nu<0 || a.nu>=.49 || a.d<=0 || a.beta<=0 || a.beta>=1.55 || a.R<=0 || a.pb<=a.R*a.d ||
                (i && (a.r<=nodes_[i-1].r || a.pb<nodes_[i-1].pb)))
                throw std::invalid_argument("invalid or nonmonotone DPC density table (beta in radians, pressure in Pa)");
        }
        parameters(r0_);
    }
    static Model fromCSV(const std::string &path, double r0)
    {
        std::ifstream f(path); if (!f) throw std::invalid_argument("cannot read DPC density table: "+path);
        std::vector<Node> nodes;
        for (std::string line; std::getline(f,line);)
        {
            line=line.substr(0,line.find('#'));
            if (line.find_first_not_of(" \t\r")==std::string::npos) continue;
            std::replace(line.begin(),line.end(),',',' ');
            std::istringstream in(line); Node a; std::string extra;
            if (!(in>>a.r>>a.E>>a.nu>>a.d>>a.beta>>a.R>>a.pb) || (in>>extra))
                throw std::invalid_argument("DPC CSV expects r,E_Pa,nu,d_Pa,beta_rad,R,pb_Pa");
            nodes.push_back(a);
        }
        return Model(nodes,r0);
    }
    double r0() const { return r0_; }
    double kappaMin() const { return std::log(nodes_.front().r/r0_); }
    double kappaMax() const { return std::log(nodes_.back().r/r0_); }
    const std::vector<Node> &nodes() const { return nodes_; }
    Parameters parameters(double r) const
    {
        if (!finite(r) || r<nodes_.front().r-1e-12 || r>nodes_.back().r+1e-12)
            throw std::runtime_error("DPC relative density outside literature table; extrapolation forbidden");
        r=std::max(nodes_.front().r,std::min(nodes_.back().r,r));
        size_t k=1; while(k+1<nodes_.size() && nodes_[k].r<r) ++k;
        const auto &a=nodes_[k-1]; const auto &z=nodes_[k];
        const double t=(r-a.r)/(z.r-a.r);
        auto mix=[&](double x,double y){return x+t*(y-x);};
        Parameters p; p.r=r; p.E=mix(a.E,z.E); p.nu=mix(a.nu,z.nu); p.d=mix(a.d,z.d);
        p.beta=mix(a.beta,z.beta); p.R=mix(a.R,z.R); p.pb=mix(a.pb,z.pb);
        p.K=p.E/(3*(1-2*p.nu)); p.G=p.E/(2*(1+p.nu)); p.b=std::tan(p.beta);
        p.pa=(p.pb-p.R*p.d)/(1+p.R*p.b); p.c=p.b*p.pa;
        return p;
    }
    static double yield(double p,double q,const Parameters &a)
    {
        return p<=a.pa ? q-a.d-p*a.b : std::hypot(p-a.pa,a.R*q)-a.R*(a.d+a.pa*a.b);
    }
    Result update(const Tensor &F, const State &old) const
    {
        const double J=F.determinant();
        if (!F.allFinite() || J<=0 || !old.inverse_Cp.allFinite() || !finite(old.kappa))
            throw std::runtime_error("invalid DPC deformation/history");
        const Tensor be=F*old.inverse_Cp*F.transpose();
        Eigen::SelfAdjointEigenSolver<Tensor> eig(.5*(be+be.transpose()));
        if(eig.info()!=Eigen::Success || eig.eigenvalues().minCoeff()<=1e-14)
            throw std::runtime_error("non-positive DPC elastic metric");
        const Eigen::Vector3d ev=.5*eig.eigenvalues().array().log().matrix();
        const double etv=ev.sum();
        const Eigen::Vector3d dev=(ev.array()-etv/3).matrix();
        const double eq=std::sqrt(2./3.)*dev.norm();
        auto stress=[&](double kap,double gam,double &p,double &q,Parameters &a)
        {
            a=parameters(r0_*std::exp(kap));
            p=-a.K/J*(etv+kap-old.kappa); q=3*a.G/J*(eq-gam);
        };
        double kap=old.kappa,gam=0,p,q; Parameters a;
        stress(kap,gam,p,q,a);
        int surface=0;
        const double scale=std::max(1.,a.pb);
        if(yield(p,q,a)>1e-9*scale)
        {
            // A two-variable backward-Euler invariant return. The potential
            // scale cancels in the flow ratio, avoiding an extra multiplier.
            auto attempt=[&](int branch)->bool
            {
                Eigen::Vector2d x(old.kappa,0);
                auto residual=[&](const Eigen::Vector2d &v)->Eigen::Vector2d
                {
                    double pp,qq; Parameters z; stress(v[0],v[1],pp,qq,z);
                    const double f=branch==1 ? qq-z.d-pp*z.b : std::hypot(pp-z.pa,z.R*qq)-z.R*(z.d+z.pa*z.b);
                    const double flow=branch==1 ?
                        (v[0]-old.kappa)*(qq+z.c)-v[1]*(pp-z.pa)*z.b*z.b :
                        (v[0]-old.kappa)*z.R*z.R*qq-v[1]*(pp-z.pa);
                    return Eigen::Vector2d(f/scale,flow/scale);
                };
                auto admissible=[&](const Eigen::Vector2d &v)
                {return v.allFinite() && v[0]>=kappaMin()-1e-13 && v[0]<=kappaMax()+1e-13 && v[1]>=-1e-13 && v[1]<=eq+1e-13;};
                for(int it=0;it<35;++it)
                {
                    const Eigen::Vector2d f=residual(x);
                    if(f.lpNorm<Eigen::Infinity>()<1e-9)
                    {
                        double pp,qq; Parameters z; stress(x[0],x[1],pp,qq,z);
                        if((branch==1 && pp>z.pa+1e-7*scale) || (branch==2 && pp<z.pa-1e-7*scale) ||
                           qq < -1e-7*scale || (branch==1 && x[0]>old.kappa+1e-10) || (branch==2 && x[0]<old.kappa-1e-10)) return false;
                        kap=x[0];gam=std::max(0.,std::min(eq,x[1]));surface=branch;return true;
                    }
                    Eigen::Matrix2d jac;
                    for(int c=0;c<2;++c)
                    {
                        double h=1e-7; Eigen::Vector2d y=x;y[c]+=h;
                        // The derivative may sample the smooth extension at
                        // q<0. Only accepted iterates obey gamma<=eq. This is
                        // essential at a hydrostatic state, where eq is zero.
                        if(c==0 && (y[0]<kappaMin() || y[0]>kappaMax())){h=-h;y=x;y[c]+=h;}
                        if(c==0 && (y[0]<kappaMin() || y[0]>kappaMax()))return false;
                        jac.col(c)=(residual(y)-f)/h;
                    }
                    if(std::abs(jac.determinant())<1e-20)return false;
                    Eigen::Vector2d dx=jac.partialPivLu().solve(-f);
                    bool moved=false;
                    for(double t=1.;t>=1./16777216;t*=.5)
                    {
                        const Eigen::Vector2d y=x+t*dx;
                        if(admissible(y) && residual(y).norm()<f.norm()){x=y;moved=true;break;}
                    }
                    if(!moved)return false;
                }
                return false;
            };
            const int first=p<=a.pa ? 1:2;
            if(!attempt(first) && !attempt(3-first))
            {
                // At the tensile apex, q=0 and the deviatoric subgradient is
                // multivalued. Solve the volumetric return and check its cone.
                auto apex=[&](double k){double pp,qq;Parameters z;stress(k,eq,pp,qq,z);return pp+z.d/z.b;};
                double lo=kappaMin(),hi=old.kappa;
                if(apex(lo)<0 || apex(hi)>0) throw std::runtime_error("DPC return failed or density domain exhausted; no clamping/extrapolation");
                for(int it=0;it<70;++it){double mid=.5*(lo+hi);if(apex(mid)>0)lo=mid;else hi=mid;}
                kap=.5*(lo+hi);gam=eq;stress(kap,gam,p,q,a);
                const double lambda=(kap-old.kappa)/((p-a.pa)*a.b*a.b);
                if(lambda<0 || gam>lambda*a.c+1e-10)throw std::runtime_error("DPC tensile apex flow condition failed");
                surface=3;
            }
            stress(kap,gam,p,q,a);
        }
        const Eigen::Vector3d ec=(dev*(eq>1e-18 ? (eq-gam)/eq : 0)).array()+(etv+kap-old.kappa)/3;
        const Tensor elastic_b=eig.eigenvectors()*(2*ec.array()).exp().matrix().asDiagonal()*eig.eigenvectors().transpose();
        Result out;out.state.kappa=kap;
        const Tensor invF=F.inverse();out.state.inverse_Cp=invF*elastic_b*invF.transpose();
        out.state.inverse_Cp=.5*(out.state.inverse_Cp+out.state.inverse_Cp.transpose()).eval();
        const Eigen::Vector3d tau=(2*a.G*(ec.array()-ec.sum()/3)+a.K*ec.sum()).matrix();
        out.sigma=eig.eigenvectors()*tau.asDiagonal()*eig.eigenvectors().transpose()/J;
        out.p=p;out.q=std::max(0.,q);out.relative_density=a.r;out.plastic_shear=gam;out.surface=surface;
        out.plastic_stress_work=J*(p*(kap-old.kappa)+q*gam);
        out.elastic_energy=.5*a.K*ec.sum()*ec.sum()+a.G*(ec.array()-ec.sum()/3).matrix().squaredNorm();
        out.yield_residual=std::max(0.,yield(p,q,a))/scale;
        if(out.plastic_stress_work < -1e-6 || out.yield_residual>2e-8 || !out.sigma.allFinite())
            throw std::runtime_error("DPC admissibility/plastic stress-work check failed");
        return out;
    }
};
}} // namespace SPH::dpc
#endif
