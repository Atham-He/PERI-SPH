/* General 3-D density-dependent Drucker-Prager/Cap material-point kernel.
 * SI units. Pressure is positive in compression; yield stresses are Cauchy.
 * Finite multiplicative plasticity with isotropic elastic Hencky strain.
 * K=1, transition alpha=0. Associated cap; nonassociated shear potential
 * g_s=eta*b^2(p-pa)^2+(q+c)^2, c=b*pa, 0<=eta<=1.
 * eta=1 recovers the prior shifted ellipse; eta=0 suppresses SHEAR dilation
 * while retaining associated cap compaction. This is a controlled assumption,
 * not a new MCC measurement. See Frenning, Powder Technol. 423, 118503.
 * Optional H>=0 adds d_eff(r,alpha)=d_table(r)+H*alpha, with accumulated
 * equivalent deviatoric plastic strain alpha. The hydrostatic cap tip pb(r)
 * is unchanged. This is a phenomenological extension, not a measured MCC law.
 * The shift is derived, not a fitted dilation parameter. Density-dependent
 * elasticity is an effective extension; plastic stress work below is NOT a
 * proof of a complete thermodynamic balance with evolving elastic moduli.
 */
#ifndef SPH_DRUCKER_PRAGER_CAP_KERNEL_H
#define SPH_DRUCKER_PRAGER_CAP_KERNEL_H
#include "dpc_point_batch.h"
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <stdexcept>
#include <vector>

namespace SPH { namespace dpc {
using Tensor = Eigen::Matrix3d;
struct Parameters : Node
{
    double K, G, b, pa, c;
};
struct State
{
    Tensor inverse_Cp = Tensor::Identity();
    double kappa = 0;
    double alpha = 0;
};
struct Result
{
    State state;
    Tensor sigma = Tensor::Zero();
    double p = 0, q = 0, relative_density = 0, plastic_shear = 0;
    double plastic_stress_work = 0, elastic_energy = 0, yield_residual = 0;
    double effective_d = 0;
    int surface = 0; // 0 elastic, 1 shear, 2 cap, 3 tensile apex
};

// This view owns no memory. The CPU model and the CUDA batch both use the
// same table interpolation, trial state, branch selection and admissibility.
class PointModel
{
    const Node *nodes_;
    size_t count_;
    double r0_, shear_dilation_, shear_hardening_;

    EIGEN_DEVICE_FUNC static bool finite(double x) { return std::isfinite(x); }

    template <int rows, int cols>
    EIGEN_DEVICE_FUNC static bool allFinite(const Eigen::Matrix<double, rows, cols> &a)
    {
        for (int i = 0; i < rows * cols; ++i)
            if (!finite(a.data()[i]))
                return false;
        return true;
    }

    EIGEN_DEVICE_FUNC static Eigen::Vector2d solveNewtonSystem(const Eigen::Matrix2d &jac, const Eigen::Vector2d &f)
    {
#if defined(__CUDA_ARCH__)
        // Eigen's dynamic LU entry is host-only. Apply the same partial pivot
        // and triangular solve directly for this fixed two-invariant system.
        const int pivot = std::abs(jac(0,0)) < std::abs(jac(1,0)) ? 1 : 0;
        const int other = 1 - pivot;
        const double multiplier = jac(other,0)/jac(pivot,0);
        const double upper = jac(other,1)-multiplier*jac(pivot,1);
        const double second = (-f[other]+multiplier*f[pivot])/upper;
        return {(-f[pivot]-jac(pivot,1)*second)/jac(pivot,0), second};
#else
        return jac.partialPivLu().solve(-f);
#endif
    }

  public:
    EIGEN_DEVICE_FUNC explicit PointModel(PointInputModel input)
        : nodes_(input.nodes), count_(input.count), r0_(input.initial_density),
          shear_dilation_(input.shear_dilation), shear_hardening_(input.shear_hardening) {}

    EIGEN_DEVICE_FUNC double kappaMin() const { return std::log(nodes_[0].r/r0_); }
    EIGEN_DEVICE_FUNC double kappaMax() const { return std::log(nodes_[count_-1].r/r0_); }
    EIGEN_DEVICE_FUNC PointFailure parameters(double r, double alpha, Parameters &p) const
    {
        if (!finite(r) || r<nodes_[0].r-1e-12 || r>nodes_[count_ - 1].r+1e-12)
            return PointFailure::DensityDomain;
        r=std::max(nodes_[0].r,std::min(nodes_[count_ - 1].r,r));
        size_t k=1; while(k+1<count_ && nodes_[k].r<r) ++k;
        const auto &a=nodes_[k-1]; const auto &z=nodes_[k];
        const double t=(r-a.r)/(z.r-a.r);
        auto mix=[&](double x,double y){return x+t*(y-x);};
        p.r=r; p.E=mix(a.E,z.E); p.nu=mix(a.nu,z.nu); p.d=mix(a.d,z.d);
        p.beta=mix(a.beta,z.beta); p.R=mix(a.R,z.R); p.pb=mix(a.pb,z.pb);
        // Retain the exact old floating-point path when the option is absent.
        if(shear_hardening_>0)
        {
            if(!finite(alpha) || alpha<0)return PointFailure::ShearHistory;
            p.d+=shear_hardening_*alpha;
            if(!finite(p.d) || p.pb<=p.R*p.d)
                return PointFailure::HardeningDomain;
        }
        p.K=p.E/(3*(1-2*p.nu)); p.G=p.E/(2*(1+p.nu)); p.b=std::tan(p.beta);
        p.pa=(p.pb-p.R*p.d)/(1+p.R*p.b); p.c=p.b*p.pa;
        return PointFailure::None;
    }
    EIGEN_DEVICE_FUNC static double yield(double p,double q,const Parameters &a)
    {
        return p<=a.pa ? q-a.d-p*a.b : std::hypot(p-a.pa,a.R*q)-a.R*(a.d+a.pa*a.b);
    }
    EIGEN_DEVICE_FUNC PointFailure integrate(const Tensor &F, const State &old, Result &out) const
    {
        const double J=F.determinant();
        if (!allFinite(F) || J<=0 || !allFinite(old.inverse_Cp) || !finite(old.kappa) || !finite(old.alpha) || old.alpha<0)
            return PointFailure::Deformation;
        const Tensor be=F*old.inverse_Cp*F.transpose();
        Eigen::SelfAdjointEigenSolver<Tensor> eig;
#if defined(__CUDA_ARCH__)
        // Eigen's fixed-size direct solver is device-callable. The iterative
        // solver uses a host-only tridiagonalization specialization.
        eig.computeDirect(.5*(be+be.transpose()));
#else
        eig.compute(.5*(be+be.transpose()));
#endif
        if(eig.info()!=Eigen::Success || eig.eigenvalues().minCoeff()<=1e-14)
            return PointFailure::ElasticMetric;
        const Eigen::Vector3d ev=.5*eig.eigenvalues().array().log().matrix();
        const double etv=ev.sum();
        const Eigen::Vector3d dev=(ev.array()-etv/3).matrix();
        const double eq=std::sqrt(2./3.)*dev.norm();
        auto stress=[&](double kap,double gam,double &p,double &q,Parameters &a)
        {
            const auto failure = parameters(r0_*std::exp(kap), old.alpha+gam, a);
            if (failure == PointFailure::None)
            {
                p=-a.K/J*(etv+kap-old.kappa);
                q=3*a.G/J*(eq-gam);
            }
            return failure;
        };
        double kap=old.kappa,gam=0,p,q; Parameters a;
        auto failure = stress(kap,gam,p,q,a);
        if (failure != PointFailure::None)
            return failure;
        int surface=0;
        const double scale=std::max(1.,a.pb);
        if(yield(p,q,a)>1e-9*scale)
        {
            // A two-variable backward-Euler invariant return. The potential
            // scale cancels in the flow ratio, avoiding an extra multiplier.
            // A converged root on the continuation of one branch is a
            // predictor for the other branch when hardening moves the junction.
            // Reuse only that converged root; all branch/yield guards remain.
            Eigen::Vector2d branch_seed(old.kappa,0);
            auto attempt=[&](int branch)->bool
            {
                Eigen::Vector2d x=branch_seed;
                auto residual=[&](const Eigen::Vector2d &v)->Eigen::Vector2d
                {
                    double pp,qq;
                    Parameters z;
                    failure = stress(v[0],v[1],pp,qq,z);
                    if (failure != PointFailure::None)
                        return Eigen::Vector2d::Zero();
                    const double f=branch==1 ? qq-z.d-pp*z.b : std::hypot(pp-z.pa,z.R*qq)-z.R*(z.d+z.pa*z.b);
                    const double flow=branch==1 ?
                        (v[0]-old.kappa)*(qq+z.c)-shear_dilation_*v[1]*(pp-z.pa)*z.b*z.b :
                        (v[0]-old.kappa)*z.R*z.R*qq-v[1]*(pp-z.pa);
                    return Eigen::Vector2d(f/scale,flow/scale);
                };
                auto admissible=[&](const Eigen::Vector2d &v)
                {
                    if(!(allFinite(v) && v[0]>=kappaMin()-1e-13 && v[0]<=kappaMax()+1e-13 && v[1]>=-1e-13 && v[1]<=eq+1e-13))return false;
                    if(shear_hardening_>0)
                    {
                        Parameters unused;
                        if (v[1]<0 || parameters(r0_*std::exp(v[0]), old.alpha+v[1], unused) != PointFailure::None)
                            return false;
                    }
                    return true;
                };
                for(int it=0;it<35;++it)
                {
                    const Eigen::Vector2d f=residual(x);
                    if (failure != PointFailure::None)
                        return false;
                    if(f.lpNorm<Eigen::Infinity>()<1e-9)
                    {
                        if(shear_hardening_>0)branch_seed=x;
                        double pp,qq;
                        Parameters z;
                        failure = stress(x[0],x[1],pp,qq,z);
                        if (failure != PointFailure::None)
                            return false;
                        if((branch==1 && pp>z.pa+1e-7*scale) || (branch==2 && pp<z.pa-1e-7*scale) ||
                           qq < -1e-7*scale || (branch==1 && x[0]>old.kappa+1e-10) || (branch==2 && x[0]<old.kappa-1e-10)) return false;
                        // At the moving shear/cap junction, validate the complete piecewise
                        // surface BEFORE accepting a branch. Otherwise a locally converged
                        // shear root can fail the unchanged final cap admissibility guard.
                        if(yield(pp,qq,z)>2e-8*scale)return false;
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
                        const Eigen::Vector2d perturbed = residual(y);
                        if (failure != PointFailure::None)
                            return false;
                        jac.col(c)=(perturbed-f)/h;
                    }
                    if(std::abs(jac.determinant())<1e-20)return false;
                    const Eigen::Vector2d dx = solveNewtonSystem(jac, f);
                    bool moved=false;
                    for(double t=1.;t>=1./16777216;t*=.5)
                    {
                        const Eigen::Vector2d y=x+t*dx;
                        if (!admissible(y))
                            continue;
                        const Eigen::Vector2d candidate = residual(y);
                        if (failure != PointFailure::None)
                            return false;
                        if(candidate.norm()<f.norm()){x=y;moved=true;break;}
                    }
                    if(!moved)return false;
                }
                return false;
            };
            // Hydrostatic compaction has an exact one-dimensional return.
            // Bracketing avoids a Newton step landing on the back of the cap
            // when an imposed volumetric increment is large.
            bool hydro_returned=false;
            if(eq<1e-14 && p>a.pb)
            {
                auto cap_tip = [&](double k) {
                    double pp=0,qq=0;
                    Parameters z{};
                    failure = stress(k,0,pp,qq,z);
                    return pp-z.pb;
                };
                double lo=old.kappa,hi=std::min(kappaMax(),old.kappa-etv);
                const double tip = cap_tip(hi);
                if (failure != PointFailure::None)
                    return failure;
                if(tip>1e-8*scale)
                    return PointFailure::HydrostaticBound;
                for(int it=0;it<70;++it)
                {
                    const double mid=.5*(lo+hi);
                    const double trial=cap_tip(mid);
                    if (failure != PointFailure::None)
                        return failure;
                    if(trial>0)lo=mid;else hi=mid;
                }
                kap=.5*(lo+hi);gam=0;surface=2;hydro_returned=true;
            }
            const int first=p<=a.pa ? 1:2;
            bool returned = hydro_returned || attempt(first);
            if (failure != PointFailure::None)
                return failure;
            if (!returned)
                returned = attempt(3-first);
            if (failure != PointFailure::None)
                return failure;
            if (!returned)
            {
                // At the tensile apex q=0 and the deviatoric subgradient is
                // multivalued. Solve the volumetric return and check its cone.
                auto apex = [&](double k) {
                    double pp=0,qq=0;
                    Parameters z{};
                    failure = stress(k,eq,pp,qq,z);
                    return failure == PointFailure::None ? pp+z.d/z.b : 0.;
                };
                double lo=kappaMin(),hi=old.kappa;
                const double lower=apex(lo);
                if (failure != PointFailure::None)
                    return failure;
                if (lower<0)
                    return PointFailure::ReturnDomain;
                const double upper=apex(hi);
                if (failure != PointFailure::None)
                    return failure;
                if (upper>0)
                    return PointFailure::ReturnDomain;
                for(int it=0;it<70;++it)
                {
                    const double mid=.5*(lo+hi);
                    const double trial=apex(mid);
                    if (failure != PointFailure::None)
                        return failure;
                    if(trial>0)lo=mid;else hi=mid;
                }
                kap=.5*(lo+hi);gam=eq;
                failure=stress(kap,gam,p,q,a);
                if (failure != PointFailure::None)
                    return failure;
                if(shear_dilation_==0)
                    return PointFailure::ZeroDilationApex;
                const double lambda=(kap-old.kappa)/(shear_dilation_*(p-a.pa)*a.b*a.b);
                if(lambda<0 || gam>lambda*a.c+1e-10)
                    return PointFailure::ApexFlow;
                surface=3;
            }
            failure=stress(kap,gam,p,q,a);
            if (failure != PointFailure::None)
                return failure;
        }
        const Eigen::Vector3d ec=(dev*(eq>1e-18 ? (eq-gam)/eq : 0)).array()+(etv+kap-old.kappa)/3;
        const Tensor elastic_b=eig.eigenvectors()*(2*ec.array()).exp().matrix().asDiagonal()*eig.eigenvectors().transpose();
        out.state.kappa=kap;out.state.alpha=old.alpha+gam;out.effective_d=a.d;
        const Tensor invF=F.inverse();out.state.inverse_Cp=invF*elastic_b*invF.transpose();
        out.state.inverse_Cp=.5*(out.state.inverse_Cp+out.state.inverse_Cp.transpose()).eval();
        const Eigen::Vector3d tau=(2*a.G*(ec.array()-ec.sum()/3)+a.K*ec.sum()).matrix();
        out.sigma=eig.eigenvectors()*tau.asDiagonal()*eig.eigenvectors().transpose()/J;
        out.p=p;out.q=std::max(0.,q);out.relative_density=a.r;out.plastic_shear=gam;out.surface=surface;
        out.plastic_stress_work=J*(p*(kap-old.kappa)+q*gam);
        out.elastic_energy=.5*a.K*ec.sum()*ec.sum()+a.G*(ec.array()-ec.sum()/3).matrix().squaredNorm();
        out.yield_residual=std::max(0.,yield(p,q,a))/scale;
        if(out.plastic_stress_work < -1e-6 || out.yield_residual>2e-8 || !allFinite(out.sigma))
        {
            out.q=q;
            return PointFailure::BadReturn;
        }
        return PointFailure::None;
    }
};

class Model
{
    std::vector<Node> nodes_;
    double r0_, shear_dilation_, shear_hardening_;
    static bool finite(double x) { return std::isfinite(x); }
  public:
    explicit Model(std::vector<Node> nodes, double r0, double shear_dilation=1, double shear_hardening=0)
        : nodes_(std::move(nodes)), r0_(r0), shear_dilation_(shear_dilation), shear_hardening_(shear_hardening)
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
        if(!finite(shear_dilation_) || shear_dilation_<0 || shear_dilation_>1) throw std::invalid_argument("DPC shear-dilation scale must be in [0,1]");
        if(!finite(shear_hardening_) || shear_hardening_<0)throw std::invalid_argument("DPC shear hardening must be finite and nonnegative");
        parameters(r0_);
    }
    static Model fromCSV(const std::string &path, double r0, double shear_dilation=1, double shear_hardening=0)
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
        return Model(nodes,r0,shear_dilation,shear_hardening);
    }
    double r0() const { return r0_; }
    double shearDilation() const { return shear_dilation_; }
    double shearHardening() const { return shear_hardening_; }
    double kappaMin() const { return std::log(nodes_.front().r/r0_); }
    double kappaMax() const { return std::log(nodes_.back().r/r0_); }
    const std::vector<Node> &nodes() const { return nodes_; }
    PointInputModel pointInput() const
    {
        return {nodes_.data(), nodes_.size(), r0_, shear_dilation_, shear_hardening_};
    }

    Parameters parameters(double r, double alpha=0) const
    {
        Parameters a;
        const auto failure = PointModel(pointInput()).parameters(r, alpha, a);
        if (failure != PointFailure::None)
            throw std::runtime_error(failureReason(failure));
        return a;
    }

    static double yield(double p,double q,const Parameters &a) { return PointModel::yield(p,q,a); }

    Result update(const Tensor &F, const State &old) const
    {
        Result result;
        const auto failure = PointModel(pointInput()).integrate(F, old, result);
        if (failure != PointFailure::None)
            throwPointError(failure, F, old, result);
        return result;
    }

    [[noreturn]] void throwPointError(PointFailure failure, const Tensor &F, const State &old, const Result &result) const
    {
        std::ostringstream reason;
        if (failure == PointFailure::BadReturn)
            reason<<std::setprecision(17)<<"BAD_RETURN p="<<result.p<<" q="<<result.q<<" dkappa="<<result.state.kappa-old.kappa
                  <<" gamma="<<result.plastic_shear<<" work="<<result.plastic_stress_work
                  <<" residual="<<result.yield_residual<<" surface="<<result.surface;
        else
            reason<<failureReason(failure);
        // Preserve an exact local failure state for offline replay.
        std::ostringstream out;
        out<<std::setprecision(17)<<"DPC_FAILURE_STATE {\"reason\":"<<std::quoted(reason.str())
           <<",\"r0\":"<<r0_<<",\"eta\":"<<shear_dilation_<<",\"kappa\":"<<old.kappa
           <<",\"alpha\":"<<old.alpha<<",\"shear_hardening_Pa\":"<<shear_hardening_<<",\"F\":[";
        for(int r=0;r<3;++r)for(int c=0;c<3;++c)out<<(r||c ? ",":"")<<F(r,c);
        out<<"],\"inverse_Cp\":[";
        for(int r=0;r<3;++r)for(int c=0;c<3;++c)out<<(r||c ? ",":"")<<old.inverse_Cp(r,c);
        out<<"]}";
        throw std::runtime_error(out.str());
    }

  private:
    static const char *failureReason(PointFailure failure)
    {
        switch (failure)
        {
        case PointFailure::DensityDomain: return "DPC relative density outside literature table; extrapolation forbidden";
        case PointFailure::ShearHistory: return "Invalid accumulated DPC plastic shear";
        case PointFailure::HardeningDomain: return "DPC shear hardening exceeds the positive-cap-transition domain; no clipping";
        case PointFailure::Deformation: return "invalid DPC deformation/history";
        case PointFailure::ElasticMetric: return "non-positive DPC elastic metric";
        case PointFailure::HydrostaticBound: return "hydrostatic cap return requires density above the table upper bound";
        case PointFailure::ReturnDomain: return "DPC return failed or density domain exhausted; no clamping/extrapolation";
        case PointFailure::ZeroDilationApex: return "zero-dilation shear cannot return a hydrostatic tensile apex; no extra tensile flow assumed";
        case PointFailure::ApexFlow: return "DPC tensile apex flow condition failed";
        default: return "DPC material-point return failed";
        }
    }
};
}} // namespace SPH::dpc
#endif
