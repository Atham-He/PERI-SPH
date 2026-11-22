/* Licensed under the Apache License, Version 2.0. */
#ifndef PD_SOLID_REFERENCE_COUPLING_H
#define PD_SOLID_REFERENCE_COUPLING_H

#include "base_particles.hpp"
#include "pd_coulomb_friction.h"
#include "pd_bonded_interface_batch.h"
#include "base_body.h"
#include <Eigen/Cholesky>
#include <Eigen/Eigenvalues>
#include <array>
#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <unordered_map>

namespace SPH::peridynamics
{
/** A force port, not a constitutive adapter. PD stores external acceleration;
 * solid SPH stores external force. Both ports receive physical forces in N.
 * The caller owns clearing and the common integration time layer. */
class SolidForcePort
{
  public:
    enum class Storage { ForcePrior, PDAccelerationPrior };
    RealBody &body;
    size_t count;
    Vecd *position, *velocity, *external;
    Real *mass, *volume;
    StdVec<Vecd> reference;
    Storage storage;
    SolidForcePort(RealBody &b, Storage s) : body(b), storage(s)
    {
        auto &p=b.getBaseParticles();count=p.TotalRealParticles();
        position=p.getVariableDataByName<Vecd>("Position");
        velocity=p.getVariableDataByName<Vecd>("Velocity");
        mass=p.getVariableDataByName<Real>("Mass");
        volume=p.getVariableDataByName<Real>("VolumetricMeasure");
        external=p.registerStateVariableData<Vecd>(s==Storage::ForcePrior ? "ForcePrior" : "PDAccelerationPrior");
        reference.assign(position,position+count);
        for(size_t i=0;i<count;++i)
            if(!std::isfinite(mass[i]) || mass[i]<=0 || !std::isfinite(volume[i]) || volume[i]<=0)
                throw std::invalid_argument("Solid coupling requires positive physical masses and volumes");
    }
    void clear() {for(size_t i=0;i<count;++i) external[i].setZero();}
    void addForce(size_t i,const Vecd &force)
    {external[i] += storage==Storage::ForcePrior ? force : force/mass[i];}
};

/** Fixed-reference polynomial-reproducing MLS surface evaluation. It supports
 * nonuniform particle volumes, including explicitly resolved thickness layers.
 * It supplies no bending law and does not turn a volume body into a shell. */
class ReferenceSolidSurface
{
  public:
    struct Row {StdVec<size_t> index;StdVec<Real> weight;};
    SolidForcePort &port;
    StdVec<Vecd> points;
    StdVec<Row> rows;
    Real affine_error=0, polynomial_error=0, largest_condition_number=0;
    virtual ~ReferenceSolidSurface()=default;
    virtual void ensureCurrent()const{}
    virtual size_t mappingVersion()const{return 0;}
    const unsigned reproduction_order;
  protected:
    struct DeferredBuildTag {};
    ReferenceSolidSurface(SolidForcePort &p,const StdVec<Vecd>&targets,unsigned order,DeferredBuildTag)
        :port(p),points(targets),reproduction_order(order){}
  public:
    ReferenceSolidSurface(SolidForcePort &p,const StdVec<Vecd> &targets,Real spacing,Real support_ratio=2.5,unsigned order=1)
        : ReferenceSolidSurface(p,targets,order,DeferredBuildTag{})
    {
        static_assert(Dimensions==3,"ReferenceSolidSurface currently implements three-dimensional affine reconstruction");
        if(spacing<=0 || support_ratio<=1 || !std::isfinite(spacing*support_ratio))
            throw std::invalid_argument("Invalid reference-surface support");
        if(order!=1 && order!=2)throw std::invalid_argument("Reference surface supports polynomial order 1 or 2");
        const Real support=spacing*support_ratio;
        using Cell=std::array<long long,3>;
        struct Hash {size_t operator()(const Cell &a) const
            {size_t h=0;for(auto x:a) h^=std::hash<long long>{}(x)+0x9e3779b9+(h<<6)+(h>>2);return h;}};
        std::unordered_map<Cell,StdVec<size_t>,Hash> cells;
        auto cell=[&](const Vecd &x){return Cell{(long long)std::floor(x[0]/support),
            (long long)std::floor(x[1]/support),(long long)std::floor(x[2]/support)};};
        for(size_t j=0;j<port.count;++j) cells[cell(port.reference[j])].push_back(j);
        rows.reserve(points.size());
        for(const auto &target:points)
        {
            Row row;StdVec<Eigen::Vector4d> basis;StdVec<Real> kernel;
            Eigen::Matrix4d moment=Eigen::Matrix4d::Zero();auto center=cell(target);
            for(int a=-1;a<=1;++a) for(int b=-1;b<=1;++b) for(int c=-1;c<=1;++c)
            {
                auto found=cells.find(Cell{center[0]+a,center[1]+b,center[2]+c});
                if(found==cells.end()) continue;
                for(size_t j:found->second)
                {
                    Vecd q=(port.reference[j]-target)/spacing;
                    if(q.squaredNorm()>=support_ratio*support_ratio) continue;
                    Eigen::Vector4d phi;phi<<1,q[0],q[1],q[2];
                    Real w=std::exp(-q.squaredNorm())*port.volume[j]/std::pow(spacing,3);
                    moment.noalias()+=w*phi*phi.transpose();
                    row.index.push_back(j);basis.push_back(phi);kernel.push_back(w);
                }
            }
            if(order==2)
            {
                // Degree two requires at least three resolved normal levels on
                // a planar volume support. Reject rank loss rather than adding
                // a diagonal regularizer that silently changes reproduction.
                using V10=Eigen::Matrix<Real,10,1>;using M10=Eigen::Matrix<Real,10,10>;
                Real total=0;Matd covariance=Matd::Zero();
                for(size_t a=0;a<row.index.size();++a)
                {Vecd q=(port.reference[row.index[a]]-target)/spacing;
                 total+=kernel[a];covariance+=kernel[a]*q*q.transpose();}
                if(row.index.size()<10 || total<=0)throw std::runtime_error("Insufficient quadratic surface support");
                covariance/=total;Eigen::SelfAdjointEigenSolver<Matd> ceig(covariance);
                if(ceig.info()!=Eigen::Success || ceig.eigenvalues().minCoeff()<=0)
                    throw std::runtime_error("Invalid quadratic surface covariance");
                Matd L=ceig.eigenvectors()*ceig.eigenvalues().cwiseSqrt().cwiseInverse().asDiagonal()*ceig.eigenvectors().transpose();
                auto basis2=[&](size_t j)->V10 {Vecd z=L*((port.reference[j]-target)/spacing);
                    V10 phi;phi<<1,z[0],z[1],z[2],z[0]*z[0],z[1]*z[1],z[2]*z[2],z[0]*z[1],z[0]*z[2],z[1]*z[2];return phi;};
                M10 m=M10::Zero();
                for(size_t a=0;a<row.index.size();++a){V10 phi=basis2(row.index[a]);m+=kernel[a]*phi*phi.transpose();}
                Eigen::SelfAdjointEigenSolver<M10> eig2(m);
                if(eig2.info()!=Eigen::Success || eig2.eigenvalues().minCoeff()<=1e-12*eig2.eigenvalues().maxCoeff())
                    throw std::runtime_error("Rank-deficient quadratic volume-to-surface mapping; increase resolved normal support");
                largest_condition_number=std::max(largest_condition_number,eig2.eigenvalues().maxCoeff()/eig2.eigenvalues().minCoeff());
                V10 e0=V10::Zero();e0[0]=1;V10 coeff=m.ldlt().solve(e0),reproduced_basis=V10::Zero();
                if(!coeff.allFinite())throw std::runtime_error("Non-finite quadratic reference surface coefficients");
                Real sum=0;Vecd reproduced=Vecd::Zero();
                for(size_t a=0;a<row.index.size();++a)
                {V10 phi=basis2(row.index[a]);Real w=kernel[a]*phi.dot(coeff);row.weight.push_back(w);
                 sum+=w;reproduced+=w*port.reference[row.index[a]];reproduced_basis+=w*phi;}
                affine_error=std::max({affine_error,std::abs(sum-1),(reproduced-target).norm()/spacing});
                polynomial_error=std::max(polynomial_error,(reproduced_basis-e0).cwiseAbs().maxCoeff());
                if(!std::isfinite(polynomial_error) || polynomial_error>1e-8 || !std::isfinite(affine_error) || affine_error>1e-8)
                    throw std::runtime_error("Quadratic reference surface reproduction failed");
                rows.push_back(std::move(row));continue;
            }
            // Preserve the legacy affine arithmetic and default behavior.
            Eigen::SelfAdjointEigenSolver<Eigen::Matrix4d> eig(moment);
            if(row.index.size()<8 || eig.info()!=Eigen::Success ||
                eig.eigenvalues().minCoeff()<=1e-12*eig.eigenvalues().maxCoeff())
                throw std::runtime_error("Rank-deficient volume-to-surface mapping; resolve physical thickness/support");
            largest_condition_number=std::max(largest_condition_number,
                eig.eigenvalues().maxCoeff()/eig.eigenvalues().minCoeff());
            Eigen::Vector4d e0(1,0,0,0),coeff=moment.ldlt().solve(e0);
            Real sum=0;Vecd reproduced=Vecd::Zero();
            for(size_t a=0;a<row.index.size();++a)
            {Real w=kernel[a]*basis[a].dot(coeff);row.weight.push_back(w);sum+=w;reproduced+=w*port.reference[row.index[a]];}
            affine_error=std::max({affine_error,std::abs(sum-1),(reproduced-target).norm()/spacing});
            if(!std::isfinite(affine_error) || affine_error>1e-8)
                throw std::runtime_error("Reference surface does not reproduce affine motion");
            rows.push_back(std::move(row));
        }
    }
    Vecd position(size_t i) const {ensureCurrent();Vecd x=Vecd::Zero();const auto &r=rows[i];
        for(size_t a=0;a<r.index.size();++a)x+=r.weight[a]*port.position[r.index[a]];return x;}
    Vecd velocity(size_t i) const {ensureCurrent();Vecd v=Vecd::Zero();const auto &r=rows[i];
        for(size_t a=0;a<r.index.size();++a)v+=r.weight[a]*port.velocity[r.index[a]];return v;}
    void kinematics(size_t i,Vecd&x,Vecd&v) const
    {
        ensureCurrent();x.setZero();v.setZero();const auto&r=rows[i];
        for(size_t a=0;a<r.index.size();++a){x+=r.weight[a]*port.position[r.index[a]];v+=r.weight[a]*port.velocity[r.index[a]];}
    }
    void scatter(size_t i,const Vecd &f) {ensureCurrent();const auto &r=rows[i];
        for(size_t a=0;a<r.index.size();++a)port.addForce(r.index[a],r.weight[a]*f);}

    /** Absolute row-sum bound of sqrt(K) B M^-1 B^T sqrt(K). */
    static Real stableBound(const ReferenceSolidSurface &a,const ReferenceSolidSurface *b,
                            const StdVec<Real> &stiffness,Real safety=.2,
                            const StdVec<size_t>* included_rows=nullptr)
    {
        a.ensureCurrent();if(b)b->ensureCurrent();
        if(stiffness.size()!=a.rows.size() || (b && b->rows.size()!=a.rows.size()))
            throw std::invalid_argument("Surface constraint dimensions disagree");
        const size_t row_count=included_rows?included_rows->size():stiffness.size();
        if(row_count==0)return std::numeric_limits<Real>::infinity();
        StdVec<Real> sa(a.port.count,0),sb(b ? b->port.count:0,0);
        for(size_t q=0;q<row_count;++q)
        {
            const size_t i=included_rows?included_rows->at(q):q;
            if(i>=stiffness.size())throw std::invalid_argument("Invalid active surface row");
            if(!std::isfinite(stiffness[i]) || stiffness[i]<=0)throw std::invalid_argument("Positive constraint stiffness required");
            Real k=std::sqrt(stiffness[i]);const auto &r=a.rows[i];
            for(size_t j=0;j<r.index.size();++j)sa[r.index[j]]+=k*std::abs(r.weight[j]);
            if(b){const auto &s=b->rows[i];for(size_t j=0;j<s.index.size();++j)sb[s.index[j]]+=k*std::abs(s.weight[j]);}
        }
        Real largest=0;
        for(size_t q=0;q<row_count;++q)
        {
            const size_t i=included_rows?(*included_rows)[q]:q;
            Real rate=0,k=std::sqrt(stiffness[i]);const auto &r=a.rows[i];
            for(size_t j=0;j<r.index.size();++j)rate+=k*std::abs(r.weight[j])*sa[r.index[j]]/a.port.mass[r.index[j]];
            if(b){const auto &s=b->rows[i];for(size_t j=0;j<s.index.size();++j)rate+=k*std::abs(s.weight[j])*sb[s.index[j]]/b->port.mass[s.index[j]];}
            largest=std::max(largest,rate);
        }
        return safety/std::sqrt(largest);
    }
};

namespace surface_execution_detail
{
inline ParticleTraceMap packTrace(const ReferenceSolidSurface &surface)
{
    ParticleTraceMap result{surface.port.count, {0}, {}};
    for (const auto &row : surface.rows)
    {
        for (size_t k = 0; k < row.index.size(); ++k)
            result.entries.push_back({row.index[k], row.weight[k]});
        result.offsets.push_back(result.entries.size());
    }
    return result;
}

inline void packState(const SolidForcePort &port, std::vector<TraceForceNode> &nodes)
{
    for (size_t i = 0; i < nodes.size(); ++i)
    {
        auto &node = nodes[i];
        node.force_divisor = port.storage == SolidForcePort::Storage::ForcePrior ? 1.0 : port.mass[i];
        for (int d = 0; d < 3; ++d)
        {
            node.position[d] = port.position[i][d];
            node.velocity[d] = port.velocity[i][d];
            node.external[d] = port.external[i][d];
        }
    }
}

inline Vecd vector(const double *values)
{
    return Vecd(values[0], values[1], values[2]);
}

} // namespace surface_execution_detail

/** One reference interface energy:U=sum(k_i|W_a*x_a-W_b*x_b|^2/2).
 * Both bodies receive the transpose pullback of the SAME physical force.
 * No fluid pressure, velocity overwrite or duplicate repulsion is added. */
class BondedSolidInterface
{
    std::unique_ptr<BondedInterfaceBatch> backend_;
    std::vector<BondedInterfaceNode> nodes_a_, nodes_b_;
    std::vector<BondedInterfacePoint> points_;

    void finishAudit(const Vecd &sum, const Vecd &torque, Real power, Real force_scale, Real power_scale)
    {
        force_residual = sum.norm() / std::max(force_scale, Real(1e-30));
        torque_residual = torque.norm() / std::max(length_scale * force_scale, Real(1e-30));
        power_residual = std::abs(power) / std::max(power_scale, Real(1e-30));
        worst_force = std::max(worst_force, force_residual);
        worst_torque = std::max(worst_torque, torque_residual);
        worst_power = std::max(worst_power, power_residual);
    }

    void applyBackend()
    {
        using namespace surface_execution_detail;
        packState(a.port, nodes_a_);
        packState(b.port, nodes_b_);
        backend_->apply(nodes_a_, nodes_b_, points_);
        if (points_.size() != stiffness.size())
            throw std::logic_error("Bonded interface backend returned the wrong point count");
        for (size_t i = 0; i < nodes_a_.size(); ++i)
            a.port.external[i] = vector(nodes_a_[i].external);
        for (size_t i = 0; i < nodes_b_.size(); ++i)
            b.port.external[i] = vector(nodes_b_[i].external);

        energy = max_gap = 0;
        Vecd sum = Vecd::Zero(), torque = Vecd::Zero();
        Real power = 0, force_scale = 0, power_scale = 0;
        for (size_t i = 0; i < points_.size(); ++i)
        {
            const auto &point = points_[i];
            x_a[i] = vector(point.position_a);
            x_b[i] = vector(point.position_b);
            force_a[i] = vector(point.force_a);
            energy += point.energy;
            max_gap = std::max(max_gap, Real(point.gap));
            sum += vector(point.sum);
            torque += vector(point.torque);
            power += point.power;
            force_scale += point.force_scale;
            power_scale += point.power_scale;
        }
        finishAudit(sum, torque, power, force_scale, power_scale);
    }

  public:
    ReferenceSolidSurface &a,&b;
    StdVec<Real> stiffness;
    StdVec<Vecd> x_a,x_b,force_a;
    Real stable_dt,energy=0,max_gap=0,force_residual=0,torque_residual=0,power_residual=0;
    Real worst_force=0,worst_torque=0,worst_power=0,length_scale;size_t a_mapping_version=~size_t(0),b_mapping_version=~size_t(0);
    BondedSolidInterface(ReferenceSolidSurface &first,ReferenceSolidSurface &second,
                         const StdVec<Real> &areas,Real stiffness_per_area,Real length,
                         std::unique_ptr<BondedInterfaceBatch> backend = nullptr)
        : backend_(std::move(backend)),a(first),b(second),length_scale(length)
    {
        if(a.rows.size()!=b.rows.size() || a.points.size()!=areas.size() || length<=0)
            throw std::invalid_argument("Invalid bonded-solid interface dimensions");
        for(size_t i=0;i<areas.size();++i)
        {if((a.points[i]-b.points[i]).norm()>1e-9*length || areas[i]<=0)
            throw std::invalid_argument("Bonded-solid reference surfaces do not coincide");
         stiffness.push_back(areas[i]*stiffness_per_area);}
        x_a.resize(areas.size());x_b.resize(areas.size());force_a.resize(areas.size());
        if (backend_)
        {
            if (a.port.external == b.port.external)
                throw std::invalid_argument("Batched bonded interface requires distinct external force arrays");
            nodes_a_.resize(a.port.count);
            nodes_b_.resize(b.port.count);
        }
        refreshTimeStepBound();
    }
    void refreshTimeStepBound()
    {
        a.ensureCurrent();
        b.ensureCurrent();
        stable_dt = ReferenceSolidSurface::stableBound(a, &b, stiffness);
        if (backend_)
            backend_->configure(surface_execution_detail::packTrace(a), surface_execution_detail::packTrace(b), std::vector<double>(stiffness.begin(), stiffness.end()));
        a_mapping_version = a.mappingVersion();
        b_mapping_version = b.mappingVersion();
    }
    void apply()
    {
        a.ensureCurrent();b.ensureCurrent();if(a_mapping_version!=a.mappingVersion()||b_mapping_version!=b.mappingVersion())refreshTimeStepBound();
        if (backend_)
        {
            applyBackend();
            return;
        }
        energy=max_gap=0;Vecd sum=Vecd::Zero(),torque=Vecd::Zero();Real power=0,fs=0,ps=0;
        for(size_t i=0;i<stiffness.size();++i)
        {
            if(a.rows[i].index.empty()||b.rows[i].index.empty())
            {x_b[i]=b.rows[i].index.empty()?Vecd::Zero():b.position(i);x_a[i]=x_b[i];force_a[i].setZero();continue;}
            Vecd xa,xb,va,vb;a.kinematics(i,xa,va);b.kinematics(i,xb,vb);
            Vecd gap=xa-xb,f=-stiffness[i]*gap;
            x_a[i]=xa;x_b[i]=xb;force_a[i]=f;
            energy+=.5*stiffness[i]*gap.squaredNorm();max_gap=std::max(max_gap,gap.norm());
            for(int side=0;side<2;++side)
            {
                const auto &map=side==0?a:b;const auto &row=map.rows[i];Real sign=side==0?1:-1;
                for(size_t j=0;j<row.index.size();++j)
                {
                    Vecd fj=sign*row.weight[j]*f;size_t id=row.index[j];
                    map.port.addForce(id,fj);
                    sum+=fj;torque+=map.port.position[id].cross(fj);
                    Real p=fj.dot(map.port.velocity[id]);power+=p;fs+=fj.norm();ps+=std::abs(p);
                }
            }
            Real energy_rate=-f.dot(va-vb);power+=energy_rate;ps+=std::abs(energy_rate);
        }
        finishAudit(sum, torque, power, fs, ps);
    }
};

/** Analytic parallel planes, with optional rate-independent Coulomb friction.
 * Plane motion is normal to z; the prescribed tangential plane velocity is zero.
 * The same surface interpolation and transpose force map serve all components. */
class SolidSurfacePlaneContact
{
    std::unique_ptr<SurfacePlaneContactBatch> backend_;
    std::vector<TraceForceNode> nodes_;
    std::vector<SurfacePlaneContactPoint> points_;
    SurfaceFrictionTransfer friction_transfer_;

    void applyBackend(Real bottom, Real top, bool commit_history)
    {
        using namespace surface_execution_detail;
        packState(surface.port, nodes_);
        if (friction>0)
        {
            auto &committed=friction_transfer_.committed;
            committed.resize(friction_history.size());
            for (size_t i=0;i<friction_history.size();++i)
            {
                const auto &h=friction_history[i];
                committed[i]={{h.anchor[0],h.anchor[1]},h.stored,h.side};
            }
        }
        backend_->apply(bottom, top, nodes_, points_, friction>0?&friction_transfer_:nullptr);
        if (points_.size() != stiffness.size())
            throw std::logic_error("Plane contact backend returned the wrong point count");
        for (size_t i = 0; i < nodes_.size(); ++i)
            surface.port.external[i] = vector(nodes_[i].external);

        Vecd net_force=Vecd::Zero(), torque=Vecd::Zero();
        Real power=0, force_scale=0, power_scale=0;
        for (size_t i = 0; i < points_.size(); ++i)
        {
            const auto &point = points_[i];
            position[i] = vector(point.position);
            force[i] = vector(point.force);
            upper += point.upper;
            lower += point.lower;
            energy += point.energy;
            max_penetration = std::max(max_penetration, Real(point.penetration));
            net_force += vector(point.sum);
            torque += vector(point.torque);
            power += point.power;
            force_scale += point.force_scale;
            power_scale += point.power_scale;
            if (friction>0)
            {
                auto &h=friction_history[i];
                const auto &trial=friction_transfer_.trials[i];
                const auto &update=trial.update;
                const Real normal=point.upper+point.lower;
                if (normal>0)
                {
                    const Real tangential=std::sqrt(point.force[0]*point.force[0]+point.force[1]*point.force[1]);
                    ++active_friction_points;
                    sliding_points+=update.sliding;
                    tangential_force_sum+=tangential;
                    max_cone_ratio=std::max(max_cone_ratio,tangential/(friction*normal));
                }
                tangent_energy+=update.stored;
                if (commit_history&&(point.penetration>0||h.side))
                {
                    FrictionHistory::Result result;
                    result.anchor=Vecd(trial.anchor[0],trial.anchor[1],0);
                    result.elastic_slip=Vecd(trial.elastic_slip[0],trial.elastic_slip[1],0);
                    result.stored=update.stored;
                    result.slip_dissipation=update.slip_dissipation;
                    result.release_energy=update.release_energy;
                    result.side=update.side;
                    result.base_revision=h.revision;
                    h.commit(result);
                    slip_dissipation+=update.slip_dissipation;
                    release_energy+=update.release_energy;
                }
            }
        }
        energy+=tangent_energy;
        if(friction>0&&commit_history)++history_commits;
        worst_force=std::max(worst_force,net_force.norm()/std::max(force_scale,Real(1e-30)));
        worst_torque=std::max(worst_torque,torque.norm()/std::max(length_scale*force_scale,Real(1e-30)));
        worst_power=std::max(worst_power,std::abs(power)/std::max(power_scale,Real(1e-30)));
    }

  public:
    ReferenceSolidSurface &surface;
    StdVec<Real> areas,stiffness;
    StdVec<Vecd> position,force;
    using FrictionHistory=CoulombTangentialHistory<Vecd>;
    StdVec<FrictionHistory> friction_history;
    Real friction,tangent_factor;
    Real stable_dt,upper=0,lower=0,energy=0,max_penetration=0;
    Real tangent_energy=0,slip_dissipation=0,release_energy=0,max_cone_ratio=0;
    Real tangential_force_sum=0;
    size_t active_friction_points=0,sliding_points=0,history_commits=0;size_t surface_mapping_version=~size_t(0);
    Real worst_force=0,worst_torque=0,worst_power=0,length_scale;
    SolidSurfacePlaneContact(ReferenceSolidSurface &s,const StdVec<Real> &area,
                             Real stiffness_per_area,Real length,Real mu=0,Real kt_over_kn=1,
                             std::unique_ptr<SurfacePlaneContactBatch> backend=nullptr)
        : backend_(std::move(backend)),surface(s),areas(area),friction(mu),tangent_factor(kt_over_kn),length_scale(length)
    {
        if(areas.size()!=s.rows.size()||!std::isfinite(mu)||mu<0||
           !std::isfinite(kt_over_kn)||kt_over_kn<=0)
            throw std::invalid_argument("Invalid contact areas or Coulomb parameters");
        for(Real a:areas){if(a<=0)throw std::invalid_argument("Contact areas must be positive");stiffness.push_back(a*stiffness_per_area);}
        position.resize(area.size());force.resize(area.size());
        if(friction>0)friction_history.resize(area.size());
        refreshTimeStepBound();
    }
    bool onDevice() const noexcept { return bool(backend_); }

    void refreshTimeStepBound()
    {
        surface.ensureCurrent();
        // Stick: Kt P_t + Kn n*n^T. Sliding also couples Ft to Fn with
        // derivative bounded by mu*Kn. The sum of these operator norms is
        // a conservative bound; mu=0 preserves the original bound exactly.
        if(friction==0)
            stable_dt=ReferenceSolidSurface::stableBound(surface,nullptr,stiffness);
        else
        {
            StdVec<Real> bound=stiffness;
            for(auto &k:bound)k*=1+tangent_factor+friction;
            stable_dt=ReferenceSolidSurface::stableBound(surface,nullptr,bound);
        }
        if (onDevice())
        {
            nodes_.resize(surface.port.count);
            backend_->configure(surface_execution_detail::packTrace(surface),
                                std::vector<double>(stiffness.begin(),stiffness.end()),friction,tangent_factor);
        }
        surface_mapping_version=surface.mappingVersion();
    }
    void apply(Real bottom,Real top,bool commit_history=false)
    {
        surface.ensureCurrent();if(surface_mapping_version!=surface.mappingVersion())refreshTimeStepBound();
        upper=lower=energy=max_penetration=tangent_energy=max_cone_ratio=tangential_force_sum=0;
        active_friction_points=sliding_points=0;
        if (onDevice())
        {
            applyBackend(bottom, top, commit_history);
            return;
        }
        Vecd ferr=Vecd::Zero(),terr=Vecd::Zero();Real perr=0,fs=0,ps=0;
        for(size_t i=0;i<stiffness.size();++i)
        {
            if(surface.rows[i].index.empty())
            {
                if(friction>0&&commit_history&&friction_history[i].side)
                {auto r=friction_history[i].trial(Vecd::Zero(),0,0,stiffness[i]*tangent_factor,friction);
                 friction_history[i].commit(r);release_energy+=r.release_energy;}
                position[i].setZero();force[i].setZero();continue;
            }
            Vecd x=surface.position(i);
            Real pb=std::max(bottom-x[2],Real(0)),pt=std::max(x[2]-top,Real(0));
            position[i]=x;Vecd ft=Vecd::Zero();
            if(friction>0)
            {
                int side=pb>0?-1:pt>0?1:0;
                if(side||friction_history[i].side)
                {
                    Vecd tangent=x;tangent[2]=0;
                    auto r=friction_history[i].trial(tangent,side,stiffness[i]*(pb+pt),stiffness[i]*tangent_factor,friction);
                    ft=r.force;tangent_energy+=r.stored;sliding_points+=r.sliding;active_friction_points+=side!=0;
                    tangential_force_sum+=ft.norm();
                    if(pb+pt>0)max_cone_ratio=std::max(max_cone_ratio,ft.norm()/(friction*stiffness[i]*(pb+pt)));
                    if(commit_history){friction_history[i].commit(r);slip_dissipation+=r.slip_dissipation;release_energy+=r.release_energy;}
                }
            }
            // The default zero-friction arithmetic follows the original path.
            if(pb==0&&pt==0){force[i].setZero();continue;}
            Vecd v=surface.velocity(i),f(0,0,stiffness[i]*(pb-pt));
            if(friction>0)f+=ft;
            force[i]=f;upper+=stiffness[i]*pt;lower+=stiffness[i]*pb;
            energy+=.5*stiffness[i]*(pb*pb+pt*pt);max_penetration=std::max({max_penetration,pb,pt});
            surface.scatter(i,f);ferr-=f;terr-=x.cross(f);perr-=f.dot(v);ps+=std::abs(f.dot(v));
            const auto &r=surface.rows[i];
            for(size_t j=0;j<r.index.size();++j)
            {size_t id=r.index[j];Vecd fj=r.weight[j]*f;ferr+=fj;terr+=surface.port.position[id].cross(fj);
             Real p=fj.dot(surface.port.velocity[id]);perr+=p;fs+=fj.norm();ps+=std::abs(p);}
        }
        energy+=tangent_energy;if(friction>0&&commit_history)++history_commits;
        worst_force=std::max(worst_force,ferr.norm()/std::max(fs,Real(1e-30)));
        worst_torque=std::max(worst_torque,terr.norm()/std::max(length_scale*fs,Real(1e-30)));
        worst_power=std::max(worst_power,std::abs(perr)/std::max(ps,Real(1e-30)));
    }
};

} // namespace SPH::peridynamics
#endif
