#include "pd_isoparametric_crack_surface.h"

#include <algorithm>
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

void CrackSurfaceContactNetwork::updateDeviceTraceMap()
{
    // The backend can be shared by more than one contact network.
    bool current = trace_generation_ != 0 && trace_generation_ == projection_batch_->contactTraceGeneration() &&
                   trace_versions_.size() == sides.size();
    for (std::size_t i = 0; current && i < sides.size(); ++i)
        current = trace_versions_[i].first == sides[i].surface &&
                  trace_versions_[i].second == sides[i].surface->mappingVersion();
    if (current)
        return;

    trace_generation_ = 0;
    ParticleTraceMap map{sides.front().surface->port.count, {0}, {}};
    trace_versions_.clear();
    trace_side_offsets_.clear();
    for (const auto &side : sides)
    {
        const auto &surface = *side.surface;
        trace_side_offsets_.push_back(map.offsets.size() - 1);
        trace_versions_.emplace_back(&surface, surface.mappingVersion());
        // The last rows reconstruct normals and never receive contact loads.
        for (size_t i = 0; i < 4 + surface.natural.size(); ++i)
        {
            const auto &row = surface.rows[i];
            for (size_t j = 0; j < row.index.size(); ++j)
                map.entries.push_back({row.index[j], row.weight[j]});
            map.offsets.push_back(map.entries.size());
        }
    }
    trace_generation_ = projection_batch_->configureContactTrace(map);
}

void CrackSurfaceContactNetwork::applyDeviceTracePullback(SolidForcePort &port)
{
    updateDeviceTraceMap();
    contact_addresses_.clear();
    for (const auto &point : selected_contacts_)
    {
        const auto slave = trace_side_offsets_[point.slave] + 4 + point.quadrature;
        const auto master = trace_side_offsets_[point.master];
        contact_addresses_.push_back({{slave, master, master + 1, master + 2, master + 3}});
    }
    trace_particle_input_.resize(port.count);
    for (size_t i = 0; i < port.count; ++i)
    {
        trace_particle_input_[i].mass = port.mass[i];
        std::copy_n(port.external[i].data(), 3, trace_particle_input_[i].external);
    }
    const bool acceleration_port = port.storage == SolidForcePort::Storage::PDAccelerationPrior;
    projection_batch_->assembleContactTrace(trace_generation_, contact_addresses_, trace_particle_input_,
                                            acceleration_port, trace_particle_output_);
    if (trace_particle_output_.size() != port.count)
        throw std::logic_error("Crack backend returned an incomplete particle force batch");
    for (size_t i = 0; i < port.count; ++i)
    {
        const auto &particle = trace_particle_output_[i];
        port.external[i] = Vecd(particle.external[0], particle.external[1], particle.external[2]);
        rates[i] = Vecd(particle.rate[0], particle.rate[1], particle.rate[2]);
    }
    device_trace_pullback = true;
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
    using TraceRows = std::array<const ReferenceSolidSurface::Row *, 5>;
    using TraceWeights = std::array<Real, 5>;
    auto traceRows = [&](const ContactPoint &point) -> TraceRows
    {
        const auto &slave = *sides[point.slave].surface;
        const auto &master = *sides[point.master].surface;
        return {&slave.rows[4 + point.quadrature], &master.rows[0],
                &master.rows[1], &master.rows[2], &master.rows[3]};
    };
    auto massWeights = [&](const TraceRows &rows) -> TraceWeights
    {
        TraceWeights result{};
        const Real *mass = sides.front().surface->port.mass;
        for (size_t i = 0; i < rows.size(); ++i)
            for (size_t j = 0; j < rows[i]->index.size(); ++j)
                result[i] += std::abs(rows[i]->weight[j]) / std::sqrt(mass[rows[i]->index[j]]);
        return result;
    };
    auto hostContact = [](const ContactPoint &point, const TraceWeights &weighted)
    {
        const auto &a = point.geometry;
        Vecd residual = Vecd::Zero();
        if (point.gap < 0)
            residual = point.gap * point.normal;
        const Real off_diagonal = a.u.dot(a.v) - residual.dot(a.uv);
        Eigen::Matrix2d metric;
        metric << a.u.squaredNorm(), off_diagonal, off_diagonal, a.v.squaredNorm();
        Eigen::SelfAdjointEigenSolver<Eigen::Matrix2d> eigenvalues(metric);
        if (eigenvalues.info() != Eigen::Success ||
            eigenvalues.eigenvalues().minCoeff() <= 1e-13 * eigenvalues.eigenvalues().maxCoeff())
            throw std::runtime_error("Closest crack contact tangent is not a local minimum");

        Eigen::Matrix<Real, 15, 2> coupling;
        const std::array<Real, 5> c{1, -a.N[0], -a.N[1], -a.N[2], -a.N[3]};
        for (size_t i = 0; i < c.size(); ++i)
        {
            Vecd u = -c[i] * a.u, v = -c[i] * a.v;
            if (i != 0)
            {
                u -= a.Nu[i - 1] * residual;
                v -= a.Nv[i - 1] * residual;
            }
            coupling.template block<3, 1>(3 * i, 0) = u;
            coupling.template block<3, 1>(3 * i, 1) = v;
        }
        Eigen::Matrix<Real, 15, 15> tangent = -coupling * metric.inverse() * coupling.transpose();
        for (size_t i = 0; i < c.size(); ++i)
            for (size_t j = 0; j < c.size(); ++j)
                tangent.template block<3, 3>(3 * i, 3 * j) += c[i] * c[j] * Matd::Identity();
        tangent *= point.stiffness;

        CrackLocalContactResult result{};
        result.valid = true;
        for (size_t i = 0; i < c.size(); ++i)
            for (size_t d = 0; d < 3; ++d)
            {
                Real row = 0;
                for (size_t j = 0; j < 15; ++j)
                    row += std::abs(tangent(3 * i + d, j)) * weighted[j / 3];
                result.row_bound[i][d] = row;
            }
        if (point.gap < 0)
        {
            const Vecd force = -point.stiffness * point.gap * point.normal;
            std::copy_n(force.data(), 3, result.force_on_slave);
        }
        return result;
    };
    auto deviceInput = [](const ContactPoint &point, const TraceWeights &weighted)
    {
        CrackLocalContactInput input;
        const auto &a = point.geometry;
        std::copy_n(a.u.data(), 3, input.tangent_u);
        std::copy_n(a.v.data(), 3, input.tangent_v);
        std::copy_n(a.uv.data(), 3, input.mixed_tangent);
        std::copy_n(point.normal.data(), 3, input.normal);
        std::copy_n(a.N.data(), 4, input.shape);
        std::copy_n(a.Nu.data(), 4, input.derivative_u);
        std::copy_n(a.Nv.data(), 4, input.derivative_v);
        std::copy_n(weighted.data(), 5, input.mass_weight);
        input.gap = point.gap;
        input.stiffness = point.stiffness;
        return input;
    };

    Vecd balance = Vecd::Zero(), moment = Vecd::Zero();
    Real magnitude = 0, length = spacing, power = 0, power_scale = 0, energy_rate = 0;
    auto accumulate = [&](const ContactPoint &point, const TraceRows &rows,
                          const CrackLocalContactResult &local, bool scatter_on_host)
    {
        auto &slave = *sides[point.slave].surface;
        auto &master = *sides[point.master].surface;
        if (scatter_on_host)
            for (size_t i = 0; i < rows.size(); ++i)
                for (size_t d = 0; d < 3; ++d)
                {
                    const Real row_bound = local.row_bound[i][d];
                    for (size_t j = 0; j < rows[i]->index.size(); ++j)
                    {
                        const auto node = rows[i]->index[j];
                        rates[node][d] += std::abs(rows[i]->weight[j]) * row_bound / std::sqrt(slave.port.mass[node]);
                    }
                }

        if (point.gap >= 0)
            return;
        if (sides[point.slave].reference_patch != sides[point.master].reference_patch)
            ++neighbor_contacts;
        const Vecd force(local.force_on_slave[0], local.force_on_slave[1], local.force_on_slave[2]);
        if (scatter_on_host)
            slave.scatter(4 + point.quadrature, force);
        Vecd relative = point.velocity;
        balance += force;
        moment += point.position.cross(force);
        magnitude += force.norm();
        length = std::max(length, point.position.norm());
        power += force.dot(point.velocity);
        power_scale += std::abs(force.dot(point.velocity));

        const auto &master_state = masters[point.master];
        for (size_t j = 0; j < 4; ++j)
        {
            const Vecd reaction = -point.geometry.N[j] * force;
            if (scatter_on_host)
                master.scatter(j, reaction);
            balance += reaction;
            moment += master_state.x[j].cross(reaction);
            magnitude += reaction.norm();
            power += reaction.dot(master_state.v[j]);
            power_scale += std::abs(reaction.dot(master_state.v[j]));
            relative -= point.geometry.N[j] * master_state.v[j];
        }
        energy_rate += point.stiffness * point.gap * point.normal.dot(relative);
        energy += .5 * point.stiffness * point.gap * point.gap;
        max_penetration = std::max(max_penetration, -point.gap);
        normal_force_sum += force.norm();
        tangential_force_sum += (force - point.normal * force.dot(point.normal)).norm();
        ++active_contacts;
    };

    const bool device_contact = projection_batch_ && projection_batch_->assemblesLocalContact();
    selected_contacts_.clear();
    contact_inputs_.clear();
    device_contact_evaluations = host_contact_fallbacks = 0;
    device_trace_pullback = false;
    for (size_t si = 0; si < sides.size(); ++si)
    {
        auto &slave = *sides[si].surface;
        for (size_t q = 0; q < slave.natural.size(); ++q)
        {
            Vecd x, v;
            slave.kinematics(4 + q, x, v);
            const Vecd tangent_u = slave.position(slave.normalU(q)) - x;
            const Vecd tangent_v = slave.position(slave.normalV(q)) - x;
            const Vecd outward = sides[si].outward_sign * tangent_u.cross(tangent_v).normalized();
            if (!outward.allFinite())
                throw std::runtime_error("Slave crack normal is unresolved");
            const auto candidates = grid.find(cell(x));
            if (candidates == grid.end())
                continue;

            size_t chosen = sides.size();
            Real distance = std::numeric_limits<Real>::infinity(), gap = 0;
            BilinearCrackGeometry::Value selected;
            Vecd normal;
            for (size_t master : candidates->second)
            {
                if (master == si)
                    continue;
                auto p = closest(master, x);
                if (p.minCoeff() < -1e-9 || p.maxCoeff() > 1 + 1e-9)
                {
                    ++outside_projections;
                    continue;
                }
                const auto a = BilinearCrackGeometry::evaluate(masters[master].x, p);
                const Vecd n = sides[master].outward_sign * a.u.cross(a.v).normalized();
                if (!n.allFinite())
                    throw std::runtime_error("Master crack normal is unresolved");
                if (n.dot(outward) >= -1e-8)
                    continue;
                const Real g = n.dot(x - a.x);
                if (std::abs(g) < distance)
                {
                    distance = std::abs(g);
                    chosen = master;
                    gap = g;
                    selected = a;
                    normal = n;
                }
            }
            if (chosen == sides.size())
                continue;

            const ContactPoint point{si, q, chosen, x, v, normal, selected, gap, .5 * penalty * sides[si].areas->at(q)};
            const TraceRows rows = traceRows(point);
            const auto weighted = massWeights(rows);
            if (device_contact)
            {
                selected_contacts_.push_back(point);
                contact_inputs_.push_back(deviceInput(point, weighted));
            }
            else
                accumulate(point, rows, hostContact(point, weighted), true);
        }
    }

    if (device_contact)
    {
        projection_batch_->assembleLocalContact(contact_inputs_, contact_results_);
        if (contact_results_.size() != selected_contacts_.size())
            throw std::logic_error("Crack contact backend returned an incomplete batch");

        auto &port = sides.front().surface->port;
        const bool one_force_port = std::all_of(sides.begin(), sides.end(), [&](const Side &side)
        {
            return side.surface->port.external == port.external && side.surface->port.storage == port.storage &&
                   side.surface->port.count == port.count;
        });
        const bool valid_batch = std::all_of(contact_results_.begin(), contact_results_.end(), [](const auto &result)
        {
            return result.valid;
        });
        if (!selected_contacts_.empty() && one_force_port && valid_batch && projection_batch_->assemblesTracePullback())
            applyDeviceTracePullback(port);

        for (size_t i = 0; i < selected_contacts_.size(); ++i)
        {
            const auto &point = selected_contacts_[i];
            const auto rows = traceRows(point);
            if (contact_results_[i].valid)
            {
                accumulate(point, rows, contact_results_[i], !device_trace_pullback);
                ++device_contact_evaluations;
            }
            else
            {
                accumulate(point, rows, hostContact(point, massWeights(rows)), true);
                ++host_contact_fallbacks;
            }
        }
    }

    Real bound=0;for(const auto&r:rates)bound=std::max(bound,r.maxCoeff());stable_dt=bound>0?.2/std::sqrt(bound):std::numeric_limits<Real>::infinity();force_residual=balance.norm()/std::max(magnitude,Real(1e-30));torque_residual=moment.norm()/std::max(magnitude*length,Real(1e-30));power_residual=std::abs(power+energy_rate)/std::max(power_scale,Real(1e-30));
}
}
