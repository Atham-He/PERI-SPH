#pragma once
/** Live-bond reconstruction for resolved-thickness NOSB-PD shells.
 * Geometry cuts remove reciprocal material interactions. No particles, masses,
 * material constants, fracture energies or interface springs are deleted here.
 * The normal derivative stays inside each material thickness column; same-level
 * tangential derivatives are rebuilt from the visible live support. This class
 * supplies kinematics/topology only, NOT a fracture initiation/evolution law.
 * Its use therefore does not by itself validate strength or Gc-controlled growth.
 */
#include "adaptation.h"
#include "pd_consistent_reconstruction.h"
#include "pd_layered_correspondence.h"
#include <Eigen/Geometry>
#include <array>
#include <cstdint>
#include <utility>
#include <numeric>
#include <algorithm>

namespace SPH::solid_dynamics
{
class PDLayeredFractureReference : public PDLayeredQuadraticReference
{
  public:
    using V9=Eigen::Matrix<Real,9,1>;using M9=Eigen::Matrix<Real,9,9>;
    /** Reference crack geometry. Disk radius infinity gives an infinite plane.
     * Triangles/disks are barriers only; no penalty or cohesive law is hidden. */
    struct CrackFacet
    {
        enum class Kind { PlaneDisk, Triangle };
        Kind kind=Kind::PlaneDisk;
        Vecd origin=Vecd::Zero(),normal=Vecd::Zero();
        std::array<Vecd,3> vertex{Vecd::Zero(),Vecd::Zero(),Vecd::Zero()};
        Real radius=std::numeric_limits<Real>::infinity();
        uint64_t id=0;
    };
    struct TopologyState
    {
        StdVec<StdVec<unsigned char>> live;
        StdVec<CrackFacet> facets;
        size_t physical_broken_bonds=0;StdVec<Vecd>transfer_positions,transfer_velocity;StdVec<Matd>transfer_gradient;
    };
    StdVec<Vecd> X,transfer_positions,transfer_velocity;StdVec<Matd>transfer_gradient;
    StdVec<int> degree;
    StdVec<size_t> changed_particles;
    size_t topology_version=0,affine_rows=0,physical_broken_bonds=0;
    Real maximum_affine_error=0,maximum_quadratic_error=0;
    Real maximum_gradient_norm=0;

    explicit PDLayeredFractureReference(BaseInnerRelation&r,const peridynamics::LayeredParticleGeometry&g)
      :PDLayeredQuadraticReference(r,g),
       broken_fraction(particles_->registerStateVariableData<Real>("BrokenBondFraction")),
       reproduction(particles_->registerStateVariableData<int>("PDReproductionOrder")),
       bond_count(particles_->getVariableDataByName<int>("BondNumber")) {}

    void initialize(){build();}
    void build() override
    {
        if(initialized_){updateAfterFailure();return;}
        // Canonicalize all coupled neighbour arrays once, before any layered
        // stencil/reverse index is built. Parallel cell insertion order must
        // not change a cross-process physical checkpoint or summation order.
        for(auto&nb:inner_configuration_){std::vector<size_t>order(nb.current_size_);std::iota(order.begin(),order.end(),0);std::sort(order.begin(),order.end(),[&](size_t a,size_t b){return nb.j_[a]<nb.j_[b];});auto permute=[&](auto&values){if(values.size()<nb.current_size_)throw std::runtime_error("Incomplete reference neighbour arrays");auto old=values;for(size_t n=0;n<order.size();++n)values[n]=old[order[n]];};permute(nb.j_);permute(nb.W_ij_);permute(nb.dW_ij_);permute(nb.r_ij_);permute(nb.e_ij_);permute(nb.bondLive_);}
        X=geometry.reference;
        PDLayeredQuadraticReference::build();
        if(X.size()!=rows.size())throw std::logic_error("Layered fracture reference geometry mismatch");
        cached_live_.resize(rows.size());degree.assign(rows.size(),2);
        size_t initial_broken_bonds=0;
        for(size_t i=0;i<rows.size();++i)
        {
            cached_live_[i].assign(rows[i].size(),1);
            const auto&nb=inner_configuration_[i];
            for(size_t n=0;n<rows[i].size();++n)
            {
                auto&e=rows[i][n];
                if(nb.bondLive_[n]!=inner_configuration_[e.j].bondLive_[e.reverse])
                    throw std::runtime_error("Nonreciprocal initial fracture topology");
                if(i<e.j&&!nb.bondLive_[n])++initial_broken_bonds;
            }
            broken_fraction[i]=0;reproduction[i]=2;bond_count[i]=int(rows[i].size());
        }
        physical_broken_bonds=initial_broken_bonds;
        initialized_=true;
        if(initial_broken_bonds)updateAfterFailure();
    }
    const Neighborhood&neighbors(size_t i)const{return inner_configuration_.at(i);}
    bool live(size_t i,size_t n)const{return inner_configuration_.at(i).bondLive_.at(n);}
    const StdVec<CrackFacet>&crackFacets()const{return facets_;}
    bool visibleReference(const Vecd&a,const Vecd&b)const
    {for(const auto&f:facets_)if(intersects(f,a,b,false))return false;return true;}
    // A trace exactly on a crack intersection has more than one material
    // limit. Resolve incident faces from the interior of the owning patch;
    // the specified inward normal distinguishes the two limits on its plane.
    StdVec<Vecd> traceSector(const Vecd &point, const Vecd &patch_interior, const Vecd &inward) const
    {
        StdVec<Vecd> sides;
        for (const auto &facet : facets_)
        {
            const Real tol=tolerance(facet,point,patch_interior);
            if (std::abs(facet.normal.dot(point-facet.origin))>tol || !contains(facet,point,tol))
                continue;
            Real side=facet.normal.dot(patch_interior-facet.origin);
            if (std::abs(side)<=tol)
                side=facet.normal.dot(inward);
            if (side!=0)
                sides.push_back(side>0 ? facet.normal : Vecd(-facet.normal));
        }
        return sides;
    }
    bool visible(size_t anchor,size_t node)const
    {
        requireInitialized();
        for(size_t n=0;n<rows.at(anchor).size();++n)if(rows[anchor][n].j==node&&!live(anchor,n))return false;
        return visibleReference(X.at(anchor),X.at(node));
    }
    bool topologyPending()const{return geometry_changed_||pending_bonds_||reconstruction_pending_;}
    void setTransferState(const StdVec<Vecd>&q,const StdVec<Matd>&F,const StdVec<Vecd>&v){
        if(topologyPending()||q.size()!=X.size()||F.size()!=X.size()||v.size()!=X.size())throw std::invalid_argument("Invalid/overlapping conservative transfer state");
        for(size_t i=0;i<q.size();++i)if(!q[i].allFinite()||!v[i].allFinite()||!F[i].allFinite()||F[i].determinant()<=0)throw std::invalid_argument("Nonfinite conservative transfer state");
        transfer_positions=q;transfer_velocity=v;transfer_gradient=F;reconstruction_pending_=true;
    }
    Real transferVelocityTimeScale(){auto&mat=DynamicCast<HughesWingetSolid>(this,sph_body_->getMatterMaterial());return sph_body_->getSPHAdaptation().ReferenceSpacing()/mat.ReferenceSoundSpeed();}
    void registerTraceSupport(Real radius){if(!std::isfinite(radius)||radius<=0)throw std::invalid_argument("Invalid crack visibility support radius");maximum_trace_support_=std::max(maximum_trace_support_,radius);}

    /** An energy/growth controller may call this after its OWN physical event
     * check. Rebuilding is explicit, so a batch can be committed atomically. */
    bool severBondPair(size_t i,size_t n)
    {
        requireInitialized();
        if(i>=rows.size()||n>=rows[i].size())throw std::out_of_range("Invalid fracture bond");
        const auto&e=rows[i][n];
        if(geometry.column_of[i]==geometry.column_of[e.j])
            throw std::runtime_error("A crack cuts a thickness column: split-column kinematics are required; no normal bond was silently retained or removed");
        auto&nb=inner_configuration_[i];auto&other=inner_configuration_[e.j];
        if(nb.bondLive_[n]!=other.bondLive_[e.reverse])throw std::runtime_error("Nonreciprocal bond before fracture event");
        if(!nb.bondLive_[n])return false;
        nb.bondLive_[n]=false;other.bondLive_[e.reverse]=false;++physical_broken_bonds;pending_bonds_=true;
        return true;
    }
    size_t cutPlane(const Vecd&origin,const Vecd&normal,
                    Real radius=std::numeric_limits<Real>::infinity(),bool rebuild=true)
    {
        CrackFacet f;f.kind=CrackFacet::Kind::PlaneDisk;f.origin=origin;
        if(!origin.allFinite()||!normal.allFinite()||normal.norm()<=0||!std::isfinite(normal.norm())||!(radius>0))
            throw std::invalid_argument("Invalid precrack plane/disk");
        f.normal=normal.normalized();f.radius=radius;
        return insertFacet(f,rebuild);
    }
    size_t cutTriangle(const Vecd&a,const Vecd&b,const Vecd&c,bool rebuild=true)
    {
        if(!a.allFinite()||!b.allFinite()||!c.allFinite())throw std::invalid_argument("Nonfinite precrack triangle");
        Vecd n=(b-a).cross(c-a);Real scale=std::max((b-a).squaredNorm(),(c-a).squaredNorm());
        if(scale<=0||n.norm()<=128*std::numeric_limits<Real>::epsilon()*scale)
            throw std::invalid_argument("Degenerate precrack triangle");
        CrackFacet f;f.kind=CrackFacet::Kind::Triangle;f.origin=a;f.normal=n.normalized();f.vertex={a,b,c};
        return insertFacet(f,rebuild);
    }
    TopologyState captureTopology()const
    {
        requireInitialized();TopologyState t;t.transfer_positions=transfer_positions;t.transfer_velocity=transfer_velocity;t.transfer_gradient=transfer_gradient;t.facets=facets_;t.physical_broken_bonds=physical_broken_bonds;t.live.resize(rows.size());
        for(size_t i=0;i<rows.size();++i){t.live[i].resize(rows[i].size());for(size_t n=0;n<rows[i].size();++n)t.live[i][n]=live(i,n);}
        return t;
    }
    /** For transactional candidate evaluation/checkpoint restoration only.
     * No fracture history or thermodynamic acceptance is inferred by restore. */
    void restoreTopology(const TopologyState&t)
    {
        requireInitialized();
        validateTopologyState(t); // Complete input validation precedes any mutation.
        auto next_facets=t.facets;
        auto before=captureTopology();
        auto old_rows=rows;
        auto old_cached_live=cached_live_;
        auto old_degree=degree;
        auto old_changed=changed_particles;
        StdVec<Real> old_broken_fraction(broken_fraction,broken_fraction+rows.size());
        StdVec<int> old_reproduction(reproduction,reproduction+rows.size());
        StdVec<int> old_bond_count(bond_count,bond_count+rows.size());
        const auto old_version=topology_version,old_affine_rows=affine_rows;
        const auto old_minimum_rank=minimum_surface_rank;
        const Real old_affine_error=maximum_affine_error,old_quadratic_error=maximum_quadratic_error;
        const Real old_gradient_norm=maximum_gradient_norm,old_surface_error=max_surface_functional_error;
        const Real old_condition=max_moment_condition,old_reproduction_error=max_reproduction_error;
        const Real old_dt=elastic_dt_bound;
        const bool old_geometry_changed=geometry_changed_,old_pending_bonds=pending_bonds_;
        const bool old_full_geometry=full_geometry_rebuild_;const size_t old_committed_facets=committed_facets_;const bool old_reconstruction_pending=reconstruction_pending_;
        try
        {
            for(size_t i=0;i<rows.size();++i)for(size_t n=0;n<rows[i].size();++n)
                inner_configuration_[i].bondLive_[n]=bool(t.live[i][n]);
            transfer_positions=t.transfer_positions;transfer_velocity=t.transfer_velocity;transfer_gradient=t.transfer_gradient;facets_.swap(next_facets);physical_broken_bonds=t.physical_broken_bonds;geometry_changed_=true;full_geometry_rebuild_=true;
            updateAfterFailure();
        }
        catch(...)
        {
            // Rank loss/allocation failure during reconstruction must also roll
            // back derived stencils, diagnostics and the topology generation.
            for(size_t i=0;i<rows.size();++i)for(size_t n=0;n<rows[i].size();++n)
                inner_configuration_[i].bondLive_[n]=bool(before.live[i][n]);
            facets_.swap(before.facets);physical_broken_bonds=before.physical_broken_bonds;
            rows.swap(old_rows);cached_live_.swap(old_cached_live);degree.swap(old_degree);changed_particles.swap(old_changed);
            std::copy(old_broken_fraction.begin(),old_broken_fraction.end(),broken_fraction);
            std::copy(old_reproduction.begin(),old_reproduction.end(),reproduction);
            std::copy(old_bond_count.begin(),old_bond_count.end(),bond_count);
            topology_version=old_version;affine_rows=old_affine_rows;minimum_surface_rank=old_minimum_rank;
            maximum_affine_error=old_affine_error;maximum_quadratic_error=old_quadratic_error;
            maximum_gradient_norm=old_gradient_norm;max_surface_functional_error=old_surface_error;
            max_moment_condition=old_condition;max_reproduction_error=old_reproduction_error;elastic_dt_bound=old_dt;
            geometry_changed_=old_geometry_changed;pending_bonds_=old_pending_bonds;full_geometry_rebuild_=old_full_geometry;committed_facets_=old_committed_facets;transfer_positions=before.transfer_positions;transfer_velocity=before.transfer_velocity;transfer_gradient=before.transfer_gradient;reconstruction_pending_=old_reconstruction_pending;
            throw;
        }
    }

    void updateAfterFailure()
    {
        using namespace peridynamics::layered_detail;
        requireInitialized();changed_particles.clear();bool changed=geometry_changed_||reconstruction_pending_;
        for(size_t i=0;i<rows.size();++i)
        {
            const auto&nb=inner_configuration_[i];bool row_changed=false;
            if(nb.current_size_!=rows[i].size())throw std::runtime_error("Fracture reference neighborhood changed");
            for(size_t n=0;n<rows[i].size();++n)
            {
                const auto&e=rows[i][n];
                if(nb.j_[n]!=e.j||nb.bondLive_[n]!=inner_configuration_[e.j].bondLive_[e.reverse])
                    throw std::runtime_error("Layered fracture topology is not reciprocal");
                if(!nb.bondLive_[n]&&geometry.column_of[i]==geometry.column_of[e.j])
                    throw std::runtime_error("Thickness-column fracture requires split-column kinematics");
                row_changed|=cached_live_[i][n]!=static_cast<unsigned char>(nb.bondLive_[n]);
            }
            if(row_changed){changed=true;changed_particles.push_back(i);}
        }
        if(!changed){pending_bonds_=false;return;}
        // A new visibility barrier can affect interface traces even if its
        // nearby bonds had already failed. Rebuild trace owners conservatively.
        if(geometry_changed_){
            // Every trace segment is shorter than twice its registered support.
            // Expanded AABBs conservatively contain all affected endpoints.
            // Restoring topology still forces a complete visibility rebuild.
            StdVec<unsigned char>marked(rows.size(),0);for(size_t j:changed_particles)marked[j]=1;
            bool all=full_geometry_rebuild_||committed_facets_>facets_.size();
            Real margin=2*std::max(maximum_trace_support_,sph_body_->getSPHAdaptation().getKernel()->CutOffRadius());
            for(size_t k=committed_facets_;k<facets_.size()&&!all;++k){const auto&f=facets_[k];
                if(f.kind==CrackFacet::Kind::PlaneDisk&&std::isinf(f.radius)){all=true;break;}
                Vecd lo,hi;if(f.kind==CrackFacet::Kind::Triangle){lo=hi=f.vertex[0];for(const auto&v:f.vertex){lo=lo.cwiseMin(v);hi=hi.cwiseMax(v);}}else{lo=f.origin-Vecd::Constant(f.radius);hi=f.origin+Vecd::Constant(f.radius);}
                lo.array()-=margin;hi.array()+=margin;
                for(size_t i=0;i<X.size();++i)if((X[i].array()>=lo.array()).all()&&(X[i].array()<=hi.array()).all())marked[i]=1;
            }
            changed_particles.clear();for(size_t i=0;i<X.size();++i)if(all||marked[i])changed_particles.push_back(i);
            if(changed_particles.empty()){changed_particles.resize(rows.size());for(size_t i=0;i<rows.size();++i)changed_particles[i]=i;}
        }
        affine_rows=0;maximum_affine_error=maximum_quadratic_error=maximum_gradient_norm=0;
        minimum_surface_rank=9;max_surface_functional_error=0;max_moment_condition=0;
        Real dp=sph_body_->getSPHAdaptation().ReferenceSpacing();
        for(size_t i=0;i<rows.size();++i)
        {
            const auto&nb=inner_configuration_[i];const size_t level=geometry.layer_of[i],column=geometry.column_of[i];
            Vecd normal=geometry.director[column];Matd T=Matd::Identity()-normal*normal.transpose();
            Real sum=0;Matd covariance=Matd::Zero();size_t live_count=0;
            for(size_t n=0;n<rows[i].size();++n)
            {
                auto&e=rows[i][n];e.gradient.setZero();
                if(!nb.bondLive_[n])continue;++live_count;
                if(geometry.layer_of[e.j]!=level)continue;
                if(!visible(i,e.j))throw std::runtime_error("A live tangent bond crosses a recorded crack: topology/geometry were not committed consistently");
                Real w=nb.W_ij_[n]*volume[e.j];Vecd d=X[e.j]-X[i];sum+=w;covariance+=w*d*d.transpose();
            }
            if(sum<=0)throw std::runtime_error("Unresolved fractured shell support: no same-layer neighbors; refine or supply verified fragment kinematics, no material was erased");
            Matd L=whitening(covariance/sum,dp);M9 M=M9::Zero();
            for(size_t n=0;n<rows[i].size();++n)if(nb.bondLive_[n]&&geometry.layer_of[rows[i][n].j]==level)
            {V9 b=basis9(L*(X[rows[i][n].j]-X[i]));M+=nb.W_ij_[n]*volume[rows[i][n].j]/sum*b*b.transpose();}
            MomentInverse<9> inverse(M);max_moment_condition=std::max(max_moment_condition,inverse.condition);minimum_surface_rank=std::min(minimum_surface_rank,inverse.rank);
            Eigen::Matrix<Real,9,3>D=Eigen::Matrix<Real,9,3>::Zero();D.topRows<3>()=L*T;
            Real error=(M*inverse.inverse*D-D).norm()/std::max(D.norm(),std::numeric_limits<Real>::min());
            bool quadratic=inverse.rank>=5&&std::isfinite(error)&&error<=1e-9;
            if(quadratic)
            {
                for(size_t n=0;n<rows[i].size();++n)if(nb.bondLive_[n]&&geometry.layer_of[rows[i][n].j]==level)
                {V9 b=basis9(L*(X[rows[i][n].j]-X[i]));rows[i][n].gradient=nb.W_ij_[n]*volume[rows[i][n].j]/sum*D.transpose()*inverse.inverse*b;}
                max_surface_functional_error=std::max(max_surface_functional_error,error);
            }
            else
            {
                // On a curved one-sided patch the normal coordinate is a
                // nuisance affine mode. Include it in the moment, but request
                // only the TWO tangent derivatives. Flat patches have rank 2.
                Matd A=Matd::Zero();
                for(size_t n=0;n<rows[i].size();++n)if(nb.bondLive_[n]&&geometry.layer_of[rows[i][n].j]==level)
                {Vecd b=L*(X[rows[i][n].j]-X[i]);A+=nb.W_ij_[n]*volume[rows[i][n].j]/sum*b*b.transpose();}
                MomentInverse<3> affine(A);max_moment_condition=std::max(max_moment_condition,affine.condition);Matd target=L*T;
                Real affine_error=(A*affine.inverse*target-target).norm()/std::max(target.norm(),std::numeric_limits<Real>::min());
                if(affine.rank<2||!std::isfinite(affine_error)||affine_error>1e-9)
                    throw std::runtime_error("Unresolved live tangent derivative after fracture; no numerical support/bonds/particles were removed");
                for(size_t n=0;n<rows[i].size();++n)if(nb.bondLive_[n]&&geometry.layer_of[rows[i][n].j]==level)
                {Vecd b=L*(X[rows[i][n].j]-X[i]);rows[i][n].gradient=nb.W_ij_[n]*volume[rows[i][n].j]/sum*target.transpose()*affine.inverse*b;}
                ++affine_rows;max_surface_functional_error=std::max(max_surface_functional_error,affine_error);
            }
            for(size_t m=0;m<geometry.columns[column].size();++m)if(m!=level)
            {
                size_t j=geometry.columns[column][m];auto it=std::find_if(rows[i].begin(),rows[i].end(),[&](const Edge&e){return e.j==j;});
                if(it==rows[i].end())throw std::runtime_error("Normal column node left physical horizon");
                size_t n=size_t(it-rows[i].begin());if(!nb.bondLive_[n])throw std::runtime_error("Broken thickness-column support is not representable");
                it->gradient+=normal*(2/geometry.thickness[column])*geometry.derivative[level][m];
            }
            if(!transfer_positions.empty()){
                Matd current=Matd::Zero();for(const auto&e:rows[i])current+=(transfer_positions[e.j]-transfer_positions[i])*e.gradient.transpose();
                Matd delta=transfer_gradient[i]-current;
                if(delta.norm()>1e-10*std::max(Real(1),transfer_gradient[i].norm())){
                    if((delta*normal).norm()>5e-9)throw std::runtime_error("Fracture transfer would alter normal-column kinematics");
                    StdVec<size_t>ids;for(size_t n=0;n<rows[i].size();++n)if(nb.bondLive_[n]&&geometry.layer_of[rows[i][n].j]==level)ids.push_back(n);
                    Eigen::MatrixXd B(ids.size(),9),Y(ids.size(),3),G(ids.size(),3);Eigen::VectorXd W(ids.size());
                    for(size_t a=0;a<ids.size();++a){const auto&e=rows[i][ids[a]];B.row(a)=basis9(L*(X[e.j]-X[i])).transpose();Y.row(a)=(transfer_positions[e.j]-transfer_positions[i]).transpose();G.row(a)=e.gradient.transpose();W[a]=nb.W_ij_[ids[a]]*volume[e.j]/sum;}
                    Eigen::MatrixXd target=Y.transpose()*G+delta;
                    try
                    {
                        try
                        {
                            peridynamics::transferReconstruction(B,W,Y,G,target,1e-9);
                        }
                        catch(const peridynamics::UnresolvedReconstruction&)
                        {
                            // Keep the same live tangent neighbors and gradient
                            // target. A sparse crack-side row may only support
                            // affine moments after conservative state transfer.
                            peridynamics::transferReconstruction(B.leftCols(3),W,Y,G,target,1e-9);
                            if(quadratic)++affine_rows;
                            quadratic=false;
                        }
                    }
                    catch(const std::exception&e)
                    {
                        throw std::runtime_error("Gradient transfer particle "+std::to_string(i)+": "+e.what());
                    }
                    for(size_t a=0;a<ids.size();++a)rows[i][ids[a]].gradient=G.row(a).transpose();
                }
                Matd check=Matd::Zero();for(const auto&e:rows[i])check+=(transfer_positions[e.j]-transfer_positions[i])*e.gradient.transpose();
                if((check-transfer_gradient[i]).norm()>1e-8)throw std::runtime_error("Conservative fracture gradient transfer failed");
            }
            Matd reproduced=Matd::Zero();
            Eigen::Matrix<Real,9,3>second=Eigen::Matrix<Real,9,3>::Zero(),exact=Eigen::Matrix<Real,9,3>::Zero();exact.topRows<3>()=Matd::Identity()/dp;
            for(const auto&e:rows[i])
            {
                if(!e.gradient.allFinite())throw std::runtime_error("Nonfinite live fracture gradient");
                Vecd d=X[e.j]-X[i];reproduced+=d*e.gradient.transpose();second+=basis9(d/dp)*e.gradient.transpose();maximum_gradient_norm=std::max(maximum_gradient_norm,e.gradient.norm());
            }
            Real affine_error=(reproduced-Matd::Identity()).norm();maximum_affine_error=std::max(maximum_affine_error,affine_error);
            if(!std::isfinite(affine_error)||affine_error>1e-8)throw std::runtime_error("Fractured layered gradient failed Cartesian affine/rigid-motion reproduction");
            if(quadratic){Real e=(second-exact).norm()/exact.norm();maximum_quadratic_error=std::max(maximum_quadratic_error,e);if(e>1e-8)throw std::runtime_error("Fractured layered gradient failed requested quadratic reproduction");}
            degree[i]=reproduction[i]=quadratic?2:1;bond_count[i]=int(live_count);broken_fraction[i]=rows[i].empty()?0:1-Real(live_count)/rows[i].size();
        }
        max_reproduction_error=std::max(maximum_affine_error,maximum_quadratic_error);
        updateBound();
        for(size_t i=0;i<rows.size();++i)for(size_t n=0;n<rows[i].size();++n)cached_live_[i][n]=inner_configuration_[i].bondLive_[n];
        geometry_changed_=false;pending_bonds_=false;full_geometry_rebuild_=false;committed_facets_=facets_.size();reconstruction_pending_=false;++topology_version;
    }
    size_t topologyVersion()const override{return topology_version;}
    void checkIntact(size_t i)const override
    {
        requireInitialized();if(topologyPending())throw std::runtime_error("Layered fracture topology is pending reconstruction");const auto&nb=inner_configuration_.at(i);
        if(nb.current_size_!=rows.at(i).size())throw std::runtime_error("Layered fracture neighborhood changed");
        for(size_t n=0;n<rows[i].size();++n)
            if(nb.j_[n]!=rows[i][n].j||static_cast<unsigned char>(nb.bondLive_[n])!=cached_live_[i][n])
                throw std::runtime_error("Fracture topology changed without rebuilding the layered gradient");
    }

  private:
    bool initialized_=false,geometry_changed_=false,pending_bonds_=false,reconstruction_pending_=false;bool full_geometry_rebuild_=false;size_t committed_facets_=0;Real maximum_trace_support_=0;
    StdVec<CrackFacet> facets_;
    StdVec<StdVec<unsigned char>> cached_live_;
    Real*broken_fraction;int*reproduction,*bond_count;
    void requireInitialized()const{if(!initialized_)throw std::logic_error("Initialize the layered fracture reference first");}
    void validateFacet(const CrackFacet&f)const
    {
        constexpr Real normal_tolerance=1e-10;
        if((f.kind!=CrackFacet::Kind::PlaneDisk&&f.kind!=CrackFacet::Kind::Triangle)||
           !f.origin.allFinite()||!f.normal.allFinite()||
           std::abs(f.normal.stableNorm()-1)>normal_tolerance||!(f.radius>0))
            throw std::invalid_argument("Invalid restored crack facet kind, origin, unit normal or radius");
        for(const auto&v:f.vertex)if(!v.allFinite())
            throw std::invalid_argument("Nonfinite restored crack facet vertex");
        if(f.kind==CrackFacet::Kind::Triangle)
        {
            const Vecd u=f.vertex[1]-f.vertex[0],v=f.vertex[2]-f.vertex[0],n=u.cross(v);
            const Vecd origin_offset=f.origin-f.vertex[0];
            const Real scale=std::max(u.squaredNorm(),v.squaredNorm()),norm=n.stableNorm();
            if(!origin_offset.allFinite()||!std::isfinite(scale)||!std::isfinite(norm)||scale<=0||
               norm<=128*std::numeric_limits<Real>::epsilon()*scale||
               std::abs(std::abs(f.normal.dot(n/norm))-1)>normal_tolerance||
               std::abs(f.normal.dot(origin_offset))>tolerance(f,f.vertex[0],f.origin))
                throw std::invalid_argument("Restored crack triangle is degenerate or disagrees with its plane");
        }
    }
    void validateTopologyState(const TopologyState&t)const
    {
        if((!t.transfer_positions.empty()||!t.transfer_gradient.empty()||!t.transfer_velocity.empty())&&(t.transfer_positions.size()!=rows.size()||t.transfer_gradient.size()!=rows.size()||t.transfer_velocity.size()!=rows.size()))throw std::invalid_argument("Restored transfer-state dimensions");
        for(size_t i=0;i<t.transfer_positions.size();++i)if(!t.transfer_positions[i].allFinite()||!t.transfer_velocity[i].allFinite()||!t.transfer_gradient[i].allFinite()||t.transfer_gradient[i].determinant()<=0)throw std::invalid_argument("Invalid restored conservative transfer state");
        if(t.live.size()!=rows.size())throw std::invalid_argument("Fracture topology restore size");
        // Check every row size before indexing a reciprocal row supplied by t.
        for(size_t i=0;i<rows.size();++i)if(t.live[i].size()!=rows[i].size())
            throw std::invalid_argument("Fracture row restore size");
        for(const auto&f:t.facets)validateFacet(f);
        size_t broken_pairs=0;
        for(size_t i=0;i<rows.size();++i)for(size_t n=0;n<rows[i].size();++n)
        {
            const auto&e=rows[i][n];
            if(t.live[i][n]>1||t.live[i][n]!=t.live[e.j][e.reverse])
                throw std::invalid_argument("Nonreciprocal restored fracture topology");
            if(!t.live[i][n]&&geometry.column_of[i]==geometry.column_of[e.j])
                throw std::invalid_argument("Restored fracture cuts a normal thickness column");
            if(i>=e.j)continue;
            if(!t.live[i][n])++broken_pairs;
            for(const auto&f:t.facets)if(intersects(f,X[i],X[e.j],true)&&t.live[i][n])
                throw std::invalid_argument("Restored live bond crosses its recorded crack geometry");
        }
        if(broken_pairs!=t.physical_broken_bonds)
            throw std::invalid_argument("Restored broken-bond pair count disagrees with its reciprocal masks");
        // No reverse inference is made: manual severBondPair precracks may
        // have broken masks without any recorded geometric facets.
    }
    Real tolerance(const CrackFacet&f,const Vecd&a,const Vecd&b)const
    {return 128*std::numeric_limits<Real>::epsilon()*std::max({sph_body_->getSPHAdaptation().ReferenceSpacing(),(b-a).norm(),(a-f.origin).norm(),(b-f.origin).norm()});}
    bool contains(const CrackFacet&f,const Vecd&p,Real tol)const
    {
        if(f.kind==CrackFacet::Kind::PlaneDisk)
        {Vecd q=p-f.origin;q-=f.normal*f.normal.dot(q);return std::isinf(f.radius)||q.norm()<=f.radius+tol;}
        Vecd u=f.vertex[1]-f.vertex[0],v=f.vertex[2]-f.vertex[0],q=p-f.vertex[0];
        Real uu=u.squaredNorm(),vv=v.squaredNorm(),uv=u.dot(v),den=uu*vv-uv*uv;
        Real b=(vv*q.dot(u)-uv*q.dot(v))/den,c=(uu*q.dot(v)-uv*q.dot(u))/den;
        Real eps=tol/std::max(std::sqrt(std::max(uu,vv)),std::numeric_limits<Real>::min());
        return b>=-eps&&c>=-eps&&b+c<=1+eps;
    }
    bool intersects(const CrackFacet&f,const Vecd&a,const Vecd&b,bool reject_node)const
    {
        Real tol=tolerance(f,a,b),da=f.normal.dot(a-f.origin),db=f.normal.dot(b-f.origin),length=(b-a).norm();
        if(reject_node&&((std::abs(da)<=tol&&contains(f,a,tol))||(std::abs(db)<=tol&&contains(f,b,tol))))
            throw std::runtime_error("Precrack surface passes through a material node; provide split nodes or a crack between nodes rather than a hidden geometric offset");
        if(length<=tol||std::abs(da-db)<=tol)return false;
        if((da>tol&&db>tol)||(da< -tol&&db< -tol))return false;
        Real t=da/(da-db),eps=tol/length;
        if(t<=eps||t>=1-eps)return false;
        return contains(f,a+t*(b-a),tol);
    }
    size_t insertFacet(CrackFacet f,bool rebuild)
    {
        requireInitialized();StdVec<std::pair<size_t,size_t>>cut;
        // Validate every geometry/topology event before changing any live mask.
        for(size_t i=0;i<rows.size();++i)for(size_t n=0;n<rows[i].size();++n)
        {
            const auto&e=rows[i][n];if(i>=e.j)continue;
            if(!intersects(f,X[i],X[e.j],true))continue;
            if(geometry.column_of[i]==geometry.column_of[e.j])
                throw std::runtime_error("Precrack intersects a thickness column; this requires explicit split-column kinematics");
            if(live(i,n))cut.push_back({i,n});
        }
        // A public single-facet insertion is transactional. Rank loss in a
        // proposed crack must not leave half-written stencils or live masks.
        TopologyState before;if(rebuild)before=captureTopology();
        try{
            f.id=facets_.empty()?0:facets_.back().id+1;facets_.push_back(f);geometry_changed_=true;
            for(auto [i,n]:cut)severBondPair(i,n);
            if(rebuild)updateAfterFailure();
        }catch(...){if(rebuild)restoreTopology(before);throw;}
        return cut.size();
    }
    void updateBound()
    {
        StdVec<Real>stiffness(rows.size(),0);auto&mat=DynamicCast<HughesWingetSolid>(this,sph_body_->getMatterMaterial());Real rho=mat.ReferenceDensity(),C=3*mat.getBulkModulusforHGC(3);
        for(size_t i=0;i<rows.size();++i)
        {
            Vecd self=Vecd::Zero();Real sum=0;for(const auto&e:rows[i]){self-=e.gradient;sum+=e.gradient.norm();}
            sum+=self.norm();Real factor=volume[i]*C*sum;stiffness[i]+=factor*self.norm();for(const auto&e:rows[i])stiffness[e.j]+=factor*e.gradient.norm();
        }
        Real bound=0;for(size_t i=0;i<rows.size();++i)bound=std::max(bound,stiffness[i]/(rho*volume[i]));
        if(!std::isfinite(bound)||bound<=0)throw std::runtime_error("Invalid live layered spectral time bound");
        elastic_dt_bound=.2/std::sqrt(bound);
    }
};
}
