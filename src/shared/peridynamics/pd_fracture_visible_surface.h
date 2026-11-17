#pragma once
/** Material-anchored, one-sided value traces for fractured layered NOSB-PD.
 * Reference support is searched once with a spatial hash. After a crack event,
 * candidate nodes are filtered by the anchor's geometric/live-bond visibility.
 * A polynomial value functional is accepted only if it lies in the represented
 * moment range. Failure never creates an empty map or silently debonds a body.
 */
#include "pd_solid_reference_coupling.h"
#include "pd_layered_fracture_reference.h"
#include <array>
#include <unordered_map>

namespace SPH::peridynamics
{
class VisibleCrackSurface : public ReferenceSolidSurface
{
  public:
    StdVec<int> row_order;
    size_t affine_points=0,inactive_points=0;
    Real maximum_value_error=0,maximum_weight_l1=0;
    VisibleCrackSurface(SolidForcePort&p,const StdVec<Vecd>&targets,
                        const StdVec<size_t>&target_anchors,
                        solid_dynamics::PDLayeredFractureReference&r,
                        Real spacing,Real support_ratio=2.5,unsigned order=2,const StdVec<Vecd>&side_normals={},const StdVec<Vecd>&transfer_values={},const StdVec<Vecd>&transfer_velocities={})
      :ReferenceSolidSurface(p,targets,order,DeferredBuildTag{}),reference_(r),
       anchors_(target_anchors),original_targets_(targets),side_normals_(side_normals),transfer_values_(transfer_values),transfer_velocity_values_(transfer_velocities),dp_(spacing),support_(spacing*support_ratio)
    {
        if(spacing<=0||support_ratio<=1||!std::isfinite(support_)||order<1||order>2)
            throw std::invalid_argument("Invalid one-sided crack trace support/order");
        if(&p.body!=&r.getSPHBody()||targets.size()!=target_anchors.size()||r.X.size()!=p.count||r.rows.size()!=p.count||r.topologyPending())
            throw std::invalid_argument("Crack trace requires initialized/current same-body topology and one material anchor per target");
        for(size_t i=0;i<p.count;++i)
            if((p.reference[i]-r.X[i]).norm()>128*std::numeric_limits<Real>::epsilon()*std::max(dp_,r.X[i].norm()))
                throw std::invalid_argument("Force-port and fracture reference coordinates differ");
        if((!transfer_values_.empty()&&transfer_values_.size()!=targets.size())||(!transfer_velocity_values_.empty()&&transfer_velocity_values_.size()!=targets.size()))throw std::invalid_argument("Crack trace transfer-value dimensions");
        if(!side_normals_.empty()&&side_normals_.size()!=targets.size())throw std::invalid_argument("Crack trace side-normal dimensions");
        for(size_t q=0;q<side_normals_.size();++q)if(!side_normals_[q].allFinite()||std::abs(side_normals_[q].norm()-1)>1e-8||(p.reference[anchors_[q]]-targets[q]).dot(side_normals_[q])<=0)throw std::invalid_argument("Crack trace side does not contain its material anchor");
        std::unordered_map<Cell,StdVec<size_t>,Hash>cells;
        for(size_t j=0;j<p.count;++j)cells[cell(p.reference[j])].push_back(j);
        candidate_.resize(targets.size());affected_.resize(p.count);
        rows.resize(targets.size());row_order.assign(targets.size(),0);fit_info_.resize(targets.size());
        for(size_t q=0;q<targets.size();++q)
        {
            if(!targets[q].allFinite()||anchors_[q]>=p.count)
                throw std::invalid_argument("Invalid crack-trace target/material anchor");
            if((p.reference[anchors_[q]]-targets[q]).norm()>=support_)
                throw std::invalid_argument("Crack-trace material anchor lies outside its physical reference support");
            Cell center=cell(targets[q]);
            for(int a=-1;a<=1;++a)for(int b=-1;b<=1;++b)for(int c=-1;c<=1;++c)
            {
                auto it=cells.find(Cell{center[0]+a,center[1]+b,center[2]+c});if(it==cells.end())continue;
                for(size_t j:it->second)if((p.reference[j]-targets[q]).squaredNorm()<support_*support_)candidate_[q].push_back(j);
            }
            std::sort(candidate_[q].begin(),candidate_[q].end());
            candidate_[q].erase(std::unique(candidate_[q].begin(),candidate_[q].end()),candidate_[q].end());
            if(candidate_[q].empty())throw std::runtime_error("No material support for a crack trace");
            for(size_t j:candidate_[q])affected_[j].push_back(q);
        }
        reference_.registerTraceSupport(support_);
        rebuildAll();
    }
    size_t materialAnchor(size_t q)const{return anchors_.at(q);}
    size_t builtTopologyVersion()const{return version_;}
    bool isCurrent()const{return version_==reference_.topology_version&&!reference_.topologyPending()&&!transfer_pending_;}
    size_t mappingVersion()const override{return generation_;}
    void ensureCurrent()const override
    {if(!isCurrent())throw std::runtime_error("One-sided crack surface was not rebuilt after topology change");}
    const StdVec<Vecd>&transferValues()const{return transfer_values_;}
    const StdVec<Vecd>&transferVelocities()const{return transfer_velocity_values_;}
    void captureTransferValues(bool include_velocity=false){ensureCurrent();StdVec<Vecd>next(points.size()),vel;if(include_velocity)vel.resize(points.size());for(size_t q=0;q<points.size();++q){next[q]=position(q);if(include_velocity)vel[q]=velocity(q);}transfer_values_=std::move(next);transfer_velocity_values_=std::move(vel);transfer_pending_=true;}
    void restoreTransferValues(const StdVec<Vecd>&values,const StdVec<Vecd>&velocities={}){if((!values.empty()&&values.size()!=points.size())||(!velocities.empty()&&velocities.size()!=points.size()))throw std::invalid_argument("Restored trace transfer dimensions");transfer_values_=values;transfer_velocity_values_=velocities;transfer_pending_=true;}
    void rebuildAll()
    {
        requireReferenceCurrent();StdVec<size_t>ids(points.size());for(size_t q=0;q<ids.size();++q)ids[q]=q;rebuild(ids);
    }
    void rebuildChanged()
    {
        requireReferenceCurrent();if(version_==reference_.topology_version&&!transfer_pending_)return;
        if(transfer_pending_){rebuildAll();return;}
        if(version_+1!=reference_.topology_version){rebuildAll();return;}
        if(reference_.changed_particles.empty())throw std::runtime_error("Crack trace topology changed without an affected-owner list");
        StdVec<unsigned char>marked(points.size(),0);StdVec<size_t>ids;
        for(size_t j:reference_.changed_particles)
        {
            if(j>=affected_.size())throw std::runtime_error("Invalid changed material index for a crack trace");
            for(size_t q:affected_[j])marked[q]=1;
        }
        for(size_t q=0;q<marked.size();++q)if(marked[q])ids.push_back(q);
        rebuild(ids);
    }

  private:
    using Cell=std::array<long long,3>;
    struct Hash
    {size_t operator()(const Cell&a)const{size_t h=0;for(long long v:a)h^=std::hash<long long>{}(v)+0x9e3779b9+(h<<6)+(h>>2);return h;}};
    struct FitInfo {int order=0;Real affine=0,polynomial=0,value=0,condition=0,l1=0;};
    solid_dynamics::PDLayeredFractureReference&reference_;
    StdVec<size_t>anchors_;
    StdVec<Vecd>original_targets_,side_normals_,transfer_values_,transfer_velocity_values_;
    Real dp_,support_;
    size_t version_=~size_t(0),generation_=0;bool transfer_pending_=false;
    StdVec<StdVec<size_t>>candidate_,affected_;
    StdVec<FitInfo>fit_info_;
    Cell cell(const Vecd&x)const
    {
        Cell result;
        for(size_t d=0;d<3;++d)
        {
            Real v=std::floor(x[d]/support_);
            if(!std::isfinite(v)||v<=Real(std::numeric_limits<long long>::lowest()+1)||v>=Real(std::numeric_limits<long long>::max()-1))
                throw std::invalid_argument("Crack trace coordinates exceed spatial-hash index range");
            result[d]=static_cast<long long>(v);
        }
        return result;
    }
    void requireReferenceCurrent()const
    {
        if(reference_.topologyPending())throw std::runtime_error("Rebuild the fracture gradient before rebuilding its one-sided surface");
        if(points.size()!=original_targets_.size())throw std::runtime_error("Crack trace targets changed without constructing a new reference map");
        for(size_t q=0;q<points.size();++q)if((points[q]-original_targets_[q]).squaredNorm()!=0)
            throw std::runtime_error("Crack trace reference targets were modified in place");
    }
    template<int N>using Vector=Eigen::Matrix<Real,N,1>;
    template<int N>static Vector<N>basis(const Vecd&z)
    {
        Vector<N>b;b.template head<4>()<<1,z[0],z[1],z[2];
        if constexpr(N==10)b.template tail<6>()<<z[0]*z[0],z[1]*z[1],z[2]*z[2],z[0]*z[1],z[0]*z[2],z[1]*z[2];
        return b;
    }
    template<int N>bool fit(size_t q,const StdVec<size_t>&nodes,const StdVec<Real>&weight,
                            const Matd&L,Row&result,FitInfo&info)const
    {
        using Matrix=Eigen::Matrix<Real,N,N>;Matrix M=Matrix::Zero();
        for(size_t a=0;a<nodes.size();++a){auto b=basis<N>(L*((port.reference[nodes[a]]-points[q])/dp_));M+=weight[a]*b*b.transpose();}
        layered_detail::MomentInverse<N>inverse(M);Vector<N>e0=Vector<N>::Zero();e0[0]=1;
        Real functional_error=(M*inverse.inverse*e0-e0).norm();
        if(!std::isfinite(functional_error)||functional_error>1e-9)return false;
        Vector<N>coefficient=inverse.inverse*e0,polynomial=Vector<N>::Zero();
        if(!coefficient.allFinite())return false;
        result.index=nodes;result.weight.clear();Real sum=0,l1=0;Vecd affine=Vecd::Zero();
        for(size_t a=0;a<nodes.size();++a)
        {
            Vecd z=(port.reference[nodes[a]]-points[q])/dp_;
            Real w=weight[a]*basis<N>(L*z).dot(coefficient);
            if(!std::isfinite(w))return false;
            result.weight.push_back(w);sum+=w;l1+=std::abs(w);affine+=w*z;polynomial+=w*basis<N>(z);
        }
        Real ae=std::max(std::abs(sum-1),affine.norm()),pe=(polynomial-e0).cwiseAbs().maxCoeff();
        if(!std::isfinite(ae)||!std::isfinite(pe)||ae>1e-8||pe>1e-8)return false;
        info.order=N==10?2:1;info.affine=ae;info.polynomial=pe;info.value=functional_error;info.condition=inverse.condition;info.l1=l1;
        return true;
    }
    bool onSide(size_t q,size_t j)const{return side_normals_.empty()||(port.reference[j]-points[q]).dot(side_normals_[q])>=-128*std::numeric_limits<Real>::epsilon()*dp_;}
    void buildRow(size_t q,Row&row,FitInfo&info)const
    {
        size_t anchor=anchors_[q];
        if(!reference_.visibleReference(reference_.X[anchor],points[q]))
            throw std::runtime_error("Crack trace target lies across a crack from its own material anchor");
        // Geometric visibility and local LIVE connectivity are both needed:
        // severBondPair() may legitimately record topology before an explicit
        // surface is available. Never average a disconnected opposite side.
        StdVec<unsigned char>member(port.count,0),connected(port.count,0);
        for(size_t j:candidate_[q])if(onSide(q,j))member[j]=1;
        StdVec<size_t>queue{anchor};connected[anchor]=1;
        for(size_t k=0;k<queue.size();++k){size_t i=queue[k];const auto&nb=reference_.neighbors(i);for(size_t n=0;n<nb.current_size_;++n){size_t j=nb.j_[n];if(reference_.live(i,n)&&member[j]&&!connected[j]&&reference_.visible(i,j)){connected[j]=1;queue.push_back(j);}}}
        StdVec<size_t>nodes;StdVec<Real>weight;Real total=0;Matd covariance=Matd::Zero();
        for(size_t j:candidate_[q])
        {
            if(!onSide(q,j)||!connected[j]||!reference_.visible(anchor,j)||!reference_.visibleReference(points[q],reference_.X[j]))continue;
            Vecd z=(port.reference[j]-points[q])/dp_;Real w=std::exp(-z.squaredNorm())*port.volume[j]/std::pow(dp_,3);
            if(!std::isfinite(w)||w<=0)throw std::runtime_error("Nonpositive/underflowed crack trace quadrature weight");
            nodes.push_back(j);weight.push_back(w);total+=w;covariance+=w*z*z.transpose();
        }
        if(nodes.empty()||!std::isfinite(total)||total<=0)
            throw std::runtime_error("No visible one-sided material support; the interface was not silently released");
        for(Real&w:weight)w/=total;
        Matd L=layered_detail::whitening(covariance/total,1);
        bool accepted=reproduction_order==2&&fit<10>(q,nodes,weight,L,row,info);
        if(!accepted)accepted=fit<4>(q,nodes,weight,L,row,info);
        if(!accepted)throw std::runtime_error("One-sided crack trace value functional is unresolved at target "+std::to_string(q)+" (anchor "+std::to_string(anchor)+"); no bonding constraint or material was removed to bypass fit failure");
        if(!transfer_values_.empty()){
            if(reference_.transfer_positions.size()!=port.count)throw std::runtime_error("Trace transfer has no matching continuous-state field");
            const size_t modes=info.order==2?10:4,fields=transfer_velocity_values_.empty()?3:6;Real tau=reference_.transferVelocityTimeScale();if(fields==6&&reference_.transfer_velocity.size()!=port.count)throw std::runtime_error("Missing trace velocity transfer state");Eigen::MatrixXd B(nodes.size(),modes),Y(nodes.size(),fields),G(nodes.size(),1);Eigen::VectorXd W(nodes.size());
            for(size_t a=0;a<nodes.size();++a){Vecd z=L*((port.reference[nodes[a]]-points[q])/dp_);if(modes==10)B.row(a)=basis<10>(z).transpose();else B.row(a)=basis<4>(z).transpose();Y.block(a,0,1,3)=reference_.transfer_positions[nodes[a]].transpose();if(fields==6)Y.block(a,3,1,3)=tau*reference_.transfer_velocity[nodes[a]].transpose();G(a,0)=row.weight[a];W[a]=weight[a];}
            Eigen::MatrixXd target(fields,1);target.topRows(3)=transfer_values_[q];if(fields==6)target.bottomRows(3)=tau*transfer_velocity_values_[q];try{transferReconstruction(B,W,Y,G,target,2e-15);}catch(const std::exception&e){throw std::runtime_error("Surface transfer point "+std::to_string(q)+": "+e.what());}
            Real sum=0;Vecd affine=Vecd::Zero();Eigen::VectorXd polynomial=Eigen::VectorXd::Zero(modes);info.l1=0;
            for(size_t a=0;a<nodes.size();++a){row.weight[a]=G(a,0);sum+=row.weight[a];affine+=row.weight[a]*(port.reference[nodes[a]]-points[q])/dp_;polynomial+=row.weight[a]*B.row(a).transpose();info.l1+=std::abs(row.weight[a]);}
            polynomial[0]-=1;info.affine=std::max(std::abs(sum-1),affine.norm());info.polynomial=polynomial.norm();
            if(info.affine>1e-8||info.polynomial>1e-8)throw std::runtime_error("State transfer lost crack-trace polynomial consistency");
        }
        // Check the committed row against current visibility, independently of
        // polynomial reproduction. No weight clipping or post-normalization.
        for(size_t j:row.index)if(!onSide(q,j)||!connected[j]||!reference_.visible(anchor,j)||!reference_.visibleReference(points[q],reference_.X[j]))
            throw std::runtime_error("A crack trace interpolation row crosses its recorded crack geometry");
    }
    void rebuild(const StdVec<size_t>&ids)
    {
        // Two-phase commit leaves every old row untouched if any new fit fails.
        StdVec<Row>next(ids.size());StdVec<FitInfo>info(ids.size());
        for(size_t k=0;k<ids.size();++k)buildRow(ids[k],next[k],info[k]);
        for(size_t k=0;k<ids.size();++k){rows[ids[k]]=std::move(next[k]);fit_info_[ids[k]]=info[k];row_order[ids[k]]=info[k].order;}
        affine_points=inactive_points=0;affine_error=polynomial_error=largest_condition_number=maximum_value_error=maximum_weight_l1=0;
        for(const auto&f:fit_info_)
        {
            if(f.order==0)throw std::logic_error("A crack trace row has no accepted material value functional");
            affine_points+=f.order==1;affine_error=std::max(affine_error,f.affine);polynomial_error=std::max(polynomial_error,f.polynomial);
            largest_condition_number=std::max(largest_condition_number,f.condition);maximum_value_error=std::max(maximum_value_error,f.value);maximum_weight_l1=std::max(maximum_weight_l1,f.l1);
        }
        version_=reference_.topology_version;transfer_pending_=false;++generation_;
    }
};
}
