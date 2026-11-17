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
                        Real spacing,Real support_ratio=2.5,unsigned order=2)
      :ReferenceSolidSurface(p,targets,order,DeferredBuildTag{}),reference_(r),
       anchors_(target_anchors),original_targets_(targets),dp_(spacing),support_(spacing*support_ratio)
    {
        if(spacing<=0||support_ratio<=1||!std::isfinite(support_)||order<1||order>2)
            throw std::invalid_argument("Invalid one-sided crack trace support/order");
        if(&p.body!=&r.getSPHBody()||targets.size()!=target_anchors.size()||r.X.size()!=p.count||r.rows.size()!=p.count||r.topologyPending())
            throw std::invalid_argument("Crack trace requires initialized/current same-body topology and one material anchor per target");
        for(size_t i=0;i<p.count;++i)
            if((p.reference[i]-r.X[i]).norm()>128*std::numeric_limits<Real>::epsilon()*std::max(dp_,r.X[i].norm()))
                throw std::invalid_argument("Force-port and fracture reference coordinates differ");
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
        rebuildAll();
    }
    size_t materialAnchor(size_t q)const{return anchors_.at(q);}
    size_t builtTopologyVersion()const{return version_;}
    bool isCurrent()const{return version_==reference_.topology_version&&!reference_.topologyPending();}
    size_t mappingVersion()const override{return version_;}
    void ensureCurrent()const override
    {if(!isCurrent())throw std::runtime_error("One-sided crack surface was not rebuilt after topology change");}
    void rebuildAll()
    {
        requireReferenceCurrent();StdVec<size_t>ids(points.size());for(size_t q=0;q<ids.size();++q)ids[q]=q;rebuild(ids);
    }
    void rebuildChanged()
    {
        requireReferenceCurrent();if(version_==reference_.topology_version)return;
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
    StdVec<Vecd>original_targets_;
    Real dp_,support_;
    size_t version_=~size_t(0);
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
    void buildRow(size_t q,Row&row,FitInfo&info)const
    {
        size_t anchor=anchors_[q];
        if(!reference_.visibleReference(reference_.X[anchor],points[q]))
            throw std::runtime_error("Crack trace target lies across a crack from its own material anchor");
        // Geometric visibility and local LIVE connectivity are both needed:
        // severBondPair() may legitimately record topology before an explicit
        // surface is available. Never average a disconnected opposite side.
        StdVec<unsigned char>member(port.count,0),connected(port.count,0);
        for(size_t j:candidate_[q])member[j]=1;
        StdVec<size_t>queue{anchor};connected[anchor]=1;
        for(size_t k=0;k<queue.size();++k){size_t i=queue[k];const auto&nb=reference_.neighbors(i);for(size_t n=0;n<nb.current_size_;++n){size_t j=nb.j_[n];if(reference_.live(i,n)&&member[j]&&!connected[j]&&reference_.visible(i,j)){connected[j]=1;queue.push_back(j);}}}
        StdVec<size_t>nodes;StdVec<Real>weight;Real total=0;Matd covariance=Matd::Zero();
        for(size_t j:candidate_[q])
        {
            if(!connected[j]||!reference_.visible(anchor,j)||!reference_.visibleReference(points[q],reference_.X[j]))continue;
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
        // Check the committed row against current visibility, independently of
        // polynomial reproduction. No weight clipping or post-normalization.
        for(size_t j:row.index)if(!connected[j]||!reference_.visible(anchor,j)||!reference_.visibleReference(points[q],reference_.X[j]))
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
        version_=reference_.topology_version;
    }
};
}
