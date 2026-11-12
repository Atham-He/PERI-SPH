#pragma once
#include "pd_quadratic_correspondence.h"
#include "pd_layered_geometry.h"

namespace SPH::solid_dynamics
{
/** Volume NOSB reference gradient for resolved thickness columns. Cartesian
 * quadratic completeness is retained. Normal polynomials are differentiated
 * along their actual quadrature columns rather than inferred from curvature
 * across a wide, thin Cartesian neighborhood. Forces use the existing EXACT
 * transpose of these rows; this does not overwrite one component of F. */
class PDLayeredQuadraticReference : public PDQuadraticReference
{
  public:
    const peridynamics::LayeredParticleGeometry &geometry;
    Real max_surface_functional_error=0;
    int minimum_surface_rank=9;
    PDLayeredQuadraticReference(BaseInnerRelation&r,const peridynamics::LayeredParticleGeometry&g):PDQuadraticReference(r),geometry(g){}
    void build() override
    {
        using namespace peridynamics::layered_detail;
        using V9=Eigen::Matrix<Real,9,1>;using M9=Eigen::Matrix<Real,9,9>;
        size_t count=particles_->TotalRealParticles();Real dp=sph_body_->getSPHAdaptation().ReferenceSpacing();rows.resize(count);
        if(count!=geometry.column_of.size())throw std::invalid_argument("Layered reference geometry count changed");
        for(size_t i=0;i<count;++i)
        {
            const auto&nb=inner_configuration_[i];size_t level=geometry.layer_of[i],column=geometry.column_of[i];Vecd normal=geometry.director[column];Matd T=Matd::Identity()-normal*normal.transpose();
            Real total=0;Matd covariance=Matd::Zero();rows[i].reserve(nb.current_size_);
            for(size_t n=0;n<nb.current_size_;++n)
            {
                size_t j=nb.j_[n];rows[i].push_back({j,0,Vecd::Zero()});
                if(geometry.layer_of[j]!=level)continue;
                Real w=nb.W_ij_[n]*volume[j];Vecd d=position[j]-position[i];total+=w;covariance+=w*d*d.transpose();
            }
            if(total<=0)throw std::runtime_error("No same-level surface support");
            Matd L=whitening(covariance/total,dp);M9 M=M9::Zero();
            for(size_t n=0;n<nb.current_size_;++n)if(geometry.layer_of[nb.j_[n]]==level)
            {V9 b=basis9(L*(position[nb.j_[n]]-position[i]));M+=nb.W_ij_[n]*volume[nb.j_[n]]/total*b*b.transpose();}
            MomentInverse<9> inverse(M);minimum_surface_rank=std::min(minimum_surface_rank,inverse.rank);max_moment_condition=std::max(max_moment_condition,inverse.condition);
            if(inverse.rank<5)throw std::runtime_error("Insufficient quadratic tangent support");
            Eigen::Matrix<Real,9,3>D=Eigen::Matrix<Real,9,3>::Zero();D.topRows<3>()=L*T;
            Real error=(M*inverse.inverse*D-D).norm()/D.norm();max_surface_functional_error=std::max(max_surface_functional_error,error);
            if(error>1e-9)throw std::runtime_error("Surface support cannot represent the requested tangential derivative");
            for(size_t n=0;n<nb.current_size_;++n)if(geometry.layer_of[nb.j_[n]]==level)
            {V9 b=basis9(L*(position[nb.j_[n]]-position[i])),q=inverse.inverse*b;rows[i][n].gradient=(nb.W_ij_[n]*volume[nb.j_[n]]/total)*T*L.transpose()*q.head<3>();}
            for(size_t m=0;m<geometry.columns[column].size();++m)if(m!=level)
            {
                size_t j=geometry.columns[column][m];auto it=std::find_if(rows[i].begin(),rows[i].end(),[&](const Edge&e){return e.j==j;});
                if(it==rows[i].end())throw std::runtime_error("Thickness-column node is outside the physical reference support");
                it->gradient+=normal*(2/geometry.thickness[column])*geometry.derivative[level][m];
            }
            Eigen::Matrix<Real,9,3>reproduced=Eigen::Matrix<Real,9,3>::Zero(),expected=Eigen::Matrix<Real,9,3>::Zero();expected.topRows<3>()=Matd::Identity()/dp;
            for(auto&e:rows[i])reproduced+=basis9((position[e.j]-position[i])/dp)*e.gradient.transpose();
            max_reproduction_error=std::max(max_reproduction_error,(reproduced-expected).norm()/expected.norm());
        }
        if(max_reproduction_error>1e-9)throw std::runtime_error("Layered Cartesian quadratic reproduction failed");
        // Same conservative graph bound as the original exact-adjoint operator.
        StdVec<Real>stiffness(count,0);auto&mat=DynamicCast<HughesWingetSolid>(this,sph_body_->getMatterMaterial());Real rho=mat.ReferenceDensity(),C=3*mat.getBulkModulusforHGC(3);
        for(size_t i=0;i<count;++i)
        {
            Vecd self=Vecd::Zero();Real sum=0;
            for(auto&e:rows[i]){auto&other=rows[e.j];auto it=std::find_if(other.begin(),other.end(),[&](const Edge&v){return v.j==i;});if(it==other.end())throw std::runtime_error("Asymmetric layered reference support");e.reverse=size_t(it-other.begin());self-=e.gradient;sum+=e.gradient.norm();}
            sum+=self.norm();Real factor=volume[i]*C*sum;stiffness[i]+=factor*self.norm();for(auto&e:rows[i])stiffness[e.j]+=factor*e.gradient.norm();
        }
        Real bound=0;for(size_t i=0;i<count;++i)bound=std::max(bound,stiffness[i]/(rho*volume[i]));elastic_dt_bound=.2/std::sqrt(bound);
    }
};
}
