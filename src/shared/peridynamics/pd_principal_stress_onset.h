#pragma once
/** Strength-based irreversible bond events, with first-crack-only default.
 * Scan the symmetric average Cauchy stress of each undirected live bond.
 * Commit reciprocal irreversible bond failure only when explicitly requested.
 * No hidden stretch fallback, parameter calibration or fracture-energy claim.
 * Repeated events require the caller to rebuild damaged gradients and surface
 * maps. The default first-crack caller must stop before an intact-only update.
 */
#include "base_particle_dynamics.h"
#include "base_particles.hpp"
#include <Eigen/Eigenvalues>
#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <stdexcept>

namespace SPH::solid_dynamics
{
class PrincipalStressBondOnset : public LocalDynamics, public DataDelegateInner
{
  public:
    struct Candidate {size_t i,n,j,reverse;Real stress;};
    StdVec<Candidate> candidates;
    Real maximum_particle_tension=0,maximum_pair_upper_bound=0,maximum_triggered_pair_tension=0,maximum_damage=0;
    bool committed=false;
    bool repeatable_failure;
    size_t total_physical_broken_bonds=0;
    PrincipalStressBondOnset(BaseInnerRelation &relation,Real threshold,bool repeated=false)
      :LocalDynamics(relation.getSPHBody()),DataDelegateInner(relation),limit(threshold),
       stress(particles_->getVariableDataByName<Matd>("PDStress")),
       position(particles_->getVariableDataByName<Vecd>("Position")),
       damage(particles_->getVariableDataByName<Real>("Damage")),
       bond_count(particles_->getVariableDataByName<int>("BondNumber")),
       particle_live(particles_->getVariableDataByName<int>("ParticleLive")),repeatable_failure(repeated)
    {if(!std::isfinite(limit)||limit<=0)throw std::invalid_argument("Positive explicit tensile fracture threshold required");}

    bool scan()
    {
        if(committed&&!repeatable_failure)throw std::logic_error("First-crack diagnostic cannot advance a fractured configuration");
        candidates.clear();maximum_particle_tension=maximum_pair_upper_bound=maximum_triggered_pair_tension=0;
        size_t count=particles_->TotalRealParticles();StdVec<Real> maximum(count);
        for(size_t i=0;i<count;++i)
        {maximum[i]=principal(stress[i]);maximum_particle_tension=std::max(maximum_particle_tension,maximum[i]);}
        for(size_t i=0;i<count;++i)
        {
            const auto &row=inner_configuration_[i];
            for(size_t n=0;n<row.current_size_;++n)
            {
                size_t j=row.j_[n];if(i>=j||!row.bondLive_[n])continue;
                Real bound=.5*(maximum[i]+maximum[j]);maximum_pair_upper_bound=std::max(maximum_pair_upper_bound,bound);
                // Convexity of lambda_max makes this an exact rejection bound.
                if(bound<=limit)continue;
                Real value=principal(.5*(stress[i]+stress[j]));if(value<=limit)continue;
                const auto &other=inner_configuration_[j];size_t reverse=other.current_size_;
                for(size_t k=0;k<other.current_size_;++k)if(other.j_[k]==i){reverse=k;break;}
                if(reverse==other.current_size_||!other.bondLive_[reverse])throw std::runtime_error("First-crack diagnostic requires reciprocal live bonds");
                candidates.push_back({i,n,j,reverse,value});maximum_triggered_pair_tension=std::max(maximum_triggered_pair_tension,value);
            }
        }
        return !candidates.empty();
    }
    void commitFirstBreak()
    {
        if((committed&&!repeatable_failure)||candidates.empty())throw std::logic_error("No uncommitted bond-failure event");
        for(const auto &c:candidates){inner_configuration_[c.i].bondLive_[c.n]=false;inner_configuration_[c.j].bondLive_[c.reverse]=false;}
        for(size_t i=0;i<particles_->TotalRealParticles();++i)
        {
            const auto &row=inner_configuration_[i];size_t live=0;for(size_t n=0;n<row.current_size_;++n)live+=row.bondLive_[n];
            bond_count[i]=int(live);damage[i]=row.current_size_?1-Real(live)/Real(row.current_size_):0;
            particle_live[i]=live?1:0;maximum_damage=std::max(maximum_damage,damage[i]);
        }
        total_physical_broken_bonds+=candidates.size();committed=true;
    }
    void writeBonds(const std::string &path) const
    {
        std::ofstream f(path);f<<std::setprecision(17)<<"i,j,average_max_principal_Pa,xi_m,yi_m,zi_m,xj_m,yj_m,zj_m\n";
        for(const auto &c:candidates)f<<c.i<<','<<c.j<<','<<c.stress<<','<<position[c.i][0]<<','<<position[c.i][1]<<','<<position[c.i][2]<<','<<position[c.j][0]<<','<<position[c.j][1]<<','<<position[c.j][2]<<'\n';
    }
    Real threshold() const {return limit;}
  private:
    Real limit;Matd *stress;Vecd *position;Real *damage;int *bond_count,*particle_live;
    static Real principal(const Matd &s)
    {
        if(!s.allFinite())throw std::runtime_error("Non-finite stress in tensile onset diagnostic");
        Eigen::SelfAdjointEigenSolver<Matd> eig;eig.computeDirect(.5*(s+s.transpose()),Eigen::EigenvaluesOnly);
        if(eig.info()!=Eigen::Success)throw std::runtime_error("Principal stress computation failed");
        return eig.eigenvalues().maxCoeff();
    }
};
}
