#pragma once
#include "pd_broken_bond_energy.h"
#include "pd_fracture_visible_surface.h"
#include "pd_cohesive_surface_evolution.h"
#include "pd_isoparametric_crack_surface.h"
#include <memory>
#include <functional>
#include <fstream>
#include <set>

namespace SPH::peridynamics {
/** Strength-triggered extrinsic cohesive surfaces on an explicit candidate bank.
 * NOSB constitutive elasticity and physical masses are unchanged. Actual cross-
 * surface bonds are removed once, then one-sided material/stabilizer/trace maps
 * are rebuilt. Surface work alone pays Gc; no point crack-band damage is mixed in.
 * Crack path is restricted by the supplied geometrical bank, NOT prescribed by
 * a time/damage schedule. Equal mixed-mode strength/toughness and a through-
 * thickness shell facet bank are explicit modelling assumptions.
 *
 * Every insertion is transactional, including external tie/platen maps via the
 * refresh callback. Its signed total energy jump is measured, never canceled
 * by offsets, velocity/force changes or extra fracture work. A jump larger than
 * the configured numerical acceptance fraction of Gc*new_area is rejected
 * for the physical residual. Removed numerical LS energy is separately
 * recorded as stabilizer transfer loss, never as physical Gc work.
 */
class ExtrinsicCrackGrowth {
 public:
    using Reference=solid_dynamics::PDLayeredFractureReference;
    using Bulk=solid_dynamics::PDBrokenBondEnergy;
    using Field=std::vector<Vecd>;using History=CohesiveInterfaceIMEX::History;
    struct Candidate {
        uint64_t id=0;std::array<Vecd,4>vertex;std::vector<std::array<Vecd,3>>barrier;
        std::vector<Vecd>point;std::vector<Real>area;
        std::vector<std::pair<size_t,size_t>>stress_pairs;
        Vecd normal=Vecd::Zero();
    };
    struct Patch {
        size_t candidate=0;
        std::unique_ptr<VisibleCrackSurface>minus_material,plus_material;
        std::unique_ptr<BilinearCrackSurface>minus,plus;
    };
    std::unique_ptr<CrackSurfaceContactNetwork>contact_network;
    Reference&reference;Bulk&bulk;SolidForcePort&port;
    Real ft,Gc,delta_f,spacing,contact_K;
    std::vector<Candidate>candidates;std::vector<Patch>patches;
    std::vector<History>history;Field last_traction;std::vector<unsigned char>active;
    Real maximum_insertion_fraction=.01,insertion_defect=0,max_insertion_fraction=0;
    Real stored=0,dissipated=0,initial_dissipation=0,contact_energy=0,stabilizer_transfer_loss=0;
    Real max_gap=0,max_penetration=0,contact_dt=std::numeric_limits<Real>::infinity();
    size_t insertion_events=0,rejected_events=0,failed_quadrature=0,unresolved_candidates=0;bool last_candidate_unresolved=false;
    std::string last_rejection;
    struct InsertionAudit {size_t candidate=0;Real area=0,old_energy=0,new_energy=0,old_bulk=0,new_bulk=0,old_cohesive=0,new_cohesive=0,old_contact=0,new_contact=0,new_dissipation=0,maximum_birth_gap=0,maximum_birth_velocity_jump=0,defect=0,fraction=0,old_stabilizer=0,new_stabilizer=0,total_defect=0;} last_insertion;
    ExtrinsicCrackGrowth(Reference&r,Bulk&b,SolidForcePort&p,Real strength,Real energy,Real dp,Real K,const std::vector<Candidate>&bank)
      :reference(r),bulk(b),port(p),ft(strength),Gc(energy),delta_f(2*energy/strength),spacing(dp),contact_K(K),candidates(bank),active(bank.size(),0){
        if(std::min({ft,Gc,spacing,contact_K})<=0)throw std::invalid_argument("Invalid extrinsic fracture inputs");
        for(const auto&c:bank){if(c.point.empty()||c.point.size()!=c.area.size()||c.barrier.empty()||c.stress_pairs.empty()||std::abs(c.normal.norm()-1)>1e-8)throw std::invalid_argument("Invalid crack candidate bank");for(auto ij:c.stress_pairs)if(ij.first>=port.count||ij.second>=port.count)throw std::invalid_argument("Crack strength sample index");}
    }
    Field positions()const{return Field(port.position,port.position+port.count);}
    void rebuildTraces(){for(auto&p:patches){p.minus->rebuildChanged();p.plus->rebuildChanged();}}
    std::vector<CohesiveSurfaceEvolution::Quadrature>quadrature()const {
        std::vector<CohesiveSurfaceEvolution::Quadrature>q;
        for(const auto&p:patches){const auto&c=candidates[p.candidate];for(size_t k=0;k<c.point.size();++k){CohesiveSurfaceEvolution::Quadrature a;a.area=c.area[k];a.trace={CohesiveSurfaceEvolution::TracePoint{p.plus.get(),4+k},{p.minus.get(),4+k},{p.minus.get(),4+k},{p.minus.get(),p.minus->normalU(k)},{p.minus.get(),p.minus->normalV(k)}};q.push_back(a);}}
        return q;
    }
    void statistics(const Field&q){stored=dissipated=max_gap=0;failed_quadrature=0;CohesiveSurfaceEvolution model(ft,Gc,quadrature());
        for(size_t k=0;k<history.size();++k){auto g=model.geometry(k,q);Real A=model.quadrature[k].area;stored+=A*model.stored(g.jump,history[k].kappa);dissipated+=A*model.dissipation(history[k].kappa);max_gap=std::max(max_gap,std::max(g.normal_gap,Real(0)));failed_quadrature+=history[k].kappa>=delta_f;}}
    Real strengthRatio(size_t k)const {const auto&c=candidates[k];Real maximum=0;for(auto ij:c.stress_pairs){Matd P=.5*(bulk.P[ij.first]+bulk.P[ij.second]),F=.5*(bulk.F[ij.first]+bulk.F[ij.second]);Vecd n=F.determinant()*F.inverse().transpose()*c.normal;if(!n.allFinite()||n.norm()==0)throw std::runtime_error("Invalid current crack initiation normal");n.normalize();Vecd T=P*c.normal;Real normal=T.dot(n);if(normal>0)maximum=std::max(maximum,T.norm()/ft);}return maximum;}
    size_t select()const {size_t best=candidates.size();Real ratio=1;for(size_t k=0;k<candidates.size();++k)if(!active[k]){Real r=strengthRatio(k);if(r>=ratio){ratio=r;best=k;}}return best;}
    Patch makePatch(size_t k,const StdVec<Vecd>&minus_values={},const StdVec<Vecd>&plus_values={},const StdVec<Vecd>&minus_velocity={},const StdVec<Vecd>&plus_velocity={}){
        const auto&c=candidates.at(k);Patch p;p.candidate=k;StdVec<Vecd>points(c.vertex.begin(),c.vertex.end());
            auto anchors=[&](Real side){StdVec<size_t>a;for(const auto&x:points){Real best=std::numeric_limits<Real>::infinity();size_t owner=port.count;
                for(size_t j=0;j<port.count;++j)if(side*(reference.X[j]-x).dot(c.normal)>64*std::numeric_limits<Real>::epsilon()*spacing){Real d=(reference.X[j]-x).squaredNorm();if(d<best){best=d;owner=j;}}
                if(owner==port.count)throw std::runtime_error("Missing physical crack-side owner");a.push_back(owner);}return a;};
            StdVec<Vecd>minus_normals(points.size(),-c.normal),plus_normals(points.size(),c.normal);
            p.minus_material=std::make_unique<VisibleCrackSurface>(port,points,anchors(-1),reference,spacing,3.5,2,minus_normals,minus_values,minus_velocity);
            p.plus_material=std::make_unique<VisibleCrackSurface>(port,points,anchors(1),reference,spacing,3.5,2,plus_normals,plus_values,plus_velocity);
            p.minus=std::make_unique<BilinearCrackSurface>(*p.minus_material,c.point);
            p.plus=std::make_unique<BilinearCrackSurface>(*p.plus_material,c.point);
        return p;
    }
    void restorePatches(const std::vector<size_t>&order,const std::vector<History>&saved,const StdVec<StdVec<Vecd>>&minus_values={},const StdVec<StdVec<Vecd>>&plus_values={},const StdVec<StdVec<Vecd>>&minus_velocity={},const StdVec<StdVec<Vecd>>&plus_velocity={}){
        std::vector<unsigned char>seen(candidates.size(),0);size_t expected=0;
        for(size_t k:order){if(k>=candidates.size()||seen[k])throw std::invalid_argument("Invalid restored crack patch ordering");seen[k]=1;expected+=candidates[k].point.size();}
        if((!minus_velocity.empty()&&minus_velocity.size()!=order.size())||(!plus_velocity.empty()&&plus_velocity.size()!=order.size()))throw std::invalid_argument("Restored trace velocity layout");
        if(saved.size()!=expected)throw std::invalid_argument("Restored cohesive history dimensions");
        if(!reference.transfer_positions.empty()&&(minus_values.size()!=order.size()||plus_values.size()!=order.size()))throw std::invalid_argument("Missing persistent crack-trace transfer state");
        std::vector<Patch>next;for(size_t i=0;i<order.size();++i)next.push_back(makePatch(order[i],minus_values.empty()?StdVec<Vecd>{}:minus_values[i],plus_values.empty()?StdVec<Vecd>{}:plus_values[i],minus_velocity.empty()?StdVec<Vecd>{}:minus_velocity[i],plus_velocity.empty()?StdVec<Vecd>{}:plus_velocity[i]));patches=std::move(next);contact_network.reset();active=std::move(seen);history=saved;last_traction.assign(saved.size(),Vecd::Zero());statistics(positions());
    }
    static std::vector<Candidate>readBank(const std::string&file,const LayeredParticleGeometry&geometry){
        std::ifstream f(file);std::string magic;size_t count=0,columns=0;f>>magic>>count>>columns;
        if(magic!="PD_COHESIVE_BANK_V2"||columns!=geometry.columns.size()||count==0||count>1000000)throw std::runtime_error("Invalid physical crack bank schema/geometry");
        std::vector<Candidate>bank;std::set<uint64_t>ids;
        for(size_t i=0;i<count;++i){Candidate c;size_t nb=0,nq=0,ns=0;f>>c.id>>nb>>nq>>ns;
            if(!ids.insert(c.id).second||nb==0||nb>4||nq==0||nq>32||ns==0||ns>32)throw std::runtime_error("Invalid crack bank patch");
            auto point=[&](Vecd&v){f>>v[0]>>v[1]>>v[2];if(!f||!v.allFinite())throw std::runtime_error("Invalid crack bank coordinate");};
            for(auto&v:c.vertex)point(v);c.normal=(c.vertex[1]-c.vertex[0]).cross(c.vertex[2]-c.vertex[0]).normalized();
            if(!c.normal.allFinite())throw std::runtime_error("Degenerate crack-bank face");
            c.barrier.resize(nb);for(auto&t:c.barrier)for(auto&v:t)point(v);c.point.resize(nq);c.area.resize(nq);
            for(size_t q=0;q<nq;++q){point(c.point[q]);f>>c.area[q];if(!f||!std::isfinite(c.area[q])||c.area[q]<=0)throw std::runtime_error("Invalid crack area");}
            for(size_t j=0;j<ns;++j){size_t minus=0,plus=0;f>>minus>>plus;if(!f||minus>=geometry.reference.size()||plus>=geometry.reference.size()||geometry.layer_of[minus]!=geometry.layer_of[plus]||(geometry.reference[plus]-geometry.reference[minus]).dot(c.normal)<=0)throw std::runtime_error("Crack bank stress-sample side/layer mismatch");c.stress_pairs.push_back({minus,plus});}
            bank.push_back(c);
        }
        std::string extra;if(f>>extra)throw std::runtime_error("Trailing crack bank input");return bank;
    }
    /** totalEnergy returns bulk+surface+contact+external constraint recoverable
     * energy. q/v/plastic histories stay fixed throughout this topology event. */
    bool insert(size_t k,bool precrack,const std::function<void()>&refresh,const std::function<Real()>&totalEnergy,
                const std::function<void()>&prepare={},const std::function<void()>&rollback={}){
        if(k>=candidates.size()||active[k])throw std::invalid_argument("Repeated/invalid physical crack insertion");
        last_candidate_unresolved=false;const auto&c=candidates[k];auto topology=reference.captureTopology();auto old_history=history;auto old_traction=last_traction;size_t old_patch_count=patches.size();
        StdVec<StdVec<Vecd>>old_minus,old_plus,old_minus_velocity,old_plus_velocity;for(const auto&p:patches){old_minus.push_back(p.minus->transferValues());old_plus.push_back(p.plus->transferValues());old_minus_velocity.push_back(p.minus->transferVelocities());old_plus_velocity.push_back(p.plus->transferVelocities());}
        Real oldEnergy=totalEnergy(),oldD=dissipated,oldHG=bulk.stabilizationEnergy();last_insertion={};last_insertion.candidate=k;last_insertion.old_energy=oldEnergy;last_insertion.old_bulk=bulk.energy();last_insertion.old_cohesive=stored;last_insertion.old_contact=contact_energy;last_insertion.old_stabilizer=oldHG;Real A=0;for(Real a:c.area)A+=a;
        try{
            StdVec<Vecd>target(c.vertex.begin(),c.vertex.end());StdVec<size_t>owner;
            for(const auto&x:target){size_t nearest=port.count;Real best=std::numeric_limits<Real>::infinity();for(size_t j=0;j<port.count;++j){Real d=(reference.X[j]-x).squaredNorm();if(d<best){best=d;nearest=j;}}owner.push_back(nearest);}
            VisibleCrackSurface continuous(port,target,owner,reference,spacing,3.5,2);StdVec<Vecd>common(target.size());for(size_t q=0;q<common.size();++q)common[q]=continuous.position(q);
            auto q=positions();StdVec<Matd>F(bulk.F,bulk.F+port.count);
            for(auto&p:patches){p.minus->captureTransferValues();p.plus->captureTransferValues();}
            if(prepare)prepare();reference.setTransferState(q,F,StdVec<Vecd>(port.velocity,port.velocity+port.count));
            size_t cut=0;for(const auto&t:c.barrier)cut+=reference.cutTriangle(t[0],t[1],t[2],false);
            if(cut==0){last_candidate_unresolved=true;throw std::runtime_error("Geometric facet has no remaining material interactions at this resolution");}
            reference.updateAfterFailure();bulk.rebuild();rebuildTraces();
            patches.push_back(makePatch(k,common,common));for(size_t i=0;i<c.point.size();++i){history.push_back({precrack?delta_f:0});last_traction.push_back(Vecd::Zero());}
            refresh();for(size_t q=0;q<4;++q)last_insertion.maximum_birth_velocity_jump=std::max(last_insertion.maximum_birth_velocity_jump,(patches.back().plus->velocity(q)-patches.back().minus->velocity(q)).norm());CohesiveSurfaceEvolution model(ft,Gc,quadrature());
            for(size_t j=0;j<history.size();++j){Real r=model.geometry(j,q).jump.norm();if(history[j].kappa<delta_f&&r>history[j].kappa+4e-14)throw std::runtime_error("New/reconstructed cohesive surface has an artificial birth opening");}
            statistics(q);Real after=totalEnergy(),newD=dissipated-oldD,newHG=bulk.stabilizationEnergy();
            Real total_defect=after-oldEnergy+(precrack?0:newD),defect=total_defect-(newHG-oldHG),fraction=std::abs(defect)/(Gc*A);
            last_insertion.area=A;last_insertion.new_energy=after;last_insertion.new_bulk=bulk.energy();last_insertion.new_cohesive=stored;last_insertion.new_contact=contact_energy;last_insertion.new_stabilizer=newHG;last_insertion.total_defect=total_defect;last_insertion.new_dissipation=newD;last_insertion.defect=defect;last_insertion.fraction=fraction;
            for(size_t j=old_history.size();j<history.size();++j)last_insertion.maximum_birth_gap=std::max(last_insertion.maximum_birth_gap,model.geometry(j,q).jump.norm());
            if(!precrack&&newHG>oldHG+1e-10*Gc*A)throw std::runtime_error("Removing live LS stabilization increased its energy");
            if(!precrack&&fraction>maximum_insertion_fraction)throw std::runtime_error("Physical crack insertion energy defect exceeds its unchanged numerical acceptance bound (fraction="+std::to_string(fraction)+")");
            active[k]=1;++insertion_events;if(precrack)initial_dissipation+=newD;else{insertion_defect+=defect;stabilizer_transfer_loss+=oldHG-newHG;max_insertion_fraction=std::max(max_insertion_fraction,fraction);}return true;
        }catch(const std::exception&e){last_rejection=e.what();if(last_candidate_unresolved){++unresolved_candidates;active[k]=2;}else ++rejected_events;patches.resize(old_patch_count);history=std::move(old_history);last_traction=std::move(old_traction);
            for(size_t i=0;i<old_patch_count;++i){patches[i].minus->restoreTransferValues(old_minus[i],old_minus_velocity[i]);patches[i].plus->restoreTransferValues(old_plus[i],old_plus_velocity[i]);}
            if(rollback)rollback();reference.restoreTopology(topology);bulk.rebuild();rebuildTraces();refresh();statistics(positions());totalEnergy();return false;}
    }
    Field cohesiveForces(const Field&q)const {
        Field f(q.size(),Vecd::Zero());CohesiveSurfaceEvolution model(ft,Gc,quadrature());
        for(size_t k=0;k<history.size();++k){Real kap=history[k].kappa;auto g=model.geometry(k,q,kap==0?model.absolute_tolerance+model.relative_tolerance*delta_f:0);Vecd t=Vecd::Zero();if(kap>0&&kap<delta_f)t=ft*(1/kap-1/delta_f)*g.jump;else if(kap==0&&k<last_traction.size())t=last_traction[k];for(const auto&e:g.linear.entries)f[e.node]-=model.quadrature[k].area*e.block().transpose()*t;}
        return f;
    }
    void applyContact(){
        contact_energy=max_penetration=0;contact_dt=std::numeric_limits<Real>::infinity();if(patches.empty()){contact_network.reset();return;}
        if(!contact_network||contact_network->sides.size()!=2*patches.size()){
            StdVec<CrackSurfaceContactNetwork::Side>sides;for(const auto&p:patches){const auto&c=candidates[p.candidate];sides.push_back({p.minus.get(),&c.area,1,p.candidate});sides.push_back({p.plus.get(),&c.area,-1,p.candidate});}
            contact_network=std::make_unique<CrackSurfaceContactNetwork>(sides,contact_K,spacing);
        }
        contact_network->apply();contact_energy=contact_network->energy;max_penetration=contact_network->max_penetration;contact_dt=contact_network->stable_dt;
    }

    CohesiveSurfaceEvolution::Result drift(const Field&qfree,const Field&qold,const Field&inverse_mass,Real dt){
        CohesiveSurfaceEvolution model(ft,Gc,quadrature());auto r=model.solve(qfree,qold,history,inverse_mass,dt);history=r.history;last_traction=r.traction;stored=r.stored;dissipated=r.dissipated;max_gap=0;failed_quadrature=0;
        for(size_t k=0;k<r.history.size();++k){failed_quadrature+=r.history[k].kappa>=delta_f;max_gap=std::max(max_gap,model.geometry(k,r.positions).normal_gap);}return r;
    }
};
}
