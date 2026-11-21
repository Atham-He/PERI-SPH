#pragma once
#include "pd_broken_bond_energy.h"
#include "pd_fracture_visible_surface.h"
#include "pd_cohesive_surface_evolution.h"
#include "pd_fracture_patch_partition.h"
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
    CrackProjectionBatch *contact_projection_;

 public:
    enum class SurfaceLaw { Cohesive, ContactOnly };
    const SurfaceLaw surface_law;
    bool contactOnly() const { return surface_law == SurfaceLaw::ContactOnly; }
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
        std::shared_ptr<const FractureIntegrationPatch> integration;
        std::shared_ptr<VisibleCrackSurface>minus_material,plus_material;
        std::shared_ptr<BilinearCrackSurface>minus,plus;
    };
    std::unique_ptr<CrackSurfaceContactNetwork>contact_network;
    Reference&reference;Bulk&bulk;SolidForcePort&port;
    Real ft,Gc,delta_f,spacing,contact_K;
    std::vector<Candidate>candidates;std::vector<Patch>patches;
    bool affine_contact_traces=false;
    // Only retained until the first subsequent topology rebuild.
    StdVec<Vecd> affine_contact_transfer_positions;
    std::vector<History>history;Field last_traction;std::vector<unsigned char>active;
    Real maximum_insertion_fraction=.01,insertion_defect=0,max_insertion_fraction=0;
    Real stored=0,dissipated=0,initial_dissipation=0,contact_energy=0,stabilizer_transfer_loss=0;
    Real max_gap=0,max_penetration=0,contact_dt=std::numeric_limits<Real>::infinity();
    size_t insertion_events=0,rejected_events=0,failed_quadrature=0,unresolved_candidates=0;bool last_candidate_unresolved=false;
    std::string last_rejection;
    struct InsertionAudit {size_t candidate=0;Real area=0,old_energy=0,new_energy=0,old_bulk=0,new_bulk=0,old_cohesive=0,new_cohesive=0,old_contact=0,new_contact=0,new_dissipation=0,maximum_birth_gap=0,maximum_birth_velocity_jump=0,defect=0,fraction=0,old_stabilizer=0,new_stabilizer=0,total_defect=0;} last_insertion;
    ExtrinsicCrackGrowth(Reference&r,Bulk&b,SolidForcePort&p,Real strength,Real energy,Real dp,Real K,const std::vector<Candidate>&bank,
                         SurfaceLaw law=SurfaceLaw::Cohesive,CrackProjectionBatch *contact_projection=nullptr)
      :contact_projection_(contact_projection),surface_law(law),reference(r),bulk(b),port(p),ft(strength),Gc(energy),delta_f(2*energy/strength),spacing(dp),contact_K(K),candidates(bank),active(bank.size(),0){
        if(std::min({ft,Gc,spacing,contact_K})<=0)throw std::invalid_argument("Invalid extrinsic fracture inputs");
        for(const auto&c:bank){if(c.point.empty()||c.point.size()!=c.area.size()||c.barrier.empty()||c.stress_pairs.empty()||std::abs(c.normal.norm()-1)>1e-8)throw std::invalid_argument("Invalid crack candidate bank");for(auto ij:c.stress_pairs)if(ij.first>=port.count||ij.second>=port.count)throw std::invalid_argument("Crack strength sample index");}
    }
    Field positions()const{return Field(port.position,port.position+port.count);}
    void rebuildTraces()
    {
        for(auto &patch:patches)
        {
            for (auto side : {-1,1})
            {
                try { (side<0 ? patch.minus : patch.plus)->rebuildChanged(); }
                catch(const std::exception &e)
                {
                    throw std::runtime_error("crack "+std::to_string(candidates[patch.candidate].id)+
                                             " side "+std::to_string(side)+": "+e.what());
                }
            }
        }
    }
    const StdVec<Real> &integrationAreas(const Patch &patch) const
    {
        return patch.integration ? patch.integration->areas : candidates[patch.candidate].area;
    }
    FractureIntegrationPatch integrationPatch(const Patch &patch) const
    {
        const auto &candidate=candidates[patch.candidate];
        return patch.integration ? *patch.integration :
            FractureIntegrationPatch{candidate.vertex,candidate.point,candidate.area};
    }
    size_t integrationPointCount() const
    {
        size_t count=0;
        for (const auto &patch : patches)
            count+=integrationAreas(patch).size();
        return count;
    }
    std::vector<CohesiveSurfaceEvolution::Quadrature>quadrature()const {
        std::vector<CohesiveSurfaceEvolution::Quadrature>q;
        for(const auto&p:patches){const auto&c=candidates[p.candidate];for(size_t k=0;k<c.point.size();++k){CohesiveSurfaceEvolution::Quadrature a;a.area=c.area[k];a.trace={CohesiveSurfaceEvolution::TracePoint{p.plus.get(),4+k},{p.minus.get(),4+k},{p.minus.get(),4+k},{p.minus.get(),p.minus->normalU(k)},{p.minus.get(),p.minus->normalV(k)}};q.push_back(a);}}
        return q;
    }
    static Real maximumNormalOpening(const Patch&patch)
    {
        Real opening=0;
        const auto&minus=*patch.minus;
        const auto&plus=*patch.plus;
        for(size_t i=0;i<minus.natural.size();++i){
            const Vecd x=minus.position(4+i);
            const Vecd normal=(minus.position(minus.normalU(i))-x).cross(minus.position(minus.normalV(i))-x);
            const Real length=normal.norm();
            if(!std::isfinite(length)||length<=0)throw std::runtime_error("Crack face normal lost geometric rank");
            opening=std::max(opening,(plus.position(4+i)-x).dot(normal)/length);
        }
        return opening;
    }
    void statistics(const Field&q){stored=dissipated=max_gap=0;failed_quadrature=0;
        if(contactOnly()){
            for(const auto&patch:patches)max_gap=std::max(max_gap,maximumNormalOpening(patch));
            failed_quadrature=history.size();
            return;
        }
        CohesiveSurfaceEvolution model(ft,Gc,quadrature());
        for(size_t k=0;k<history.size();++k){auto g=model.kinematics(k,q);Real A=model.quadrature[k].area;stored+=A*model.stored(g.jump,history[k].kappa);dissipated+=A*model.dissipation(history[k].kappa);max_gap=std::max(max_gap,std::max(g.normal_gap,Real(0)));failed_quadrature+=history[k].kappa>=delta_f;}}
    Real strengthRatio(size_t k)const {const auto&c=candidates[k];Real maximum=0;for(auto ij:c.stress_pairs){Matd P=.5*(bulk.P[ij.first]+bulk.P[ij.second]),F=.5*(bulk.F[ij.first]+bulk.F[ij.second]);Vecd n=F.determinant()*F.inverse().transpose()*c.normal;if(!n.allFinite()||n.norm()==0)throw std::runtime_error("Invalid current crack initiation normal");n.normalize();Vecd T=P*c.normal;Real normal=T.dot(n);if(normal>0)maximum=std::max(maximum,T.norm()/ft);}return maximum;}
    size_t select()const {size_t best=candidates.size();Real ratio=1;for(size_t k=0;k<candidates.size();++k)if(!active[k]){Real r=strengthRatio(k);if(r>=ratio){ratio=r;best=k;}}return best;}
    Patch makePatch(size_t k,const StdVec<Vecd>&minus_values={},const StdVec<Vecd>&plus_values={},
                    const StdVec<Vecd>&minus_velocity={},const StdVec<Vecd>&plus_velocity={},
                    std::shared_ptr<const FractureIntegrationPatch> integration={}){
        const auto &c=candidates.at(k);
        Patch p;
        p.candidate=k;
        p.integration=std::move(integration);
        const auto &vertices=p.integration ? p.integration->vertices : c.vertex;
        const auto &quadrature_points=p.integration ? p.integration->points : c.point;
        const auto &areas=integrationAreas(p);
        if (quadrature_points.empty() || quadrature_points.size()!=areas.size() ||
            std::any_of(areas.begin(),areas.end(),[](Real area){return !std::isfinite(area)||area<=0;}))
            throw std::runtime_error("Invalid contact patch integration layout");
        StdVec<Vecd>points(vertices.begin(),vertices.end());
            const Vecd interior=(points[0]+points[1]+points[2]+points[3])*.25;
            auto anchors=[&](Real side){StdVec<size_t>a;for(const auto&x:points){Real best=std::numeric_limits<Real>::infinity();size_t owner=port.count;
                for(size_t j=0;j<port.count;++j)if(side*(reference.X[j]-x).dot(c.normal)>64*std::numeric_limits<Real>::epsilon()*spacing){Real d=(reference.X[j]-x).squaredNorm();if(d<best){best=d;owner=j;}}
                if(owner==port.count)throw std::runtime_error("Missing physical crack-side owner");a.push_back(owner);}return a;};
            StdVec<Vecd>minus_normals(points.size(),-c.normal),plus_normals(points.size(),c.normal);
            const unsigned order=affine_contact_traces?1:2;
            const auto*positions=affine_contact_traces?&affine_contact_transfer_positions:nullptr;
            using Ownership=VisibleCrackSurface::AnchorOwnership;
            using Transfer=VisibleCrackSurface::TransferOrder;
            p.minus_material=std::make_shared<VisibleCrackSurface>(port,points,anchors(-1),reference,spacing,3.5,order,minus_normals,minus_values,minus_velocity,Ownership::MaterialSide,Transfer::PreservePolynomial,positions,interior,bool(p.integration));
            p.plus_material=std::make_shared<VisibleCrackSurface>(port,points,anchors(1),reference,spacing,3.5,order,plus_normals,plus_values,plus_velocity,Ownership::MaterialSide,Transfer::PreservePolynomial,positions,interior,bool(p.integration));
            p.minus=std::make_shared<BilinearCrackSurface>(*p.minus_material,quadrature_points);
            p.plus=std::make_shared<BilinearCrackSurface>(*p.plus_material,quadrature_points);
        return p;
    }
    void restorePatches(const std::vector<size_t>&order,const std::vector<History>&saved,
                        const StdVec<StdVec<Vecd>>&minus_values={},const StdVec<StdVec<Vecd>>&plus_values={},
                        const StdVec<StdVec<Vecd>>&minus_velocity={},const StdVec<StdVec<Vecd>>&plus_velocity={},
                        const StdVec<FractureIntegrationPatch>&integrations={})
    {
        if(!integrations.empty() && (!contactOnly() || integrations.size()!=order.size()))
            throw std::invalid_argument("Invalid restored contact integration layout");
        if((!minus_velocity.empty()&&minus_velocity.size()!=order.size())||(!plus_velocity.empty()&&plus_velocity.size()!=order.size()))
            throw std::invalid_argument("Restored trace velocity layout");
        if ((!minus_values.empty()&&minus_values.size()!=order.size()) || (!plus_values.empty()&&plus_values.size()!=order.size()) ||
            (!reference.transfer_positions.empty()&&(minus_values.size()!=order.size()||plus_values.size()!=order.size())))
            throw std::invalid_argument("Missing or invalid persistent crack-trace transfer state");

        std::vector<unsigned char>seen(candidates.size(),0);
        std::vector<Patch>next;
        size_t expected=0;
        for(size_t i=0;i<order.size();++i)
        {
            const size_t k=order[i];
            const bool split=!integrations.empty()&&!integrations[i].areas.empty();
            if(k>=candidates.size()||(seen[k]&&(!split||seen[k]!=2)))
                throw std::invalid_argument("Invalid restored crack patch ordering");
            seen[k]=split?2:1;
            auto layout=split?std::make_shared<const FractureIntegrationPatch>(integrations[i]):nullptr;
            next.push_back(makePatch(k,minus_values.empty()?Field{}:minus_values[i],plus_values.empty()?Field{}:plus_values[i],
                minus_velocity.empty()?Field{}:minus_velocity[i],plus_velocity.empty()?Field{}:plus_velocity[i],std::move(layout)));
            expected+=integrationAreas(next.back()).size();
        }
        if(saved.size()!=expected)
            throw std::invalid_argument("Restored cohesive history dimensions");
        for (auto &flag : seen)
            flag=(flag!=0);
        patches=std::move(next);
        contact_network.reset();
        active=std::move(seen);
        history=saved;
        last_traction.assign(saved.size(),Vecd::Zero());
        statistics(positions());
    }

    static Vecd patchInterior(const std::array<Vecd,4> &vertices)
    {
        return .25*(vertices[0]+vertices[1]+vertices[2]+vertices[3]);
    }

    static Field transferPatchCorners(const std::array<Vecd,4> &parent,const Field &values,
                                     const std::array<Vecd,4> &child)
    {
        if (values.size()!=4)
            throw std::logic_error("Patch subdivision requires its four pre-cut trace values");
        const std::array<Vecd,4> state{values[0],values[1],values[2],values[3]};
        Field result;
        for (const auto &vertex : child)
        {
            const auto uv=BilinearCrackGeometry::closest(parent,vertex);
            result.push_back(BilinearCrackGeometry::evaluate(state,uv).x);
        }
        return result;
    }

    void splitExistingContactPatches(const StdVec<Reference::CrackFacet> &facets)
    {
        std::vector<Patch> next;
        next.reserve(patches.size());
        for (const auto &patch : patches)
        {
            const auto geometry=integrationPatch(patch);
            auto children=fracture_partition::split(geometry,facets);
            if (children.empty())
            {
                next.push_back(patch);
                continue;
            }
            for (auto &child : children)
            {
                const auto minus=transferPatchCorners(geometry.vertices,patch.minus->transferValues(),child.vertices);
                const auto plus=transferPatchCorners(geometry.vertices,patch.plus->transferValues(),child.vertices);
                next.push_back(makePatch(patch.candidate,minus,plus,{},{},
                    std::make_shared<const FractureIntegrationPatch>(std::move(child))));
            }
        }
        if (next.size()!=patches.size())
        {
            contact_network.reset();
            patches=std::move(next);
        }
    }

    void transferContactTracesToAffine()
    {
        if(!contactOnly()||affine_contact_traces)throw std::logic_error("Affine contact handoff requires contact-only legacy traces");
        // Match the current face corners, rather than rewinding to the last
        // topology event. Leave the bulk and all other trace targets intact.
        affine_contact_transfer_positions=positions();
        affine_contact_traces=true;
        std::vector<Patch>rebuilt;
        rebuilt.reserve(patches.size());
        for(const auto&patch:patches){
            Field minus(4),plus(4);
            for(size_t q=0;q<4;++q){minus[q]=patch.minus->position(q);plus[q]=patch.plus->position(q);}
            rebuilt.push_back(makePatch(patch.candidate,minus,plus,{},{},patch.integration));
        }
        contact_network.reset();
        patches.swap(rebuilt);
        statistics(affine_contact_transfer_positions);
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
        last_candidate_unresolved=false;const auto&c=candidates[k];auto topology=reference.captureTopology();auto old_history=history;auto old_traction=last_traction;const auto old_patches=patches;size_t old_patch_count=patches.size();
        StdVec<StdVec<Vecd>>old_minus,old_plus,old_minus_velocity,old_plus_velocity;for(const auto&p:patches){old_minus.push_back(p.minus->transferValues());old_plus.push_back(p.plus->transferValues());old_minus_velocity.push_back(p.minus->transferVelocities());old_plus_velocity.push_back(p.plus->transferVelocities());}
        Real oldEnergy=totalEnergy(),oldD=dissipated,oldHG=bulk.stabilizationEnergy();last_insertion={};last_insertion.candidate=k;last_insertion.old_energy=oldEnergy;last_insertion.old_bulk=bulk.energy();last_insertion.old_cohesive=stored;last_insertion.old_contact=contact_energy;last_insertion.old_stabilizer=oldHG;Real A=0;for(Real a:c.area)A+=a;
        auto old_contact_positions=affine_contact_transfer_positions;
        const char *phase="pre-insertion trace";
        try{
            struct NewPatchTrace
            {
                std::shared_ptr<const FractureIntegrationPatch> layout;
                Field common;
            };
            std::vector<NewPatchTrace> new_traces;
            auto capture=[&](const std::array<Vecd,4> &vertices,std::shared_ptr<const FractureIntegrationPatch> layout={})
            {
                const Field target(vertices.begin(),vertices.end());
                std::vector<size_t> owner;
                for (const auto &x : target)
                {
                    size_t nearest=port.count;
                    Real best=std::numeric_limits<Real>::infinity();
                    for (size_t j=0;j<port.count;++j)
                    {
                        const Real distance=(reference.X[j]-x).squaredNorm();
                        if (distance<best) { best=distance;nearest=j; }
                    }
                    owner.push_back(nearest);
                }
                VisibleCrackSurface continuous(port,target,owner,reference,spacing,3.5,2,{},{},{},
                    VisibleCrackSurface::AnchorOwnership::ReferencePoint,
                    VisibleCrackSurface::TransferOrder::PreservePolynomial,nullptr,
                    layout?std::optional<Vecd>(patchInterior(vertices)):std::nullopt,bool(layout));
                NewPatchTrace trace;
                trace.layout=std::move(layout);
                for (size_t corner=0;corner<target.size();++corner)
                    trace.common.push_back(continuous.position(corner));
                new_traces.push_back(std::move(trace));
            };
            StdVec<FractureIntegrationPatch> subdivisions;
            if (contactOnly())
                subdivisions=fracture_partition::split({c.vertex,c.point,c.area},reference.crackFacets());
            if (subdivisions.empty())
                capture(c.vertex);
            else
                for (auto &subdivision : subdivisions)
                {
                    auto layout=std::make_shared<const FractureIntegrationPatch>(std::move(subdivision));
                    capture(layout->vertices,layout);
                }
            auto q=positions();StdVec<Matd>F(bulk.F,bulk.F+port.count);
            for(auto&p:patches){p.minus->captureTransferValues();p.plus->captureTransferValues();}
            if(prepare)prepare();
            affine_contact_transfer_positions.clear();
            reference.setTransferState(q,F,StdVec<Vecd>(port.velocity,port.velocity+port.count));
            const auto prior_facets=reference.crackFacets().size();
            size_t cut=0;for(const auto&t:c.barrier)cut+=reference.cutTriangle(t[0],t[1],t[2],false);
            if(cut==0){last_candidate_unresolved=true;throw std::runtime_error("Geometric facet has no remaining material interactions at this resolution");}
            phase="material reconstruction";
            reference.updateAfterFailure();bulk.rebuild();
            phase="existing crack traces";
            if (contactOnly())
            {
                const auto &facets=reference.crackFacets();
                splitExistingContactPatches(StdVec<Reference::CrackFacet>(facets.begin()+prior_facets,facets.end()));
            }
            rebuildTraces();
            phase="new crack traces";
            const size_t first_new_patch=patches.size();
            for (auto &trace : new_traces)
                patches.push_back(makePatch(k,trace.common,trace.common,{},{},std::move(trace.layout)));
            if (contactOnly())
            {
                history.resize(integrationPointCount());
                last_traction.resize(history.size(),Vecd::Zero());
            }
            else
                for(size_t i=0;i<c.point.size();++i)
                {
                    history.push_back({precrack?delta_f:0});
                    last_traction.push_back(Vecd::Zero());
                }
            phase="external surfaces and constraints";
            refresh();
            phase="crack insertion energy and opening";
            for (size_t i=first_new_patch;i<patches.size();++i)
                for(size_t q=0;q<4;++q)
                    last_insertion.maximum_birth_velocity_jump=std::max(last_insertion.maximum_birth_velocity_jump,
                        (patches[i].plus->velocity(q)-patches[i].minus->velocity(q)).norm());
            if(contactOnly()){
                for (size_t i=first_new_patch;i<patches.size();++i)
                    for(size_t j=0;j<integrationAreas(patches[i]).size();++j)
                        last_insertion.maximum_birth_gap=std::max(last_insertion.maximum_birth_gap,
                            (patches[i].plus->position(4+j)-patches[i].minus->position(4+j)).norm());
                if(last_insertion.maximum_birth_gap>4e-14)throw std::runtime_error("New crack surface has an artificial birth opening");
            }else{
                CohesiveSurfaceEvolution model(ft,Gc,quadrature());
                for(size_t j=0;j<history.size();++j){Real r=model.kinematics(j,q).jump.norm();if(history[j].kappa<delta_f&&r>history[j].kappa+4e-14)throw std::runtime_error("New/reconstructed cohesive surface has an artificial birth opening");}
                for(size_t j=old_history.size();j<history.size();++j)last_insertion.maximum_birth_gap=std::max(last_insertion.maximum_birth_gap,model.kinematics(j,q).jump.norm());
            }
            statistics(q);Real after=totalEnergy(),newD=dissipated-oldD,newHG=bulk.stabilizationEnergy();
            Real total_defect=after-oldEnergy+(precrack?0:newD),defect=total_defect-(newHG-oldHG),fraction=std::abs(defect)/(Gc*A);
            last_insertion.area=A;last_insertion.new_energy=after;last_insertion.new_bulk=bulk.energy();last_insertion.new_cohesive=stored;last_insertion.new_contact=contact_energy;last_insertion.new_stabilizer=newHG;last_insertion.total_defect=total_defect;last_insertion.new_dissipation=newD;last_insertion.defect=defect;last_insertion.fraction=fraction;
            // Gc is an audit scale only when the crack carries no cohesive traction.
            if(!precrack&&newHG>oldHG+1e-10*Gc*A)throw std::runtime_error("Removing live LS stabilization increased its energy");
            if(!precrack&&fraction>maximum_insertion_fraction)throw std::runtime_error("Physical crack insertion energy defect exceeds its unchanged numerical acceptance bound (fraction="+std::to_string(fraction)+")");
            active[k]=1;++insertion_events;if(precrack)initial_dissipation+=newD;else{insertion_defect+=defect;stabilizer_transfer_loss+=oldHG-newHG;max_insertion_fraction=std::max(max_insertion_fraction,fraction);}return true;
        }catch(const std::exception&e){last_rejection="candidate "+std::to_string(c.id)+" "+phase+": "+e.what();if(last_candidate_unresolved){++unresolved_candidates;active[k]=2;}else ++rejected_events;patches=old_patches;contact_network.reset();history=std::move(old_history);last_traction=std::move(old_traction);
            for(size_t i=0;i<old_patch_count;++i){patches[i].minus->restoreTransferValues(old_minus[i],old_minus_velocity[i]);patches[i].plus->restoreTransferValues(old_plus[i],old_plus_velocity[i]);}
            affine_contact_transfer_positions=std::move(old_contact_positions);
            if(rollback)rollback();reference.restoreTopology(topology);bulk.rebuild();rebuildTraces();refresh();statistics(positions());totalEnergy();return false;}
    }
    Field cohesiveForces(const Field&q)const {
        Field f(q.size(),Vecd::Zero());
        if(contactOnly())return f;
        CohesiveSurfaceEvolution model(ft,Gc,quadrature());
        for(size_t k=0;k<history.size();++k){Real kap=history[k].kappa;auto g=model.geometry(k,q,kap==0?model.absolute_tolerance+model.relative_tolerance*delta_f:0);Vecd t=Vecd::Zero();if(kap>0&&kap<delta_f)t=ft*(1/kap-1/delta_f)*g.jump;else if(kap==0&&k<last_traction.size())t=last_traction[k];for(const auto&e:g.linear.entries)f[e.node]-=model.quadrature[k].area*e.block().transpose()*t;}
        return f;
    }
    void applyContact(){
        contact_energy=max_penetration=0;contact_dt=std::numeric_limits<Real>::infinity();if(patches.empty()){contact_network.reset();return;}
        if(!contact_network||contact_network->sides.size()!=2*patches.size()){
            StdVec<CrackSurfaceContactNetwork::Side>sides;
            for (const auto &patch : patches)
            {
                const auto &areas=integrationAreas(patch);
                sides.push_back({patch.minus.get(),&areas,1,patch.candidate});
                sides.push_back({patch.plus.get(),&areas,-1,patch.candidate});
            }
            contact_network=std::make_unique<CrackSurfaceContactNetwork>(sides,contact_K,spacing,contact_projection_);
        }
        contact_network->apply();contact_energy=contact_network->energy;max_penetration=contact_network->max_penetration;contact_dt=contact_network->stable_dt;
    }

    CohesiveSurfaceEvolution::Result trialDrift(const Field&qfree,const Field&qold,const Field&inverse_mass,Real dt)const
    {
        if(contactOnly())throw std::logic_error("Contact-only cracks use the explicit bulk drift");
        CohesiveSurfaceEvolution model(ft,Gc,quadrature());
        return model.solve(qfree,qold,history,inverse_mass,dt);
    }
    void acceptDrift(const CohesiveSurfaceEvolution::Result&r)
    {
        history=r.history;last_traction=r.traction;stored=r.stored;dissipated=r.dissipated;max_gap=0;failed_quadrature=0;
        CohesiveSurfaceEvolution model(ft,Gc,quadrature());
        for(size_t k=0;k<r.history.size();++k){
            failed_quadrature+=r.history[k].kappa>=delta_f;
            max_gap=std::max(max_gap,model.kinematics(k,r.positions).normal_gap);
        }
    }
    CohesiveSurfaceEvolution::Result drift(const Field&qfree,const Field&qold,const Field&inverse_mass,Real dt)
    {
        auto result=trialDrift(qfree,qold,inverse_mass,dt);
        acceptDrift(result);
        return result;
    }
};
}
