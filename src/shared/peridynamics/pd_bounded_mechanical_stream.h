/* Licensed under the Apache License, Version 2.0. */
#ifndef PD_BOUNDED_MECHANICAL_STREAM_H
#define PD_BOUNDED_MECHANICAL_STREAM_H
#include "pd_bounded_session_archive.h"
#include <optional>

namespace SPH::peridynamics
{
/** Borrow stored arrays, never construct a fresh reference or initial scene. */
struct BoundedStoredSceneView
{
    struct {int mode=1;} config;
    const BoundedGelModel &model;
    const ReferenceSurfaceTransfer &surface,&dense_surface;
    const Matd &frame;
    const std::array<Vecd,2> &normal;
    explicit BoundedStoredSceneView(const BoundedArchiveBody &body):model(body.model),surface(body.surface),
        dense_surface(body.dense_surface),frame(body.frame),normal(body.normal){}
};

/** Every accepted physical output quantity, separate from ephemeral Newton
 * trace/process metadata. Full point material history is stored in body.
 */
struct BoundedMechanicalObservation
{
    StdVec<Vecd> contact_impulse,normalized_momentum_residual;
    std::array<Vecd,2> plate_impulse{Vecd::Zero(),Vecd::Zero()},normal_impulse{Vecd::Zero(),Vecd::Zero()},
        tangent_impulse{Vecd::Zero(),Vecd::Zero()},reaction{Vecd::Zero(),Vecd::Zero()};
    std::array<Real,2> normal_work{0,0},tangent_work{0,0},whole_gap{0,0},endpoint_gap{0,0};
    Real kinetic_old=0,kinetic_new=0,material_old=0,material_new=0,nonaffine_old=0,nonaffine_new=0,
        Maxwell_loss=0,normal_loss=0,friction_loss=0,algorithmic_loss=0,energy_residual=0,
        material_work_residual=0,nonaffine_work_residual=0;
    Vecd total_momentum_residual=Vecd::Zero(),internal_force_sum=Vecd::Zero(),
        angular_drift_minus_external=Vecd::Zero(),internal_old_torque=Vecd::Zero();

    template<class Registry>void registerCheckpoint(Registry &cp)
    {
        cp.addArray("ContactImpulse",contact_impulse);cp.addArray("NormalizedMomentumResidual",normalized_momentum_residual);
        for(size_t p=0;p<2;++p)
        {
            auto plate=cp.scope(p==0?"lower":"upper");
            plate.addScalar("Impulse",plate_impulse[p]);plate.addScalar("NormalImpulse",normal_impulse[p]);
            plate.addScalar("TwoTangentImpulse",tangent_impulse[p]);plate.addScalar("IntervalMeanReaction",reaction[p]);
            plate.addScalar("NormalWork",normal_work[p]);plate.addScalar("TangentWork",tangent_work[p]);
            plate.addScalar("WholeCubicMinimumGap",whole_gap[p]);plate.addScalar("EndpointMinimumGap",endpoint_gap[p]);
        }
        cp.addScalar("KineticOld",kinetic_old);cp.addScalar("KineticNew",kinetic_new);
        cp.addScalar("MaterialStoredOld",material_old);cp.addScalar("MaterialStoredNew",material_new);
        cp.addScalar("NonaffineStoredOld",nonaffine_old);cp.addScalar("NonaffineStoredNew",nonaffine_new);
        cp.addScalar("MaxwellLoss",Maxwell_loss);cp.addScalar("NormalContactLoss",normal_loss);
        cp.addScalar("FrictionLoss",friction_loss);cp.addScalar("AlgorithmicLoss",algorithmic_loss);
        cp.addScalar("EnergyResidual",energy_residual);cp.addScalar("MaterialTransposeWorkResidual",material_work_residual);
        cp.addScalar("NonaffineWorkResidual",nonaffine_work_residual);
        cp.addScalar("TotalMomentumResidual",total_momentum_residual);cp.addScalar("InternalForceSum",internal_force_sum);
        cp.addScalar("AngularDriftMinusExternal",angular_drift_minus_external);cp.addScalar("OldInternalTorque",internal_old_torque);
    }
};

struct BoundedMechanicalRecord
{
    BoundedArchiveBody body;
    BoundedMechanicalObservation observation;
    uint64_t kind=0,interval_phase=0;
    Real old_time=0;
    template<class Registry>void registerCheckpoint(Registry &cp)
    {
        auto b=cp.scope("body");body.registerCheckpoint(b);auto o=cp.scope("observation");observation.registerCheckpoint(o);
        cp.addScalar("RecordKind",kind);cp.addScalar("IntervalPhase",interval_phase);cp.addScalar("OldTime",old_time);
    }
};

inline canonical_archive::Bytes encodeBoundedMechanicalRecord(BoundedMechanicalRecord record)
{
    canonical_archive::Builder builder(record.body.accepted.gel.steps+1);canonical_archive::FieldWriter fields(builder);
    record.registerCheckpoint(fields);return builder.finish();
}
inline BoundedMechanicalRecord decodeBoundedMechanicalRecord(const canonical_archive::Bytes &bytes,uint64_t expected_step)
{
    const canonical_archive::Candidate candidate(bytes,expected_step+1);
    return canonical_archive::decodeDetached<BoundedMechanicalRecord>(candidate,[](const auto &record){record.body.validateStored();});
}

inline BoundedMechanicalObservation boundedMechanicalObservation(const BoundedGelFrictionStepResult &r,
    Real dt,const std::array<Real,2> &whole,const std::array<Real,2> &endpoint)
{
    if(!r.accepted || !(dt>0))throw std::invalid_argument("mechanical stream: unaccepted proposal");
    BoundedMechanicalObservation o;o.contact_impulse=r.contact_impulse;o.normalized_momentum_residual=r.normalized_momentum_residual;
    o.plate_impulse=r.particle_plate_impulse;o.normal_impulse=r.normal_plate_impulse;o.tangent_impulse=r.tangent_plate_impulse;
    o.normal_work=r.normal_plate_work;o.tangent_work=r.tangent_plate_work;o.whole_gap=whole;o.endpoint_gap=endpoint;
    for(size_t p=0;p<2;++p)o.reaction[p]=-r.particle_plate_impulse[p]/dt;
    o.kinetic_old=r.kinetic_old;o.kinetic_new=r.kinetic_new;
    o.material_old=r.evaluation.material_stored_old;o.material_new=r.evaluation.material_stored_new;
    o.nonaffine_old=r.evaluation.stabilization_stored_old;o.nonaffine_new=r.evaluation.stabilization_stored_new;
    o.Maxwell_loss=r.evaluation.physical_material_loss;o.normal_loss=r.normal_contact_loss;o.friction_loss=r.friction_loss;
    o.algorithmic_loss=r.algorithmic_loss;o.energy_residual=r.energy_residual;
    o.material_work_residual=r.evaluation.transpose_work_residual;o.nonaffine_work_residual=r.evaluation.stabilization_work_residual;
    o.total_momentum_residual=r.total_momentum_residual;o.internal_force_sum=r.internal_force_sum;
    o.angular_drift_minus_external=r.angular_drift_minus_external;o.internal_old_torque=r.internal_old_torque;return o;
}

inline canonical_archive::Bytes boundedPhysicalStreamHeader(const canonical_archive::Digest &trusted_contract)
{
    canonical_archive::Bytes bytes{'P','E','R','I','P','H','Y','1'};canonical_archive::appendU64(bytes,1);
    bytes.insert(bytes.end(),trusted_contract.begin(),trusted_contract.end());return bytes;
}

namespace bounded_mechanical_detail
{
using namespace bounded_session_detail;
inline void requireBytes(const canonical_archive::Bytes &expected,const canonical_archive::Bytes &actual,const char *reason)
{
    if(expected==actual)return;size_t offset=0;
    while(offset<std::min(expected.size(),actual.size()) && expected[offset]==actual[offset])++offset;
    throw std::invalid_argument(std::string(reason)+" at canonical byte "+std::to_string(offset)+
        ", expected/actual sizes "+std::to_string(expected.size())+"/"+std::to_string(actual.size()));
}
inline void validateStatic(const BoundedArchiveBody &body,const BoundedSessionContract &contract)
{require(staticBodyDigest(body)==digestArray(contract.static_body_sha256),"mechanical stream: static body identity");}

inline void validateInitial(const BoundedMechanicalRecord &record)
{
    const auto &s=record.body.accepted.gel;const auto &reference=record.body.model.reference;
    require(record.kind==0 && record.interval_phase==0 && sameBits(record.old_time,0) && s.steps==0 &&
        sameBits(s.time,0) && sameBits(s.last_dt,0),"mechanical stream: initial clock");
    for(size_t i=0;i<reference.size();++i)
    {
        for(size_t a=0;a<3;++a)require(sameBits(s.position[i][a],reference.position[i][a]) && s.velocity[i][a]==0,
            "mechanical stream: initial kinematics");
        require(s.viscous[i].squaredNorm()==0 && s.dissipation[i]==0,"mechanical stream: initial material history");
    }
    for(const auto *ledger:{&record.body.accepted.ledger_sum,&record.body.accepted.ledger_correction})
        for(Real x:*ledger)require(x==0,"mechanical stream: initial cumulative ledger");
    BoundedMechanicalRecord expected;expected.body=record.body;
    requireBytes(encodeBoundedMechanicalRecord(expected),encodeBoundedMechanicalRecord(record),
        "mechanical stream: initial physical observation");
}

/** Re-evaluate ONE saved accepted endpoint in scratch; no solveCoulombContact
 * call, reference construction, reset or update of either stored body.
 */
inline void validateAdjacent(const BoundedMechanicalRecord &previous,const BoundedMechanicalRecord &record,
    const BoundedSessionContract &contract)
{
    const auto &old=previous.body;const auto &state=record.body.accepted.gel;
    const auto interval=contract.protocol.interval(old.accepted.gel.steps,old.accepted.gel.time);
    require(record.kind==1 && record.interval_phase==interval.phase && sameBits(record.old_time,interval.old_time) &&
        state.steps==old.accepted.gel.steps+1 && sameBits(state.time,interval.new_time) && sameBits(state.last_dt,interval.actual_dt),
        "mechanical stream: adjacent exact clock");
    const BoundedStoredSceneView scene(old);
    const auto rows=boundedGelFrictionContactsAtOffsets(scene,interval.old_offset,interval.new_offset,old.mu,interval.wall_tangent_velocity);
    const BoundedGelStepContext context(old.model,old.accepted.gel,interval.actual_dt);
    const auto velocity=flattenPD(state.velocity);
    require(context.inDomain(velocity),"mechanical stream: original material path domain");
    const CoulombContactSystem system(flattenPD(old.accepted.gel.position),interval.actual_dt,context.total_mass,rows,contract.solver);
    StdVec<Real> impulse(3*rows.size());
    for(size_t k=0;k<rows.size();++k)
    {impulse[3*k]=state.last_normal_impulse[k];impulse[3*k+1]=record.body.accepted.last_tangent_impulse[k][0];impulse[3*k+2]=record.body.accepted.last_tangent_impulse[k][1];}
    const auto original=system.evaluate(velocity,impulse,[&](const auto &v){return context.evaluate(v).equation;},true);
    const Real KKT=std::max(original.merit,original.KKT);
    require(original.precision_budget_ok && KKT<=contract.solver.relative_tolerance,
        "mechanical stream: original nonsymmetric equation/KKT or precision gate");
    require(sameBits(state.last_KKT,KKT),"mechanical stream: original KKT cursor");
    CoulombContactResult final;final.converged=true;final.velocity=velocity;final.impulses=std::move(impulse);final.normalized_KKT=KKT;
    final.trace.resize(1);final.trace.back().iteration=state.last_newton_iterations;
    const auto proposal=bounded_friction_detail::finishProposal(old.model,old.accepted,interval.actual_dt,rows,old.plate_ids,
        context,system,std::move(final));
    require(proposal.accepted,"mechanical stream: original physical proposal gates");
    auto expected_body=old;expected_body.accepted=proposal.next;
    requireBytes(encodeBoundedArchiveBody(expected_body,1),encodeBoundedArchiveBody(record.body,1),
        "mechanical stream: accepted history byte mismatch");
    const auto whole=boundedFrictionIntervalGaps(scene,old.accepted.gel.position,proposal.next.gel.position,interval);
    const auto endpoint=boundedFrictionIntervalGaps(scene,old.accepted.gel.position,proposal.next.gel.position,interval,false);
    require(std::min(whole[0],whole[1])>=-contract.penetration_stroke_gate*contract.protocol.normal_stroke,
        "mechanical stream: whole-cubic penetration gate");
    auto expected=record;expected.observation=boundedMechanicalObservation(proposal,interval.actual_dt,whole,endpoint);
    requireBytes(encodeBoundedMechanicalRecord(expected),encodeBoundedMechanicalRecord(record),
        "mechanical stream: physical observation byte mismatch");
}
}

/** Proof of a selected committed physical prefix, bound to its exact candidate
 * archive bytes. This proof still has no filesystem/live-session commit API.
 */
class VerifiedBoundedMechanicalPrefix
{
    canonical_archive::VerifiedOutputPrefix prefix_;
    canonical_archive::Digest candidate_digest_;
    static canonical_archive::Digest candidateDigest(const BoundedSessionCandidate &candidate,uint64_t generation)
    {
        const auto bytes=encodeBoundedSessionCandidate(candidate,generation);return canonical_archive::sha256(bytes,bytes.size());
    }
  public:
    VerifiedBoundedMechanicalPrefix(const std::string &path,const BoundedSessionCandidate &candidate,uint64_t generation,
        const canonical_archive::Digest &external_trusted_contract,bool exact_snapshot)
        :prefix_(path,candidate.output.bytes,bounded_session_detail::digestArray(candidate.output.prefix_sha256),exact_snapshot),
         candidate_digest_(candidateDigest(candidate,generation))
    {
        using namespace bounded_session_detail;using namespace canonical_archive;
        candidate.validateCandidate(external_trusted_contract);
        const auto header=boundedPhysicalStreamHeader(external_trusted_contract);Bytes stored(header.size());prefix_.readAt(0,stored.data(),stored.size());
        require(stored==header && sha256(header,header.size())==digestArray(candidate.output.header_sha256),
            "mechanical stream: trusted physical header");
        uint64_t offset=header.size();std::optional<BoundedMechanicalRecord> previous;
        for(uint64_t step=0;step<candidate.output.records;++step)
        {
            std::array<uint8_t,8> length_bytes{};prefix_.readAt(offset,length_bytes.data(),length_bytes.size());offset+=8;
            uint64_t length=0;for(size_t k=0;k<8;++k)length|=uint64_t(length_bytes[k])<<(8*k);
            require(length>=72 && length<=maximum_bytes && offset<=prefix_.bytes() && length<=prefix_.bytes()-offset,
                "mechanical stream: record length");
            Bytes bytes(size_t(length),0);prefix_.readAt(offset,bytes.data(),bytes.size());offset+=length;
            try
            {
                auto record=decodeBoundedMechanicalRecord(bytes,step);
                require(record.body.accepted.gel.steps==step,"mechanical stream: record accepted index");
                bounded_mechanical_detail::validateStatic(record.body,candidate.contract);
                if(previous)bounded_mechanical_detail::validateAdjacent(*previous,record,candidate.contract);
                else bounded_mechanical_detail::validateInitial(record);
                previous=std::move(record);
            }
            catch(const std::exception &error)
            {throw std::invalid_argument("mechanical stream: saved record step "+std::to_string(step)+": "+error.what());}
        }
        require(offset==prefix_.bytes() && previous.has_value(),"mechanical stream: exact committed byte cursor");
        require(encodeBoundedArchiveBody(previous->body,1)==encodeBoundedArchiveBody(candidate.body,1),
            "mechanical stream: final archive/physical output body mismatch");
    }
    void requireBoundCandidate(const BoundedSessionCandidate &candidate,uint64_t generation)const
    {
        bounded_session_detail::require(candidateDigest(candidate,generation)==candidate_digest_,"mechanical stream: candidate changed after verification");
        prefix_.assertUnchanged();
    }
    void stageExclusiveCopy(const BoundedSessionCandidate &candidate,uint64_t generation,int private_directory,const std::string &leaf)const
    {
        requireBoundCandidate(candidate,generation);
        prefix_.stageExclusiveCopy(private_directory,leaf);
    }
};
}
#endif
