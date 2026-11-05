/* Licensed under the Apache License, Version 2.0. */
#ifndef PD_BOUNDED_SESSION_ARCHIVE_H
#define PD_BOUNDED_SESSION_ARCHIVE_H
#include "pd_bounded_archive_body.h"
#include "pd_bounded_physical_protocol.h"
#include "pd_output_prefix.h"
#include <cfenv>
#if defined(__x86_64__)
#include <fpu_control.h>
#include <xmmintrin.h>
#endif

namespace SPH::peridynamics
{
namespace bounded_session_detail
{
using namespace canonical_archive;
inline void require(bool condition,const char *reason)
{if(!condition)throw std::invalid_argument(reason);}
inline bool sameBits(Real a,Real b){return doubleBits(a)==doubleBits(b);}
inline void digestField(const StdVec<uint8_t> &value)
{
    require(value.size()==32 && std::any_of(value.begin(),value.end(),[](uint8_t x){return x!=0;}),
            "bounded session: digest profile");
}
inline Digest digestArray(const StdVec<uint8_t> &value)
{digestField(value);Digest result;std::copy(value.begin(),value.end(),result.begin());return result;}
inline void currentFloatingPointProfile()
{
#if defined(__FAST_MATH__) || (defined(__FINITE_MATH_ONLY__) && __FINITE_MATH_ONLY__)
    throw std::invalid_argument("bounded session: unsupported fast floating point mode");
#endif
#if defined(__x86_64__)
    fpu_control_t control;_FPU_GETCW(control);
    require(std::numeric_limits<Real>::digits==53 && std::numeric_limits<long double>::digits==64 &&
        sizeof(long double)==16 && sizeof(size_t)==8 && std::fegetround()==FE_TONEAREST &&
        (control&0x0300)==0x0300 && (control&0x003f)==0x003f &&
        (_mm_getcsr()&0xffc0)==0x1f80,"bounded session: unsupported current floating point profile");
#else
    throw std::invalid_argument("bounded session: unqualified architecture profile");
#endif
}

// Only the stored immutable model/contact fields; no clock, state, PID or path.
// Generation1 is the fingerprint encoding version, not a selected CP cursor.
inline Bytes staticBodyBytes(BoundedArchiveBody body)
{
    Builder builder(1);FieldWriter writer(builder);
    auto model=writer.scope("model");body.model.registerCheckpoint(model);
    auto surface=writer.scope("surface");body.surface.registerCheckpoint(surface);
    auto dense=writer.scope("dense_surface");body.dense_surface.registerCheckpoint(dense);
    auto contact=writer.scope("contact");
    contact.addScalar("Frame",body.frame);contact.addScalar("LowerNormal",body.normal[0]);
    contact.addScalar("UpperNormal",body.normal[1]);contact.addScalar("Tangent0",body.tangent[0]);
    contact.addScalar("Tangent1",body.tangent[1]);contact.addScalar("LowerMu",body.mu[0]);
    contact.addScalar("UpperMu",body.mu[1]);contact.addArray("PlateIDs",body.plate_ids);
    return builder.finish();
}
inline Digest staticBodyDigest(const BoundedArchiveBody &body)
{const auto bytes=staticBodyBytes(body);return sha256(bytes,bytes.size());}
}

/** The external frozen run plan supplies a trusted fingerprint of this entire
 * contract. Taking a fingerprint FROM an untrusted archive is not admission.
 * Version1 is the current33-node, mode1, fixed16/8/16 friction protocol only.
 */
struct BoundedSessionContract
{
    uint64_t schema=1,equation_version=1,material_path_version=1,transfer_version=1,
        compensation_version=1,floating_point_profile=1,physical_stream_schema=1;
    uint64_t mode=1,normal_base_steps=16,rotation=0,azimuths=8;
    Real cap_width=1./64,penetration_stroke_gate=.01;
    BoundedFrictionProtocol protocol; // explicit schema1 only
    BoundedPhysicalTimeProtocol physical_protocol; // explicit schema2 only
    ManyContactSettings solver;
    StdVec<uint64_t> checkpoint_steps; // schema2 explicit publication cursors
    StdVec<uint8_t> source_manifest_sha256,ELF_sha256,build_runtime_manifest_sha256,
        input_sha256,static_body_sha256;

    template<class Registry>void registerCheckpoint(Registry &cp)
    {
        cp.addScalar("Schema",schema);cp.addScalar("OriginalEquationVersion",equation_version);
        cp.addScalar("AffineFMaterialPathVersion",material_path_version);cp.addScalar("SignedTransferVersion",transfer_version);
        cp.addScalar("CompensationVersion",compensation_version);cp.addScalar("FloatingPointProfile",floating_point_profile);
        cp.addScalar("PhysicalStreamSchema",physical_stream_schema);cp.addScalar("SceneMode",mode);
        if(schema==1)cp.addScalar("NormalBaseSteps",normal_base_steps);
        cp.addScalar("Rotation",rotation);cp.addScalar("Azimuths",azimuths);
        cp.addScalar("CapWidth",cap_width);cp.addScalar("PenetrationStrokeGate",penetration_stroke_gate);
        if(schema==1)
        {
            cp.addScalar("ProtocolKind",protocol.kind);cp.addScalar("NominalDt",protocol.nominal_dt);
            cp.addScalar("NormalStroke",protocol.normal_stroke);cp.addScalar("TwoTangentStrokes",protocol.tangent_stroke);
        }
        else
        {
            bounded_session_detail::require(schema==2,"bounded session: unsupported explicit schema");
            auto &p=physical_protocol;cp.addArray("CheckpointSteps",checkpoint_steps);cp.addScalar("TimeDenominator",p.physical.denominator);
            cp.addScalar("NominalTicks",p.nominal_ticks);
            for(size_t k=0;k<3;++k)
            {
                cp.addScalar("NormalEvent"+std::to_string(k),p.physical.normal_ends[k]);
                cp.addScalar("TangentEvent"+std::to_string(k),p.physical.tangent_ends[k]);
            }
            cp.addScalar("PhysicalNormalStroke",p.physical.normal_stroke);
            for(size_t k=0;k<2;++k)cp.addScalar("PhysicalTangentStroke"+std::to_string(k),p.physical.tangent_stroke[k]);
        }
        cp.addScalar("SolverVelocityScale",solver.velocity_scale);cp.addScalar("SolverRelativeTolerance",solver.relative_tolerance);
        cp.addScalar("MaxNewton",solver.max_newton_steps);cp.addScalar("MaxLineSearch",solver.max_line_search);
        cp.addArray("SourceManifestSHA256",source_manifest_sha256);cp.addArray("ELFSHA256",ELF_sha256);
        cp.addArray("BuildRuntimeManifestSHA256",build_runtime_manifest_sha256);cp.addArray("InputSHA256",input_sha256);
        cp.addArray("StaticBodySHA256",static_body_sha256);
    }
    uint64_t steps()const{return schema==1?protocol.steps():physical_protocol.steps();}
    Real normalStroke()const{return schema==1?protocol.normal_stroke:physical_protocol.physical.normal_stroke;}
    uint64_t phaseAtStep(uint64_t step)const
    {return schema==1?(step==protocol.steps()?3:(step<16?0:(step<24?1:2))):physical_protocol.phaseAtStep(step);}
    BoundedSessionInterval interval(uint64_t step,Real time)const
    {
        if(schema==2)return physical_protocol.interval(step,time);
        BoundedSessionInterval r;static_cast<BoundedFrictionInterval&>(r)=protocol.interval(step,time);return r;
    }
    void validateProfile()const
    {
        using namespace bounded_session_detail;
        currentFloatingPointProfile();
        require(schema==1 || schema==2,"bounded session: unsupported explicit schema");
        if(schema==1)protocol.validate();else
        {
            physical_protocol.validate();
            require(!checkpoint_steps.empty() && checkpoint_steps.size()<=64 && checkpoint_steps.back()==physical_protocol.steps(),
                "bounded session: explicit checkpoint cursor profile");
            uint64_t previous=0;
            for(uint64_t step:checkpoint_steps)
            {require(step>previous && step<=physical_protocol.steps(),"bounded session: checkpoint cursor order");previous=step;}
        }
        require(equation_version==1 && material_path_version==1 && transfer_version==1 && compensation_version==1 &&
            floating_point_profile==1 && physical_stream_schema==schema && mode==1 && (schema==2 || normal_base_steps==16),
            "bounded session: incompatible algorithm version");
        require(rotation<=1 && azimuths>=8 && azimuths<=16 && std::isfinite(cap_width) && cap_width>0 && cap_width<.5,
            "bounded session: scene profile");
        require(sameBits(penetration_stroke_gate,.01) && sameBits(solver.velocity_scale,1) &&
            sameBits(solver.relative_tolerance,1e-12) && solver.max_newton_steps==48 && solver.max_line_search==32,
            "bounded session: incompatible frozen solver or geometry gate");
        for(const auto *value:{&source_manifest_sha256,&ELF_sha256,&build_runtime_manifest_sha256,&input_sha256,&static_body_sha256})
            digestField(*value);
    }
};

inline canonical_archive::Bytes boundedContractBytes(BoundedSessionContract contract)
{
    canonical_archive::Builder builder(1);canonical_archive::FieldWriter fields(builder);
    contract.registerCheckpoint(fields);return builder.finish();
}
inline canonical_archive::Digest boundedContractDigest(const BoundedSessionContract &contract)
{const auto bytes=boundedContractBytes(contract);return canonical_archive::sha256(bytes,bytes.size());}

/** Explicit cursor, not a path or automatic latest-generation lookup.
 * This structural profile does NOT verify the actual stream or its semantics.
 * A loader must check the whole selected prefix before private staging/commit.
 */
struct BoundedPhysicalOutputCursor
{
    uint64_t schema=1,records=0,bytes=0,last_step=0;
    Real last_time=0;
    StdVec<uint8_t> header_sha256,prefix_sha256;
    template<class Registry>void registerCheckpoint(Registry &cp)
    {
        cp.addScalar("Schema",schema);cp.addScalar("Records",records);cp.addScalar("Bytes",bytes);
        cp.addScalar("LastStep",last_step);cp.addScalar("LastTime",last_time);
        cp.addArray("HeaderSHA256",header_sha256);cp.addArray("PrefixSHA256",prefix_sha256);
    }
};

/** Detached archive candidate only. Deliberately no resume/commit method and no
 * constructor which takes a live owner. A structurally valid cursor is not a
 * proof of the physical output prefix, transaction or restart qualification.
 */
struct BoundedSessionCandidate
{
    BoundedArchiveBody body;
    BoundedSessionContract contract;
    BoundedPhysicalOutputCursor output;
    uint64_t accepted_boundary=1,next_phase=0,accepted_tick=0;

    template<class Registry>void registerCheckpoint(Registry &cp)
    {
        auto b=cp.scope("body");body.registerCheckpoint(b);
        auto c=cp.scope("contract");contract.registerCheckpoint(c);
        auto o=cp.scope("output");output.registerCheckpoint(o);
        auto algorithm=cp.scope("accepted_algorithm");algorithm.addScalar("NoInflightTrial",accepted_boundary);
        algorithm.addScalar("NextPhase",next_phase);
        if(contract.schema==2)algorithm.addScalar("AcceptedTick",accepted_tick);
    }
    void validateCandidate(const canonical_archive::Digest &external_trusted_contract)const
    {
        using namespace bounded_session_detail;
        contract.validateProfile();
        require(boundedContractDigest(contract)==external_trusted_contract,"bounded session: external trusted contract mismatch");
        body.validateStored();
        require(staticBodyDigest(body)==digestArray(contract.static_body_sha256),"bounded session: static body fingerprint mismatch");
        Matd expected_frame=Matd::Identity();if(contract.rotation)expected_frame<<0,-1,0,1,0,0,0,0,1;
        for(size_t i=0;i<3;++i)for(size_t j=0;j<3;++j)
            require(sameBits(body.frame(i,j),expected_frame(i,j)),"bounded session: rotation/frame mismatch");
        require(body.surface.reference.size()==7*contract.azimuths+2,"bounded session: azimuth count mismatch");
        const auto &state=body.accepted.gel;
        require(accepted_boundary==1 && state.steps<=contract.steps(),"bounded session: accepted-only boundary");
        Real time=0,last_dt=0;
        if(contract.schema==1)
        {
            // Original v1 binary clock arithmetic and byte schema unchanged.
            for(size_t k=0;k<state.steps;++k)
            {
                const auto interval=contract.protocol.interval(k,time);last_dt=interval.actual_dt;time=interval.new_time;
            }
        }
        else
        {
            const auto &p=contract.physical_protocol;
            require(accepted_tick==p.tickAtStep(state.steps),"bounded session: accepted rational tick/history mismatch");
            time=p.timeAtStep(state.steps);
            if(state.steps)last_dt=p.interval(state.steps-1,p.timeAtStep(state.steps-1)).actual_dt;
        }
        const uint64_t phase=contract.phaseAtStep(state.steps);
        require(next_phase==phase && sameBits(state.time,time) && sameBits(state.last_dt,last_dt),
            "bounded session: exact accepted protocol clock");
        require(output.schema==1 && output.records==state.steps+1 && output.last_step==state.steps &&
            sameBits(output.last_time,time) && output.bytes>=48 && output.bytes<=canonical_archive::maximum_output_prefix,
            "bounded session: structural output cursor");
        digestField(output.header_sha256);digestField(output.prefix_sha256);
    }
};

inline canonical_archive::Bytes encodeBoundedSessionCandidate(BoundedSessionCandidate candidate,uint64_t explicit_generation)
{
    canonical_archive::Builder builder(explicit_generation);canonical_archive::FieldWriter fields(builder);
    candidate.registerCheckpoint(fields);return builder.finish();
}
inline BoundedSessionCandidate decodeBoundedSessionCandidate(const canonical_archive::Bytes &bytes,
    uint64_t explicit_generation,const canonical_archive::Digest &external_trusted_contract)
{
    const canonical_archive::Candidate archive(bytes,explicit_generation);
    return canonical_archive::decodeDetached<BoundedSessionCandidate>(archive,
        [&](const auto &candidate){candidate.validateCandidate(external_trusted_contract);});
}
}
#endif
