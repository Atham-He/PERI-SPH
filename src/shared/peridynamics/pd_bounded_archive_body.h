/* Licensed under the Apache License, Version 2.0. */
#ifndef PD_BOUNDED_ARCHIVE_BODY_H
#define PD_BOUNDED_ARCHIVE_BODY_H
#include "pd_archive_fields.h"
#include "pd_bounded_gel_friction.h"

namespace SPH::peridynamics
{
/** Value-owned body of the bounded33-node checkpoint. Default construction
 * leaves reference/material/contact packets empty, never rebuilds a scene.
 * Protocol, algorithm contract and committed output transaction are separate
 * mandatory session data; a body roundtrip alone is not restart qualification.
 */
struct BoundedArchiveBody
{
    BoundedGelModel model;
    ReferenceSurfaceTransfer surface,dense_surface;
    Matd frame=Matd::Identity();
    std::array<Vecd,2> normal{Vecd::Zero(),Vecd::Zero()};
    std::array<Vecd,2> tangent{Vecd::Zero(),Vecd::Zero()};
    std::array<Real,2> mu{0,0};
    StdVec<size_t> plate_ids;
    BoundedGelFrictionState accepted;

    template<class Registry>void registerCheckpoint(Registry &registry)
    {
        auto body=registry.scope("model");model.registerCheckpoint(body);
        auto control=registry.scope("surface");surface.registerCheckpoint(control);
        auto dense=registry.scope("dense_surface");dense_surface.registerCheckpoint(dense);
        auto state=registry.scope("accepted");accepted.registerCheckpoint(state);
        auto contact=registry.scope("contact");
        contact.addScalar("Frame",frame);
        contact.addScalar("LowerNormal",normal[0]);contact.addScalar("UpperNormal",normal[1]);
        contact.addScalar("Tangent0",tangent[0]);contact.addScalar("Tangent1",tangent[1]);
        contact.addScalar("LowerMu",mu[0]);contact.addScalar("UpperMu",mu[1]);contact.addArray("PlateIDs",plate_ids);
    }

    /** Readonly consistency of the stored body, no cache repair or solve.
     * Full session validation must additionally pin the trusted static input,
     * exact protocol/solver contract and complete committed output histories.
     */
    void validateStored()const
    {
        model.validate();const auto &ref=model.reference;const size_t n=ref.size();
        if(n!=33 || !frame.allFinite() || (frame.transpose()*frame-Matd::Identity()).norm()>1e-12 ||
            std::abs(frame.determinant()-1)>1e-12)
            throw std::invalid_argument("bounded archive body: dimensions/frame");
        if(!normal[0].allFinite() || !normal[1].allFinite() || !tangent[0].allFinite() || !tangent[1].allFinite() ||
            (normal[0]-frame*Vecd::UnitY()).norm()>1e-12 || (normal[1]+normal[0]).norm()>1e-12 ||
            (tangent[0]-frame*Vecd::UnitX()).norm()>1e-12 || (tangent[1]-frame*Vecd::UnitZ()).norm()>1e-12)
            throw std::invalid_argument("bounded archive body: stored contact frame");
        for(Real value:mu)if(!std::isfinite(value) || value<0)
            throw std::invalid_argument("bounded archive body: friction coefficient");
        surface.validateStored(ref.position.data(),n);dense_surface.validateStored(ref.position.data(),n);
        const size_t m=2*surface.reference.size();
        if(m>256 || plate_ids.size()!=m)throw std::invalid_argument("bounded archive body: contact dimensions");
        for(size_t k=0;k<m;++k)if(plate_ids[k]!=k%2)
            throw std::invalid_argument("bounded archive body: contact order");
        const auto &s=accepted.gel;
        if(s.position.size()!=n || s.velocity.size()!=n || s.green.size()!=n || s.viscous.size()!=n ||
            s.accepted_stress.size()!=n || s.dissipation.size()!=n || s.density.size()!=n || s.steps>98304 ||
            !std::isfinite(s.time) || s.time<0 || !std::isfinite(s.last_dt) || s.last_dt<0 ||
            !std::isfinite(s.last_KKT) || s.last_KKT<0 || s.last_KKT>1e-12 || s.last_newton_iterations>48)
            throw std::invalid_argument("bounded archive body: accepted dimensions/scalars");
        if((s.steps==0 && (s.time!=0 || s.last_dt!=0 || s.last_KKT!=0 || s.last_newton_iterations!=0)) ||
            (s.steps>0 && (!(s.time>0) || !(s.last_dt>0) || s.last_dt>s.time)))
            throw std::invalid_argument("bounded archive body: basic accepted clock");
        Real mass=0,volume=0;long double physical_loss=0;
        for(size_t i=0;i<n;++i)
        {
            if(!s.position[i].allFinite() || !s.velocity[i].allFinite())
                throw std::invalid_argument("bounded archive body: nonfinite kinematics");
            const Matd f=ref.deformation(i,s.position);
            const auto response=gelResponse(f,model.material,s.history(i));
            if((s.green[i]-gelDeviatoricGreen(f)).norm()>1e-12*(1+s.green[i].norm()) ||
                !std::isfinite(s.density[i]) || std::abs(s.density[i]-model.reference_density/f.determinant())>1e-10*model.reference_density ||
                !s.accepted_stress[i].allFinite() || (s.accepted_stress[i]-response.cauchy).norm()>1e-10*(1+response.cauchy.norm()))
                throw std::invalid_argument("bounded archive body: accepted material cursor");
            mass+=model.node_mass[i];volume+=ref.volume[i];physical_loss+=static_cast<long double>(ref.volume[i])*s.dissipation[i];
        }
        if(s.last_normal_impulse.size()!=m || accepted.last_tangent_impulse.size()!=m ||
            accepted.tangent_slip_sum.size()!=m || accepted.tangent_slip_correction.size()!=m ||
            accepted.ledger_sum.size()!=size_t(ContactLedger::Count) || accepted.ledger_correction.size()!=size_t(ContactLedger::Count))
            throw std::invalid_argument("bounded archive body: history dimensions");
        for(size_t k=0;k<m;++k)
        {
            const Real lambda=s.last_normal_impulse[k];
            if(!std::isfinite(lambda) || !accepted.last_tangent_impulse[k].allFinite() ||
                !accepted.tangent_slip_sum[k].allFinite() || !accepted.tangent_slip_correction[k].allFinite())
                throw std::invalid_argument("bounded archive body: nonfinite contact history");
            if(lambda < -1e-10*mass || accepted.last_tangent_impulse[k].norm()>mu[plate_ids[k]]*std::max(Real(0),lambda)+1e-10*mass)
                throw std::invalid_argument("bounded archive body: stored dual cone");
            if(s.steps==0 && (lambda!=0 || accepted.last_tangent_impulse[k].squaredNorm()!=0 ||
                accepted.tangent_slip_sum[k].squaredNorm()!=0 || accepted.tangent_slip_correction[k].squaredNorm()!=0))
                throw std::invalid_argument("bounded archive body: initial contact history");
        }
        for(size_t k=0;k<accepted.ledger_sum.size();++k)
            if(!std::isfinite(accepted.ledger_sum[k]) || !std::isfinite(accepted.ledger_correction[k]))
                throw std::invalid_argument("bounded archive body: nonfinite cumulative history");
        const Real energy=std::max(mass,volume*std::max(model.material.bulk,model.material.equilibrium_shear));
        if(std::abs(physical_loss-accepted.ledger(ContactLedger::MaxwellLoss))>1e-10*energy)
            throw std::invalid_argument("bounded archive body: Maxwell cumulative history");
        // Friction/algorithmic cumulative values and travel compensations need
        // the committed output prefix, not an invented terminal energy rule.
    }
};

inline canonical_archive::Bytes encodeBoundedArchiveBody(BoundedArchiveBody body,uint64_t generation)
{
    canonical_archive::Builder archive(generation);canonical_archive::FieldWriter fields(archive);
    body.registerCheckpoint(fields);return archive.finish();
}

inline BoundedArchiveBody decodeBoundedArchiveBody(const canonical_archive::Bytes &bytes,uint64_t generation)
{
    const canonical_archive::Candidate archive(bytes,generation);
    return canonical_archive::decodeDetached<BoundedArchiveBody>(archive,
        [](const BoundedArchiveBody &candidate){candidate.validateStored();});
}
}
#endif
