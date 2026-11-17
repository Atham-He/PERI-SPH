#pragma once
/* Licensed under the Apache License, Version 2.0. */
#include "pd_layered_fracture_reference.h"
#include "pd_cohesive_crack_law.h"
#include "pd_cohesive_interface_imex.h"
#include <algorithm>
#include <array>
#include <cstdint>
#include <cstring>
#include <limits>
#include <set>
#include <string>
#include <vector>

namespace SPH::peridynamics
{
/** Independent persistent physical-crack state. Version 1 has NO conversion
 * from a legacy point-damage field or cached reconstructed/visible stencils.
 * The byte codec can be included as a named attachment in DiskCheckpoint's
 * generation-atomic publication. It deliberately does not create another file
 * publication protocol. FNV-1a64 and a commit footer follow pd_checkpoint.h;
 * this detects accidental corruption and is NOT an adversarial signature.
 * Integers and IEEE binary64 scalars use explicit little-endian encoding.
 *
 * Source identity requires an explicit nonempty compatibility/source tag and
 * ordered, unique original material-point IDs. makeIdentity hashes the EXACT
 * reference coordinates, layered geometry, ordered neighbor connections and
 * reference kernel quadrature; changing point order/geometry is not tolerated.
 * The caller's tag must additionally identify source/model/material inputs
 * outside that geometry. Current deformed positions are owned by the ordinary
 * particle checkpoint and are not substituted for reference geometry here.
 *
 * IMEX initially-rigid kappa and finite-K bilinear (kappa,d,D) are DIFFERENT
 * named typed blocks. Each carries its own unchanged ft/Gc[/K] and ordered,
 * unique surface-quadrature history IDs. Layout must already be known by the
 * receiver; a file cannot silently create/reorder/convert a constitutive law.
 * History IDs belong to the caller's physical surface quadrature, not to a
 * cache row number; include that quadrature definition in source_tag.
 *
 * All bytes, checksum, identity and history values are decoded/validated into
 * temporaries BEFORE mutation. restore then calls ref.restoreTopology(), whose
 * full topology validation and reconstruction are transactional; histories are
 * staged before that call and installed only by noexcept swaps after success.
 * Trace caches are NOT serialized. Existing VisibleCrackSurface owners must
 * ensureCurrent() against the new topologyVersion before their next use.
 * This archive alone is not a full position/velocity/material restart.
 */
class FractureTopologyCheckpoint
{
  public:
    using Reference=solid_dynamics::PDLayeredFractureReference;
    using Bytes=std::vector<uint8_t>;
    struct Identity
    {
        uint64_t dimensions=Dimensions,real_bytes=sizeof(Real),reference_geometry_hash=0;
        std::string source_tag;
        std::vector<uint64_t> source_point_ids;
        std::vector<size_t> neighbor_counts;
    };
    struct IMEXBlock
    {
        std::string name;Real ft=0,Gc=0;
        std::vector<uint64_t> history_ids;
        std::vector<CohesiveInterfaceIMEX::History> history;
    };
    struct FiniteKBlock
    {
        std::string name;Real ft=0,Gc=0,K=0;
        std::vector<uint64_t> history_ids;
        std::vector<BilinearCohesiveCrackLaw::History> history;
    };
    struct Histories
    {
        std::vector<IMEXBlock> imex;
        std::vector<FiniteKBlock> finite_k;
        void swap(Histories&other)noexcept{imex.swap(other.imex);finite_k.swap(other.finite_k);}
    };
    struct State {Identity identity;Reference::TopologyState topology;Histories histories;};
    struct Limits
    {
        uint64_t maximum_bytes=uint64_t(16)<<30;
        uint64_t maximum_facets=10000000;
        uint64_t maximum_history_entries=100000000;
    };
    static constexpr uint64_t schema_version=1;

    static Identity makeIdentity(const Reference&ref,const std::string&source_tag,
                                 const std::vector<uint64_t>&source_point_ids)
    {
        // captureTopology also verifies initialization; pending reconstruction
        // must never be disguised as a committed checkpoint.
        if(ref.topologyPending())throw std::runtime_error("crack checkpoint: pending topology reconstruction");
        if(ref.X.size()!=ref.rows.size()||ref.geometry.reference.size()!=ref.X.size()||source_point_ids.size()!=ref.X.size())
            throw std::invalid_argument("crack checkpoint: source/reference point count mismatch");
        Identity id;id.source_tag=source_tag;id.source_point_ids=source_point_ids;id.neighbor_counts.resize(ref.rows.size());
        Writer geometry(false);geometry.label("NOSB-LAYERED-REFERENCE-1");geometry.u64(ref.X.size());
        for(size_t i=0;i<ref.X.size();++i)
        {
            geometry.vec(ref.X[i]);geometry.vec(ref.geometry.reference[i]);
            const auto&nb=ref.neighbors(i);id.neighbor_counts[i]=ref.rows[i].size();
            if(nb.current_size_!=ref.rows[i].size())throw std::invalid_argument("crack checkpoint: changed reference neighborhood");
            geometry.u64(ref.rows[i].size());
            for(size_t n=0;n<ref.rows[i].size();++n)
            {
                const auto&e=ref.rows[i][n];if(nb.j_.at(n)!=e.j)throw std::invalid_argument("crack checkpoint: neighbor identity mismatch");
                geometry.u64(e.j);geometry.u64(e.reverse);geometry.real(nb.W_ij_.at(n));geometry.real(nb.dW_ij_.at(n));
                geometry.real(nb.r_ij_.at(n));geometry.vec(nb.e_ij_.at(n));
            }
        }
        auto integers=[&](const auto&v){geometry.u64(v.size());for(auto n:v)geometry.u64(n);};
        auto scalars=[&](const auto&v){geometry.u64(v.size());for(Real x:v)geometry.real(x);};
        integers(ref.geometry.column_of);integers(ref.geometry.layer_of);geometry.u64(ref.geometry.columns.size());
        for(const auto&column:ref.geometry.columns)integers(column);
        geometry.u64(ref.geometry.director.size());for(const Vecd&v:ref.geometry.director)geometry.vec(v);
        scalars(ref.geometry.xi);scalars(ref.geometry.thickness);geometry.u64(ref.geometry.derivative.size());
        for(const auto&row:ref.geometry.derivative)scalars(row);
        id.reference_geometry_hash=geometry.hash;validateIdentity(id);return id;
    }
    static State capture(const Reference&ref,const std::string&source_tag,
                         const std::vector<uint64_t>&source_point_ids,const Histories&histories)
    {
        State state;state.topology=ref.captureTopology();state.identity=makeIdentity(ref,source_tag,source_point_ids);
        state.histories=histories;validateState(state);return state;
    }
    static Bytes encode(const Reference&ref,const std::string&source_tag,
                        const std::vector<uint64_t>&source_point_ids,const Histories&histories)
    {return encode(capture(ref,source_tag,source_point_ids,histories));}
    static Bytes encode(const State&state)
    {
        validateState(state);Writer w;w.label("PERI-PHYSICAL-CRACK-CHECKPOINT");w.u64(schema_version);
        w.u64(state.identity.dimensions);w.u64(state.identity.real_bytes);w.u64(64); // scalar wire width
        w.label(state.identity.source_tag);w.u64(state.identity.source_point_ids.size());
        for(auto id:state.identity.source_point_ids)w.u64(id);w.u64(state.identity.reference_geometry_hash);
        w.label("RECIPROCAL-LIVE-MASKS");w.u64(state.topology.live.size());
        for(const auto&row:state.topology.live){w.u64(row.size());for(auto live:row)w.byte(live);}
        w.u64(state.topology.physical_broken_bonds);w.label("REFERENCE-CRACK-FACETS");w.u64(state.topology.facets.size());
        for(const auto&f:state.topology.facets)
        {w.u64(f.kind==Reference::CrackFacet::Kind::PlaneDisk?0:1);w.u64(f.id);w.vec(f.origin);w.vec(f.normal);for(const auto&v:f.vertex)w.vec(v);w.u64(std::isinf(f.radius)?1:0);if(std::isfinite(f.radius))w.real(f.radius);}
        w.label("IMEX-INITIALLY-RIGID-KAPPA-BLOCKS");w.u64(state.histories.imex.size());
        for(const auto&b:state.histories.imex)
        {w.label(b.name);w.real(b.ft);w.real(b.Gc);w.u64(b.history.size());for(size_t i=0;i<b.history.size();++i){w.u64(b.history_ids[i]);w.real(b.history[i].kappa);}}
        w.label("FINITE-K-BILINEAR-KAPPA-D-DISSIPATION-BLOCKS");w.u64(state.histories.finite_k.size());
        for(const auto&b:state.histories.finite_k)
        {w.label(b.name);w.real(b.ft);w.real(b.Gc);w.real(b.K);w.u64(b.history.size());for(size_t i=0;i<b.history.size();++i){w.u64(b.history_ids[i]);w.real(b.history[i].kappa);w.real(b.history[i].damage);w.real(b.history[i].dissipation);}}
        const uint64_t hash=w.hash;w.u64(hash);w.u64(footer);return std::move(w.bytes);
    }
    static State decode(const Bytes&bytes,const Identity&expected,const Histories&expected_layout)
    {return decode(bytes,expected,expected_layout,Limits{});}
    static State decode(const Bytes&bytes,const Identity&expected,const Histories&expected_layout,const Limits&limits)
    {
        validateIdentity(expected);validateHistories(expected_layout);
        if(bytes.size()>limits.maximum_bytes||bytes.size()<16)throw std::runtime_error("crack checkpoint: archive size limit/truncation");
        Reader trailer(bytes,bytes.size()-16,bytes.size());uint64_t stored_hash=trailer.u64();
        if(trailer.u64()!=footer||hashBytes(bytes.data(),bytes.size()-16)!=stored_hash)
            throw std::runtime_error("crack checkpoint: checksum/commit footer mismatch");
        Reader r(bytes,0,bytes.size()-16);r.expectLabel("PERI-PHYSICAL-CRACK-CHECKPOINT");r.expect(schema_version);
        State state;state.identity=expected;r.expect(expected.dimensions);r.expect(expected.real_bytes);r.expect(64);
        r.expectLabel(expected.source_tag);r.expect(expected.source_point_ids.size());for(auto id:expected.source_point_ids)r.expect(id);
        r.expect(expected.reference_geometry_hash);r.expectLabel("RECIPROCAL-LIVE-MASKS");r.expect(expected.source_point_ids.size());
        state.topology.live.resize(expected.neighbor_counts.size());
        for(size_t i=0;i<state.topology.live.size();++i)
        {r.expect(expected.neighbor_counts[i]);r.need(expected.neighbor_counts[i]);auto&row=state.topology.live[i];row.resize(expected.neighbor_counts[i]);for(auto&live:row){live=r.byte();if(live>1)throw std::invalid_argument("crack checkpoint: invalid live mask byte");}}
        state.topology.physical_broken_bonds=r.index();r.expectLabel("REFERENCE-CRACK-FACETS");
        size_t facet_count=r.count(limits.maximum_facets,144);state.topology.facets.resize(facet_count);
        for(auto&f:state.topology.facets)
        {uint64_t kind=r.u64();if(kind>1)throw std::invalid_argument("crack checkpoint: unknown facet kind");f.kind=kind==0?Reference::CrackFacet::Kind::PlaneDisk:Reference::CrackFacet::Kind::Triangle;f.id=r.u64();f.origin=r.vec();f.normal=r.vec();for(auto&v:f.vertex)v=r.vec();uint64_t infinite=r.u64();if(infinite>1)throw std::invalid_argument("crack checkpoint: invalid radius tag");f.radius=infinite?std::numeric_limits<Real>::infinity():r.real();}
        r.expectLabel("IMEX-INITIALLY-RIGID-KAPPA-BLOCKS");r.expect(expected_layout.imex.size());state.histories.imex=expected_layout.imex;
        for(auto&b:state.histories.imex)
        {r.expectLabel(b.name);r.expectReal(b.ft);r.expectReal(b.Gc);r.expect(b.history_ids.size());if(b.history_ids.size()>limits.maximum_history_entries)throw std::invalid_argument("crack checkpoint: history entry limit");r.needProduct(b.history_ids.size(),16);for(size_t i=0;i<b.history.size();++i){r.expect(b.history_ids[i]);b.history[i].kappa=r.real();}}
        r.expectLabel("FINITE-K-BILINEAR-KAPPA-D-DISSIPATION-BLOCKS");r.expect(expected_layout.finite_k.size());state.histories.finite_k=expected_layout.finite_k;
        for(auto&b:state.histories.finite_k)
        {r.expectLabel(b.name);r.expectReal(b.ft);r.expectReal(b.Gc);r.expectReal(b.K);r.expect(b.history_ids.size());if(b.history_ids.size()>limits.maximum_history_entries)throw std::invalid_argument("crack checkpoint: history entry limit");r.needProduct(b.history_ids.size(),32);for(size_t i=0;i<b.history.size();++i){r.expect(b.history_ids[i]);b.history[i].kappa=r.real();b.history[i].damage=r.real();b.history[i].dissipation=r.real();}}
        if(!r.atEnd())throw std::invalid_argument("crack checkpoint: trailing or unrecognized schema data");
        validateState(state);return state;
    }
    /** Strong physical-state transaction: allocations/history validation occur
     * before restoreTopology; that method rolls back masks/facets/stencils if
     * reconstruction fails. No throwing work follows successful restoration. */
    static void restore(Reference&ref,const State&state,const std::string&source_tag,
                        const std::vector<uint64_t>&source_point_ids,Histories&histories)
    {
        Identity current=makeIdentity(ref,source_tag,source_point_ids);expectIdentity(state.identity,current);
        validateState(state);expectLayout(state.histories,histories);Histories staged=state.histories;
        ref.restoreTopology(state.topology);histories.swap(staged);
    }
    static void load(Reference&ref,const Bytes&bytes,const std::string&source_tag,
                     const std::vector<uint64_t>&source_point_ids,Histories&histories)
    {auto current=makeIdentity(ref,source_tag,source_point_ids);auto state=decode(bytes,current,histories);restore(ref,state,source_tag,source_point_ids,histories);}

  private:
    static constexpr uint64_t hash_offset=14695981039346656037ULL,hash_prime=1099511628211ULL;
    static constexpr uint64_t footer=0x5044435241434b31ULL;
    static uint64_t hashBytes(const uint8_t*data,size_t n)
    {uint64_t hash=hash_offset;for(size_t i=0;i<n;++i){hash^=data[i];hash*=hash_prime;}return hash;}
    struct Writer
    {
        Bytes bytes;uint64_t hash=hash_offset;bool recording=true;
        explicit Writer(bool record=true):recording(record){}
        void byte(uint8_t x){if(recording)bytes.push_back(x);hash^=x;hash*=hash_prime;}
        void u64(uint64_t x){for(size_t i=0;i<8;++i)byte(uint8_t(x>>(8*i)));}
        void real(Real x){static_assert(sizeof(double)==8&&std::numeric_limits<double>::is_iec559,"Checkpoint requires IEEE binary64");if(!std::isfinite(x))throw std::invalid_argument("crack checkpoint: nonfinite scalar");double d=x;uint64_t bits;std::memcpy(&bits,&d,8);u64(bits);}
        void vec(const Vecd&x){for(size_t d=0;d<3;++d)real(x[d]);}
        void label(const std::string&s){if(s.empty()||s.size()>4096)throw std::invalid_argument("crack checkpoint: invalid label");u64(s.size());for(unsigned char c:s)byte(c);}
    };
    struct Reader
    {
        const Bytes&bytes;size_t at,end;
        Reader(const Bytes&b,size_t begin,size_t stop):bytes(b),at(begin),end(stop){}
        void need(size_t n)const{if(n>end-at)throw std::runtime_error("crack checkpoint: truncated data");}
        void needProduct(size_t n,size_t width)const{if(n>(end-at)/width)throw std::runtime_error("crack checkpoint: truncated counted data");}
        uint8_t byte(){need(1);return bytes[at++];}
        uint64_t u64(){need(8);uint64_t x=0;for(size_t i=0;i<8;++i)x|=uint64_t(bytes[at++])<<(8*i);return x;}
        size_t index(){auto n=u64();if(n>std::numeric_limits<size_t>::max())throw std::invalid_argument("crack checkpoint: index overflow");return size_t(n);}
        size_t count(uint64_t maximum,size_t width){auto n=u64();if(n>maximum||n>std::numeric_limits<size_t>::max()||n>(end-at)/width)throw std::invalid_argument("crack checkpoint: counted allocation limit");return size_t(n);}
        Real real(){uint64_t bits=u64();double d;std::memcpy(&d,&bits,8);Real x=Real(d);if(!std::isfinite(d)||!std::isfinite(x)||double(x)!=d)throw std::invalid_argument("crack checkpoint: nonfinite/unrepresentable scalar");return x;}
        Vecd vec(){Vecd x;for(size_t d=0;d<3;++d)x[d]=real();return x;}
        void expect(uint64_t value){if(u64()!=value)throw std::invalid_argument("crack checkpoint: incompatible schema/dimension/count/identity");}
        void expectReal(Real value){Real x=real();if(x!=value||std::signbit(x)!=std::signbit(value))throw std::invalid_argument("crack checkpoint: material input mismatch");}
        void expectLabel(const std::string&value){expect(value.size());need(value.size());for(unsigned char c:value)if(byte()!=c)throw std::invalid_argument("crack checkpoint: source/typed block label mismatch");}
        bool atEnd()const{return at==end;}
    };
    static void uniqueIds(const std::vector<uint64_t>&ids)
    {auto sorted=ids;std::sort(sorted.begin(),sorted.end());if(std::adjacent_find(sorted.begin(),sorted.end())!=sorted.end())throw std::invalid_argument("crack checkpoint: duplicate source/history ID");}
    static void label(const std::string&s){if(s.empty()||s.size()>4096)throw std::invalid_argument("crack checkpoint: empty/oversized compatibility or history name");}
    static void validateIdentity(const Identity&id)
    {if(id.dimensions!=3||Dimensions!=3||id.real_bytes!=sizeof(Real)||id.source_point_ids.empty()||id.neighbor_counts.size()!=id.source_point_ids.size())throw std::invalid_argument("crack checkpoint: incompatible identity");label(id.source_tag);uniqueIds(id.source_point_ids);}
    static void expectIdentity(const Identity&a,const Identity&b)
    {if(a.dimensions!=b.dimensions||a.real_bytes!=b.real_bytes||a.source_tag!=b.source_tag||a.source_point_ids!=b.source_point_ids||a.neighbor_counts!=b.neighbor_counts||a.reference_geometry_hash!=b.reference_geometry_hash)throw std::invalid_argument("crack checkpoint: source point/reference geometry identity mismatch");}
    static void validateHistories(const Histories&histories)
    {
        std::set<std::string>names;
        auto common=[&](const auto&b){label(b.name);if(!names.insert(b.name).second||b.history_ids.size()!=b.history.size())throw std::invalid_argument("crack checkpoint: duplicate history name or history count mismatch");uniqueIds(b.history_ids);if(!std::isfinite(b.ft)||!std::isfinite(b.Gc)||b.ft<=0||b.Gc<=0||!std::isfinite(2*b.Gc/b.ft)||2*b.Gc/b.ft<=0)throw std::invalid_argument("crack checkpoint: invalid cohesive ft/Gc");};
        for(const auto&b:histories.imex){common(b);for(const auto&h:b.history)if(!std::isfinite(h.kappa)||h.kappa<0)throw std::invalid_argument("crack checkpoint: invalid IMEX kappa");}
        for(const auto&b:histories.finite_k)
        {
            common(b);BilinearCohesiveCrackLaw law(b.ft,b.Gc,b.K);Real tol=256*std::numeric_limits<Real>::epsilon();
            for(const auto&h:b.history)
            {if(!std::isfinite(h.kappa)||h.kappa<0||!std::isfinite(h.damage)||!std::isfinite(h.dissipation)||h.damage<0||h.damage>1||h.dissipation<0)throw std::invalid_argument("crack checkpoint: invalid finite-K history");auto expected=law.historyForKappa(h.kappa);if(std::abs(h.damage-expected.damage)>tol||std::abs(h.dissipation-expected.dissipation)>tol*std::max(b.Gc,Real(1)))throw std::invalid_argument("crack checkpoint: finite-K history inconsistent with its kappa");}
        }
    }
    static void expectLayout(const Histories&a,const Histories&b)
    {
        validateHistories(b);if(a.imex.size()!=b.imex.size()||a.finite_k.size()!=b.finite_k.size())throw std::invalid_argument("crack checkpoint: cohesive block layout mismatch");
        for(size_t k=0;k<a.imex.size();++k){const auto&x=a.imex[k];const auto&y=b.imex[k];if(x.name!=y.name||x.ft!=y.ft||x.Gc!=y.Gc||x.history_ids!=y.history_ids)throw std::invalid_argument("crack checkpoint: IMEX block identity mismatch");}
        for(size_t k=0;k<a.finite_k.size();++k){const auto&x=a.finite_k[k];const auto&y=b.finite_k[k];if(x.name!=y.name||x.ft!=y.ft||x.Gc!=y.Gc||x.K!=y.K||x.history_ids!=y.history_ids)throw std::invalid_argument("crack checkpoint: finite-K block identity mismatch");}
    }
    static void validateState(const State&state)
    {
        validateIdentity(state.identity);validateHistories(state.histories);
        if(state.topology.live.size()!=state.identity.neighbor_counts.size())throw std::invalid_argument("crack checkpoint: live-mask point count mismatch");
        size_t broken_directed=0;
        for(size_t i=0;i<state.topology.live.size();++i)
        {const auto&row=state.topology.live[i];if(row.size()!=state.identity.neighbor_counts[i])throw std::invalid_argument("crack checkpoint: live-mask neighbor count mismatch");for(auto live:row){if(live>1)throw std::invalid_argument("crack checkpoint: invalid live mask");if(!live){if(broken_directed==std::numeric_limits<size_t>::max())throw std::invalid_argument("crack checkpoint: broken count overflow");++broken_directed;}}}
        if(broken_directed%2||broken_directed/2!=state.topology.physical_broken_bonds)throw std::invalid_argument("crack checkpoint: physical broken-pair count mismatch");
        std::set<uint64_t>facet_ids;
        for(const auto&f:state.topology.facets)
        {
            if((f.kind!=Reference::CrackFacet::Kind::PlaneDisk&&f.kind!=Reference::CrackFacet::Kind::Triangle)||!f.origin.allFinite()||!f.normal.allFinite()||std::abs(f.normal.stableNorm()-1)>1e-10||!(f.radius>0)||!facet_ids.insert(f.id).second)throw std::invalid_argument("crack checkpoint: invalid/duplicate physical facet");
            for(const Vecd&v:f.vertex)if(!v.allFinite())throw std::invalid_argument("crack checkpoint: nonfinite facet vertex");
            if(f.kind==Reference::CrackFacet::Kind::Triangle)
            {Vecd u=f.vertex[1]-f.vertex[0],v=f.vertex[2]-f.vertex[0],n=u.cross(v);Real scale=std::max(u.squaredNorm(),v.squaredNorm()),length=n.stableNorm();if(!std::isfinite(scale)||!std::isfinite(length)||scale<=0||length<=128*std::numeric_limits<Real>::epsilon()*scale||std::abs(std::abs(f.normal.dot(n/length))-1)>1e-10)throw std::invalid_argument("crack checkpoint: degenerate/inconsistent triangle");}
        }
    }
};
} // namespace SPH::peridynamics
