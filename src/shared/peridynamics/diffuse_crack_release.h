#pragma once
#include "pd_energy_crack_band.h"
#include <fstream>

namespace SPH::peridynamics
{
/** Transport through a damage-only crack band, not an explicit free surface.
 * Grassl, arXiv:0809.2758 Eq.7: w = band * d * strain (plastic strain is zero).
 * Adaptations: current normal plane spacing from F^-T and layers in SERIES.
 * Open area is aperture/CURRENT band width, not aperture/reference width.
 * Fick conductance of the open area is D*A*open_fraction/current_thickness.
 * A closed/intact layer blocks the crack contribution exactly. Base EC
 * permeability is an independent physical input, not a numerical leak.
 */
struct DiffuseBarrierResult
{
    Real open_fraction=0,transport_length=0,max_equivalent_aperture=0,min_damage=1;
};

inline DiffuseBarrierResult diffuseLayerBarrier(const StdVec<Matd>&F,
    const StdVec<Real>&damage,const StdVec<Vecd>&normal,const StdVec<Real>&weights,
    const Vecd&thickness_normal,Real thickness,Real band_width)
{
    if(F.empty()||F.size()!=damage.size()||F.size()!=normal.size()||F.size()!=weights.size()||
       thickness<=0||band_width<=0||std::abs(thickness_normal.norm()-1)>1e-8)
        throw std::invalid_argument("Invalid layered crack barrier");
    Real sum=0;for(Real w:weights){if(w<=0)throw std::invalid_argument("Invalid layer weights");sum+=w;}
    DiffuseBarrierResult r;Real resistance_length=0;bool blocked=false;
    for(size_t l=0;l<F.size();++l)
    {
        if(!F[l].allFinite()||F[l].determinant()<=0||damage[l]<0||damage[l]>1)
            throw std::runtime_error("Invalid crack barrier material state");
        Matd inverse_transpose=F[l].inverse().transpose();
        // Current normal spacing of reference parallel planes (Nanson).
        Real radial_stretch=1/(inverse_transpose*thickness_normal).norm();
        Real length=thickness*weights[l]/sum*radial_stretch;
        r.transport_length+=length;r.min_damage=std::min(r.min_damage,damage[l]);
        Real fraction=0;
        if(damage[l]>0)
        {
            if(std::abs(normal[l].norm()-1)>1e-7)throw std::runtime_error("Missing crack normal for damaged shell");
            // A reference slab of width h and normal N has current normal
            // spacing h/|F^-T N|. Use that SAME current width in the open-area
            // fraction; reference-width/current-area mixing is inconsistent.
            Real normal_stretch=1/(inverse_transpose*normal[l]).norm();
            Real tensile=std::max(Real(0),normal_stretch-1);
            Real aperture=band_width*damage[l]*tensile;
            r.max_equivalent_aperture=std::max(r.max_equivalent_aperture,aperture);
            // Area averaging of a slot across its represented crack-band width.
            fraction=aperture/(band_width*normal_stretch);
            if(fraction>1)throw std::runtime_error("Crack-band open area exceeds represented band; no clipping");
        }
        if(fraction==0)blocked=true;else resistance_length+=length/fraction;
    }
    if(!blocked)r.open_fraction=r.transport_length/resistance_length;
    return r;
}

template<class Geometry,class Shell>
void writeLayeredReleaseBarrier(const std::string&path,const Geometry&g,const Shell&shell,
    const Matd*F,const StdVec<Vecd>&directions,const StdVec<Real>&weights)
{
    std::ofstream f(path);if(!f)throw std::runtime_error("Cannot write crack barrier");
    f<<std::setprecision(17)<<"column,nx,ny,nz,open_fraction,transport_length_m,equivalent_aperture_m,min_layer_damage\n";
    for(size_t c=0;c<g.columns.size();++c)
    {
        StdVec<Matd> state;StdVec<Real>d;StdVec<Vecd>n;
        for(size_t j:g.columns[c]){state.push_back(F[j]);d.push_back(shell.damage[j]);n.push_back(shell.crack_normal[j]);}
        auto r=diffuseLayerBarrier(state,d,n,weights,directions[c],g.thickness[c],shell.fracture.width);
        f<<c<<','<<directions[c][0]<<','<<directions[c][1]<<','<<directions[c][2]<<','<<r.open_fraction<<','<<r.transport_length<<','<<r.max_equivalent_aperture<<','<<r.min_damage<<'\n';
    }
    if(!f)throw std::runtime_error("Crack barrier output failed");
}
}
