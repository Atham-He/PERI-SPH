#pragma once
#include "pd_isoparametric_crack_surface.h"

namespace SPH::peridynamics
{
/** Reference geometry and area measure of a contact integration patch. */
struct FractureIntegrationPatch
{
    std::array<Vecd,4> vertices;
    StdVec<Vecd> points;
    StdVec<Real> areas;
};

namespace fracture_partition
{
using Facet=solid_dynamics::PDLayeredFractureReference::CrackFacet;
using Polygon=StdVec<Vecd>;

inline Vecd polygonCenter(const Polygon &poly)
{
    Vecd center=Vecd::Zero();
    for (const auto &p : poly)
        center+=p;
    return center/Real(poly.size());
}

inline void appendDistinct(Polygon &poly, const Vecd &p, Real tolerance)
{
    if (poly.empty() || (poly.back()-p).norm()>tolerance)
        poly.push_back(p);
}

inline bool intersectsChord(const Facet &f, const Vecd &a, const Vecd &b, Real tolerance)
{
    const Vecd direction=b-a;
    const Real length=direction.norm();
    if (length<=2*tolerance)
        return false;
    if (f.kind==Facet::Kind::PlaneDisk)
    {
        if (std::isinf(f.radius))
            return true;
        const Real t=std::clamp(direction.dot(f.origin-a)/direction.squaredNorm(),
                                tolerance/length, 1-tolerance/length);
        return (a+t*direction-f.origin).norm()<f.radius;
    }

    // Clip the seam segment against the actual finite triangle. A supporting
    // plane on its own cannot introduce a partition on a distant face.
    Real lower=tolerance/length, upper=1-lower;
    for (size_t i=0;i<3 && lower<upper;++i)
    {
        const Vecd &origin=f.vertex[i], &end=f.vertex[(i+1)%3];
        const Vecd inward=f.normal.cross(end-origin).normalized();
        const Real value=inward.dot(a-origin), slope=inward.dot(direction);
        if (std::abs(slope)<=tolerance)
        {
            if (value < -tolerance)
                return false;
        }
        else if (slope>0)
            lower=std::max(lower,(-tolerance-value)/slope);
        else
            upper=std::min(upper,(-tolerance-value)/slope);
    }
    return (upper-lower)*length>2*tolerance;
}

inline bool splitPolygon(const Polygon &poly,const Facet &facet,Real tolerance,Polygon &minus,Polygon &plus)
{
    StdVec<Real> signed_distance;
    bool positive=false,negative=false;
    for (const auto &p : poly)
    {
        const Real d=facet.normal.dot(p-facet.origin);
        const Real snapped=std::abs(d)<=tolerance ? 0 : d;
        signed_distance.push_back(snapped);
        positive|=snapped>0;
        negative|=snapped<0;
    }
    if (!positive || !negative)
        return false;

    Polygon seam;
    for (size_t i=0;i<poly.size();++i)
    {
        const auto next=(i+1)%poly.size();
        const Real a=signed_distance[i], b=signed_distance[next];
        if (a<=0) appendDistinct(minus,poly[i],tolerance);
        if (a>=0) appendDistinct(plus,poly[i],tolerance);
        if (a==0) appendDistinct(seam,poly[i],tolerance);
        if ((a<0 && b>0) || (a>0 && b<0))
        {
            const Vecd intersection=poly[i]+(a/(a-b))*(poly[next]-poly[i]);
            appendDistinct(minus,intersection,tolerance);
            appendDistinct(plus,intersection,tolerance);
            appendDistinct(seam,intersection,tolerance);
        }
    }
    for (auto *side : {&minus,&plus})
        if (side->size()>1 && (side->front()-side->back()).norm()<=tolerance)
            side->pop_back();
    if (minus.size()<3 || plus.size()<3 || seam.size()<2)
        return false;
    auto farthest=std::max_element(seam.begin()+1,seam.end(),[&](const Vecd &u,const Vecd &v)
    {return (u-seam.front()).squaredNorm()<(v-seam.front()).squaredNorm();});
    return intersectsChord(facet,seam.front(),*farthest,tolerance);
}

inline void quadrulate(const Polygon &poly, StdVec<FractureIntegrationPatch> &patches,Real tolerance)
{
    auto add=[&](const std::array<Vecd,4> &quad)
    {
        FractureIntegrationPatch child;
        child.vertices=quad;
        for (Real u : {.5-.5/std::sqrt(3.),.5+.5/std::sqrt(3.)})
            for (Real v : {.5-.5/std::sqrt(3.),.5+.5/std::sqrt(3.)})
            {
                const auto p=BilinearCrackGeometry::evaluate(quad,Vec2d(u,v));
                const Real measure=.25*p.u.cross(p.v).norm();
                if (!(measure>tolerance*tolerance))
                    throw std::runtime_error("Split crack quadrature has a degenerate child");
                child.points.push_back(p.x);
                child.areas.push_back(measure);
            }
        patches.push_back(std::move(child));
    };
    if (poly.size()==4)
    {
        add({poly[0],poly[1],poly[2],poly[3]});
        return;
    }

    // The centroid fan partitions a convex polygon; each triangle becomes
    // three bilinear quadrilaterals without adding gaps or overlapping area.
    const Vecd centroid=polygonCenter(poly);
    for (size_t i=0;i<poly.size();++i)
    {
        const Vecd &a=poly[i], &b=poly[(i+1)%poly.size()];
        if ((b-a).cross(centroid-a).norm()<=tolerance*tolerance)
            continue;
        const Vecd ab=.5*(a+b), bc=.5*(b+centroid), ca=.5*(centroid+a), center=(a+b+centroid)/3;
        add({a,ab,center,ca});
        add({b,bc,center,ab});
        add({centroid,ca,center,bc});
    }
}

inline StdVec<FractureIntegrationPatch> split(const FractureIntegrationPatch &parent,const StdVec<Facet> &facets)
{
    Polygon original(parent.vertices.begin(),parent.vertices.end());
    const Vecd origin=original.front();
    const Vecd normal=(original[1]-origin).cross(original[3]-origin);
    Real scale=0;
    for (const auto &p : original)
        scale=std::max(scale,(p-origin).norm());
    const Real tolerance=256*std::numeric_limits<Real>::epsilon()*scale;
    if (!(normal.norm()>tolerance*tolerance))
        throw std::runtime_error("Crack integration patch has no reference area");
    const Vecd unit_normal=normal.normalized();
    std::vector<Polygon> polygons{original};
    for (const auto &facet : facets)
    {
        if (std::abs(facet.normal.dot(unit_normal))>1-64*std::numeric_limits<Real>::epsilon())
            continue;
        std::vector<Polygon> next;
        for (auto &poly : polygons)
        {
            Polygon minus,plus;
            if (splitPolygon(poly,facet,tolerance,minus,plus))
            {
                next.push_back(std::move(minus));
                next.push_back(std::move(plus));
            }
            else
                next.push_back(std::move(poly));
        }
        polygons=std::move(next);
    }
    if (polygons.size()==1)
        return {};
    for (const auto &p : original)
        if (std::abs(unit_normal.dot(p-origin))>tolerance)
            throw std::runtime_error("Crack integration subdivision requires a planar reference patch");

    StdVec<FractureIntegrationPatch> result;
    for (const auto &poly : polygons)
        quadrulate(poly,result,tolerance);
    if (result.empty())
        throw std::runtime_error("Empty crack integration subdivision");

    // Keep the authored area and first moment. In a curved shell bank the
    // supplied areas measure the annular section rather than its chord.
    const Vecd e0=(original[1]-origin).normalized(), e1=unit_normal.cross(e0);
    auto moment=[&](const Vecd &p)
    {
        const Vecd d=(p-origin)/scale;
        return Vecd(1,d.dot(e0),d.dot(e1));
    };
    Vecd target=Vecd::Zero();
    for (size_t i=0;i<parent.points.size();++i)
        target+=parent.areas[i]*moment(parent.points[i]);
    Vecd current=Vecd::Zero();
    Matd gram=Matd::Zero();
    for (const auto &patch : result)
        for (size_t i=0;i<patch.points.size();++i)
        {
            const Vecd m=moment(patch.points[i]);
            current+=patch.areas[i]*m;
            gram+=patch.areas[i]*m*m.transpose();
        }
    const Vecd correction=gram.ldlt().solve(target-current);
    Vecd integrated=Vecd::Zero();
    for (auto &patch : result)
        for (size_t i=0;i<patch.points.size();++i)
        {
            const Vecd m=moment(patch.points[i]);
            patch.areas[i]*=1+m.dot(correction);
            if (!(patch.areas[i]>0) || !std::isfinite(patch.areas[i]))
                throw std::runtime_error("Split crack quadrature cannot retain a positive area measure");
            integrated+=patch.areas[i]*m;
        }
    if ((integrated-target).norm()>1e-12*target[0])
        throw std::runtime_error("Split crack quadrature did not conserve area and first moment");
    return result;
}
}
}
