#pragma once
#include "pd_layered_geometry.h"
#include "pd_solid_reference_coupling.h"

namespace SPH::peridynamics
{
/** One material trace per thickness column at the chosen physical shell face.
 * The normal polynomial and its transpose retain the surface lever arm. A
 * column does not use other columns as support and a through-thickness crack
 * between columns therefore cannot change this map. */
class ColumnFaceSolidSurface : public ReferenceSolidSurface
{
  public:
    ColumnFaceSolidSurface(SolidForcePort &port, const LayeredParticleGeometry &geometry,
                           Real face_coordinate)
        : ReferenceSolidSurface(port, {}, 2, DeferredBuildTag{})
    {
        const auto weights = geometry.trace(face_coordinate);
        for (const auto &column : geometry.columns)
        {
            Vecd point = Vecd::Zero();
            for (size_t level = 0; level < column.size(); ++level)
                point += weights[level] * port.reference.at(column[level]);
            points.push_back(point);
            rows.push_back({column, weights});
        }
    }
};
}
