/** Compile-time type coverage only, no executable or restart qualification. */
#include <Eigen/Core>
#include "../src/shared/peridynamics/pd_archive_fields.h"
using namespace SPH::peridynamics::canonical_archive;
namespace
{
struct Components
{
    uint64_t generation=1;
    double time=0;
    std::vector<size_t> offsets;
    std::vector<uint8_t> masks;
    std::vector<double> values;
    std::vector<Eigen::Vector3d> positions;
    std::vector<Eigen::Matrix3d> tensors;
    std::vector<Eigen::Vector2d> tangent;
    template<class Registry>void registerCheckpoint(Registry &registry)
    {
        registry.addScalar("generation",generation);registry.addScalar("time",time);
        registry.addArray("offsets",offsets);registry.addArray("masks",masks);registry.addArray("values",values);
        registry.addArray("positions",positions);registry.addArray("tensors",tensors);registry.addArray("tangent",tangent);
    }
};
}
void instantiateArchiveFieldsTypes()
{
    Components components;Builder builder(7);FieldWriter writer(builder);
    components.registerCheckpoint(writer);
    Candidate archive(builder.finish(),7);
    auto candidate=decodeDetached<Components>(archive,[](const Components&){});
    (void)candidate;
}
