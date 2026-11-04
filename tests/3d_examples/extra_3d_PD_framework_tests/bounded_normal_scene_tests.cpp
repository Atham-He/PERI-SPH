#include "pd_bounded_normal_scene.h"
#include <gtest/gtest.h>
using namespace SPH;using namespace SPH::peridynamics;

TEST(BoundedNormalScene, CompleteMapsAndPlateGeometryUseTheSamePhysicalClock)
{
    BoundedNormalConfig c;c.steps=16;c.dt=.01;const BoundedNormalScene scene(c);
    EXPECT_EQ(scene.surface.reference.size(),114u);EXPECT_EQ(scene.dense_surface.reference.size(),640u);
    EXPECT_EQ(scene.contacts(0,.01).size(),228u);EXPECT_EQ(scene.plateIDs().size(),228u);
    EXPECT_EQ(c.endTime(),Real(.16));EXPECT_EQ(c.upperDisplacement(0),0);EXPECT_EQ(c.upperDisplacement(c.endTime()),c.stroke);
    const auto rows=scene.contacts(.15,.16);
    for(size_t s=0;s<rows.size();++s)
    {
        EXPECT_DOUBLE_EQ(rows[s].old_offset,-1+(s%2?c.upperDisplacement(.15):0));
        EXPECT_DOUBLE_EQ(rows[s].new_offset,-1+(s%2?c.stroke:0));
    }
}

TEST(BoundedNormalScene, EntireIntervalDetectsPenetrationHiddenByExactEndpointContact)
{
    for(size_t count:{size_t(8),size_t(16)})
    {
        BoundedNormalConfig c;c.steps=count;c.dt=.16/count;const BoundedNormalScene scene(c);
        const Real t0=Real(count-1)*c.dt,t1=c.endTime();auto q0=scene.model.reference.position,q1=q0;
        for(size_t i=0;i<q0.size();++i)
        {
            q0[i][1]=(1-.5*c.upperDisplacement(t0))*q0[i][1]-.5*c.upperDisplacement(t0);
            q1[i][1]=(1-.5*c.upperDisplacement(t1))*q1[i][1]-.5*c.upperDisplacement(t1);
        }
        const auto endpoints=scene.denseGaps(q1,t1),interval=scene.denseIntervalGaps(q0,q1,t0,t1);
        EXPECT_GE(endpoints[1],-1e-12);EXPECT_LT(interval[1],0);
        if(count==8)EXPECT_NEAR(-interval[1]/c.stroke,.010259712548221647,1e-11);
        else EXPECT_NEAR(-interval[1]/c.stroke,.002746921031773374,1e-11);
    }
}

TEST(BoundedNormalScene, CoordinateFrameRotationRetainsEverySignedReferenceWeight)
{
    BoundedNormalConfig c;const BoundedNormalScene base(c);c.rotation=1;const BoundedNormalScene rotated(c);
    EXPECT_EQ(base.surface.targets,rotated.surface.targets);ASSERT_EQ(base.surface.weights.size(),rotated.surface.weights.size());
    for(size_t i=0;i<base.surface.weights.size();++i)EXPECT_NEAR(base.surface.weights[i],rotated.surface.weights[i],1e-12);
    EXPECT_LT((rotated.normal[0]-Vecd(-1,0,0)).norm(),1e-15);
    c.initial_affine(0,0)=1.1;
    EXPECT_THROW(c.validate(),std::invalid_argument);
}
