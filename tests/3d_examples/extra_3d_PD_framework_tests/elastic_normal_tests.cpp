#include "pd_elastic_normal_step.h"
#include <gtest/gtest.h>
using namespace SPH;
using namespace SPH::peridynamics;

TEST(ElasticNormal, EffectiveMetricIsNotPhysicalMass)
{
    const auto r=advanceElasticNormal({2},{8},{0},{0},0.5,{{{0},{1},1,0,0.1}});
    ASSERT_TRUE(r.projection.converged);
    EXPECT_NEAR(r.projection.velocity[0],0.2,1e-14);
    EXPECT_NEAR(r.projection.normal_impulse[0],0.6,1e-14);
    EXPECT_NEAR(r.kinetic_new,0.04,1e-14);
    EXPECT_NEAR(r.projection.accepted_quadratic,0.06,1e-14);
    EXPECT_NEAR(r.elastic_new,0.04,1e-14);
    EXPECT_NEAR(r.algorithmic_loss,0.04,1e-14);
    EXPECT_NEAR(r.wall_work,0.12,1e-14);
    EXPECT_NEAR(r.energy_residual,0,1e-14);
    EXPECT_NEAR(r.maximum_momentum_residual,0,1e-14);
}

TEST(ElasticNormal, FreeSpringPhysicalBalanceAndLoss)
{
    const auto r=advanceElasticNormal({2},{8},{0.3},{-0.2},0.1,{});
    ASSERT_TRUE(r.projection.converged);
    EXPECT_NEAR(r.maximum_momentum_residual,0,1e-14);
    EXPECT_NEAR(r.energy_residual,0,1e-14);
    EXPECT_GT(r.algorithmic_loss,0);
    EXPECT_DOUBLE_EQ(r.wall_work,0);
}

TEST(ElasticNormal, SeparatingWallDoesNotPull)
{
    const auto r=advanceElasticNormal({1},{0},{0},{0},0.1,{{{0},{1},1,0,-1}});
    EXPECT_TRUE(r.projection.converged);
    EXPECT_DOUBLE_EQ(r.projection.normal_impulse[0],0);
    EXPECT_DOUBLE_EQ(r.projection.position[0],0);
}

TEST(ElasticNormal, ZeroStiffnessPreservesMassContactAPI)
{
    const StdVec<NormalContactConstraint> c={{{0},{1},1,0,0.2}};
    const auto a=advanceElasticNormal({3},{0},{0.1},{-0.4},0.1,c);
    const auto b=projectNormalContact({3},{0.1},{-0.4},0.1,c);
    EXPECT_EQ(a.projection.velocity,b.velocity);
    EXPECT_EQ(a.projection.normal_impulse,b.normal_impulse);
    EXPECT_NEAR(a.algorithmic_loss,b.contact_step_dissipation,1e-14);
}

TEST(ElasticNormal, RejectsInvalidPhysicalInputs)
{
    EXPECT_THROW(advanceElasticNormal({1},{-1},{0},{0},0.1,{}),std::invalid_argument);
    EXPECT_THROW(advanceElasticNormal({0},{1},{0},{0},0.1,{}),std::invalid_argument);
    EXPECT_THROW(advanceElasticNormal({1},{1},{0},{0},0,{}),std::invalid_argument);
}
