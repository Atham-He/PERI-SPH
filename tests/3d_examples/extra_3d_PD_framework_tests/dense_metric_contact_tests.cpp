#include "pd_dense_metric_contact.h"
#include <gtest/gtest.h>
using namespace SPH;
using namespace SPH::peridynamics;

TEST(DenseMetric, OffDiagonalSolveAndNonSPDRejection)
{
    const auto r = DensePositiveMetric({{2, 1}, {1, 2}}).solve({1, 0});
    EXPECT_NEAR(r.value[0], 2.0/3, 1e-15); EXPECT_NEAR(r.value[1], -1.0/3, 1e-15);
    EXPECT_LT(r.relative_residual, 1e-15);
    EXPECT_THROW(DensePositiveMetric({{1, 2}, {2, 1}}), std::invalid_argument);
    EXPECT_THROW(DensePositiveMetric({{1, 1}, {1, 1}}), std::invalid_argument);
    EXPECT_THROW(DensePositiveMetric({{.1, .1}, {.1, .1}}), std::invalid_argument);
    EXPECT_THROW(DensePositiveMetric({{2, 1}, {0, 2}}), std::invalid_argument);
}

TEST(DenseMetric, ConstraintPropagatesThroughFullMetric)
{
    const auto r = projectDenseMetricContact({{2, 1}, {1, 2}}, {0, 0}, {0, 0}, 1, {{{1, 0}, 0, 1}});
    ASSERT_TRUE(r.converged);
    EXPECT_NEAR(r.velocity[0], 1, 1e-14); EXPECT_NEAR(r.velocity[1], -0.5, 1e-14);
    EXPECT_NEAR(r.normal_impulse[0], 1.5, 1e-14);
    // A diagonal substitution would give v2=0 and impulse=2; it cannot pass.
}

TEST(DenseMetric, SignedGatherHasIdenticalTransposeScatter)
{
    const DenseContactMatrix metric{{2, 0.5}, {0.5, 1}};
    const auto r = projectDenseMetricContact(metric, {0, 0}, {0, 0}, 0.1, {{{1.5, -0.5}, 0, 0.02}});
    ASSERT_TRUE(r.converged);
    EXPECT_NEAR(1.5*r.velocity[0]-0.5*r.velocity[1], 0.2, 1e-14);
    EXPECT_NEAR(2*r.velocity[0]+.5*r.velocity[1], 1.5*r.normal_impulse[0], 1e-14);
    EXPECT_NEAR(.5*r.velocity[0]+r.velocity[1], -.5*r.normal_impulse[0], 1e-14);
}

TEST(DenseMetric, InfeasibleRowsCannotConverge)
{
    NormalProjectionSettings settings; settings.max_sweeps=32;
    const auto r = projectDenseMetricContact({{1}}, {0}, {0}, 1, {{{1}, 0, 1}, {{-1}, 0, 1}}, settings);
    EXPECT_FALSE(r.converged); EXPECT_EQ(r.sweeps, 32u);
}

TEST(DenseElastic, ActualMassEnergyAndConvexStiffness)
{
    const auto r = advanceDenseElasticContact({.01, .02}, {{2, -1}, {-1, 2}}, {0, 0}, {0, 0}, .1, {{{1, 0}, 0, .01}});
    ASSERT_TRUE(r.projection.converged);
    EXPECT_NEAR(r.energy_residual, 0, 1e-16);
    EXPECT_NEAR(r.momentum_residual[0], 0, 1e-16); EXPECT_NEAR(r.momentum_residual[1], 0, 1e-16);
    EXPECT_NEAR(r.kinetic_new, (.01*std::pow(r.projection.velocity[0],2)+.02*std::pow(r.projection.velocity[1],2))/2, 1e-16);
    EXPECT_GT(r.algorithmic_loss, 0);
    EXPECT_NO_THROW(validateDenseStiffness({{1, -1}, {-1, 1}}));
    EXPECT_THROW(validateDenseStiffness({{1, 2}, {2, 1}}), std::invalid_argument);
}
