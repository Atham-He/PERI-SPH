#include "pd_normal_projection.h"
#include <gtest/gtest.h>

using namespace SPH;
using namespace SPH::peridynamics;

namespace
{
NormalContactConstraint lower(size_t i = 0, Real old_face = 0, Real new_face = 0)
{ return {{i}, {1}, 1, old_face, new_face}; }
NormalContactConstraint upper(size_t i = 0, Real old_face = 1, Real new_face = 1)
{ return {{i}, {1}, -1, -old_face, -new_face}; }
void checkBalances(const NormalProjectionResult &r)
{
    ASSERT_TRUE(r.converged);
    EXPECT_LT(std::abs(r.momentum_residual), 1e-11);
    EXPECT_LT(std::abs(r.energy_residual), 1e-11);
    EXPECT_GE(r.contact_step_dissipation, -1e-13);
}
}

TEST(PDNormalProjection, FreeFlightIsUnchanged)
{
    const auto r = projectNormalContact({1, 2}, {1, 1}, {-0.1, 0.2}, 0.1, {});
    checkBalances(r);
    EXPECT_EQ(r.velocity, (StdVec<Real>{-0.1, 0.2}));
    EXPECT_NEAR(r.position[0], 0.99, 1e-15); EXPECT_NEAR(r.position[1], 1.02, 1e-15);
    EXPECT_DOUBLE_EQ(r.contact_step_dissipation, 0);
}

TEST(PDNormalProjection, SingleFixedWallMatchesAnalyticImpulseAndDissipation)
{
    const auto r = projectNormalContact({2}, {0.1}, {-1}, 0.2, {lower()});
    checkBalances(r);
    EXPECT_NEAR(r.velocity[0], -0.5, 1e-14); EXPECT_NEAR(r.position[0], 0, 1e-14);
    EXPECT_NEAR(r.normal_impulse[0], 1, 1e-14);
    EXPECT_NEAR(r.contact_step_dissipation, 0.75, 1e-14);
}

TEST(PDNormalProjection, MovingWallWorkHasExplicitNonnegativeContactStepDissipation)
{
    const auto r = projectNormalContact({2}, {0}, {0}, 0.2, {lower(0, 0, 0.1)});
    checkBalances(r);
    EXPECT_NEAR(r.velocity[0], 0.5, 1e-14); EXPECT_NEAR(r.position[0], 0.1, 1e-14);
    EXPECT_NEAR(r.wall_work, 0.5, 1e-14); EXPECT_NEAR(r.accepted_kinetic, 0.25, 1e-14);
    EXPECT_NEAR(r.contact_step_dissipation, 0.25, 1e-14);
}

TEST(PDNormalProjection, SignedMLSTransferUsesTheSameMassMetricTranspose)
{
    NormalContactConstraint row{{0, 1}, {1.25, -0.25}, 1, 0, 0};
    const auto r = projectNormalContact({2, 1}, {0.1, 0.1}, {-1, -1}, 0.2, {row});
    checkBalances(r);
    const Real impulse = 0.5 / (1.25 * 1.25 / 2 + 0.25 * 0.25);
    EXPECT_NEAR(r.normal_impulse[0], impulse, 1e-13);
    EXPECT_NEAR(1.25 * r.position[0] - 0.25 * r.position[1], 0, 1e-13);
    EXPECT_NEAR(r.velocity[0], -1 + 1.25 * impulse / 2, 1e-13);
    EXPECT_NEAR(r.velocity[1], -1 - 0.25 * impulse, 1e-13);
}

TEST(PDNormalProjection, BothWallImpulsesAreNonzeroWhenNetImpulseIsZero)
{
    const auto r = projectNormalContact({1, 1}, {0, 1}, {-1, 1}, 0.1, {lower(0), upper(1)});
    checkBalances(r);
    EXPECT_NEAR(r.normal_impulse[0], 1, 1e-14); EXPECT_NEAR(r.normal_impulse[1], 1, 1e-14);
    EXPECT_NEAR(r.contact_step_dissipation, 1, 1e-14);
    EXPECT_EQ(r.velocity, (StdVec<Real>{0, 0}));
}

TEST(PDNormalProjection, RedundantRowsPreserveUniqueVelocityWithoutAssumingUniqueDuals)
{
    const auto r = projectNormalContact({1}, {0}, {-2}, 0.1, {lower(), lower()});
    checkBalances(r);
    EXPECT_NEAR(r.velocity[0], 0, 1e-14);
    EXPECT_NEAR(r.normal_impulse[0] + r.normal_impulse[1], 2, 1e-14);
}

TEST(PDNormalProjection, CoupledRowsMatchTheIndependentTwoConstraintSolution)
{
    const StdVec<NormalContactConstraint> rows{{{0, 1}, {0.7, 0.3}, 1, 0, 0}, {{0, 1}, {0.2, 0.8}, 1, 0, 0}};
    const auto r = projectNormalContact({2, 3}, {0, 0}, {-1, -0.5}, 0.1, rows);
    checkBalances(r);
    EXPECT_NEAR(r.velocity[0], 0, 1e-10); EXPECT_NEAR(r.velocity[1], 0, 1e-10);
    EXPECT_NEAR(r.normal_impulse[0], 2.6, 1e-10); EXPECT_NEAR(r.normal_impulse[1], 0.9, 1e-10);
}

TEST(PDNormalProjection, InfeasibleMotionAndInvalidInitialPenetrationCannotPass)
{
    NormalProjectionSettings settings; settings.max_sweeps = 64;
    const auto r = projectNormalContact({1}, {0.5}, {0}, 0.1, {lower(0, 0, 0.7), upper(0, 1, 0.3)}, settings);
    EXPECT_FALSE(r.converged); EXPECT_EQ(r.sweeps, 64); EXPECT_GT(r.primal_velocity_error, 1);
    EXPECT_THROW(projectNormalContact({1}, {-0.1}, {0}, 0.1, {lower()}), std::invalid_argument);
    EXPECT_THROW(projectNormalContact({0}, {0}, {0}, 0.1, {lower()}), std::invalid_argument);
}

TEST(PDNormalProjection, ImpactRefinementKeepsPositionFeasibleAndTotalImpulseEnergyExact)
{
    for (Real candidate : {0.36, 0.18, 0.09, 0.045})
    {
        Real t = 0, q = 1, v = -1, impulse = 0, loss = 0;
        while (t < 1.7)
        {
            const Real dt = std::min(candidate, 1.7 - t);
            const auto r = projectNormalContact({1}, {q}, {v}, dt, {lower()});
            checkBalances(r);
            EXPECT_GE(r.position[0], -1e-13);
            impulse += r.normal_impulse[0]; loss += r.contact_step_dissipation;
            q = r.position[0]; v = r.velocity[0]; t = dt == 1.7 - t ? 1.7 : t + dt;
        }
        EXPECT_NEAR(q, 0, 1e-12); EXPECT_NEAR(v, 0, 1e-12);
        EXPECT_NEAR(impulse, 1, 1e-12); EXPECT_NEAR(loss, 0.5, 1e-12);
    }
}
