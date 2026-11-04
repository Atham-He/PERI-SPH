#include "pd_visco_friction_step.h"
#include "pd_elastic_normal_step.h"
#include <gtest/gtest.h>
using namespace SPH;
using namespace SPH::peridynamics;

TEST(ScalarMaxwell, HeldRelaxationAndExactPhysicalEnergy)
{
    const ScalarMaxwellParameters p{1, 0.5, 0.3};
    const auto r = advanceScalarMaxwell(0.2, 0.2, 0.2, 0.15, p);
    EXPECT_NEAR(r.extension_new, 0.2 * std::exp(-0.5), 1e-15);
    EXPECT_NEAR(r.physical_dissipation, 0.01 * (1 - std::exp(-1.0)), 1e-15);
    EXPECT_NEAR(r.energy_residual, 0, 1e-15);
}

TEST(ScalarMaxwell, SmallIntervalsDoNotLoseRampDissipation)
{
    const ScalarMaxwellParameters p{0, 2, 1};
    for (Real dt : {1e-10, 1e-7, 1e-4})
    {
        const auto r = advanceScalarMaxwell(0, 0, dt, dt, p);
        EXPECT_GT(r.physical_dissipation, 0);
        EXPECT_NEAR(r.physical_dissipation / (2 * dt * dt * dt / 3), 1, dt);
        EXPECT_NEAR(r.extension_new / dt, 1, dt);
        EXPECT_NEAR(r.energy_residual, 0, 1e-16 * dt * dt);
    }
}

TEST(ScalarMaxwell, LinearPathSplitHasExactHistoryWorkAndDissipation)
{
    const ScalarMaxwellParameters p{1, 0.5, 0.3};
    const auto all = advanceScalarMaxwell(0.1, -0.02, 0.3, 0.4, p);
    const auto first = advanceScalarMaxwell(0.1, -0.02, 0.2, 0.2, p);
    const auto second = advanceScalarMaxwell(0.2, first.extension_new, 0.3, 0.2, p);
    EXPECT_NEAR(all.extension_new, second.extension_new, 1e-15);
    EXPECT_NEAR(all.physical_dissipation, first.physical_dissipation + second.physical_dissipation, 1e-15);
    EXPECT_NEAR(all.average_force * 0.2, (first.average_force + second.average_force) * 0.1, 1e-15);
}

TEST(ScalarCoulomb, StickSlipSeparationAndCommonVelocity)
{
    const auto stick = solveScalarCoulomb(0.01, 0, 0.05, 0.3, 0.01);
    EXPECT_EQ(stick.mode, 1);
    EXPECT_NEAR(stick.velocity, 0.05, 1e-15);
    const auto slip = solveScalarCoulomb(0.01, 0, 0.05, 0.3, 0.001);
    const auto boost = solveScalarCoulomb(0.01, 3.25, 3.3, 0.3, 0.001);
    EXPECT_EQ(slip.mode, 3);
    EXPECT_NEAR(slip.impulse, 0.0003, 1e-15);
    EXPECT_GT(slip.physical_dissipation, 0);
    EXPECT_NEAR(slip.impulse, boost.impulse, 1e-15);
    EXPECT_NEAR(slip.physical_dissipation, boost.physical_dissipation, 1e-15);
    const auto separated = solveScalarCoulomb(1, 0, 0.1, 0.3, 0);
    EXPECT_EQ(separated.mode, 0);
    EXPECT_DOUBLE_EQ(separated.impulse, 0);
}

TEST(ViscoFriction, ActualMassPhysicalLedgerAndZeroBranchElasticLimit)
{
    const ViscoFrictionState old{};
    const auto r = advanceViscoFriction(old, 0.01, 0.01, {1, 0.5, 0.3}, 0.01, 0, 0.001, 0.05, 0.3);
    EXPECT_NEAR(r.energy_residual, 0, 1e-15);
    EXPECT_NEAR(r.momentum_normal, 0, 1e-15);
    EXPECT_NEAR(r.momentum_tangent, 0, 1e-15);
    EXPECT_GT(r.material.physical_dissipation, 0);
    EXPECT_GT(r.algorithmic_loss, 0);
    EXPECT_GT(r.friction_loss, 0);
    const auto zero = advanceViscoFriction(old, 0.01, 0.01, {1, 0, 0.3}, 0.01, 0, 0.001, 0, 0);
    const auto elastic = advanceElasticNormal({0.01}, {1}, {0}, {0}, 0.01, {{{0}, {1}, 1, 0, 0.001}});
    EXPECT_NEAR(zero.state.normal_velocity, elastic.projection.velocity[0], 1e-15);
    EXPECT_NEAR(zero.normal.normal_impulse[0], elastic.projection.normal_impulse[0], 1e-15);
}

TEST(ViscoFriction, InvalidInputsCannotCommit)
{
    EXPECT_THROW(advanceScalarMaxwell(0, 0, 0, 0, {1, 1, 1}), std::invalid_argument);
    EXPECT_THROW(solveScalarCoulomb(1, 0, 0, 0.3, -1), std::invalid_argument);
    EXPECT_THROW(solveScalarCoulomb(1, 0, 0, -0.3, 1), std::invalid_argument);
}
