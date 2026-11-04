#include "pd_gel_affine_path.h"
#include <gtest/gtest.h>
using namespace SPH;
using namespace SPH::peridynamics;
namespace
{
GelParameters parameters{5000, 1000, 500, .3};
Matd mixedEndpoint()
{
    Matd f; f << 1.08, .2, -.03, .01, .93, .04, 0, -.02, 1.02; return f;
}
}

TEST(GelAffinePath, QuadraticGreenHistoryAndWorkBalance)
{
    Matd final = Matd::Identity(); final(0, 1) = .2;
    const GelHistory old{};
    const auto r = advanceGelAffinePath(Matd::Identity(), final, .15, parameters, old);
    EXPECT_GT(r.physical_dissipation, 0);
    EXPECT_LT(std::abs(r.energy_residual) / parameters.bulk, 1e-12);
    const auto linear_green = advanceGel(final, .15, parameters, old);
    EXPECT_GT((linear_green.viscous - r.history.viscous).norm(), 1e-5);
    EXPECT_NEAR(r.history.dissipation, r.physical_dissipation, 1e-14);
}

TEST(GelAffinePath, RelaxationHasAnalyticDissipation)
{
    Matd f = Matd::Identity(); f(0, 0) = 1.1; f(1, 1) = .95;
    GelHistory old; old.green = gelDeviatoricGreen(f);
    const auto r = advanceGelAffinePath(f, f, .15, parameters, old);
    EXPECT_NEAR(r.physical_dissipation, parameters.branch_shear * (1 - std::exp(-1.0)) * old.green.squaredNorm(), 1e-12);
    EXPECT_LT((r.history.green - r.history.viscous - std::exp(-.5) * old.green).norm(), 1e-14);
    EXPECT_DOUBLE_EQ(r.work, 0);
    EXPECT_LT(std::abs(r.energy_residual) / parameters.bulk, 1e-12);
}

TEST(GelAffinePath, ConstantSuperposedRotationPreservesHistoryWorkAndStress)
{
    Matd rotation; rotation << 0, -1, 0, 1, 0, 0, 0, 0, 1;
    const auto base = advanceGelAffinePath(Matd::Identity(), mixedEndpoint(), .15, parameters, {});
    const auto rotated = advanceGelAffinePath(rotation, rotation * mixedEndpoint(), .15, parameters, {});
    EXPECT_LT((rotated.average_pk1 - rotation * base.average_pk1).norm() / parameters.bulk, 1e-13);
    EXPECT_LT((rotated.history.viscous - base.history.viscous).norm(), 1e-14);
    EXPECT_NEAR(rotated.physical_dissipation, base.physical_dissipation, 1e-12);
    EXPECT_NEAR(rotated.work, base.work, 1e-12);
}

TEST(GelAffinePath, AnalyticEndpointTangentMatchesDirectionalDifference)
{
    const Matd final = mixedEndpoint();
    const auto base = advanceGelAffinePath(Matd::Identity(), final, .15, parameters, {});
    for (size_t j = 0; j < 9; ++j)
    {
        Matd direction = Matd::Zero(); direction(j / 3, j % 3) = 1;
        const auto plus = advanceGelAffinePath(Matd::Identity(), final + 1e-6 * direction, .15, parameters, {});
        const auto minus = advanceGelAffinePath(Matd::Identity(), final - 1e-6 * direction, .15, parameters, {});
        const Matd difference = (plus.average_pk1 - minus.average_pk1) / 2e-6;
        EXPECT_LT((difference - base.endpoint_tangent[j]).norm() / parameters.bulk, 1e-8);
    }
    EXPECT_GT(base.tangent_asymmetry, 1e-8); // Do not silently force an SPD material tangent.
    auto elastic = parameters; elastic.branch_shear = 0;
    const auto conservative = advanceGelAffinePath(Matd::Identity(), final, .15, elastic, {});
    EXPECT_LT(conservative.tangent_asymmetry, 1e-13);
}

TEST(GelAffinePath, NoAdvancePreservesCompleteAcceptedHistory)
{
    GelHistory old; old.viscous(0, 0) = .01; old.viscous(1, 1) = -.01; old.dissipation = 3.25;
    const auto r = advanceGelAffinePath(Matd::Identity(), Matd::Identity(), 0, parameters, old);
    EXPECT_EQ(r.history.green, old.green); EXPECT_EQ(r.history.viscous, old.viscous);
    EXPECT_DOUBLE_EQ(r.history.dissipation, old.dissipation);
    EXPECT_FALSE(r.tangent_available); EXPECT_DOUBLE_EQ(r.work, 0); EXPECT_DOUBLE_EQ(r.physical_dissipation, 0);
}

TEST(GelAffinePath, FullPathDeterminantAndHistoryMismatchAreRejected)
{
    Matd final = Matd::Identity(); final(0, 0) = final(1, 1) = -1;
    EXPECT_THROW(advanceGelAffinePath(Matd::Identity(), final, .1, parameters, {}), std::runtime_error);
    GelHistory wrong; wrong.green(0, 0) = .01; wrong.green(1, 1) = -.01;
    EXPECT_THROW(advanceGelAffinePath(Matd::Identity(), Matd::Identity(), .1, parameters, wrong), std::invalid_argument);
}

TEST(GelAffinePath, SmallRelaxationIncrementRetainsNonnegativeLoss)
{
    const auto r = advanceGelAffinePath(Matd::Identity(), mixedEndpoint(), 1e-10, parameters, {});
    EXPECT_GT(r.physical_dissipation, 0);
    EXPECT_LT(std::abs(r.energy_residual) / parameters.bulk, 1e-12);
}
