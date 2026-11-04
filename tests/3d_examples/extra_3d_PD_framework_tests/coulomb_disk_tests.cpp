#include "pd_coulomb_disk.h"
#include <gtest/gtest.h>
using namespace SPH;
using namespace SPH::peridynamics;

TEST(CoulombDisk, StickInteriorHasOriginalVelocityEquation)
{
    const auto r = coulombDiskNaturalMap(1, {.1,.05}, {0,0}, .3);
    ASSERT_TRUE(r.within_budget);
    EXPECT_EQ(r.branch, CoulombDiskBranch::StickInterior);
    EXPECT_EQ(r.residual.norm(), 0);
    EXPECT_EQ(r.impulse_derivative.norm(), 0);
    EXPECT_EQ(r.normal_derivative.norm(), 0);
    EXPECT_LT((r.velocity_derivative-ContactTangentMatrix::Identity()).norm(), 1e-15);
}

TEST(CoulombDisk, BothSlipComponentsGiveOpposingImpulseAndPositivePhysicalLoss)
{
    for (Real sign : {Real(-1), Real(1)})
    {
        const ContactTangentVector t(-sign*.24,-sign*.18), w(sign*.4,sign*.3);
        const auto r = coulombDiskNaturalMap(1, t, w, .3);
        ASSERT_TRUE(r.within_budget);
        EXPECT_LT(r.residual.norm(), 1e-15);
        EXPECT_NEAR(r.friction_loss, .15, 1e-15);
        EXPECT_LT(r.cone_excess_upper, 1e-15);
        EXPECT_LT(r.maximum_dissipation_error_upper, 1e-15);
    }
}

TEST(CoulombDisk, OriginalGeneralizedDerivativeIncludesNormalRadiusCoupling)
{
    const StdVec<Real> x{.7,.2,-.35,.5,.1}; const Real epsilon=1e-6;
    const auto r = coulombDiskNaturalMap(x[0], {x[1],x[2]}, {x[3],x[4]}, .4);
    ASSERT_TRUE(r.within_budget);
    EXPECT_GT(r.normal_derivative.norm(), .39);
    for (size_t column=0;column<5;++column)
    {
        auto plus=x,minus=x;plus[column]+=epsilon;minus[column]-=epsilon;
        const auto a=coulombDiskNaturalMap(plus[0],{plus[1],plus[2]},{plus[3],plus[4]},.4);
        const auto b=coulombDiskNaturalMap(minus[0],{minus[1],minus[2]},{minus[3],minus[4]},.4);
        const ContactTangentVector numerical=(a.residual-b.residual)/(2*epsilon);
        const ContactTangentVector expected=column==0 ? r.normal_derivative :
            (column<3 ? ContactTangentVector(r.impulse_derivative.col(column-1)) :
                        ContactTangentVector(r.velocity_derivative.col(column-3)));
        EXPECT_LT((numerical-expected).norm(), 2e-9);
    }
}

TEST(CoulombDisk, TangentialBasisRotationPreservesMapAndWork)
{
    ContactTangentMatrix q; q<<0,-1,1,0;
    const ContactTangentVector t(.2,-.35),w(.5,.1);
    const auto a=coulombDiskNaturalMap(.7,t,w,.4);
    const auto b=coulombDiskNaturalMap(.7,q*t,q*w,.4);
    ASSERT_TRUE(a.within_budget && b.within_budget);
    EXPECT_LT((b.residual-q*a.residual).norm(),1e-15);
    EXPECT_LT((b.impulse_derivative-q*a.impulse_derivative*q.transpose()).norm(),1e-15);
    EXPECT_LT((b.velocity_derivative-q*a.velocity_derivative*q.transpose()).norm(),1e-15);
    EXPECT_LT((b.normal_derivative-q*a.normal_derivative).norm(),1e-15);
    EXPECT_NEAR(a.friction_loss,b.friction_loss,1e-15);
}

TEST(CoulombDisk, ZeroCoefficientAndSeparationRetainBothTangentialCoordinates)
{
    const auto zero=coulombDiskNaturalMap(1,{.2,.1},{.4,.3},0);
    ASSERT_TRUE(zero.within_budget);
    EXPECT_LT((zero.residual-ContactTangentVector(.2,.1)).norm(),1e-15);
    EXPECT_LT((zero.impulse_derivative-ContactTangentMatrix::Identity()).norm(),1e-15);
    const auto separated=coulombDiskNaturalMap(0,{0,0},{.4,.3},.3);
    EXPECT_TRUE(separated.within_budget);
    EXPECT_EQ(separated.residual.norm(),0);
    EXPECT_EQ(separated.friction_loss,0);
}

TEST(CoulombDisk, ExactBoundaryAndOriginUseFiniteLimitingDerivatives)
{
    const auto boundary=coulombDiskNaturalMap(1,{.5,0},{0,0},.5);
    ASSERT_TRUE(boundary.within_budget);
    EXPECT_EQ(boundary.branch,CoulombDiskBranch::SlipExterior);
    EXPECT_EQ(boundary.residual.norm(),0);
    EXPECT_NEAR(boundary.impulse_derivative(0,0),1,1e-15);
    EXPECT_NEAR(boundary.impulse_derivative(1,1),0,1e-15);
    EXPECT_NEAR(boundary.normal_derivative[0],-.5,1e-15);
    const auto origin=coulombDiskNaturalMap(0,{0,0},{0,0},.5);
    EXPECT_TRUE(origin.within_budget);
    EXPECT_EQ(origin.branch,CoulombDiskBranch::ZeroRadiusOrigin);
    EXPECT_LT((origin.impulse_derivative-ContactTangentMatrix::Identity()).norm(),1e-15);
    EXPECT_EQ(origin.velocity_derivative.norm(),0);
}

TEST(CoulombDisk, LargeImpulseSmallRelativeVelocityDoesNotLoseInteriorResidual)
{
    const auto r=coulombDiskNaturalMap(1e8,{1e8,0},{1e-8,0},1);
    ASSERT_TRUE(r.within_budget);
    EXPECT_EQ(r.branch,CoulombDiskBranch::StickInterior);
    EXPECT_DOUBLE_EQ(r.residual[0],1e-8);
    EXPECT_LT(r.residual_roundoff_bound[0],1e-12);
    // Binary64 subtraction followed by projection loses this small residual.
    const Real naive=Real(1e8)-std::min(Real(1e8),Real(1e8)-Real(1e-8));
    EXPECT_GT(std::abs(naive-r.residual[0]),1e-10);
}

TEST(CoulombDisk, SquareBoxCannotReplaceIsotropicDisk)
{
    const auto r=coulombDiskNaturalMap(1,{-.3,-.3},{1,1},.3);
    ASSERT_TRUE(r.within_budget);
    EXPECT_GT(r.residual.norm(),.12);
    EXPECT_GT(r.cone_excess_upper,.12);
}

TEST(CoulombDisk, PrecisionBudgetRejectsUnsupportedRangeWithoutInfCandidate)
{
    const auto large=coulombDiskNaturalMap(1e150,{-.5e150,0},{1e-150,0},.5);
    EXPECT_FALSE(large.within_budget);
    EXPECT_GT(large.residual_roundoff_bound.maxCoeff(),1e-12);
    const Real maximum=std::numeric_limits<Real>::max();
    const auto impossible=coulombDiskNaturalMap(maximum,{maximum,0},{-maximum,0},1);
    EXPECT_FALSE(impossible.within_budget);
    EXPECT_FALSE(impossible.representable);
    EXPECT_TRUE(impossible.residual.allFinite());
    EXPECT_TRUE(std::isfinite(impossible.friction_loss));
}

TEST(CoulombDisk, InvalidCoefficientAndNonfiniteInputsReject)
{
    EXPECT_THROW(coulombDiskNaturalMap(1,{0,0},{0,0},-.1),std::invalid_argument);
    EXPECT_THROW(coulombDiskNaturalMap(std::numeric_limits<Real>::infinity(),{0,0},{0,0},.3),std::invalid_argument);
    EXPECT_THROW(coulombDiskNaturalMap(1,{std::numeric_limits<Real>::quiet_NaN(),0},{0,0},.3),std::invalid_argument);
}
