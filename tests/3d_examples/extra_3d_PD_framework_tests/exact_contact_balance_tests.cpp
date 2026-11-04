#include "pd_many_contact.h"
#include <gtest/gtest.h>
#include <cfenv>
using namespace SPH;
using namespace SPH::peridynamics;

TEST(ExactContactBalance, HugeCancellationRetainsSmallMomentumAndProductLowBits)
{
    const Real huge=std::ldexp(Real(1),60);
    const auto r=exactContactBalance(.125,{1,1,1},{huge,.25,-huge},1);
    EXPECT_EQ(r.normalized_value,-.125);EXPECT_EQ(r.normalized_roundoff_bound,0);EXPECT_TRUE(r.within_budget);
    const Real delta=std::ldexp(Real(1),-27);
    const auto product=exactContactBalance(1,{1+delta},{1-delta},1);
    EXPECT_EQ(product.normalized_value,std::ldexp(Real(1),-54));EXPECT_EQ(product.normalized_roundoff_bound,0);
}

TEST(ExactContactBalance, IntermediateProductOverflowIsNotAFalseInputFailure)
{
    const Real largest=std::numeric_limits<Real>::max();
    const auto r=exactContactBalance(.5,{largest,-largest},{2,2},1);
    EXPECT_EQ(r.normalized_value,.5);EXPECT_EQ(r.normalized_roundoff_bound,0);
    EXPECT_THROW(exactContactBalance(largest,{}, {},std::numeric_limits<Real>::denorm_min()),std::overflow_error);
}

TEST(ExactContactBalance, SubnormalAndTiesUseExactIntegerRounding)
{
    const Real tiny=std::numeric_limits<Real>::denorm_min(),half_ulp=std::ldexp(Real(1),-53);
    const auto underflow=exactContactBalance(0,{tiny},{tiny},1);
    EXPECT_EQ(underflow.normalized_value,0);EXPECT_EQ(underflow.normalized_roundoff_bound,tiny);
    const auto even=exactContactBalance(1,{-half_ulp},{1},1);
    EXPECT_EQ(even.normalized_value,1);EXPECT_EQ(even.normalized_roundoff_bound,half_ulp);
    const auto odd=exactContactBalance(1+2*half_ulp,{-half_ulp},{1},1);
    EXPECT_EQ(odd.normalized_value,1+4*half_ulp);EXPECT_EQ(odd.normalized_roundoff_bound,half_ulp);
}

TEST(ExactContactBalance, PrecisionBudgetAndInvalidityAreExplicit)
{
    const auto precise=exactContactBalance(1,{}, {},3), insufficient=exactContactBalance(1e9,{}, {},3);
    EXPECT_EQ(precise.normalized_value,Real(1)/3);EXPECT_TRUE(precise.within_budget);
    EXPECT_FALSE(insufficient.within_budget);EXPECT_GT(insufficient.normalized_roundoff_bound,1e-12);
    EXPECT_THROW(exactContactBalance(1,{}, {},0),std::invalid_argument);
    EXPECT_THROW(exactContactBalance(1,{1}, {},1),std::invalid_argument);
    EXPECT_THROW(exactContactBalance(std::numeric_limits<Real>::infinity(),{}, {},1),std::invalid_argument);
}

TEST(ExactContactBalance, RoundingModeDoesNotChangeExactRatioDecision)
{
    const int previous=std::fegetround();
    const auto nearest=exactContactBalance(1,{}, {},3);
    ASSERT_EQ(std::fesetround(FE_UPWARD),0);
    const auto upward=exactContactBalance(1,{}, {},3);
    EXPECT_EQ(std::fesetround(previous),0);
    EXPECT_EQ(nearest.normalized_value,upward.normalized_value);
    EXPECT_EQ(nearest.normalized_roundoff_bound,upward.normalized_roundoff_bound);
}

TEST(ExactContactBalance, WholeSolverRejectsInitialAndTrialPrecisionLossWithoutCommit)
{
    auto initial=[](const StdVec<Real> &v){return NonlinearContactEquation{{1e9+v[0]},{{1}}};};
    auto domain=[](const auto &){return true;};
    const auto a=solveManyNonsymmetricContact({0},1,3,{},initial,domain);
    EXPECT_FALSE(a.converged);EXPECT_EQ(a.status,"contact_accumulation_precision_budget");EXPECT_TRUE(a.velocity.empty());
    ASSERT_TRUE(a.precision_rejection.present);EXPECT_FALSE(a.precision_rejection.during_line_search);
    auto trial=[](const StdVec<Real> &v){return NonlinearContactEquation{{v[0]*v[0]*v[0]+1e-12*v[0]+.125},{{3*v[0]*v[0]+1e-12}}};};
    const auto b=solveManyNonsymmetricContact({0},1,3,{},trial,domain);
    EXPECT_FALSE(b.converged);EXPECT_EQ(b.status,"contact_accumulation_precision_budget");EXPECT_TRUE(b.velocity.empty());
    ASSERT_TRUE(b.precision_rejection.present);EXPECT_TRUE(b.precision_rejection.during_line_search);
    EXPECT_EQ(b.precision_rejection.trial_alpha,1);EXPECT_GT(b.precision_rejection.roundoff_bounds[0],1e-12);
}
