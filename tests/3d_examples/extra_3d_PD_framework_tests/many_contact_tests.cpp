#include "pd_many_contact.h"
#include <gtest/gtest.h>
using namespace SPH;
using namespace SPH::peridynamics;

TEST(ManyContact, RedundantLinearSolveRetainsOriginalEquationAndMinimumNorm)
{
    const auto r = solveConsistentContactSystem({{1,1},{2,2}}, {3,6});
    ASSERT_EQ(r.rank, 1u);
    EXPECT_NEAR(r.value[0], 1.5, 1e-13); EXPECT_NEAR(r.value[1], 1.5, 1e-13);
    EXPECT_LT(r.relative_residual, 1e-12);
    EXPECT_THROW(solveConsistentContactSystem({{1,1},{2,2}}, {3,7}), std::runtime_error);
}

TEST(ManyContact, ComplementarityRetainsSmallProductsAndRejectsWrongSigns)
{
    EXPECT_NEAR(contactFischerBurmeister(1,1e-20).value, -1e-20, 1e-34);
    EXPECT_GT(contactFischerBurmeister(-1,0).value, 0);
    EXPECT_GT(contactFischerBurmeister(0,-1).value, 0);
    EXPECT_EQ(contactFischerBurmeister(0,1).value, 0);
    EXPECT_NEAR(contactFischerBurmeister(0,0).gap_derivative, std::sqrt(.5)-1, 1e-16);
}

namespace
{
NonlinearContactEquation manufactured(const StdVec<Real> &v)
{
    NonlinearContactEquation e;
    e.momentum = {2*(v[0]-.2)+.4*(v[1]+.1)+.5*(v[0]*v[0]*v[0]-.008)+1,
        -.1*(v[0]-.2)+3*(v[1]+.1)+.25*(v[1]*v[1]*v[1]+.001)};
    e.jacobian = {{2+1.5*v[0]*v[0],.4},{-.1,3+.75*v[1]*v[1]}};
    return e;
}
StdVec<DenseContactConstraint> redundantRows()
{
    StdVec<DenseContactConstraint> rows;
    for (size_t i=0;i<16;++i) rows.push_back({{1,0},0,.2});
    rows.push_back({{-1,0},0,-.4}); rows.push_back({{0,1},0,-.1});
    return rows;
}
}

TEST(ManyContact, NonlinearNonsymmetricRootWithRedundantOpenAndBiactiveRows)
{
    const auto rows = redundantRows();
    const auto r = solveManyNonsymmetricContact({0,0},1,1,rows,manufactured,[](const auto &){return true;});
    ASSERT_TRUE(r.converged) << r.status;
    EXPECT_NEAR(r.velocity[0],.2,1e-12); EXPECT_NEAR(r.velocity[1],-.1,1e-12);
    Real impulse = 0;
    for(size_t i=0;i<16;++i) { impulse += r.normal_impulse[i]; EXPECT_GE(r.normal_impulse[i],-1e-12); }
    EXPECT_NEAR(impulse,1,1e-12); EXPECT_NEAR(r.normal_impulse[16],0,1e-12); EXPECT_NEAR(r.normal_impulse[17],0,1e-12);
    EXPECT_LT(r.normalized_KKT,1e-12);
    // The path may stop before the near-root derivative loses numerical
    // rank. Test the exactly redundant root system itself, not that timing.
    DenseContactMatrix h(20,StdVec<Real>(20,0));
    const auto root=manufactured({.2,-.1});
    for(size_t i=0;i<2;++i)for(size_t j=0;j<2;++j)h[i][j]=root.jacobian[i][j];
    for(size_t s=0;s<18;++s)
    {
        const auto fb=contactFischerBurmeister(s<16?Real(1)/16:0,s==16?.2:0);
        for(size_t i=0;i<2;++i){h[i][2+s]=-rows[s].coefficients[i];h[2+s][i]=fb.gap_derivative*rows[s].coefficients[i];}
        h[2+s][2+s]=fb.impulse_derivative;
    }
    StdVec<Real> perturbation(20,.03),rhs(20);perturbation[0]=.1;perturbation[1]=.2;
    for(size_t i=0;i<20;++i)rhs[i]=denseContactDot(h[i],perturbation);
    const auto redundant=solveConsistentContactSystem(h,rhs);
    EXPECT_EQ(redundant.rank,5u); EXPECT_LT(redundant.relative_residual,1e-12);
}

TEST(ManyContact, RowPermutationPreservesPhysicalScatterWithoutClaimingUniqueDual)
{
    auto rows = redundantRows(); std::reverse(rows.begin(),rows.end());
    const auto r = solveManyNonsymmetricContact({0,0},1,1,rows,manufactured,[](const auto &){return true;});
    ASSERT_TRUE(r.converged) << r.status;
    StdVec<Real> scatter(2,0);
    for(size_t s=0;s<rows.size();++s)for(size_t i=0;i<2;++i)scatter[i]+=rows[s].coefficients[i]*r.normal_impulse[s];
    EXPECT_NEAR(scatter[0],1,1e-12); EXPECT_NEAR(scatter[1],0,1e-12);
    EXPECT_NEAR(r.velocity[0],.2,1e-12); EXPECT_NEAR(r.velocity[1],-.1,1e-12);
}

TEST(ManyContact, NewtonAndDomainFailureReturnNoAcceptedPhysicalState)
{
    ManyContactSettings s; s.max_newton_steps=1;
    const auto r = solveManyNonsymmetricContact({0,0},1,1,redundantRows(),manufactured,[](const auto &){return true;},s);
    ASSERT_FALSE(r.converged); EXPECT_TRUE(r.velocity.empty()); EXPECT_TRUE(r.normal_impulse.empty());
    EXPECT_EQ(r.status,"Newton_iteration_limit"); ASSERT_EQ(r.trace.back().iteration,1u);
    EXPECT_GT(std::max(r.trace.back().merit,r.trace.back().KKT),s.relative_tolerance);
    const auto d = solveManyNonsymmetricContact({0,0},1,1,redundantRows(),manufactured,
        [](const auto &v){return v[0]==0 && v[1]==0;});
    EXPECT_FALSE(d.converged); EXPECT_EQ(d.status,"line_search_limit");
    EXPECT_TRUE(d.velocity.empty()); EXPECT_EQ(d.trace.back().domain_rejections,32u);
}

TEST(ManyContact, InvalidDataIsRejectedBeforePhysicalCommit)
{
    auto rows=redundantRows(); rows[0].coefficients={0,0};
    EXPECT_THROW(solveManyNonsymmetricContact({0,0},1,1,rows,manufactured,[](const auto &){return true;}),std::invalid_argument);
    EXPECT_THROW(solveManyNonsymmetricContact({0,0},1,0,{},manufactured,[](const auto &){return true;}),std::invalid_argument);
    EXPECT_THROW(solveManyNonsymmetricContact({0,0},1,1,{},manufactured,[](const auto &){return false;}),std::invalid_argument);
}
