#include "pd_tetra_gel_constraints.h"
#include <gtest/gtest.h>
using namespace SPH;
using namespace SPH::peridynamics;
namespace
{
NonlinearContactEquation quadratic(const StdVec<Real> &v)
{
    return {{.02*v[0]-.005*v[1], -.005*v[0]+.03*v[1]}, {{.02,-.005},{-.005,.03}}};
}
NonlinearContactEquation cubic(const StdVec<Real> &v)
{
    return {{v[0]+.2*v[0]*v[0]*v[0]+.1*v[1], 2*v[1]+.2*v[0]}, {{1+.6*v[0]*v[0],.1},{.2,2}}};
}
bool allDomain(const StdVec<Real> &) { return true; }
}

TEST(GeneralContactLinear, NonsymmetricAndIndefiniteAlgebraWithoutSymmetrization)
{
    const auto r=solveGeneralContactSystem({{2,3},{-1,1}},{5,0});
    EXPECT_NEAR(r.value[0],1,1e-15); EXPECT_NEAR(r.value[1],1,1e-15);
    EXPECT_LT(r.relative_residual,1e-15);
    const auto indefinite=solveGeneralContactSystem({{0,1},{1,0}},{2,3});
    EXPECT_NEAR(indefinite.value[0],3,1e-15); EXPECT_NEAR(indefinite.value[1],2,1e-15);
    EXPECT_THROW(solveGeneralContactSystem({{1,1},{1,1}},{1,1}),std::runtime_error);
}

TEST(PrimalDualContact, QuadraticLimitMatchesFullMetricControl)
{
    const auto r=solveNonsymmetricContact({0,0},.1,.03,{{{1,0},0,.01}},quadratic,allDomain);
    ASSERT_TRUE(r.converged);
    EXPECT_NEAR(r.velocity[0],.1,1e-14); EXPECT_NEAR(r.velocity[1],1.0/60,1e-14);
    EXPECT_NEAR(r.normal_impulse[0],.0019166666666666667,1e-14);
    EXPECT_LT(r.normalized_KKT,1e-12);
}

TEST(PrimalDualContact, NonlinearNonsymmetricEqualityRoot)
{
    const auto r=solveNonsymmetricContact({0,0},1,1,{{{1,0},0,.2}},cubic,allDomain);
    ASSERT_TRUE(r.converged);
    EXPECT_NEAR(r.velocity[0],.2,1e-14); EXPECT_NEAR(r.velocity[1],-.02,1e-14);
    EXPECT_NEAR(r.normal_impulse[0],.1996,1e-14);
}

TEST(TetraGel, RigidTranslationHasNoInternalForceOrHistoryChange)
{
    TetraGelProblem p;
    for(size_t node=0;node<4;++node) { p.velocity[3*node]=.1;p.velocity[3*node+1]=-.2;p.velocity[3*node+2]=.3; }
    const auto r=advanceTetraGel(p);
    ASSERT_TRUE(r.solve.converged);
    for(size_t i=0;i<12;++i) EXPECT_NEAR(r.velocity[i],p.velocity[i],1e-13);
    EXPECT_LT(r.evaluation.material.average_pk1.norm(),1e-12);
    EXPECT_LT(r.accepted_history.viscous.norm(),1e-14);
    EXPECT_LT(std::abs(r.energy_residual),1e-13);
}

TEST(TetraGel, TwoWallCompressionCommitsOnlyConvergedPhysicalState)
{
    TetraGelProblem p; StdVec<Real> lower(12,0),upper(12,0); lower[1]=1;upper[7]=-1;
    p.contacts={{lower,0,0},{upper,-1,-.99}};
    const auto old=p.history; const auto r=advanceTetraGel(p);
    ASSERT_TRUE(r.solve.converged);
    EXPECT_NEAR(r.position[1],0,1e-13); EXPECT_NEAR(r.position[7],.99,1e-13);
    EXPECT_GE(r.solve.normal_impulse[0],0); EXPECT_GE(r.solve.normal_impulse[1],0);
    for(Real residual:r.momentum_residual) EXPECT_LT(std::abs(residual),1e-12);
    EXPECT_LT(std::abs(r.energy_residual),1e-12);
    EXPECT_LT(std::abs(r.evaluation.transpose_work_residual),1e-12);
    EXPECT_LT(r.total_momentum_residual.norm(),1e-12);
    EXPECT_LT(r.internal_force_sum.norm(),1e-12);
    EXPECT_LT(r.angular_balance_old_cursor.norm(),1e-12);
    EXPECT_GT(r.algorithmic_loss,0);
    EXPECT_EQ(p.history.green,old.green); EXPECT_EQ(p.history.viscous,old.viscous);
    EXPECT_DOUBLE_EQ(p.history.dissipation,old.dissipation);
    auto invalid=p; invalid.node_mass[0]=0;
    EXPECT_THROW(advanceTetraGel(invalid),std::invalid_argument);
}

TEST(TetraGel, FullGlobalJacobianMatchesResidualDifferenceAndRetainsSkew)
{
    TetraGelProblem p; p.dt=.15;
    Matd f;f<<1.08,.2,-.03,.01,.93,.04,0,-.02,1.02;
    const Matd change=f-Matd::Identity(); StdVec<Real> v(12,0);
    for(size_t node=1;node<4;++node) for(size_t i=0;i<3;++i) v[3*node+i]=change(i,node-1)/p.dt;
    const auto base=evaluateTetraGel(p,v);
    EXPECT_GT(base.jacobian_asymmetry,1e-8);
    for(size_t j=0;j<12;++j)
    {
        auto plus=v,minus=v;plus[j]+=1e-6;minus[j]-=1e-6;
        const auto a=evaluateTetraGel(p,plus),b=evaluateTetraGel(p,minus);
        for(size_t i=0;i<12;++i) EXPECT_NEAR((a.equation.momentum[i]-b.equation.momentum[i])/2e-6,
                                           base.equation.jacobian[i][j],1e-8);
    }
    EXPECT_LT(std::abs(base.transpose_work_residual),1e-12);
}

TEST(PrimalDualContact, InfeasibleRowsDoNotProduceAcceptedState)
{
    const auto equation=[](const StdVec<Real> &v) { return NonlinearContactEquation{{v[0]},{{1}}}; };
    const auto r=solveNonsymmetricContact({0},1,1,{{{1},0,.1},{{-1},0,.1}},equation,allDomain);
    EXPECT_FALSE(r.converged); EXPECT_EQ(r.status,"no_KKT_feasible_active_set");
    EXPECT_TRUE(r.velocity.empty()); EXPECT_EQ(r.dependent_sets_skipped,1u);
}

TEST(PrimalDualContact, NewtonLimitIsTerminalWithoutAnotherAlgorithm)
{
    PrimalDualSettings settings;settings.max_newton_steps=1;
    const auto r=solveNonsymmetricContact({0,0},1,1,{{{1,0},0,.2}},cubic,allDomain,settings);
    EXPECT_FALSE(r.converged); EXPECT_EQ(r.status,"Newton_iteration_limit");
    EXPECT_TRUE(r.velocity.empty());
}
