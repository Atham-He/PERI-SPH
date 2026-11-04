#include "pd_coulomb_contact.h"
#include <gtest/gtest.h>
using namespace SPH;
using namespace SPH::peridynamics;
namespace
{
CoulombContactConstraint plane(Real mu)
{
    CoulombContactConstraint c;c.weights={1};c.mu=mu;return c;
}
std::function<NonlinearContactEquation(const StdVec<Real> &)> particle(StdVec<Real> free_velocity,Real mass=1)
{
    return [=](const StdVec<Real> &v)
    {
        NonlinearContactEquation e;e.momentum.resize(3);e.jacobian.assign(3,StdVec<Real>(3));
        for(size_t i=0;i<3;++i){e.momentum[i]=mass*(v[i]-free_velocity[i]);e.jacobian[i][i]=mass;}
        return e;
    };
}
bool domain(const StdVec<Real> &) { return true; }
}

TEST(CoulombContact, ExactRationalBranchSurvivesNormalizationBeyondBinaryPrecision)
{
    using Q=ContactExactScalar;
    exact_contact_detail::Integer denominator=1;
    for(size_t i=0;i<100;++i)denominator*=10;
    const Q small=Q(1)/Q(denominator);
    const auto r=coulombDiskNaturalMapExact(Q(1),{Q(1),Q(0)},{small,Q(0)},Q(1));
    ASSERT_TRUE(r.within_budget);
    EXPECT_EQ(r.branch,CoulombDiskBranch::StickInterior);
    EXPECT_NEAR(r.residual[0],1e-100,1e-115);
    EXPECT_EQ(r.impulse_derivative.norm(),0);
    EXPECT_LT((r.velocity_derivative-ContactTangentMatrix::Identity()).norm(),1e-15);
}

TEST(CoulombContact, SharedSignedPartitionConstructsAllThreeExactTransposeRows)
{
    auto c=plane(.3);c.weights={-.2,1.2};
    const CoulombContactSystem system(StdVec<Real>(6),1,1,{c});
    const auto &g=system.coefficients();
    ASSERT_EQ(g.size(),3u);
    for(size_t a=0;a<3;++a)
    {
        EXPECT_DOUBLE_EQ(g[a][a],-.2);EXPECT_DOUBLE_EQ(g[a][3+a],1.2);
        for(size_t b=0;b<3;++b)if(a!=b){EXPECT_EQ(g[a][b],0);EXPECT_EQ(g[a][3+b],0);}
    }
    const StdVec<Real> v{.1,.2,.3,-.4,.5,-.6},p{.7,-.8,.9};
    long double gathered=0,scattered=0;
    for(size_t a=0;a<3;++a)gathered+=contactExactToWide(
        coulomb_contact_detail::exactDot(g[a],v)).convert_to<long double>()*p[a];
    for(size_t i=0;i<6;++i)for(size_t a=0;a<3;++a)scattered+=static_cast<long double>(v[i])*g[a][i]*p[a];
    EXPECT_NEAR(gathered,scattered,1e-15);
}

TEST(CoulombContact, AnalyticActualMassImpactHasTwoComponentSlideAndSeparateLoss)
{
    const CoulombContactSystem system({0,0,0},1,1,{plane(.3)});
    const auto r=solveCoulombContact(system,particle({-1,2,1}),domain);
    ASSERT_TRUE(r.converged)<<r.status;
    const Real magnitude=std::sqrt(Real(5));
    EXPECT_NEAR(r.velocity[0],0,1e-11);
    EXPECT_NEAR(r.velocity[1],2*(1-.3/magnitude),1e-11);
    EXPECT_NEAR(r.velocity[2],1-.3/magnitude,1e-11);
    EXPECT_NEAR(r.impulses[0],1,1e-11);
    EXPECT_NEAR(r.impulses[1],-.6/magnitude,1e-11);
    EXPECT_NEAR(r.impulses[2],-.3/magnitude,1e-11);
    const auto &last=r.trace.back().original;
    EXPECT_LT(last.KKT,1e-12);
    const Real kinetic=.5*(r.velocity[0]*r.velocity[0]+r.velocity[1]*r.velocity[1]+r.velocity[2]*r.velocity[2]);
    const Real algorithmic=.5*((r.velocity[0]+1)*(r.velocity[0]+1)+(r.velocity[1]-2)*(r.velocity[1]-2)+(r.velocity[2]-1)*(r.velocity[2]-1));
    const Real friction=last.disk[0].friction_loss;
    EXPECT_NEAR(3-kinetic,algorithmic+friction,1e-10);
    EXPECT_GT(friction,0);
}

TEST(CoulombContact, StrictStickAndZeroFrictionHaveDistinctPhysicalSolutions)
{
    const auto material=particle({-1,2,1});
    const auto stick=solveCoulombContact(CoulombContactSystem({0,0,0},1,1,{plane(3)}),material,domain);
    ASSERT_TRUE(stick.converged)<<stick.status;
    EXPECT_NEAR(stick.velocity[1],0,1e-11);EXPECT_NEAR(stick.velocity[2],0,1e-11);
    EXPECT_NEAR(stick.impulses[1],-2,1e-11);EXPECT_NEAR(stick.impulses[2],-1,1e-11);
    EXPECT_NEAR(stick.trace.back().original.disk[0].friction_loss,0,1e-11);
    const auto free=solveCoulombContact(CoulombContactSystem({0,0,0},1,1,{plane(0)}),material,domain);
    ASSERT_TRUE(free.converged)<<free.status;
    EXPECT_NEAR(free.velocity[0],0,1e-11);EXPECT_NEAR(free.velocity[1],2,1e-11);EXPECT_NEAR(free.velocity[2],1,1e-11);
    EXPECT_NEAR(free.impulses[1],0,1e-11);EXPECT_NEAR(free.impulses[2],0,1e-11);
}

TEST(CoulombContact, TangentBasisRotationPreservesPhysicalVelocityAndScatter)
{
    auto c=plane(.3),rotated=c;
    rotated.tangent[0]=ContactSpatialVector::UnitZ();rotated.tangent[1]=-ContactSpatialVector::UnitY();
    const auto material=particle({-1,2,1});
    const auto a=solveCoulombContact(CoulombContactSystem({0,0,0},1,1,{c}),material,domain);
    const auto b=solveCoulombContact(CoulombContactSystem({0,0,0},1,1,{rotated}),material,domain);
    ASSERT_TRUE(a.converged && b.converged)<<a.status<<" "<<b.status;
    for(size_t i=0;i<3;++i)EXPECT_NEAR(a.velocity[i],b.velocity[i],1e-11);
    EXPECT_NEAR(b.impulses[0],a.impulses[0],1e-11);
    EXPECT_NEAR(b.impulses[1],a.impulses[2],1e-11);
    EXPECT_NEAR(b.impulses[2],-a.impulses[1],1e-11);
}

TEST(CoulombContact, NormalizationAndSharedWallVelocityDoNotScalePhysicalMass)
{
    ManyContactSettings settings;settings.velocity_scale=2.5;
    const auto material=particle({-1,2,1});
    const auto a=solveCoulombContact(CoulombContactSystem({0,0,0},1,1,{plane(.3)}),material,domain);
    const auto b=solveCoulombContact(CoulombContactSystem({0,0,0},1,7,{plane(.3)},settings),material,domain,settings);
    ASSERT_TRUE(a.converged && b.converged)<<a.status<<" "<<b.status;
    for(size_t i=0;i<3;++i){EXPECT_NEAR(a.velocity[i],b.velocity[i],1e-10);EXPECT_NEAR(a.impulses[i],b.impulses[i],1e-10);}
    auto wall=plane(.3);wall.wall_tangent_velocity={.4,-.2};
    const auto common=solveCoulombContact(CoulombContactSystem({0,0,0},1,1,{wall}),particle({-1,.4,-.2}),domain);
    ASSERT_TRUE(common.converged)<<common.status;
    EXPECT_NEAR(common.velocity[1],.4,1e-11);EXPECT_NEAR(common.velocity[2],-.2,1e-11);
    EXPECT_NEAR(common.trace.back().original.disk[0].friction_loss,0,1e-11);
}

TEST(CoulombContact, FullOriginalNonsymmetricBlockMatchesJointFiniteDifferences)
{
    const DenseContactMatrix a{{2,.2,-.1},{.05,1.4,.3},{-.12,.1,1.2}};
    const auto material=[&](const StdVec<Real> &v)
    {
        NonlinearContactEquation e;e.jacobian=a;e.momentum.resize(3);
        for(size_t i=0;i<3;++i)e.momentum[i]=denseContactDot(a[i],v)-(.2+.1*i);
        return e;
    };
    const CoulombContactSystem system({0,0,0},1,1,{plane(.4)});
    const StdVec<Real> velocity{.1,.3,-.2},impulses{.5,-.1,.04};
    const auto r=system.evaluate(velocity,impulses,material);
    ASSERT_TRUE(r.precision_budget_ok);
    EXPECT_NEAR(r.jacobian[0][1],.2,1e-15);EXPECT_NEAR(r.jacobian[1][0],.05,1e-15);
    const Real epsilon=1e-6;
    for(size_t column=0;column<6;++column)
    {
        auto vp=velocity,vm=velocity,pp=impulses,pm=impulses;
        if(column<3){vp[column]+=epsilon;vm[column]-=epsilon;}
        else{pp[column-3]+=epsilon;pm[column-3]-=epsilon;}
        const auto plus=system.evaluate(vp,pp,material,false),minus=system.evaluate(vm,pm,material,false);
        for(size_t row=0;row<6;++row)
            EXPECT_NEAR((plus.residual[row]-minus.residual[row])/(2*epsilon),r.jacobian[row][column],3e-8);
    }
}

TEST(CoulombContact, CoupledNonsymmetricManufacturedRootRetainsPhysicalTranspose)
{
    const DenseContactMatrix a{{2,.2,-.1},{.05,1.4,.3},{-.12,.1,1.2}};
    const StdVec<Real> target{0,.2,-.3};const Real norm=std::sqrt(.13);
    const StdVec<Real> p{.6,-.24*.2/norm,.24*.3/norm};
    StdVec<Real> b(3);
    for(size_t i=0;i<3;++i)b[i]=denseContactDot(a[i],target)-p[i];
    const auto material=[&](const StdVec<Real> &v)
    {
        NonlinearContactEquation e;e.jacobian=a;e.momentum.resize(3);
        for(size_t i=0;i<3;++i)e.momentum[i]=denseContactDot(a[i],v)-b[i];return e;
    };
    const auto r=solveCoulombContact(CoulombContactSystem({0,0,0},1,1,{plane(.4)}),material,domain);
    ASSERT_TRUE(r.converged)<<r.status;
    for(size_t i=0;i<3;++i){EXPECT_NEAR(r.velocity[i],target[i],1e-10);EXPECT_NEAR(r.impulses[i],p[i],1e-10);}
}

TEST(CoulombContact, RedundantContactsPreserveUniqueTotalImpulse)
{
    const auto c=plane(.3);const auto material=particle({-1,2,1});
    const auto one=solveCoulombContact(CoulombContactSystem({0,0,0},1,1,{c}),material,domain);
    const auto two=solveCoulombContact(CoulombContactSystem({0,0,0},1,1,{c,c}),material,domain);
    ASSERT_TRUE(one.converged && two.converged)<<one.status<<" "<<two.status;
    ASSERT_EQ(two.trace.back().original.residual.size(),9u);
    for(size_t i=0;i<3;++i)
    {EXPECT_NEAR(one.velocity[i],two.velocity[i],1e-10);EXPECT_NEAR(one.impulses[i],two.impulses[i]+two.impulses[3+i],1e-10);}
}

TEST(CoulombContact, UnsupportedPrecisionReturnsNoRootAndOriginalWitness)
{
    const auto r=solveCoulombContact(CoulombContactSystem({0,0,0},1,1,{plane(.3)}),particle({-1e150,1e150,0}),domain);
    EXPECT_FALSE(r.converged);EXPECT_TRUE(r.velocity.empty());EXPECT_TRUE(r.impulses.empty());
    EXPECT_TRUE(r.precision_rejection.present);EXPECT_EQ(r.status,"original_contact_arithmetic_precision_budget");
}

TEST(CoulombContact, InvalidFramePartitionAndInitialGeometryReject)
{
    auto bad=plane(.3);bad.tangent[0]=bad.normal;
    EXPECT_THROW(CoulombContactSystem({0,0,0},1,1,{bad}),std::invalid_argument);
    bad=plane(.3);bad.weights={.9};
    EXPECT_THROW(CoulombContactSystem({0,0,0},1,1,{bad}),std::invalid_argument);
    bad=plane(-.3);EXPECT_THROW(CoulombContactSystem({0,0,0},1,1,{bad}),std::invalid_argument);
    EXPECT_THROW(CoulombContactSystem({-.1,0,0},1,1,{plane(.3)}),std::invalid_argument);
    ManyContactSettings loose;loose.relative_tolerance=1e-9;
    EXPECT_THROW(CoulombContactSystem({0,0,0},1,1,{plane(.3)},loose),std::invalid_argument);
}

TEST(CoulombContact, ExplicitLargerLinearEntryKeepsOldBoundAndOriginalResidual)
{
    for(size_t n:{size_t(447),size_t(783)})
    {
        DenseContactMatrix a(n,StdVec<Real>(n));StdVec<Real> rhs(n,1);
        for(size_t i=0;i<n;++i)a[i][i]=1;
        a[0][1]=.2;a[1][0]=-.1;
        EXPECT_THROW(solveConsistentContactSystem(a,rhs),std::invalid_argument);
        const auto r=solveBoundedConsistentContactSystem(a,rhs,1024);
        EXPECT_EQ(r.rank,n);EXPECT_LT(r.relative_residual,1e-12);
        for(size_t i=0;i<n;++i)EXPECT_NEAR(denseContactDot(a[i],r.value),rhs[i],1e-11);
    }
    EXPECT_THROW(solveBoundedConsistentContactSystem({{1,1},{2,2}},{3,7},1024),std::runtime_error);
}
