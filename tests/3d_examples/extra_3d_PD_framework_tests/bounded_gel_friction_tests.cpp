#include "pd_bounded_gel_friction.h"
#include <gtest/gtest.h>
using namespace SPH;
using namespace SPH::peridynamics;

TEST(BoundedGelFriction, ActualSphereUsesEverySignedThreeComponentContact)
{
    BoundedNormalConfig c;c.azimuths=8;const BoundedNormalScene scene(c);
    const auto rows=boundedGelFrictionContacts(scene,0,c.dt,{.3,.3},
        {ContactTangentVector::Zero(),ContactTangentVector(.02,.01)});
    const auto normal=scene.contacts(0,c.dt);
    const CoulombContactSystem system(flattenPD(scene.model.reference.position),c.dt,4*std::acos(Real(-1))/3,rows);
    EXPECT_EQ(system.dofs,99u);EXPECT_EQ(system.contacts_count,116u);EXPECT_EQ(system.unknowns,447u);
    for(size_t s=0;s<rows.size();++s)for(size_t i=0;i<99;++i)
        EXPECT_DOUBLE_EQ(system.coefficients()[3*s][i],normal[s].coefficients[i]);
    bool negative=false;
    for(const auto &row:rows)for(Real w:row.weights)negative|=w<0;
    EXPECT_TRUE(negative);
}

TEST(BoundedGelFriction, TwoAcceptedActualPDStepsCarryHistoryAndSeparatedWorkLedgers)
{
    BoundedNormalConfig c;c.azimuths=8;const BoundedNormalScene scene(c);
    auto old=initialBoundedGelFrictionState(scene);
    const std::array<ContactTangentVector,2> wall{ContactTangentVector::Zero(),ContactTangentVector(.02,.01)};
    const auto original=old;
    Real energy=0,friction=0;Vecd plate_impulse=Vecd::Zero();
    for(size_t step=0;step<2;++step)
    {
        const auto rows=boundedGelFrictionContacts(scene,old.gel.time,old.gel.time+c.dt,{.3,.3},wall);
        const auto r=advanceBoundedGelFriction(scene.model,old,c.dt,rows,scene.plateIDs());
        ASSERT_TRUE(r.accepted)<<r.status<<" "<<r.solve.status;
        EXPECT_EQ(r.next.gel.steps,step+1);EXPECT_EQ(r.next.gel.time,Real(step+1)*c.dt);
        EXPECT_LT(std::abs(r.energy_residual),1e-10);
        EXPECT_GT(r.normal_plate_impulse[1].norm(),0);
        EXPECT_GT(std::abs(r.tangent_plate_impulse[1][0]),0);
        EXPECT_GT(std::abs(r.tangent_plate_impulse[1][2]),0);
        EXPECT_GE(r.friction_loss,-1e-12);
        const auto gap=scene.denseIntervalGaps(old.gel.position,r.next.gel.position,old.gel.time,r.next.gel.time);
        EXPECT_GE(std::min(gap[0],gap[1]),-.01*c.stroke);
        energy+=r.normal_plate_work[0]+r.normal_plate_work[1]+r.tangent_plate_work[0]+r.tangent_plate_work[1];
        friction+=r.friction_loss;plate_impulse+=r.particle_plate_impulse[1];
        old=r.next;
    }
    EXPECT_NEAR(old.ledger(ContactLedger::LowerWork)+old.ledger(ContactLedger::UpperWork),energy,1e-15);
    EXPECT_NEAR(old.ledger(ContactLedger::FrictionLoss),friction,1e-15);
    EXPECT_NEAR(old.ledger(ContactLedger::UpperImpulseX),plate_impulse[0],1e-15);
    EXPECT_GT(old.ledger(ContactLedger::MaxwellLoss),0);
    EXPECT_GT(old.ledger(ContactLedger::AlgorithmicLoss),0);
    EXPECT_GT(old.tangent_slip_sum[0].norm()+old.tangent_slip_sum[1].norm(),0);
    EXPECT_EQ(original.gel.steps,0u);EXPECT_EQ(original.gel.time,0);
    EXPECT_EQ(original.ledger_sum,StdVec<Real>(size_t(ContactLedger::Count),0));
}

TEST(BoundedGelFriction, FreeTranslationWithNoContactsUsesActualMassAndHistory)
{
    BoundedNormalConfig c;c.mode=2;c.translation={.1,-.2,.3};const BoundedNormalScene scene(c);
    const auto old=initialBoundedGelFrictionState(scene);
    const auto r=advanceBoundedGelFriction(scene.model,old,c.dt,{},scene.plateIDs());
    ASSERT_TRUE(r.accepted)<<r.status;
    EXPECT_EQ(r.next.gel.steps,1u);
    for(size_t i=0;i<33;++i)
    {
        EXPECT_LT((r.next.gel.velocity[i]-old.gel.velocity[i]).norm(),1e-11);
        EXPECT_LT((r.next.gel.position[i]-old.gel.position[i]-c.dt*c.translation).norm(),1e-14);
    }
    EXPECT_LT(std::abs(r.energy_residual),1e-12);
    EXPECT_NEAR(r.friction_loss,0,1e-15);
}

TEST(BoundedGelFriction, FailedRootDoesNotMutateAcceptedHistoryOrLedgers)
{
    BoundedNormalConfig c;c.azimuths=8;const BoundedNormalScene scene(c);
    const auto old=initialBoundedGelFrictionState(scene);
    const auto rows=boundedGelFrictionContacts(scene,0,c.dt,{.3,.3},
        {ContactTangentVector::Zero(),ContactTangentVector(.02,.01)});
    ManyContactSettings settings;settings.max_newton_steps=1;
    const auto r=advanceBoundedGelFriction(scene.model,old,c.dt,rows,scene.plateIDs(),settings);
    EXPECT_FALSE(r.accepted);EXPECT_FALSE(r.solve.converged);
    EXPECT_TRUE(r.solve.velocity.empty());EXPECT_TRUE(r.solve.impulses.empty());
    EXPECT_EQ(old.gel.steps,0u);EXPECT_EQ(old.gel.time,0);
    EXPECT_EQ(old.ledger_sum,StdVec<Real>(size_t(ContactLedger::Count),0));
    for(const auto &v:old.last_tangent_impulse)EXPECT_EQ(v.norm(),0);
}

TEST(BoundedGelFriction, NonfiniteDiagnosticHistoryRejectsBeforeSolve)
{
    BoundedNormalConfig c;c.azimuths=8;const BoundedNormalScene scene(c);
    auto old=initialBoundedGelFrictionState(scene);
    const auto rows=boundedGelFrictionContacts(scene,0,c.dt,{.3,.3},
        {ContactTangentVector::Zero(),ContactTangentVector::Zero()});
    old.gel.last_normal_impulse[0]=std::numeric_limits<Real>::quiet_NaN();
    EXPECT_THROW(advanceBoundedGelFriction(scene.model,old,c.dt,rows,scene.plateIDs()),std::invalid_argument);
    old=initialBoundedGelFrictionState(scene);
    old.tangent_slip_correction[0][0]=std::numeric_limits<Real>::infinity();
    EXPECT_THROW(advanceBoundedGelFriction(scene.model,old,c.dt,rows,scene.plateIDs()),std::invalid_argument);
}

TEST(BoundedGelFriction, ExplicitCompensationRetainsSmallDiagnosticIncrements)
{
    Real sum=0,correction=0;
    addCompensatedContactValue(1e16,sum,correction);
    addCompensatedContactValue(1,sum,correction);
    Real copied_sum=sum,copied_correction=correction;
    addCompensatedContactValue(-1e16,sum,correction);
    addCompensatedContactValue(-1e16,copied_sum,copied_correction);
    EXPECT_EQ(sum+correction,1);
    EXPECT_DOUBLE_EQ(sum,copied_sum);EXPECT_DOUBLE_EQ(correction,copied_correction);
}
