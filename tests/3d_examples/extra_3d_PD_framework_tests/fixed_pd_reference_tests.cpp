#include "pd_bounded_gel.h"
#include "pd_surface_quadrature.h"
#include <gtest/gtest.h>
using namespace SPH;
using namespace SPH::peridynamics;
namespace
{
FixedPDReference reference()
{
    auto q=boundedSphereReference();return {q,StdVec<Real>(q.size(),4*std::acos(Real(-1))/(3*q.size())),1.25};
}
BoundedGelModel model(){return {reference(),{5,1,.5,.3},1,.1};}
Eigen::MatrixXd elastic(const FixedPDReference &ref,bool collapsed=false)
{
    const size_t n=ref.size();Eigen::MatrixXd k=Eigen::MatrixXd::Zero(3*n,3*n);
    Eigen::Matrix<Real,9,9> c;
    for(size_t a=0;a<3;++a)for(size_t b=0;b<3;++b)for(size_t d=0;d<3;++d)for(size_t e=0;e<3;++e)
        c(3*a+b,3*d+e)=Real(a==d && b==e)+Real(a==e && b==d)+(5-2./3)*Real(a==b && d==e);
    Matd global=Matd::Zero();for(size_t i=0;i<n;++i)global+=ref.volume[i]*ref.position[i]*ref.position[i].transpose();
    for(size_t i=0;i<n;++i)
    {
        auto b=pdGradientMatrix(ref,i);
        if(collapsed)
        {
            b.setZero();for(size_t j=0;j<n;++j)for(size_t a=0;a<3;++a)for(size_t d=0;d<3;++d)
                b(3*a+d,3*j+a)=ref.volume[j]*(global.inverse()*ref.position[j])[d];
        }
        k.noalias()+=ref.volume[i]*b.transpose()*c*b;
    }
    if(!collapsed)for(size_t i=0;i<n;++i)for(size_t j=0;j<n;++j)for(size_t a=0;a<3;++a)k(3*i+a,3*j+a)+=.1*ref.nonaffine_unit[i*n+j];
    return k;
}
}

TEST(FixedPDReference, FullSphereHasActualLocalShapesAndAffineReproduction)
{
    const auto r=reference();ASSERT_EQ(r.size(),33u);EXPECT_GT(r.target.size(),200u);
    EXPECT_EQ(r.offsets.size(),34u);EXPECT_GT(*std::min_element(r.shape_eigenvalue_ratio.begin(),r.shape_eigenvalue_ratio.end()),1e-3);
    Matd f;f<<1.1,.2,-.1,0,.9,.03,.05,0,1.02;StdVec<Vecd> q;
    for(const auto &x:r.position)q.push_back(f*x+Vecd(.13,-.27,.33));
    for(size_t i=0;i<r.size();++i)EXPECT_LT((r.deformation(i,q)-f).norm(),1e-12);
    EXPECT_LT(nonaffineEnergy(r,q,.1),1e-25);
    for(const auto &force:nonaffineForce(r,q,.1))EXPECT_LT(force.norm(),1e-12);
}

TEST(FixedPDReference, StabilizationIsPositiveAndItsForceIsTheEnergyGradient)
{
    const auto r=reference();auto q=r.position;q[16]+=Vecd(.02,-.03,.01);
    const auto force=nonaffineForce(r,q,.1);EXPECT_GT(nonaffineEnergy(r,q,.1),0);
    for(size_t i:{size_t(0),size_t(16),size_t(32)})for(size_t a=0;a<3;++a)
    {
        auto plus=q,minus=q;plus[i][a]+=1e-6;minus[i][a]-=1e-6;
        const Real derivative=(nonaffineEnergy(r,plus,.1)-nonaffineEnergy(r,minus,.1))/2e-6;
        EXPECT_NEAR(derivative,-force[i][a],1e-10);
    }
    Vecd net=Vecd::Zero(),torque=Vecd::Zero();for(size_t i=0;i<q.size();++i){net+=force[i];torque+=q[i].cross(force[i]);}
    EXPECT_LT(net.norm(),1e-12);EXPECT_LT(torque.norm(),1e-12);
}

TEST(FixedPDReference, ActualLinearizedSpectrumHasOnlySixRigidModes)
{
    const auto r=reference();const auto k=elastic(r);Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> spectrum(k);
    ASSERT_EQ(spectrum.info(),Eigen::Success);const auto &e=spectrum.eigenvalues();const Real scale=e[e.size()-1];
    EXPECT_GT(scale,0);EXPECT_LT(std::abs(e[0])/scale,1e-10);EXPECT_LT(std::abs(e[5])/scale,1e-10);EXPECT_GT(e[6]/scale,1e-9);
    // Replacing all local operators by one globally affine average preserves
    // affine tests but introduces unphysical nonaffine zero modes. Reject it.
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> bad(elastic(r,true));
    EXPECT_LT(std::abs(bad.eigenvalues()[6])/scale,1e-10);
}

TEST(FixedPDReference, CachedReferenceCorruptionIsRejectedWithoutRepair)
{
    const auto r=reference();auto bad=r;std::fill(bad.nonaffine_unit.begin(),bad.nonaffine_unit.end(),0);
    EXPECT_THROW(bad.validate(),std::invalid_argument);
    bad=r;bad.inverse_shape[0](0,0)*=1.01;
    EXPECT_THROW(bad.validate(),std::invalid_argument);
    bad=r;bad.active[0]=0;
    EXPECT_THROW(bad.validate(),std::invalid_argument);
    bad=r;++bad.generation;
    EXPECT_THROW(bad.validate(),std::invalid_argument);
    EXPECT_NO_THROW(r.validate());
}

TEST(FixedPDReference, SingularAndDisconnectedGeometriesAreRejected)
{
    EXPECT_THROW(FixedPDReference({Vecd(0,0,0),Vecd(1,0,0),Vecd(0,1,0),Vecd(1,1,0)},StdVec<Real>(4,1),2),std::invalid_argument);
    StdVec<Vecd> q{Vecd(0,0,0),Vecd(.1,0,0),Vecd(0,.1,0),Vecd(0,0,.1)};
    for(size_t i=0;i<4;++i)q.push_back(q[i]+Vecd(10,0,0));
    EXPECT_THROW(FixedPDReference(q,StdVec<Real>(8,1),.3),std::invalid_argument);
}

TEST(BoundedPDGel, GlobalNonlinearAssemblyHasTheTrueNonsymmetricTangent)
{
    const auto m=model();const auto old=m.initialState();const BoundedGelStepContext ctx(m,old,.04);
    StdVec<Vecd> v(old.position.size());Matd rate;rate<<.1,.3,-.1,.02,-.15,.03,.01,.08,.04;
    for(size_t i=0;i<v.size();++i)v[i]=rate*old.position[i];v[16]+=Vecd(.2,-.1,.15);
    const auto flat=flattenPD(v);const auto before=old.green;const auto r=ctx.evaluate(flat);
    EXPECT_GT(r.jacobian_asymmetry,1e-9);EXPECT_LT(std::abs(r.transpose_work_residual),1e-10);EXPECT_LT(std::abs(r.stabilization_work_residual),1e-10);
    for(size_t column:{size_t(0),size_t(8),size_t(49),size_t(98)})
    {
        auto plus=flat,minus=flat;plus[column]+=1e-6;minus[column]-=1e-6;
        const auto p=ctx.evaluate(plus),q=ctx.evaluate(minus);
        for(size_t row=0;row<flat.size();++row)EXPECT_NEAR((p.equation.momentum[row]-q.equation.momentum[row])/2e-6,r.equation.jacobian[row][column],2e-8);
    }
    EXPECT_EQ(old.green,before);EXPECT_EQ(old.steps,0u);
}

TEST(BoundedPDGel, ActualManyPointFreeTranslationAcceptsWithoutMaterialReset)
{
    const auto m=model();auto old=m.initialState();for(auto &v:old.velocity)v=Vecd(.1,-.2,.3);
    const auto r=advanceBoundedGel(m,old,.01,{},{});ASSERT_TRUE(r.accepted)<<r.solve.status;
    for(size_t i=0;i<old.position.size();++i)
    {
        EXPECT_LT((r.next.velocity[i]-old.velocity[i]).norm(),1e-10);
        EXPECT_LT((r.next.position[i]-old.position[i]-.01*old.velocity[i]).norm(),1e-10);
    }
    EXPECT_LT(std::abs(r.energy_residual),1e-10);EXPECT_EQ(old.steps,0u);EXPECT_EQ(r.next.steps,1u);
}

TEST(BoundedPDGel, PhysicalMassAndAcceptedHistoryMismatchCannotAdvance)
{
    auto m=model();auto old=m.initialState();m.node_mass[0]*=2;
    EXPECT_THROW(advanceBoundedGel(m,old,.01,{},{}),std::invalid_argument);
    m=model();old.viscous[0](0,0)=.2;old.viscous[0](1,1)=-.2;
    EXPECT_THROW(advanceBoundedGel(m,old,.01,{},{}),std::invalid_argument);
}

TEST(BoundedPDGel, PoleInclusiveFullSurfaceHasPositiveAreaAndSignedAffineMLS)
{
    const auto r=reference();const PolePanelSphereRule rule(Vecd::Zero(),1);ASSERT_EQ(rule.position.size(),58u);
    Real area=0;Matd moment=Matd::Zero();size_t poles=0;
    for(size_t s=0;s<rule.position.size();++s)
    {
        EXPECT_GT(rule.area[s],0);EXPECT_NEAR(rule.position[s].norm(),1,1e-14);
        area+=rule.area[s];moment+=rule.area[s]*rule.position[s]*rule.position[s].transpose();
        if(std::abs(rule.position[s][1])==1)++poles;
    }
    EXPECT_EQ(poles,2u);EXPECT_NEAR(area,4*std::acos(Real(-1)),1e-12);
    EXPECT_LT((moment-area/3*Matd::Identity()).norm(),1e-12);
    ReferenceSurfaceTransfer transfer(r.position.data(),r.size(),rule,1.25);
    EXPECT_LT(*std::min_element(transfer.weights.begin(),transfer.weights.end()),0);
    EXPECT_NO_THROW(transfer.validateStored(r.position.data(),r.size()));
    auto bad=transfer;for(auto &w:bad.weights)if(w<0)w=0;
    EXPECT_THROW(bad.validateStored(r.position.data(),r.size()),std::invalid_argument);
}

TEST(BoundedPDGel, SphereCompressionRetainsEveryFullSurfaceConstraintAndImpulse)
{
    const auto m=model();const auto old=m.initialState();const PolePanelSphereRule rule(Vecd::Zero(),1);
    const ReferenceSurfaceTransfer transfer(m.reference.position.data(),m.reference.size(),rule,1.25);
    StdVec<DenseContactConstraint> rows;StdVec<size_t> plates;
    for(size_t s=0;s<rule.position.size();++s)for(size_t plate=0;plate<2;++plate)
    {
        DenseContactConstraint row;row.coefficients.assign(3*old.position.size(),0);row.old_offset=-1;row.new_offset=-1+(plate?.001:0);
        for(size_t p=transfer.offsets[s];p<transfer.offsets[s+1];++p)row.coefficients[3*transfer.targets[p]+1]=(plate?-1:1)*transfer.weights[p];
        rows.push_back(row);plates.push_back(plate);
    }
    const auto r=advanceBoundedGel(m,old,.02,rows,plates);ASSERT_TRUE(r.accepted)<<r.solve.status;
    ASSERT_EQ(rows.size(),116u);EXPECT_EQ(r.solve.normal_impulse.size(),rows.size());
    const auto q=flattenPD(r.next.position);
    for(const auto &row:rows)EXPECT_GE(denseContactDot(row.coefficients,q)-row.new_offset,-1e-10);
    EXPECT_GT(r.particle_plate_impulse[0].norm(),0);EXPECT_GT(r.particle_plate_impulse[1].norm(),0);
    EXPECT_LT(std::abs(r.energy_residual),1e-10);EXPECT_GT(r.algorithmic_loss,0);
    EXPECT_EQ(old.steps,0u);EXPECT_EQ(r.next.steps,1u);EXPECT_EQ(r.next.time,.02);
}
