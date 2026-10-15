#include "sphinxsys.h"
#include <gtest/gtest.h>
#include <limits>
#include "io_vtk.hpp"
#include <fstream>
#include <sstream>

using namespace SPH;

class PDTestBody
{
  public:
    SPHSystem system;
    PDBody body;
    explicit PDTestBody(Real scale = 1)
        : system(BoundingBoxd(Vecd(-scale, -scale, -scale), Vecd(scale, scale, scale)), 0.1 * scale),
          body(system, makeShared<GeometricShapeBox>(Transform(), Vecd(0.3 * scale, 0.3 * scale, 0.3 * scale), "TestSolid")) {}
    UniquePtr<InnerRelation> inner;
    template <class Material, class... Args> void generate(Args... args)
    {
        body.defineMatterMaterial<Material>(args...);
        body.generateParticles<BaseParticles, Lattice>();
        inner = makeUnique<InnerRelation>(body);
        system.initializeSystemCellLinkedLists();
        system.initializeSystemConfigurations();
    }
};

TEST(PDTime, RejectsInvalidAndClipsToEveryBoundary)
{
    EXPECT_DOUBLE_EQ(peridynamics::boundedTimeStep(0.7, {1.0, 0.2, 0.3}), 0.2);
    EXPECT_THROW(peridynamics::boundedTimeStep(0.0, {1.0}), std::runtime_error);
    EXPECT_THROW(peridynamics::boundedTimeStep(std::numeric_limits<Real>::quiet_NaN(), {1.0}), std::runtime_error);
    EXPECT_THROW(peridynamics::boundedTimeStep(1.0, {0.0}), std::runtime_error);
}

TEST(PDMotion, ExactDisplacementAndBothWallVelocityFields)
{
    PDTestBody f;
    f.generate<HughesWingetSolid>(1000.0, 1e6, 0.25);
    auto &p = f.body.getBaseParticles();
    Vecd initial = p.getVariableDataByName<Vecd>("Position")[0];
    SimpleDynamics<solid_dynamics::PrescribedTranslation> move(f.body, 0.0, 0.3, Vecd(0.0, -0.6, 0.0));
    Real &time = *f.system.getSystemVariableDataByName<Real>("PhysicalTime");
    for (Real t : {0.0, 0.02, 0.17, 0.29, 0.3, 0.5})
    {
        time = t;
        move.exec();
        Vecd expected = initial + Vecd(0, -2.0 * std::min(t, Real(0.3)), 0);
        EXPECT_LT((p.getVariableDataByName<Vecd>("Position")[0] - expected).norm(), 1e-13);
        Vecd velocity = t < 0.3 ? Vecd(0, -2, 0) : Vecd::Zero();
        EXPECT_LT((p.getVariableDataByName<Vecd>("Velocity")[0] - velocity).norm(), 1e-13);
        EXPECT_LT((p.getVariableDataByName<Vecd>("AverageVelocity")[0] - velocity).norm(), 1e-13);
    }
}

TEST(PDElasticity, HughesWingetPureRotationPreservesStress)
{
    HughesWingetSolid material(1000, 1e6, 0.25);
    Real theta = 0.23;
    Matd G = Matd::Zero();
    G(0, 1) = -2.0 * std::tan(theta / 2.0);
    G(1, 0) = -G(0, 1);
    Matd stress = Vecd(100, 200, 300).asDiagonal();
    Matd R = Matd::Identity();
    R(0, 0) = R(1, 1) = std::cos(theta);
    R(0, 1) = -std::sin(theta); R(1, 0) = std::sin(theta);
    EXPECT_LT((material.StressHW(G, stress, 0) - R * stress * R.transpose()).norm(), 1e-10);
}

TEST(PDGeometry, UniformDeformationReproducedAtSurfaceAndInterior)
{
    PDTestBody f;
    f.generate<HughesWingetSolid>(1000, 1e6, 0.25);
    InteractionWithUpdate<solid_dynamics::NosbPDShapeMatrix> shape(*f.inner);
    InteractionWithUpdate<solid_dynamics::NosbPDSecondStep> second(*f.inner);
    shape.exec();
    Matd A = Matd::Identity(); A(0, 0) = 1.07; A(1, 1) = 0.94; A(0, 2) = 0.03;
    auto &p = f.body.getBaseParticles();
    auto *position = p.getVariableDataByName<Vecd>("Position");
    for (size_t i = 0; i < p.TotalRealParticles(); ++i) position[i] = A * position[i] + Vecd(0.1, 0.2, -0.1);
    second.exec(1e-5);
    auto *F = p.getVariableDataByName<Matd>("DeformationGradient");
    for (size_t i = 0; i < p.TotalRealParticles(); ++i) EXPECT_LT((F[i] - A).norm(), 1e-10);
}

TEST(PDFracture, DerivedShearAndPressureCriteriaDispatchThroughBase)
{
    PDTestBody f;
    f.generate<HughesWingetSolid>(1000, 1e6, 0.25);
    solid_dynamics::BondBreakBySigma1andSigma3 shear(*f.inner, 1000, 10);
    solid_dynamics::BondBreakBySigma1andNorm1 pressure(*f.inner, 1000, 10);
    solid_dynamics::BondBreakByPrinStress *base = &shear;
    Matd stress = Matd::Zero(); stress(0, 1) = stress(1, 0) = 20;
    Real stretch = 0;
    EXPECT_FALSE(base->checkBondLive(stress, stretch));
    base = &pressure;
    stress = -20 * Matd::Identity();
    EXPECT_FALSE(base->checkBondLive(stress, stretch));
}

TEST(PDJ2, NonlinearHardeningReturnsToYieldSurface)
{
    PDTestBody f;
    f.generate<J2PlasticityforPD>(1000, 1e6, 0.25, 1000, 10000, 2000);
    auto &material = DynamicCast<J2PlasticityforPD>(&f, f.body.getMatterMaterial());
    auto &p = f.body.getBaseParticles();
    Matd G = Matd::Zero(); G(0, 0) = 0.02; G(1, 1) = G(2, 2) = -0.01;
    Matd epsilon_p = Matd::Zero();
    Matd stress = material.PlasticConstitutiveRelation(G, Matd::Zero(), epsilon_p, 0, 1e-3);
    Real alpha = p.getVariableDataByName<Real>("IsotropicHardeningParam")[0];
    Matd back = material.KinHardeningFunc(p.getVariableDataByName<Matd>("KinematicHardeningParam")[0]);
    EXPECT_GT(alpha, 0);
    EXPECT_NEAR(material.YieldFunc(stress, material.IsoHardeningFunc(alpha), back), 0, 1e-6);
    EXPECT_GT((stress.array() * epsilon_p.array()).sum(), 0);
}

static Real viscousStressNorm(Real viscosity)
{
    PDTestBody f;
    f.generate<ViscoPlasticSolidforPD>(1000, 1e6, 0.25, 1000, viscosity, 10000);
    auto &material = DynamicCast<ViscoPlasticSolidforPD>(&f, f.body.getMatterMaterial());
    Matd G = Matd::Zero(); G(0, 0) = 0.02; G(1, 1) = G(2, 2) = -0.01;
    Matd ep = Matd::Zero();
    return material.PlasticConstitutiveRelation(G, Matd::Zero(), ep, 0, 1e-3).norm();
}

TEST(PDViscoplasticity, ViscosityIncreasesOverstressAndInvalidDtFails)
{
    EXPECT_GT(viscousStressNorm(1000), viscousStressNorm(0));
    PDTestBody f;
    f.generate<ViscoPlasticSolidforPD>(1000, 1e6, 0.25, 1000, 1000);
    auto &material = DynamicCast<ViscoPlasticSolidforPD>(&f, f.body.getMatterMaterial());
    Matd G = Matd::Zero(); G(0, 0) = 0.02; G(1, 1) = G(2, 2) = -0.01;
    Matd ep = Matd::Zero();
    EXPECT_THROW(material.PlasticConstitutiveRelation(G, Matd::Zero(), ep, 0, 0), std::invalid_argument);
}

TEST(PDDruckerPrager, HydrostaticApexAndNearbyPathStayFiniteAndConsistent)
{
    for (Real shear : {0.0, 1e-8, 0.001})
    {
        PDTestBody f;
        f.generate<DruckerPragerPlasticityforPD>(1000, 1e6, 0.25, 1000, 0.2, 10000);
        auto &material = DynamicCast<DruckerPragerPlasticityforPD>(&f, f.body.getMatterMaterial());
        Matd G = 0.01 * Matd::Identity(); G(0, 1) = G(1, 0) = shear;
        Matd ep = Matd::Zero();
        Matd stress = material.PlasticConstitutiveRelation(G, Matd::Zero(), ep, 0, 1e-3);
        Real alpha = f.body.getBaseParticles().getVariableDataByName<Real>("IsotropicHardeningParam")[0];
        EXPECT_TRUE(stress.allFinite()); EXPECT_TRUE(ep.allFinite()); EXPECT_GT(alpha, 0);
        EXPECT_NEAR(material.YieldFunc(stress, material.IsoHardeningFunc(alpha)), 0, 1e-6);
        EXPECT_LT((stress - stress.trace() / 3.0 * Matd::Identity()).norm(), 1e-8);
    }
}

TEST(PDContact, FrictionOpposesSlipAndIsInvariantToCommonVelocity)
{
    Vecd fn(0, 10, 0), slip(2, 0, 0), common(3, 4, 5);
    Vecd force = solid_dynamics::pdContactPairForce(fn, slip, 0.3, 0.01);
    EXPECT_NEAR(force[0], -3, 1e-12); EXPECT_DOUBLE_EQ(force[1], 10);
    EXPECT_LT((force - fn).dot(slip), 0);
    EXPECT_LT((force - solid_dynamics::pdContactPairForce(fn, (slip + common) - common, 0.3, 0.01)).norm(), 1e-12);
}

TEST(PDContact, AppliedForceMatchesRecordedForceAndWallReaction)
{
    PDTestBody f;
    f.generate<HughesWingetSolid>(1000, 1e6, 0.25);
    SolidBody wall(f.system, makeShared<GeometricShapeBox>(Transform(Vecd(0, -0.4, 0)), Vecd(0.5, 0.15, 0.5), "Wall"));
    wall.defineMatterMaterial<Solid>(1000, 1e6);
    wall.generateParticles<BaseParticles, Lattice>();
    SimpleDynamics<NormalDirectionFromBodyShape> wall_normals(wall);
    wall_normals.exec();
    auto &wall_particles = wall.getBaseParticles();
    auto *wall_position = wall_particles.getVariableDataByName<Vecd>("Position");
    // Introduce a controlled sub-spacing contact overlap, independently of lattice clipping.
    for (size_t i = 0; i < wall_particles.TotalRealParticles(); ++i) wall_position[i][1] += 0.025;
    ContactRelation contact(f.body, {&wall}), reverse(wall, {&f.body});
    InteractionDynamics<solid_dynamics::ContactRepulsionFactorforPD> factor(contact);
    InteractionDynamics<solid_dynamics::ContactRepulsionFromWallforPD> force(contact, 0.3, 0.01);
    InteractionDynamics<solid_dynamics::ContactRepulsionToWallforPD> reaction(reverse, 0.3, 0.01);
    f.system.initializeSystemCellLinkedLists(); f.system.initializeSystemConfigurations();
    auto &p = f.body.getBaseParticles();
    auto *v = p.getVariableDataByName<Vecd>("Velocity");
    auto *a = p.getVariableDataByName<Vecd>("PDAccelerationPrior");
    auto *volume = p.getVariableDataByName<Real>("VolumetricMeasure");
    for (size_t i = 0; i < p.TotalRealParticles(); ++i) {v[i] = Vecd(1, 0, 0); a[i] = Vecd(0, -9.8, 0);}
    factor.exec(); force.exec(); reaction.exec();
    auto *recorded = p.getVariableDataByName<Vecd>("ContactRepulsionFromWall");
    Vecd sum = Vecd::Zero();
    for (size_t i = 0; i < p.TotalRealParticles(); ++i)
    {
        EXPECT_LT(((a[i] - Vecd(0, -9.8, 0)) * 1000 * volume[i] - recorded[i]).norm(), 1e-9);
        sum += recorded[i];
    }
    EXPECT_GT(sum.norm(), 0);
    auto &wp = wall.getBaseParticles();
    auto *r = wp.getVariableDataByName<Vecd>("ContactRepulsionToWall");
    for (size_t i = 0; i < wp.TotalRealParticles(); ++i) sum += r[i];
    EXPECT_LT(sum.norm(), 1e-8);
}

TEST(PDAcousticDamping, ConsumesForceOncePreservesAccumulatorAndDissipates)
{
    SPHSystem system(BoundingBoxd(Vecd(-1, -1, -1), Vecd(1, 1, 1)), 0.1);
    FluidBody body(system, makeShared<GeometricShapeBox>(Transform(), Vecd(0.3, 0.3, 0.3), "Fluid"));
    body.defineMatterMaterial<WeaklyCompressibleFluid>(1000, 10);
    body.generateParticles<BaseParticles, Lattice>();
    InnerRelation inner(body);
    Dynamics1Level<fluid_dynamics::Integration1stHalfInnerRiemann> integration(inner);
    InteractionDynamics<fluid_dynamics::VelocityDivergenceInner> divergence(inner);
    InteractionDynamics<fluid_dynamics::AcousticDamperInnerforPD> damping(inner);
    damping.setCoeffAcousticDamper(1000, 10, 0.1);
    system.initializeSystemCellLinkedLists(); system.initializeSystemConfigurations();
    auto &p = body.getBaseParticles();
    auto *v = p.getVariableDataByName<Vecd>("Velocity");
    auto *x = p.getVariableDataByName<Vecd>("Position");
    auto *mass = p.getVariableDataByName<Real>("Mass");
    auto *upstream = p.getVariableDataByName<Vecd>("Force");
    StdVec<Vecd> old(p.TotalRealParticles());
    Real before = 0, after = 0; Vecd total_force = Vecd::Zero();
    for (size_t i = 0; i < p.TotalRealParticles(); ++i)
    {v[i] = old[i] = x[i]; upstream[i] = Vecd(1, 2, 3); before += mass[i] * v[i].squaredNorm();}
    divergence.exec(); damping.exec(0.001);
    auto *diagnostic = p.getVariableDataByName<Vecd>("AcousticDampingForceInner");
    for (size_t i = 0; i < p.TotalRealParticles(); ++i)
    {
        EXPECT_LT((mass[i] * (v[i] - old[i]) / 0.001 - diagnostic[i]).norm(), 1e-10);
        EXPECT_EQ(upstream[i], Vecd(1, 2, 3));
        after += mass[i] * v[i].squaredNorm(); total_force += diagnostic[i];
    }
    EXPECT_LT(after, before); EXPECT_LT(total_force.norm(), 1e-9);
}

TEST(PDDiagnostics, RejectsInjectedNaNAndNegativeDensity)
{
    PDTestBody f;
    f.generate<HughesWingetSolid>(1000, 1e6, 0.25);
    SimpleDynamics<solid_dynamics::NosbPDFirstStep> first(f.body);
    SimpleDynamics<peridynamics::ParticleStateCheck> check(f.body);
    check.addScalar("Density", 0.0);
    auto &p = f.body.getBaseParticles();
    auto *x = p.getVariableDataByName<Vecd>("Position");
    x[0][2] = std::numeric_limits<Real>::quiet_NaN();
    EXPECT_THROW(check.exec(), std::runtime_error);
    x[0][2] = 0;
    p.getVariableDataByName<Real>("Density")[0] = -1;
    EXPECT_THROW(check.exec(), std::runtime_error);
}

TEST(PDTime, AccelerationAndSpeedReduceTheSolidIncrement)
{
    PDTestBody f;
    f.generate<HughesWingetSolid>(1000, 1e6, 0.25);
    ReduceDynamics<solid_dynamics::PDTimeStep> step(f.body, 0.1);
    Real initial = step.exec();
    auto &p = f.body.getBaseParticles();
    p.getVariableDataByName<Vecd>("Velocity")[0] = Vecd(1000, 0, 0);
    EXPECT_LT(step.exec(), initial);
    p.getVariableDataByName<Vecd>("Velocity")[0].setZero();
    p.getVariableDataByName<Vecd>("PDAcceleration")[0] = Vecd(1e12, 0, 0);
    EXPECT_LT(step.exec(), initial);
}

TEST(PDMaterials, InvalidElasticAndViscousParametersAreRejected)
{
    EXPECT_THROW(HughesWingetSolid(1000, -1, 0.25), std::invalid_argument);
    EXPECT_THROW(HughesWingetSolid(1000, 1e6, 0.5), std::invalid_argument);
    EXPECT_THROW(ViscoPlasticSolidforPD(1000, 1e6, 0.25, 1000, -1), std::invalid_argument);
}

TEST(PDContact, NormalWallMotionDoesNotCreateTangentialFriction)
{
    Vecd force(3, 10, 1), relative(0, -2, 0), wall_normal(0, 1, 0);
    EXPECT_LT((solid_dynamics::pdContactPairForce(force, relative, wall_normal, 0.3, 0.01) - force).norm(), 1e-12);
}

static StdVec<Vecd> hourglassAcceleration(Real scale, bool perturb, Vecd &net_force)
{
    PDTestBody f(scale);
    f.generate<HughesWingetSolid>(1000, 1e6, 0.25);
    InteractionWithUpdate<solid_dynamics::NosbPDShapeMatrix> shape(*f.inner);
    InteractionWithUpdate<solid_dynamics::NosbPDSecondStep> second(*f.inner);
    InteractionDynamics<solid_dynamics::LittleWoodHourGlassControl> hourglass(*f.inner);
    shape.exec();
    auto &p = f.body.getBaseParticles();
    auto *x = p.getVariableDataByName<Vecd>("Position");
    if (perturb)
    {
        size_t center = 0;
        for (size_t i = 1; i < p.TotalRealParticles(); ++i)
            if (x[i].squaredNorm() < x[center].squaredNorm()) center = i;
        x[center] += scale * Vecd(0.001, 0.0004, -0.0003);
    }
    else
    {
        Matd A = Matd::Identity(); A(0,0) = 1.02; A(1,1) = 0.99; A(0,2) = 0.01;
        for (size_t i = 0; i < p.TotalRealParticles(); ++i) x[i] = A * x[i] + scale * Vecd(0.1, 0.2, -0.1);
    }
    second.exec(1e-5); hourglass.exec();
    auto *a = p.getVariableDataByName<Vecd>("PDAcceleration");
    auto *mass = p.getVariableDataByName<Real>("Mass");
    StdVec<Vecd> result(p.TotalRealParticles());
    net_force.setZero();
    for (size_t i = 0; i < p.TotalRealParticles(); ++i) {result[i] = a[i]; net_force += mass[i] * a[i];}
    return result;
}

TEST(PDHourglass, GeometricScalingHasInverseLengthAcceleration)
{
    Vecd force;
    auto small = hourglassAcceleration(1, true, force);
    auto large = hourglassAcceleration(2, true, force);
    ASSERT_EQ(small.size(), large.size());
    Real norm = 0, error = 0;
    for (size_t i = 0; i < small.size(); ++i) {norm += small[i].squaredNorm(); error += (small[i] - 2 * large[i]).squaredNorm();}
    ASSERT_GT(norm, 1e-20);
    EXPECT_LT(std::sqrt(error / norm), 1e-9);
}

TEST(PDHourglass, VanishesForAffineMotionAndBalancesMomentum)
{
    Vecd net;
    auto affine = hourglassAcceleration(1, false, net);
    for (const Vecd &a : affine) EXPECT_LT(a.norm(), 1e-8);
    auto perturbed = hourglassAcceleration(1, true, net);
    Real scale = 0;
    for (const Vecd &a : perturbed) scale += a.norm();
    EXPECT_LT(net.norm(), 1e-9 * std::max(scale, Real(1)));
}


TEST(PDEnergy, LinearEOSStorageMatchesPressureWorkAndIsPositive)
{
    Real rho0=935, c=25;
    EXPECT_DOUBLE_EQ(peridynamics::linearEOSSpecificEnergy(rho0, rho0, c), 0);
    for (Real ratio : {0.9, 0.99999999, 1.00000001, 1.1})
    {
        Real rho=rho0*ratio;
        EXPECT_GE(peridynamics::linearEOSSpecificEnergy(rho,rho0,c),0);
        Real h=1e-3*rho0;
        if (std::abs(ratio-1)>0.01)
        {
            h=1e-6*rho0;
            Real derivative=(peridynamics::linearEOSSpecificEnergy(rho+h,rho0,c)-peridynamics::linearEOSSpecificEnergy(rho-h,rho0,c))/(2*h);
            EXPECT_NEAR(derivative,c*c*(rho-rho0)/(rho*rho),1e-9);
        }
    }
    EXPECT_THROW(peridynamics::linearEOSSpecificEnergy(0,rho0,c),std::invalid_argument);
}

TEST(PDEnergy, SolidKickWorkPartitionsKineticIncrementExactly)
{
    PDTestBody f; f.generate<HughesWingetSolid>(1000,1e6,0.25);
    SimpleDynamics<solid_dynamics::NosbPDFirstStep> first(f.body);
    SimpleDynamics<solid_dynamics::NosbPDFourthStep> fourth(f.body);
    auto &p=f.body.getBaseParticles();
    auto *v=p.getVariableDataByName<Vecd>("Velocity");
    auto *a=p.getVariableDataByName<Vecd>("PDAcceleration");
    auto *prior=p.getVariableDataByName<Vecd>("PDAccelerationPrior");
    auto *m=p.getVariableDataByName<Real>("Mass");
    const Vecd ai(2,-3,1), ah(-1,0.5,0.2), ad(-0.2,-0.1,0);
    Real before=0,after=0,work=0;
    for(size_t i=0;i<p.TotalRealParticles();++i)
    {
        v[i]=Vecd(0.1,0.2,-0.1); prior[i]=Vecd(0,-9.8,0); a[i]=prior[i]+ai+ah+ad;
        p.getVariableDataByName<Vecd>("PDInternalForce")[i]=m[i]*ai;
        p.getVariableDataByName<Vecd>("HourglassForce")[i]=m[i]*ah;
        p.getVariableDataByName<Vecd>("PDNumericalDampingForce")[i]=m[i]*ad;
        before+=0.5*m[i]*v[i].squaredNorm();
    }
    fourth.exec(0.001);
    for(size_t i=0;i<p.TotalRealParticles();++i)
    {
        after+=0.5*m[i]*v[i].squaredNorm();
        for(const char *name : {"PDInternalWork","PDExternalWork","PDHourglassWork","PDDampingWork"})
            work+=p.getVariableDataByName<Real>(name)[i];
    }
    EXPECT_NEAR(after-before,work,1e-12);
}

TEST(PDEnergy, FluidKickAndAcousticDamperWorkAreMeasuredIndependently)
{
    PDTestBody f; f.generate<HughesWingetSolid>(1000,1e6,0.25);
    auto &p=f.body.getBaseParticles();
    p.registerStateVariableData<Vecd>("Velocity");
    p.registerStateVariableData<Vecd>("ViscousForce");
    p.registerStateVariableData<Real>("Density",1000.0);
    SimpleDynamics<peridynamics::FluidKineticWork> ledger(f.body);
    InteractionDynamics<fluid_dynamics::VelocityDivergenceInner> divergence(*f.inner);
    InteractionDynamics<fluid_dynamics::AcousticDamperInnerforPD> damping(*f.inner);
    damping.setCoeffAcousticDamper(1000,10,0.1);
    auto *v=p.getVariableDataByName<Vecd>("Velocity");
    auto *m=p.getVariableDataByName<Real>("Mass");
    Real before=0,after=0,work=0,dt=0.001;
    for(size_t i=0;i<p.TotalRealParticles();++i)
    {
        v[i]=p.getVariableDataByName<Vecd>("Position")[i];
        before+=0.5*m[i]*v[i].squaredNorm();
        p.getVariableDataByName<Vecd>("Force")[i]=m[i]*Vecd(0.2,0.3,0.4);
        p.getVariableDataByName<Vecd>("ForcePrior")[i]=m[i]*Vecd(-0.1,0.2,0.1);
        v[i]+=dt*Vecd(0.1,0.5,0.5);
    }
    ledger.exec(dt); divergence.exec(); damping.exec(dt);
    for(size_t i=0;i<p.TotalRealParticles();++i)
    {
        after+=0.5*m[i]*v[i].squaredNorm();
        for(const char *name : {"FluidPressureWork","FluidPriorWork","AcousticDampingWorkInner"})
            work+=p.getVariableDataByName<Real>(name)[i];
    }
    EXPECT_NEAR(after-before,work,1e-12);
}


TEST(PDContact, FiniteRigidBoxForceIsNegativePotentialGradientAndReactionBalances)
{
    PDTestBody f; f.generate<HughesWingetSolid>(1000,1e6,0.25);
    SimpleDynamics<solid_dynamics::RigidBoxContactforPD> contact(f.body,0,0.01);
    contact.addBox("Plate",Vecd(0,-0.4,0),Vecd(0.5,0.15,0.5),Vecd(0,-0.1,0));
    auto &p=f.body.getBaseParticles();
    auto *x=p.getVariableDataByName<Vecd>("Position");
    auto *a=p.getVariableDataByName<Vecd>("PDAccelerationPrior");
    for(size_t i=0;i<p.TotalRealParticles();++i) x[i]=Vecd(0,0.5,0);
    auto evaluate=[&](Vecd point)
    {
        x[0]=point; a[0].setZero(); contact.exec();
        return p.getVariableDataByName<Real>("AnalyticContactEnergy")[0];
    };
    for(Vecd point : {Vecd(0,-0.21,0),Vecd(0.52,-0.22,0),Vecd(0,-0.27,0)})
    {
        evaluate(point);
        Vecd force=p.getVariableDataByName<Vecd>("ContactRepulsionFromWall")[0];
        Vecd reaction=p.getVariableDataByName<Vecd>("PlateReaction")[0];
        EXPECT_LT((force+reaction).norm(),1e-12);
        EXPECT_NEAR(p.getVariableDataByName<Real>("ContactBoundaryPower")[0],force.dot(Vecd(0,-0.1,0)),1e-12);
        for(int d=0;d<Dimensions;++d)
        {
            Vecd delta=Vecd::Zero(); delta[d]=1e-7;
            Real derivative=(evaluate(point+delta)-evaluate(point-delta))/(2e-7);
            EXPECT_NEAR(force[d],-derivative,1e-4*std::max(force.norm(),Real(1)));
        }
    }
    EXPECT_DOUBLE_EQ(evaluate(Vecd(0,-0.15,0)),0);
    evaluate(Vecd(0,-0.21,0)); Vecd previous=a[0];
    contact.applyAcceleration(false); contact.exec();
    EXPECT_LT((a[0]-previous).norm(),1e-12);
}


TEST(PDFSI, PressureTractionAndDensityFluxExchangeEnergyAndMomentumAtZeroExtrapolation)
{
    PDTestBody f; f.generate<HughesWingetSolid>(1000,1e6,0.25);
    auto &sp=f.body.getBaseParticles();
    sp.registerStateVariableData<Vecd>("AverageVelocity");
    sp.registerStateVariableData<Vecd>("AverageAcceleration");
    FluidBody fluid(f.system,makeShared<GeometricShapeBox>(Transform(Vecd(0,-0.6,0)),Vecd(0.3,0.3,0.3),"Fluid"));
    fluid.defineMatterMaterial<WeaklyCompressibleFluid>(935,25);
    fluid.generateParticles<BaseParticles,Lattice>();
    fluid.getBaseParticles().registerStateVariableData<Real>("VelocityDivergence");
    ContactRelationSelectKernel walls(fluid,{&f.body},true),reverse(f.body,{&fluid},false);
    fluid_dynamics::PressureWallTermsforPD fluid_force(walls,0);
    fluid_dynamics::Integration2ndHalf<Contact<Wall>,NoRiemannSolver> flux(walls);
    InteractionDynamics<solid_dynamics::PressureForceFromFluidforPD<NoRiemannSolver>> shell_force(reverse);
    shell_force.setAccelerationPressureFactor(0);
    f.system.initializeSystemCellLinkedLists(); f.system.initializeSystemConfigurations();
    auto &fp=fluid.getBaseParticles();
    Vecd vf(0.3,0.2,-0.1),vs(-0.2,-0.1,0.1);
    for(size_t i=0;i<fp.TotalRealParticles();++i)
    {
        fp.getVariableDataByName<Vecd>("Velocity")[i]=vf;
        fp.getVariableDataByName<Real>("Pressure")[i]=1000;
        fp.getVariableDataByName<Vecd>("ForcePrior")[i]=Vecd(0,30,0);
    }
    for(size_t i=0;i<sp.TotalRealParticles();++i)
    {
        sp.getVariableDataByName<Vecd>("AverageVelocity")[i]=vs;
        sp.getVariableDataByName<Vecd>("AverageAcceleration")[i]=Vecd(0,-10,0);
        sp.getVariableDataByName<Vecd>("NormalDirection")[i]=Vecd(0,-1,0);
    }
    for(size_t i=0;i<fp.TotalRealParticles();++i) fluid_force.interaction(i);
    shell_force.exec();
    for(size_t i=0;i<fp.TotalRealParticles();++i) flux.interaction(i);
    Vecd resultant=Vecd::Zero();Real power=0,scale=0;
    for(size_t i=0;i<fp.TotalRealParticles();++i)
    {
        Vecd force=fp.getVariableDataByName<Vecd>("Force")[i];
        Real rho=fp.getVariableDataByName<Real>("Density")[i];
        Real compression=fp.getVariableDataByName<Real>("Mass")[i]*1000/(rho*rho)*fp.getVariableDataByName<Real>("DensityChangeRate")[i];
        resultant+=force; power+=force.dot(vf)+compression; scale+=force.norm();
    }
    for(size_t i=0;i<sp.TotalRealParticles();++i)
    {
        Vecd force=sp.getVariableDataByName<Vecd>("PressureForceFromFluid")[i];
        resultant+=force;power+=force.dot(vs);
    }
    ASSERT_GT(scale,1);
    EXPECT_LT(resultant.norm(),1e-11*scale);
    EXPECT_LT(std::abs(power),1e-11*scale);
}


class VtpProbe : public BodyStatesRecordingToVtp
{
  public:
    using BodyStatesRecordingToVtp::BodyStatesRecordingToVtp;
    using BodyStatesRecordingToVtp::writeParticlesToVtk;
};
static std::pair<std::string,StdVec<Real>> vtkField(const std::string &xml,const std::string &name)
{
    size_t name_at=xml.find("Name=\""+name+"\"");
    if(name_at==std::string::npos) throw std::runtime_error("missing VTK field "+name);
    size_t begin=xml.rfind("<DataArray",name_at),end=xml.find('>',name_at),close=xml.find("</DataArray>",end);
    std::istringstream values(xml.substr(end+1,close-end-1));
    StdVec<Real> data; Real value;
    while(values>>value) data.push_back(value);
    return {xml.substr(begin,end-begin),data};
}
TEST(PDOutput, FullPrecisionRoundTripsSmallScalarVectorMatrixPositionAndTime)
{
    PDTestBody f(1e-3); f.generate<HughesWingetSolid>(1000,1e6,0.25);
    auto &p=f.body.getBaseParticles();
    auto *scalar=p.registerStateVariableData<Real>("TinyScalar"); scalar[0]=9.350000000000002e-13;
    auto *vector=p.registerStateVariableData<Vecd>("TinyVector"); vector[0]=Vecd(1.234567890123456e-19,-5.43210987654321e-24,2.345678901234567e-11);
    auto *matrix=p.registerStateVariableData<Matd>("TinyMatrix"); matrix[0]=Matd::Zero(); matrix[0](0,1)=-1.234567890123456e-25;
    p.addVariableToWrite<Real>("TinyScalar"); p.addVariableToWrite<Vecd>("TinyVector"); p.addVariableToWrite<Matd>("TinyMatrix");
    VtpProbe legacy(f.body); std::ostringstream rounded; legacy.writeParticlesToVtk(rounded,p);
    EXPECT_EQ(vtkField(rounded.str(),"TinyScalar").second[0],0);
    auto *position=p.getVariableDataByName<Vecd>("Position"); position[0][0]=1.234567890123456e-9;
    Real &time=*f.system.getSystemVariableDataByName<Real>("PhysicalTime"); time=1.234567890123456e-11;
    f.system.setStateRecording(true); f.body.setNewlyUpdated();
    VtpProbe precise(f.body); precise.setFullPrecisionOutput(true); precise.writeToFile(9876543);
    std::ifstream stream("output/TestSolid_ite_0009876543.vtp");
    ASSERT_TRUE(stream.good()); std::string xml((std::istreambuf_iterator<char>(stream)),{});
    for(const auto &name : {"TinyScalar","TinyVector","TinyMatrix","Position","TimeValue"})
        EXPECT_NE(vtkField(xml,name).first.find("Float64"),std::string::npos);
    auto sv=vtkField(xml,"TinyScalar").second; ASSERT_GT(sv.size(),0); EXPECT_DOUBLE_EQ(sv[0],scalar[0]);
    auto vv=vtkField(xml,"TinyVector").second; ASSERT_GE(vv.size(),3);
    for(int d=0;d<3;++d) EXPECT_DOUBLE_EQ(vv[d],vector[0][d]);
    auto mv=vtkField(xml,"TinyMatrix").second; ASSERT_GE(mv.size(),9); EXPECT_DOUBLE_EQ(mv[3],matrix[0](0,1));
    EXPECT_DOUBLE_EQ(vtkField(xml,"Position").second[0],position[0][0]);
    EXPECT_DOUBLE_EQ(vtkField(xml,"TimeValue").second[0],time);
}


TEST(PDFSI, ViscosityPairsBalanceDissipateAndIgnoreLaggedAverageVelocity)
{
    PDTestBody f; f.generate<HughesWingetSolid>(1000,1e6,0.25);
    FluidBody fluid(f.system,makeShared<GeometricShapeBox>(Transform(Vecd(0,-0.6,0)),Vecd(0.3,0.3,0.3),"ViscousFluid"));
    fluid.defineMatterMaterial<WeaklyCompressibleFluid>(935,25);
    fluid.addMaterialProperty<Viscosity>(0.25); fluid.generateParticles<BaseParticles,Lattice>();
    fluid.getBaseParticles().registerStateVariableData<Vecd>("Velocity");
    ContactRelationSelectKernel walls(fluid,{&f.body},true),reverse(f.body,{&fluid},false);
    InteractionWithUpdate<fluid_dynamics::ViscousWallTermsforPD> fluid_force(walls);
    InteractionWithUpdate<solid_dynamics::ViscousForceFromFluidforPD> shell_force(reverse);
    f.system.initializeSystemCellLinkedLists(); f.system.initializeSystemConfigurations();
    auto &fp=fluid.getBaseParticles();auto &sp=f.body.getBaseParticles();
    Vecd vf(0.3,0.2,-0.1),vs(-0.2,-0.1,0.1);
    for(size_t i=0;i<fp.TotalRealParticles();++i) fp.getVariableDataByName<Vecd>("Velocity")[i]=vf;
    for(size_t i=0;i<sp.TotalRealParticles();++i)
    {
        sp.getVariableDataByName<Vecd>("Velocity")[i]=vs;
        sp.getVariableDataByName<Vecd>("AverageVelocity")[i]=Vecd(10,-30,50);
    }
    fluid_force.exec();shell_force.exec(); Vecd net=Vecd::Zero();Real power=0,scale=0;
    for(size_t i=0;i<fp.TotalRealParticles();++i)
    {Vecd force=fp.getVariableDataByName<Vecd>("ViscousForce")[i];net+=force;power+=force.dot(vf);scale+=force.norm();}
    for(size_t i=0;i<sp.TotalRealParticles();++i)
    {Vecd force=sp.getVariableDataByName<Vecd>("ViscousForceFromFluid")[i];net+=force;power+=force.dot(vs);}
    ASSERT_GT(scale,1e-8);EXPECT_LT(net.norm(),1e-11*scale);EXPECT_LT(power,0);
    // A common instantaneous translation must have no interface viscosity,
    // even when the displacement-averaged wall velocity is one half step old.
    for(size_t i=0;i<fp.TotalRealParticles();++i)
    {
        fp.getVariableDataByName<Vecd>("Velocity")[i]=vs;
        fp.getVariableDataByName<Vecd>("ViscousForce")[i].setZero();
    }
    fluid_force.exec();shell_force.exec();
    for(size_t i=0;i<fp.TotalRealParticles();++i) EXPECT_LT(fp.getVariableDataByName<Vecd>("ViscousForce")[i].norm(),1e-15);
    for(size_t i=0;i<sp.TotalRealParticles();++i) EXPECT_LT(sp.getVariableDataByName<Vecd>("ViscousForceFromFluid")[i].norm(),1e-15);
}
