#include "sphinxsys.h"
#include <gtest/gtest.h>
#include <limits>

using namespace SPH;

class PDTestBody
{
  public:
    SPHSystem system{BoundingBoxd(Vecd(-1, -1, -1), Vecd(1, 1, 1)), 0.1};
    PDBody body{system, makeShared<GeometricShapeBox>(Transform(), Vecd(0.3, 0.3, 0.3), "TestSolid")};
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
