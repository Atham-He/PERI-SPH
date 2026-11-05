#include "sphinxsys.h"
#include <gtest/gtest.h>
#include <limits>

using namespace SPH;

namespace
{
constexpr Real E = 1e6, nu = 0.25, yield = 1000;

/** A real material/particle registry, without contacts, fracture or PD forces. */
class MaterialPoint
{
  public:
    SPHSystem system;
    PDBody body;
    MaterialPoint(const PDIsotropicHardening &law, bool viscous = false, Real eta = 0,
                  Real kinematic = 0, Real spacing = 0.1)
        : system(BoundingBoxd(Vecd(-1, -1, -1), Vecd(1, 1, 1)), spacing),
          body(system, makeShared<GeometricShapeBox>(Transform(), Vecd(0.15, 0.15, 0.15), "MaterialPoint"))
    {
        if (viscous)
            body.defineMatterMaterial<ViscoPlasticSolidforPD>(1000, E, nu, yield, eta, law, kinematic);
        else
            body.defineMatterMaterial<J2PlasticityforPD>(1000, E, nu, yield, law, kinematic);
        body.generateParticles<BaseParticles, Lattice>();
    }
    J2PlasticityforPD &material() { return DynamicCast<J2PlasticityforPD>(&body, body.getMatterMaterial()); }
    BaseParticles &particles() { return body.getBaseParticles(); }
    Matd &stress() { return particles().getVariableDataByName<Matd>("PDStress")[0]; }
    Matd &plastic() { return particles().getVariableDataByName<Matd>("PlasticStrain")[0]; }
    Real alpha() { return particles().getVariableDataByName<Real>("IsotropicHardeningParam")[0]; }
    void step(const Matd &increment, Real dt = 0.01)
    {
        stress() = material().PlasticConstitutiveRelation(increment, stress(), plastic(), 0, dt);
    }
};

// Independent scalar uniaxial-stress backward-Euler reference. It solves for
// accumulated plastic strain using axial compatibility, not the tensor code's
// return residual/YieldFunc. Monotonic tension and tensile holds only.
struct UniaxialReference
{
    long double strain = 0, lateral = 0, alpha = 0, stress = 0;
    Matd advance(long double next_strain, long double dt, long double H = 0,
                 long double Q = 0, long double b = 0, long double eta = 0)
    {
        const long double old_alpha = alpha;
        auto resistance = [&](long double a)
        { return yield + H * a - Q * std::expm1(-b * a) + eta * (a - old_alpha) / dt; };
        const long double trial = E * (next_strain - old_alpha);
        if (trial > resistance(old_alpha))
        {
            long double lo = old_alpha, hi = next_strain;
            for (int k = 0; k < 100; ++k)
            {
                const long double mid = (lo + hi) / 2;
                if (E * (next_strain - mid) > resistance(mid)) lo = mid;
                else hi = mid;
            }
            alpha = (lo + hi) / 2;
        }
        stress = E * (next_strain - alpha);
        const long double next_lateral = -nu * stress / E - alpha / 2;
        Matd increment = Matd::Zero();
        increment(0, 0) = Real(next_strain - strain);
        increment(1, 1) = increment(2, 2) = Real(next_lateral - lateral);
        strain = next_strain;
        lateral = next_lateral;
        return increment;
    }
};

void expectUniaxial(MaterialPoint &point, const UniaxialReference &reference)
{
    EXPECT_NEAR(point.stress()(0, 0), Real(reference.stress), 2e-6);
    EXPECT_NEAR(point.stress()(1, 1), 0, 2e-6);
    EXPECT_NEAR(point.stress()(2, 2), 0, 2e-6);
    EXPECT_NEAR(point.alpha(), Real(reference.alpha), 2e-12);
    EXPECT_NEAR(point.plastic()(0, 0), Real(reference.alpha), 2e-12);
    EXPECT_NEAR(point.plastic().trace(), 0, 1e-13);
}
} // namespace

TEST(PDHardening, ExplicitLawsExposeParametersAndRetainSmallSaturationIncrements)
{
    const auto perfect = PDIsotropicHardening::Perfect();
    EXPECT_EQ(perfect.law(), PDIsotropicHardening::Law::Perfect);
    EXPECT_DOUBLE_EQ(perfect.value(0.7), 0);
    const auto linear = PDIsotropicHardening::Linear(2000);
    EXPECT_DOUBLE_EQ(linear.value(0.5), 1000);
    const auto saturation = PDIsotropicHardening::Saturation(2000, 300, 4);
    EXPECT_DOUBLE_EQ(saturation.linearModulus(), 2000);
    EXPECT_DOUBLE_EQ(saturation.saturationIncrement(), 300);
    EXPECT_DOUBLE_EQ(saturation.saturationRate(), 4);
    const auto tiny = PDIsotropicHardening::Saturation(0, 300, 4);
    EXPECT_NEAR(tiny.value(1e-20) / 1e-20, 1200, 1e-10);
}

TEST(PDHardening, InvalidParametersAndHistoryReject)
{
    const Real nan = std::numeric_limits<Real>::quiet_NaN();
    const Real inf = std::numeric_limits<Real>::infinity();
    EXPECT_THROW(PDIsotropicHardening::Linear(-1), std::invalid_argument);
    EXPECT_THROW(PDIsotropicHardening::Saturation(0, -1, 1), std::invalid_argument);
    EXPECT_THROW(PDIsotropicHardening::Saturation(0, 1, -1), std::invalid_argument);
    EXPECT_THROW(PDIsotropicHardening::Saturation(nan, 1, 1), std::invalid_argument);
    EXPECT_THROW(PDIsotropicHardening::Saturation(0, inf, 1), std::invalid_argument);
    EXPECT_THROW(PDIsotropicHardening::Legacy(0), std::invalid_argument);
    EXPECT_THROW(PDIsotropicHardening::Perfect().value(-1), std::invalid_argument);
    EXPECT_THROW(PDIsotropicHardening::Perfect().value(nan), std::invalid_argument);
}

TEST(PDJ2Material, LegacyConstructorPreservesOriginalHardeningFormula)
{
    J2PlasticityforPD legacy(1000, E, nu, yield, 2000);
    EXPECT_EQ(legacy.HardeningParameters().law(), PDIsotropicHardening::Law::Legacy);
    for (Real alpha : {0.0, 1e-14, 0.01, 1.0})
        EXPECT_DOUBLE_EQ(legacy.IsoHardeningFunc(alpha),
            2000 * alpha + (1.589 * yield - yield) * (1 - std::exp(-16.93 * alpha)));
    EXPECT_GT(legacy.IsoHardeningFunc(0.01), 0);
}

TEST(PDJ2Material, PerfectPlasticUniaxialPlateauElasticUnloadAndResidual)
{
    MaterialPoint point(PDIsotropicHardening::Perfect());
    UniaxialReference ref;
    for (int k = 1; k <= 40; ++k)
    {
        point.step(ref.advance(0.0002L * k, 0.01L));
        expectUniaxial(point, ref);
        EXPECT_NEAR(ref.stress, std::min(Real(E * ref.strain), yield), 1e-10);
    }
    const Real residual = point.alpha();
    Real previous_stress = point.stress()(0, 0);
    const Real decrement = 0.0001;
    for (int k = 0; k < 10; ++k)
    {
        Matd increment = Matd::Zero();
        increment(0, 0) = -decrement;
        increment(1, 1) = increment(2, 2) = nu * decrement;
        point.step(increment);
        EXPECT_NEAR(point.stress()(0, 0) - previous_stress, -E * decrement, 1e-7);
        EXPECT_NEAR(point.alpha(), residual, 2e-12);
        previous_stress = point.stress()(0, 0);
    }
    EXPECT_NEAR(point.stress().norm(), 0, 2e-6);
    EXPECT_NEAR(Real(ref.strain) - yield / E, residual, 2e-12);
}

TEST(PDJ2Material, LinearHardeningMatchesUniaxialAnalyticTangent)
{
    const Real H = 20000;
    MaterialPoint point(PDIsotropicHardening::Linear(H));
    UniaxialReference ref;
    Real previous = 0;
    for (int k = 1; k <= 40; ++k)
    {
        point.step(ref.advance(0.0002L * k, 0.01L, H));
        expectUniaxial(point, ref);
        if (k > 6)
            EXPECT_NEAR((point.stress()(0, 0) - previous) / 0.0002, E * H / (E + H), 0.01);
        previous = point.stress()(0, 0);
    }
}

TEST(PDJ2Material, SaturationHardeningMatchesIndependentScalarPath)
{
    MaterialPoint point(PDIsotropicHardening::Saturation(20000, 700, 80));
    UniaxialReference ref;
    for (int k = 1; k <= 60; ++k)
    {
        point.step(ref.advance(0.0002L * k, 0.01L, 20000, 700, 80));
        expectUniaxial(point, ref);
    }
}

TEST(PDJ2Material, PureShearYieldAndHydrostaticNoFlow)
{
    MaterialPoint shear(PDIsotropicHardening::Perfect());
    Matd increment = Matd::Zero(); increment(0, 1) = increment(1, 0) = 0.005;
    shear.step(increment);
    EXPECT_NEAR(shear.stress()(0, 1), yield / std::sqrt(3.0), 1e-6);
    EXPECT_GT(shear.alpha(), 0);
    MaterialPoint hydro(PDIsotropicHardening::Perfect());
    hydro.step(0.01 * Matd::Identity());
    EXPECT_DOUBLE_EQ(hydro.alpha(), 0);
    EXPECT_DOUBLE_EQ(hydro.plastic().norm(), 0);
    EXPECT_NEAR(hydro.stress()(0, 0), E / (1 - 2 * nu) * 0.01, 1e-8);
}

TEST(PDJ2Material, HardeningZeroLimitsAndRateIndependentPhysicalDuration)
{
    MaterialPoint perfect(PDIsotropicHardening::Perfect());
    MaterialPoint zero(PDIsotropicHardening::Saturation(0, 0, 20));
    MaterialPoint linear(PDIsotropicHardening::Linear(2000));
    MaterialPoint zero_rate(PDIsotropicHardening::Saturation(2000, 700, 0));
    Matd increment = Matd::Zero(); increment(0, 0) = 0.0004;
    increment(1, 1) = increment(2, 2) = -0.0002;
    for (int k = 0; k < 40; ++k)
    {
        perfect.step(increment, 0.001); zero.step(increment, 10);
        linear.step(increment, 0.001); zero_rate.step(increment, 10);
        EXPECT_LT((perfect.stress() - zero.stress()).norm(), 1e-7);
        EXPECT_LT((linear.stress() - zero_rate.stress()).norm(), 1e-7);
        EXPECT_NEAR(perfect.alpha(), zero.alpha(), 1e-12);
    }
}

TEST(PDJ2Material, FiniteRotationCarriesPlasticAndKinematicHistoryWithoutFlow)
{
    MaterialPoint point(PDIsotropicHardening::Linear(2000), false, 0, 3000);
    Matd increment = Matd::Zero(); increment(0, 0) = 0.01;
    increment(1, 1) = increment(2, 2) = -0.005;
    point.step(increment);
    const Matd stress = point.stress(), plastic = point.plastic();
    auto *back = point.particles().getVariableDataByName<Matd>("KinematicHardeningParam");
    const Matd back_old = back[0];
    const Real alpha = point.alpha(), angle = 0.43;
    Matd R = Matd::Identity();
    R(0, 0) = R(1, 1) = std::cos(angle); R(0, 1) = -std::sin(angle); R(1, 0) = std::sin(angle);
    Matd spin = Matd::Zero(); spin(0, 1) = -2 * std::tan(angle / 2); spin(1, 0) = -spin(0, 1);
    point.step(spin);
    EXPECT_LT((point.stress() - R * stress * R.transpose()).norm(), 1e-6);
    EXPECT_LT((point.plastic() - R * plastic * R.transpose()).norm(), 1e-12);
    EXPECT_LT((back[0] - R * back_old * R.transpose()).norm(), 1e-12);
    EXPECT_NEAR(point.alpha(), alpha, 1e-12);
}

TEST(PDViscoplasticMaterial, ZeroViscosityRecoversSameConfiguredJ2)
{
    const auto law = PDIsotropicHardening::Saturation(2000, 700, 80);
    MaterialPoint j2(law), viscous(law, true, 0);
    Matd increment = Matd::Zero(); increment(0, 0) = 0.0004;
    increment(1, 1) = increment(2, 2) = -0.0002;
    for (int k = 0; k < 40; ++k)
    {
        j2.step(increment); viscous.step(increment);
        EXPECT_LT((j2.stress() - viscous.stress()).norm(), 1e-10);
        EXPECT_NEAR(j2.alpha(), viscous.alpha(), 1e-14);
    }
}

TEST(PDViscoplasticMaterial, SaturationLoadAndHoldMatchIndependentBackwardEuler)
{
    MaterialPoint point(PDIsotropicHardening::Saturation(2000, 700, 80), true, 20000);
    UniaxialReference ref;
    for (int k = 1; k <= 80; ++k)
    {
        const long double strain = 0.0002L * std::min(k, 40);
        point.step(ref.advance(strain, 0.01L, 2000, 700, 80, 20000));
        expectUniaxial(point, ref);
    }
}

TEST(PDViscoplasticMaterial, PhysicalRateAndTimeRefinementMatchContinuousUniaxialSolution)
{
    const Real eta = 20000, end = 0.4;
    for (Real rate : {0.01, 0.02})
    {
        const Real exact = yield + eta * rate * (1 - std::exp(-(end - yield / (E * rate)) * E / eta));
        Real previous_error = std::numeric_limits<Real>::infinity();
        for (int steps : {20, 40, 80})
        {
            MaterialPoint point(PDIsotropicHardening::Perfect(), true, eta);
            UniaxialReference ref;
            const long double dt = static_cast<long double>(end) / steps;
            for (int k = 1; k <= steps; ++k)
            {
                point.step(ref.advance(rate * dt * k, dt, 0, 0, 0, eta), Real(dt));
                expectUniaxial(point, ref);
            }
            const Real error = std::abs(point.stress()(0, 0) - exact);
            EXPECT_LT(error, previous_error);
            previous_error = error;
        }
        EXPECT_LT(previous_error, 0.001);
    }
}

TEST(PDPlasticPatch, AffineDeviatoricIncrementMatchesAnalyticHardeningAtEveryParticle)
{
    // Exercise the actual PD shape reconstruction and Plastic SecondStep, not
    // just the point material. Cayley geometry gives exactly G=diag(e,-e/2,-e/2).
    const Real e = 0.01, dt = 0.01, shear = E / (2 * (1 + nu));
    for (Real H : {0.0, 20000.0})
        for (Real spacing : {0.1, 0.075, 0.05})
        {
            const auto law = H == 0 ? PDIsotropicHardening::Perfect() : PDIsotropicHardening::Linear(H);
            MaterialPoint point(law, false, 0, 0, spacing);
            InnerRelation reference(point.body);
            point.system.initializeSystemCellLinkedLists();
            point.system.initializeSystemConfigurations();
            InteractionWithUpdate<solid_dynamics::NosbPDShapeMatrix> shape(reference);
            InteractionWithUpdate<solid_dynamics::NosbPDSecondStepPlastic> second(reference);
            shape.exec();
            Matd D = Matd::Zero(); D(0, 0) = e; D(1, 1) = D(2, 2) = -e / 2;
            const Matd F = (Matd::Identity() + 0.5 * D) * (Matd::Identity() - 0.5 * D).inverse();
            auto &particles = point.particles();
            auto *position = particles.getVariableDataByName<Vecd>("Position");
            auto *velocity = particles.getVariableDataByName<Vecd>("Velocity");
            for (size_t i = 0; i < particles.TotalRealParticles(); ++i)
            {
                const Vecd initial = position[i];
                position[i] = F * initial;
                velocity[i] = (F - Matd::Identity()) * initial / dt;
            }
            second.exec(dt);
            const Real alpha = (3 * shear * e - yield) / (3 * shear + H);
            const Real equivalent = yield + H * alpha;
            const Matd expected = Vecd(2 * equivalent / 3, -equivalent / 3, -equivalent / 3).asDiagonal();
            const auto *stress = particles.getVariableDataByName<Matd>("PDStress");
            const auto *deformation = particles.getVariableDataByName<Matd>("DeformationGradient");
            const auto *history = particles.getVariableDataByName<Real>("IsotropicHardeningParam");
            ASSERT_GT(particles.TotalRealParticles(), 0);
            for (size_t i = 0; i < particles.TotalRealParticles(); ++i)
            {
                EXPECT_LT((deformation[i] - F).norm(), 1e-10);
                EXPECT_LT((stress[i] - expected).norm(), 2e-6);
                EXPECT_NEAR(history[i], alpha, 2e-12);
            }
        }
}
