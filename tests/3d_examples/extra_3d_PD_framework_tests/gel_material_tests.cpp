#include "sphinxsys.h"
#include "pd_checkpoint.h"
#include "pd_gel_material.h"
#include "pd_gel_dynamics.h"
#include <gtest/gtest.h>
#include <filesystem>
#include <limits>

using namespace SPH;
using namespace SPH::peridynamics;

namespace
{
const GelParameters parameters = GelParameters::fromYoungPoisson(1e6, 0.35, 1.2, 0.7);
Matd deformation()
{
    Matd F = Matd::Identity();
    F(0, 0) = 1.17; F(1, 1) = 0.91; F(0, 1) = 0.08;
    return F;
}

Matd rotation()
{
    Matd R = Matd::Identity();
    const Real angle = 0.43;
    R(0, 0) = R(1, 1) = std::cos(angle);
    R(0, 1) = -std::sin(angle); R(1, 0) = std::sin(angle);
    return R;
}
} // namespace

TEST(PDGelMaterial, FiniteDeformationStressIsFreeEnergyGradientAtFixedHistory)
{
    const Matd F = deformation();
    GelHistory history;
    history.viscous = 0.2 * gelDeviatoricGreen(F);
    const auto response = gelResponse(F, parameters, history);
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
        {
            Matd delta = Matd::Zero(); delta(i, j) = 1e-7;
            const Real derivative = (gelResponse(F + delta, parameters, history).free_energy -
                                     gelResponse(F - delta, parameters, history).free_energy) / 2e-7;
            EXPECT_NEAR(response.pk1(i, j), derivative, 1e-2);
        }
}

TEST(PDGelMaterial, SuperposedFiniteRotationRotatesStressAndPreservesMemoryEnergy)
{
    const Matd F = deformation(), R = rotation();
    const auto a = advanceGel(F, 0.13, parameters, GelHistory{});
    const auto b = advanceGel(R * F, 0.13, parameters, GelHistory{});
    const auto before = gelResponse(F, parameters, a), after = gelResponse(R * F, parameters, b);
    EXPECT_LT((a.viscous - b.viscous).norm(), 1e-15);
    EXPECT_NEAR(a.dissipation, b.dissipation, 1e-9);
    EXPECT_NEAR(before.free_energy, after.free_energy, 1e-9);
    EXPECT_LT((after.pk1 - R * before.pk1).norm(), 1e-8);
    EXPECT_LT((after.cauchy - R * before.cauchy * R.transpose()).norm(), 1e-8);
    EXPECT_LT(gelResponse(R, parameters, GelHistory{}).pk1.norm(), 1e-8);
}

TEST(PDGelMaterial, FixedStrainRelaxationMatchesExponentialAndDissipationBalancesEnergy)
{
    const Matd F = deformation();
    GelHistory held;
    held.green = gelDeviatoricGreen(F); // State immediately after an ideal instantaneous strain step.
    const auto initial = gelResponse(F, parameters, held);
    const Real dt = 0.83;
    const auto relaxed = advanceGel(F, dt, parameters, held);
    EXPECT_LT((held.green - relaxed.viscous - std::exp(-dt / parameters.relaxation_time) * held.green).norm(), 1e-15);
    const auto final = gelResponse(F, parameters, relaxed);
    EXPECT_NEAR(initial.free_energy - final.free_energy, relaxed.dissipation, 1e-9);
    EXPECT_GT(relaxed.dissipation, 0);
}

TEST(PDGelMaterial, ExactLinearGreenStrainSegmentsAreIndependentOfSubdivision)
{
    const Real duration = 0.6;
    auto F = [&](Real fraction)
    {
        Matd value = Matd::Identity();
        value(0, 0) = std::sqrt(1 + 0.4 * fraction);
        value(1, 1) = std::sqrt(1 - 0.12 * fraction);
        return value;
    };
    const auto one = advanceGel(F(1), duration, parameters, GelHistory{});
    for (int steps : {2, 7, 64})
    {
        GelHistory many;
        for (int k = 1; k <= steps; ++k)
            many = advanceGel(F(Real(k) / steps), duration / steps, parameters, many);
        EXPECT_LT((many.viscous - one.viscous).norm(), 1e-14);
        EXPECT_NEAR(many.dissipation, one.dissipation, 1e-8);
    }
}

TEST(PDGelMaterial, NonlinearHistoryShowsSecondOrderTimeRefinementAgainstAnalyticSolution)
{
    // dev(E)=A*t^2, q=A*(2*tau*t-2*tau^2*(1-exp(-t/tau))).
    Matd A = Matd::Zero(); A(0, 0) = 0.1; A(1, 1) = -0.05; A(2, 2) = -0.05;
    const Real duration = 0.6, tau = parameters.relaxation_time;
    const Matd q_exact = A * (2 * tau * duration - 2 * tau * tau * (1 - std::exp(-duration / tau)));
    StdVec<Real> errors;
    for (int steps : {8, 16, 32})
    {
        GelHistory state;
        for (int k = 1; k <= steps; ++k)
        {
            const Real t = duration * k / steps;
            Matd F = Matd::Identity();
            for (int d = 0; d < 3; ++d) F(d, d) = std::sqrt(1 + 2 * A(d, d) * t * t);
            state = advanceGel(F, duration / steps, parameters, state);
        }
        errors.push_back((state.green - state.viscous - q_exact).norm());
    }
    EXPECT_GT(errors[0] / errors[1], 3.9);
    EXPECT_GT(errors[1] / errors[2], 3.9);
    EXPECT_LT(errors[1] / errors[2], 4.1);
}

TEST(PDGelMaterial, StableTinyAndLongRelaxationIncrementsDoNotClipDissipation)
{
    for (Real ratio : {1e-10, 1e-6, 1e-3, 1., 1000.})
    {
        const auto state = advanceGel(deformation(), ratio * parameters.relaxation_time, parameters, GelHistory{});
        EXPECT_GE(state.dissipation, 0);
        EXPECT_TRUE(state.viscous.allFinite());
    }
    const auto unchanged = advanceGel(deformation(), 0, parameters, GelHistory{});
    EXPECT_EQ(unchanged.green.squaredNorm(), 0);
    EXPECT_EQ(unchanged.viscous.squaredNorm(), 0);
    EXPECT_EQ(unchanged.dissipation, 0);
}

TEST(PDGelMaterial, InversionBadMemoryAndInvalidParametersAreRejected)
{
    Matd inverted = Matd::Identity(); inverted(0, 0) = -1;
    EXPECT_THROW(gelResponse(inverted, parameters, GelHistory{}), std::runtime_error);
    GelHistory bad; bad.viscous(0, 0) = std::numeric_limits<Real>::quiet_NaN();
    EXPECT_THROW(advanceGel(deformation(), 0.1, parameters, bad), std::runtime_error);
    EXPECT_THROW(GelParameters::fromYoungPoisson(1e6, 0.5, 1, 1), std::invalid_argument);
    EXPECT_THROW(GelParameters::fromYoungPoisson(1e6, 0.3, 1, 0), std::invalid_argument);
    EXPECT_THROW(advanceGel(deformation(), -1, parameters, GelHistory{}), std::invalid_argument);
}

TEST(PDGelMaterial, MidpointPredictionLeavesAcceptedMemoryUnchangedAndCommitMatchesPureKernel)
{
    SPHSystem system(BoundingBoxd(Vecd(-1, -1, -1), Vecd(1, 1, 1)), 0.1);
    PDBody body(system, makeShared<GeometricShapeBox>(Transform(), Vecd(0.3, 0.3, 0.3), "GelPatch"));
    body.defineMatterMaterial<GelViscoelasticSolid>(1000, 1e6, 0.35, 1.2, 0.7);
    body.generateParticles<CheckpointParticles, Lattice>();
    InnerRelation reference(body);
    InteractionWithUpdate<solid_dynamics::NosbPDShapeMatrix> shape(reference);
    using Update = solid_dynamics::GelConstitutiveUpdate;
    InteractionDynamics<Update> predict(reference, Update::Stage::MidpointForce), commit(reference, Update::Stage::AcceptedHistory);
    system.initializeSystemCellLinkedLists(); system.initializeSystemConfigurations(); shape.exec();
    auto &p = body.getBaseParticles();
    auto *x = p.getVariableDataByName<Vecd>("Position");
    StdVec<Vecd> initial(x, x + p.TotalRealParticles());
    const Matd Fmid = deformation();
    Matd Fend = Fmid; Fend(0, 0) += 0.02;
    for (size_t i = 0; i < p.TotalRealParticles(); ++i) x[i] = Fmid * initial[i];
    predict.exec(0.2);
    for (size_t i = 0; i < p.TotalRealParticles(); ++i)
    {
        EXPECT_EQ(p.getVariableDataByName<Matd>("GelViscousStrain")[i].squaredNorm(), 0);
        EXPECT_EQ(p.getVariableDataByName<Real>("GelHistoryTime")[i], 0);
        x[i] = Fend * initial[i];
    }
    commit.exec(0.2);
    const auto expected = advanceGel(Fend, 0.1, parameters, advanceGel(Fmid, 0.1, parameters, GelHistory{}));
    for (size_t i = 0; i < p.TotalRealParticles(); ++i)
    {
        EXPECT_LT((p.getVariableDataByName<Matd>("GelViscousStrain")[i] - expected.viscous).norm(), 1e-13);
        EXPECT_NEAR(p.getVariableDataByName<Real>("GelDissipationDensity")[i], expected.dissipation, 1e-7);
        EXPECT_DOUBLE_EQ(p.getVariableDataByName<Real>("GelHistoryTime")[i], 0.2);
        EXPECT_LT((p.getVariableDataByName<Matd>("DeformationGradient")[i] - Fend).norm(), 1e-13);
    }
}

TEST(PDGelMaterial, AffineFiniteDeformationIsReproducedAtBoundaryAcrossThreeParticleSpacings)
{
    const Matd imposed = deformation();
    for (Real spacing : {0.2, 0.15, 0.1})
    {
        SPHSystem system(BoundingBoxd(Vecd(-1, -1, -1), Vecd(1, 1, 1)), spacing);
        PDBody body(system, makeShared<GeometricShapeBox>(Transform(), Vecd(0.3, 0.3, 0.3), "GelSpatialPatch"));
        body.defineMatterMaterial<GelViscoelasticSolid>(1000, 1e6, 0.35, 1.2, 0.7);
        body.generateParticles<BaseParticles, Lattice>();
        InnerRelation reference(body);
        InteractionWithUpdate<solid_dynamics::NosbPDShapeMatrix> shape(reference);
        InteractionDynamics<solid_dynamics::GelConstitutiveUpdate> predict(
            reference, solid_dynamics::GelConstitutiveUpdate::Stage::MidpointForce);
        system.initializeSystemCellLinkedLists(); system.initializeSystemConfigurations(); shape.exec();
        auto &p = body.getBaseParticles();
        auto *x = p.getVariableDataByName<Vecd>("Position");
        for (size_t i = 0; i < p.TotalRealParticles(); ++i) x[i] = imposed * x[i] + Vecd(0.1, -0.07, 0.03);
        predict.exec(0.2);
        const auto expected = gelResponse(imposed, parameters, advanceGel(imposed, 0.1, parameters, GelHistory{}));
        for (size_t i = 0; i < p.TotalRealParticles(); ++i)
        {
            EXPECT_LT((p.getVariableDataByName<Matd>("DeformationGradient")[i] - imposed).norm(), 1e-12);
            EXPECT_LT((p.getVariableDataByName<Matd>("PDPK1")[i] - expected.pk1).norm(), 1e-6);
        }
    }
}

TEST(PDGelCheckpoint, NonzeroViscousMemoryAndDissipationRestoreWithoutAdvancing)
{
    SPHSystem system(BoundingBoxd(Vecd(-1, -1, -1), Vecd(1, 1, 1)), 0.1);
    PDBody body(system, makeShared<GeometricShapeBox>(Transform(), Vecd(0.3, 0.3, 0.3), "GelCheckpoint"));
    body.defineMatterMaterial<GelViscoelasticSolid>(1000, 1e6, 0.35, 1.2, 0.7);
    body.generateParticles<CheckpointParticles, Lattice>();
    InnerRelation reference(body);
    InteractionWithUpdate<solid_dynamics::NosbPDShapeMatrix> shape(reference);
    using Update = solid_dynamics::GelConstitutiveUpdate;
    InteractionDynamics<Update> predict(reference, Update::Stage::MidpointForce), commit(reference, Update::Stage::AcceptedHistory);
    system.initializeSystemCellLinkedLists(); system.initializeSystemConfigurations(); shape.exec();
    auto &p = body.getBaseParticles();
    auto *x = p.getVariableDataByName<Vecd>("Position");
    for (size_t i = 0; i < p.TotalRealParticles(); ++i) x[i] = deformation() * x[i];
    predict.exec(0.2); commit.exec(0.2);
    auto *viscous = p.getVariableDataByName<Matd>("GelViscousStrain");
    auto *loss = p.getVariableDataByName<Real>("GelDissipationDensity");
    const Matd saved_viscous = viscous[0]; const Real saved_loss = loss[0];
    ASSERT_GT(saved_viscous.norm(), 0); ASSERT_GT(saved_loss, 0);
    const auto root = std::filesystem::temp_directory_path() / ("gel-memory-unit-" +
        std::to_string(std::chrono::steady_clock::now().time_since_epoch().count()));
    DiskCheckpoint checkpoint(root, "nonzero-gel-history");
    checkpoint.addBody(dynamic_cast<CheckpointParticles &>(p));
    checkpoint.addConfiguration("Reference", reference.inner_configuration_, p.ParticlesBound());
    const auto first = checkpoint.save({0.2, 0, 1, 1, 0});
    viscous[0].setZero(); loss[0] = -1;
    checkpoint.restore(root);
    EXPECT_EQ((viscous[0] - saved_viscous).squaredNorm(), 0);
    EXPECT_DOUBLE_EQ(loss[0], saved_loss);
    const auto second = checkpoint.save({0.2, 0, 1, 1, 0});
    EXPECT_EQ(checkpoint_detail::fileHash(first / "state.bin"), checkpoint_detail::fileHash(second / "state.bin"));
    std::filesystem::remove_all(root);
}
