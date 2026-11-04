#include "sphinxsys.h"
#include "pd_particle_compression.h"
#include <gtest/gtest.h>
#include <limits>

using namespace SPH;

namespace
{
struct SolidSphere
{
    SPHSystem system{BoundingBoxd(Vecd(-1, -1, -1), Vecd(1, 1, 1)), 0.1};
    PDBody body{system, makeShared<GeometricShapeBall>(Vecd::Zero(), 0.3, "SolidGelControl")};

    SolidSphere()
    {
        body.defineMatterMaterial<HughesWingetSolid>(1000, 1e6, 0.25);
        body.generateParticles<BaseParticles, Lattice>();
    }
};
} // namespace

TEST(PDGelGeometry, SphereContainsItsInteriorInsteadOfAHollowShellOrDummyFluid)
{
    SolidSphere solid;
    auto &particles = solid.body.getBaseParticles();
    const auto *position = particles.getVariableDataByName<Vecd>("Position");
    Real nearest = 1;
    ASSERT_GT(particles.TotalRealParticles(), 0u);
    for (size_t i = 0; i < particles.TotalRealParticles(); ++i)
    {
        nearest = std::min(nearest, position[i].norm());
        EXPECT_LE(position[i].norm(), 0.3 + 1e-12);
    }
    EXPECT_LT(nearest, 0.1);
}

TEST(PDGelOptions, GeometryMaterialContactAndTimeRemainIndependentlyControllable)
{
    peridynamics::ParticleCompressionOptions options;
    options.validate();
    const Real duration = options.duration(), horizon = options.horizon();
    options.contact.friction_mu = 0.2;
    options.material.youngs_modulus *= 2;
    options.validate();
    EXPECT_DOUBLE_EQ(options.duration(), duration);
    EXPECT_DOUBLE_EQ(options.horizon(), horizon);
    const auto serialized = options.bindings().serialize();
    EXPECT_NE(serialized.find("friction_mu="), std::string::npos);
    EXPECT_NE(serialized.find("youngs_modulus="), std::string::npos);
    EXPECT_EQ(serialized.find("shell_thickness="), std::string::npos);
}

TEST(PDGelOptions, RejectsSingularMaterialUnresolvedSphereAndNonFiniteInputs)
{
    peridynamics::ParticleCompressionOptions options;
    options.material.poisson_ratio = 0.5;
    EXPECT_THROW(options.validate(), std::invalid_argument);
    options.material.poisson_ratio = 0.35;
    options.geometry.spacing = options.geometry.radius;
    EXPECT_THROW(options.validate(), std::invalid_argument);
    options.geometry.spacing = 2e-5; // Exactly three spacings in the radius.
    EXPECT_NO_THROW(options.validate());
    options.geometry.spacing = 1e-5;
    options.contact.friction_mu = std::numeric_limits<Real>::quiet_NaN();
    EXPECT_THROW(options.validate(), std::invalid_argument);
}

TEST(PDGelMotion, PlateGeometryVelocityAndStopUseTheSameExplicitTime)
{
    const peridynamics::CompressionPlateMotion motion{0.5, 0.1, 0.2, 2};
    for (Real time : {-1., 0., 0.7, 1.9, 2., 3.})
    {
        EXPECT_NEAR(motion.center(time)[1], 0.55 - 0.2 * std::clamp(time / 2, Real(0), Real(1)), 1e-15);
        EXPECT_DOUBLE_EQ(motion.velocity(time)[1], time >= 0 && time < 2 ? -0.1 : 0);
    }
    EXPECT_THROW((peridynamics::CompressionPlateMotion{0.5, 0.1, 0.2, 0}), std::invalid_argument);
    EXPECT_THROW(motion.center(std::numeric_limits<Real>::quiet_NaN()), std::invalid_argument);
}

TEST(PDContactObservation, AcceptedStateSamplingDoesNotOverwriteAppliedLoadsOrAccelerations)
{
    SolidSphere solid;
    auto &particles = solid.body.getBaseParticles();
    SimpleDynamics<solid_dynamics::RigidBoxContactforPD> contact(solid.body, 0.2, 0.01);
    contact.addBox("MovingPlate", Vecd(0, -0.4, 0), Vecd(0.5, 0.15, 0.5), Vecd(0, -0.1, 0));
    auto *position = particles.getVariableDataByName<Vecd>("Position");
    auto *velocity = particles.getVariableDataByName<Vecd>("Velocity");
    auto *prior = particles.getVariableDataByName<Vecd>("PDAccelerationPrior");
    for (size_t i = 0; i < particles.TotalRealParticles(); ++i)
    {
        position[i] = Vecd(0, 0.5, 0);
        velocity[i] = Vecd(0.03, 0, 0);
        prior[i] = Vecd(1, 2, 3);
    }
    position[0] = Vecd(0, -0.21, 0);
    contact.exec();
    // Move the observation point after the applied load was stored. This makes
    // an accidental re-execution detectable, including force/power arrays.
    position[0][1] += 0.005;
    const auto *applied = particles.getVariableDataByName<Vecd>("ContactRepulsionFromWall");
    const Vecd last_applied = applied[0];
    const char *vectors[] = {"Position", "Velocity", "PDAccelerationPrior", "PDAcceleration",
                            "ContactRepulsionFromWall", "ContactFrictionForce", "MovingPlateReaction"};
    StdVec<StdVec<Vecd>> before_vectors;
    for (const auto *name : vectors)
    {
        auto *values = particles.getVariableDataByName<Vecd>(name);
        before_vectors.emplace_back(values, values + particles.TotalRealParticles());
    }
    const char *scalars[] = {"AnalyticContactEnergy", "ContactBoundaryPower", "ContactFrictionBoundaryPower"};
    StdVec<StdVec<Real>> before_scalars;
    for (const auto *name : scalars)
    {
        auto *values = particles.getVariableDataByName<Real>(name);
        before_scalars.emplace_back(values, values + particles.TotalRealParticles());
    }
    const auto expected = contact.evaluate(0, 0);
    const auto observed = contact.observe();
    ASSERT_EQ(observed.size(), 1u);
    EXPECT_EQ(observed[0].contacting_particles, 1u);
    EXPECT_LT((observed[0].reaction + expected.force).norm(), 1e-12);
    EXPECT_GT((expected.force - last_applied).norm(), 1e-6);
    EXPECT_DOUBLE_EQ(observed[0].potential_energy, expected.potential_energy);
    EXPECT_DOUBLE_EQ(observed[0].max_penetration, expected.penetration);
    for (size_t k = 0; k < std::size(vectors); ++k)
        for (size_t i = 0; i < particles.TotalRealParticles(); ++i)
            EXPECT_EQ((particles.getVariableDataByName<Vecd>(vectors[k])[i] - before_vectors[k][i]).squaredNorm(), 0);
    for (size_t k = 0; k < std::size(scalars); ++k)
        for (size_t i = 0; i < particles.TotalRealParticles(); ++i)
            EXPECT_DOUBLE_EQ(particles.getVariableDataByName<Real>(scalars[k])[i], before_scalars[k][i]);
}

TEST(PDContactObservation, FrictionDissipatesSlipAndSamplingIsGalileanInvariant)
{
    SolidSphere solid;
    auto &particles = solid.body.getBaseParticles();
    SimpleDynamics<solid_dynamics::RigidBoxContactforPD> contact(solid.body, 0.3, 0.01);
    const Vecd wall_velocity(0, -0.1, 0), particle_velocity(0.03, -0.02, 0.01), common(1, 2, 3);
    const Vecd center(0, -0.4, 0);
    contact.addBox("Plate", center, Vecd(0.5, 0.15, 0.5), wall_velocity);
    particles.getVariableDataByName<Vecd>("Position")[0] = Vecd(0, -0.21, 0);
    auto *velocity = particles.getVariableDataByName<Vecd>("Velocity");
    velocity[0] = particle_velocity;
    const auto before = contact.evaluate(0, 0);
    EXPECT_LT(before.friction.dot(particle_velocity - wall_velocity), 0);
    EXPECT_DOUBLE_EQ(before.boundary_power, before.force.dot(wall_velocity));
    velocity[0] += common;
    contact.setKinematics(0, center, wall_velocity + common);
    const auto after = contact.evaluate(0, 0);
    EXPECT_LT((before.force - after.force).norm(), 1e-10);
    EXPECT_DOUBLE_EQ(before.potential_energy, after.potential_energy);
    EXPECT_NEAR(after.boundary_power - before.boundary_power, before.force.dot(common), 1e-10);
}

TEST(PDContactObservation, IntegrationVelocityPredictorDoesNotChangeAcceptedStateObservation)
{
    SolidSphere solid;
    auto &p = solid.body.getBaseParticles();
    SimpleDynamics<solid_dynamics::RigidBoxContactforPD> contact(solid.body, 0.3, 0.01);
    contact.addBox("Plate", Vecd(0, -0.4, 0), Vecd(0.5, 0.15, 0.5));
    auto *position = p.getVariableDataByName<Vecd>("Position");
    auto *velocity = p.getVariableDataByName<Vecd>("Velocity");
    StdVec<Vecd> predicted(p.TotalRealParticles(), Vecd(-0.03, 0, 0));
    for (size_t i = 0; i < p.TotalRealParticles(); ++i)
    { position[i] = Vecd(0, 0.5, 0); velocity[i] = Vecd(0.03, 0, 0); }
    position[0] = Vecd(0, -0.21, 0);
    contact.setForceSamplingVelocityData(predicted.data());
    contact.exec();
    EXPECT_GT(p.getVariableDataByName<Vecd>("ContactFrictionForce")[0][0], 0);
    EXPECT_LT(contact.evaluate(0, 0).friction[0], 0);
    EXPECT_GT(contact.observe()[0].reaction[0], 0);
    EXPECT_DOUBLE_EQ(velocity[0][0], 0.03);
}
