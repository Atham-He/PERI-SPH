#include "sphinxsys.h"
#include "pd_checkpoint.h"
#include "pd_surface_quadrature.h"
#include <gtest/gtest.h>

using namespace SPH;
using namespace SPH::peridynamics;

namespace
{
StdVec<SurfaceVector> cloud()
{
    StdVec<SurfaceVector> points;
    for (int i = -5; i <= 5; ++i)
        for (int j = -5; j <= 5; ++j)
            for (int k = -5; k <= 5; ++k)
            {
                SurfaceVector x(0.2 * i, 0.2 * j, 0.2 * k);
                if (x.squaredNorm() < 1) points.push_back(x);
            }
    return points;
}
}

TEST(PDSurfaceContact, PositiveAreaAndRigidSpherePlaneHaveAnalyticIntegral)
{
    const Real pi = std::acos(Real(-1));
    StdVec<Real> error;
    for (size_t rings : {32, 64, 128})
    {
        SphereSurfaceRule rule(SurfaceVector(0, 1, 0), 1, rings, 16);
        Real area = 0, maximum = 0;
        for (Real value : rule.area) { EXPECT_GT(value, 0); area += value; }
        EXPECT_NEAR(area, 4 * pi, 2e-12);
        for (Real indentation : {0.001, 0.003, 0.008, 0.015, 0.03, 0.1})
        {
            Real force = 0, energy = 0;
            for (size_t s = 0; s < rule.area.size(); ++s)
            {
                const auto response = surfacePlaneResponse(rule.position[s] - SurfaceVector(0, indentation, 0),
                    SurfaceVector::Zero(), SurfaceVector(0, 1, 0), 0, SurfaceVector::Zero(), rule.area[s], 1, 0, 0.01);
                force += response.force[1]; energy += response.energy;
            }
            maximum = std::max(maximum, std::abs(force / (pi * indentation * indentation) - 1));
            if (rings == 128)
                EXPECT_NEAR(energy, pi * std::pow(indentation, 3) / 3, 0.01 * pi * std::pow(indentation, 3) / 3);
        }
        error.push_back(maximum);
    }
    EXPECT_LT(error.back(), 0.01);
    EXPECT_LT(error.back(), error.front());
}

TEST(PDSurfaceContact, SignedReferenceTransferReproducesAffineFieldsWithoutClipping)
{
    auto points = cloud();
    SphereSurfaceRule rule(SurfaceVector::Zero(), 1, 8, 8);
    ReferenceSurfaceTransfer transfer(points.data(), points.size(), rule, 0.8);
    Matd F = Matd::Identity(); F(0, 1) = 0.2; F(1, 1) = 0.8;
    const SurfaceVector shift(0.12, -0.23, 0.07);
    auto deformed = points;
    for (auto &x : deformed) x = F * x + shift;
    EXPECT_LT(*std::min_element(transfer.weights.begin(), transfer.weights.end()), 0);
    for (size_t s = 0; s < rule.area.size(); ++s)
        EXPECT_LT((transfer.gather(s, deformed.data()) - F * rule.position[s] - shift).norm(), 1e-12);
    EXPECT_THROW(ReferenceSurfaceTransfer(points.data(), points.size(), rule, 0.01), std::runtime_error);
}

TEST(PDSurfaceContact, FixedGradedRuleResolvesTinyContactCapsWithoutReferenceAdaptation)
{
    const Real pi = std::acos(Real(-1));
    StdVec<Real> force_errors, energy_errors;
    for (size_t rings : {8, 16, 32})
    {
        SphereSurfaceRule rule(SurfaceVector(0, 1, 0), 1, rings, 8, 14);
        ASSERT_EQ(rule.area.size(), 2 * rings * 8 * 15);
        Real area = 0, max_force = 0, max_energy = 0;
        for (Real value : rule.area) { EXPECT_GT(value, 0); area += value; }
        EXPECT_NEAR(area, 4 * pi, 2e-12);
        for (int k = 0; k <= 70; ++k)
        {
            const Real indentation = std::pow(Real(10), -8 + 0.1 * k);
            Real force = 0, energy = 0;
            for (size_t s = 0; s < rule.area.size(); ++s)
            {
                const auto response = surfacePlaneResponse(rule.position[s] - SurfaceVector(0, indentation, 0),
                    SurfaceVector::Zero(), SurfaceVector(0, 1, 0), 0, SurfaceVector::Zero(), rule.area[s], 1, 0, 0.01);
                force += response.force[1]; energy += response.energy;
            }
            max_force = std::max(max_force, std::abs(force / (pi * indentation * indentation) - 1));
            max_energy = std::max(max_energy, std::abs(energy / (pi * std::pow(indentation, 3) / 3) - 1));
        }
        force_errors.push_back(max_force); energy_errors.push_back(max_energy);
    }
    EXPECT_LT(force_errors[1], 0.005); EXPECT_LT(energy_errors[1], 0.005);
    EXPECT_LT(force_errors[2], 0.002); EXPECT_LT(energy_errors[2], 0.002);
    EXPECT_LT(force_errors[2], force_errors[0]);
}

TEST(PDSurfaceContact, TransposeForceTransferConservesForceTorqueAndPotentialGradient)
{
    auto points = cloud();
    SphereSurfaceRule rule(SurfaceVector::Zero(), 1, 8, 8);
    ReferenceSurfaceTransfer transfer(points.data(), points.size(), rule, 0.8);
    StdVec<SurfaceVector> force(points.size(), SurfaceVector::Zero());
    SurfaceVector sample_force = SurfaceVector::Zero(), sample_torque = SurfaceVector::Zero();
    auto potential = [&]()
    {
        Real energy = 0;
        for (size_t s = 0; s < rule.area.size(); ++s)
            energy += surfacePlaneResponse(transfer.gather(s, points.data()), SurfaceVector::Zero(),
                SurfaceVector(0, 1, 0), -0.8, SurfaceVector::Zero(), rule.area[s], 2, 0, 0.01).energy;
        return energy;
    };
    for (size_t s = 0; s < rule.area.size(); ++s)
    {
        const auto x = transfer.gather(s, points.data());
        const auto response = surfacePlaneResponse(x, SurfaceVector::Zero(), SurfaceVector(0, 1, 0), -0.8,
            SurfaceVector::Zero(), rule.area[s], 2, 0, 0.01);
        transfer.scatter(s, response.force, force.data());
        sample_force += response.force; sample_torque += x.cross(response.force);
    }
    SurfaceVector total = SurfaceVector::Zero(), torque = SurfaceVector::Zero();
    size_t active = 0;
    for (size_t i = 0; i < points.size(); ++i)
    {
        total += force[i]; torque += points[i].cross(force[i]);
        if (force[i].norm() > force[active].norm()) active = i;
    }
    EXPECT_LT((total - sample_force).norm(), 1e-12);
    EXPECT_LT((torque - sample_torque).norm(), 1e-12);
    for (int d = 0; d < 3; ++d)
    {
        const Real original = points[active][d], epsilon = 1e-6;
        points[active][d] = original + epsilon; const Real plus = potential();
        points[active][d] = original - epsilon; const Real minus = potential();
        points[active][d] = original;
        EXPECT_NEAR(force[active][d], -(plus - minus) / (2 * epsilon), 1e-9);
    }
}

TEST(PDSurfaceContact, FrictionDissipatesAndRejectsInvalidGeometry)
{
    const auto response = surfacePlaneResponse(SurfaceVector(0, -0.1, 0), SurfaceVector(0.2, 0, 0),
        SurfaceVector(0, 1, 0), 0, SurfaceVector(0.1, 0, 0), 1, 2, 0.3, 0.01);
    EXPECT_LT(response.friction.dot(SurfaceVector(0.1, 0, 0)), 0);
    EXPECT_LE(response.friction.norm(), 0.3 * 0.2);
    EXPECT_DOUBLE_EQ(response.boundary_power, response.force[0] * 0.1);
    EXPECT_THROW(surfacePlaneResponse(SurfaceVector::Zero(), SurfaceVector::Zero(), SurfaceVector(0, 2, 0),
        0, SurfaceVector::Zero(), 1, 2, 0, 0.01), std::invalid_argument);
}

TEST(PDSurfaceContact, CompleteReferenceCacheRestoresBitwiseWithoutRefresh)
{
    auto points = cloud();
    SphereSurfaceRule rule(SurfaceVector::Zero(), 1, 8, 8);
    ReferenceSurfaceTransfer transfer(points.data(), points.size(), rule, 0.8);
    const auto expected = transfer.weights;
    const auto root = std::filesystem::temp_directory_path() / ("surface-cache-unit-" +
        std::to_string(std::chrono::steady_clock::now().time_since_epoch().count()));
    std::filesystem::create_directories(root / "output");
    DiskCheckpoint checkpoint(root, "surface-cache-control", root / "output");
    transfer.registerCheckpoint(checkpoint);
    const auto first = checkpoint.save({0.2, 0, 1, 1, 0});
    std::fill(transfer.weights.begin(), transfer.weights.end(), 0);
    transfer.generation = 99;
    checkpoint.restore(root);
    EXPECT_EQ(transfer.weights, expected); EXPECT_EQ(transfer.generation, 1);
    const auto second = checkpoint.save({0.2, 0, 1, 1, 0});
    EXPECT_EQ(checkpoint_detail::fileHash(first / "state.bin"), checkpoint_detail::fileHash(second / "state.bin"));
    std::filesystem::remove_all(root);
}
