/* Licensed under the Apache License, Version 2.0. */
#ifndef PD_PARTICLE_COMPRESSION_H
#define PD_PARTICLE_COMPRESSION_H

#include "pd_case_config.h"
#include <algorithm>

namespace SPH::peridynamics
{
/** Geometry/material/contact inputs have separate ownership. These defaults are
 * small numerical controls, not inferred experimental specimen properties. */
struct ParticleCompressionOptions
{
    struct Geometry
    {
        Real radius = 60e-6;
        Real spacing = 10e-6;
        Real plate_half_width = 100e-6;
        Real plate_thickness = 20e-6;
        Real initial_gap = 0;
    } geometry;

    struct Material
    {
        Real density = 1000;
        Real youngs_modulus = 1e6;
        Real poisson_ratio = 0.35;
        bool viscoplastic = false;
        bool viscoelastic = true;
        Real relaxing_shear_ratio = 1;
        Real relaxation_time = 1e-6;
        Real yield_stress = 5e5;
        Real hardening_modulus = 5e6;
        Real plastic_viscosity = 500;
    } material;

    struct Contact
    {
        Real friction_mu = 0;
        Real friction_velocity = 0.005;
        bool reference_surface = false;
        Real surface_support = 40e-6;
        Real penalty_length = 40e-6;
        int surface_rings = 128;
        int surface_azimuths = 16;
        int surface_polar_levels = 0;
    } contact;

    struct Loading
    {
        Real speed = 0.1;
        Real stroke = 0.5e-6;
        int output_frames = 10;
    } loading;

    struct Numerics
    {
        Real cfl = 0.1;
        Real horizon_ratio = 3.015;
        Real hourglass_coefficient = 0.5;
        Real maximum_dt = 0; // Zero disables this optional temporal refinement cap.
    } numerics;

    CaseParameters bindings()
    {
        CaseParameters inputs;
        inputs.bind("radius", geometry.radius);
        inputs.bind("spacing", geometry.spacing);
        inputs.bind("plate_half_width", geometry.plate_half_width);
        inputs.bind("plate_thickness", geometry.plate_thickness);
        inputs.bind("initial_gap", geometry.initial_gap);
        inputs.bind("density", material.density);
        inputs.bind("youngs_modulus", material.youngs_modulus);
        inputs.bind("poisson_ratio", material.poisson_ratio);
        inputs.bind("viscoplastic", material.viscoplastic);
        inputs.bind("gel_viscoelastic", material.viscoelastic);
        inputs.bind("relaxing_shear_ratio", material.relaxing_shear_ratio);
        inputs.bind("relaxation_time", material.relaxation_time);
        inputs.bind("yield_stress", material.yield_stress);
        inputs.bind("hardening_modulus", material.hardening_modulus);
        inputs.bind("plastic_viscosity", material.plastic_viscosity);
        inputs.bind("friction_mu", contact.friction_mu);
        inputs.bind("friction_velocity", contact.friction_velocity);
        inputs.bind("reference_surface_contact", contact.reference_surface);
        inputs.bind("surface_support", contact.surface_support);
        inputs.bind("contact_penalty_length", contact.penalty_length);
        inputs.bind("surface_rings", contact.surface_rings);
        inputs.bind("surface_azimuths", contact.surface_azimuths);
        inputs.bind("surface_polar_levels", contact.surface_polar_levels);
        inputs.bind("compression_speed", loading.speed);
        inputs.bind("compression_stroke", loading.stroke);
        inputs.bind("output_frames", loading.output_frames);
        inputs.bind("pd_cfl", numerics.cfl);
        inputs.bind("horizon_ratio", numerics.horizon_ratio);
        inputs.bind("hourglass_coefficient", numerics.hourglass_coefficient);
        inputs.bind("maximum_dt", numerics.maximum_dt);
        return inputs;
    }

    void validate() const
    {
        // Programmatic callers receive the same finite-value checks as the CLI.
        for (Real value : {geometry.radius, geometry.spacing, geometry.plate_half_width,
                           geometry.plate_thickness, geometry.initial_gap, material.density,
                           material.youngs_modulus, material.poisson_ratio, material.yield_stress,
                           material.hardening_modulus, material.plastic_viscosity,
                           material.relaxing_shear_ratio, material.relaxation_time,
                           contact.friction_mu, contact.friction_velocity, contact.surface_support, contact.penalty_length, loading.speed,
                           loading.stroke, numerics.cfl, numerics.horizon_ratio,
                           numerics.hourglass_coefficient, numerics.maximum_dt})
            if (!std::isfinite(value))
                throw std::invalid_argument("particle compression requires finite inputs");

        if (geometry.spacing <= 0 || geometry.radius / geometry.spacing < 3 ||
            geometry.plate_half_width <= geometry.radius || geometry.plate_thickness <= 0 ||
            geometry.initial_gap < 0 || material.density <= 0 || material.youngs_modulus <= 0 ||
            material.poisson_ratio <= -1 || material.poisson_ratio >= 0.5 ||
            material.yield_stress <= 0 || material.hardening_modulus < 0 || material.plastic_viscosity < 0 ||
            (material.viscoplastic && material.viscoelastic) || material.relaxing_shear_ratio < 0 || material.relaxation_time <= 0 ||
            contact.friction_mu < 0 || contact.friction_velocity <= 0 || loading.speed <= 0 ||
            contact.surface_support <= 0 || contact.penalty_length <= 0 ||
            contact.surface_rings < 4 || contact.surface_rings > 256 || contact.surface_azimuths < 8 || contact.surface_azimuths > 256 ||
            contact.surface_polar_levels < 0 || contact.surface_polar_levels > 20 ||
            loading.stroke <= 2 * geometry.initial_gap || loading.stroke >= 2 * geometry.radius ||
            loading.output_frames < 1 || numerics.cfl <= 0 || numerics.cfl > 0.2 ||
            numerics.horizon_ratio < 2 || numerics.hourglass_coefficient < 0 || numerics.maximum_dt < 0)
            throw std::invalid_argument("invalid particle geometry, material, contact or loading inputs");

        if (!std::isfinite(duration()) || duration() <= 0 || !std::isfinite(horizon()) || horizon() <= 0 ||
            !std::isfinite(initialUpperFace()) || !std::isfinite(6 * geometry.radius))
            throw std::invalid_argument("non-finite derived particle controls");
    }

    Real duration() const { return loading.stroke / loading.speed; }
    Real horizon() const { return geometry.spacing * numerics.horizon_ratio; }
    Real centerHeight() const { return geometry.radius + geometry.initial_gap; }
    Real initialUpperFace() const { return 2 * centerHeight(); }
};

/** One source of translating-plate geometry and velocity, including exact stop.
 * Caller supplies the actual integration/observation time; no global time read. */
struct CompressionPlateMotion
{
    const Real initial_face, thickness, stroke, duration;

    CompressionPlateMotion(Real face, Real plate_thickness, Real distance, Real time)
        : initial_face(face), thickness(plate_thickness), stroke(distance), duration(time)
    {
        if (!std::isfinite(face) || !std::isfinite(plate_thickness) || !std::isfinite(distance) ||
            !std::isfinite(time) || plate_thickness <= 0 || distance <= 0 || time <= 0)
            throw std::invalid_argument("invalid compression plate motion");
    }

    Vecd center(Real time) const
    {
        return Vecd(0, initial_face + 0.5 * thickness - displacement(time), 0);
    }
    Vecd velocity(Real time) const
    {
        requireFinite(time);
        return time >= 0 && time < duration ? Vecd(0, -stroke / duration, 0) : Vecd::Zero();
    }
    Real displacement(Real time) const
    {
        requireFinite(time);
        return stroke * std::clamp(time / duration, Real(0), Real(1));
    }

  private:
    static void requireFinite(Real time)
    {
        if (!std::isfinite(time))
            throw std::invalid_argument("non-finite compression plate time");
    }
};
} // namespace SPH::peridynamics
#endif
