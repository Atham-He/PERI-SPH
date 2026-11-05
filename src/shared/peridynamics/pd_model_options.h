/* Licensed under the Apache License, Version 2.0. */
#ifndef PD_MODEL_OPTIONS_H
#define PD_MODEL_OPTIONS_H
#include "pd_case_config.h"
#include "pd_hardening.h"

namespace SPH::peridynamics
{
/** Geometry-independent material inputs. No contact or failure threshold lives
 * here. "legacy" is accepted only by the explicit input migration function. */
struct PDMaterialOptions
{
    std::string model = "elastic", hardening_model = "none";
    Real density = 1000, youngs_modulus = 1e6, poisson_ratio = 0.3;
    Real yield_stress = 1000, hardening_modulus = 0, kinematic_hardening_modulus = 0;
    Real saturation_increment = 0, saturation_rate = 0, plastic_viscosity = 0;
    Real relaxing_shear_ratio = 1, relaxation_time = 1, pressure_sensitivity = 0;

    bool isPlastic() const { return model == "j2" || model == "viscoplastic" || model == "drucker_prager"; }
    bool isViscoelastic() const { return model == "viscoelastic"; }
    void bind(CaseParameters &p)
    {
        p.bind("material_model", model); p.bind("hardening_model", hardening_model);
        p.bind("density", density); p.bind("youngs_modulus", youngs_modulus); p.bind("poisson_ratio", poisson_ratio);
        p.bind("yield_stress", yield_stress); p.bind("hardening_modulus", hardening_modulus);
        p.bind("kinematic_hardening_modulus", kinematic_hardening_modulus);
        p.bind("saturation_increment", saturation_increment); p.bind("saturation_rate", saturation_rate);
        p.bind("plastic_viscosity", plastic_viscosity); p.bind("relaxing_shear_ratio", relaxing_shear_ratio);
        p.bind("relaxation_time", relaxation_time); p.bind("pressure_sensitivity", pressure_sensitivity);
    }
    PDIsotropicHardening hardening() const
    {
        if (hardening_model == "perfect") return PDIsotropicHardening::Perfect();
        if (hardening_model == "linear") return PDIsotropicHardening::Linear(hardening_modulus);
        if (hardening_model == "saturation")
            return PDIsotropicHardening::Saturation(hardening_modulus, saturation_increment, saturation_rate);
        if (hardening_model == "legacy") return PDIsotropicHardening::Legacy(yield_stress, hardening_modulus);
        throw std::invalid_argument("material requires a supported isotropic hardening law");
    }
    void validateScalars() const
    {
        for (Real x : {density, youngs_modulus, poisson_ratio, yield_stress, hardening_modulus,
                       kinematic_hardening_modulus, saturation_increment, saturation_rate, plastic_viscosity,
                       relaxing_shear_ratio, relaxation_time, pressure_sensitivity})
            if (!std::isfinite(x)) throw std::invalid_argument("PD material requires finite parameters");
        if (density <= 0 || youngs_modulus <= 0 || poisson_ratio <= -1 || poisson_ratio >= 0.5 ||
            yield_stress <= 0 || hardening_modulus < 0 || kinematic_hardening_modulus < 0 ||
            saturation_increment < 0 || saturation_rate < 0 || plastic_viscosity < 0 ||
            relaxing_shear_ratio < 0 || relaxation_time <= 0 || pressure_sensitivity < 0)
            throw std::invalid_argument("invalid PD material parameter range");
    }
    void validate() const
    {
        validateScalars();
        if (model != "elastic" && model != "j2" && model != "viscoplastic" &&
            model != "viscoelastic" && model != "drucker_prager")
            throw std::invalid_argument("unsupported material_model: " + model);
        if (isPlastic())
        {
            (void)hardening();
            if (hardening_model == "perfect" && hardening_modulus != 0)
                throw std::invalid_argument("perfect isotropic hardening requires hardening_modulus=0");
            if (hardening_model != "saturation" && (saturation_increment != 0 || saturation_rate != 0))
                throw std::invalid_argument("saturation parameters require hardening_model=saturation");
        }
        else if (hardening_model != "none" || hardening_modulus != 0 || kinematic_hardening_modulus != 0 ||
                 saturation_increment != 0 || saturation_rate != 0)
            throw std::invalid_argument("elastic/viscoelastic models require hardening_model=none and zero hardening parameters");
        if (model != "viscoplastic" && plastic_viscosity != 0)
            throw std::invalid_argument("plastic_viscosity requires material_model=viscoplastic");
        if (model == "drucker_prager" && (hardening_model != "legacy" || kinematic_hardening_modulus != 0))
            throw std::invalid_argument("current Drucker-Prager supports legacy isotropic hardening and zero kinematic hardening only");
    }
};

/** Resolve old input defaults once, before validation, allocation or output.
 * New model selections use perfect J2 hardening by default, not the legacy
 * specimen defaults. Explicit incompatible parameters are never discarded. */
inline void resolvePDMaterialInputs(PDMaterialOptions &m, const CaseParameters &p, const std::string &legacy_model)
{
    m.validateScalars();
    const bool old = m.model == "legacy";
    if (old) m.model = legacy_model;
    if (!p.wasProvided("hardening_model"))
        m.hardening_model = m.isPlastic() ? (old || m.model == "drucker_prager" ? "legacy" : "perfect") : "none";
    if ((!old || !m.isPlastic() || m.hardening_model == "perfect") && !p.wasProvided("hardening_modulus"))
        m.hardening_modulus = 0;
    if (m.model != "viscoplastic" && !p.wasProvided("plastic_viscosity")) m.plastic_viscosity = 0;
    // Old nonplastic inputs used to carry inactive plastic parameters.
    if (old && !m.isPlastic())
    {
        if (p.wasProvided("hardening_model") && m.hardening_model != "none")
            throw std::invalid_argument("nonplastic material cannot select a hardening law");
        m.hardening_model = "none"; m.hardening_modulus = 0; m.plastic_viscosity = 0;
    }
    m.validate();
}

struct PDFractureOptions
{
    std::string model = "none";
    Real stress_limit = 1e5, plastic_strain_limit = 1, shear_limit = 1e5, pressure_limit = 1e5;
    void bind(CaseParameters &p)
    {
        p.bind("fracture_model", model); p.bind("fracture_stress", stress_limit);
        p.bind("fracture_plastic_strain", plastic_strain_limit); p.bind("fracture_shear", shear_limit);
        p.bind("fracture_pressure", pressure_limit);
    }
    void validate(const PDMaterialOptions &material) const
    {
        material.validate();
        for (Real value : {stress_limit, plastic_strain_limit, shear_limit, pressure_limit})
            if (!std::isfinite(value) || value <= 0) throw std::invalid_argument("positive finite fracture thresholds required");
        if (model != "none" && model != "principal_stress" && model != "principal_plastic_strain" &&
            model != "stress_shear" && model != "stress_pressure")
            throw std::invalid_argument("unsupported fracture_model: " + model);
        if (model != "none" && material.isViscoelastic())
            throw std::invalid_argument("current viscoelastic updater requires intact reference bonds; fracture_model must be none");
        if (model == "principal_plastic_strain" && !material.isPlastic())
            throw std::invalid_argument("principal_plastic_strain fracture requires a plastic material");
    }
};
} // namespace SPH::peridynamics
#endif
