/* Licensed under the Apache License, Version 2.0. */
#ifndef PD_MODEL_ASSEMBLY_H
#define PD_MODEL_ASSEMBLY_H
#include "pd_model_options.h"
#include "pd_body.h"
#include "base_body.hpp"
#include "base_particles.hpp"
#include "base_particle_dynamics.h"
#include "pd_gel_dynamics.h"
#include "nosbpd_dynamics.h"

namespace SPH::peridynamics
{
inline void definePDMaterial(PDBody &body, const PDMaterialOptions &m)
{
    m.validate();
    if (m.model == "elastic")
        body.defineMatterMaterial<HughesWingetSolid>(m.density, m.youngs_modulus, m.poisson_ratio);
    else if (m.model == "j2")
        body.defineMatterMaterial<J2PlasticityforPD>(m.density, m.youngs_modulus, m.poisson_ratio,
            m.yield_stress, m.hardening(), m.kinematic_hardening_modulus);
    else if (m.model == "viscoplastic")
        body.defineMatterMaterial<ViscoPlasticSolidforPD>(m.density, m.youngs_modulus, m.poisson_ratio,
            m.yield_stress, m.plastic_viscosity, m.hardening(), m.kinematic_hardening_modulus);
    else if (m.model == "drucker_prager")
        body.defineMatterMaterial<DruckerPragerPlasticityforPD>(m.density, m.youngs_modulus, m.poisson_ratio,
            m.yield_stress, m.pressure_sensitivity, m.hardening_modulus, 0);
    else
        body.defineMatterMaterial<GelViscoelasticSolid>(m.density, m.youngs_modulus, m.poisson_ratio,
            m.relaxing_shear_ratio, m.relaxation_time);
}

struct PDConstitutiveDynamics
{
    UniquePtr<BaseDynamics<void>> force, commit;
};
/** Select by the actual material type. The optional commit runs once after
 * the final drift, never during observation. No geometry-specific rule. */
inline PDConstitutiveDynamics makePDConstitutiveDynamics(
    BaseInnerRelation &relation, solid_dynamics::NosbPDBondWorkspace *bond_workspace = nullptr,
    std::unique_ptr<NosbPDElasticBatch> elastic_backend = nullptr,
    std::unique_ptr<NosbPDPlasticBatch> plastic_backend = nullptr,
    GelPointBatch *gel_backend = nullptr)
{
    auto &material = relation.getSPHBody().getMatterMaterial();
    PDConstitutiveDynamics dynamics;
    if (dynamic_cast<GelViscoelasticSolid *>(&material))
    {
        using Update = solid_dynamics::GelConstitutiveWithExecution;
        dynamics.force = makeUnique<Update>(relation, Update::Stage::MidpointForce, gel_backend, bond_workspace);
        dynamics.commit = makeUnique<Update>(relation, Update::Stage::AcceptedHistory, gel_backend, bond_workspace);
    }
    else if (dynamic_cast<PlasticSolidforPD *>(&material))
        dynamics.force = makeUnique<solid_dynamics::NosbPDPlasticWithExecution>(
            relation, bond_workspace, std::move(plastic_backend));
    else
    {
        if (!dynamic_cast<HughesWingetSolid *>(&material))
            throw std::invalid_argument("PD constitutive assembly requires a supported PD material");
        dynamics.force = makeUnique<solid_dynamics::NosbPDElasticWithExecution>(
            relation, bond_workspace, std::move(elastic_backend));
    }
    return dynamics;
}

inline UniquePtr<BaseDynamics<void>> makePDFracture(BaseInnerRelation &relation,
    const PDMaterialOptions &material, const PDFractureOptions &failure,
    solid_dynamics::NosbPDBondWorkspace *bond_workspace = nullptr)
{
    failure.validate(material);
    auto &actual = relation.getSPHBody().getMatterMaterial();
    if (!dynamic_cast<HughesWingetSolid *>(&actual))
        throw std::invalid_argument("fracture assembly requires an actual PD material");
    const bool actual_plastic = dynamic_cast<PlasticSolidforPD *>(&actual) != nullptr;
    const bool actual_gel = dynamic_cast<GelViscoelasticSolid *>(&actual) != nullptr;
    if (actual_plastic != material.isPlastic() || actual_gel != material.isViscoelastic())
        throw std::invalid_argument("declared material does not match actual PD material history");
    if (failure.model == "none") return nullptr; // Skip every threshold, including stretch fallbacks.
    if (dynamic_cast<GelViscoelasticSolid *>(&actual))
        throw std::invalid_argument("fracture assembly requires intact-bond compatible material dynamics");
    using namespace solid_dynamics;
    using Law = NosbPDBondFailureLaw;
    if (failure.model == "principal_plastic_strain")
    {
        if (!dynamic_cast<PlasticSolidforPD *>(&actual))
            throw std::invalid_argument("plastic-strain fracture requires actual plastic material history");
        return makeUnique<NosbPDBondFractureWithExecution<BondBreakByPlasticStrain>>(relation, bond_workspace,
            NosbPDBondFailureParameters{Law::PrincipalPlasticStrain, failure.plastic_strain_limit, 0, Dimensions},
            failure.plastic_strain_limit);
    }
    if (failure.model == "principal_stress")
        return makeUnique<NosbPDBondFractureWithExecution<BondBreakByPrinStress>>(relation, bond_workspace,
            NosbPDBondFailureParameters{Law::PrincipalStress, failure.stress_limit, 0, Dimensions}, failure.stress_limit);
    if (failure.model == "stress_shear")
        return makeUnique<NosbPDBondFractureWithExecution<BondBreakBySigma1andSigma3>>(relation, bond_workspace,
            NosbPDBondFailureParameters{Law::StressShear, failure.stress_limit, failure.shear_limit, Dimensions},
            failure.stress_limit, failure.shear_limit);
    return makeUnique<NosbPDBondFractureWithExecution<BondBreakBySigma1andNorm1>>(relation, bond_workspace,
        NosbPDBondFailureParameters{Law::StressPressure, failure.stress_limit, failure.pressure_limit, Dimensions},
        failure.stress_limit, failure.pressure_limit);
}
} // namespace SPH::peridynamics
#endif
