#include "sphinxsys.h"
#include "pd_particle_compression.h"
#include <gtest/gtest.h>
#include <chrono>
#include <filesystem>
#include <fstream>
#include <limits>

using namespace SPH;
using namespace SPH::peridynamics;

namespace
{
class InputFile
{
    std::filesystem::path path_;
  public:
    explicit InputFile(const std::string &text)
        : path_(std::filesystem::temp_directory_path() /
                ("pd-model-options-" + std::to_string(std::chrono::steady_clock::now().time_since_epoch().count()) + ".cfg"))
    { std::ofstream out(path_); out << text; }
    ~InputFile() { std::error_code error; std::filesystem::remove(path_, error); }
    void parse(CaseParameters &inputs)
    {
        std::string executable = "test", flag = "--case-config", file = path_.string();
        char *argv[] = {executable.data(), flag.data(), file.data(), nullptr};
        int argc = 3; inputs.parse(argc, argv);
        EXPECT_EQ(argc, 1);
    }
};

class ConfiguredBody
{
    static SharedPtr<Shape> shape(bool sphere)
    {
        if (sphere) return makeShared<GeometricShapeBall>(Vecd::Zero(), 0.2, "Sphere");
        return makeShared<GeometricShapeBox>(Transform(), Vecd(0.2,0.2,0.2), "Box");
    }
  public:
    SPHSystem system;
    PDBody body;
    UniquePtr<InnerRelation> reference;
    ConfiguredBody(const PDMaterialOptions &material, bool sphere = false)
        : system(BoundingBoxd(Vecd(-1,-1,-1), Vecd(1,1,1)), 0.1),
          body(system, shape(sphere))
    {
        definePDMaterial(body, material);
        body.generateParticles<BaseParticles, Lattice>();
        reference = makeUnique<InnerRelation>(body);
        system.initializeSystemCellLinkedLists(); system.initializeSystemConfigurations();
        InteractionWithUpdate<solid_dynamics::NosbPDShapeMatrix> shape(*reference); shape.exec();
    }
};

PDMaterialOptions plastic(const std::string &model = "j2")
{
    PDMaterialOptions m; m.model = model; m.hardening_model = "perfect"; return m;
}
} // namespace

TEST(PDModelConfig, StringOptionsAndNumericAliasesRoundTrip)
{
    PDMaterialOptions material;
    CaseParameters input; material.bind(input); input.bindAlias("youngs_shell", "youngs_modulus");
    InputFile original("material_model=j2\nhardening_model=linear\nhardening_modulus=2000\nyoungs_shell=1e6\nyoungs_modulus=1000000\n");
    original.parse(input);
    material.validate();
    EXPECT_TRUE(input.wasProvided("youngs_shell")); EXPECT_TRUE(input.wasProvided("youngs_modulus"));
    PDMaterialOptions restored; CaseParameters second; restored.bind(second); second.bindAlias("youngs_shell", "youngs_modulus");
    InputFile saved(input.serialize()); saved.parse(second); restored.validate();
    EXPECT_EQ(second.serialize(), input.serialize());
}

TEST(PDModelConfig, ConflictingAliasesRejectInEitherOrder)
{
    for (const std::string &text : {"youngs_modulus=1e6\nyoungs_shell=2e6\n", "youngs_shell=2e6\nyoungs_modulus=1e6\n"})
    {
        PDMaterialOptions m; CaseParameters p; m.bind(p); p.bindAlias("youngs_shell", "youngs_modulus");
        InputFile file(text); EXPECT_THROW(file.parse(p), std::invalid_argument);
    }
}

TEST(PDModelConfig, UnknownDuplicateAndUnsupportedSelectionsReject)
{
    for (const std::string &text : {"material_model=unknown\n", "material_model=j2\nmaterial_model=j2\n", "unknown_key=1\n"})
    {
        PDMaterialOptions m; CaseParameters p; m.bind(p); InputFile file(text);
        EXPECT_THROW({file.parse(p); m.validate();}, std::invalid_argument);
    }
}

TEST(PDModelConfig, LegacyGelAndPlasticControlsResolveWithoutChangingTheirLaw)
{
    ParticleCompressionOptions gel; auto gp = gel.bindings(); gel.resolveModels(gp); gel.validate();
    EXPECT_EQ(gel.material.model, "viscoelastic"); EXPECT_EQ(gel.fracture.model, "none");
    ParticleCompressionOptions old; auto p = old.bindings();
    InputFile file("viscoplastic=true\ngel_viscoelastic=false\nhardening_modulus=5000000\nplastic_viscosity=500\n");
    file.parse(p); old.resolveModels(p); old.validate();
    EXPECT_EQ(old.material.model, "viscoplastic"); EXPECT_EQ(old.material.hardening_model, "legacy");
    EXPECT_DOUBLE_EQ(old.material.hardening().value(0.01),
        5e6*0.01 + (1.589*5e5-5e5)*(1-std::exp(-16.93*0.01)));
}

TEST(PDModelConfig, NewPerfectSelectionResetsOnlyUnspecifiedLegacyDefaults)
{
    ParticleCompressionOptions options; auto p = options.bindings();
    InputFile file("material_model=j2\nhardening_model=perfect\nfracture_model=none\n");
    file.parse(p); options.resolveModels(p); options.validate();
    EXPECT_DOUBLE_EQ(options.material.hardening_modulus, 0);
    EXPECT_DOUBLE_EQ(options.material.plastic_viscosity, 0);
    EXPECT_FALSE(options.material.viscoplastic); EXPECT_FALSE(options.material.viscoelastic);
    ParticleCompressionOptions restored; auto q = restored.bindings(); InputFile saved(p.serialize());
    saved.parse(q); restored.resolveModels(q); restored.validate();
    EXPECT_EQ(q.serialize(), p.serialize());
}

TEST(PDModelConfig, ExplicitContradictionsAreNotSilentlyDiscarded)
{
    for (const std::string &text : {
        "material_model=j2\nplastic_viscosity=1\n",
        "material_model=j2\nhardening_model=perfect\nhardening_modulus=10\n",
        "material_model=j2\nviscoplastic=true\n",
        "material_model=viscoelastic\nfracture_model=principal_stress\n",
        "material_model=elastic\nfracture_model=principal_plastic_strain\n"})
    {
        ParticleCompressionOptions options; auto p=options.bindings(); InputFile file(text);
        EXPECT_THROW({file.parse(p); options.resolveModels(p); options.validate();}, std::invalid_argument);
    }
}

TEST(PDModelAssembly, AllMaterialTypesSelectTheirActualMatchingUpdater)
{
    for (const std::string &model : {"elastic", "j2", "viscoplastic", "drucker_prager", "viscoelastic"})
    {
        PDMaterialOptions m; m.model = model;
        if (m.isPlastic()) m.hardening_model = model == "drucker_prager" ? "legacy" : "perfect";
        if (model == "viscoplastic") m.plastic_viscosity = 100;
        ConfiguredBody fixture(m); auto update = makePDConstitutiveDynamics(*fixture.reference);
        EXPECT_NE(update.force, nullptr);
        EXPECT_EQ(update.commit != nullptr, model == "viscoelastic");
        if (m.isPlastic())
            EXPECT_NE(dynamic_cast<InteractionWithUpdate<solid_dynamics::NosbPDSecondStepPlastic> *>(update.force.get()), nullptr);
        else if (m.isViscoelastic())
            EXPECT_NE(dynamic_cast<InteractionDynamics<solid_dynamics::GelConstitutiveUpdate> *>(update.force.get()), nullptr);
        else
            EXPECT_NE(dynamic_cast<InteractionWithUpdate<solid_dynamics::NosbPDSecondStep> *>(update.force.get()), nullptr);
        auto &actual = fixture.body.getMatterMaterial();
        EXPECT_EQ(dynamic_cast<PlasticSolidforPD *>(&actual) != nullptr, m.isPlastic());
        if (model == "j2") EXPECT_EQ(dynamic_cast<ViscoPlasticSolidforPD *>(&actual), nullptr);
        if (model == "viscoplastic") EXPECT_NE(dynamic_cast<ViscoPlasticSolidforPD *>(&actual), nullptr);
        if (model == "drucker_prager") EXPECT_NE(dynamic_cast<DruckerPragerPlasticityforPD *>(&actual), nullptr);
        if (model == "viscoelastic") EXPECT_NE(dynamic_cast<GelViscoelasticSolid *>(&actual), nullptr);
    }
}

TEST(PDModelAssembly, SameConfiguredHardeningWorksForDifferentGeometries)
{
    auto m=plastic(); m.hardening_model="saturation"; m.hardening_modulus=2000;
    m.saturation_increment=700; m.saturation_rate=80;
    for (bool sphere : {false,true})
    {
        ConfiguredBody fixture(m,sphere);
        auto &actual=DynamicCast<J2PlasticityforPD>(&fixture.body,fixture.body.getMatterMaterial());
        EXPECT_NEAR(actual.IsoHardeningFunc(0.01), 20+700*(1-std::exp(-0.8)),1e-10);
        EXPECT_NE(makePDConstitutiveDynamics(*fixture.reference).force,nullptr);
    }
}

TEST(PDModelAssembly, FractureSelectionIsIndependentAndNoneCannotBreakBonds)
{
    auto m=plastic(); ConfiguredBody fixture(m);
    auto &particles=fixture.body.getBaseParticles();
    auto *stress=particles.getVariableDataByName<Matd>("PDStress");
    for (size_t i=0;i<particles.TotalRealParticles();++i) stress[i]=1e6*Matd::Identity();
    PDFractureOptions failure;
    auto off=makePDFracture(*fixture.reference,m,failure);
    EXPECT_EQ(off,nullptr);
    for (const auto &neighbors:fixture.reference->inner_configuration_)
        for(size_t n=0;n<neighbors.current_size_;++n) EXPECT_TRUE(neighbors.bondLive_[n]);
    failure.model="principal_stress";
    auto on=makePDFracture(*fixture.reference,m,failure); ASSERT_NE(on,nullptr); on->exec();
    size_t broken=0;
    for(const auto &neighbors:fixture.reference->inner_configuration_)
        for(size_t n=0;n<neighbors.current_size_;++n) broken+=!neighbors.bondLive_[n];
    EXPECT_GT(broken,0);
}

TEST(PDModelAssembly, PlasticStrainCriterionRequiresActualPlasticHistory)
{
    PDMaterialOptions elastic; ConfiguredBody fixture(elastic);
    PDFractureOptions failure; failure.model="principal_plastic_strain";
    EXPECT_THROW(makePDFracture(*fixture.reference,elastic,failure),std::invalid_argument);
    const auto declared=plastic();
    EXPECT_THROW(makePDFracture(*fixture.reference,declared,failure),std::invalid_argument);
}

TEST(PDModelConfig, ModelAndFractureChoicesArePartOfSerializedIdentity)
{
    auto m=plastic(); PDFractureOptions failure; CaseParameters p; m.bind(p); failure.bind(p);
    const auto intact=p.serialize(); failure.model="principal_stress";
    EXPECT_NE(p.serialize(),intact);
    const auto fractured=p.serialize(); m.hardening_model="linear";m.hardening_modulus=100;
    EXPECT_NE(p.serialize(),fractured);
}

TEST(PDModelConfig, DruckerPragerRejectsUnsupportedHardeningCombinations)
{
    PDMaterialOptions m; m.model="drucker_prager";m.hardening_model="perfect";
    EXPECT_THROW(m.validate(),std::invalid_argument);
    m.hardening_model="legacy";m.kinematic_hardening_modulus=1;
    EXPECT_THROW(m.validate(),std::invalid_argument);
}

TEST(PDModelConfig, InvalidLegacyScalarsCannotBeHiddenByNormalization)
{
    PDMaterialOptions m; m.model="legacy";
    m.hardening_modulus=-1; CaseParameters empty;
    EXPECT_THROW(resolvePDMaterialInputs(m,empty,"elastic"),std::invalid_argument);
    m.hardening_modulus=std::numeric_limits<Real>::quiet_NaN();
    EXPECT_THROW(resolvePDMaterialInputs(m,empty,"elastic"),std::invalid_argument);
}

TEST(PDModelAssembly, MaterialFractureMatrixConstructsOnlySupportedCombinations)
{
    size_t supported=0,rejected=0;
    for(const std::string &name:{"elastic","j2","viscoplastic","drucker_prager","viscoelastic"})
    {
        PDMaterialOptions material; material.model=name;
        if(material.isPlastic()) material.hardening_model=name=="drucker_prager" ? "legacy" : "perfect";
        ConfiguredBody fixture(material);
        for(const std::string &mode:{"none","principal_stress","principal_plastic_strain","stress_shear","stress_pressure"})
        {
            PDFractureOptions failure; failure.model=mode;
            const bool allowed=(name!="viscoelastic" || mode=="none") &&
                (mode!="principal_plastic_strain" || material.isPlastic());
            if(allowed)
            {
                EXPECT_NO_THROW({auto rule=makePDFracture(*fixture.reference,material,failure);
                    EXPECT_EQ(rule!=nullptr,mode!="none");});
                ++supported;
            }
            else
            {
                EXPECT_THROW(makePDFracture(*fixture.reference,material,failure),std::invalid_argument);
                ++rejected;
            }
        }
    }
    EXPECT_EQ(supported,20); EXPECT_EQ(rejected,5);
}
