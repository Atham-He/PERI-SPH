# Common PD material and fracture configuration

`pd_model_options.h` owns geometry-independent material and fracture inputs.
`pd_model_assembly.h` builds materials and selects constitutive dynamics from
the actual material type. The particle-compression and liquid-shell drivers
reuse these functions; neither public header reads specimen dimensions, plate
IDs or experimental fitting data.

## Configuration

```ini
material_model=j2
density=1000
youngs_modulus=1000000
poisson_ratio=0.3
yield_stress=1000
hardening_model=perfect
hardening_modulus=0
kinematic_hardening_modulus=0
plastic_viscosity=0
fracture_model=none
```

The same material/failure file can be passed to either driver. Geometry,
contact and loading stay with the calling example. Reusable numerical
controls are in `peri/configs/materials`; their values are not experimental
material fits.

Material choices:

| `material_model` | Constitutive update | Hardening |
|---|---|---|
| `elastic` | Hughes–Winget incremental elasticity | `none` |
| `j2` | rate-independent J2 | `perfect`, `linear`, `saturation`, `legacy` |
| `viscoplastic` | Perzyna-type J2 | same choices as J2 |
| `drucker_prager` | existing associated pressure-sensitive plasticity | currently `legacy` only, zero kinematic hardening |
| `viscoelastic` | existing equilibrium neo-Hookean/Green-Maxwell model | `none` |

Linear/saturation parameters are `hardening_modulus`, `saturation_increment`
and `saturation_rate`. `kinematic_hardening_modulus` is separate, as is
`plastic_viscosity`. `Perfect` here names the isotropic law: strict ideal
plasticity additionally needs zero kinematic hardening and zero plastic
viscosity. Viscoelastic parameters remain `relaxing_shear_ratio` and
`relaxation_time`; DP uses `pressure_sensitivity`.

Fracture is independently selected:

| `fracture_model` | Thresholds | Compatible materials |
|---|---|---|
| `none` | no dynamics object is created or executed | all listed materials |
| `principal_stress` | `fracture_stress` | elastic and plastic materials |
| `principal_plastic_strain` | `fracture_plastic_strain` | plastic materials only |
| `stress_shear` | `fracture_stress`, `fracture_shear` | elastic and plastic materials |
| `stress_pressure` | `fracture_stress`, `fracture_pressure` | elastic and plastic materials |

These are the existing bond-breaking criteria, including their existing
stretch fallback rules. They are not a fracture-energy regularization.
The plastic-strain criterion uses the maximum principal plastic-strain tensor
component, not accumulated equivalent plastic strain or total Hencky strain.
`none` skips the whole criterion, including those fallbacks, rather than
attempting to disable failure by increasing one threshold.

The current gel constitutive updater requires intact reference bonds, so
all its fracture combinations except `none` are rejected. The liquid-shell
driver additionally rejects `viscoelastic`: its FSI integration does not yet
implement that model's staged history commit. That material remains available
through the shared factory and the solid-particle driver. Assembly compatibility
does not establish physical validation of every material/fracture combination.

## Old input compatibility and effective parameters

`CaseParameters` now supports string values and registered numeric aliases.
Both aliases may appear in one file only if their parsed, serialized values
agree. Different values are rejected in either input order. Exact duplicate
keys and unknown keys still reject. Strings and aliases round-trip through
the saved effective configuration.

The liquid-shell aliases include:

| Old name | Canonical name |
|---|---|
| `rho0_shell` | `density` |
| `youngs_shell` | `youngs_modulus` |
| `poisson_shell` | `poisson_ratio` |
| `yield_stress_shell` | `yield_stress` |
| `hardening_shell` | `hardening_modulus` |
| `viscosity_shell` | `plastic_viscosity` |
| `max_tension_stress` | `fracture_stress` |
| `critical_plastic_strain` | `fracture_plastic_strain` |

Without a new `material_model`, old booleans resolve to their previous
material law. The old liquid-shell defaults also retain their previous
fracture selection. New explicit material selections default to no fracture;
an explicit `fracture_model` always controls the criterion. Saved inputs
contain the resolved choices, so subsequent loading does not infer failure
from the material type.

For new explicit J2 selections, omitted hardening uses the perfect law with
zero linear modulus. A non-viscoplastic selection resets only an *unspecified*
legacy viscosity default to zero; explicitly supplied incompatible viscosity
or hardening values reject. If old booleans are supplied alongside the new
model, they must agree with it. Resolved compatibility booleans are retained
in the effective configuration for existing readers.

Normalization validates scalar inputs before removing old inactive defaults.
Malformed numeric data cannot be hidden by changing the material selection.
Legacy hardening still means the versioned original saturation constants;
its independent saturation fields are inactive, not fitted replacements.

All model, hardening, viscosity and fracture fields participate in the
serialized checkpoint identity. The changed example identities are versioned;
old executable/state qualifications are not inherited by the new driver.
The existing checkpoint storage algorithm is unchanged.

The gel description now identifies its already-existing
`solid-particle-spatial-v1` resource profile (8192 particles, 98304 steps,
eight threads). The old preflight still expected 1000 particles and initially
failed on this pre-existing mismatch. It now checks this named profile with
exact limits while retaining the unlabeled legacy 1000-particle check.
Unknown profiles, unlabeled enlargement and changed limits reject. Solver
resource limits, R8 guards and physical acceptance thresholds are unchanged.

## Verification

The dedicated C++ target covers parser/alias round-trips and conflicts,
legacy migration, incompatible explicit options, actual material/updater
selection, use on box and sphere geometries, disabled/enabled failure, and
the 25 material/fracture construction combinations (20 allowed, 5 rejected).
These are constructor/update-path and local criterion tests, not 20 separate
fracture benchmark qualifications.

`peri_model_config_check.py` checks both executable interfaces, every retained
input config, reusable new material files, early refusal without output, and
one new tiny dry integration with a handful of PD steps. The latter verifies
that the serialized perfect-plastic selection runs with zero damage and that
its effective parameters can be read back. No SPH evolution, experimental fit,
old S4 trajectory, or independent high-precision review is dispatched.

```sh
cmake -S peri -B build-model-config -G Ninja -DCMAKE_BUILD_TYPE=Release
cmake --build build-model-config --target extra_3d_PD_model_config_tests \
  extra_3d_PD_material_tests extra_3d_PD_framework_tests \
  extra_3d_PD_capsule_compression extra_3d_PD_gel_compression -j 4
ctest --test-dir build-model-config \
  -R 'extra_3d_PD_(model_config|material|framework)_tests|peri_model_config_cli|peri_gel_model_description' \
  --output-on-failure
```

Recorded result (2026-11-05): seven selected CTest groups passed, including
156 retained framework tests, 13 material tests, 14 configuration/assembly
tests, 41 CLI/integration controls and four profile unit tests. The tiny dry
integration accepted one PD step with zero damage and no numerical failure;
its saved parameters reloaded to the same selections. All five retained
application targets compiled. The initial stale-preflight failure is retained
with the later passing records. This change does not modify the other
session's S4 controller, high-precision reviewer, frozen studies or resource
limits, and does not qualify new physical checkpoint continuations.
