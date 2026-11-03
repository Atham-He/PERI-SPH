# Effective-gel baseline development

This change implements the next code-development stage following
`CAPSULE_CODE_ASSESSMENT_20261103.md`. The order is: develop and verify the code,
qualify the zero-oil baseline, calibrate identifiable baseline material parameters,
then extend/calibrate the oil-bearing particle. Contact is checked separately and
can be refined between material-fit stages while the material is frozen.

## Model boundaries

`extra_3d_PD_gel_compression` is a small solid-sphere PD numerical control. It has
one deformable body and two translating analytic rigid boxes. There is no hollow
cavity, dummy liquid allocation, fracture, oil phase or fluid transport. It is an
effective-gel geometry baseline, not proof of a homogeneous experimental interior.

The default 60 micrometre radius, spacing, loading speed and elastic properties
are numerical-control inputs. They are not the specimen diameter, rate or fitted
properties of `Micro-H_0TO-3SA-H`. Neither a group D50 nor the nanoindentation
5/10/5-second protocol supplies those missing microscale inputs.

The available constitutive paths remain Hughes-Winget elasticity and optional
Perzyna-type J2 viscoplasticity. Neither is presented as a calibrated hydrated-gel
viscoelastic/poroelastic law. Verifying the appropriate finite-strain gel material,
including recovery/relaxation and possible water transport, remains development
work before full-curve material calibration.

The existing `extra_3d_PD_capsule_compression` remains the idealized single-core
shell reference. Its fluid-side observer has not been certified as a common-state
instantaneous-reaction measurement by this change. Frozen R8 results retain their
original executable/source provenance and limitations.

## Code organization

| Component | Responsibility |
|---|---|
| `src/shared/peridynamics/pd_particle_compression.h` | Separate geometry, material, contact, loading and numerical input ownership; strict validation; one explicit-time plate motion definition |
| `src/shared/peridynamics/pd_contact.h/.cpp` | Shared finite-box contact force/potential; pure per-particle evaluation and accepted-state observation; integration application remains a separate operation |
| `tests/3d_examples/extra_3d_PD_gel_compression/gel_compression.cpp` | Solid-sphere setup, PD stepping, checkpoint cursors and diagnostics |
| `tests/3d_examples/extra_3d_PD_gel_compression/execution_limits.h` | Fixed prospective small-control cost/resource limits; host-volume check precedes particle/output creation |
| `tools/peri_calibration.py` | Source-verified development/calibration sequence; separate varying-parameter blocks; no optimizer or simulation dispatch |
| `tools/peri_gel_preflight.py` | CLI/model validation and a guaranteed pre-allocation resource-stop negative control |

The gel control has no global mutable model parameters. Geometry is explicitly a
ball rather than the difference of two balls. Contact/material controls are
serialized independently for later sensitivity studies. There is no fitted force
multiplier or automatic density/position/reference reset.

## Reaction and state semantics

`RigidBoxContactforPD::evaluate(particle, box) const` returns force, friction,
penetration, potential and boundary power without touching particle arrays.
`observe() const` reduces the same law into each box's reaction, contact count and
maximum penetration. `update()` applies the evaluated force and writes the
integration fields. These paths share the law, not mutable observer scratch.

The gel diagnostics distinguish:

- `upper_observed_Fy_N`: contact law evaluated at current accepted positions,
  velocities and analytic plate geometry. This is a contact observation, not a
  declaration of quasistatic equilibrium or a total coupled-fluid reaction.
- `upper_last_applied_Fy_N` and `last_contact_sample_time_s`: the midpoint load
  actually used by the preceding PD integration step.
- `observed_contact_potential_J`: accepted-state penalty potential, separate from
  the accumulated work of the applied integration loads.
- `kick_identity_residual_J`: kinetic increment/work arithmetic. It is not an
  independently established complete constitutive/transport energy balance.

Observation never initializes accelerations, reapplies a force, advances time or
adds a terminal kick. Current contact count and maximum penetration are available
for contact qualification, instead of hiding contact error in material fitting.

The inherited normal penalty uses the material's `ContactStiffness()`, derived
from its reference density/sound speed. Therefore fixed friction/regularization
settings do not make normal contact independent of E/nu. The driver records the
actual penalty modulus and penetration/spacing. Contact penetration and reaction
sensitivity must be rechecked over the material-fit candidate range. No free
normal-force multiplier is introduced to conceal that dependency.

## Checkpoints and execution limits

The new driver registers all PD particle state, fixed reference neighborhoods,
initial masses, accepted time/step count, last-contact sample time, completed
output frame and committed diagnostic prefixes with the existing atomic PD
checkpoint implementation. Executable/config/state-recording identity is part of
the compatibility key. Snapshot tails beyond a restored cursor are retained in an
orphan directory. Existing run directories are rejected for a fresh run.

This is checkpoint implementation, not new driver-level restart qualification.
No-advance, output rollback and split-continuous tests of this driver are still
required before release, and must be repeated when gel material memories or new
coupling/topology caches are added.

The control is bounded to 1000 particles, eight threads, 98304 accepted steps,
8 GiB process RSS and 3600 seconds total recorded wall time. Storage reserves
remain 120 GiB on the actual host volume and 150 GiB in the run filesystem.
`PERI_HOST_STORAGE` identifies the actual host volume; the WSL project default is
`/mnt/e`. An unknown host volume is rejected. Resource stops return 125, wall-cost
stops return 124. Neither is evidence of density instability. A stopped run is not
automatically retried or resumed.

## Checks that do not dispatch a solver simulation

```sh
cmake --build build/peri-latest-20261103 --parallel 8 --target \
  extra_3d_PD_gel_compression extra_3d_PD_capsule_compression extra_3d_PD_framework_tests
ctest --test-dir build/peri-latest-20261103 -R extra_3d_PD_framework_tests --output-on-failure
ctest --test-dir build/peri-latest-20261103 -R peri_analysis_unit_tests --output-on-failure
ctest --test-dir build/peri-latest-20261103 -R peri_gel_preflight --output-on-failure
```

The new CTest registrations only describe/validate inputs and force an absent-host
preflight rejection. No gel/capsule trajectory, numerical convergence or material
fit is performed by these tests. The framework tests additionally check force as
the normal-potential gradient, action/reaction, friction dissipation, common-
velocity invariance and exact non-mutation of applied load/acceleration arrays by
observation.

Prepare the source-verified calibration sequence from the existing read-only
experimental export using `tools/peri_calibration.py --manifest <manifest.json>
--output <new-plan.json>`. It retains unknown geometry/time and leaves all
scientific release gates false. It refuses mixed contact/material varying blocks
and force rescaling. It does not treat preparation as a calibration pass.

## Remaining release work

1. Meet the existing storage reserves, then perform a versioned small-control
   integration and restart regression; do not retry the old stopped 48us state.
2. Qualify time/space/contact sensitivity and recoverability of the material law.
3. Resolve or explicitly propagate individual geometry/loading/contact-zero
   uncertainty before fitting parameters it controls. Missing time prevents a
   unique rate-dependent material fit.
4. Freeze qualified contact, calibrate the zero-oil material, and independently
   validate it. Contact refinement uses the frozen material and repeated numerical
   checks before another material-fit stage.
5. Add justified multi-chamber/homogenized oil geometry and transport/coupling,
   then qualify their mass/impulse/work and complete history before oil-case fits.
