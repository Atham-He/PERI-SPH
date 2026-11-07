# All-solid SPH core-shell compression

## Configured whole-capsule calculation

The fixed-input first case is `configs/EC0SS_1mm_first_case.cfg`. It selects the
rate-independent J2 elastic-perfectly-plastic core and elastic plane-stress shell,
not an elastic-only trial. Paper geometry and assumed material values are marked
separately in the file. Run the complete 0--5% compression case with:

```sh
../all-sph-build/capsule/test_3d_solid_core_shell_compression \
  --config=tests/3d_examples/test_3d_solid_core_shell_compression/configs/EC0SS_1mm_first_case.cfg \
  --output=/absolute/path/whole-capsule-01
```

Each of 20 displacement levels must converge, with three successive checks of
nodal-force residual, plate imbalance and total kinetic energy (including shell
rotation). No success summary is written if a level fails to converge. The
dynamic-relaxation clock is numerical, not experimental time. Both constitutive
models are rate-independent; `equilibria.csv` also lists displacement divided by
the paper speed 1 micrometer/s as an experimental reference coordinate.

For this spherical, isotropic, gravity-free, frictionless problem the static
solution is computed in a centered frame to eliminate a slow rigid translation.
Convert to the fixed-lower-plate frame by `z_fixed = z_centered - displacement/2`.
VTP files use the centered frame. This equivalence is for static equilibria, not
for a physical dynamic test with an accelerating frame.

The three milestones are fixed inputs/constitutive choices, completion of all
equilibrium levels, and interpretation of the whole-case force/deformation data.
Material identification, fracture and full spatial/contact convergence remain
outside this first configured calculation. The earlier accelerated modes below
remain available for development but are not run as additional small prototypes.

This **uncalibrated demonstrator** uses a native SphinXsys volume-SPH core and
native SphinXsys midsurface-SPH shell. No PD, liquid or fracture is involved.
The two bodies retain separate masses and native stress integration. A reference
MLS penalty tie, with affine reproduction and transpose force scatter, approximates
a fully bonded interface. It is case-local and is not presented as an upstream
validated core-shell implementation. Both normal and tangential displacement
are tied; debonding and slip are unsupported.

## Build in a separate directory

From the repository root:

```sh
cmake -S solid_capsule -B ../all-sph-build -G Ninja \
  -DCMAKE_BUILD_TYPE=Release -DCMAKE_CXX_FLAGS_RELEASE=-O1
cmake --build ../all-sph-build --target test_3d_solid_core_shell_compression -j2
ctest --test-dir ../all-sph-build --output-on-failure
```

The ordinary root build also discovers this example when 3D examples are enabled.
The dedicated entry avoids building unrelated examples. It does not alter the
existing `peri` entry or the other chat's data. GCC 13 optimization encountered an
internal compiler error at `-O3` in an existing library source on this workstation;
`-O1` is used for the isolated validation build.
The same compiler also ICEd in the existing CK `derived_solid_state.cpp`; the
dedicated build entry compiles that one source at `-O0` for GCC 13 only.

## First runs

Use a **new or empty** absolute output directory each time; outputs are never erased.

```sh
EXE=../all-sph-build/capsule/test_3d_solid_core_shell_compression
"$EXE" --case=bare --output=/absolute/path/bare
"$EXE" --case=coated --output=/absolute/path/coated
"$EXE" --case=bare --material=j2 --output=/absolute/path/bare-j2
"$EXE" --patch-test
```

For a repeatable five-run check (elastic/J2 bare/coated, plus a half-step coated
comparison), use a new output directory:

```sh
python3 tests/3d_examples/test_3d_solid_core_shell_compression/verify_runs.py \
  --exe ../all-sph-build/capsule/test_3d_solid_core_shell_compression \
  --output /absolute/path/verification
```

The verifier records source/executable hashes and confirms that the J2 histories
actually evolve. Its completion is a runtime check, not quasi-static or experimental
acceptance. The CTest suite covers affine interface force/torque/power, unloaded
equilibrium, loaded integration and invalid input. Test outputs use new temporary
directories so repeated CTest runs do not overwrite earlier evidence.

All quantities are SI. Defaults: diameter 1 mm, shell thickness 10 micrometers,
particle spacing 0.1 mm, core density 1355 kg/m3 (provisional paper reference),
shell density 1000 kg/m3 (demonstration assumption), core E=1 MPa, shell E=10 MPa,
both Poisson ratios 0.3. **Neither modulus is an identified material property.**
The core radius is 0.49 mm in the coated case and 0.5 mm in the bare case; the
shell midsurface radius is 0.495 mm. The samples are separate nominal geometries,
not a claim of identical experimental pre-coating and post-coating cores.

The default 1% compression uses a cosine ramp over 0.3 ms followed by a 0.3 ms
hold. **This is accelerated numerical loading, not the paper's 1 micrometer/s
test.** Both speeds are written in `parameters.txt`. No fitted experimental
prediction or quasi-static qualification is claimed merely because a run completes.

Useful options (all `--key=value`):

- `case=bare|coated`, `material=elastic|j2` (core material only)
- `diameter`, `thickness`, `spacing`
- `core-E`, `shell-E`, `core-nu`, `shell-nu`, `core-density`, `shell-density`
- `yield-stress` (default 20 kPa, demonstration only), `hardening` (default zero)
- `strain`, `load-time`, `hold-time`, `dt-scale`, `cfl`
- `contact-factor`, `tie-factor`, `drag` (artificial mass drag, default zero)
- `max-steps`, `output`

The time step is the minimum of the native core/shell acoustic limits, conservative
interface stiffness limit and plane-contact stiffness limit. It is not a supplied
fixed step. `dt-scale=0.5` halves this limit for a time-step comparison.

## Outputs and interpretation

- `parameters.txt`: all inputs and assumptions, actual masses, particle counts,
  loading speed and stiffness-based time-step bounds.
- `history.csv`: both plate reactions (positive compressive magnitudes), displacement,
  translational kinetic energy, elastic core energy (only meaningful in elastic
  mode), tie/contact energy, tie gap and penalty penetration, deformation determinant,
  interface force/power residual, external work and optional drag loss.
- `vtp/`: native core and shell snapshots for ParaView.
- `summary.json`: completion status and extrema. A nonzero exit means failure;
  no success summary is written after a caught exception.

Shell strain/rotational energy and plastic dissipation are **not yet audited**.
Consequently these columns do not constitute a complete system energy balance.
The interface is compliant, not an exact kinematic constraint. Force/power patch
checks do not replace penalty sensitivity, spatial convergence or experimental
validation. The coarse core particle volume is normalized to the continuum sphere
volume; the resulting boundary quadrature needs refinement studies before Hertz
or material identification can be claimed. Core contact uses a reference offset
from the outermost axial particle center to the specified sphere pole (half a
spacing for the default lattice); shell contact uses a half-thickness offset.
This is a discrete plane-penalty approximation, not an independently resolved
continuum contact patch.

The immediate milestone is a stable native all-SPH run with nonzero plate loading
and demonstrably reciprocal interface forces. Subsequent qualification includes
particle refinement, step/load/contact/tie sensitivity and the selected 0SS/EC0SS
experimental comparison. Original pressure-fluid coupling is not reused.

## Whole-curve calibration and frozen-core comparison

`calibrate_paper_curves.py` performs two main acceptance stages, without launching
the legacy prototype tests. It calibrates native J2 perfect-plastic core E and
yield stress to digitized **experimental markers** for the nominal 1 mm bare
particle in Figure 20(a) of DOI 10.1007/s11440-026-03253-6, then freezes those
parameters for the coated calculation against Figure 20(b). Model-fit lines in
the paper are not calibration targets. `configs/paper_curve_calibration.cfg`
contains the shared settings and fixed shell assumptions.

Calibration uses 0.5-8% compression; 8.5-10% is held out. Each material-search
evaluation solves the same full 0-10% history. The supported input compression
limit is now 15% for these pre-fracture calculations; no rupture prediction is
implied. The driver exploits force scaling with E at fixed yield/E, and verifies
the selected actual E/yield by a final full bare solve before the coated solve.
It writes the frozen core parameters and hash before starting the coated solve.
The material search stops once fit NRMSE is at most 5% and maximum fit error is
at most 10% of the reference peak; holdout data do not select parameters. It does
not claim a global optimum or a completed identifiability profile. A control-only
resume may reuse completed histories after checking unchanged solver, input and
reference-data hashes; the previous launch receipt is preserved.
Acceptance thresholds and failure states are explicit in the driver's result.
Fixed nu, shell properties, contact/tie factors and discretization are not fitted.
Parameters remain effective, discretization-dependent estimates; identifiability,
experimental variability, spatial convergence and full energy accounting are not
inferred from curve agreement.

## Unified exterior contact (2026-11-07)

`plate-contact=outer_surface` uses the same antipodal spherical quadrature at the
nominal OUTER radius, the same reference-area weights, and the same frictionless
plane penalty energy in both cases: U = sum_q k_area A_q p_q^2 / 2. The outer
reference geometry is independent of whether a shell exists. No Hertz force law
or paper-fit force curve is inserted. `legacy` remains available for old runs.

For the bare body, virtual surface coordinates are obtained with a fixed-reference
affine-complete MLS map of core particle positions; the transpose map scatters
the contact load. Virtual points carry no mass or material. For the shell, each
point is the actual outer offset x_mid + t*n_director/2. Its force supplies both
the midsurface force and the director load t*f/2, which is included in the native
shell rotational-acceleration equation by a local dynamics subclass. Constitutive
classes and their stress updates are unchanged. Contact stability bounds include
the mapped core stiffness or shell translational/director inertia, respectively.

The shell and bare contact quadratures use identical angular samples and outer
reference areas. Shell material mass still integrates the midsurface area times
thickness; it is not replaced with contact quadrature mass. The bonded-interface
law remains unchanged. At this nominal 1 mm / dp=62.5 um resolution, both have
806 exterior points and one initial highest point per platen.

The complete main histories record contact force/torque/power transfer residuals
and `contact_surface/surface_<frame>.csv` at the initial state and every accepted
load point. The calibration driver checks affine reproduction, identical initial
exterior geometry and mapping conservation within the two main acceptance stages.
It calibrates the bare curve, freezes the actual verified core parameters, then
runs the coated case. Fit-only adequacy stops further search at NRMSE <=3% and
maximum error/reference peak <=8%; holdout does not choose parameters. Acceptance
still uses the declared 10% fit / 15% holdout / 15% coated thresholds.
Use `configs/unified_surface_calibration.cfg`. The current complete-energy and
spatial-convergence limitations continue to apply.
