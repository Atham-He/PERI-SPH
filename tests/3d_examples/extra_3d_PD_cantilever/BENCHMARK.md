# Cantilever verification

The original `extra_3d_PD_cantilever` executable remains the large-deflection
demonstration. `pd_cantilever_benchmark` is a separate, intact, small-strain
verification using the current Hughes–Winget NOSB-PD implementation. It does
not build, run, or modify the closed-fluid-cavity prototype.

## Build and run

From the repository root (C++17, CMake, Eigen, Simbody, Boost, spdlog and TBB):

```sh
cmake -S tests/3d_examples/extra_3d_PD_cantilever \
  -B build/cantilever-verification -G Ninja -DCMAKE_BUILD_TYPE=Release
cmake --build build/cantilever-verification --target pd_cantilever_benchmark -j 4
ctest --test-dir build/cantilever-verification --output-on-failure
python3 tools/peri_cantilever_study.py \
  --executable build/cantilever-verification/bin/pd_cantilever_benchmark \
  --output build/peri_validation/cantilever-study --suite full
python3 tools/analyze_cantilever.py build/peri_validation/cantilever-study
```

The runner only requires Python's standard library. Analysis requires NumPy,
SciPy and Matplotlib; plots use Arial when installed. `--suite pilot` measures
cost over a quarter period, `comparison` selects the three operators, and
`--only name ...` selects named cases from the full matrix. Failed cases are
retained, reported and excluded from accuracy metrics. Reuse is allowed only
for identical parameters and executable hash. To rerun failures or change the
binary, use a fresh output directory. CTest runs a bounded smoke test, not the
full research matrix. Long studies must be launched explicitly.

On a sufficiently large machine, `--threads 8 --jobs 2` runs two independent
cases with eight TBB threads each. The coordinator alone writes the aggregate
manifest. The effective thread count is recorded per case and cannot silently
change for a reused result. The `performance` suite runs a short n=6 thread
comparison. The initial study's 1/4/8-thread normalized displacement histories
differed by less than 4e-13 in relative L2; this is not a bitwise reproducibility
claim. Interrupted attempts are retained separately and never used as results.

## Geometry, material, and boundary

- Baseline square section H=W=1 m, free length L=12H; E=10 MPa,
  rho=1000 kg/m³, nu=0.3. These are verification parameters, not calibrated
  capsule or ice properties.
- Integer-indexed cell-centre lattice, dx=H/n, with n=3,4,6,8. Cells span
  x in (0,L), y,z in (-H/2,H/2). The support is at least H and one horizon
  thick. The reference configuration is constructed once and never sorted
  or rebuilt.
- All support particles x<0 are fixed. Their acceleration is zeroed before
  the second half drift and work accumulation. Their unconstrained force
  is recorded as the reaction. Free-body load is not applied to support.
- Default horizon delta=3.015 dx. A separate fixed-delta series uses
  delta/H=0.75375. The production shape tensor and boundary truncation are
  retained; there is no benchmark-specific surface correction or effective
  length fitting.
- No fracture. The production pair numerical damper is enabled by default,
  matching the existing large-deflection demonstration. Use the runner option
  `--no-numerical-damping` for an explicitly labelled undamped diagnostic.
  It is not safe to assume the undamped operator is stable: the initial study
  found strong high-frequency growth at n=4 and n=6.

## Analytical references and loading

Let A=H², I=H⁴/12, beta=1.875104068711961 and
sigma=(cosh beta + cos beta)/(sinh beta + sin beta). The first bending mode is

```
phi(x) = [cosh(beta x/L)-cos(beta x/L)
          -sigma (sinh(beta x/L)-sin(beta x/L))] / phi_raw(L)
omega_EB = beta² sqrt(E I/(rho A)) / L²
T_EB = 2 pi / omega_EB
```

The free-vibration run starts with zero displacement and prescribed velocity
`v_y = a omega_EB phi(x)`, `v_x = -y a omega_EB dphi/dx`, `v_z=0`, where
a/H=0.001. The axial velocity initializes beam-section rotation. The reference
is `u_y = a phi(x) sin(omega_EB t)`. The mass-weighted transverse modal
projection is fit to one exponentially modulated sinusoid plus a constant;
frequency, decay, amplitude and fit residual are reported. The derived
undamped frequency sqrt(omega_d²+decay_rate²) assumes a single viscously damped
mode and is labelled separately; it is not a new exact analytical solution.
The EB frequency ignores 3-D shear, rotary inertia, local clamping details and
nonlocality. Therefore it is a slender-beam reference, not an exact solution
of the discrete or 3-D PD problem. No fitted frequency is substituted into the
reference history error. A half-amplitude run tests linear-response sensitivity.
Fits are rejected when the optimizer fails, reaches a search boundary, or the
relative single-mode residual exceeds 5%. A bounded multimode response may still
be physically meaningful, but is not reported as a resolved fundamental mode.

The static run applies uniform body acceleration g=q/(rho A), where
`q=8 E I a/L⁴`. The EB and shear-corrected static deflections are

```
w_EB(x) = q x² (6L² - 4Lx + x²)/(24 E I)
w_shear(x) = q (Lx-x²/2)/(kappa G A), kappa=5/6, G=E/[2(1+nu)]
w_EB(L) = a
```

Static equilibrium is approached using split exact velocity drag
`dv/dt=-2 omega_EB v` (half-step before/after the physical update). This
benchmark-only relaxation is zero in all vibration runs. Its kinetic-energy
loss is recorded separately. Before interpreting static accuracy require:
last-quarter-period tip drift <0.1%, final RMS physical acceleration/g <1%,
and RMS velocity/(a omega_EB) <0.1%. A failed equilibrium criterion calls for
longer relaxation or investigation, not an accuracy claim.
In addition, check the clamp reaction against -qL with a 1% relative target.
A motionless discrete solution can still fail this global balance when its
internal stabilization forces do not cancel. The analyzer reports this
separately as `stationary_but_reaction_imbalanced`; a small displacement error
does not override that failure.

The reported tip is the mean of the final particle cross section
`x=L-dx/2`. Every analytical value is evaluated at that same position. Profiles
are section means and standard deviations, not a noisy single surface particle.

The analyzer explicitly reports provisional engineering targets: 5% finest-grid
frequency/static-deflection error, 0.1% time-history change, 1% finest-pair
frequency change for grid independence, and 1e-6 relative similarity mismatch.
These are declared assessment targets, not literature-derived universal limits.

Reference derivations:
- TU Graz, [Cantilever beam](https://lampz.tugraz.at/~hadley/memm/mechanics/cbeam.php),
  eigenvalue equation, mode shape and frequency.
- Clemson University, [Strength of Materials and Design](https://www.clemson.edu/cecas/departments/me/academics/graduate/qualifying_exams/media/solid_mechanics/solid-3.pdf),
  uniform-load cantilever deflection table.
- The shear correction above follows direct integration of
  `w'_shear=V/(kappa G A)`, `V=q(L-x)`, with w_shear(0)=0.

## Original / new stabilization

The benchmark uses the same current constitutive update, time integration,
particle set, clamp, material and loading for all three selections:

| Value | Operator |
|---|---|
| 0 | No hourglass stabilization |
| 1 | Original one-endpoint F_i and 12 K/(pi delta³) transcription |
| 2 | Production endpoint-averaged projected pair operator and 12 K/(pi delta⁵) |

The original is transcribed from Icebreaking_PD-SPH
`SPHINXsys/src/shared/particle_dynamics/solid_dynamics/elastic_dynamics.cpp`,
LittleWoodHourGlassControl, lines 488–537. This is a controlled operator
comparison, not a rerun of the entire historical solver or ice-jet paper.
Both use coefficient 0.5 without geometry-dependent calibration. In 3-D, the
original delta^-3 formula does not have the same dimensional scaling as the new
delta^-5 force-density form. Equal numerical coefficients are deliberate for
the implementation comparison, not a claim of equal physical penalty strength.
Operator 1 exists only in this benchmark translation unit.

Net hourglass force divided by sum of force magnitudes, and total torque
divided by sum of |x||force|, diagnose pair cancellation. The entire body,
including support particles before constraint removal, enters these sums.
They are not free-beam momentum conservation claims: the clamp supplies reaction.

## Scale and convergence interpretation

Similarity uses s=0.001,1,10, with all lengths including dx and delta scaled by
s, t by s, displacement by s, initial velocity unchanged, body acceleration
by 1/s and static q by s. E, rho, nu, a/H, CFL and stabilization coefficient
remain fixed. Compare u/H at matching t/T, frequency ratios, particle counts
and bond counts. No per-scale retuning is permitted.

- Time refinement: CFL=0.2,0.1,0.05 at n=4. Compare complete modal histories
  and frequencies. Output intervals clip time steps consistently.
- Coupled spacing/horizon refinement: n=3,4,6,8 at delta/dx=3.015. This changes
  both quadrature and nonlocal length; do not call it fixed-horizon convergence.
  Both original and new operators use this same four-grid sequence, so a
  coarse-grid advantage is not extrapolated to finer grids.
- Fixed-horizon quadrature refinement: n=3,4,6 at delta/H=0.75375. This does
  not establish the local continuum limit.
- Compare successive changes and analytical discrepancies separately.
  Do not infer an order or extrapolated limit from non-monotone differences.

`history.csv`, `profiles.csv`, effective parameters, completion or failure
markers, command logs, code snapshots, hashes and source revision are retained.
`summary.json` and `metrics.csv` contain computed results, including failures.
`elastic_small_strain` is a diagnostic approximation that excludes stabilization
energy and support strain energy; it is not evidence of full energy closure.
The analysis additionally flags free vibration with max kinetic energy over
twice its initial kinetic energy, even if it remains within the small-strain
envelope. This is a boundedness diagnostic, not an energy-conservation proof.
Runs fail on dead/damaged particles, numerical flags, non-finite values, clamp
motion or |F-I|>0.05. Completion alone does not establish beam accuracy.

The full matrix contains 23 cases, including static original/new comparisons
at n=4 and n=8. The original n=8 static run initially lasts 2 T_EB; all other
static runs last 4 T_EB. Every static case must independently pass the same
equilibrium criteria, regardless of its relaxation duration.
