# Explicit isotropic hardening for the general PD material library

The J2 and Perzyna material constructors now accept an explicit isotropic
hardening descriptor. This supplies a genuine zero-hardening option without
changing the existing scalar-modulus constructors' legacy law. It is independent
of specimen geometry, contact, fluid coupling and fracture.

```cpp
auto perfect = PDIsotropicHardening::Perfect();
auto linear = PDIsotropicHardening::Linear(2.0e4);
auto saturation = PDIsotropicHardening::Saturation(2.0e4, 700.0, 80.0);

body.defineMatterMaterial<J2PlasticityforPD>(rho, E, nu, initial_yield, perfect);
// The final optional argument is the tensor kinematic-hardening coefficient.
body.defineMatterMaterial<ViscoPlasticSolidforPD>(
    rho, E, nu, initial_yield, eta, saturation, H_kin);
```

The material response uses

```
sigma_y(alpha) = initial_yield + q(alpha)
q(alpha) = H * alpha + Q * (1 - exp(-b * alpha)).
```

`alpha` is accumulated equivalent plastic strain. `H` and `Q` have stress
units; `b` is inverse strain. `Perfect()` sets `H=Q=0`, `Linear(H)` sets
`Q=0`, and `Saturation(H,Q,b)` exposes all three parameters. Strict
rate-independent perfect plasticity additionally requires `H_kin=0` and
`eta=0` (or the rate-independent J2 class). `Perfect()` alone does not
disable the separate kinematic or rate contributions.

The descriptor validates finite, nonnegative parameters and rejects negative
or nonfinite equivalent plastic strain. New saturation configurations use
`expm1` to retain small saturation increments. The original constructors
explicitly select `Legacy(initial_yield,H)`, retaining `Q=1.589*y-y`,
`b=16.93` and the original `1-exp` floating-point evaluation. Their `H=0`
continues to include saturation hardening; it is not silently reinterpreted.
`HardeningParameters()` exposes the selected law and actual parameters.

The nonlinear corrector counts the complete isotropic resistance through
`q(alpha)` exactly once. The existing bracketed solve, kinematic convention
and viscosity term are retained. Drucker–Prager is unchanged by this update.

## Verification

The separate `extra_3d_PD_material_tests` target exercises actual material
objects with registered particle histories. Its reference uses uniaxial axial
compatibility and a long-double scalar solve, independently of the tensor
return residual and `YieldFunc`:

```
sigma = E * (epsilon - alpha)
sigma = y + H*alpha + Q*(1-exp(-b*alpha))
        + eta*(alpha-alpha_old)/dt                 # flowing tensile branch
epsilon_22 = epsilon_33 = -nu*sigma/E - alpha/2.
```

The prescribed lateral increments enforce a uniaxial **stress** path. A
fixed lateral-strain test would be a different stress state and must not be
compared with the same uniaxial analytical plateau.

Coverage includes perfect-plastic loading/unloading/residual strain, linear
post-yield tangent, nonlinear saturation, pure shear yield, hydrostatic
no-flow, parameter zero limits, physical-duration independence of J2,
rotation of plastic/kinematic history, zero-viscosity equivalence, and
viscoplastic loading/holding against backward Euler. A fixed physical-rate
time-refinement test separately compares with the continuous perfect-plastic
Perzyna uniaxial solution. Its strain increments and dt are both refined;
the material API's `G` already contains the time increment.

A separate affine PD patch test calls the real shape reconstruction and
`NosbPDSecondStepPlastic` on three particle spacings. It compares every
particle, including surface particles, with the analytical proportional
deviatoric J2 solution for perfect and linear hardening. This is a patch
test, not a convergence certificate for a free-boundary structural problem.

These are constitutive verification controls. They do not establish general
finite-strain plastic energy consistency, a calibrated polymer law, crack
energy regularization, or full SPH–PD coupling accuracy. No new history field
is introduced. Existing case CLI options keep their old meanings; application
configuration/assembly can explicitly opt into this API in a later change.

Run in a separate build/output directory:

```sh
cmake -S peri -B build-materials -G Ninja -DCMAKE_BUILD_TYPE=Release
cmake --build build-materials --target extra_3d_PD_material_tests extra_3d_PD_framework_tests -j 4
ctest --test-dir build-materials -R 'extra_3d_PD_(material|framework)_tests' --output-on-failure
```

Frozen production executables, prior physical trajectories, independent
mathematical qualifications and the other session's S4 work are not rebuilt
or rerun by this material-library task.

## Recorded engineering result (2026-11-05)

A clean Release build in a separate worktree based on `b81069f86` passed
all 13 new material/PD-patch tests and all 156 existing framework tests.
The capsule, gel, ice-jet, brittle-disc and cantilever application targets
also compiled; their physical trajectories were not dispatched.

The initial parallel application build encountered a GCC 13.3 internal
compiler segmentation fault in the gel translation unit's checkpoint
function. The unchanged original source compiled in a separate control.
The changed gel target then compiled and linked with the same `-O3` settings
using one build job. The initial failure and serial-build logs are retained;
this does not establish the cause of that compiler failure. No compiler
workaround or checkpoint-source modification is included in this change.
