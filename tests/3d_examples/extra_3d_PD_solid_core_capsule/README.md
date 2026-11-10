# NOSB-PD shell / solid-SPH core compression

This example uses the native NOSB-PD volume dynamics for the shell and native
total-Lagrangian solid-SPH dynamics with the DPC material adapter for the core.
There is no WCSPH fluid. Both bodies retain their physical mass and inertia.

`pd_solid_reference_coupling.h` is a reusable 3-D, fixed-reference, bonded-solid
interface. Its force ports handle PD external acceleration and SPH external
force explicitly. Affine MLS reconstructs coincident physical interface points;
one penalty energy supplies opposite forces through transpose scatter. Analytic
frictionless plates act on a separate reconstruction of the physical outer face.
No fluid coupling, prescribed particle velocities, or second repulsion is added.

## Fixed working input

- Core diameter 1 mm; shell thickness 10 um; outer diameter 1.02 mm.
- Core spacing 25 um, 33,552 volume particles.
- Shell angular spacing 25 um, 5,230 directions and 3 Gaussian thickness layers,
  giving 15,690 NOSB-PD volume integration points. Physical shell mass is retained.
- Core DPC input is copied exactly from the beta=62 degree bare-core case.
  Cohesion is paired with beta to retain the earlier model's initial uniaxial
  yield; this is a controlled numerical input, not measured strength.
- EC shell: Hughes-Winget elastic update, E=2.151 GPa, nu=0.3, density=1200 kg/m3.
  This is the working reference, not an additional fitted modulus. Fracture off.
- Bonded interface; absolute contact K=2.8491556e15 Pa/m and
  interface K=1.4245778e15 Pa/m. Plate quadrature spacing 12.5 um.
- Twenty equilibrated stages to 10% outer-diameter shortening. A shared step is
  limited by both bodies, DPC acoustic bound, physical radial spacing and the
  combined interface/contact mass-stiffness bound. Relaxation is numerical.

## Run

Build this directory directly with CMake, or build its target from the 3-D
project. The dated study also provides a build receipt and `build_hybrid.py`
for a verified existing native library. The driver never overwrites a history.

From this example directory, after building:

```sh
/absolute/path/to/pd_solid_capsule --config=input/capsule.cfg \
  --dpc-table="$PWD/input/mcc_beta62.csv" --preflight-only \
  --output=/absolute/path/to/new_preflight_directory
/absolute/path/to/pd_solid_capsule --config=input/capsule.cfg \
  --dpc-table="$PWD/input/mcc_beta62.csv" \
  --output=/absolute/path/to/new_compression_directory
```

The absolute table path is required because output generation uses an isolated
runtime working directory. Source identities produced by CMake use a documented
semicolon-joined hash; the dated build script uses a sorted file/hash manifest.

## Acceptance and present limits

Preflight checks physical mass, PD shape-tensor rank, initial and manufactured
affine deformation, force/torque transfer, and a loaded finite-difference energy
gradient against particle-force work. The full solve requires 20 actual accepted
points, force residual <=0.02, plate imbalance <=0.01, kinetic ratio <=1e-5,
valid DPC domain and deformation gradients, and conservative mapped loads.

The shell is a curved **3-D volume quadrature**, not a validated reduced-order
PD shell theory. Three Gaussian layers do not demonstrate convergence in bending
or horizon. The horizon is 75.375 um, larger than the 10 um thickness; quantitative
shell-mechanics claims need thin-layer/horizon verification. A passed preflight
or completed curve is not an experimental fit, interface identification, full
stored-energy audit, or validation of fracture. These are reported separately.
