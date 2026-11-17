# Physical fracture in the layered NOSB-PD capsule

The current mechanical route is one S-NSNI/RKPM DPC solid core and one layered
NOSB-PD shell. The shell still uses finite-volume Biot elasticity. Physical
particle masses, the material table, tensile strength and fracture energy are
unchanged. This is an implementation change, not a material calibration.

## Topology and surface law

A tensile nominal-traction test activates an initially rigid extrinsic cohesive
surface. Its actual reciprocal material bonds are severed. Constitutive gradients,
least-squares stabilization and one-sided material traces rebuild from live
support. No second point-damage multiplier pays fracture energy again.

For the current equal-mode law, effective separation includes the tangential
jump and only the positive normal jump. Its irreversible history is the largest
effective separation. The failure opening is 2 Gc / ft. Compression is carried
by unilateral frictionless crack-face contact. At zero opening, the cohesive
subgradient is a tensile half-ball: t dot n >= 0 and |t| <= ft. Its proximal
projection is solved inside the coupled dual problem, without postprocessing
forces. The zero-gap active set uses the existing nonlinear convergence tolerance
to avoid branch chatter caused by round-off; there is no new physical gap or
material parameter. Contact includes neighbouring
patches after tangential sliding. Cohesion and contact use the same isoparametric
surface geometry. Forces use the exact transpose of the current reconstruction.

Particle positions, velocities and masses do not change during insertion.
Constrained GMLS transfers the current deformation gradient and the current
interface position/velocity functionals while preserving polynomial moments.
The transferred operator is frozen between topology events. These operator targets
are history and are included in checkpoints; they are not discarded as caches.
An unresolved constraint is an error, not permission to erase material.

The signed physical insertion-energy residual is measured against Gc times new
physical area. Numerical least-squares energy released by removing support is
reported separately as stabilizer transfer loss. It is never counted as Gc.
There is no additive energy offset, velocity kick, force scaling or stress reset.

Bulk forces use explicit half kicks. The surface drift solves coupled cohesive
constraints with physical inverse masses and moving-normal geometry. The budget
reports cohesive force work, irreversible Gc work and time-integration loss
separately. Contact retains the original absolute penalty stiffness and contributes
to the physical-mass time-step bound. This is not implicit integration of the bulk.

## Candidate geometry and limits

Capsule cracks are restricted to an explicitly declared radial surface bank;
their activation times and softening are not prescribed. The bank derives macro
Voronoi width from twice the existing PD horizon so separated regions retain
several material columns. Angular subdivision resolves the curved face geometry
and does not introduce crack planes inside those regions. Gc is charged once on
the exact physical annular section; visibility barriers extending into the empty
core volume are not charged. Material spacing and horizon are unchanged.

The bank is a crack-path discretization assumption. The current layered route
does not support fracture through a thickness column, unrestricted branching,
or a general edge-edge collision guarantee. Equal shear and opening toughness
is an assumption because a separate experimental mixed-mode law is unavailable.

## Acceptance sequence

1. The same 768-point plate checks an actual precrack opening, closure, sliding,
   rotation, stress-free separated fragments, force/energy derivatives, physical
   explicit time integration and checkpoint rejection/roundtrip.
2. Reset that specimen to its intact topology, introduce one precrack and let
   three remaining facets activate under the fixed ft criterion. Integrate the
   actual cohesive force work through complete softening. Check the real gap,
   Gc work, insertion residual, raw/accounted energy budgets, contact across
   neighbouring failed patches and restart of the fractured operator state.
3. Check the capsule driver without advancing physical material history:
   unchanged DPC update, affine shell/core kinematics and the derivative of
   tie/platen constraint energy. Bind each receipt to source and executable hashes.
4. Only after both benchmark and driver gates pass, start one capsule history.
   The first accepted point writes a full checkpoint and the same history resumes
   from it, verifying the force before further loading. A failed gate or native
   failure stops the owner; no parameter search or automatic retry is allowed.

The old diffuse directional-damage history remains a diagnostic result only.
It cannot be resumed into this physical-topology model. The old source is in Git
history; there is only one active source tree and one active capsule entry.

## Reproducibility

Build `pd_physical_crack_benchmark` and `pd_snsni_capsule`. Generate the capsule
bank with `tools/verification/generate_capsule_crack_bank.py` using the frozen
inner radius, thickness, PD spacing and existing horizon. Run the formal benchmark
before `pd_snsni_capsule --preflight-only`. The production owner is
`tools/peri_run_physical_capsule.py`; it requires PASSED receipts and matching
source, input and executable hashes. It performs no parameter modification.

The constrained transfer is an adaptation of conservation-aware state transfer,
not a literal implementation of Powell-Sabin finite-element remeshing. A useful
reference is Chen, Li and de Borst, Energy conservation during remeshing in dynamic
fracture, IJNME 120 (2019), 433-446, doi:10.1002/nme.6142. The inner projected
surface solver follows the initially rigid cohesive-contact direction in Doyen,
Ern and Piperno (2013), doi:10.1007/s00466-012-0819-2; moving normals and nonlinear
softening are solved in an outer iteration here.
