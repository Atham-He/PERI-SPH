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
interface position functionals while preserving polynomial moments.
Particle velocities themselves are unchanged, so physical kinetic energy and
linear/angular momentum are preserved. Trace velocities are derived observations,
not independent massive degrees of freedom. The physical route therefore does
not impose their independent equality as an extra reconstruction constraint:
near-dependent position and velocity fields can otherwise produce extremely
large weights without adding an energetic invariant. The resulting trace-velocity
jump at birth is recorded, never removed by changing particle velocities. Full
force-work and raw energy checks remain required.

The transferred operator is frozen between topology events. Full capsule
checkpoints retain its exact gradient rows and spectral bound, including the
intact initial operator, so algebraic reassociation cannot alter restored forces. These operator targets
are history and are included in checkpoints; they are not discarded as caches.
An unresolved constraint is an error, not permission to erase material.

Boundary quadrature points and explicit crack-side nodes have different ownership.
An inner/outer boundary point is attached to its immutable reference coordinate;
the initially nearest PD node is only a search anchor. If a new crack separates
that node from the point, its owner is resolved to the nearest geometrically
visible node within the original support. Live connectivity and visibility still
filter every interpolation entry. Explicitly sided crack traces retain fixed
material-node ownership and reject a crossing; their side is never reassigned.
Point coordinates, quadrature areas, penalty stiffness and physical states are
unchanged. The existing constrained transfer preserves the boundary value and
the force scatter remains the exact transpose of the rebuilt map.

The same slab verification deliberately puts reference boundary points opposite
their initial search anchors, then cuts the slab. It checks zero crossing entries,
independent fragment motion, physical/constraint energy transfer, force work and
deterministic reconstruction. This is part of the existing specimen's acceptance,
not a new parameter or loading matrix.

The signed physical insertion-energy residual is measured against Gc times new
physical area. Numerical least-squares energy released by removing support is
reported separately as stabilizer transfer loss. It is never counted as Gc.
There is no additive energy offset, velocity kick, force scaling or stress reset.

The coupled convex dual uses a scalar Gershgorin majorizer for each facet's
three-vector block. For its symmetric Hessian H, diag(L_k I)-H is positive
semidefinite by diagonal dominance. This is FISTA under a fixed linear change
of variables, with the same exact half-ball projection and physical KKT/force
checks. It prevents the large compliance of an old open facet from setting the
iteration step for a newly inserted facet. No physical time step, mass or law
parameter is changed. The zero-opening subgradient includes the strength
boundary; a feasible boundary traction alone is not damage. When an emerging
opening requires tighter force accuracy, the same inner problem is refined
before any history is committed.

Bulk forces use explicit half kicks. The surface drift solves coupled cohesive
constraints with physical inverse masses and moving-normal geometry. The budget
reports cohesive force work, irreversible Gc work and time-integration loss
separately. Contact retains the original absolute penalty stiffness and contributes
to the physical-mass time-step bound. This is not implicit integration of the bulk.

## Contact-only fracture option

Use `--fracture-contact-only` on the capsule, slab verification runner and
capsule owner. The same tensile strength test and candidate geometry sever the
reciprocal bonds. The surviving NOSB gradient, stabilizer and one-sided traces
rebuild with the same topology update. The physical-mass integrator advances
explicitly and uses the existing frictionless compression-only crack contact.
Fracture parameters stay in the input. Gc scales the numerical topology energy
check in this option. No cohesive traction, softening or prescribed Gc work
enters the mechanical history.

A cohesive checkpoint can initialize this option only while the shell is
intact. Any cut or active cohesive history is rejected by the native reader.
Checkpoints written in the contact-only mode can resume their own history.
`MaximumNormalOpening_m` measures the reconstructed gap between two faces.

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
   dynamic force work and the raw/accounted time-integration budgets. A subcritical
   ligament may remain after that free tensile history. Save the complete dynamic
   state, then finish a prescribed rigid-fragment opening of 2 delta_f along the
   pending faces' normals. History advances only from actual separation. This
   second phase checks spatial virtual work and complete Gc expenditure; its load
   parameter is displacement, not dynamic time. Keep its energy ledger separate.
   Finally check contact across neighbouring failed patches and state restoration.
   Neither the material inputs nor the original dynamic loading are changed.
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

## Reading the capsule trajectory

`capsule.pvd` contains shell, core, platens and the two actual faces of every
activated crack facet. Shell `BrokenBondFraction` counts removed material
interactions. On crack cells, `FractureWorkFraction` is the area-weighted spent
Gc fraction, `FullyFailedAreaFraction` is the fraction with zero remaining
cohesive traction, and `MaximumNormalOpening_m` is the current positive gap at
that facet's quadrature. A completely decohered face can close under contact;
a colour alone is not an opening measurement. Geometry is not magnified. Same-time insertion cascades publish one final
trajectory frame; every individual strength/energy event remains in the CSV.

The formal runner `tools/verification/run_physical_fracture.py` executes one
physical dynamic history, saves its complete state and resumes it in a fresh
process for the spatial opening-work/contact checks. `--completed-state` can
repeat the remaining checks with the same pinned physics. Coupled neighbour
arrays are sorted together by material-point ID before layered stencils are
built, preserving each bond/kernel value and making cross-process identity
independent of parallel cell insertion order.

## Reproducibility

Build `pd_physical_crack_benchmark` and `pd_snsni_capsule`. Generate the capsule
bank with `tools/verification/generate_capsule_crack_bank.py` using the frozen
inner radius, thickness, PD spacing and existing horizon. Run the formal benchmark
before `pd_snsni_capsule --preflight-only`. The production owner is
`tools/peri_run_physical_capsule.py`; it requires PASSED receipts and matching
source, input and executable hashes. It performs no parameter modification.

`--resume-checkpoint` together with `--resume-manifest` continues a saved capsule
history into a fresh result root after an implementation repair. The owner verifies
the original schema/checkpoint identity and exact config, DPC table, geometry cache
and crack-bank hashes against the new inputs. The original history fingerprint
remains attached to unchanged checkpoint bytes; the current execution has a
separate source/binary fingerprint. Native checksum and restored-force validation
must pass before the next topology event. Prior CSV and trajectory segments remain
part of the same history. A repaired executable still needs fresh slab and coupled
preflight acceptance; there is no automatic retry or initial-state restart.

The constrained transfer is an adaptation of conservation-aware state transfer,
not a literal implementation of Powell-Sabin finite-element remeshing. A useful
reference is Chen, Li and de Borst, Energy conservation during remeshing in dynamic
fracture, IJNME 120 (2019), 433-446, doi:10.1002/nme.6142. The inner projected
surface solver follows the initially rigid cohesive-contact direction in Doyen,
Ern and Piperno (2013), doi:10.1007/s00466-012-0819-2; moving normals and nonlinear
softening are solved in an outer iteration here.
