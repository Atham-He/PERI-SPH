# Original simultaneous normal/two-tangent contact and bounded actual PD

This is an implemented candidate with engineering checks, not a frozen
independent global/friction trajectory or complete restart qualification.
The earlier v15 disk qualification remains the original frozen kernel scope.

Each contact stores one signed nodal partition, a unit normal and two
orthogonal unit tangents. All three gather rows are constructed from that
same partition; their identical transposes scatter physical impulses.
Normal geometry bounds, original physical impulse/I and (Gv-wall)/speed
are evaluated as exact rational quantities. Binary64 products are first
summed as integers at a common binary exponent, then divided exactly.
The squared-norm branch comparison for the disk is exact, including a
stick/slip distinction smaller than the 256-bit working precision.
The exact differences are converted directly to 256-bit binary arithmetic
for the disk value/selected derivative. Original inputs are never replaced
by intermediate binary64 normalized impulses or gathered velocities.

The top-left global block is the supplied actual nonsymmetric dR/dv.
Top-right blocks are the original negative transposes. Normal FB and the
two disk natural-map equations form the lower blocks, including normal
radius coupling and full2x2 tangent derivatives. Conversion enclosures
are propagated through derivative-times-gather assembly. Each original
momentum/contact residual, normal complementarity product, disk excess,
maximum-dissipation work error and linear backward residual is retained.
Precision rejection carries an uncommitted original-input witness.
There is one semismooth Newton/Armijo/JacobiSVD method: no shifts,
symmetrization, row selection, physical mass substitution or fallback.

The old normal-only consistent linear API retains its384-unknown bound.
A separately requested bounded entry permits up to1024 unknowns and
checks the original matrix backward residual. Actual447 and783 matrices
are engineering-tested; a size guard is not a scalability qualification.
The contact API limits primal coordinates to128, contacts to256 and uses
whole3D nodes, so the current actual body has99 coordinates, not a
production8192-node solver. Solver tolerance cannot be loosened above
1e-12 through this new entry.

pd_bounded_gel_friction.h couples the new system to the actual33-node
fixed-reference gel context. Every material trial reads the same immutable
accepted point histories; the full local B/transpose and true material
tangent remain unchanged. Two normal/tangent wall work terms per plate,
physical M kinetic energy, Maxwell loss, normal contact loss, Coulomb loss
and algorithmic loss are accounted separately. Rational physical ledgers
are rounded once using the exact integer-ratio conversion; impulse,
energy and length quantities use their own declared scales. Full interval
geometry and durable output commit are additional caller gates.

BoundedGelFrictionState records point material state, last normal/tangent
impulses, compensated cumulative slip and12 compensated plate/work/loss
ledgers. Registration hooks alone are NOT the complete checkpoint.
The old state remains unchanged after every failed solve or post-solve
gate; a complete uncommitted candidate is explicitly distinguished.
The caller must never commit solely because the nonlinear solver converges.

Engineering checks cover analytic actual-M impacts/sliding/sticking,
zero friction, common body/wall speed, tangent-basis covariance,
nonsymmetric manufactured roots and joint original-matrix differences,
redundant contacts by total scatter, invalid histories/frames/partitions,
precision rejection and explicit larger linear bounds. The actual sphere
checks use every58 full-sphere sample,116 contacts and348 impulse
components (447 unknowns), including two successive coupled material steps,
both nonzero tangent components, separated energy ledgers, frozen dense
whole-interval gap observations, free translation and no mutation on failure.
This is engineering evidence, not a physical convergence claim.

Before any new native study, implement full original-matrix/native history
observation and independent exact-rational/95-digit readers, then freeze
plan/inputs/source/ELFs/compiler dependencies/support/resources/gates.
Readers must reconstruct branch signs from exact original fractions,
stream large records to remain within budget, and independently check all
point histories/tangents, all Newton matrices, work/impulse/energy/KKT and
the chosen derivative, not merely printed norm fields. Zero-mu comparison
against the sealed normal raw requires its SAME114-sample228-contact map
(783 unknowns), not a58-sample substitute. Do not rerun the old normal solver.

Full S3 still requires independent simultaneous global/actual PD friction
trajectories and new complete transactional state/output restart. Archive
all reference/neighbors/shape/HG/surface/cache/material/contact/dual/algorithm/
compensation and committed output cursors, restore to a temporary receiver
with private output staging, validate all semantics, then swap. Explicit
generation selection, preservation of qualification evidence, full-byte
noadvance/split/tailrollback and semantic-corruption no-live-mutation proof
are required. Do not borrow the old archive's fallback/retention behavior.
S4 time/space/horizon/large deformation/inertia/rate and the original
2% force/1% penetration/1% algorithmic-loss targets precede0TO calibration.
