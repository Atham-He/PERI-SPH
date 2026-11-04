# Actual bounded3D PD gel integration under development

This candidate uses33volume nodes of a full3D unit sphere at spacing.5,
99physical coordinates. It is a normalized numerical benchmark, not a
measured capsule diameter, calibrated material or scalable production mesh.
Production gel ELF and old frozen studies are unchanged.

FixedPDReference constructs all intact directed neighbors inside horizon1.25,
with canonical target order and w=(1-|xi|²/delta²)^3. It stores shape inverses,
all gradients g_ij=wVj K_i^-1 xi, center partition gradients, volumes,
neighbor measures/weights/active flags and generation. F_i=B_i q uses the
actual local neighborhood, not a single global affine surrogate. All
persistent caches can register for a NEW complete checkpoint. Stored-cache
validation recomputes expected quantities in read-only scratch and never
repairs/replaces an evolved reference. Disconnected/rank-deficient/duplicate
geometry and inconsistent stored caches reject.

The positive nonaffine energy is
eta*mu/(2delta²)*sum_i Vi*sum_j wVj*|(qj-qi)-Fi*xi|²/(sum_j wVj).
Its fixed H^T H scalar node cache is shared across3physical components.
Dividing by neighbor measure preserves physical dimensions. Positive
residual-square evaluation avoids subtractive near-affine energy loss;
force and Hessian are the exact mathematical derivatives of that quadratic
form. eta=.1 is a declared discretization parameter, not material loss or a
force fit. Its sensitivity belongs to S4. The actual combined equilibrium
linearized99x99operator has only6rigid modes in engineering checks; no
general horizon/mesh stability theorem is inferred from this one spectrum.

BoundedGelStepContext holds immutable accepted material history and physical
masses rho0*Vi. Each nonlinear trial evaluates all33qualified affine-F
Maxwell paths. Material force is -sum_i Vi B_i^T Pbar_i, and the full true
Jacobian is M+dt²sum_i Vi B_i^T T_i B_i plus the consistent nonaffine term.
No diagonal replacement, symmetrization or mass shift is performed.

Normal contact retains all fixed signed MLS sample constraints. The new
PolePanelSphereRule has58positive-area samples including two actual poles,
using fixed composite three-point Lobatto panels in mu with cap width1/64.
It integrates constant and second spherical moments exactly to roundoff;
this is not surface-convergence evidence. ReferenceSurfaceTransfer preserves
all signed weights and stores support/generation, with read-only restored
cache validation. Neither clipping nor nodal-only substitution is used.
Both plates give116constraint rows. Whole-plate ledgers may have more than
256original signed products, so the exact arithmetic primitive now accepts
up to32768products; the nonlinear solver still has at most256rows. New
engineering bounds checks exercise the enlarged arithmetic budget; oldv13
qualification remains its original scope until new independent validation.

advanceBoundedGel returns a new candidate only after original all-row KKT,
coordinate momentum, actual-M global balance, material/nonaffine work and
energy gates. Physical Maxwell loss, nonaffine stored energy, first-order
algorithmic loss and normal gap loss remain separate. Plate impulse is
distributed via the same original signed rows, with exact accumulation.
Accepted Cauchy stress/density/history and applied interval force have
different cursors. Angular drift is retained; the old-position torque
identity is bookkeeping, not angular momentum conservation.

Current engineering regressions include local affine reproduction,
nonaffine energy/force derivative, true linearized rigid modes and rejection
of a collapsed global-average operator, cache/shape/geometry rejection,
sampled global nonlinear finite differences, free translation, and one
actual33node/116row compression step with nonzero plate impulses and energy
ledger checks. They do NOT constitute frozen independent trajectory,
contact/time/space/horizon/large-strain/inertia or full restart qualification.

Next: freeze a complete native observer/runner and independent reader for
all-reference/MLS reconstruction, actual material/global operator samples
and multistep normal sphere trajectories. Then qualify3D tangential
friction and a NEW complete material/reference/neighbor/cache/contact/dual/
algorithm/committed-output cursor checkpoint with noadvance/split/tail
rollback/corruption rejection. Full S3 is still active. S4/0TO release is
false; no zero-friction or one-step result substitutes for the full goal.
