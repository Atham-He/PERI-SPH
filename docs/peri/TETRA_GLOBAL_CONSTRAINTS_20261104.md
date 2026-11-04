# Bounded nonsymmetric global contact before PD sphere integration

This is a four-node3D affine reference-gradient control, not a sphere or a
scalable/nonlocal PD implementation. Reference vertices are0,e1,e2,e3; volume
is1/6, gradients(-1,-1,-1),e1,e2,e3. B maps12nodal coordinates to9deformation
entries. Its exact transpose distributes material forces. Each gradient sums
to zero for translation, and B applied to reference vertices is identity.
The control moduli K=5,G=1,Gv=.5,tau=.3 and nodal masses.01 are prescribed
normalized numerical inputs, not a target experiment or calibrated material.

Momentum residual is M(v1-v0)+dt*V0*B^T*Pbar-C^T*lambda. The actual analytic
Jacobian M+dt²*V0*B^T*T_path*B retains every nonsymmetric entry. Material paths
and histories are those independently qualified in the affine-F control.
Each nonlinear trial re-evaluates from immutable accepted history; accepted
history is returned only after KKT and precommit physical ledger guards pass.

The versioned algorithm enumerates at most16active subsets for<=4rows. Dependent
row subsets are identified by rank, while every original row remains in the
final KKT check. Each independent equality root uses cold zero iteration
variables (not a physical reset), analytic true Jacobian and scaled partial-
pivot nonsymmetric Gaussian elimination. No Cholesky, matrix symmetrization,
diagonal substitution, mass shift or alternate fallback is used. Fixed limits:
32Newton steps (<=64allowed by API),24line-search halvings (<=32by API),Armijo
1e-4,normalized root tolerance1e-12,linear residual1e-12 and relative pivot1e-14.
Domain guards evaluate the full cubic determinant path before trial evaluation.
Other material/arithmetic errors propagate; they are not masked as domain errors.
A failed global linear solve/Newton/line search aborts the case. An equality root
with wrong dual sign or inactive feasibility is rejected as that active set;
all sets exhausted yields an explicit no-KKT-root/no-commit status.

Every accepted coordinate's actual-M momentum, both wall impulses, storage,
physical material loss, normal gap loss and first-order algorithmic loss are
recorded. Wall work equals Delta(kinetic+stored)+materialloss+normalgaploss+
algorithmicloss. Work -dt*fint.v equals V0*Pbar:dF via B/B^T. No force/work
correction enforces these equalities. Native precommit arithmetic gates1e-10
and a separate95-digit independent reader check each coordinate and full12x12
Jacobian by numerical differences of independent48point path integrals.

Angular momentum is explicitly diagnostic: Delta L-external impulse torque at
old-position cursor equals dt*internal torque at that cursor. This bookkeeping
identity is verified, but the internal old-cursor torque need not vanish for
this first-order path-average step. The actual drift is retained. Local frame
objectivity and a small old-cursor residual are not angular-conservation proof.
Time refinement/angular qualification,1%algorithmic-loss and2%reaction-curve
goals remain open for trajectories; these are one-step algebra controls.

Freeze plan/source/ELF/libs/all inputs and independent support before native
dispatch. Native5s/256MiB address space/RSS/32MiB output; Fhost/output reserves
120/150GiB. Unexpected applicable execution/integrity/resource/accuracy failure
stops subsequent cases and preserves precise raw/returncode/traceback/SHA.
Existing S3a quadratic raw is only a reference for the new algorithm's limit;
the old solver/state/review is never redispatched. Negative controls must reject
wrong transpose force, symmetrized full Jacobian and using effective J as
physical kinetic mass. A separate single immutable review is mandatory.

Production PD sphere/neighbor-cache integration, scalability, full new history
checkpoint, time/space/horizon/larger-strain/inertia qualification and0TO fitting
remain gated. No production gel ELF or original sealed study is changed.

Stopped v11 is preserved: native one-step-limit negative control correctly
returned2/unchanged history, but the original reader wrongly expected a terminal
residual trace entry (the native guard returns before writing that entry).
A separate95-digit independent Jacobian reconstruction proves the recorded full
Newton update leaves normalized merit1.39365e-5>original1e-12. A new assessor
freezes this stronger mathematical guard witness, keeps all physical/input/ELF/
threshold values unchanged, reuses11existingraw anchors and executes only the
3unattempted invalid controls. Old state/terminal/raw/reader remain stopped,
not passed. Complete termination/trial trace instrumentation is still required
in the later trajectory/restart version, before production qualification.

Post-qualification diagnostic-only source change now records the terminal
Newton-limit residual before returning. Its regression checks that last entry.
Qualified v11/v11b native binaries/raw remain unchanged and are not rerun.
No accepted root equation, Jacobian, history or physical gate was altered.
