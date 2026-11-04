# S3a: full coupled metric before nonlinear gel integration

The bounded dense research operator retains every off-diagonal entry. It is
not a scalable PD solver, a production gel tangent, or a completed S3 stage.
A positive metric is factored with a long-double normalized Cholesky solve.
Symmetry tolerance is 1e-14 of matrix infinity norm; scaled pivots must exceed
1e-14 of that norm to reject numerical rank deficiency. The operator is never
shifted, clipped or replaced by its diagonal. Linear residual must be <=1e-12.
Native bounds remain128 coordinates/1024 rows; exact audit cases use <=6 DOF.

For each arbitrary dense row C, compliance is the full A^-1 C^T. The projected
dual iteration uses its actual diagonal C A^-1 C^T and applies the identical
row for observation/gather and transpose impulse scatter. Signed interpolation
and two-node flattened3D normals are included. Redundant rows need not have
unique impulses; independent KKT and unique primal comparisons do not assert
a unique dual distribution. Initially penetrating and infeasible geometry
must fail explicitly without a fallback or candidate-state commit.

The quadratic coupled-elastic control uses A=M+dt²K/2, rhs=M*v0-dt*K*q0,
q1=q0+dt*v1, internal force=-K*(q0+q1)/2. K is checked for convexity; null rigid
modes are allowed. Physical kinetic energy and momentum always use actual M.
The discrete ledger is

    Wwall = Delta(0.5*v^T M v + 0.5*q^T K q)
              +0.5*(v1-v0)^T M (v1-v0)+sum(Jn*g0/dt),

up to normal complementarity. These one-step algebra controls intentionally
include impact loss; they do not claim the future trajectory's1% algorithmic-
loss/time accuracy gate or qualify nonlinear/finite-strain material physics.

The independent Fraction reader decides SPD via exact leading principal
determinants and stiffness PSD via all exact principal determinants. It uses
rational Gaussian elimination and active-set enumeration rather than native
factorization or projected sweeps. It reconstructs every coordinate's KKT and
physical transfer rather than relying on a zero resultant. Arithmetic gate is
1e-10 in frozen scales. Negative controls must reject both diagonal operator
substitution and using the effective A quadratic as physical kinetic energy.

Before dispatch freeze the new plan, source commit/ELF/libraries/support/inputs.
Each native5s/256MiB/32MiB; actual F host/output reserves120/150GiB. Preserve
return codes and all raw before parsing; stop subsequent cases on unexpected
execution, resource, integrity or applicable arithmetic/accuracy failure.
S3b nonlinear gel operator/history and S3c-e3D integration/complete new restart
remain required. No original study, old state, production ELF or material fit
is launched or changed by these algebra probes.
