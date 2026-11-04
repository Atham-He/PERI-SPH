# Bounded manufactured many-contact control, before PD integration

The four-node global control has at most four contact rows and enumerates
subsets. This candidate instead keeps all original rows in a single normalized
nonlinear complementarity system. Its measured scale will be stated from
actual controls, never inferred from API limits or called production PD.

Momentum is unchanged: R(v)-C^T lambda=0. Its true nonsymmetric derivative J
is retained. With velocity scale V and impulse scale I=M_total*V, use
a=lambda/I and b=(Cv-bound)/V. Each contact equation is
phi(a,b)=hypot(a,b)-a-b. It is zero exactly for a>=0,b>=0,ab=0.
Away from the origin the derivatives are a/hypot(a,b)-1 and b/hypot(a,b)-1;
at the origin this candidate chooses both sqrt(.5)-1. The normalized Newton
matrix has top blocks J*V/I and -C^T, bottom blocks diag(phi_b)*C and
diag(phi_a). No symmetrization or physical mass change is used.

The complementarity formula is described in the official PETSc TAO manual:
https://petsc.gitlab.io/petsc/release/manual/tao/ . That source describes
additional globalization variants. This candidate deliberately has ONE
predeclared linear method and no gradient/runtime algorithm fallback. No
general nonlinear/PD convergence theorem or empirical production claim is
borrowed from PETSc.

Every Newton linear system uses Eigen JacobiSVD, a minimum-norm rank-revealing
solve with relative singular cutoff1e-13. The ORIGINAL unshifted linear
backward residual must be<=1e-12, including rank-deficient cases; an
inconsistent system returns explicit no-commit failure. This truncation
selects an algebraic solution, not a regularized physical metric. Rank,
residual, complete original matrix, residual vector and direction are saved.
All rows remain in final primal, dual and summed complementarity checks.

Fixed limits:128primal coordinates,256contact rows,48Newton updates,
32Armijo halvings, squared-norm decrease1e-4, root/KKT1e-12. Each native
control has5seconds/256MiB address space/RSS/32MiB output; Fhost/output
reserves120/150GiB. Root failure, domain line-search failure or inconsistent
linear solve returns no accepted state. A caller must keep old history
immutable through all trials and commit only a completely accepted root.
This generic manufactured driver has no gel history or physical commit.

The independent95digit reader checks every recorded original Newton linear
equation, all generalized-Jacobian entries and updates, actual nonlinear
residual, every original contact KKT row, and unique primal/scatter against
prescribed exact roots. The manufactured operator is A*v+beta*v^3-load.
An exact rational positive diagonal-dominance bound on Sym(A), plus beta>=0,
proves strong monotonicity. For two KKT roots, monotonicity gives a strictly
positive pairing for unequal velocities, while complementarity/feasibility
give a nonpositive pairing. Thus v and C^Tlambda are unique, but redundant
individual multipliers need not be. The reader never demands arbitrary
multiplier uniqueness. Consistent rank deficiency and inconsistent RHS have
separate exact rank certificates.

Near a biactive corner the selected generalized derivative is sensitive to
roundoff of a signed gap. The native trace therefore records that actual
rounded argument, full original generalized matrix and residual. Independent
high-precision geometry first checks that argument at the original1e-10
normalized accuracy gate, then checks the selected scalar derivatives and
all matrix entries. The independent linear backward residual uses the
recorded ORIGINAL unmodified matrix and residual and retains1e-12. It does
not reconstruct a different corner derivative from truncated coordinates.
Raw floating values are decoded to exact binary fractions before decimal
arithmetic. Physical roots/Jacobians retain the independent1e-10 gate.

Declared controls include mixed active/open/biactive rows, duplicated rows,
row permutation/positive rescaling, signed axis rotation and changed equation
normalizers. A128DOF no-contact and a3DOF/256row control test the two bounds
separately, NOT their combined maximum or scalable sparse performance.
Explicit Newton/domain/infeasible/invalid controls must reject; wrong full
operator symmetrization, nonzero impulse and generalized transpose matrix
are independent negative controls. Any unexpected applicable gate stops
remaining dispatch and preserves exact raw/returncode/traceback/size/SHA.

This is manufactured algebra only. Actual path-consistent gel/PD reference
operators and neighbor/shape caches, nonaffine/hourglass stability, sphere
geometry, material/algorithmic/normal/friction energy and angular behavior,
trajectory force/penetration accuracy and complete new restart remain
unqualified. S3, S4 and0TO release stay false. Original frozen studies and
production gel ELF must remain unchanged; no old state is redispatched.
