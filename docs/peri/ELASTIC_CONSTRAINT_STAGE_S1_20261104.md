# Stage S1: elastic coupling before production PD integration

This is a scalar spring/mass and moving plane control, not a capsule simulation.
Mass M remains physical. A=M+dt²K/2 is only the positive diagonal metric of the
velocity projection. The mass-specialized existing API retains its physical
kinetic ledger; the new generic API explicitly labels its quadratic forms.

From M(v1-v0)=-dt*K*(q0+q1)/2+C'lambda and q1=q0+dt*v1,
multiplication by v1 gives

    Delta(0.5*M*v²+0.5*K*q²)
      = lambda*u_wall - lambda*g0/dt - 0.5*M*(v1-v0)²
        + lambda*(C*v1-b).

The final term is zero at complementarity. Thus the explicit algorithmic loss
is nonnegative from feasible initial geometry. It includes first-order velocity
loss in free flight; it is not gel viscosity or a calibrated restitution law.
No effective-metric quadratic energy is used as physical kinetic energy.

Frozen prospective control: M=.01, K=1, q0=v0=0, D=.2, T=2,
w=D*(10s³-15s⁴+6s⁵), dt=.05/.025/.0125. Exact maintained-contact solution is
q=w, v=w', force=K*w+M*w''. Reaction positivity follows from w>=D/2 for s>=.5
and |w''|<6D/T² (the extrema are at (3±sqrt(3))/6); the earlier half has
nonnegative w and w''. K/2>6M/T² proves admissibility throughout.

Each native trajectory advances only its own state. The assessor uses a separate
exact Fraction scalar active/inactive solution and the continuum derivatives.
Endpoint velocity is compared to w'(t1); interval impulse/dt to force at its
midpoint, explicitly not a fabricated endpoint force. Nonzero reference L2 norms
are used. Every step and physical energy/momentum arithmetic is checked.

Before native dispatch, source commit, source files, ELF, libraries, support,
prospective plan and inputs are frozen by size/SHA. The runtime stops on resource,
execution or arithmetic failure. Review is a separate immutable process/receipt.
Predeclared gates: normalized arithmetic/KKT/nonpenetration 1e-10; finest
force/velocity L2 2%; algorithmic loss/positive wall work 1%. Negative controls
must reject the use of A as kinetic mass and a corrupted nonzero reaction.

Passing this stage permits further reduced-model viscoelastic/friction research.
It does not release production PD, full-history restart, spatial/horizon or
large-strain convergence, original dimensions or 0TO calibration.
