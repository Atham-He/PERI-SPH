# Normal-constraint candidate: bounded mathematical qualification

The v4 penalty study remains failed (2.254% penetration, 106.134% force change).
The v5 observation/restart control remains qualified only for recording and
complete restart. No existing data, failed receipt or production solver changes
are required for this new normal-contact prototype.

## Candidate and scope

The literature supports investigating impulse/velocity complementarity for
stiff contact and explicitly assessing discrete contact energy:
[Anitescu and Potra, 2002, DOI10.1002/nme.512](https://www.mcs.anl.gov/papers/P884.pdf),
[Acary, energy properties of nonsmooth elastodynamics](https://arxiv.org/abs/1410.2499).
The implementation here is our **normal-only, frictionless, first-order,
perfectly inelastic** specialization, not a reproduction or qualification of
either full published frictional scheme. Its contact-step dissipation is
algorithmic impact dissipation, not the gel's calibrated viscous dissipation.
An impact law, friction, deformable-body coupling, time/spatial/horizon effects
and restart of future warm starts all remain to be qualified.

For one normal coordinate per material node, retain reference surface weights
phi (including negative MLS weights), sum(phi)=1. Each row C_s=n_s*phi_s has
n_s=+1 or -1 and demands C_s*y_new>=a_s,new. With y_new=y_old+dt*v_new,
b_s=(a_s,new-C_s*y_old)/dt. Solve the convex mass-metric velocity problem

    min_v 0.5*(v-v_free)^T M (v-v_free), subject to C v >= b.

Its KKT conditions are M(v-v_free)=C^T*lambda, lambda>=0, r=Cv-b>=0,
lambda*r=0. A projected dual coordinate iteration uses H_ss=C_s*M^-1*C_s^T.
Primal, natural, stationarity and energy-complementarity residuals must all
converge; an exhausted iteration budget does not permit a commit. Redundant
rows can have nonunique dual distributions; the primal velocity is unique.
Initial penetration is rejected rather than silently repaired. Positions are
the declared time-integration result, not a particle-position clamp/reset.

The same transpose transfers body impulse. Since sum(phi)=1, total body and
wall impulses balance even for signed interpolation weights. With old gap g,
wall velocity u=(a_new-a_old)/dt and lambda in physical impulse units,

    wall_work = sum(lambda*u)
    D_step = sum(lambda*g/dt) + 0.5*||v-v_free||_M^2
    wall_work - (K_new-K_free) - D_step = 0   (up to KKT tolerance).

Nonnegative D_step for initially feasible geometry is assessed independently;
it is not defined by setting the energy residual to zero. No contact penalty,
artificial viscosity, mass scaling or physical gel parameter is used.

## Independent exact qualification

`peri_normal_projection.py` freezes the header, native probe, inputs, controller,
independent reader, ELF and runtime-library SHAs in a new F-disk directory.
The native executable reads bounded input and constructs no PD/SPH body. It
captures its actual /proc argv, PID/parent/start ticks/boot/executable, allowing
fast probes to retain native OS identity evidence even after they exit.

The independent Python reader enumerates active subsets and solves rational
equality systems by exact Gaussian elimination. It does not call or reimplement
the projected iteration. It checks the unique primal, work, energy, dissipation,
native dual KKT and momentum balance, retaining nonuniqueness of redundant duals.
The predeclared independent relative tolerance is 1e-10 in physical mass/velocity
scales, including an SI-scale check. Tests cover free flight, fixed/moving walls,
signed MLS weights, two opposite walls with individually nonzero impulses,
redundant/coupled rows, infeasible geometry and invalid initial penetration.

Four mass-one perfectly inelastic impact chains use dt=0.36/0.18/0.09/0.045
through T=1.7 from y=1,v=-1 against a fixed wall. Each step uses the previous
actual native output, not an expected result. Final impulse=1 and energy
loss=0.5 are analytic checks. The impulse/dt peak is reported separately:
a rigid impact need not have a finite pointwise-force limit. This does not
change the failed gel force criterion or the eventual 2% quasi-static force
target; it distinguishes the mathematical observables of two different tests.

The bounded prototype permits <=128 coordinates, <=1024 rows, <=256 native
calls, <=5s/native, <=256MiB native RSS and <=32MiB raw outputs, with the actual
F reserve >=150GiB. Infeasible/invalid cases are explicit negative controls,
not solver retries. Other failures stop subsequent probes and preserve raw.
No old R8, gel solver, queue, material fit or reference cache is dispatched.
