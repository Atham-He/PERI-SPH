# S2: reduced Maxwell / normal-constraint / Coulomb coupling

This two-coordinate research candidate is separate from the production finite-
strain gel/PD driver. Parameters are prescribed controls, not experimental fits.
Normal coordinate q has U=K0*q²/2+Kv*x²/2, xdot=qdot-x/tau. For a linear accepted
path q=q0+t*v, h=dt/tau, a=(1-exp(-h))/h:

    x1=exp(-h)*x0+a*dt*v
    Favg=K0*(q0+dt*v/2)+Kv*(a*x0+tau*(1-a)*v).

Analytic integration of Kv*x(t)²/tau gives physical Maxwell dissipation.
Small-h series for the quadratic integrals avoid exponential cancellation.
The actual normal mass remains M; the normal effective metric is
M+dt*(K0*dt/2+Kv*tau*(1-a)). It is never used as physical kinetic mass.

At fixed normal impulse Jn>=0 the tangent solution is
Jt=clip(Mt*(Uwall-vfree),-mu*Jn,+mu*Jn). Friction loss is
-Jt*(vt1-Uwall)>=0; the box variational inequality checks maximum dissipation.
There is no tensile contact, fitted force multiplier, penalty, density clamp,
mass scaling, warm-start fallback or reference reset.

With body force f and feasible old gap g0, the full discrete balance is

    Wwall+f*(q1-q0) = Delta(K_actual_mass+U)+D_Maxwell+D_friction
                       +Jn*g0/dt+0.5*||v1-v0||_M²,

up to normal complementarity. Normal inelastic loss Jn*g0/dt and first-order
algorithmic velocity loss are separate columns. Prescribed standalone material
controls test only constitutive work/storage/dissipation; their actuator has no
inertial model. Coupled controls retain normal and tangent actual-mass inertia.

Native records every accepted state, exact material history, each impulse,
interval reaction, all energy components, momentum/KKT/friction residuals and
slip mode. Checkpoints serialize all 15 immutable inputs, accepted state, step
and output byte cursors, seven compensated ledgers, last impulses/wall cursor
and mode without C++ padding. New journals may roll back only to their explicitly
supplied committed checkpoint cursor. Restore validates immutable inputs.

Before dispatch the complete new plan, source commit, relevant sources/ELF/
linked libraries/assessor/inputs are frozen. Each native has 5 s, 256 MiB address
space/RSS, <=10000 steps, and <=32 MiB outputs. F host and output reserves remain
120/150 GiB. A failed applicable arithmetic, resource, execution or finest
accuracy gate stops later cases and leaves immutable terminal evidence.

Independent assessment uses 65-digit Decimal closed-form active/box solutions,
dimensional exponential antiderivatives (not native series), polynomial-plus-
exponential continuum convolution and exact continuous Coulomb interval impulses.
Velocity is an endpoint quantity; normal reaction is interval impulse/dt compared
to the continuum midpoint. Tangential reaction compares exact continuum interval
impulse/dt, including stick transitions. Zero reference channels remain unresolved
and are never used as a zero-net-force proof. Fine nonzero channel L2<=2%,
algorithmic loss/positive input work<=1%, normalized arithmetic/cone<=1e-10.
Signed work and every physical dissipation remain separate.

Controls include relaxation, ramp/hold/release, normal SLS refinement, friction
with fixed preload, combined refinement, zero Maxwell branch against sealed S1
raw (without rerunning/reviewing S1), zero mu, Galilean common translation,
separation, full-state no-advance/split-continuous/tail-rollback byte equality.
The assessor must reject energy-label substitution that leaves the total loss
unchanged and corruption of a nonzero friction force. Passing does not qualify
3D coupled stiffness, finite-strain gel/PD integration, spatial/horizon/larger-
strain convergence, original dimensions, production history restart or calibration.
