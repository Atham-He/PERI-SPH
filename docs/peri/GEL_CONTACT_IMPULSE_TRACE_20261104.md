# Accepted-step contact and impulse observations

The preserved v4 penalty study failed rigid-contact accuracy; its anchored
adaptive-time diagnosis passed the 2% stored-support force criterion at
unchanged physical loading (0.201046% then 0.0517021%). The sparse anchor has
only 41 full frames and cannot resolve all fast contact reactions. Neither
result releases normal penalty, quasi-static response or calibration.

The new v5 driver emits each accepted step's actual midpoint reaction on each
plate separately, each plate's accumulated impulse, accepted body momentum,
momentum-balance residual, kinetic energy, gel free energy and physical viscous
dissipation. Actual sample and accepted time cursors remain explicit. These
are observations only; no force averaging, velocity/mass modification,
artificial damping or physical material change is made.

Wall impulse is the integral of the **last applied** force over its accepted
kick interval. It is not an instantaneous accepted-terminal force. Plate
reactions have the opposite sign to the body's contact force, so the reported
momentum residual is P-P_initial+I_lower+I_upper. A nonzero residual is retained,
not normalized away. Each plate is checked separately; zero net force is not
a test of either individual plate's accuracy.

Compensated sums, their compensation memories, initial momentum, accepted
time and step count are part of the complete checkpoint. The loaded impulse
cursor must equal the driver cursor before any advance. Unit checks include
rejected duplicate/wrong-phase updates, nonzero individual plate impulses with
zero net force, no-advance byte identity and bitwise split continuation.

The separate Python reviewer independently integrates the emitted per-step
reactions using Neumaier summation (the writer uses Kahan summation), checks
reported impulse/momentum arithmetic, and compares stored-frame particle
momentum/energy with the corresponding step aggregates. A negative control
with correct zero net force but missing individual wall impulses is rejected.
The executable describes its required trace format; a frozen plan requiring
that format cannot silently accept a missing trace. These are aggregate
arithmetic and sparse raw anchors, not independent per-pair reconstruction or
independent energy reconstruction from every particle at every step.

A new prospective `--contact-trace-restart` control uses a distinct 0.3um/3us
trajectory, 12 output intervals, 100nm penalty length, unchanged 0.1m/s speed,
the fixed graded full-sphere surface and 4ns candidate cap. It qualifies only
the new recording/history and complete nonzero restart/output rollback. It
does not repeat a prior study or claim physical contact convergence. Existing
resource guards and finite/mass/penetration-spacing checks remain active.
Full normal-constraint, horizon, spatial, larger-deformation and rate-aware
quasi-static qualification still precede 0TO calibration.
