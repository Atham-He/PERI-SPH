# Prospective normal-contact qualification

This study follows the separately reviewed soft-penalty surface control at
968c4eabb. That review passed numerical time/quadrature/particle/restart checks;
it did not qualify rigid plate contact: upper penetration consumed about 93%
of the imposed stroke. Its immutable data and receipt are prerequisites.

## Fixed surface and actual time control

The full reference sphere now supports a fixed geometric partition of the
polar integration variable z. With 14 levels the panels are [0,2^-14],
[2^-14,2^-13], ..., [1/2,1]. Each has its own Gauss rule; the full reference
area and signed MLS transfer are retained. Nothing adapts or refreshes during
deformation or restart. Analytic rigid-sphere/plane force and potential are
checked at 71 indentation/radius ratios from 1e-8 to 0.1. The 16-point rule per
panel must have both errors <0.5%; the 32-point rule must have errors <0.2%.

Temporal refinement uses explicit maximum accepted-step candidates of
4/2/1 ns at the same CFL=0.1. Wave, material relaxation, configured cap,
output-time clipping and mapped contact stability remain active. Initial and
median **actual** dt must both halve. This avoids treating a CFL edit as proof
of temporal refinement when contact stiffness controls the step.

The output records lower and upper maximum accepted penetration separately,
their sum divided by the total prescribed stroke, and all timestep candidates.
The sum is a conservative bound even if the maxima occur at different samples.
Pure accepted observation remains distinct from the actual midpoint load.

## Predeclared penalty ladder and gates

Independent initial states use penalty lengths 100, 25, 6.25, 1.5625 and
0.390625 nm. The contact modulus, material, mass, geometry, physical motion,
MLS support and reference horizon remain fixed. These lengths change numerical
constraint compliance; they are not fitted force multipliers.

The first four points are prospective coarse constraint controls. They are not
claimed as rigid-contact-accurate individually. The final point is eligible
only when its two-wall penetration is <=1% of the total physical stroke and
its force difference from the preceding point is <=2%. This target was recorded
before dispatch; it is additional to, and does not relax, the existing mass,
penetration/spacing, finite-state or execution limits.

Only after the penalty gates pass does the sequence run time refinement,
complete nonzero-history/cache restart and output-tail rollback, 8/16/32-point
panel quadrature, and dx=10/7.5/5 micrometre particle refinement at fixed
40 micrometre nonlocal horizon/support. Any applicable execution, integrity,
resource or accuracy failure stops subsequent cases. Existing state is never
redispatched and failed receipts are never changed to passed.

The separate frozen reviewer reconstructs all stored material/contact fields,
both wall penetrations and force scatter, and compares complete checkpoint
physical bytes. It verifies prerequisite receipt/state/plan SHAs without
re-reviewing or modifying the preceding study.

## Prospective resource budget and remaining release gates

Eight threads, 8192 particles, 98304 accepted steps, 3600 s per case, 8 GiB RSS,
4 GiB output per case, actual host reserve 120 GiB and output reserve 150 GiB
remain the engineering limits of this new gel study on F. The preceding
7208-particle/1889-step surface case took 139.825 s with a 589.476 MiB checkpoint;
the approximately 1.875x larger surface rule and 5000-step fine temporal control
predict roughly 694 s of work, leaving margin within the existing time budget.
No original R8 guard or failed study is modified.

Even successful completion is a small-control result. Nonlocal horizon,
larger deformation, inertia, nominal experimental geometry and conditional
loading-rate sensitivity must still be qualified before 0TO fitting. No solver
speed, mass or material stiffness is changed to conceal a numerical error.
