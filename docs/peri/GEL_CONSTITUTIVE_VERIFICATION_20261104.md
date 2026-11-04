# Finite-deformation gel material and verification protocol

This version adds a recoverable effective-gel material and a frozen verification
runner. Parameters are numerical examples, not fitted values. Solvent transport,
permanent plastic flow, damage and oil chambers remain outside this model.

## Constitutive definition

Let `C=F^T F`, `J=det(F)`, `Edev=dev((C-I)/2)`, with symmetric deviatoric internal
strain `Ev`. Reference-volume free energy is

```
Psi(F,Ev) = Ginf/2 * (J^(-2/3)*tr(C)-3)
          + K/4 * (J^2-1-2*ln(J))
          + Gv * ||Edev-Ev||^2.
dEv/dt = (Edev-Ev)/tau.
P = dPsi/dF at fixed Ev; sigma = P*F^T/J.
Ddot = 2*Gv/tau * ||Edev-Ev||^2 >= 0.
```

The equilibrium potential follows the compressible neo-Hookean formulation in
[Sandia's LAME manual, equations 4.8–4.10](https://www.sandia.gov/files/sierra/SM_LAME_5_26/main/models/neo_hookean.html).
The reference-Green-strain Maxwell branch and evolution above define our chosen
effective model. They are not a claim that a particular published multiplicative
viscoelastic model was reproduced. Energy/internal-variable thermodynamic
formulations are discussed, for example, by
[Kunzemann et al.](https://doi.org/10.1177/1045389X241283074); their fitted PDMS
parameters are not transferred to alginate.

The branch is deviatoric in reference Green strain; it is not a strictly
isochoric multiplicative branch at arbitrary finite deformation. Its admissible
use range still requires numerical/material assessment. The existing determinant
guard `J>=0.05` is retained and is not a proof of material validity over that
entire deformation range.

Inputs `youngs_modulus` and `poisson_ratio` specify equilibrium `Ginf,K`.
`relaxing_shear_ratio=Gv/Ginf` and `relaxation_time=tau` define the transient
response. The initial longitudinal modulus `K+4*(Ginf+Gv)/3` sets the new gel
wave-speed bound and inherited contact penalty. Candidate material changes
therefore require renewed contact checks.

## Exact history segment and integration stage

For a linear segment in `Edev`, let `z=dt/tau`, `a=exp(-z)`, `q=Edev-Ev`:

```
q1 = a*q0 + (1-a)/z * (Edev1-Edev0).
```

The dissipated energy is integrated analytically along that segment, including
the initial-memory/cross/increment terms. `expm1` and a small-z series avoid
cancellation. No negative-dissipation clipping or history reset is used. A
zero material-point increment leaves history unchanged.

In the PD driver, the first drift locates midpoint geometry. A predicted half
segment supplies midpoint stress and force while leaving accepted memory intact.
After the kick and final drift, a second half segment commits accepted history
at accepted geometry. Last-applied `PDStress/PDPK1/PDForceStateT0` remain distinct
from `GelAcceptedStress`. Standard `DeformationGradient` and density describe
accepted geometry after commit.

Original artificial pair viscosity is zero for this new material; physical
Maxwell dissipation is explicit. Hourglass work is still recorded separately.
Other existing material paths retain their original damping implementation.

The new driver samples contact friction with the explicit midpoint predictor
`v_n + dt/2*a_previous`, stored separately from actual velocity. Integration
binds that sampling array without changing the particle velocity state; pure
observation uses the actual accepted velocity. Raw snapshots retain base velocity,
previous acceleration and sampling velocity so the phase can be reconstructed.

Checkpoint state includes accepted and predicted internal strains, dissipated
energy densities, accepted material time, force-stage arrays, contact sample
positions/velocities and all original PD state/reference configurations. The
driver also commits the output prefix of the accepted-step trace.

## Independent verification

Framework checks cover energy-gradient stress, superposed rotation, exact fixed-
strain relaxation/energy loss, subdivision of linear strain segments, second-
order refinement of a nonlinear prescribed strain history against an analytic
solution, small/large relaxation increments, finite-deformation spatial patch
reproduction at boundary particles, and nonzero-memory checkpoint restoration.

`tools/peri_gel_validation.py` prepares a new immutable source archive, solver,
dynamic-library hashes, support tools and case inputs. It rejects reuse of an
existing dispatched state. The case sequence is:

1. CFL `0.1/0.05/0.025`, checking the *actual* initial dt ratios and complete
   accepted-step trace at the same physical endpoint.
2. Nonzero-state split/restart, no-advance state-digest equality, CSV and snapshot
   tail rollback, and bytewise split/continuous diagnostics/trace comparison.
3. Fixed material with friction regularization velocity divided by two/four.
4. Three particle spacings, fixed physical horizon 40 micrometres and fixed
   geometry, material and loading. Time resolution uses CFL 0.025.

Prospective acceptance: nonzero supported reaction curves differ by <=2%, mass
and per-particle mass errors <=1e-10, penetration/spacing <=0.1, actual complete
time/stroke support, nonnegative accumulated dissipation, consistent memory/time
cursors, unchanged frozen data, and exact restart state/output comparisons.
Space error must decrease and the finest pair must meet the force criterion.
No interpolation outside actual common stroke support and no zero-force/noise
denominator passes are allowed. Mechanical residuals are reported; no new R8
integrated-energy physical threshold is introduced.

Any execution, accuracy, integrity or resource failure stops subsequent cases.
The old stopped 48us case is not part of this sequence. Resource budgets remain
8 GiB RSS, 98304 steps, eight threads and 1000 particles for these small controls,
with a prospective 3600 s/case and 2 GiB/case output budget. New output resides on
F with the same 120 GiB host and 150 GiB output-filesystem reserves; E-hosted old
studies remain sealed. No existing source tree or old evidence is relocated.

The separate `review` process rehashes all recorded raw sizes/SHA, source archive,
ELF and dependencies. `peri_gel_audit.py` independently reconstructs energy,
accepted Cauchy stress, Green strain, kinematic density, actual last-applied
particle/box force and accepted-state contact from every available full-precision
VTP frame using Python scalar algebra and `fsum`. These sparse full frames do not
constitute independent raw energy verification of every accepted step.

## Release scope

Passing these checks qualifies the tested *small numerical control*, not an
original-size experimental capsule, an oil-bearing model or a parameter fit.
The selected gel law must still be assessed against the experimental deformation
range and mechanisms. A unique relaxation-time/viscosity fit needs a relevant
microcompression loading history; the nanoindentation time protocol and group D50
do not provide missing per-curve microcompression time or specimen diameter.
0TO calibration remains gated until the necessary numerical and experimental
mapping evidence is available.
