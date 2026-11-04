# Gel sphere contact and spatial verification

The first independently reviewed gel study at commit 824b7b6ed passed time
refinement and full-state restart, but the 15-to-10 micrometre particle-spacing
force difference was 41.2119489%, above the predefined 2% limit. Its failed
state, source archive and review remain immutable. No calibration is released.

The old dp/2 particle-radius contact starts on a flat lattice layer. Initial
volume quadrature errors at spacings 20, 15 and 10 micrometres were respectively
20.2504%, 4.44543% and 0.798131%. This motivates a separate contact discretization;
it does not establish that contact alone explains the spatial error.

## Opt-in reference surface

`reference_surface_contact=1` enables massless samples of the full reference
sphere. Gauss-Legendre polar integration uses mu=1-z^2 on each hemisphere, and
uniform azimuths. Sample reference areas are positive and sum to 4*pi*R^2.
They add neither particle mass nor a dummy liquid.

Fixed compact reference MLS weights reproduce constants and affine fields.
Signed weights are retained. Position and velocity use x_s=sum(phi_si*x_i).
Contact force uses the identical transpose, f_i=sum(phi_si*f_s), preserving
virtual work, global force and torque. Rank-deficient support is a hard error.
The map is not refreshed on a deformed body or after checkpoint restore.

For a plane with inward unit normal n and permitted region n.x>=a,
delta=max(a-n.x_s,0), U_s=0.5*(K_contact/ell)*A_s*delta^2. Normal force is the
negative position gradient of U. Regularized Coulomb friction acts on the
same surface sample. `ell` is a declared numerical penalty length, not a fitted
force multiplier. It needs its own refinement before experimental fitting.
Reference-area integration is explicitly total Lagrangian; finite-strain
contact behavior still needs qualification. A sample escaping the actual
finite plate support is rejected rather than treated as contact with an
infinite plane.

The integration uses midpoint velocity prediction. Pure accepted-state
observation does not update forces, material history or the reference cache.
A conservative mapped active stiffness/friction bound supplements the PD
wave and relaxation-time step limits. Every actual accepted dt is recorded.

Checkpoint state includes reference IDs, offsets, weights, areas, positions,
generation, active wall masks, last sample positions/velocities/reactions and
the applied-time cursor. The reference CSV is verified against the restored
cache; it is never silently regenerated. Independent Python arithmetic gathers
and scatters all stored samples using fsum and verifies force/potential against
full-precision particle fields. Sparse frames remain insufficient for an
independent full-step energy claim.

## New prospective cost budget

This new solid gel spatial study permits 8192 particles. The measured
912-particle, 920-step case took 2.764 s with a 9.94 MiB checkpoint. At fixed
40 micrometre nonlocal horizon, halving particle spacing predicts roughly
128 times the work and 64 times neighbor storage. The case budgets remain
3600 s, 8 GiB RSS, 98304 steps and eight threads; output budget is 4 GiB.
Actual host/output reserves remain 120/150 GiB on F, where new data reside.
The sealed original R8 N<=1000 limit and all old studies remain unchanged.

The new frozen study compares time steps, friction regularization, 32/64/128
surface rules and dx=10/7.5/5 micrometre particle spacing at fixed 40 micrometre
horizon. Force criterion remains 2%. Failure stops the sequence and retains
all bytes. Fixed-horizon particle convergence is not horizon convergence.

## Experimental assumptions

The user has no original loading-speed or individual-diameter record. The
prospective speed is 0.5 mm/s, with 0.25/0.5/1 mm/s sensitivity, based on the
single alginate microbead protocol of Essifi et al., DOI
[10.1016/j.matpr.2020.05.747](https://www.sciencedirect.com/science/article/pii/S2214785320343601).
It is not attributed to the target experiment. Its 2.34 mm group D50 is not an
individual specimen diameter. See the versioned assumptions JSON. Relaxation
times inferred later will be conditional on assumed speed/geometry. The current
small numerical controls at 0.1 m/s are not experimental calibration runs.
