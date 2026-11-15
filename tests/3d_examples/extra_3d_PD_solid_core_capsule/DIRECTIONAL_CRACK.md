# Directional crack-band closure for the layered NOSB-PD shell

The optional `fracture-kinematics=directional_contact` mode adapts Definition 2.1
of Hakimzadeh et al., JMPS 167 (2022) 104994,
https://doi.org/10.1016/j.jmps.2022.104994 . This is an effective finite-width
crack-volume contact model, not a sharp fragment-to-fragment search algorithm.
The default `volume_deviatoric` mode is unchanged for existing callers.

The intact finite-volume Biot energy and its E and nu are unchanged. In a frozen
reference crack frame N, QR normal-column shear modes are relaxed; the normal
height is relaxed only above its traction-free elastic value. Tangential metric
and normal closure energy remain. The blend is (1-d) W(F) + d Wd(F,N). Exact local
first and second derivatives yield PK1 and a conservative absolute-Hessian row
bound for the explicit time step. Crack normals are frozen on first damage.
There is no residual strength, mass scaling, force clipping or particle erosion.
The existing scalar irreversible Gc crack-band history remains. Its initiation
energy is derived from the new release W-Wd at the same prescribed strength;
E, nu, ft and Gc are not fitted or rescaled.

The quadratic residual stabilizer statically condenses a point-centered
Heaviside jump mode, retaining the orthogonal residual energy. This is our
energy-consistent adaptation, not code copied from the JMPS paper. Stabilizer
release enters the same damage budget. A globally sharp jump is not exactly
represented by all point-centered planes: the benchmark reports its nonzero
residual rather than claiming exact sharp-fragment separation. The closure law
is frictionless. General sharp-fragment collision and geometric nonpenetration
have not been validated by this test. The intact Biot energy is not globally
rank-one convex; leaving the stable zero-shear relaxation branch raises an
explicit error instead of introducing an undocumented constitutive fix.

One 480-point layered crack specimen is checked in opening, finite slip, closure,
closed slip and reopening. Checks include analytic closure energy, zero opening
and slip tractions, frame indifference, exact transpose force/energy and moment
balance, local jump projection, retained non-crack stabilization, prescribed Gc
work, and 200 undamped physical-mass velocity-Verlet steps. Zero-response balance
errors are normalized by the physical E-volume scale to avoid dividing roundoff
by roundoff. The test is not a spatial convergence or capsule validation study.

Build target: pd_directional_crack_benchmark.
CTest: pd_directional_crack_open_close_slip.
The production acceptance and build receipts are archived under
framework_work/2026-11-15_SS30_release_reaction/directional_fracture_20261115.

All available full SS30 checkpoints already contain old-law damage histories.
The corrected same-input capsule therefore starts from the initial state in a
new directory; no old damage or release/reaction histories are converted.
