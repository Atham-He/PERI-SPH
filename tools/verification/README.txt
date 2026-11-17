Physical fracture qualification (17 November 2026)
=================================================
Canonical source: /home/hezhy58/SPHinXsys
Canonical build: /home/hezhy58/SPHinXsys-build

Build:
  cmake --build /home/hezhy58/SPHinXsys-build --target pd_physical_crack_benchmark -j4
Run:
  /home/hezhy58/SPHinXsys-build/pd_physical_crack_benchmark <output-directory>
CTest registration: pd_physical_precrack_open_close_slip

This is ONE 768-point layered NOSB-PD specimen with a prescribed precrack.
The same geometry is opened, closed, slid and rigidly rotated; a physical-mass
undamped Verlet hold checks force/energy consistency. No material fitting,
particle erosion, displacement amplification, or altered mass is used.

Validated implementation scope:
* Reciprocal material bonds are actually removed, and broken gradients are zero.
* Layered live-support material gradients and live polynomial stabilization
  are reconstructed; insufficient quadratic support uses an explicit qualified
  affine fit, never a hidden residual spring across the crack.
* Crack-face trace interpolation uses local live connectivity plus geometric
  visibility; topology epochs invalidate material, stabilization and trace caches.
* Known opposite interior faces have frictionless unilateral contact. The
  contact force is the exact transpose/energy derivative including moving normals.
* The physical precrack PVD contains unscaled positions and explicit face polygons.
  PVD time values label prescribed modes, not a physical loading-time history.
* A separate finite-K cohesive law checks the retained ft/Gc work budget. This
  does NOT validate dynamic insertion or growing fractures in a capsule.

Limits that must remain explicit:
* Dynamic initiation/growth and a work-consistent continuous-to-discrete handoff
  have NOT been qualified in the capsule. Precrack passing cannot clear that gate.
* Damage from the legacy directional continuum law is not BrokenBondFraction.
  It is invalid to convert point Damage=1 directly into removed material bonds.
* The contact qualification covers known interior face pairs, not arbitrary
  fragment edge-edge collision or a global remeshing/contact algorithm.
* Cuts through a thickness column require split-column kinematics and currently
  fail explicitly; no column or particle is silently deleted to hide this limit.
* The standalone fracture checkpoint is not a replacement for the capsule's
  complete mechanical/material/time-integrator checkpoint.

The current capsule is deliberately stopped. Do not restart the old diffuse
run from this precrack test or migrate its post-damage checkpoint to broken
bonds by relabeling fields, bypassing source identities, or paying Gc twice.

Literature motivating the checks (not claims of literal implementation):
Lu et al., Theor. Appl. Fract. Mech. 118 (2022), 103082,
  doi:10.1016/j.tafmec.2021.103082 (crack contact and precracks).
Foster, Silling & Chen, Int. J. Multiscale Comput. Eng. 9 (2011), 675-687,
  doi:10.1615/IntJMultCompEng.2011002407 (state-based critical work).
Doyen, Ern & Piperno, Comput. Mech. 52 (2013), 401-416,
  doi:10.1007/s00466-012-0819-2 (quasi-explicit cohesive integration).
The latter two do not by themselves justify using spherical 3-D critical-work
normalizations for this quadrature-resolved curved thin layer.
