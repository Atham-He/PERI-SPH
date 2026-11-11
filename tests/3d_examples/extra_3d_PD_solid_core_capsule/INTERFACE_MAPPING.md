# Optional quadratic interface trace

The reference bonded interface evaluates the displacement of each volume body
on the same physical interface, then scatters one penalty force through the
transpose of each evaluation map. The default remains affine MLS.

`interface-order=2` enables complete degree-two reconstruction on **both**
interface traces. `interface-core-support-ratio` and
`interface-shell-support-ratio` specify support radii relative to each body's
spacing; both default to 2.5. This does not change the outer plate-contact map,
the constitutive laws, reference volumes, shell thickness, or penalty stiffness.

On the 25 micrometre spherical core grid, a 2.5-spacing support near an axial
pole includes only two normal coordinate levels. It is rank deficient for a
complete quadratic polynomial. Use a resolved support, such as 3.0 for this
core; the three-layer PD shell supports the fit at 2.5. The implementation
whitens reference coordinates and rejects rank loss. It never silently
regularizes an incomplete fit. This support choice is geometric and is not a
material parameter inferred from experimental reaction force.

The actual signed weights enter the existing transpose force scatter and the
absolute-row-sum mass/stiffness time-step bound. Do not retain a time-step bound
from a different reconstruction. A complete polynomial patch check, a loaded
virtual-work check, and a physical loaded comparison serve different purposes:
passing the first two does not validate the complete capsule response.

This option assumes an intact fixed reference support. Fracture, interface
debonding, slip and particle redistribution need separate support/state logic.
