# Two-component Coulomb disk kernel, bounded S3 precursor

The normal Signorini law and two tangential maximum-dissipation conditions
remain separate simultaneous equations. The dimensionless tangential
residual is t - projection_disk(mu*max(l,0))(t-w), with l=lambda_n/I,
t=tau/I, w=(G_t*v-v_wall,t)/speed. A full associated cone projection would
alter the normal law; independent component clipping would replace the
isotropic disk by a square. Neither is this candidate's constitutive law.

pd_coulomb_disk.h computes the original disk map, normal-radius coupling,
2x2 impulse and relative-velocity derivatives, physical friction loss
-tau dot relative_velocity, cone violation and maximum-dissipation error.
At the circle it selects the exterior limiting derivative. At the
zero-radius origin it selects a limiting derivative from negative normal
trial impulse. The final normal solve must still enforce lambda_n>=0.
There is no numerical smoothing or penalty spring in this kernel.

Binary64 inputs are exact in the 256-bit binary reference arithmetic.
Projection, derivative and work expressions use a conservative 2^-240
absolute-magnitude envelope, leaving 16 guard bits over fewer than 64
operations per expression. Positive squared norms contain no cancellation;
near binary64 subtraction is exact at this precision, distant subtraction
has a small relative error, and projection remains continuous across its
branch boundary. Conversion back to binary64 adds the explicitly measured
conversion difference and rounds the enclosure outward with nextafter.
Derivative enclosures concern the selected branch; they are not a theorem
of global nonlinear convergence at every nonsmooth input. Frozen independent
95-digit value/derivative/conversion tests are still required. No valid
arithmetic candidate may be consumed when its representation or 1e-12
roundoff budget fails; finite placeholder fields in that witness are not
physical or accepted values.

The new native probe only evaluates a fixed input. It saves actual
PID/parent/boot/startticks/argv, inputs, original values/derivatives/guards,
and zero accepted physical steps. A resource, unexpected execution,
integrity or accuracy failure stops remaining inputs. The expected invalid
and arithmetic-range rejections are declared before dispatch. An independent
Decimal reader checks exact binary inputs, radial maximum dissipation,
conversion enclosures, independent centered-difference derivatives away
from corners, selected limiting branches, basis rotation, and deliberate
corruption of cancellation handling, disk shape, normal coupling, loss sign
and precision acceptance. Existing physical/raw studies are not rerun.

This is only a required kernel precursor. It does not qualify a coupled
normal/friction solve, a true PD friction trajectory, or complete restart.
The next global system must retain all normal/tangent coefficients from the
same signed MLS map and exact transpose, true nonsymmetric material J and
physical M. Current normal-only linear API has a 384-unknown bound; 33 actual
nodes and 116 three-component contacts require 447 unknowns, and a same-map
zero-mu comparison against the sealed114-sample normal raw would require
783. Introduce a separately versioned bounded API with original-system
linear residual checks before those cases; do not silently reuse or enlarge
the old qualification, symmetrize J, or replace it by physical mass.

After global/actual PD friction qualification, a new complete transactional
checkpoint must validate material/reference/neighbors/shape/HG/surface/
contact/normal-tangent dual/algorithm/slip-loss and committed output cursors
in temporary state plus staged outputs before swapping live state/output.
The old DiskCheckpoint automatic fallback/retention behavior is not the
new qualification: use explicit committed generation selection, preserve
all qualification evidence, and reject corruption without live mutation.
Full-byte noadvance/split/tailrollback and semantic-corruption tests are
required. S3 remains active; S4 force/time/space/horizon/large-deformation/
inertia/rate and the existing 1% algorithmic-loss target precede 0TO fitting.
