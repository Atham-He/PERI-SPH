# Exact contact impulse balance, guarded before physical commit

The preserved v12 infeasible control returned no accepted state, but its
opposing duals reached4.21585e8. Real-by-Real subtraction lost2.470536e-8
of normalized momentum accuracy and stopped the batch at the original1e-10
gate. This is a real trace-accuracy failure, not a passed negative control.
The old state/raw/ELF/source/reader remain stopped and unchanged.

The new candidate replaces only the arithmetic of supplied
R_i-sum_s C_si*lambda_s, followed by normalization by the supplied finite
positive impulse scale. It does not change R,C,lambda,M, material history,
the true nonsymmetric tangent, contact constraints or physical thresholds.

Every finite IEEE binary64 number is a signed integer mantissa times2^e,
with e>=-1074. Every exact product therefore has e>=-2148. A cpp_int
accumulator with unit2^-2148 holds all contact products and the supplied R
exactly, even if individual floating products would overflow. Each final
normalization is an exact integer ratio times a power of two. Integer
quotient/remainder implements nearest/ties-to-even binary64 rounding once;
no floating multiply-add or summation precision assumption is used.

The outward normalized error bound is zero for an exact result, otherwise
half the target bit unit. At the subnormal bottom it rounds outward to one
minimum subnormal. Normalized overflow is rejected. Static assertions require
IEEE binary64 Real. Current ambient rounding mode does not alter integer
decisions; the final conversion/scaling is exactly representable.

A result with normalized rounding enclosure>1e-12 is rejected before any
physical commit. At an initial or line-search point the solver returns
contact_accumulation_precision_budget and empty accepted vectors. The exact
supplied momentum, candidate velocity/duals, rounded residual and enclosure
are preserved as an explicitly UNACCEPTED precision witness, not a valid
full nonlinear iterate. A precision failure does not switch method or try
a different alpha. Every valid trace carries the per-coordinate enclosure;
root merit includes |normalized_momentum|+enclosure. Original physical
independent1e-10 and original linear backward1e-12 gates remain unchanged.

This bound covers only products, summation and normalization of SUPPLIED
finite floating inputs. It does NOT certify material evaluation, shape
gradients, tangent construction or complete physical energy. Those retain
independent checks and cannot be waived by exact contact accumulation.

New frozen primitive controls use binary bit-pattern wire inputs. An
independent Python Fraction oracle forms exact products/sums/ratios, uses
CPython's independent rational-to-float rounding and verifies the enclosure
and bit pattern. Controls include huge cancellation, product low bits,
individually overflowing products, subnormals, normal/subnormal boundary,
positive/negative even/odd ties, nonbinary normalization, changed ambient
rounding,256 shuffled extreme cancelling pairs, invalidity and overflow.
New six/twelve-coordinate nonlinear roots, row permutation, NEW infeasible
geometry and initial/trial precision rejection exercise the integrated
solver. The four previously unattempted v12 invalid inputs are allowed;
preparation rejects any SHA matching a v12 ATTEMPTED solver input.

All compiler-reported headers, including actual Boost/Eigen dependencies,
are copied and hashed before dispatch, along with source/ELF/support/plan/
inputs/linked libraries. Preparation copy budget512MiB, each native control
5s/256MiB/32MiB, single thread, Fhost/output reserves120/150GiB. Unexpected
applicable failure stops remaining dispatch and preserves exact evidence.

This qualifies arithmetic and bounded manufactured solver behavior only.
No actual PD sphere, material trajectory, energy/angular/inertia convergence
or full new history checkpoint qualification follows from it. S3 remains
active until actual PD integration and new complete restart pass; S4 and
0TO/oil release stay gated. Production gel ELF and sealed studies stay intact.
