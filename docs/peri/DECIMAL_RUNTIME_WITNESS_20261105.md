# Failed-reader arithmetic observations

The v19 native trajectory completed, but its first full independent reader
raised `RuntimeError: internal error in flags_as_exception`. The failed review
and terminal remain unchanged. No physical tail is dispatched. Instrumented
read-only replay preserves the original 95-digit operations and all gates;
it does not retry an operation, clear traps, change precision or call a solver.
It records the step, material point, finite-difference direction, exact operands,
bytecode/source columns and context if the original exception recurs. Absence
of recurrence is not an automatic positive study verdict.

The separate minimal-witness utility uses bounded integer/Fraction arithmetic
for nearest-even significant-digit rounding and can compare a captured scalar
operation with C and pure-Python decimal implementations in NEW declared clean
contexts. This is an independent scalar witness, not an alternate solver or
automatic runtime fallback. It cannot explain the original error merely by
obtaining a correct result in a clean context. It reports private unknown
signal bits explicitly; it never writes or repairs the original context.

Context-word interpretation is restricted to the documented 64-bit
little-endian 96-byte layout, guarded by agreement with all observed public
fields and signal mappings. An unsupported or inconsistent layout is rejected,
not guessed. Arithmetic source is parsed with a small AST whitelist, never
evaluated as code. Nonfinite, excessive-size and exponent-boundary operands
are outside this bounded witness profile. Whole-trajectory, material, contact,
energy, negative-control and restart qualification remain separate requirements.
