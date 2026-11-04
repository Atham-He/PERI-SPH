# Nonzero bounded body archive component

This layer does not qualify a complete session restart. It stores the actual
33-node model, reference neighborhoods, cached shape/gradient/nonaffine operator,
control and dense signed MLS maps, contact frame and ordering, accepted material
and contact histories, cumulative sums AND corrections. Decoding constructs an
empty owned receiver. Validation reads its stored values and computes expected
values in scratch; it never repairs caches or reconstructs a fresh scene.

The independent byte oracle starts from the preserved v19 40-step final raw,
not an empty initial fixture. Python struct/hashlib writes every field using
explicit little-endian values and dimensions. Native decoding and re-encoding
must match the whole file bit for bit. Every listed semantic corruption has a
new valid SHA. The native decoder must refuse the specific semantic cause and
leave the original owned baseline bytes unchanged. Python independently parses
each corrupted envelope, checks its valid digest and exact changed-field set.
No material/contact advance or physical solver call occurs.

Engineering gates before execution: all bytes/fields/count companions match;
18 cause-specific negative controls; nonzero viscous/tangent/travel-correction
history; no reference refresh; no physical advance. Build cap180 CPU/240 wall,
2GiB AS/RSS; component checker30 CPU/60 wall,1GiB AS/RSS,one thread,256MiB
file/1GiB total outputs. Host F>=120GiB, output reserve>=150GiB, available
memory>=8GiB. Freeze sources, actual compiler dependencies, ELF/libraries,
input, plan and checker before the single engineering execution. An applicable
failure preserves its evidence and stops that state.

Full session work still must pin protocol, solver/algorithm contract, complete
committed output prefix and generation identity; validate in detached private
storage and atomically publish afterward. New physical noadvance,
split-continuous and tail rollback must compare all bytes, including outputs.
This component result cannot substitute for those tests, full friction law
qualification, S4 convergence or material calibration. The original v19 C
decimal failure and separate fixed pure-Python raw audit remain unchanged.
