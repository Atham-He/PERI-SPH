# Bounded actual PD checkpoint driver: implementation and prospective qualification

This driver calls the original nonsymmetric Coulomb solve exactly once per
physical proposal. It does not symmetrize the tangent or change true mass,
material, prescribed normal speed, tolerances or contact penetration gate.
The immutable original probe ELF and production gel ELF are not rebuilt.

`--initialize ROOT` constructs one independent 33-node initial reference and
emits the canonical initial body, contract and physical record. The externally
computed full contract digest and every artifact must be frozen before any
physical mode. The zero-advance preflight is not a physical qualification.
The shared observer extraction preserves the original full Newton/point tangent
record. The fresh physical mode emits its original t0 header using stored data.
Checkpoint publication metadata goes to its separate immutable selection files,
not extra JSON rows in that original identity/header/steps/end debug protocol.
Resume has an explicit resume_header and requires the preregistered branch trace
comparison/import policy; it never pretends to be a fresh t0 full trajectory.

Commands use absolute ROOT/STORE paths. ROOT supplies frozen original.input,
source_manifest.json, build_runtime_manifest.json and initialized files.

* `--fresh ROOT STORE TRUST TARGET LABEL`: STORE must have no selected HEAD;
  loads the frozen independent t0 packet, publishes it, and advances to TARGET.
* `--resume ROOT STORE TRUST SELECTION_LEAF GENERATION EXPECTED_HEAD_SHA TARGET LABEL`:
  loads only the explicit generation with complete prefix proof before staging.
* `--validate ROOT STORE TRUST SELECTION_LEAF GENERATION EXPECTED_HEAD_SHA`:
  same detached validation, zero new steps, no branch or live owner publication.
  It uses the cooperating store lock, and may create its empty lock file.
* `--restore-noadvance ROOT STORE TRUST SELECTION_LEAF GENERATION EXPECTED_HEAD_SHA LABEL`:
  actually installs the restored owned body, then saves a new complete branch
  without any solve or clock advance. The explicit logical generation is kept,
  even if it differs from the physical step. Archive and entire physical prefix
  must remain byte-identical. This is the required positive noadvance operation;
  the read-only validator alone is insufficient to qualify owned restoration.

The new zero-advance preflight uses explicit logicalG7 at physicalstep0 to check
the installed/re-saved path and prevent implicit generation=step retagging. It
does not replace required NONZERO press/hold/release physical restart tests.

Resume never constructs a fresh scene, recomputes reference or initial material
history, or guesses another generation. The entire old body is owned across
the solve/commit/debug path; no borrowed view can dangle after pointer swaps.
Each accepted physical record is checked against the original endpoint equation,
all histories and observations in scratch, then appended/read back/fsynced before
the noexcept owned-state swap. That endpoint check does not call a second Newton.
Actual mechanical durable cursor is updated before full debug formatting.
Failures print the mechanical cursor, separate current-process full debug-row
cursor, actual solve count and first append failure witness while owner lives.
Checkpoint publication errors also retain the explicit error number and stage,
including the selected new owner after a post-rename directory-sync failure.
`--diagnostic-fixture` is a zero-advance controlled formatter check, not a disk
fault or physical failure; its distinct mechanical/debug cursors, numeric error
and escaped quote/backslash/newline are independently checked.
No subsequent solve occurs after an exception, rejected gate or uncertain write.
CP sealing occurs at accepted steps8/20/32/40 and any declared early terminal.
All branch prefixes, selected pointers and unselected later tails are retained.

Before physical dispatch: the current fixed-reference full raw audit must pass;
then the two never-attempted friction tails each require native plus independent
review. Freeze a NEW CP campaign and actual resource/cost budget. A proposed
new two-axis input is dt=.01, normal stroke=.01, tangent(.0006,-.0004), mu=.22;
it is not a replay of old attempted native input. Continuous40 must independently
pass all original pointT/R/J/trueM/history/geometry gates and negatives before
dependent noadvance at8/20/32, split8+32 and old20->40 tail rollback commands.
Full accepted archives and committed physical streams must match byte-for-byte.
Correct-hash semantic corruptions must reject before staging/live mutation.
Retain all full Newton trial raw. New physical S3/S4/calibration remain gated.
The old 1% trajectory algorithm-loss target remains for S4; it is not proved by
the archive or by this new driver implementation.
