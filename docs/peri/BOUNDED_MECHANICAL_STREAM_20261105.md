# Complete physical prefix and scratch proposal validation

New codec:48byte magic/version/trusted-contract header, then initial state and
one explicit-length canonical138field record per accepted step. Each contains
all stored model/reference/cache/material/contact/dual/compensated history,
both plate normal/tangent/total impulse, actual-interval mean reaction and
work, trueM kinetic/material/nonaffine energies, separate physical Maxwell,
normal contact/friction and algorithmic losses, physical residuals and full
cubic-wall interval/end minimum gaps. Initial kind0 has no interval; its
observation fields are explicit zero sentinels. No PID/path/walltime or split
end marker occurs in this deterministic physical stream.

The same verified file descriptor supplies bounded semantic reads and later
exclusive private staging. Reopening a path between hash and parsing is not
used. Explicit selected bytes may exclude a preserved later tail; immutable
snapshots must have exact length. Per record<=16MiB,total<=4GiB. All
record schemas/generations and final archive/prefix body bits must match.

Refactor only: native solve and readonly validation share finishProposal,
which is the original accepted-state/history/ledger arithmetic in the same
order. Native advance still runs the original nonsymmetric solver. Prefix
validation does NOT run that solver: it builds context from stored arrays,
checks saved endpoint original R/KKT/precision, then recomputes the proposal
in scratch. Full encoded body/observation equality checks material history,
all contact sums AND compensation, each plate reaction/work and all losses.
No accepted arrays or reference caches are repaired/replaced. Contact and
geometry helpers now accept a readonly borrowed stored-scene view; their
numerical expressions remain unchanged.

Engineering input is converted EXISTING raw, not a new native physical case.
All41records/40 accepted observations are independently encoded/parsed by
Python struct/hashlib, with JSON -0 token sign preserved. Previous static/
accepted/friction/reaction packets were inspected: no negative zeros in
those saved fields; the prior body-component result is not invalidated.
Frozen checker source/ELF identities describe this engineering reader, not
the old physical producer. No old verdict or independent review is changed.

Prospective single test:2positive selected-prefix cases,20valid-SHA semantic
refusals+2integrity+3external-bound/staging controls; baseline/copy bytes and
preserved later tail identical. Native budget900CPU/1200wall/1GiBAS/1thread;
build180CPU/240wall/2GiBAS. Cost basis: old native40step205.775s includes all
Newton; this checks3full readonly walks plus early controls without Newton.
256MiB/file,1GiB all outputs, F reserve150GiB/host120GiB,memory8GiB. Plan,
source/actual includes/libraries/ELF/inputs freeze before the unique execution.
An applicable failure stops that state with exact exception/record cursor/
returncode/evidence, no solver rerun/retry/fallback or gate relaxation.

This is still NOT a live generation transaction or physical restart proof.
Next qualify actual F fsync/atomic publication, detached candidate/private
branch ownership and explicit commit with uncertainty handling. Then NEW
physical noadvance/split/tailrollback full state AND physical output byte
tests, gated by complete independent friction evidence and untouched tails.
S3 remains active; S4 convergence and0TO remain gated.
