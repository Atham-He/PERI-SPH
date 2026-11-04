# Explicit generation publication and owned checkpoint state

The store owns exactly one checkpoint state and its working output handle.
Callers cannot supply or swap another live owner. A nonblocking exclusive
cooperating-writer lock is retained for its lifetime. Root/branch/file access
is directory-FD relative, validated safe leaves and O_NOFOLLOW; root ancestors
and trusted contract must come from the frozen driver, never archive paths.

Only explicit generation selections are parsed. HEAD bytes are an expected-
current compare guard, not an instruction to search latest or fallback.
Selectors canonically contain branch, archive/prefix/contract digests and
generation, independently checksum protected. Selected immutable state and
prefix are fully verified before publication. Archives never rebuild a scene.

Candidate validation, exact prefix-proof binding and source unchanged checks
precede private branch creation. New exclusive branch retains immutable
state.bin/prefix.phy plus a separate active.phy copy. All allocations/open/
hash/readback/fsync complete first. Save previous pointer bytes,fsync pending
pointer and directories,rename HEAD atomically,then swap unique_ptr ownership
without throw. Previous branches and pending failure evidence are retained.

Pre-rename failure throws without HEAD/live change. A post-rename directory
sync failure returns an explicit published-but-durability-uncertain result:
new HEAD/new owner are known,old/new branches preserved,no false unchanged or
silent rollback claim. That owner blocks all further publication until an
external explicit recovery creates a new process/store. Fault injection is
only a declared engineering control and is reported as injected; it is not
an actual filesystem/power-loss incident.

This class has NO accepted-step advance/write API yet. Physical execution
must later bind the ACTUAL running solver ELF/manifests to the contract,
prepare candidate+record+hash before durable append,and perform a noexcept
accepted cursor/state move only after fsync. Converted-data store ownership
tests cannot substitute for that or full physical restart qualification.

New engineering data reuses the already-converted valid41record prefix,
not a physical solver. First process checks19specified refusal controls,
2actual durable publications and1controlled post-rename uncertainty.
Second distinct process explicitly restores37 while current HEAD38,retains
the newer branch and a SYNTHETIC uncommitted abandoned-working-file tail,
publishes another37 branch with whole archive/prefix bytes identical.
No actual accepted PD tail, nonlinear solve or material advance is inferred.
Python struct/hashlib independently checks every pointer/archive/prefix and
the old evidence before and after the second process. Symlink targets and
pending/previous pointer evidence stay intact. No recursive delete/move.

Prospective budgets:build180CPU/240wall/2GiBAS;each of2native phases
180CPU/240wall/1GiBAS/1thread;256MiB/file/1GiB all outputs;F output150GiB/
host120GiB,memory8GiB. Cost reference:prior3prefix walks+25controls22.487s.
Freeze source,actual compiler dependencies,ELF,inputs and plan before first
execution;freeze concrete phase1 outputs before phase2. Unexpected failure
preserves exact evidence/stops later phases,no retry/fallback/gate relaxation.
Actual F fsync/lock/rename support is required,never silently weakened. This
does not claim a power-loss test,physical restart or whole S3 completion.
