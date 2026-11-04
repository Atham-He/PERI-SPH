# Output-prefix byte foundation

The bounded PD campaign remains frozen. These standalone files are not in its
compiler or runtime dependencies and do not change a physical trajectory.

The canonical archive and output-prefix reader share a streaming AF_ALG
SHA256 implementation. A 64 KiB buffer permits an explicit selected prefix
up to the existing 4 GiB output bound without allocating the entire physical
stream. A later tail does not select a newer checkpoint. An exact snapshot
must have the selected byte length and hash; regular-file, no-symlink and
unchanged-source checks apply. Private staging creates only a new exclusive
leaf and fsyncs that file. It cannot accept or overwrite a live output file
descriptor. Partial staged evidence is retained on failure.

Independent Python hashlib/struct checks are required before using this byte
foundation. Tests distinguish an old selected prefix from a later tail,
corruption/truncation/wrong cursors, existing or symlink destinations, unsafe
names, source changes after proof, and digest lifecycle errors. Refusal must
have its expected cause, not merely any exception. Protected sentinel files
and original source snapshots must remain byte-identical.

This layer does not parse a physical record or validate a material, contact,
dual, cache, algorithmic history or accepted output cursor. It does not switch
an active branch pointer or implement the live owner transaction. Full CP
integration still requires empty temporary model receivers, complete semantic
validation, explicit-generation selection, private output staging and tested
commit semantics; noadvance, split-continuous, tail rollback and valid-hash
semantic corruption refusal must be qualified on the actual new PD algorithm.
Zero physical steps in these byte tests cannot qualify S3, S4 or calibration.
