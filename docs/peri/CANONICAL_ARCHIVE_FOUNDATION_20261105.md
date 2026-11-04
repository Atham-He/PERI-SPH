# Detached canonical archive byte layer

pd_canonical_archive.h is a NEW unused byte-format component. It is not included
by the active v18 executable or its readers, and does not alter any v18 frozen
dependency. It is not integrated into the model, output stream or restore
transaction. No complete physical checkpoint qualification is claimed.

The component stores explicit typed bytes, little-endian uint64 and original
IEEE-binary64 bits. Shape/type/name/count/size bounds precede candidate use.
Canonical field order is lexicographic, independent of registration order.
It preserves negative zero, minimum subnormal and extreme finite scalars.
It never serializes Eigen/vector object memory or binds decoded fields to live
simulation pointers. Required schema names/types/shapes are explicit.

The supported digest profile is Linux AF_ALG hash/sha256, using the installed
kernel provider. It has no alternative algorithm/provider fallback and fails
explicitly if unavailable. Known abc/empty vectors and streamed1MiB content are
checked; a separately implemented Python struct/hashlib parser verifies the
canonical fixture. No cryptographic hash algorithm is reimplemented here.

The native engineering test covers deterministic registration order, exact
scalar bits, detached noadvance bytes, explicit-generation mismatch, exact
schema/type/dimension/name refusal, duplicate fields, nonfinite data, corrupted
hash/footer, and valid-checksum structural/nonfinite corruption. These are
byte-format checks only, zero physical steps. Domain-invalid but byte-valid
material/reference/cache/contact/history values require the future full
session adapter's semantic tests.

Still required after actualPD friction qualification: complete source/ELF/
config/model/reference/neighbor/cache/material/contact/dual/algorithm/
compensation/acceptedclock/output-cursor binding; empty temporary receiver;
all read-only semantic validation; private output branch and atomic commit;
full noadvance/split/rollback bytes and valid-checksum semantic-invalid
no-live-mutation. Old DiskCheckpoint auto-fallback/retention is not reused.

Do not add this header to or rebuild the active v18 targets. The standalone
byte engineering test is compiled separately with only standard/Linux headers,
with its own bounded artifacts. CMake integration and actual checkpoint
trajectories wait for v18 completion. S3,S4 and0TO remain gated.
