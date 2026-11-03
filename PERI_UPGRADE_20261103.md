# SPHinXsys official update merged into the existing Peri checkout

The active source directory is `/home/hezhy58/SPHinXsys`, branch `capsule`.
The official baseline is `1612330e44c8947406f24b6756449e80dcd7158e`.
The existing Peri history, local cantilever changes, PD module, application
examples, checkpoint implementation and audit tools have been retained.
The original Git directory and registered frozen validation worktree remain.

## Recovery record

Before the in-place merge, the original source tree was saved in:

`/mnt/e/PERI-SPH/upstream_upgrade_20261103/inplace-backup-20261103T145631Z/original-source.tar.gz`

Original commit: `9d9ccd5e68f1b53315b90502467915028513dd94`.
Recovery branch: `backup/capsule-before-latest-20261103T145631Z`.
Retained local-change stash: `e75bfac7a0702301b7c69c5ecd0f4d06ea575851`.
The archive includes tracked source and retained untracked/ignored Peri sources;
the original build/results/checkpoint tree remains in place and was not archived
or replaced. Do not resume the sealed failed studies as part of source recovery.

The Windows-visible import record directory is
`E:/PERI-SPH/upstream_upgrade_20261103`.
Read `inplace_merge_receipt.json` for final identities and checks.

## Build in this existing directory

```sh
cd /home/hezhy58/SPHinXsys
cmake -S peri -B build/peri-latest-20261103 -G Ninja -DCMAKE_BUILD_TYPE=Release
cmake --build build/peri-latest-20261103 --parallel 12
ctest --test-dir build/peri-latest-20261103 \
  -R '^(extra_3d_PD_framework_tests|peri_analysis_unit_tests)$' --output-on-failure
```

The isolated build-output subdirectory prevents old validated binaries from
being replaced. The source is the existing checkout, not an independent source
copy. Long capsule, brittle-disc and ice-jet CTests are disabled by default.
The imported cantilever analysis additionally requires NumPy/SciPy/Matplotlib;
see `peri/requirements-analysis.txt`. Its five unit checks passed using the
existing Windows monitoring Python environment. The system Ubuntu Python does
not have NumPy installed.

Building and passing unit checks does not transfer the previous scientific
validation results to the new executable. Original frozen R8 source/binaries,
failure receipts, checkpoint data and the stopped queue retain their locations.
No production simulation or original queue is started by this migration.
