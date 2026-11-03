# Upstream sync policy (avoid fork drift)

## Independent upstream import on 2026-11-03

The current checkout imports official master `1612330e` while retaining the
local Peri working tree from `9d9ccd5e`. See `PERI_UPGRADE_20261103.md` and
the `E:/PERI-SPH/upstream_upgrade_20261103` evidence directory. The new `peri/`
build entry and capsule/framework tool path adjustments are build integration
changes. Original validated experiments retain their original source/binaries.

This fork follows a simple two-branch model:

- `master` — tracks the upstream `Xiangyu-Hu/SPHinXsys` master as closely as
  possible. Never commit our own work here.
- `capsule` — our development line: the peridynamics module, its core patches
  and the PD application examples.

## Why the discipline

Industrial forks of SPHinXsys commonly drift thousands of commits behind
upstream (e.g. 17 ahead / 5473 behind), making later merges painful. Our core
patch surface is intentionally tiny — keep it that way:

1. `src/shared/particle_neighborhood/neighborhood.{h,cpp}` — per-bond
   `bondLive_` state (+10 lines).
2. `src/shared/include/sphinxsys.h` — module registration (+1 line).

**Rule:** new physics goes into `src/shared/peridynamics/` or example
directories. Touching core files outside the two spots above requires an
explicit note in this file.

## Routine sync (suggested every 1–2 months)

```bash
git remote add upstream https://github.com/Xiangyu-Hu/SPHinXsys.git   # once
git checkout master
git fetch upstream
git merge upstream/master        # fast-forward expected
git checkout capsule
git merge master                 # conflicts only possible in the 11 core-patch lines
ninja -C build extra_3d_PD_brittle_disc extra_3d_PD_capsule_compression   # smoke build
```

After merging, run one short case (e.g. the brittle-disc smoke case) before
resuming production runs.

## Upstream PR candidates

When the module matures (paper published), consider contributing back, in
order of increasing scope:

1. `bondLive_` per-bond state in `neighborhood` (generic, small, reviewable).
2. `ContactRepulsionFromWallforPD` smooth-Coulomb friction (self-contained).
3. The whole `peridynamics` module + one validation example
   (the official contribution path: feature + test case + PR review).

Contributing back shrinks our patch surface and lets the paper cite the
upstream repository directly.

## Explicit core I/O extension (2026-10-15)

`src/shared/io_system/io_vtk.{h,hpp,cpp}` adds the opt-in
`setFullPrecisionOutput(true)` setting. Fixed nine-decimal ASCII output rounded
the 1 mm capsule's particle masses (~9.35e-13 kg) and small forces to zero.
The setting writes scientific notation with `Real::max_digits10`, Float64
declarations and double VTK arrays/points. It changes no dynamics or default
upstream output behavior. Capsule snapshots opt in; scalar, vector, matrix,
coordinate and physical-time serialization have round-trip regression tests.
Keep this I/O fix separate from the PD physics when syncing upstream. The
current workstation validates the ASCII build; VTK-library builds require
their own build check.
