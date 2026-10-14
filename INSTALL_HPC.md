# INSTALL (HPC) — SPHinXsys with the NOSB-PD module

This repository is a fork of [Xiangyu-Hu/SPHinXsys](https://github.com/Xiangyu-Hu/SPHinXsys)
extended with a non-ordinary state-based peridynamics (NOSB-PD) module
(`src/shared/peridynamics/`) and PD application examples
(`tests/3d_examples/extra_3d_PD_*`).

**The installation procedure is identical to vanilla SPHinXsys** — the module
is compiled into `libsphinxsys_3d` automatically by CMake globbing; no extra
build flags are required.

## 1. Dependencies (same as vanilla SPHinXsys)

- CMake ≥ 3.16 (< 4.0), C++17 compiler (GCC ≥ 8 on Linux), Python3, Git
- Libraries: Eigen3, Intel TBB, boost-program-options, boost-geometry,
  Simbody, GoogleTest, spdlog

On a bare Ubuntu (20.04+) login node:

```bash
sudo apt update
sudo apt install -y build-essential curl zip unzip tar pkg-config git cmake ccache python3-dev gfortran
sudo apt-get install -y autoconf automake autoconf-archive libtool

cd $HOME
git clone -b 2026.04.27 https://www.github.com/microsoft/vcpkg.git
cd vcpkg
./bootstrap-vcpkg.sh
./vcpkg install --clean-after-build \
    eigen3 tbb boost-program-options boost-geometry simbody gtest spdlog
```

On an HPC cluster with environment modules, the equivalents are usually
`module load eigen tbb boost cmake gcc` — Simbody may still need a manual
install or vcpkg in `$HOME`.

## 2. Get the source and build

Clone the repository:

```bash
cd $HOME
git clone <this-repo-url> sphinxsys-pd
cd sphinxsys-pd
```

or unpack the source tarball (e.g. `SPHinXsys-PD-capsule-src_YYYY-MM-DD.tar.gz`)
delivered with this fork:

```bash
cd $HOME
tar -xzf SPHinXsys-PD-capsule-src_2026-10-12.tar.gz
cd SPHinXsys
```

Then configure and build (works in either case):

```bash
cmake -G "Unix Makefiles" \
      -D CMAKE_BUILD_TYPE=Release \
      -D CMAKE_C_COMPILER=gcc -D CMAKE_CXX_COMPILER=g++ \
      -D CMAKE_TOOLCHAIN_FILE="$HOME/vcpkg/scripts/buildsystems/vcpkg.cmake" \
      -D CMAKE_C_COMPILER_LAUNCHER=ccache -D CMAKE_CXX_COMPILER_LAUNCHER=ccache \
      -S . -B ./build
cmake --build build/ --target extra_3d_PD_capsule_compression
```

Useful build targets (one executable per example):

| Target | Content |
|---|---|
| `extra_3d_PD_cantilever` | elasto-plastic cantilever (J2 plasticity, validation) |
| `extra_3d_PD_brittle_disc` | brittle disc splitting (elastic + principal-stress fracture) |
| `extra_3d_PD_capsule_compression` | liquid-filled capsule: PD shell + SPH fluid + FSI + plate contact |

Each builds into `tests/3d_examples/<case>/bin/`.

## 3. Run

```bash
cd tests/3d_examples/extra_3d_PD_capsule_compression/bin
./extra_3d_PD_capsule_compression
```

Outputs (body states `.vtp` time series + plate reaction `.dat`) are written
to `./output`. Command-line options (`--help`): `--state_recording`,
`--restart_step`, `--regression`, `--log_level`, etc.

## 4. What is inside, beyond vanilla

| Path | Content |
|---|---|
| `src/shared/peridynamics/` | NOSB-PD module: integration chain (shape matrix → first/second/third/fourth step), Hughes-Winget / J2 / visco-plastic / Drucker-Prager materials, bond-break criteria, Littlewood hourglass control, pair numerical damping, PD-wall contact (with smooth-Coulomb friction), PD–fluid coupling |
| `src/shared/particle_neighborhood/neighborhood.{h,cpp}` | +10 lines: per-bond `bondLive_` state (required by the module) |
| `src/shared/include/sphinxsys.h` | +1 line: `#include "all_peridynamics.h"` |
| `tests/3d_examples/extra_3d_PD_*/` | application examples (one `main` cpp each) |

Note: `tests/3d_examples/extra_3d_PD_ice_jet/` is a third-party case
(SJTU ice-breaking group) kept only for local regression and is **not part of
this distribution** (see `.gitignore`).
