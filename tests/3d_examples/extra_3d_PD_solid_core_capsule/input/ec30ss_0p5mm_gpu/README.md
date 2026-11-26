# GPU single-capsule compression

The solver source for the production run is commit
`770a1746857c6eed151670fa07efa309afba24c9`. This case uses a 30% SS + 70% MCC
solid core with S-NSNI/RKPM and DPC, an EC shell with layered NOSB-PD, and a
bonded solid interface. The model outer diameter is 0.510 mm and the shell
thickness is 5 μm. The numerical end point is 40.5% platen closure relative
to the initial diameter, in 81 loading stages.

`material.csv`, `reference_cells_h25.bin`, and `crack_bank.txt` are the exact
production material table, reference geometry, and candidate crack faces.
The comments and material-table path in `mechanics.cfg` are portable.
Its numerical and model settings match the production input. Verify the
included files with `sha256sum -c SHA256SUMS` from this directory.

## Build on Linux or WSL

Install CMake, a CUDA-compatible C++17 compiler, the CUDA toolkit, Simbody,
Eigen3, oneTBB, and Boost. Run from the repository root.
The example below uses CUDA 12.8 and architecture 89 for the RTX 4070 SUPER.
Set the compiler path and architecture for the available GPU.

```sh
example=tests/3d_examples/extra_3d_PD_solid_core_capsule
build="$PWD/build/capsule-gpu"
case_dir="$PWD/$example/input/ec30ss_0p5mm_gpu"
cmake -S "$example" -B "$build" \
  -DCMAKE_BUILD_TYPE=Release -DPERI_SPH_ENABLE_CUDA=ON \
  -DCMAKE_CUDA_COMPILER=/usr/local/cuda-12.8/bin/nvcc \
  -DCMAKE_CUDA_ARCHITECTURES=89
cmake --build "$build" --target pd_snsni_capsule pd_physical_crack_benchmark --parallel 2
```

## Check the production configuration

Use a new output directory for each invocation. The production solver uses an
implicit ramp with explicit fallback and explicit holds. The original plate
reaction settling window is 10 ns. Shell fracture deletes bonds and uses
contact on closure. The Gc input is retained in the configuration, with
cohesive softening disabled by the contact-only run option.

```sh
checks="$(mktemp -d)"
"$build/pd_physical_crack_benchmark" "$checks/benchmark" \
  --execution=gpu --fracture-contact-only
"$build/pd_snsni_capsule" \
  --config="$case_dir/mechanics.cfg" --dpc-table="$case_dir/material.csv" \
  --snsni-cells="$case_dir/reference_cells_h25.bin" \
  --fracture-facets="$case_dir/crack_bank.txt" --output="$checks/preflight" \
  --execution=gpu --fracture-contact-only --affine-crack-contact \
  --affine-inner-interface --local-platen-step-bound --tighter-shell-step-bound \
  --reaction-settling-window-s=1e-08 --hybrid-loading --preflight-only
```

The [production launcher](../../../../../tools/peri_run_physical_capsule.py)
takes the matching benchmark, preflight, and hybrid qualification receipts.
For the launcher, prepare a working copy of the configuration with the
absolute path to `material.csv`. The launcher checks receipts against that
configuration, the three inputs, and the executable. A continuation requires
the original run manifest, checkpoint, and matching historical input bytes.
Checkpoints, numerical results, and machine-specific receipts remain with
the original run.
