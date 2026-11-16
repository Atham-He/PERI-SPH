# Development snapshot — 2026-11-16

This snapshot contains the isolated liquid-core branch. The solid MCC/SS branch
has not been merged. No completed 0–20% experimental validation is claimed.

## Numerical implementation

- Add a 240 μm PUF-shell / paraffin-oil-core compression driver.
- Derive `IdealJ2` to disable both linear/kinematic hardening and the upstream
  J2 class's hard-coded saturation hardening. Preserve upstream class defaults.
- Use ordinary radial volume layers and exact spherical-sector particle volumes
  inside the actual 2 μm shell; no artificial shell thickening or force offsets.
- Retain the original PD Winfunc influence function when refining the body.
- Use continuous density evolution for the closed WCSPH core instead of copying
  the water-jet free-surface density reset.
- Reuse original virtual-fluid pressure coupling. Schedule both phases at the
  smaller stable time step rather than claiming original multirate subcycling.
- Match the contact-density offset to `NeighborBuilderSurfaceContact`'s
  average-smoothing-length Wendland kernel. Retain the original repulsion law.
- Add force balance, kinetic/recoverable-energy, pressure-volume and resolved
  inner-surface crossing checks. A failed gate does not trigger material tuning.
- Preserve full state checkpoints for reviewed same-model continuation.

## Library changes

1. `large_data_containers.h`: use C++17 standard allocation for neighborhood
   arrays instead of padding every small allocation to a cache line.
2. `cell_linked_list_supplementary.cpp`: remove eager 36-entry allocation in every
   empty grid cell. Search cutoffs, stencils and force formulas are unchanged.
3. `fluid_dynamics_multi_phase.hpp`: correct a CompressibleFluid pointer conversion
   in an upstream template so the library compiles with the current toolchain.

## Reproducible packaging

- Add a standalone CMake target and the existing strict ideal-J2 material audit
  as a CTest test. Full compression is not an automatic test.
- Preserve upstream license texts, file notices, references and source provenance.
- Provide a local-output launcher without developer-machine absolute paths.
- Keep paper-derived inputs separate from external assumptions. The JSON manifest
  documents inputs; a general native configuration parser is not implemented here.
