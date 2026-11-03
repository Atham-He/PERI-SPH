# Single-particle compression: code and experimental assessment

Source reviewed: merged `capsule` commit `364f12f1a`, based on official master
`1612330e`. Experimental source: Yang, Chen and Lourenco,
Powder Technology 469 (2026) 121906, DOI `10.1016/j.powtec.2025.121906`,
and the supplied `Nano-micro mechanics.xlsx`. All data imports are read-only.

## Findings that change the experimental mapping

1. The tested materials have an alginate matrix with distributed oil chambers,
   rather than a demonstrated single uniform liquid core and thin shell.
   Sections 2.1/3.1 and figure 7 establish this distinction. Zero TO means no
   tung oil, while hydrated gel still contains water. Neither the current hollow
   dry-shell control nor a pure-tung-oil-core capsule is an experimentally matched
   `0TO-3SA-H` model. Keep the idealized core-shell model as a numerical control.
2. Section 3.1 reports group particle-size distributions: the 2.34 mm value for
   `0TO-3SA-H` is a D50, not the measured diameter of its exported force curve.
   The paper does not establish the current model's 20 micrometre uniform shell.
3. The microscale test uses parallel plates (section 2.4). Its speed and a
   per-point time axis are not specified in the supplied workbook/paper.
   The 5/10/5-second protocol belongs to Berkovich nanoindentation only.
4. Micro sheets use micrometres/millinewtons; Nano sheets use
   nanometres/micronewtons. Nano loading/unloading order must be retained.
   The original Micro-H duplicate pairs are also equal in the source workbook;
   they are not an artefact of our CSV export and are not independent replicates.
5. Local nanoindentation moduli and effective whole-particle Hertz moduli probe
   different scales. Neither can be substituted blindly as a thin-shell modulus.
   The initial Hertz interpretation is restricted to approximately d/D < 0.1.
6. Section 3.3 reports large deformation and fluid exudation, including cases
   without macroscopic fracture. Force changes/leakage are not sufficient evidence
   for calibrating a brittle bond-breaking threshold.

The machine-readable source facts and unknowns are in
`peri/experiments/powtec_2025_121906.json`. It is not a runnable case config.

## Read-only data preparation

```sh
python tools/peri_experiment.py --workbook /path/to/Nano-micro\ mechanics.xlsx \
  --output /path/to/a/new/experiment-export
```

This optional command requires `openpyxl` from `peri/requirements-analysis.txt`.
The supplied workbook produced 36 curves and 230812 numeric pairs. Every exported
pair, original row index and SI conversion was independently checked against the
workbook. All 18 prior Micro CSV hashes remain unchanged. The Micro sheets have
no decreasing-displacement transitions in source row order; Nano sheets do. Neither
observation supplies missing physical timestamps or a microscale loading rate.

## Existing code and remaining development

| Area | Existing implementation | Assessment and next action |
|---|---|---|
| Liquid-PD interface | `pd_fsi` and the capsule's pressure/density/viscous operators | R8 corrected gradients remain in an independent prescribed-wall prototype restricted to <=1000 fluid particles and one wall. Porting requires explicit paired operators, changing topology, phase cursors and independent conservation checks. Do not remove the prototype guard to claim scalability. |
| Capsule geometry | Lattice shell from the difference of two spheres | Add an explicit effective-gel solid-sphere numerical baseline for zero TO, then assess a homogenized/multi-chamber model for oil-bearing particles. Geometry and loading inputs must be justified before experimental fitting. |
| Material | Objective elastic update and Perzyna-type J2 viscoplasticity | These are retained implementations, not calibrated gel material laws. Verify reversible relaxation and permanent deformation separately, and assess whether a poro/visco-elasto-plastic model is needed for exudation. Do not identify all rate parameters from one monotonic curve. |
| Force observation | `measure` evaluates contact and fluid force fields at output points | The observer writes acceleration/force-related arrays. A dedicated read-only observer and explicit pressure/geometry/velocity phase metadata are required before claiming common-state instantaneous reaction. Do not add a terminal dynamics kick. |
| Conservation | Mass inventory, kick-work identities and mechanical residual outputs | Existing diagnostics are useful. Complete independent interface work/impulse checks; shell storage is still a proxy and fracture energy is unresolved. Sparse R8 samples do not establish full-step energy closure. |
| PD spatial accuracy | Shape matrix, stabilization, fixed material neighborhoods | Preserve objectivity and pair symmetry. Separately refine spacing, horizon and shell/solid resolution, including boundary support and damage. The current nominal three shell layers are not a convergence result. |
| Performance | Sparse cached small prototype; reusable CPU PD code | Develop scalable operators/preconditioning before original-size corrected coupling. Avoid raw-output/archive growth. GPU migration follows a validated CPU baseline. |
| Restart | Existing PD full-state registration and atomic generations | Keep old checkpoint evidence. Extend the schema for new pair IDs, correction weights, active masks, reference/cache generation, material memory and committed output cursors; verify rollback and split-continuous execution again. |

## Optimizations completed in this pass

- Added `tools/peri_experiment.py`: a single worksheet traversal extracts all
  labelled pairs, preserves row order/repeated displacements and both duplicate
  labels, converts known units to SI, records hashes and keeps time/geometry
  unknown. Invalid pairs/units/labels fail explicitly. The old workbook and CSVs
  are not overwritten. This is data preparation, not constitutive calibration.
- The launcher now includes the standalone `peri/` build/profile inputs in its
  source receipt. It prunes generated output/checkpoint/build trees before
  descending into them, retaining authored input/reload geometry. Stale-source
  preflight fails before creating a misleading run directory.
- Analytic rigid-box contact no longer refreshes the three particle-contact
  configurations that its force path does not consume. The particle-contact
  branch keeps the original refreshes. Bounds checks and fluid/shell coupling
  refreshes remain intact. No constitutive/EOS/mass/time-step threshold changes
  are introduced by this optimization.

## Sequence and acceptance

1. Keep the existing R8 24/36us evidence attached to its frozen old executable.
   The 48us attempt remains resource-stopped and is not a physical pass.
2. Bridge the new version with small numerical controls and complete corrected
   closed-cavity operators. Confirm the effective-gel geometry baseline and
   interface assumptions before treating a full-size curve as experimental.
3. Establish time, EOS, spacing/horizon, reaction, mass and work accuracy at the
   intended dimensions. Preserve all existing numerical/resource limits.
4. Build and verify the relevant finite-strain material law, then calibrate only
   identifiable parameters. Retain `0TO-3SA-H` as the first experimental target.
5. Add calibrated damage/multi-chamber transport only when required by the
   observed failure/release mechanism; account for released water/oil by flux.

This pass does not run the full capsule, resume the original queue, fit material
parameters, change a frozen R8 artifact, or declare the new solver converged.
Coupled numerical equivalence of the contact-refresh optimization still needs
a versioned solver regression after the existing disk-reserve requirement is met.
