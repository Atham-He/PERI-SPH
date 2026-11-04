# Contact time diagnosis, anchored to preserved negative evidence

The ea2273184 normal-penalty study stopped after the finest force comparison:
1.061340806 relative curve difference exceeded 0.02. Its finest two-wall
penetration fraction was 0.0225442, also above 0.01. Subsequent time, restart,
quadrature and spatial cases were not dispatched. Separate frozen review
verified all five completed cases' raw arithmetic, inputs and size/SHA evidence;
the normal-contact study remains unqualified.

The accepted output shows intermittent upper contact: the finest sampled
reaction reaches 82.4 micronewtons but is zero at several common times. Final
kinetic energy is 1.995e-12 J versus gel energy 3.697e-12 J. This motivates
distinguishing integration error from transient inertia/contact impacts.
Small energy balance residual alone does not prove reaction convergence.

`peri_gel_contact_time.py` makes two new independent initial states from the
same **existing frozen ELF and physical source archive**, reducing CFL and
maximum-dt cap by 1/2 and 1/4. It never dispatches the anchor or repeats its raw
physical review. Physical speed, stroke, material, mass, geometry, contact
penalty, friction and reference support are checked unchanged. The current
verification-controller commit is recorded separately from the physical-source
commit. Every actual observed process argv, parent/start ticks, boot and ELF is
checked; budgets and stop-on-failure policy remain in force.

Initial actual dt must halve; adaptive median dt must decrease and its actual
ratio is reported. Exact median halving is deliberately **not claimed** for
this adaptive diagnostic. This does not change or pass the original study's
unreached strict temporal gate. Stored-curve relative-force target remains 2%.
A failed first refinement stops the second one. Complete accepted-step force
history was not emitted by the anchor: the diagnosis compares the same 41
stored support points and must not claim a resolved high-frequency force curve
or independent full-step energy.

The independent reviewer reconstructs only the new full-precision particle
frames and verifies source/ELF/dependency/plan/raw SHA evidence. The anchor's
already reviewed receipt and curves are authenticated through its terminal
artifact index; its state and review are immutable. Passing this diagnosis
would provide adaptive-time evidence at the sampled support, not normal-penalty,
rigid-contact, quasi-static, original-size or calibration release.
