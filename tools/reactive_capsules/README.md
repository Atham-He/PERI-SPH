# SS30 capsule release and reaction

This isolated development retains S-NSNI/RKPM-DPC for the core and layered
NOSB-PD crack-band damage for the shell. The original running source, binaries,
configuration and pure-MCC evidence are not modified.

The first release representation is **fracture-assisted dissolution/diffusion
through a diffuse crack band**. Exposure means hydraulic access through the
shell, not a resolved free crack face, ejected solid fragments or a liquid core.
The model must never equate a shell damage value with a released mass fraction.

## Evidence and choices

* Target experiment: Yang et al., Acta Geotechnica,
  DOI 10.1007/s11440-026-03253-6, sections 2.1, 2.4 and 3.3.1.
  EC30SS-0.5mm contains 30% SS / 70% MCC. SS modulus (SiO2/Na2O) is 3.
  Shell thickness for the 0.5-mm group is 5 micrometres. The mixture is exposed
  to 1 mol/L CaCl2; the shear test follows 15 minutes of immersion. The crushed
  experimental group is pre-milled; it is not a measured in-shear release movie.
* Crack opening measure: Grassl, "A lattice approach to model flow in cracked
  concrete", arXiv:0809.2758, equation 7. For a damage-only solid, the equivalent
  inelastic opening is `band_width * damage * positive_normal_strain`.
  Our objective finite-strain Biot normal strain and layered series barrier are
  explicit adaptations, not equations claimed verbatim from that paper.
* Transport: finite-volume Fick fluxes with nonnegative conductances and paired
  inventory transfer. The original conforming core cells are reused. A crack
  crossing the complete shell thickness is required for hydraulic access.
* Chemistry: Parkhurst and Appelo (2013), USGS PHREEQC v3,
  https://pubs.usgs.gov/tm/06/a43/; Lothenbach et al. (2019), Cemdata18,
  DOI 10.1016/j.cemconres.2018.04.018. The unmodified Cemdata18 PHREEQC database
  is pinned by SHA256. Ideal CSHQ endmembers are used as a first equilibrium
  assemblage, with finite Ca and released Si/Na inventories.

## Explicit limits

Local reaction equilibrium assumes chemistry is fast relative to the transport
step; it does not supply reaction times or a nucleation law. The initial CSHQ
selection omits alkali uptake by C-S-H. Cemdata18 Appendix C states an approximate
ionic-strength range up to 1 molal for its extended Debye-Hueckel correction;
1 M CaCl2 exceeds that range. Reaction outputs must carry this extrapolation
flag and cannot yet be treated as quantitative predictions for the target bath.
No precipitation amount is converted to sand-bond strength without a separate
constitutive relationship and evidence. Dissolution-to-mechanics feedback and
solid fragmentation are not silently supplied.

The mechanical interface remains the inherited bonded trace constraint. This
first model describes hydraulic exposure through a diffuse shell crack band;
it does not predict free shell fragments, mechanical debonding or solid-core
extrusion. Pores are assumed to imbibe rapidly once connected to external water.
Finite-volume diffusion uses the normal projection of the Piola-pulled-back
isotropic diffusivity, with deformed boundary areas. Its two-point face flux is
an explicit approximation for sheared/nonorthogonal cells. Pore volume is mapped
from the assumed DPC relative-density state; that mapping is not measured for
the target mixture. The external dilute-sink approximation is checked against
the bath equilibrium concentration. None of these assumptions is a fitted
force correction.

## Acceptance nodes

1. An intact shell blocks release; a through-thickness tensile crack permits
   access. Closure blocks the crack contribution. Access is not a fitted switch
   at a chosen capsule strain.
2. One linked material/barrier/transport/reaction gate verifies nonnegative
   inventories, conserved Si/Na/Ca, and zero C-S-H for zero SS or zero calcium.
3. One EC30SS-0.5mm capsule loading history exports actual damage, exposure and
   release histories. No parameter matrix and no claim of experimental release
   validation without release measurements.
