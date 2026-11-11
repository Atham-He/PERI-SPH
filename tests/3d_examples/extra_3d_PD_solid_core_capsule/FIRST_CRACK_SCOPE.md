# First-crack diagnostic and independent shell material selection

`shell-material=elastic` remains the default. Optional `shell-material=j2`
requires an explicit positive `shell-yield-stress`; the default hardening modulus
is zero, with optional `shell-hardening-modulus` for linear isotropic hardening.
The quadratic plastic path uses the native Hughes–Winget J2 return mapping.
Selecting a material never selects or identifies a fracture threshold.

The supported diagnostic is `fracture-model=principal_stress` together with
`fracture-scope=first_crack` and an explicit `fracture-stress` in Pa. It detects
the largest principal value of the symmetric average stress of each live bond's
two endpoints, commits both directed states irreversibly, and outputs `Damage`,
`BondNumber`, the broken-bond list, and `first_crack.json`. There is no hidden
50% stretch fallback. Material-point maximum stress alone is not the criterion.

**The diagnostic stops immediately after the first sampled bond failure.** It
does not evaluate an intact-only gradient or surface map on the broken topology.
Its recorded reaction is before bond failure, not a post-fracture equilibrium.
Full fracture growth is deliberately rejected at input validation until live-bond
gradient reconstruction, crack-side surface transfer, fragment contact, and an
appropriate fracture-energy treatment are qualified. Do not remove this guard
to obtain an apparent post-peak curve.

A tensile strength from another EC formulation is a prior, not this shell's
measured strength, and it must not be reinterpreted as J2 yield stress. A first
local bond-failure event and the experimentally observed macroscopic rupture
peak are different events. Existing intact-field screens cannot be continued
as damaged histories once the criterion has triggered.
