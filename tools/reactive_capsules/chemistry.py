"""Finite-inventory SS/CaCl2 equilibrium using unmodified PHREEQC/Cemdata18.

The input SS is Na2O.(SiO2)_modulus on a dry, anhydrous basis. Its hydration
state is an explicit assumption. No hand-written precipitation rate or fixed
Ca/Si ratio is imposed. This batch operator does NOT predict reaction kinetics
or the mechanical strength of a sand-grain bond.
"""
from pathlib import Path
import csv
import hashlib
import json
import math
import subprocess


class CemdataReaction:
    COMPONENTS = ("CSHQ-TobH", "CSHQ-TobD", "CSHQ-JenH", "CSHQ-JenD")

    def __init__(self, binary, database, workdir, water_kg, calcium_mol,
                 temperature_C=25.0, modulus=3.0):
        self.binary, self.database = Path(binary), Path(database)
        self.workdir = Path(workdir)
        self.workdir.mkdir(parents=True, exist_ok=True)
        self.water_kg, self.calcium_mol = water_kg, calcium_mol
        self.temperature_C, self.modulus = temperature_C, modulus
        if not self.binary.is_file() or not self.database.is_file():
            raise ValueError("PHREEQC executable and pinned Cemdata database are required")
        if min(water_kg, modulus) <= 0 or calcium_mol < 0:
            raise ValueError("Invalid finite bath")
        self.calls = 0

    def equilibrate(self, released_ss_formula_mol, label):
        if not math.isfinite(released_ss_formula_mol) or released_ss_formula_mol < 0:
            raise ValueError("Negative/non-finite released SS inventory")
        self.calls += 1
        stem = self.workdir / label
        if stem.with_suffix(".pqi").exists():
            raise FileExistsError("Reaction evidence must not be overwritten: " + str(stem))
        ca = self.calcium_mol / self.water_kg
        ss = released_ss_formula_mol / self.water_kg
        comp = "\n".join(f"    -comp {name} 0" for name in self.COMPONENTS)
        # Normalize to one kg of water for solver tolerances, then scale moles
        # back. Ca and chloride inventories are finite; neither is replenished.
        source = f"""TITLE Finite SS release into calcium chloride; local equilibrium hypothesis
SOLUTION 1
    temp {self.temperature_C:.17g}
    units mol/kgw
    pH 7 charge
    Ca {ca:.17g}
    Cl {2*ca:.17g}
    -water 1
SOLID_SOLUTIONS 1
  CSHQ
{comp}
REACTION 1
    Na2O 1
    SiO2 {self.modulus:.17g}
    {ss:.17g} moles
SELECTED_OUTPUT 1
    -file {stem.with_suffix('.selected')}
    -reset false
    -high_precision true
    -user_punch true
USER_PUNCH 1
    -headings pH ionic_strength aq_Si aq_Ca aq_Na aq_Cl csh_Si csh_Ca csh_Na csh_Cl sys_Si sys_Ca sys_Na sys_Cl
    -start
10 PUNCH -LA("H+"), MU, TOTMOLE("Si"), TOTMOLE("Ca"), TOTMOLE("Na"), TOTMOLE("Cl")
20 PUNCH SUM_S_S("CSHQ", "Si"), SUM_S_S("CSHQ", "Ca"), SUM_S_S("CSHQ", "Na"), SUM_S_S("CSHQ", "Cl")
30 PUNCH SYS("Si"), SYS("Ca"), SYS("Na"), SYS("Cl")
    -end
END
"""
        stem.with_suffix(".pqi").write_text(source)
        r = subprocess.run([str(self.binary), str(stem.with_suffix('.pqi')),
                            str(stem.with_suffix('.out')), str(self.database)],
                           cwd=self.workdir, capture_output=True, text=True, timeout=120)
        stem.with_suffix(".log").write_text(r.stdout + r.stderr)
        if r.returncode:
            raise RuntimeError(f"PHREEQC failed ({r.returncode}); see {stem}.out")
        lines = [x.split() for x in stem.with_suffix('.selected').read_text().splitlines() if x.strip()]
        if len(lines) < 2:
            raise RuntimeError("PHREEQC produced no selected output")
        raw = dict(zip(lines[0], map(float, lines[-1])))
        result = {k: (v if k in ('pH', 'ionic_strength') else v*self.water_kg)
                  for k, v in raw.items()}
        expected = dict(Si=self.modulus*released_ss_formula_mol,
                        Ca=self.calcium_mol, Na=2*released_ss_formula_mol,
                        Cl=2*self.calcium_mol)
        residual = {}
        for element, total in expected.items():
            recovered = result['aq_'+element] + result['csh_'+element]
            # Normalize with actual finite inventory, with a tiny absolute floor.
            residual[element] = abs(recovered-total)/max(total, 1e-20)
            if abs(recovered-total) > max(1e-18, 1e-7*total):
                raise RuntimeError(f"Reaction violates {element} mass balance: {recovered} vs {total}")
        result.update(released_ss_formula_mol=released_ss_formula_mol,
                      relative_element_balance_residual=residual,
                      equilibrium_not_kinetics=True,
                      csh_model='Cemdata18 four-endmember ideal CSHQ',
                      high_ionic_strength_extrapolation=raw['ionic_strength'] > 1,
                      database_sha256=hashlib.sha256(self.database.read_bytes()).hexdigest())
        stem.with_suffix('.json').write_text(json.dumps(result, indent=2)+'\n')
        return result
