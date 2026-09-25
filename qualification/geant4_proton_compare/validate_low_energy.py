"""Recalculate the exploratory low-energy electronic-CSDA acceptance receipt."""
import csv
import json
from pathlib import Path

root = Path(__file__).resolve().parent

def read(name):
    with (root/name).open(newline="") as handle:
        return list(csv.DictReader(handle))

g4 = {float(read(f"geant4_low_{label}.csv")[0]["energy_MeV"]):
      read(f"geant4_low_{label}.csv")[0] for label in ("1p2","2","5")}
half = {float(read(f"geant4_low_{label}_halfstep.csv")[0]["energy_MeV"]):
        read(f"geant4_low_{label}_halfstep.csv")[0] for label in ("1p2","2")}
pstar = {float(row["energy_MeV"]): row for row in read("pstar_comparison_128x16.csv")}
pstar_coarse = {float(row["energy_MeV"]): row for row in read("pstar_comparison_64x8.csv")}
char = {float(row["energy_MeV"]): row for row in read("g4_characteristic_low_128x16.csv")}
char_coarse = {float(row["energy_MeV"]): row for row in read("g4_characteristic_low_64x8.csv")}
radiant = {}
for row in read("radiant_low_energy_scores.csv"):
    energy = float(row["energy_MeV"])
    if energy not in radiant or int(row["lower_energy_groups"]) > int(radiant[energy]["lower_energy_groups"]):
        radiant[energy] = row
pstar_input_radiant = {float(row["energy_MeV"]):row for row in read("radiant_pstar_input_scores.csv")}

checks = []
def check(name, value, limit, energy, material=""):
    checks.append(dict(name=name, energy_MeV=energy, material=material,
                       observed=value, limit=limit, pass_check=value <= limit))

for E in (1.2,2.0,5.0):
    check("PSTAR min exit above Radiant lower boundary", 1/float(pstar[E]["pstar_min_exit_MeV"]), 1.0, E)
    check("G4 characteristic min exit above lower boundary", 1/float(char[E]["minimum_exit_MeV"]), 1.0, E)
    for m in ("Al","Cu"):
        G = float(g4[E][f"{m}_edep_MeV"])
        SE = float(g4[E][f"{m}_se_MeV"])
        R = float(radiant[E][f"{m}_edep_MeV"])
        P = float(pstar[E][f"pstar_{m}_edep_MeV"])
        C = float(char[E][f"{m}_edep_MeV"])
        check("Radiant versus shared-table characteristic percent", abs(R/C-1)*100, .25, E, m)
        check("Radiant versus independent PSTAR percent", abs(R/P-1)*100, .5, E, m)
        if E in pstar_input_radiant:
            RP = float(pstar_input_radiant[E][f"{m}_edep_MeV"])
            check("PSTAR-fed Radiant versus PSTAR characteristic percent",
                  abs(RP/P-1)*100, .25, E, m)
        check("Radiant versus Geant4 combined percent/SE", abs(R-G)/(0.005*G+3*SE), 1, E, m)
        check("G4 stopping versus PSTAR percent", abs(float(pstar[E][f"stopping_G4_minus_PSTAR_pct_{m}"])), .5, E, m)
        check("PSTAR quadrature change percent", abs(float(pstar_coarse[E][f"pstar_{m}_edep_MeV"])/P-1)*100, .01, E, m)
        check("G4 characteristic quadrature change percent", abs(float(char_coarse[E][f"{m}_edep_MeV"])/C-1)*100, .01, E, m)
        if E in half:
            H = float(half[E][f"{m}_edep_MeV"])
            check("G4 half-step change percent", abs(H/G-1)*100, .05, E, m)

for row in read("pstar_range_consistency.csv"):
    check("PSTAR total-stop/CSDA-range internal consistency percent",
          abs(float(row["difference_percent"])), .2,
          f"{row['low_MeV']}-{row['high_MeV']}", row["material"])

receipt = {"scope":"Al/Cu, 1.2/2/5 MeV, electronic stopping only",
           "criteria":"exploratory, model-relative thresholds; not measurement validation",
           "radiant_selected_grid": {str(E): [int(radiant[E]["voxels_per_material"]),
                                              int(radiant[E]["lower_energy_groups"])] for E in g4},
           "checks":checks,
           "passed":all(item["pass_check"] for item in checks)}
(root/"low_energy_validation_receipt.json").write_text(json.dumps(receipt,indent=2)+"\n")
for item in checks:
    if not item["pass_check"]:
        print("FAIL",item)
print(sum(item["pass_check"] for item in checks),"/",len(checks),"exploratory checks passed")
if not receipt["passed"]:
    raise SystemExit(1)
