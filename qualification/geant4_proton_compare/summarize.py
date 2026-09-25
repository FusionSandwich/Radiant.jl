"""Compare per-source deposited MeV from matched Geant4 and Radiant slab cases."""
import csv
import math
from pathlib import Path

root = Path(__file__).resolve().parent
with (root / "geant4_scores.csv").open(newline="") as handle:
    geant4 = {float(row["energy_MeV"]): row for row in csv.DictReader(handle)}
with (root / "continuous_reference.csv").open(newline="") as handle:
    continuous = {float(row["energy_MeV"]): row for row in csv.DictReader(handle)}

rows = []
with (root / "radiant_scores_mass_corrected.csv").open(newline="") as handle:
    for radiant in csv.DictReader(handle):
        energy = float(radiant["energy_MeV"])
        g4 = geant4[energy]
        for material in ("Al", "Cu"):
            value = float(radiant[f"{material}_edep_MeV"])
            reference = float(g4[f"{material}_edep_MeV"])
            characteristic = float(continuous[energy][f"{material}_edep_MeV"])
            se = float(g4[f"{material}_se_MeV"])
            assert all(math.isfinite(x) for x in (value, reference, se)) and reference > 0 and se > 0
            rows.append({
                "energy_MeV": energy,
                "material": material,
                "voxels_per_material": int(radiant["voxels_per_material"]),
                "lower_energy_groups": int(radiant["lower_energy_groups"]),
                "geant4_MeV_per_source": reference,
                "geant4_standard_error_MeV": se,
                "radiant_MeV_per_source": value,
                "characteristic_MeV_per_source": characteristic,
                "radiant_minus_geant4_pct": 100 * (value / reference - 1),
                "radiant_minus_characteristic_pct": 100 * (value / characteristic - 1),
                "geant4_minus_characteristic_pct": 100 * (reference / characteristic - 1),
                "difference_over_geant4_SE": (value - reference) / se,
            })

out = root / "deposition_comparison.csv"
with out.open("w", newline="") as handle:
    writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
    writer.writeheader()
    writer.writerows(rows)
for row in rows:
    print(f"{row['energy_MeV']:7.2f} MeV {row['material']} "
          f"{row['voxels_per_material']:2d}x{row['lower_energy_groups']:3d}: "
          f"{row['radiant_minus_geant4_pct']:+8.3f}% "
          f"({row['difference_over_geant4_SE']:+7.2f} G4 SE)")
