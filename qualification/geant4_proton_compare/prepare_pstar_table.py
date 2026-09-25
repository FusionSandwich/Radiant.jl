"""Convert official PSTAR mass electronic stopping into Radiant's MeV/cm input.

Use the matched Geant4 elemental densities. Preserve the original NIST tables
separately; this derived candidate does not include nuclear stopping.
"""
import csv
from pathlib import Path

root = Path(__file__).resolve().parent
data = {}
for material in ("Al", "Cu"):
    with (root / f"pstar_{material}.csv").open(newline="") as handle:
        rows = list(csv.DictReader(handle))
    data[material] = {float(row["energy_MeV"]): float(row["electronic_MeV_cm2_g"])
                      for row in rows if 1 <= float(row["energy_MeV"]) <= 500}
grid = sorted(data["Al"].keys() & data["Cu"].keys())
if len(grid) < 50 or grid[0] != 1 or grid[-1] != 500:
    raise RuntimeError("PSTAR Al/Cu grids do not cover Radiant's 1–500 MeV contract")
with (root / "pstar_electronic_stopping.csv").open("w", newline="") as handle:
    writer = csv.writer(handle)
    writer.writerow(("energy_MeV", "Al_MeV_cm", "Cu_MeV_cm"))
    for e in grid:
        writer.writerow((f"{e:.8g}", f"{data['Al'][e]*2.699:.12g}",
                         f"{data['Cu'][e]*8.96:.12g}"))
print(len(grid), "PSTAR points from", grid[0], "to", grid[-1], "MeV")
