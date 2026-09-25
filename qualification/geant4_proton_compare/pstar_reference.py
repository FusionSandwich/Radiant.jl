"""Independent PSTAR electronic-CSDA slab benchmark for the matched source.

PSTAR's collision stopping data are independent of the Geant4 hIoni table.
Nuclear stopping is excluded because the controlled Geant4 and Radiant cases
transport electronic loss only. This is a characteristic solver, not matRad.
"""
import bisect
import csv
import math
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent
nx = int(sys.argv[1]) if len(sys.argv) > 1 else 64
ne = int(sys.argv[2]) if len(sys.argv) > 2 else 8
max_step_cm = float(sys.argv[3]) if len(sys.argv) > 3 else 1e-5
if nx <= 0 or ne <= 0 or nx % 2 or ne % 2 or max_step_cm <= 0:
    raise ValueError("Quadratures must be positive even counts; path step positive")

def read_csv(path):
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle))

tables = {m: read_csv(ROOT / f"pstar_{m}.csv") for m in ("Al", "Cu")}
g4table = read_csv(ROOT / "geant4_stopping.csv")
grids = {id(rows): [float(row["energy_MeV"]) for row in rows]
         for rows in (*tables.values(), g4table)}

def interpolate(rows, field, e):
    grid = grids[id(rows)]
    if not grid[0] <= e <= grid[-1]:
        raise ValueError(f"Energy {e} escaped table [{grid[0]}, {grid[-1]}]")
    i = min(bisect.bisect_right(grid, e) - 1, len(grid) - 2)
    a, b = float(rows[i][field]), float(rows[i+1][field])
    return a + (b-a)*(e-grid[i])/(grid[i+1]-grid[i])

def electronic_dedx(material, e, density):
    return density*interpolate(tables[material], "electronic_MeV_cm2_g", e)

def propagate(material, e, density, path):
    steps = max(1, math.ceil(path/max_step_cm))
    h = path/steps
    for _ in range(steps):
        f = lambda value: -electronic_dedx(material, value, density)
        k1 = f(e)
        k2 = f(e+h*k1/2)
        k3 = f(e+h*k2/2)
        k4 = f(e+h*k3)
        e += h*(k1+2*k2+2*k3+k4)/6
    return e

def simpson_weights(n):
    return [1 if i in (0,n) else 4 if i % 2 else 2 for i in range(n+1)]

cases = []
for name in ("geant4_low_1p2.csv", "geant4_low_2.csv", "geant4_low_5.csv", "geant4_scores.csv"):
    for row in read_csv(ROOT / name):
        E = float(row["energy_MeV"])
        L = float(row.get("slab_cm", .01))
        rho = {m: float(row[f"{m}_density_g_cm3"]) for m in ("Al", "Cu")}
        result = {"energy_MeV": E, "slab_cm": L}
        result["pstar_Al_edep_MeV"] = result["pstar_Cu_edep_MeV"] = 0.0
        min_exit = math.inf
        for ix, wx in enumerate(simpson_weights(nx)):
            x = L*ix/nx
            for ie, we in enumerate(simpson_weights(ne)):
                e0 = E-.005+.01*ie/ne
                weight = wx*we/(3*nx*3*ne*2)
                for sign in (-1,1):
                    al_path = (x if sign < 0 else L-x)*math.sqrt(3)
                    eal = propagate("Al",e0,rho["Al"],al_path)
                    ecu = propagate("Cu",eal,rho["Cu"],L*math.sqrt(3)) if sign > 0 else eal
                    result["pstar_Al_edep_MeV"] += weight*(e0-eal)
                    result["pstar_Cu_edep_MeV"] += weight*(eal-ecu)
                    min_exit = min(min_exit,ecu)
        result["pstar_min_exit_MeV"] = min_exit
        for m in ("Al", "Cu"):
            pstar = interpolate(tables[m], "electronic_MeV_cm2_g", E)
            geant4 = interpolate(g4table, f"{m}_MeV_cm", E)/rho[m]
            result[f"pstar_{m}_stopping_MeV_cm2_g"] = pstar
            result[f"geant4_{m}_stopping_MeV_cm2_g"] = geant4
            result[f"stopping_G4_minus_PSTAR_pct_{m}"] = 100*(geant4/pstar-1)
            result[f"geant4_{m}_edep_MeV"] = float(row[f"{m}_edep_MeV"])
            result[f"geant4_{m}_se_MeV"] = float(row[f"{m}_se_MeV"])
            result[f"deposition_G4_minus_PSTAR_pct_{m}"] = 100*(result[f"geant4_{m}_edep_MeV"]/result[f"pstar_{m}_edep_MeV"]-1)
        cases.append(result)

out = ROOT / f"pstar_comparison_{nx}x{ne}.csv"
with out.open("w", newline="") as handle:
    writer = csv.DictWriter(handle, fieldnames=list(cases[0]))
    writer.writeheader()
    writer.writerows(cases)
for row in cases:
    print(f"E={row['energy_MeV']:g} MeV L={row['slab_cm']:g} cm "
          f"PSTAR Al={row['pstar_Al_edep_MeV']:.7g} Cu={row['pstar_Cu_edep_MeV']:.7g} "
          f"G4-PSTAR Al={row['deposition_G4_minus_PSTAR_pct_Al']:+.3f}% "
          f"Cu={row['deposition_G4_minus_PSTAR_pct_Cu']:+.3f}% "
          f"min_exit={row['pstar_min_exit_MeV']:.5g}")
