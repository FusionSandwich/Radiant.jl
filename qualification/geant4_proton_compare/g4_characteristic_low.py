"""Characteristic limit of the *shared* Geant4 stopping table at low energy.

This separates deterministic transport discretization from the independent
PSTAR stopping-physics comparison. It is not an independent physics model.
"""
import bisect
import csv
import math
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent
nx = int(sys.argv[1]) if len(sys.argv) > 1 else 64
ne = int(sys.argv[2]) if len(sys.argv) > 2 else 8
hmax = float(sys.argv[3]) if len(sys.argv) > 3 else 1e-6
if nx <= 0 or ne <= 0 or nx % 2 or ne % 2 or hmax <= 0:
    raise ValueError("Expected positive even quadrature counts and positive path step")

def read(path):
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle))

table = read(ROOT / "geant4_stopping.csv")
energies = [float(row["energy_MeV"]) for row in table]

def stop(m, e):
    if e < energies[0] or e > energies[-1]:
        raise ValueError(f"Energy {e} MeV escaped Geant4 table")
    i = min(bisect.bisect_right(energies,e)-1,len(energies)-2)
    y0 = float(table[i][f"{m}_MeV_cm"])
    y1 = float(table[i+1][f"{m}_MeV_cm"])
    return y0+(y1-y0)*(e-energies[i])/(energies[i+1]-energies[i])

def path(m,e,d):
    n = max(1,math.ceil(d/hmax))
    h = d/n
    for _ in range(n):
        k1 = -stop(m,e)
        k2 = -stop(m,e+h*k1/2)
        k3 = -stop(m,e+h*k2/2)
        k4 = -stop(m,e+h*k3)
        e += h*(k1+2*k2+2*k3+k4)/6
    return e

def weights(n):
    return [1 if i in (0,n) else 4 if i%2 else 2 for i in range(n+1)]

out = []
for label in ("1p2", "2", "5"):
    case = read(ROOT / f"geant4_low_{label}.csv")[0]
    E,L = float(case["energy_MeV"]),float(case["slab_cm"])
    al = cu = 0.0
    emin = math.inf
    for ix,wx in enumerate(weights(nx)):
        x=L*ix/nx
        for ie,we in enumerate(weights(ne)):
            e0=E-.005+.01*ie/ne
            w=wx*we/(3*nx*3*ne*2)
            for sign in (-1,1):
                eal=path("Al",e0,(x if sign<0 else L-x)*math.sqrt(3))
                ecu=path("Cu",eal,L*math.sqrt(3)) if sign>0 else eal
                al+=w*(e0-eal)
                cu+=w*(eal-ecu)
                emin=min(emin,ecu)
    out.append(dict(energy_MeV=E,slab_cm=L,Al_edep_MeV=al,Cu_edep_MeV=cu,
                    minimum_exit_MeV=emin))
with (ROOT/f"g4_characteristic_low_{nx}x{ne}.csv").open("w",newline="") as handle:
    writer=csv.DictWriter(handle,fieldnames=list(out[0]))
    writer.writeheader()
    writer.writerows(out)
for row in out:
    print(row)
