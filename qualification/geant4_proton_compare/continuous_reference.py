"""Independent 1D characteristic integration using the exported Geant4 dE/dx table.

This checks the deterministic transport limit, not the correctness of the shared
electronic stopping data. Source position and energy are uniform; directions
are the same two S2 ordinates used in the Geant4 and Radiant cases.
"""
import bisect
import csv
import math
import sys
from pathlib import Path

root = Path(__file__).resolve().parent
nx = int(sys.argv[1]) if len(sys.argv) > 1 else 100
max_step_cm = float(sys.argv[2]) if len(sys.argv) > 2 else 2e-5
if nx <= 0 or nx % 2 or max_step_cm <= 0:
    raise ValueError("Source-position intervals must be positive and even; step must be positive")
with (root / "geant4_stopping.csv").open(newline="") as handle:
    rows = list(csv.DictReader(handle))
energies = [float(row["energy_MeV"]) for row in rows]
stopping = {
    material: [float(row[f"{material}_MeV_cm"]) for row in rows]
    for material in ("Al", "Cu")
}


def dedx(material, energy):
    if not energies[0] <= energy <= energies[-1]:
        raise ValueError(f"Energy {energy} MeV escaped the stopping table")
    index = min(bisect.bisect_right(energies, energy) - 1, len(energies) - 2)
    fraction = (energy - energies[index]) / (energies[index + 1] - energies[index])
    values = stopping[material]
    return values[index] + fraction * (values[index + 1] - values[index])


def propagate(material, energy, path_cm):
    steps = max(1, math.ceil(path_cm / max_step_cm))
    h = path_cm / steps
    value = energy
    for _ in range(steps):
        k1 = -dedx(material, value)
        k2 = -dedx(material, value + h*k1/2)
        k3 = -dedx(material, value + h*k2/2)
        k4 = -dedx(material, value + h*k3)
        value += h*(k1 + 2*k2 + 2*k3 + k4)/6
    return value


def weights(n):
    assert n > 0 and n % 2 == 0
    return [1 if i in (0,n) else 4 if i % 2 else 2 for i in range(n + 1)]


L = 0.01
mu = 1/math.sqrt(3)
ne = 10
out = []
for center in (10.0,100.0,499.99):
    al_mean = cu_mean = 0.0
    min_exit = math.inf
    for ix, wx in enumerate(weights(nx)):
        x = L * ix/nx
        for ie, we in enumerate(weights(ne)):
            e0 = center - 0.005 + 0.01*ie/ne
            integration_weight = wx*we/(3*nx*3*ne*2)
            for sign in (-1,1):
                al_path = x/mu if sign < 0 else (L-x)/mu
                e_after_al = propagate("Al",e0,al_path)
                e_after_cu = propagate("Cu",e_after_al,L/mu) if sign > 0 else e_after_al
                al_mean += integration_weight*(e0-e_after_al)
                cu_mean += integration_weight*(e_after_al-e_after_cu)
                min_exit = min(min_exit,e_after_cu)
    out.append({"energy_MeV":center,"Al_edep_MeV":al_mean,"Cu_edep_MeV":cu_mean,"minimum_exit_MeV":min_exit})

with (root / "continuous_reference.csv").open("w",newline="") as handle:
    writer = csv.DictWriter(handle,fieldnames=list(out[0]))
    writer.writeheader()
    writer.writerows(out)
for row in out:
    print(row)
