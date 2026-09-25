"""Assemble matched low-energy results into one reproducible comparison table."""
import csv
from pathlib import Path

root=Path(__file__).resolve().parent
def read(name):
    with (root/name).open(newline="") as handle:
        return list(csv.DictReader(handle))

g4={float(read(f"geant4_low_{label}.csv")[0]["energy_MeV"]):
    read(f"geant4_low_{label}.csv")[0] for label in ("1p2","2","5")}
char={float(r["energy_MeV"]):r for r in read("g4_characteristic_low_128x16.csv")}
pstar={float(r["energy_MeV"]):r for r in read("pstar_comparison_128x16.csv")}
native={float(r["energy_MeV"]):r for r in read("radiant_pstar_input_scores.csv")}
radiant={}
for row in read("radiant_low_energy_scores.csv"):
    e=float(row["energy_MeV"])
    if e not in radiant or int(row["lower_energy_groups"])>int(radiant[e]["lower_energy_groups"]):
        radiant[e]=row

out=[]
for e in (1.2,2.0,5.0):
    for m in ("Al","Cu"):
        g=float(g4[e][f"{m}_edep_MeV"])
        r=float(radiant[e][f"{m}_edep_MeV"])
        c=float(char[e][f"{m}_edep_MeV"])
        p=float(pstar[e][f"pstar_{m}_edep_MeV"])
        n=float(native[e][f"{m}_edep_MeV"])
        out.append(dict(energy_MeV=e,material=m,slab_cm=g4[e]["slab_cm"],
                        geant4_MeV=g,geant4_SE_MeV=g4[e][f"{m}_se_MeV"],
                        radiant_G4_table_MeV=r,G4_characteristic_MeV=c,
                        PSTAR_characteristic_MeV=p,radiant_PSTAR_table_MeV=n,
                        radiant_minus_G4_percent=100*(r/g-1),
                        radiant_minus_G4_characteristic_percent=100*(r/c-1),
                        radiant_minus_PSTAR_percent=100*(r/p-1),
                        radiant_PSTAR_table_minus_PSTAR_percent=100*(n/p-1)))

with (root/"low_energy_summary.csv").open("w",newline="") as handle:
    writer=csv.DictWriter(handle,fieldnames=list(out[0]))
    writer.writeheader()
    writer.writerows(out)
for row in out:
    print(f"{row['energy_MeV']:4g} {row['material']:2s} "
          f"G4={row['geant4_MeV']:.8f}±{float(row['geant4_SE_MeV']):.7f} "
          f"Radiant={row['radiant_G4_table_MeV']:.8f} "
          f"G4-char={row['G4_characteristic_MeV']:.8f} "
          f"PSTAR={row['PSTAR_characteristic_MeV']:.8f} "
          f"R-G4={row['radiant_minus_G4_percent']:+.3f}% "
          f"R-PSTAR={row['radiant_minus_PSTAR_percent']:+.3f}%")
