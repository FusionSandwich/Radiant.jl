"""Check parsed PSTAR CSDA ranges against its total stopping column.

NIST defines these ranges using electronic plus nuclear stopping; this
consistency check is separate from the electronic-only slab benchmark.
"""
import bisect
import csv
from pathlib import Path

root = Path(__file__).resolve().parent
out = []
for material in ("Al","Cu"):
    with (root/f"pstar_{material}.csv").open(newline="") as handle:
        rows = list(csv.DictReader(handle))
    grid = [float(row["energy_MeV"]) for row in rows]
    stop = [float(row["total_MeV_cm2_g"]) for row in rows]
    ranges = {float(row["energy_MeV"]):float(row["csda_range_g_cm2"]) for row in rows}
    def reciprocal(e):
        i = min(bisect.bisect_right(grid,e)-1,len(grid)-2)
        s = stop[i]+(stop[i+1]-stop[i])*(e-grid[i])/(grid[i+1]-grid[i])
        return 1/s
    for low,high in ((1,2),(2,5),(5,10)):
        n = int((high-low)/.001)
        h = (high-low)/n
        integral = h*(reciprocal(low)+reciprocal(high)+
                      2*sum(reciprocal(low+i*h) for i in range(1,n)))/2
        listed = ranges[high]-ranges[low]
        difference = 100*(integral/listed-1)
        out.append(dict(material=material,low_MeV=low,high_MeV=high,
                        integrated_range_increment_g_cm2=integral,
                        pstar_range_increment_g_cm2=listed,difference_percent=difference))

with (root/"pstar_range_consistency.csv").open("w",newline="") as handle:
    writer=csv.DictWriter(handle,fieldnames=list(out[0]))
    writer.writeheader()
    writer.writerows(out)
for row in out:
    print(row)
if any(abs(row["difference_percent"])>.5 for row in out):
    raise SystemExit("PSTAR range/stopping consistency check failed")
