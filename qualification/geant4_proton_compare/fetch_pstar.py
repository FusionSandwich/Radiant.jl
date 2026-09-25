"""Fetch and preserve independent NIST PSTAR Al/Cu proton tables.

Only the official NIST text interface is queried. No Python dependencies.
"""
import csv
import hashlib
import html
import re
import urllib.parse
import urllib.request
from pathlib import Path

ROOT = Path(__file__).resolve().parent
URL = "https://physics.nist.gov/cgi-bin/Star/ap_table-t.pl"
GRID = sorted(set([*(round(1.0 + i * .05, 3) for i in range(1, 21)),
                   *(round(2.1 + i * .1, 3) for i in range(30)),
                   *(round(5.5 + i * .5, 3) for i in range(30)),
                   20.0, 100.0, 500.0]))
ROW = re.compile(r"^\s*([\d.]+E[+-]\d+)\s+" +
                 r"\s+".join([r"([\d.]+E[+-]\d+)"] * 5 + [r"([\d.]+)"]) + r"\s*$")
fields = ["energy_MeV", "electronic_MeV_cm2_g", "nuclear_MeV_cm2_g",
          "total_MeV_cm2_g", "csda_range_g_cm2", "projected_range_g_cm2",
          "detour_factor"]
for material, number in (("Al", "013"), ("Cu", "029")):
    payload = urllib.parse.urlencode({"prog": "PSTAR", "matno": number,
                                      "ShowDefault": "on",
                                      "Energies": "\n".join(map(str, GRID))}).encode()
    request = urllib.request.Request(URL, data=payload,
                                     headers={"User-Agent": "Radiant-proton-qualification/1.0"})
    with urllib.request.urlopen(request, timeout=30) as response:
        raw = response.read()
    if len(raw) > 200_000:
        raise RuntimeError("PSTAR response exceeded planned acquisition size")
    text = html.unescape(re.sub(r"<[^>]*>", "\n", raw.decode("utf-8", errors="replace")))
    rows = [match.groups() for line in text.splitlines()
            if (match := ROW.fullmatch(line))]
    energies = [float(row[0]) for row in rows]
    if len(rows) < 130 or any(b <= a for a, b in zip(energies, energies[1:])):
        bad = [(a, b) for a, b in zip(energies, energies[1:]) if b <= a]
        raise RuntimeError(f"Invalid {material} PSTAR table: {len(rows)} rows; non-increasing {bad}")
    missing = [target for target in (1.0, 1.2, 2.0, 5.0, 10.0, 100.0, 500.0)
               if not any(abs(e - target) < 1e-9 for e in energies)]
    if missing:
        raise RuntimeError(f"{material} PSTAR table is missing benchmark energies {missing}; "
                           f"near 1 MeV: {[e for e in energies if .9 <= e <= 2.1]}")
    (ROOT / f"pstar_{material}_response.html").write_bytes(raw)
    with (ROOT / f"pstar_{material}.csv").open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(fields)
        writer.writerows(rows)
    print(material, len(rows), "rows", len(raw), "bytes", hashlib.sha256(raw).hexdigest())
