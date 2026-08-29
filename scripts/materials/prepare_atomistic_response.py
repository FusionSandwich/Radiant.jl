#!/usr/bin/env python3
"""Validate and hash DFT/MD response-table products before RadiantHTS ingestion.

The script intentionally uses only Python's standard library. It does not run DFT or MD; it checks
that exported tables are complete, energy-conserving where applicable, and bound to an immutable
manifest. Radiant's Julia adapters perform the final schema conversion.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from pathlib import Path
from typing import Iterable

SUBKEV_CHANNELS = (
    "ionization",
    "electronic_excitation",
    "prompt_phonon",
    "athermal_phonon",
    "quasiparticle",
    "defect_storage",
    "optical_emission",
    "electron_escape",
)


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as stream:
        lines = (line for line in stream if line.strip() and not line.lstrip().startswith("#"))
        reader = csv.DictReader(lines)
        if not reader.fieldnames:
            raise ValueError(f"{path} has no CSV header")
        rows = list(reader)
    if not rows:
        raise ValueError(f"{path} has no data rows")
    return rows


def numeric(row: dict[str, str], key: str) -> float:
    if key not in row or row[key] is None or row[key].strip() == "":
        raise ValueError(f"missing numeric field {key}")
    value = float(row[key])
    if not math.isfinite(value):
        raise ValueError(f"non-finite numeric field {key}")
    return value


def require_columns(rows: list[dict[str, str]], columns: Iterable[str]) -> None:
    missing = [column for column in columns if column not in rows[0]]
    if missing:
        raise ValueError(f"missing required columns: {', '.join(missing)}")


def validate_subkev(rows: list[dict[str, str]], tolerance: float) -> dict[str, object]:
    require_columns(rows, ("energy_eV", "temperature_K", *SUBKEV_CHANNELS))
    grid: set[tuple[float, float]] = set()
    for line, row in enumerate(rows, start=2):
        energy = numeric(row, "energy_eV")
        temperature = numeric(row, "temperature_K")
        if energy <= 0 or temperature <= 0:
            raise ValueError(f"row {line}: energy and temperature must be positive")
        key = (energy, temperature)
        if key in grid:
            raise ValueError(f"row {line}: duplicate energy-temperature point {key}")
        grid.add(key)
        fractions = [numeric(row, channel) for channel in SUBKEV_CHANNELS]
        if any(value < 0 or value > 1 for value in fractions):
            raise ValueError(f"row {line}: partition fraction outside [0,1]")
        if abs(sum(fractions) - 1.0) > tolerance:
            raise ValueError(f"row {line}: partition fractions sum to {sum(fractions)}")
    energies = sorted({value[0] for value in grid})
    temperatures = sorted({value[1] for value in grid})
    if len(grid) != len(energies) * len(temperatures):
        raise ValueError("sub-keV table is not a complete energy-temperature grid")
    return {
        "table_type": "subkev_partition",
        "row_count": len(rows),
        "energy_count": len(energies),
        "temperature_count": len(temperatures),
        "energy_closure_tolerance": tolerance,
    }


def validate_green_kubo(rows: list[dict[str, str]]) -> dict[str, object]:
    columns = ("temperature_K", "replica", "kxx_W_m_K", "kyy_W_m_K", "kzz_W_m_K")
    require_columns(rows, columns)
    seen: set[tuple[float, str]] = set()
    counts: dict[float, int] = {}
    for line, row in enumerate(rows, start=2):
        temperature = numeric(row, "temperature_K")
        replica = row["replica"].strip()
        if temperature <= 0 or not replica:
            raise ValueError(f"row {line}: invalid temperature or replica")
        key = (temperature, replica)
        if key in seen:
            raise ValueError(f"row {line}: duplicate Green-Kubo replica {key}")
        seen.add(key)
        counts[temperature] = counts.get(temperature, 0) + 1
        conductivities = [numeric(row, column) for column in columns[2:]]
        if any(value < 0 for value in conductivities):
            raise ValueError(f"row {line}: negative thermal conductivity")
    return {
        "table_type": "green_kubo_conductivity",
        "row_count": len(rows),
        "temperatures_K": sorted(counts),
        "replicas_per_temperature": {str(key): value for key, value in sorted(counts.items())},
    }


def validate_cascade_md(rows: list[dict[str, str]], responses: list[str]) -> dict[str, object]:
    required = ("recoil_energy_eV", "orientation_id", "replica", *responses)
    require_columns(rows, required)
    seen: set[tuple[float, str, str]] = set()
    grid_counts: dict[tuple[float, str], int] = {}
    for line, row in enumerate(rows, start=2):
        energy = numeric(row, "recoil_energy_eV")
        orientation = row["orientation_id"].strip()
        replica = row["replica"].strip()
        if energy <= 0 or not orientation or not replica:
            raise ValueError(f"row {line}: invalid cascade grid coordinates")
        key = (energy, orientation, replica)
        if key in seen:
            raise ValueError(f"row {line}: duplicate cascade replica {key}")
        seen.add(key)
        grid_counts[(energy, orientation)] = grid_counts.get((energy, orientation), 0) + 1
        for response in responses:
            numeric(row, response)
    counts = set(grid_counts.values())
    if len(counts) != 1:
        raise ValueError("cascade grid has unequal replica counts")
    return {
        "table_type": "cascade_md",
        "row_count": len(rows),
        "response_columns": responses,
        "replicas_per_grid_point": next(iter(counts)),
        "grid_point_count": len(grid_counts),
    }


def validate_phonopy(path: Path) -> dict[str, object]:
    required = {"temperature", "heat_capacity", "free_energy", "entropy"}
    current: set[str] = set()
    count = 0
    active = False
    for raw_line in path.read_text(encoding="utf-8").splitlines():
        line = raw_line.strip()
        if line.startswith("thermal_properties:"):
            active = True
            continue
        if not active:
            continue
        if line.startswith("-"):
            if current:
                if not required.issubset(current):
                    raise ValueError(f"incomplete Phonopy entry: {sorted(current)}")
                count += 1
            current = set()
            line = line[1:].strip()
        if ":" in line:
            key = line.split(":", 1)[0].strip()
            if key in required:
                current.add(key)
    if current:
        if not required.issubset(current):
            raise ValueError(f"incomplete Phonopy entry: {sorted(current)}")
        count += 1
    if count == 0:
        raise ValueError("no Phonopy thermal-property entries found")
    return {"table_type": "phonopy_thermal_properties", "entry_count": count}


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("kind", choices=("subkev", "green-kubo", "cascade-md", "phonopy"))
    parser.add_argument("input", type=Path)
    parser.add_argument("output_manifest", type=Path)
    parser.add_argument("--material", required=True)
    parser.add_argument("--material-state-hash", required=True)
    parser.add_argument("--model-id", required=True)
    parser.add_argument("--model-hash", required=True)
    parser.add_argument("--responses", default="surviving_defects,stored_energy_eV")
    parser.add_argument("--closure-tolerance", type=float, default=1.0e-8)
    args = parser.parse_args()

    if not args.input.is_file():
        raise SystemExit(f"input file does not exist: {args.input}")
    if args.kind == "phonopy":
        details = validate_phonopy(args.input)
    else:
        rows = read_rows(args.input)
        if args.kind == "subkev":
            details = validate_subkev(rows, args.closure_tolerance)
        elif args.kind == "green-kubo":
            details = validate_green_kubo(rows)
        else:
            responses = [value.strip() for value in args.responses.split(",") if value.strip()]
            details = validate_cascade_md(rows, responses)

    manifest = {
        "schema": "radiant.hts.atomistic_response_input/v1",
        "classification": "candidate-until-independent-reference-pass",
        "kind": args.kind,
        "material": args.material,
        "material_state_hash": args.material_state_hash,
        "model_id": args.model_id,
        "model_hash": args.model_hash,
        "input_path": str(args.input.resolve()),
        "input_sha256": sha256_file(args.input),
        "details": details,
    }
    args.output_manifest.parent.mkdir(parents=True, exist_ok=True)
    args.output_manifest.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n",
                                    encoding="utf-8")
    print(args.output_manifest)


if __name__ == "__main__":
    main()
