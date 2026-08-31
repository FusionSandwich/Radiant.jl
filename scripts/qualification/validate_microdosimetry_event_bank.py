#!/usr/bin/env python3
"""Validate a Julia-written weighted microdosimetry bank through Python/h5py."""

from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path
from typing import Any

import h5py
import numpy as np

SCHEMA = "radiant.hts.weighted_microdosimetry_event_bank/v2"
PARTITION_CHANNELS = (
    "ionization_eV",
    "electronic_excitation_eV",
    "prompt_lattice_heat_eV",
    "nuclear_recoil_handoff_eV",
    "defect_stored_eV",
    "optical_emission_eV",
    "escaping_particle_eV",
    "cutoff_handoff_eV",
    "unresolved_eV",
)


def _text(value: Any) -> str:
    if isinstance(value, bytes):
        return value.decode("utf-8")
    return str(value)


def _text_array(dataset: h5py.Dataset) -> list[str]:
    return [_text(value) for value in dataset[...].tolist()]


def _require(condition: bool, message: str) -> None:
    if not condition:
        raise ValueError(message)


def validate(path: Path) -> dict[str, Any]:
    with h5py.File(path, "r") as handle:
        _require(_text(handle.attrs.get("schema_id", "")) == SCHEMA, "schema mismatch")
        _require(
            bool(handle.attrs.get("representative_not_analog", False)),
            "bank must be labelled representative, not analog",
        )
        source_hash = _text(handle.attrs.get("source_hash", ""))
        kernel_hash = _text(handle.attrs.get("kernel_hash", ""))
        _require(
            len(source_hash) == 64
            and len(kernel_hash) == 64
            and all(character in "0123456789abcdefABCDEF" for character in source_hash)
            and all(character in "0123456789abcdefABCDEF" for character in kernel_hash),
            "lineage hashes must be SHA-256 digests",
        )

        events = handle["events"]
        secondaries = handle["secondaries"]
        count = int(events.attrs["event_count"])
        _require(count > 0, "event bank cannot be empty")

        required_event_vectors = ("position_cm", "direction")
        for name in required_event_vectors:
            _require(
                events[name].shape == (count, 3),
                f"events/{name} must be event,component",
            )
            _require(
                np.isfinite(events[name][...]).all(), f"events/{name} must be finite"
            )
        directions = events["direction"][...]
        _require(
            np.allclose(np.linalg.norm(directions, axis=1), 1.0, rtol=0.0, atol=1e-12),
            "parent directions must be normalized",
        )

        required_event_scalars = (
            "sample_id",
            "prototype_id",
            "particle_tag",
            "process_id",
            "material_tag",
            "layer_id",
            "correlation_id",
            "time_s",
            "direction_model",
            "incident_energy_eV",
            "deposited_energy_eV",
            "statistical_weight_events",
        )
        for name in required_event_scalars:
            _require(events[name].shape == (count,), f"events/{name} length mismatch")
        weights = events["statistical_weight_events"][...]
        _require(
            np.isfinite(weights).all() and np.all(weights > 0.0),
            "statistical weights must be finite and positive",
        )
        incident = events["incident_energy_eV"][...]
        deposited = events["deposited_energy_eV"][...]
        _require(
            np.isfinite(incident).all() and np.all(incident > 0.0),
            "incident energies must be finite and positive",
        )
        _require(
            np.isfinite(deposited).all()
            and np.all(deposited >= 0.0)
            and np.all(deposited <= incident),
            "deposited energies must be physical",
        )

        partition = events["partition"]
        available = partition["available_energy_eV"][...]
        channel_sum = sum(
            (partition[name][...] for name in PARTITION_CHANNELS),
            np.zeros(count, dtype=np.float64),
        )
        _require(
            np.allclose(available, deposited, rtol=1e-12, atol=1e-12),
            "partition available energy must equal deposited energy",
        )
        _require(
            np.allclose(channel_sum, available, rtol=1e-12, atol=1e-12),
            "energy partition does not close",
        )

        parent_index = secondaries["parent_event_index_1based"][...].astype(np.int64)
        secondary_count = len(parent_index)
        _require(
            np.all((parent_index >= 1) & (parent_index <= count)),
            "secondary parent index is out of bounds",
        )
        required_secondary_scalars = (
            "particle_tag",
            "kinetic_energy_eV",
            "direction_model",
            "correlation_id",
        )
        for name in required_secondary_scalars:
            _require(
                secondaries[name].shape == (secondary_count,),
                f"secondaries/{name} length mismatch",
            )
        _require(
            secondaries["direction"].shape == (secondary_count, 3),
            "secondaries/direction must be secondary,component",
        )
        secondary_directions = secondaries["direction"][...]
        _require(
            np.isfinite(secondary_directions).all(),
            "secondary directions must be finite",
        )
        if secondary_count:
            _require(
                np.allclose(
                    np.linalg.norm(secondary_directions, axis=1),
                    1.0,
                    rtol=0.0,
                    atol=1e-12,
                ),
                "secondary directions must be normalized",
            )
        secondary_energy = secondaries["kinetic_energy_eV"][...]
        _require(
            np.isfinite(secondary_energy).all() and np.all(secondary_energy >= 0.0),
            "secondary energies must be finite and nonnegative",
        )

        parent_correlation = _text_array(events["correlation_id"])
        secondary_correlation = _text_array(secondaries["correlation_id"])
        secondary_models = _text_array(secondaries["direction_model"])
        for index, model in enumerate(secondary_models):
            if model == "parent_correlated":
                parent = int(parent_index[index]) - 1
                _require(
                    secondary_correlation[index] == parent_correlation[parent],
                    "parent-correlated lineage identifier mismatch",
                )
                _require(
                    np.allclose(
                        secondary_directions[index],
                        directions[parent],
                        rtol=0.0,
                        atol=1e-12,
                    ),
                    "parent-correlated actual vector mismatch",
                )

    digest = hashlib.sha256(path.read_bytes()).hexdigest()
    return {
        "schema": "radiant.hts.microdosimetry_cross_language_replay/v1",
        "passed": True,
        "artifact": str(path.resolve()),
        "artifact_sha256": digest,
        "artifact_bytes": path.stat().st_size,
        "event_count": count,
        "secondary_count": secondary_count,
        "statistical_weight_sum": math.fsum(float(value) for value in weights),
        "source_hash": source_hash,
        "kernel_hash": kernel_hash,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("artifact", type=Path)
    parser.add_argument("--receipt", type=Path)
    args = parser.parse_args()
    result = validate(args.artifact)
    if args.receipt is not None:
        args.receipt.parent.mkdir(parents=True, exist_ok=True)
        args.receipt.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(result, sort_keys=True))


if __name__ == "__main__":
    main()
