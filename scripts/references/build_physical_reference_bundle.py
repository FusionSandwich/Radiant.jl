#!/usr/bin/env python3
"""Build `radiant.hts.physical_reference_bundle/v1` from independent response arrays.

The JSON specification names immutable source/geometry/material/model hashes and one or more
response fields. Array files may be `.npy` or numeric CSV. h5py and NumPy are imported lazily so the
script can be syntax-checked without optional scientific packages.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
from typing import Any

SCHEMA = "radiant.hts.physical_reference_bundle/v1"
STORAGE_ORDER = "canonical-c-row-major-flat/v1"
ALLOWED_CLASSIFICATIONS = {
    "physical",
    "continuous_energy_monte_carlo",
    "deterministic_reference",
    "analytic",
    "software_fixture",
}


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def load_array(path: Path):
    import numpy as np

    if path.suffix.lower() == ".npy":
        array = np.load(path, allow_pickle=False)
    else:
        array = np.loadtxt(path, delimiter=",", comments="#")
    array = np.asarray(array, dtype=np.float64)
    if array.size == 0 or not np.all(np.isfinite(array)):
        raise ValueError(f"response array is empty or non-finite: {path}")
    return array


def write_string_dict(parent, name: str, values: dict[str, Any]) -> None:
    group = parent.create_group(name)
    keys = sorted(str(key) for key in values)
    group.create_dataset("count", data=len(keys))
    for index, key in enumerate(keys, start=1):
        entry = group.create_group(f"entry_{index:06d}")
        entry.create_dataset("key", data=key)
        entry.create_dataset("value", data=str(values[key]))


def require_hash(spec: dict[str, Any], name: str) -> str:
    value = str(spec.get(name, ""))
    if not value or value == "unbound":
        raise ValueError(f"{name} must be an immutable non-unbound identifier")
    return value


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("specification", type=Path)
    parser.add_argument("output", type=Path)
    parser.add_argument("--overwrite", action="store_true")
    args = parser.parse_args()

    if args.output.exists() and not args.overwrite:
        raise SystemExit(f"refusing to overwrite {args.output}")
    spec = json.loads(args.specification.read_text(encoding="utf-8"))
    required_text = ("case_id", "producer_code", "producer_version")
    for key in required_text:
        if not str(spec.get(key, "")):
            raise ValueError(f"missing required specification field {key}")
    lineage = {
        "source_hash": require_hash(spec, "source_hash"),
        "geometry_hash": require_hash(spec, "geometry_hash"),
        "material_state_hash": require_hash(spec, "material_state_hash"),
        "model_hash": require_hash(spec, "model_hash"),
    }
    responses = spec.get("responses", [])
    if not responses:
        raise ValueError("reference specification contains no responses")

    import h5py
    import numpy as np

    args.output.parent.mkdir(parents=True, exist_ok=True)
    with h5py.File(args.output, "w", libver="earliest") as handle:
        meta = handle.create_group("meta")
        meta.create_dataset("schema", data=SCHEMA)
        meta.create_dataset("storage_order", data=STORAGE_ORDER)
        meta.create_dataset("case_id", data=str(spec["case_id"]))
        meta.create_dataset("producer_code", data=str(spec["producer_code"]))
        meta.create_dataset("producer_version", data=str(spec["producer_version"]))
        for key, value in lineage.items():
            meta.create_dataset(key, data=value)
        write_string_dict(meta, "metadata", spec.get("metadata", {}))

        response_group = handle.create_group("responses")
        response_group.create_dataset("count", data=len(responses))
        identifiers: set[str] = set()
        for index, response in enumerate(responses, start=1):
            response_id = str(response.get("response_id", ""))
            if not response_id or response_id in identifiers:
                raise ValueError(f"invalid or duplicate response_id {response_id!r}")
            identifiers.add(response_id)
            classification = str(response.get("classification", ""))
            if classification not in ALLOWED_CLASSIFICATIONS:
                raise ValueError(f"unsupported response classification {classification}")
            values_path = Path(response["values_path"])
            values = load_array(values_path)
            uncertainty_path = response.get("uncertainty_path")
            uncertainty = load_array(Path(uncertainty_path)) if uncertainty_path else None
            if uncertainty is not None:
                if uncertainty.shape != values.shape or np.any(uncertainty < 0):
                    raise ValueError(f"invalid uncertainty array for {response_id}")
            group = response_group.create_group(f"response_{index:06d}")
            group.create_dataset("response_id", data=response_id)
            group.create_dataset("units", data=str(response["units"]))
            group.create_dataset("classification", data=classification)
            group.create_dataset("scoring_semantics", data=str(response["scoring_semantics"]))
            group.create_dataset("values_flat", data=values.reshape(-1, order="C"))
            group.create_dataset("values_shape", data=np.asarray(values.shape, dtype=np.int64))
            group.create_dataset("uncertainty_present", data=np.int8(uncertainty is not None))
            if uncertainty is not None:
                group.create_dataset("uncertainty_flat", data=uncertainty.reshape(-1, order="C"))
                group.create_dataset("uncertainty_shape",
                                     data=np.asarray(uncertainty.shape, dtype=np.int64))
            metadata = dict(response.get("metadata", {}))
            metadata.update(lineage)
            metadata["values_sha256"] = sha256_file(values_path)
            if uncertainty_path:
                metadata["uncertainty_sha256"] = sha256_file(Path(uncertainty_path))
            write_string_dict(group, "metadata", metadata)

    print(json.dumps({
        "path": str(args.output.resolve()),
        "sha256": sha256_file(args.output),
        "schema": SCHEMA,
        "response_count": len(responses),
    }, sort_keys=True))


if __name__ == "__main__":
    main()
