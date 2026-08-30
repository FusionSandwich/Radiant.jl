#!/usr/bin/env python3
"""Build a hash-bound Radiant HTS physical-reference HDF5 bundle.

The input JSON specification names immutable source, result, geometry, material-state, and
normalization artifacts. Numerical arrays may be NumPy ``.npy`` files or numeric CSV files. The
output schema exactly matches ``radiant.hts.physical_reference_bundle/v2`` consumed by Julia.

This adapter never infers units, scoring semantics, classifications, normalizations, or
uncertainties. Every response must provide them explicitly, either directly or through the
specification's ``defaults`` object.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import re
from pathlib import Path
from typing import Any

SCHEMA = "radiant.hts.physical_reference_bundle/v2"
SPEC_SCHEMA = "radiant.hts.physical_reference_builder_spec/v2"
STORAGE_ORDER = "canonical-c-row-major-flat/v1"
ALLOWED_CLASSIFICATIONS = {
    "analytic",
    "synthetic",
    "continuous_energy_openmc",
    "geant4",
    "opensn",
    "experiment",
}
PHYSICAL_CLASSIFICATIONS = {
    "continuous_energy_openmc",
    "geant4",
    "opensn",
    "experiment",
}
SHA256_PATTERN = re.compile(r"^[0-9a-f]{64}$")


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def require_text(mapping: dict[str, Any], key: str, *, context: str) -> str:
    value = str(mapping.get(key, "")).strip()
    if not value:
        raise ValueError(f"{context}: missing required field {key}")
    return value


def require_sha256(mapping: dict[str, Any], key: str, *, context: str) -> str:
    value = require_text(mapping, key, context=context).lower()
    if not SHA256_PATTERN.fullmatch(value):
        raise ValueError(f"{context}: {key} must be a 64-character SHA-256 digest")
    return value


def resolve_input_path(specification: Path, value: Any, *, field: str) -> Path:
    text = str(value).strip()
    if not text:
        raise ValueError(f"missing required path field {field}")
    path = Path(text).expanduser()
    if not path.is_absolute():
        path = specification.parent / path
    path = path.resolve()
    if not path.is_file():
        raise FileNotFoundError(f"{field} does not exist: {path}")
    return path


def load_array(path: Path):
    import numpy as np

    if path.suffix.lower() == ".npy":
        array = np.load(path, allow_pickle=False)
    else:
        array = np.loadtxt(path, delimiter=",", comments="#", ndmin=1)
    array = np.asarray(array, dtype=np.float64)
    if array.size == 0 or not np.all(np.isfinite(array)):
        raise ValueError(f"response array is empty or non-finite: {path}")
    return array


def write_string(group: Any, name: str, value: str, string_dtype: Any) -> None:
    import numpy as np

    group.create_dataset(name, data=np.asarray(value, dtype=string_dtype))


def write_string_dict(
    parent: Any,
    name: str,
    values: dict[str, Any],
    string_dtype: Any,
) -> None:
    group = parent.create_group(name)
    entries = sorted((str(key), str(value)) for key, value in values.items())
    group.create_dataset("count", data=len(entries))
    for index, (key, value) in enumerate(entries, start=1):
        entry = group.create_group(f"entry_{index:06d}")
        write_string(entry, "key", key, string_dtype)
        write_string(entry, "value", value, string_dtype)


def merged_field(
    response: dict[str, Any],
    defaults: dict[str, Any],
    key: str,
    *,
    context: str,
) -> str:
    if key in response:
        return require_text(response, key, context=context)
    return require_text(defaults, key, context=f"defaults for {context}")


def merged_sha256(
    response: dict[str, Any],
    defaults: dict[str, Any],
    key: str,
    *,
    context: str,
) -> str:
    if key in response:
        return require_sha256(response, key, context=context)
    return require_sha256(defaults, key, context=f"defaults for {context}")


def build_bundle(specification: Path, output: Path, *, overwrite: bool) -> dict[str, Any]:
    if output.exists() and not overwrite:
        raise FileExistsError(f"refusing to overwrite {output}")
    if output.exists() and not output.is_file():
        raise ValueError(f"output path exists and is not a file: {output}")

    spec = json.loads(specification.read_text(encoding="utf-8"))
    if spec.get("schema") != SPEC_SCHEMA:
        raise ValueError(f"unsupported specification schema: {spec.get('schema')!r}")
    defaults = dict(spec.get("defaults", {}))
    bundle_metadata = dict(spec.get("metadata", {}))
    responses = spec.get("responses", [])
    if not isinstance(responses, list) or not responses:
        raise ValueError("reference specification contains no responses")

    import h5py
    import numpy as np

    string_dtype = h5py.string_dtype(encoding="utf-8")
    prepared: list[dict[str, Any]] = []
    identifiers: set[str] = set()
    process_keys: set[str] = set()

    for index, response_raw in enumerate(responses, start=1):
        if not isinstance(response_raw, dict):
            raise TypeError(f"response {index} must be a JSON object")
        response = dict(response_raw)
        context = f"response {index}"
        response_id = require_text(response, "response_id", context=context)
        if response_id in identifiers:
            raise ValueError(f"duplicate response_id {response_id!r}")
        identifiers.add(response_id)
        process_key = require_text(response, "process_key", context=context)
        if bool(spec.get("require_unique_process_keys", False)) and process_key in process_keys:
            raise ValueError(f"duplicate process_key {process_key!r}")
        process_keys.add(process_key)

        classification = merged_field(
            response, defaults, "classification", context=context
        ).lower()
        if classification not in ALLOWED_CLASSIFICATIONS:
            raise ValueError(f"{context}: unsupported classification {classification!r}")
        particle_tag = merged_field(response, defaults, "particle_tag", context=context)
        units = merged_field(response, defaults, "units", context=context)
        producer = merged_field(response, defaults, "producer", context=context)

        lineage = {
            key: merged_sha256(response, defaults, key, context=context)
            for key in (
                "source_artifact_hash",
                "result_artifact_hash",
                "geometry_hash",
                "material_state_hash",
                "normalization_hash",
            )
        }
        values_path = resolve_input_path(
            specification, response.get("values_path"), field=f"{context}.values_path"
        )
        uncertainty_path = resolve_input_path(
            specification,
            response.get("uncertainty_path"),
            field=f"{context}.uncertainty_path",
        )
        values = load_array(values_path)
        uncertainty = load_array(uncertainty_path)
        if uncertainty.shape != values.shape:
            raise ValueError(
                f"{context}: uncertainty shape {uncertainty.shape} does not match "
                f"values shape {values.shape}"
            )
        if np.any(uncertainty < 0.0):
            raise ValueError(f"{context}: standard uncertainty cannot be negative")

        metadata = dict(response.get("metadata", {}))
        metadata.update(
            {
                "values_path_basename": values_path.name,
                "values_sha256": sha256_file(values_path),
                "uncertainty_path_basename": uncertainty_path.name,
                "uncertainty_sha256": sha256_file(uncertainty_path),
                "scoring_semantics": require_text(
                    response, "scoring_semantics", context=context
                ),
                "physical_reference": str(
                    classification in PHYSICAL_CLASSIFICATIONS
                ).lower(),
            }
        )
        prepared.append(
            {
                "response_id": response_id,
                "process_key": process_key,
                "particle_tag": particle_tag,
                "units": units,
                "classification": classification,
                "producer": producer,
                "lineage": lineage,
                "values": values,
                "uncertainty": uncertainty,
                "metadata": metadata,
            }
        )

    output.parent.mkdir(parents=True, exist_ok=True)
    if output.exists():
        output.unlink()
    with h5py.File(output, "w", libver="earliest") as handle:
        meta = handle.create_group("meta")
        response_group = handle.create_group("responses")
        write_string(meta, "schema", SCHEMA, string_dtype)
        write_string(meta, "storage_order", STORAGE_ORDER, string_dtype)
        meta.create_dataset("response_count", data=np.int64(len(prepared)))
        bundle_metadata.update(
            {
                "builder": "scripts/references/build_physical_reference_bundle.py",
                "builder_spec_sha256": sha256_file(specification),
                "response_count": str(len(prepared)),
            }
        )
        write_string_dict(meta, "metadata", bundle_metadata, string_dtype)

        for index, response in enumerate(
            sorted(prepared, key=lambda item: item["response_id"]), start=1
        ):
            entry = response_group.create_group(f"entry_{index:06d}")
            for key in (
                "response_id",
                "process_key",
                "particle_tag",
                "units",
                "classification",
                "producer",
            ):
                write_string(entry, key, str(response[key]), string_dtype)
            for key, value in response["lineage"].items():
                write_string(entry, key, value, string_dtype)
            values = response["values"]
            uncertainty = response["uncertainty"]
            entry.create_dataset("shape", data=np.asarray(values.shape, dtype=np.int64))
            entry.create_dataset("values_flat", data=values.reshape(-1, order="C"))
            entry.create_dataset(
                "standard_uncertainty_flat",
                data=uncertainty.reshape(-1, order="C"),
            )
            write_string_dict(entry, "metadata", response["metadata"], string_dtype)

    result = {
        "schema": SCHEMA,
        "specification": str(specification.resolve()),
        "specification_sha256": sha256_file(specification),
        "output": str(output.resolve()),
        "output_sha256": sha256_file(output),
        "response_count": len(prepared),
        "response_ids": sorted(identifiers),
        "classifications": sorted({item["classification"] for item in prepared}),
    }
    return result


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("specification", type=Path)
    parser.add_argument("output", type=Path)
    parser.add_argument("--overwrite", action="store_true")
    parser.add_argument("--receipt", type=Path)
    args = parser.parse_args()

    result = build_bundle(
        args.specification.resolve(), args.output.resolve(), overwrite=args.overwrite
    )
    if args.receipt is not None:
        args.receipt.parent.mkdir(parents=True, exist_ok=True)
        args.receipt.write_text(
            json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8"
        )
    print(json.dumps(result, sort_keys=True))


if __name__ == "__main__":
    main()
