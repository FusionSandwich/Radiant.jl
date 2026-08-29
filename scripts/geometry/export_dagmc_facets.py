#!/usr/bin/env python3
"""Export DAGMC/MOAB triangle facets to Radiant's normalized faceted HDF5 schema.

This script intentionally uses pymoab rather than reverse-engineering H5M internals. It preserves
MOAB GLOBAL_ID, GEOM_DIMENSION, CATEGORY, and NAME information where available. A produced file is
an interchange artifact; it does not certify that the DAGMC model is watertight or overlap-free.
Radiant performs independent topology checks when reading it.
"""

from __future__ import annotations

import argparse
import hashlib
from pathlib import Path
from typing import Iterable

import h5py
import numpy as np

SCHEMA = "radiant.faceted_surface_mesh/v1"
STORAGE_ORDER = "external-row-major"
UTF8 = h5py.string_dtype(encoding="utf-8")


def _tag(mb, name: str):
    try:
        return mb.tag_get_handle(name)
    except Exception:
        return None


def _tag_scalar(mb, tag, entity, default: int) -> int:
    if tag is None:
        return default
    try:
        value = np.asarray(mb.tag_get_data(tag, entity, flat=True)).reshape(-1)
        return int(value[0]) if value.size else default
    except Exception:
        return default


def _tag_text(mb, tag, entity, default: str) -> str:
    if tag is None:
        return default
    try:
        value = mb.tag_get_data(tag, entity, flat=True)
        raw = np.asarray(value).reshape(-1)[0]
        if isinstance(raw, bytes):
            return raw.decode(errors="replace").replace("\x00", "").strip() or default
        return str(raw).replace("\x00", "").strip() or default
    except Exception:
        return default


def _material_name(name: str, volume_id: int) -> str:
    cleaned = name.strip()
    if cleaned.startswith("mat:"):
        cleaned = cleaned[4:]
    if cleaned.endswith("_comp"):
        cleaned = cleaned[:-5]
    return cleaned or f"volume:{volume_id}"


def export_h5m(source: Path, output: Path, length_multiplier_to_cm: float) -> None:
    try:
        from pymoab import core, types
    except ImportError as exc:
        raise SystemExit(
            "pymoab is required to read DAGMC .h5m files. Install a MOAB/pymoab build "
            "matching the DAGMC environment, then rerun this adapter."
        ) from exc

    mb = core.Core()
    mb.load_file(str(source))
    root = mb.get_root_set()

    vertices = list(mb.get_entities_by_dimension(root, 0))
    triangles = list(mb.get_entities_by_type(root, types.MBTRI))
    if not vertices or not triangles:
        raise RuntimeError("The H5M file contains no vertices or triangle entities.")

    vertex_index = {entity: index + 1 for index, entity in enumerate(vertices)}
    coordinates = np.asarray(mb.get_coords(vertices), dtype=np.float64).reshape(-1, 3)
    coordinates *= float(length_multiplier_to_cm)
    connectivity = np.empty((len(triangles), 3), dtype=np.int64)
    for index, triangle in enumerate(triangles):
        nodes = list(mb.get_connectivity(triangle))
        if len(nodes) != 3:
            raise RuntimeError(f"MOAB entity {triangle} is not a three-node triangle.")
        connectivity[index, :] = [vertex_index[node] for node in nodes]

    geom_dimension = _tag(mb, "GEOM_DIMENSION")
    global_id = _tag(mb, "GLOBAL_ID")
    name_tag = _tag(mb, "NAME")

    surface_sets = list(
        mb.get_entities_by_type_and_tag(
            root,
            types.MBENTITYSET,
            [geom_dimension] if geom_dimension is not None else [],
            [np.asarray([2], dtype=np.int32)] if geom_dimension is not None else [],
        )
    ) if geom_dimension is not None else []
    volume_sets = list(
        mb.get_entities_by_type_and_tag(
            root,
            types.MBENTITYSET,
            [geom_dimension] if geom_dimension is not None else [],
            [np.asarray([3], dtype=np.int32)] if geom_dimension is not None else [],
        )
    ) if geom_dimension is not None else []

    triangle_lookup = {entity: index for index, entity in enumerate(triangles)}
    surface_ids = np.zeros(len(triangles), dtype=np.int64)
    volume_ids = np.zeros(len(triangles), dtype=np.int64)
    material_tags = np.full(len(triangles), "unbound", dtype=object)

    volume_set_id = {
        entity: _tag_scalar(mb, global_id, entity, index + 1)
        for index, entity in enumerate(volume_sets)
    }
    volume_set_material = {
        entity: _material_name(
            _tag_text(mb, name_tag, entity, f"volume:{volume_set_id[entity]}"),
            volume_set_id[entity],
        )
        for entity in volume_sets
    }

    for fallback_surface_id, surface_set in enumerate(surface_sets, start=1):
        surface_id = _tag_scalar(mb, global_id, surface_set, fallback_surface_id)
        surface_triangles = list(mb.get_entities_by_type(surface_set, types.MBTRI))
        parents = list(mb.get_parent_meshsets(surface_set))
        parent_volumes = [parent for parent in parents if parent in volume_set_id]
        if len(parent_volumes) > 2:
            raise RuntimeError(
                f"Surface {surface_id} has more than two adjacent DAGMC volumes."
            )
        chosen_volume = parent_volumes[0] if parent_volumes else None
        chosen_volume_id = volume_set_id.get(chosen_volume, surface_id)
        chosen_material = volume_set_material.get(chosen_volume, f"surface:{surface_id}")
        for triangle in surface_triangles:
            if triangle not in triangle_lookup:
                continue
            index = triangle_lookup[triangle]
            if surface_ids[index] != 0 and surface_ids[index] != surface_id:
                raise RuntimeError(
                    f"Triangle {triangle} belongs to multiple geometric surfaces."
                )
            surface_ids[index] = surface_id
            volume_ids[index] = chosen_volume_id
            material_tags[index] = chosen_material

    # Files without DAGMC geometry sets remain importable but explicitly unbound.
    surface_ids[surface_ids == 0] = np.arange(1, len(triangles) + 1)[surface_ids == 0]
    volume_ids[volume_ids == 0] = 1

    source_sha256 = hashlib.sha256(source.read_bytes()).hexdigest()
    output.parent.mkdir(parents=True, exist_ok=True)
    with h5py.File(output, "w") as handle:
        meta = handle.create_group("meta")
        geometry = handle.create_group("geometry")
        provenance = handle.create_group("provenance")
        meta.create_dataset("schema", data=np.asarray(SCHEMA, dtype=UTF8))
        meta.create_dataset("storage_order", data=np.asarray(STORAGE_ORDER, dtype=UTF8))
        meta.create_dataset("length_unit", data=np.asarray("cm", dtype=UTF8))
        meta.create_dataset("geometry_hash", data=np.asarray(source_sha256, dtype=UTF8))
        geometry.create_dataset("vertices_cm", data=coordinates)
        geometry.create_dataset("triangles", data=connectivity)
        geometry.create_dataset("surface_ids", data=surface_ids)
        geometry.create_dataset("volume_ids", data=volume_ids)
        geometry.create_dataset(
            "material_tags", data=np.asarray(material_tags.tolist(), dtype=UTF8)
        )
        keys = [
            "source_format",
            "source_path",
            "source_sha256",
            "length_multiplier_to_cm",
            "adapter",
        ]
        values = [
            "DAGMC-H5M",
            str(source.resolve()),
            source_sha256,
            repr(float(length_multiplier_to_cm)),
            "scripts/geometry/export_dagmc_facets.py",
        ]
        provenance.create_dataset("keys", data=np.asarray(keys, dtype=UTF8))
        provenance.create_dataset("values", data=np.asarray(values, dtype=UTF8))

    digest = hashlib.sha256(output.read_bytes()).hexdigest()
    print(f"output={output}")
    print(f"source_sha256={source_sha256}")
    print(f"output_sha256={digest}")
    print(f"vertices={len(vertices)}")
    print(f"triangles={len(triangles)}")
    print(f"surfaces={len(set(surface_ids.tolist()))}")
    print(f"volumes={len(set(volume_ids.tolist()))}")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("source_h5m", type=Path)
    parser.add_argument("output_hdf5", type=Path)
    parser.add_argument(
        "--length-multiplier-to-cm",
        type=float,
        default=1.0,
        help="Multiply MOAB coordinates by this value to obtain centimetres.",
    )
    args = parser.parse_args()
    if not np.isfinite(args.length_multiplier_to_cm) or args.length_multiplier_to_cm <= 0:
        parser.error("--length-multiplier-to-cm must be finite and positive")
    export_h5m(args.source_h5m, args.output_hdf5, args.length_multiplier_to_cm)


if __name__ == "__main__":
    main()
