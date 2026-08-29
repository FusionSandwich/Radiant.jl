#!/usr/bin/env python3
"""Create deterministic Python/h5py fixtures for Radiant source interchange."""

from __future__ import annotations

import argparse
from pathlib import Path

import h5py
import numpy as np

BOUNDARY_SCHEMA = "radiant.boundary_angular_current/v1"
VOLUME_SCHEMA = "radiant.anisotropic_volume_source/v1"
STORAGE_ORDER = "external-row-major"
UTF8 = h5py.string_dtype(encoding="utf-8")


def _write_string(group: h5py.Group, name: str, value: str) -> None:
    group.create_dataset(name, data=np.asarray(value, dtype=UTF8))


def _write_normalization(
    meta: h5py.Group,
    *,
    basis: str,
    source_rate_per_s: float,
    symmetry_factor: float,
    time_interval_s: tuple[float, float],
    time_class: str,
    source_hash: str,
) -> None:
    _write_string(meta, "normalization_basis", basis)
    meta.create_dataset("source_rate_per_s", data=np.float64(source_rate_per_s))
    meta.create_dataset("symmetry_factor", data=np.float64(symmetry_factor))
    meta.create_dataset("time_interval_s", data=np.asarray(time_interval_s, dtype=np.float64))
    _write_string(meta, "time_class", time_class)
    _write_string(meta, "source_hash", source_hash)


def _write_provenance(handle: h5py.File, values: dict[str, str]) -> None:
    group = handle.create_group("provenance")
    entries = sorted((str(key), str(value)) for key, value in values.items())
    group.create_dataset("keys", data=np.asarray([key for key, _ in entries], dtype=UTF8))
    group.create_dataset("values", data=np.asarray([value for _, value in entries], dtype=UTF8))


def write_boundary_fixture(path: Path) -> None:
    directions = np.asarray(
        [[-1.0, 0.0, 0.0], [-0.5, np.sqrt(3.0) / 2.0, 0.0]], dtype=np.float64
    )
    directional_current_density = np.asarray(
        [[[1.0, 3.0], [0.5, 1.5]]], dtype=np.float64
    )
    variance = np.asarray([[[0.01, 0.04], [0.0025, 0.01]]], dtype=np.float64)

    with h5py.File(path, "w", libver="earliest") as handle:
        meta = handle.create_group("meta")
        patch = handle.create_group("patch")
        energy = handle.create_group("energy")
        angle = handle.create_group("angle")
        source = handle.create_group("source")

        _write_string(meta, "schema", BOUNDARY_SCHEMA)
        _write_string(meta, "particle", "photon")
        _write_string(meta, "source_representation", "directional_current_density")
        _write_string(meta, "storage_order", STORAGE_ORDER)
        _write_normalization(
            meta,
            basis="per_history",
            source_rate_per_s=2.0,
            symmetry_factor=4.0,
            time_interval_s=(0.0, 1.0),
            time_class="prompt",
            source_hash="python-boundary-fixture",
        )

        patch.create_dataset("id", data=np.asarray([11], dtype=np.int64))
        patch.create_dataset("centroid_cm", data=np.asarray([[0.0, 0.0, 0.0]]))
        patch.create_dataset("area_cm2", data=np.asarray([2.0]))
        patch.create_dataset("normal", data=np.asarray([[1.0, 0.0, 0.0]]))
        patch.create_dataset("tangent_1", data=np.asarray([[0.0, 1.0, 0.0]]))
        patch.create_dataset("tangent_2", data=np.asarray([[0.0, 0.0, 1.0]]))

        energy.create_dataset("edges_eV", data=np.asarray([1.0, 10.0, 100.0]))
        angle.create_dataset("direction_cosines", data=directions)
        angle.create_dataset("quadrature_weights", data=np.asarray([0.25, 0.75]))
        source.create_dataset("directional_current_density", data=directional_current_density)
        source.create_dataset("variance", data=variance)
        source.create_dataset("incoming_current", data=np.asarray([[8.0, 4.0]]))
        _write_provenance(
            handle,
            {"producer": "python-h5py", "surface_semantics": "incoming-current"},
        )


def write_volume_fixture(path: Path) -> None:
    values = np.asarray([[[1.0], [2.0]], [[3.0], [4.0]]], dtype=np.float64)
    with h5py.File(path, "w", libver="earliest") as handle:
        meta = handle.create_group("meta")
        voxel = handle.create_group("voxel")
        energy = handle.create_group("energy")
        handle.create_group("angle")
        source = handle.create_group("source")

        _write_string(meta, "schema", VOLUME_SCHEMA)
        _write_string(meta, "particle", "photon")
        _write_string(meta, "angular_representation", "isotropic")
        _write_string(meta, "storage_order", STORAGE_ORDER)
        _write_string(meta, "parent_reaction", "delayed-photon")
        _write_normalization(
            meta,
            basis="per_second",
            source_rate_per_s=1.0,
            symmetry_factor=1.0,
            time_interval_s=(100.0, 200.0),
            time_class="delayed",
            source_hash="python-volume-fixture",
        )

        voxel.create_dataset("id", data=np.asarray([1, 2], dtype=np.int64))
        voxel.create_dataset("volume_cm3", data=np.asarray([1.0, 2.0]))
        energy.create_dataset("edges_eV", data=np.asarray([1.0, 10.0, 100.0]))
        source.create_dataset("values", data=values)
        source.create_dataset("integrated_rate", data=np.asarray([7.0, 10.0]))
        _write_provenance(
            handle,
            {"producer": "python-h5py", "cooling_time_s": "100.0"},
        )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("output_directory", type=Path)
    args = parser.parse_args()
    args.output_directory.mkdir(parents=True, exist_ok=True)
    write_boundary_fixture(args.output_directory / "python_boundary_source.h5")
    write_volume_fixture(args.output_directory / "python_volume_source.h5")


if __name__ == "__main__":
    main()
