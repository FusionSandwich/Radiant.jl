#!/usr/bin/env python3
"""Create deterministic Python/h5py fixtures for the Radiant HDF5 source contract."""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Mapping

import h5py
import numpy as np

BOUNDARY_SCHEMA = "radiant.boundary_angular_current/v1"
VOLUME_SCHEMA = "radiant.anisotropic_volume_source/v1"
STORAGE_ORDER = "canonical-c-row-major-flat/v1"
UTF8 = h5py.string_dtype(encoding="utf-8")


def _write_string(group: h5py.Group, name: str, value: str) -> None:
    group.create_dataset(name, data=np.asarray(value, dtype=UTF8))


def _write_tensor(group: h5py.Group, name: str, value: np.ndarray) -> None:
    array = np.asarray(value)
    group.create_dataset(f"{name}_flat", data=array.reshape(-1, order="C"))
    group.create_dataset(f"{name}_shape", data=np.asarray(array.shape, dtype=np.int64))


def _write_optional_tensor(group: h5py.Group, name: str, value: np.ndarray | None) -> None:
    group.create_dataset(f"{name}_present", data=np.int8(value is not None))
    if value is not None:
        _write_tensor(group, name, value)


def _write_dictionary(group: h5py.Group, name: str, values: Mapping[str, str]) -> None:
    target = group.create_group(name)
    entries = sorted((str(key), str(value)) for key, value in values.items())
    target.create_dataset("count", data=np.int64(len(entries)))
    for index, (key, value) in enumerate(entries, start=1):
        entry = target.create_group(f"entry_{index:06d}")
        _write_string(entry, "key", key)
        _write_string(entry, "value", value)


def _write_particle(meta: h5py.Group, tag: str, mass: float, charge: float) -> None:
    particle = meta.create_group("particle")
    _write_string(particle, "tag", tag)
    particle.create_dataset("mass_MeV_c2", data=np.float64(mass))
    particle.create_dataset("charge_e", data=np.float64(charge))


def _write_normalization(
    meta: h5py.Group,
    *,
    basis: str,
    source_rate_per_s: float,
    symmetry_factor: float,
    time_interval_s: tuple[float, float],
    time_class: str,
    source_hash: str,
    provenance: Mapping[str, str],
) -> None:
    normalization = meta.create_group("normalization")
    _write_string(normalization, "basis", basis)
    normalization.create_dataset("source_rate_per_s", data=np.float64(source_rate_per_s))
    normalization.create_dataset("symmetry_factor", data=np.float64(symmetry_factor))
    normalization.create_dataset("time_interval_s", data=np.asarray(time_interval_s, dtype=np.float64))
    _write_string(normalization, "time_class", time_class)
    _write_string(normalization, "source_hash", source_hash)
    _write_dictionary(normalization, "provenance", provenance)


def write_boundary_fixture(path: Path) -> None:
    directions = np.asarray(
        [[-1.0, 0.0, 0.0], [-0.5, np.sqrt(3.0) / 2.0, 0.0]], dtype=np.float64
    )
    directional_current_density = np.asarray([[[1.0, 3.0], [0.5, 1.5]]], dtype=np.float64)
    variance = np.asarray([[[0.01, 0.04], [0.0025, 0.01]]], dtype=np.float64)

    with h5py.File(path, "w", libver="earliest") as handle:
        meta = handle.create_group("meta")
        _write_string(meta, "schema", BOUNDARY_SCHEMA)
        _write_string(meta, "storage_order", STORAGE_ORDER)
        _write_string(meta, "representation", "directional_current_density")
        _write_particle(meta, "photon", 0.0, 0.0)
        _write_normalization(
            meta,
            basis="per_history",
            source_rate_per_s=2.0,
            symmetry_factor=4.0,
            time_interval_s=(0.0, 1.0),
            time_class="prompt",
            source_hash="python-boundary-fixture",
            provenance={"producer": "python-h5py", "fixture": "analytic"},
        )
        _write_dictionary(meta, "provenance", {"surface_semantics": "incoming-current"})

        patch = handle.create_group("patch")
        patch.create_dataset("id", data=np.asarray([11], dtype=np.int64))
        _write_tensor(patch, "centroid_cm", np.asarray([[0.0, 0.0, 0.0]]))
        patch.create_dataset("area_cm2", data=np.asarray([2.0]))
        _write_tensor(patch, "normal", np.asarray([[1.0, 0.0, 0.0]]))
        _write_tensor(patch, "tangent_1", np.asarray([[0.0, 1.0, 0.0]]))
        _write_tensor(patch, "tangent_2", np.asarray([[0.0, 0.0, 1.0]]))

        phase = handle.create_group("phase_space")
        phase.create_dataset("energy_edges_eV", data=np.asarray([1.0, 10.0, 100.0]))
        _write_tensor(phase, "direction", directions)
        phase.create_dataset("quadrature_weight", data=np.asarray([0.25, 0.75]))
        _write_tensor(phase, "value", directional_current_density)
        _write_optional_tensor(phase, "variance", variance)


def write_volume_fixture(path: Path) -> None:
    values = np.asarray([[[1.0], [2.0]], [[3.0], [4.0]]], dtype=np.float64)
    with h5py.File(path, "w", libver="earliest") as handle:
        meta = handle.create_group("meta")
        _write_string(meta, "schema", VOLUME_SCHEMA)
        _write_string(meta, "storage_order", STORAGE_ORDER)
        _write_string(meta, "angular_representation", "isotropic")
        _write_string(meta, "parent_reaction", "delayed-photon")
        _write_particle(meta, "photon", 0.0, 0.0)
        _write_normalization(
            meta,
            basis="per_second",
            source_rate_per_s=1.0,
            symmetry_factor=1.0,
            time_interval_s=(100.0, 200.0),
            time_class="delayed",
            source_hash="python-volume-fixture",
            provenance={"producer": "python-h5py", "fixture": "analytic"},
        )
        _write_dictionary(meta, "provenance", {"cooling_time_s": "100.0"})

        volume = handle.create_group("volume")
        volume.create_dataset("voxel_id", data=np.asarray([1, 2], dtype=np.int64))
        volume.create_dataset("voxel_volume_cm3", data=np.asarray([1.0, 2.0]))

        phase = handle.create_group("phase_space")
        phase.create_dataset("energy_edges_eV", data=np.asarray([1.0, 10.0, 100.0]))
        _write_tensor(phase, "value", values)
        _write_optional_tensor(phase, "variance", None)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("output_directory", type=Path)
    args = parser.parse_args()
    args.output_directory.mkdir(parents=True, exist_ok=True)
    write_boundary_fixture(args.output_directory / "python_boundary_source.h5")
    write_volume_fixture(args.output_directory / "python_volume_source.h5")


if __name__ == "__main__":
    main()
