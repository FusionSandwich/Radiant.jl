#!/usr/bin/env python3
"""Convert OpenMC surface-source banks to a Radiant angular-current source.

The adapter reads OpenMC's ``source_bank`` records directly with h5py. Each crossing weight is
interpreted as an integrated surface-current contribution per originating source history. It is
therefore binned as current before conversion to directional current density; it is never multiplied
by another ``|u·n|`` factor.

The JSON input schema is ``radiant.openmc_surface_bank_adapter/v1``. The output is the existing
``radiant.boundary_angular_current/v1`` HDF5 contract with
``source_representation=directional_current_density``.

No OpenMC installation is required. The adapter supports both the legacy particle codes
0/1/2/3 and the PDG-backed codes used by newer OpenMC releases. It does not infer physical source
rate, patch area, outward normal, source histories, or independent-batch identity.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Iterable

import h5py
import numpy as np

SPEC_SCHEMA = "radiant.openmc_surface_bank_adapter/v1"
OUTPUT_SCHEMA = "radiant.boundary_angular_current/v1"
STORAGE_ORDER = "external-row-major"
UTF8 = h5py.string_dtype(encoding="utf-8")
SHA256_PATTERN = re.compile(r"^[0-9a-f]{64}$")
PARTICLE_CODES: dict[str, set[int]] = {
    "neutron": {0, 2112},
    "photon": {1, 22},
    "electron": {2, 11},
    "positron": {3, -11},
}
ALLOWED_CLASSIFICATIONS = {
    "software-verification",
    "physical-data-candidate",
    "physical-data-qualified",
}


@dataclass(frozen=True)
class Patch:
    patch_id: int
    surface_ids: tuple[int, ...]
    centroid_cm: np.ndarray
    area_cm2: float
    normal: np.ndarray
    tangent_1: np.ndarray
    tangent_2: np.ndarray
    selection_type: str
    selection: dict[str, Any]

    def contains(self, position_cm: np.ndarray, tolerance_cm: float) -> bool:
        offset = position_cm - self.centroid_cm
        coordinate_1 = float(np.dot(offset, self.tangent_1))
        coordinate_2 = float(np.dot(offset, self.tangent_2))
        normal_offset = abs(float(np.dot(offset, self.normal)))
        if normal_offset > tolerance_cm:
            return False
        if self.selection_type == "all":
            return True
        if self.selection_type == "local_rectangle":
            u_min, u_max = self.selection["u_bounds_cm"]
            v_min, v_max = self.selection["v_bounds_cm"]
            return (
                u_min - tolerance_cm <= coordinate_1 <= u_max + tolerance_cm
                and v_min - tolerance_cm <= coordinate_2 <= v_max + tolerance_cm
            )
        if self.selection_type == "local_circle":
            radius_cm = float(self.selection["radius_cm"])
            return coordinate_1**2 + coordinate_2**2 <= (
                radius_cm + tolerance_cm
            ) ** 2
        if self.selection_type == "local_triangle":
            vertices = np.asarray(self.selection["vertices_uv_cm"], dtype=np.float64)
            point = np.asarray([coordinate_1, coordinate_2], dtype=np.float64)
            return point_in_triangle(point, vertices, tolerance_cm)
        raise ValueError(f"unsupported patch selection type {self.selection_type!r}")


@dataclass
class BankResult:
    batch_id: str
    path: Path
    file_sha256: str
    source_histories: int
    independent: bool
    integrated_current_per_history: np.ndarray
    accepted_contribution_weights: list[float]
    angular_error_weighted_sum: float
    angular_error_weight_sum: float
    maximum_angular_error_rad: float
    counts: dict[str, int]


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def canonical_sha256(value: Any) -> str:
    encoded = json.dumps(value, sort_keys=True, separators=(",", ":")).encode("utf-8")
    return hashlib.sha256(encoded).hexdigest()


def require_text(mapping: dict[str, Any], key: str, *, context: str) -> str:
    value = str(mapping.get(key, "")).strip()
    if not value:
        raise ValueError(f"{context}: missing required field {key}")
    return value


def require_positive_float(mapping: dict[str, Any], key: str, *, context: str) -> float:
    value = float(mapping.get(key, math.nan))
    if not math.isfinite(value) or value <= 0.0:
        raise ValueError(f"{context}: {key} must be finite and positive")
    return value


def normalized_vector(value: Any, *, label: str) -> np.ndarray:
    vector = np.asarray(value, dtype=np.float64)
    if vector.shape != (3,) or not np.all(np.isfinite(vector)):
        raise ValueError(f"{label} must be a finite three-vector")
    magnitude = float(np.linalg.norm(vector))
    if magnitude <= 0.0:
        raise ValueError(f"{label} cannot be zero")
    return vector / magnitude


def point_in_triangle(point: np.ndarray, vertices: np.ndarray, tolerance: float) -> bool:
    if vertices.shape != (3, 2) or not np.all(np.isfinite(vertices)):
        raise ValueError("local_triangle vertices_uv_cm must have shape (3,2)")
    a, b, c = vertices
    v0 = c - a
    v1 = b - a
    v2 = point - a
    denominator = v0[0] * v1[1] - v1[0] * v0[1]
    if abs(denominator) <= np.finfo(np.float64).eps:
        raise ValueError("local_triangle selection is degenerate")
    u = (v2[0] * v1[1] - v1[0] * v2[1]) / denominator
    v = (v0[0] * v2[1] - v2[0] * v0[1]) / denominator
    scale = max(float(np.max(np.abs(vertices))), 1.0)
    epsilon = tolerance / scale
    return u >= -epsilon and v >= -epsilon and u + v <= 1.0 + epsilon


def parse_patch(raw: dict[str, Any], *, frame_tolerance: float) -> Patch:
    patch_id = int(raw.get("id", 0))
    if patch_id <= 0:
        raise ValueError("patch id must be positive")
    surface_ids = tuple(int(value) for value in raw.get("surface_ids", []))
    if not surface_ids or any(value <= 0 for value in surface_ids):
        raise ValueError(f"patch {patch_id}: surface_ids must be positive and nonempty")
    if len(set(surface_ids)) != len(surface_ids):
        raise ValueError(f"patch {patch_id}: surface_ids must be unique")
    centroid = np.asarray(raw.get("centroid_cm"), dtype=np.float64)
    if centroid.shape != (3,) or not np.all(np.isfinite(centroid)):
        raise ValueError(f"patch {patch_id}: centroid_cm must be a finite three-vector")
    area = require_positive_float(raw, "area_cm2", context=f"patch {patch_id}")
    normal = normalized_vector(raw.get("normal"), label=f"patch {patch_id} normal")
    tangent_1 = normalized_vector(
        raw.get("tangent_1"), label=f"patch {patch_id} tangent_1"
    )
    tangent_2 = normalized_vector(
        raw.get("tangent_2"), label=f"patch {patch_id} tangent_2"
    )
    gram = np.asarray(
        [
            [np.dot(normal, normal), np.dot(normal, tangent_1), np.dot(normal, tangent_2)],
            [np.dot(tangent_1, normal), np.dot(tangent_1, tangent_1), np.dot(tangent_1, tangent_2)],
            [np.dot(tangent_2, normal), np.dot(tangent_2, tangent_1), np.dot(tangent_2, tangent_2)],
        ]
    )
    if not np.allclose(gram, np.eye(3), rtol=frame_tolerance, atol=frame_tolerance):
        raise ValueError(f"patch {patch_id}: normal and tangents must be orthonormal")
    if float(np.dot(np.cross(tangent_1, tangent_2), normal)) < 1.0 - 10.0 * frame_tolerance:
        raise ValueError(
            f"patch {patch_id}: frame must be right-handed (tangent_1 x tangent_2 = normal)"
        )

    selection = dict(raw.get("selection", {"type": "all"}))
    selection_type = str(selection.get("type", "all")).strip().lower()
    if selection_type == "local_rectangle":
        for key in ("u_bounds_cm", "v_bounds_cm"):
            bounds = np.asarray(selection.get(key), dtype=np.float64)
            if bounds.shape != (2,) or not np.all(np.isfinite(bounds)) or bounds[1] <= bounds[0]:
                raise ValueError(f"patch {patch_id}: {key} must be increasing finite bounds")
            selection[key] = bounds.tolist()
        expected_area = (
            selection["u_bounds_cm"][1] - selection["u_bounds_cm"][0]
        ) * (selection["v_bounds_cm"][1] - selection["v_bounds_cm"][0])
        if not math.isclose(expected_area, area, rel_tol=1.0e-8, abs_tol=1.0e-14):
            raise ValueError(f"patch {patch_id}: rectangle bounds do not reproduce area_cm2")
    elif selection_type == "local_circle":
        radius = float(selection.get("radius_cm", math.nan))
        if not math.isfinite(radius) or radius <= 0.0:
            raise ValueError(f"patch {patch_id}: circle radius must be positive")
        expected_area = math.pi * radius**2
        if not math.isclose(expected_area, area, rel_tol=1.0e-8, abs_tol=1.0e-14):
            raise ValueError(f"patch {patch_id}: circle radius does not reproduce area_cm2")
    elif selection_type == "local_triangle":
        vertices = np.asarray(selection.get("vertices_uv_cm"), dtype=np.float64)
        point_in_triangle(vertices[0], vertices, 0.0)  # validates shape and degeneracy
        expected_area = 0.5 * abs(
            np.cross(vertices[1] - vertices[0], vertices[2] - vertices[0])
        )
        if not math.isclose(float(expected_area), area, rel_tol=1.0e-8, abs_tol=1.0e-14):
            raise ValueError(f"patch {patch_id}: triangle vertices do not reproduce area_cm2")
        selection["vertices_uv_cm"] = vertices.tolist()
    elif selection_type != "all":
        raise ValueError(f"patch {patch_id}: unsupported selection type {selection_type!r}")
    return Patch(
        patch_id,
        surface_ids,
        centroid,
        area,
        normal,
        tangent_1,
        tangent_2,
        selection_type,
        selection,
    )


def resolve_path(specification_path: Path, value: Any, *, label: str) -> Path:
    path = Path(str(value)).expanduser()
    if not path.is_absolute():
        path = specification_path.parent / path
    path = path.resolve()
    if not path.is_file():
        raise FileNotFoundError(f"{label} does not exist: {path}")
    return path


def _structured_component(array: np.ndarray, names: Iterable[str]) -> np.ndarray | None:
    dtype_names = array.dtype.names or ()
    for name in names:
        if name in dtype_names:
            return np.asarray(array[name], dtype=np.float64)
    return None


def extract_vector(records: np.ndarray, field_name: str) -> np.ndarray:
    names = records.dtype.names or ()
    if field_name in names:
        field = records[field_name]
        nested_names = field.dtype.names or ()
        if nested_names:
            components = []
            for options in (("x", "u"), ("y", "v"), ("z", "w")):
                component = _structured_component(field, options)
                if component is None:
                    raise ValueError(
                        f"source_bank field {field_name!r} does not contain x/y/z components"
                    )
                components.append(component)
            return np.column_stack(components)
        numeric = np.asarray(field, dtype=np.float64)
        if numeric.ndim == 2 and numeric.shape[1] == 3:
            return numeric
        raise ValueError(f"source_bank field {field_name!r} is not a three-vector")

    fallback = ("x", "y", "z") if field_name == "r" else ("u", "v", "w")
    if all(name in names for name in fallback):
        return np.column_stack([np.asarray(records[name], dtype=np.float64) for name in fallback])
    raise ValueError(f"source_bank is missing vector field {field_name!r}")


def extract_scalar(records: np.ndarray, candidates: Iterable[str], *, required: bool = True):
    names = records.dtype.names or ()
    for name in candidates:
        if name in names:
            return np.asarray(records[name])
    if required:
        raise ValueError(f"source_bank is missing all scalar fields {tuple(candidates)!r}")
    return None


def particle_mask(raw_codes: np.ndarray, particle: str) -> np.ndarray:
    if raw_codes.dtype.kind in {"S", "U", "O"}:
        normalized = np.asarray(
            [
                value.decode("utf-8") if isinstance(value, bytes) else str(value)
                for value in raw_codes
            ]
        )
        return np.char.lower(normalized.astype(str)) == particle
    codes = np.asarray(raw_codes, dtype=np.int64)
    return np.isin(codes, sorted(PARTICLE_CODES[particle]))


def energy_group_index(energy_eV: float, edges_eV: np.ndarray) -> int | None:
    if energy_eV < edges_eV[0] or energy_eV > edges_eV[-1]:
        return None
    if energy_eV == edges_eV[-1]:
        return len(edges_eV) - 2
    index = int(np.searchsorted(edges_eV, energy_eV, side="right") - 1)
    return index if 0 <= index < len(edges_eV) - 1 else None


def choose_patch(
    surface_id: int,
    position_cm: np.ndarray,
    patches: list[Patch],
    tolerance_cm: float,
) -> int | None:
    matches = [
        index
        for index, patch in enumerate(patches)
        if surface_id in patch.surface_ids and patch.contains(position_cm, tolerance_cm)
    ]
    if len(matches) > 1:
        identifiers = [patches[index].patch_id for index in matches]
        raise ValueError(
            f"source record maps to multiple patches {identifiers} on surface {surface_id}"
        )
    return matches[0] if matches else None


def choose_direction(
    direction: np.ndarray,
    patch: Patch,
    target_directions: np.ndarray,
    maximum_angle_rad: float,
) -> tuple[int, float] | None:
    incoming = np.flatnonzero(target_directions @ patch.normal < -1.0e-14)
    if incoming.size == 0:
        raise ValueError(f"patch {patch.patch_id} has no incoming target ordinate")
    dots = np.clip(target_directions[incoming] @ direction, -1.0, 1.0)
    local = int(np.argmax(dots))
    target_index = int(incoming[local])
    angle = float(math.acos(float(dots[local])))
    if angle > maximum_angle_rad:
        return None
    return target_index, angle


def policy_action(policies: dict[str, str], key: str) -> str:
    value = str(policies.get(key, "error")).lower()
    if value not in {"error", "discard"}:
        raise ValueError(f"policy {key!r} must be 'error' or 'discard'")
    return value


def handle_rejection(
    policies: dict[str, str],
    key: str,
    counts: dict[str, int],
    message: str,
) -> None:
    counts[key] = counts.get(key, 0) + 1
    if policy_action(policies, key) == "error":
        raise ValueError(message)


def read_bank(
    specification_path: Path,
    bank_spec: dict[str, Any],
    particle: str,
    patches: list[Patch],
    energy_edges_eV: np.ndarray,
    target_directions: np.ndarray,
    maximum_angle_rad: float,
    position_tolerance_cm: float,
    sense_tolerance: float,
    policies: dict[str, str],
    event_time_window_s: tuple[float, float] | None,
) -> BankResult:
    path = resolve_path(specification_path, bank_spec.get("path"), label="bank path")
    expected_hash = require_text(bank_spec, "sha256", context=f"bank {path.name}").lower()
    if not SHA256_PATTERN.fullmatch(expected_hash):
        raise ValueError(f"bank {path.name}: sha256 must be a 64-character digest")
    actual_hash = sha256_file(path)
    if actual_hash != expected_hash:
        raise ValueError(
            f"bank {path.name}: SHA-256 mismatch; expected {expected_hash}, calculated {actual_hash}"
        )
    histories = int(bank_spec.get("source_histories", 0))
    if histories <= 0:
        raise ValueError(f"bank {path.name}: source_histories must be positive")
    batch_id = require_text(bank_spec, "batch_id", context=f"bank {path.name}")
    independent = bool(bank_spec.get("independent", False))
    dataset_path = str(bank_spec.get("dataset", "source_bank"))

    with h5py.File(path, "r") as handle:
        if dataset_path not in handle:
            raise ValueError(f"bank {path.name}: dataset {dataset_path!r} is missing")
        records = handle[dataset_path][...]
    if records.dtype.names is None:
        raise ValueError(f"bank {path.name}: source_bank must be a compound dataset")

    positions = extract_vector(records, "r")
    directions = extract_vector(records, "u")
    energies = np.asarray(extract_scalar(records, ("E", "energy")), dtype=np.float64)
    weights = np.asarray(extract_scalar(records, ("wgt", "weight")), dtype=np.float64)
    surface_ids = np.asarray(
        extract_scalar(records, ("surf_id", "surface_id")), dtype=np.int64
    )
    particle_codes = extract_scalar(records, ("particle", "particle_type"))
    times = extract_scalar(records, ("time",), required=False)
    if times is None:
        times = np.zeros(len(records), dtype=np.float64)
    else:
        times = np.asarray(times, dtype=np.float64)

    count = len(records)
    for name, array in (
        ("positions", positions),
        ("directions", directions),
        ("energies", energies),
        ("weights", weights),
        ("surface_ids", surface_ids),
        ("particle", particle_codes),
        ("time", times),
    ):
        if len(array) != count:
            raise ValueError(f"bank {path.name}: {name} length does not match source_bank")
    if not np.all(np.isfinite(positions)) or not np.all(np.isfinite(directions)):
        raise ValueError(f"bank {path.name}: positions and directions must be finite")
    direction_norm = np.linalg.norm(directions, axis=1)
    if not np.allclose(direction_norm, 1.0, rtol=1.0e-8, atol=1.0e-10):
        raise ValueError(f"bank {path.name}: source directions must be unit vectors")
    if not np.all(np.isfinite(energies)) or np.any(energies < 0.0):
        raise ValueError(f"bank {path.name}: energies must be finite and nonnegative")
    if not np.all(np.isfinite(weights)) or np.any(weights <= 0.0):
        raise ValueError(f"bank {path.name}: weights must be finite and positive")
    if not np.all(np.isfinite(times)):
        raise ValueError(f"bank {path.name}: times must be finite")

    desired_particle = particle_mask(particle_codes, particle)
    current = np.zeros(
        (len(patches), len(energy_edges_eV) - 1, len(target_directions)),
        dtype=np.float64,
    )
    counts = {
        "records_total": count,
        "particle_mismatch": 0,
        "time_out_of_range": 0,
        "unmapped": 0,
        "outgoing": 0,
        "energy_out_of_range": 0,
        "angular_out_of_range": 0,
        "accepted": 0,
    }
    contributions: list[float] = []
    angular_weighted_sum = 0.0
    angular_weight_sum = 0.0
    maximum_angular_error = 0.0

    for index in range(count):
        if not bool(desired_particle[index]):
            counts["particle_mismatch"] += 1
            continue
        if event_time_window_s is not None and not (
            event_time_window_s[0] <= times[index] <= event_time_window_s[1]
        ):
            handle_rejection(
                policies,
                "time_out_of_range",
                counts,
                f"bank {path.name}: event time lies outside the requested interval",
            )
            continue
        patch_index = choose_patch(
            int(surface_ids[index]), positions[index], patches, position_tolerance_cm
        )
        if patch_index is None:
            handle_rejection(
                policies,
                "unmapped",
                counts,
                f"bank {path.name}: source record is not assigned to exactly one patch",
            )
            continue
        patch = patches[patch_index]
        normal_cosine = float(np.dot(directions[index], patch.normal))
        if normal_cosine >= -sense_tolerance:
            handle_rejection(
                policies,
                "outgoing",
                counts,
                f"bank {path.name}: selected record is outgoing or tangential for patch {patch.patch_id}",
            )
            continue
        group = energy_group_index(float(energies[index]), energy_edges_eV)
        if group is None:
            handle_rejection(
                policies,
                "energy_out_of_range",
                counts,
                f"bank {path.name}: event energy lies outside the requested group structure",
            )
            continue
        direction_mapping = choose_direction(
            directions[index], patch, target_directions, maximum_angle_rad
        )
        if direction_mapping is None:
            handle_rejection(
                policies,
                "angular_out_of_range",
                counts,
                f"bank {path.name}: event direction exceeds the maximum projection angle",
            )
            continue
        direction_index, angle_error = direction_mapping
        contribution = float(weights[index]) / histories
        current[patch_index, group, direction_index] += contribution
        contributions.append(contribution)
        angular_weighted_sum += contribution * angle_error
        angular_weight_sum += contribution
        maximum_angular_error = max(maximum_angular_error, angle_error)
        counts["accepted"] += 1

    if counts["accepted"] == 0:
        raise ValueError(f"bank {path.name}: no records survived the declared selection")
    return BankResult(
        batch_id=batch_id,
        path=path,
        file_sha256=actual_hash,
        source_histories=histories,
        independent=independent,
        integrated_current_per_history=current,
        accepted_contribution_weights=contributions,
        angular_error_weighted_sum=angular_weighted_sum,
        angular_error_weight_sum=angular_weight_sum,
        maximum_angular_error_rad=maximum_angular_error,
        counts=counts,
    )


def write_string(group: h5py.Group, name: str, value: str) -> None:
    group.create_dataset(name, data=np.asarray(value, dtype=UTF8))


def write_provenance(handle: h5py.File, values: dict[str, Any]) -> None:
    group = handle.create_group("provenance")
    entries = sorted((str(key), str(value)) for key, value in values.items())
    group.create_dataset("keys", data=np.asarray([key for key, _ in entries], dtype=UTF8))
    group.create_dataset("values", data=np.asarray([value for _, value in entries], dtype=UTF8))


def effective_sample_size(weights: Iterable[float]) -> float:
    values = np.asarray(list(weights), dtype=np.float64)
    if values.size == 0:
        return 0.0
    denominator = float(np.sum(values**2))
    return 0.0 if denominator == 0.0 else float(np.sum(values) ** 2 / denominator)


def build_source(specification_path: Path, output_path: Path, *, overwrite: bool) -> dict[str, Any]:
    spec = json.loads(specification_path.read_text(encoding="utf-8"))
    if spec.get("schema") != SPEC_SCHEMA:
        raise ValueError(f"unsupported OpenMC adapter schema: {spec.get('schema')!r}")
    classification = require_text(spec, "classification", context="adapter specification")
    if classification not in ALLOWED_CLASSIFICATIONS:
        raise ValueError(f"unsupported adapter classification {classification!r}")
    particle = require_text(spec, "particle", context="adapter specification").lower()
    if particle not in PARTICLE_CODES:
        raise ValueError(f"unsupported particle {particle!r}")
    source_rate_per_s = require_positive_float(
        spec, "source_rate_per_s", context="adapter specification"
    )
    symmetry_factor = require_positive_float(
        spec, "symmetry_factor", context="adapter specification"
    )
    time_class = require_text(spec, "time_class", context="adapter specification").lower()
    if time_class not in {"prompt", "delayed"}:
        raise ValueError("time_class must be prompt or delayed")
    normalization_basis = str(spec.get("normalization_basis", "per_history"))
    if normalization_basis != "per_history":
        raise ValueError("OpenMC surface-bank conversion requires normalization_basis=per_history")
    time_interval = np.asarray(spec.get("time_interval_s", [0.0, 0.0]), dtype=np.float64)
    if time_interval.shape != (2,) or not np.all(np.isfinite(time_interval)) or time_interval[1] < time_interval[0]:
        raise ValueError("time_interval_s must contain finite increasing start/stop values")
    event_window_raw = spec.get("event_time_window_s")
    event_time_window = None
    if event_window_raw is not None:
        event_window = tuple(float(value) for value in event_window_raw)
        if len(event_window) != 2 or not all(math.isfinite(value) for value in event_window) or event_window[1] < event_window[0]:
            raise ValueError("event_time_window_s must contain finite increasing values")

    tolerances = dict(spec.get("tolerances", {}))
    frame_tolerance = float(tolerances.get("frame", 1.0e-9))
    position_tolerance_cm = float(tolerances.get("position_cm", 1.0e-8))
    sense_tolerance = float(tolerances.get("sense_cosine", 1.0e-12))
    maximum_angle_rad = float(tolerances.get("maximum_projection_angle_rad", math.pi))
    current_rtol = float(tolerances.get("current_rtol", 1.0e-10))
    current_atol = float(tolerances.get("current_atol", 1.0e-14))
    if not all(math.isfinite(value) and value >= 0.0 for value in (
        frame_tolerance,
        position_tolerance_cm,
        sense_tolerance,
        maximum_angle_rad,
        current_rtol,
        current_atol,
    )):
        raise ValueError("adapter tolerances must be finite and nonnegative")
    if maximum_angle_rad > math.pi:
        raise ValueError("maximum_projection_angle_rad cannot exceed pi")

    patch_rows = spec.get("patches", [])
    if not isinstance(patch_rows, list) or not patch_rows:
        raise ValueError("adapter specification requires patches")
    patches = [parse_patch(dict(raw), frame_tolerance=frame_tolerance) for raw in patch_rows]
    patch_ids = [patch.patch_id for patch in patches]
    if len(set(patch_ids)) != len(patch_ids):
        raise ValueError("patch identifiers must be unique")

    energy_edges_eV = np.asarray(spec.get("energy_edges_eV"), dtype=np.float64)
    if energy_edges_eV.ndim != 1 or len(energy_edges_eV) < 2 or not np.all(np.isfinite(energy_edges_eV)) or np.any(energy_edges_eV < 0.0) or np.any(np.diff(energy_edges_eV) <= 0.0):
        raise ValueError("energy_edges_eV must be finite, nonnegative, and strictly increasing")
    angular = dict(spec.get("angular", {}))
    target_directions = np.asarray(angular.get("direction_cosines"), dtype=np.float64)
    quadrature_weights = np.asarray(angular.get("quadrature_weights"), dtype=np.float64)
    if target_directions.ndim != 2 or target_directions.shape[1] != 3 or len(target_directions) == 0:
        raise ValueError("angular.direction_cosines must have shape (direction,3)")
    if quadrature_weights.shape != (len(target_directions),) or not np.all(np.isfinite(quadrature_weights)) or np.any(quadrature_weights <= 0.0):
        raise ValueError("angular.quadrature_weights must be finite and positive")
    if not np.allclose(np.linalg.norm(target_directions, axis=1), 1.0, rtol=1.0e-9, atol=1.0e-10):
        raise ValueError("target direction cosines must be unit vectors")

    policies = {
        key: str(value).lower() for key, value in dict(spec.get("policies", {})).items()
    }
    for key in (
        "time_out_of_range",
        "unmapped",
        "outgoing",
        "energy_out_of_range",
        "angular_out_of_range",
    ):
        policy_action(policies, key)

    bank_specs = spec.get("bank_files", [])
    if not isinstance(bank_specs, list) or not bank_specs:
        raise ValueError("adapter specification requires bank_files")
    results = [
        read_bank(
            specification_path,
            dict(bank_spec),
            particle,
            patches,
            energy_edges_eV,
            target_directions,
            maximum_angle_rad,
            position_tolerance_cm,
            sense_tolerance,
            policies,
            event_time_window,
        )
        for bank_spec in bank_specs
    ]
    batch_ids = [result.batch_id for result in results]
    if len(set(batch_ids)) != len(batch_ids):
        raise ValueError("OpenMC bank batch_id values must be unique")

    total_histories = sum(result.source_histories for result in results)
    combined_integrated_current = sum(
        result.integrated_current_per_history * result.source_histories
        for result in results
    ) / total_histories
    areas = np.asarray([patch.area_cm2 for patch in patches], dtype=np.float64)
    directional_current_density = combined_integrated_current / areas[:, None, None]
    incoming_current = np.sum(combined_integrated_current, axis=2)

    variance = None
    variance_status = "not-estimated"
    if len(results) >= 2:
        independent = all(result.independent for result in results)
        equal_histories = len({result.source_histories for result in results}) == 1
        if independent and equal_histories:
            bank_density = np.stack(
                [result.integrated_current_per_history / areas[:, None, None] for result in results]
            )
            variance = np.var(bank_density, axis=0, ddof=1) / len(results)
            variance_status = "independent-equal-history-batch-mean"
        else:
            variance_status = (
                "not-estimated-batches-not-independent"
                if not independent
                else "not-estimated-unequal-history-batches"
            )

    expected_current = spec.get("expected_incoming_current_per_history")
    independent_current_closure = None
    independent_current_max_relative_error = None
    if expected_current is not None:
        expected = np.asarray(expected_current, dtype=np.float64)
        if expected.shape != incoming_current.shape or not np.all(np.isfinite(expected)):
            raise ValueError(
                "expected_incoming_current_per_history must match (patch,energy-group)"
            )
        difference = np.abs(incoming_current - expected)
        tolerance = current_atol + current_rtol * np.abs(expected)
        independent_current_closure = bool(np.all(difference <= tolerance))
        denominator = np.maximum(np.abs(expected), current_atol)
        independent_current_max_relative_error = float(
            np.max(difference / denominator)
        )
        if not independent_current_closure:
            raise ValueError(
                "converted OpenMC bank does not close the independent incoming-current ledger"
            )

    bank_hashes = [result.file_sha256 for result in results]
    specification_sha256 = sha256_file(specification_path)
    source_hash = canonical_sha256(
        {
            "schema": SPEC_SCHEMA,
            "specification_sha256": specification_sha256,
            "bank_sha256": bank_hashes,
            "particle": particle,
            "patch_ids": patch_ids,
            "energy_edges_eV": energy_edges_eV.tolist(),
            "direction_cosines": target_directions.tolist(),
            "quadrature_weights": quadrature_weights.tolist(),
            "normalization_basis": normalization_basis,
        }
    )

    if output_path.exists() and not overwrite:
        raise FileExistsError(f"refusing to overwrite {output_path}")
    if output_path.exists() and not output_path.is_file():
        raise ValueError(f"output path exists and is not a file: {output_path}")
    output_path.parent.mkdir(parents=True, exist_ok=True)
    if output_path.exists():
        output_path.unlink()

    with h5py.File(output_path, "w", libver="earliest") as handle:
        meta = handle.create_group("meta")
        patch_group = handle.create_group("patch")
        energy_group = handle.create_group("energy")
        angle_group = handle.create_group("angle")
        source_group = handle.create_group("source")
        write_string(meta, "schema", OUTPUT_SCHEMA)
        write_string(meta, "particle", particle)
        write_string(meta, "source_representation", "directional_current_density")
        write_string(meta, "storage_order", STORAGE_ORDER)
        write_string(meta, "normalization_basis", normalization_basis)
        meta.create_dataset("source_rate_per_s", data=np.float64(source_rate_per_s))
        meta.create_dataset("symmetry_factor", data=np.float64(symmetry_factor))
        meta.create_dataset("time_interval_s", data=time_interval)
        write_string(meta, "time_class", time_class)
        write_string(meta, "source_hash", source_hash)

        patch_group.create_dataset("id", data=np.asarray(patch_ids, dtype=np.int64))
        patch_group.create_dataset(
            "centroid_cm", data=np.stack([patch.centroid_cm for patch in patches])
        )
        patch_group.create_dataset("area_cm2", data=areas)
        patch_group.create_dataset(
            "normal", data=np.stack([patch.normal for patch in patches])
        )
        patch_group.create_dataset(
            "tangent_1", data=np.stack([patch.tangent_1 for patch in patches])
        )
        patch_group.create_dataset(
            "tangent_2", data=np.stack([patch.tangent_2 for patch in patches])
        )
        energy_group.create_dataset("edges_eV", data=energy_edges_eV)
        angle_group.create_dataset("direction_cosines", data=target_directions)
        angle_group.create_dataset("quadrature_weights", data=quadrature_weights)
        source_group.create_dataset(
            "directional_current_density", data=directional_current_density
        )
        source_group.create_dataset("incoming_current", data=incoming_current)
        if variance is not None:
            source_group.create_dataset("variance", data=variance)
        write_provenance(
            handle,
            {
                "adapter": "scripts/adapters/openmc_surface_bank_to_radiant.py",
                "adapter_schema": SPEC_SCHEMA,
                "classification": classification,
                "specification_sha256": specification_sha256,
                "bank_sha256": ",".join(bank_hashes),
                "bank_batch_ids": ",".join(batch_ids),
                "source_histories": total_histories,
                "crossing_weight_semantics": "integrated-current-per-originating-history",
                "angular_projection": "nearest-incoming-target-ordinate-current-preserving",
                "variance_status": variance_status,
                "physical_replay_pass": "false",
                "independent_current_ledger_present": str(
                    expected_current is not None
                ).lower(),
            },
        )

    accepted_contributions = [
        contribution
        for result in results
        for contribution in result.accepted_contribution_weights
    ]
    angular_weight = sum(result.angular_error_weight_sum for result in results)
    weighted_angle = (
        sum(result.angular_error_weighted_sum for result in results) / angular_weight
        if angular_weight > 0.0
        else 0.0
    )
    counts_total = {
        key: sum(result.counts.get(key, 0) for result in results)
        for key in results[0].counts
    }
    output_sha256 = sha256_file(output_path)
    receipt = {
        "schema": "radiant.openmc_surface_bank_adapter_receipt/v1",
        "classification": classification,
        "specification_sha256": specification_sha256,
        "output_path": str(output_path.resolve()),
        "output_sha256": output_sha256,
        "source_hash": source_hash,
        "particle": particle,
        "patch_ids": patch_ids,
        "energy_edges_eV": energy_edges_eV.tolist(),
        "direction_count": len(target_directions),
        "bank_count": len(results),
        "total_source_histories": total_histories,
        "variance_status": variance_status,
        "variance_present": variance is not None,
        "counts": counts_total,
        "effective_sample_size": effective_sample_size(accepted_contributions),
        "weighted_mean_projection_angle_rad": weighted_angle,
        "maximum_projection_angle_rad": max(
            result.maximum_angular_error_rad for result in results
        ),
        "internal_current_closure_pass": True,
        "independent_current_ledger_present": expected_current is not None,
        "independent_current_closure_pass": independent_current_closure,
        "independent_current_max_relative_error": independent_current_max_relative_error,
        "physical_replay_pass": False,
        "banks": [
            {
                "batch_id": result.batch_id,
                "basename": result.path.name,
                "sha256": result.file_sha256,
                "source_histories": result.source_histories,
                "independent": result.independent,
                "counts": result.counts,
                "effective_sample_size": effective_sample_size(
                    result.accepted_contribution_weights
                ),
                "weighted_mean_projection_angle_rad": (
                    result.angular_error_weighted_sum / result.angular_error_weight_sum
                    if result.angular_error_weight_sum > 0.0
                    else 0.0
                ),
                "maximum_projection_angle_rad": result.maximum_angular_error_rad,
            }
            for result in results
        ],
    }
    return receipt


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("specification", type=Path)
    parser.add_argument("output", type=Path)
    parser.add_argument("--receipt", type=Path)
    parser.add_argument("--overwrite", action="store_true")
    args = parser.parse_args()

    receipt = build_source(
        args.specification.resolve(), args.output.resolve(), overwrite=args.overwrite
    )
    if args.receipt is not None:
        args.receipt.parent.mkdir(parents=True, exist_ok=True)
        args.receipt.write_text(
            json.dumps(receipt, indent=2, sort_keys=True) + "\n", encoding="utf-8"
        )
    print(json.dumps(receipt, sort_keys=True))


if __name__ == "__main__":
    main()
