#!/usr/bin/env python3
"""Fail-closed preflight for the first physical OpenMC/OpenSn/Radiant YBCO replay.

The preflight validates artifacts, schemas, checksums, lineage declarations, ownership, and
software/runtime pins. It never runs transport and never sets a physical-science pass. A successful
result means only that the requested physical run is ready to start.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import re
import tempfile
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Iterable

try:
    import h5py  # type: ignore
except ImportError:  # pragma: no cover - handled explicitly at runtime
    h5py = None

SCHEMA = "radiant.hts.first_physical_replay_manifest/v1"
RECEIPT_SCHEMA = "radiant.hts.first_physical_replay_preflight/v1"
SHA256 = re.compile(r"^[0-9a-f]{64}$")
FORBIDDEN_PHYSICAL_CLASSIFICATIONS = {
    "analytic",
    "synthetic",
    "software",
    "software-verification",
    "verification-only",
    "representative",
    "unbound",
}

REQUIRED_ARTIFACTS = {
    "openmc_boundary_source": "radiant.boundary_angular_current/v1",
    "openmc_reference_bundle": "radiant.hts.physical_reference_bundle/v1",
    "local_geometry": "radiant.faceted_surface_mesh/v1",
    "neutron_multigroup_data": "radiant.hts.multigroup_neutron/v1",
    "photon_multigroup_data": "radiant.hts.multigroup_photon/v1",
    "neutron_to_photon_data": "radiant.hts.neutron_to_photon/v1",
    "material_response_bundle": "radiant.hts.material_response_bundle/v1",
    "transport_ownership_map": "parastell_damage.transport_ownership_map/v1",
}

OPTIONAL_ARTIFACTS = {
    "activation_inventory": "radiant.hts.activation_inventory/v1",
    "delayed_source_bundle": "radiant.hts.activation_delayed_source_bundle/v1",
    "gdbco_cascade_data": "radiant.hts.evaluated_gd_cascade/v1",
    "curved_reference_bundle": "radiant.hts.physical_reference_bundle/v1",
}

REQUIRED_RUNTIME_PINS = {
    "radiant_git_sha",
    "opensn_git_sha",
    "openmc_version",
    "julia_version",
    "python_version",
}

REQUIRED_PROTECTED_RESPONSES = {
    "layer_scalar_flux",
    "layer_total_deposition",
    "layer_prompt_lattice_heat",
    "peak_sublayer_deposition",
    "interface_incoming_current",
    "interface_outgoing_current",
    "particle_balance",
    "energy_balance",
}


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _read_hdf5_scalar(value: Any) -> str:
    if hasattr(value, "shape") and value.shape == ():
        value = value[()]
    if isinstance(value, bytes):
        return value.decode("utf-8")
    if hasattr(value, "item"):
        value = value.item()
        if isinstance(value, bytes):
            return value.decode("utf-8")
    return str(value)


def artifact_schema(path: Path) -> str:
    suffix = path.suffix.lower()
    if suffix in {".h5", ".hdf5", ".h5m"}:
        if h5py is None:
            raise RuntimeError("h5py is required to inspect HDF5 replay artifacts")
        with h5py.File(path, "r") as handle:
            candidates = ("meta/schema", "schema")
            for candidate in candidates:
                if candidate in handle:
                    return _read_hdf5_scalar(handle[candidate])
        raise ValueError(f"HDF5 artifact has no schema dataset: {path}")
    if suffix == ".json":
        document = json.loads(path.read_text(encoding="utf-8"))
        return str(document.get("schema", ""))
    if suffix == ".toml":
        try:
            import tomllib
        except ImportError as exc:  # pragma: no cover
            raise RuntimeError("Python 3.11+ is required for TOML preflight artifacts") from exc
        document = tomllib.loads(path.read_text(encoding="utf-8"))
        return str(document.get("schema", ""))
    raise ValueError(f"Unsupported replay artifact format: {path}")


def _inside_repository(path: Path, repository_root: Path | None) -> bool:
    if repository_root is None:
        return False
    try:
        path.resolve().relative_to(repository_root.resolve())
        return True
    except ValueError:
        return False


def _classification_is_physical(value: str) -> bool:
    normalized = value.strip().lower().replace("_", "-")
    return bool(normalized) and normalized not in FORBIDDEN_PHYSICAL_CLASSIFICATIONS


@dataclass(frozen=True)
class ArtifactCheck:
    artifact_id: str
    path: str
    expected_schema: str
    observed_schema: str | None
    expected_sha256: str
    observed_sha256: str | None
    classification: str
    exists: bool
    passed: bool
    reasons: tuple[str, ...]

    def as_mapping(self) -> dict[str, Any]:
        return {
            "artifact_id": self.artifact_id,
            "path": self.path,
            "expected_schema": self.expected_schema,
            "observed_schema": self.observed_schema,
            "expected_sha256": self.expected_sha256,
            "observed_sha256": self.observed_sha256,
            "classification": self.classification,
            "exists": self.exists,
            "passed": self.passed,
            "reasons": list(self.reasons),
        }


def check_artifact(
    artifact_id: str,
    declaration: dict[str, Any],
    expected_schema: str,
    *,
    repository_root: Path | None,
    require_external_large_artifacts: bool,
) -> ArtifactCheck:
    reasons: list[str] = []
    raw_path = str(declaration.get("path", ""))
    path = Path(raw_path).expanduser() if raw_path else Path(".")
    expected_hash = str(declaration.get("sha256", "")).lower()
    classification = str(declaration.get("classification", ""))
    exists = bool(raw_path) and path.is_file()
    observed_hash: str | None = None
    observed_schema: str | None = None

    if not raw_path:
        reasons.append("path is missing")
    elif not exists:
        reasons.append("artifact file does not exist")
    if not SHA256.fullmatch(expected_hash):
        reasons.append("expected SHA-256 is missing or malformed")
    if not _classification_is_physical(classification):
        reasons.append("artifact classification is not physical")
    if require_external_large_artifacts and raw_path and _inside_repository(path, repository_root):
        reasons.append("physical/large artifact is stored inside the source repository")

    if exists:
        observed_hash = sha256_file(path)
        if SHA256.fullmatch(expected_hash) and observed_hash != expected_hash:
            reasons.append("artifact SHA-256 does not match")
        try:
            observed_schema = artifact_schema(path)
        except Exception as exc:  # noqa: BLE001 - recorded in fail-closed receipt
            reasons.append(f"schema read failed: {type(exc).__name__}: {exc}")
        else:
            declared_schema = str(declaration.get("schema", expected_schema))
            if declared_schema != expected_schema:
                reasons.append(
                    f"manifest schema declaration {declared_schema!r} does not equal required "
                    f"schema {expected_schema!r}"
                )
            if observed_schema != expected_schema:
                reasons.append(
                    f"observed artifact schema {observed_schema!r} does not equal required "
                    f"schema {expected_schema!r}"
                )

    return ArtifactCheck(
        artifact_id=artifact_id,
        path=str(path.resolve()) if raw_path else "",
        expected_schema=expected_schema,
        observed_schema=observed_schema,
        expected_sha256=expected_hash,
        observed_sha256=observed_hash,
        classification=classification,
        exists=exists,
        passed=not reasons,
        reasons=tuple(reasons),
    )


def validate_manifest(
    manifest: dict[str, Any],
    *,
    manifest_path: Path,
    repository_root: Path | None = None,
) -> dict[str, Any]:
    blockers: list[str] = []
    warnings: list[str] = []
    if manifest.get("schema") != SCHEMA:
        blockers.append(f"manifest schema must be {SCHEMA}")
    case_id = str(manifest.get("case_id", ""))
    if not case_id:
        blockers.append("case_id is missing")
    if str(manifest.get("material_family", "")) not in {"YBCO", "GdBCO"}:
        blockers.append("material_family must be YBCO or GdBCO")
    if not _classification_is_physical(str(manifest.get("classification", ""))):
        blockers.append("manifest classification is not physical")

    runtime = manifest.get("runtime_pins", {})
    missing_runtime = sorted(REQUIRED_RUNTIME_PINS - set(runtime))
    if missing_runtime:
        blockers.append(f"missing runtime pins: {missing_runtime}")
    for key in ("radiant_git_sha", "opensn_git_sha"):
        value = str(runtime.get(key, ""))
        if value and not re.fullmatch(r"[0-9a-f]{40}", value):
            blockers.append(f"runtime pin {key} is not a 40-character Git SHA")

    responses = set(map(str, manifest.get("protected_responses", [])))
    missing_responses = sorted(REQUIRED_PROTECTED_RESPONSES - responses)
    if missing_responses:
        blockers.append(f"missing protected responses: {missing_responses}")

    ownership = manifest.get("ownership", {})
    if ownership.get("one_additive_owner_per_population_response") is not True:
        blockers.append("ownership does not require one additive owner per population/response")
    if ownership.get("comparison_responses_excluded_from_production_sum") is not True:
        blockers.append("comparison responses are not excluded from production heating")
    if ownership.get("prompt_delayed_separation") is not True:
        blockers.append("prompt and delayed source populations are not separated")

    convergence_axes = set(map(str, manifest.get("required_convergence_axes", [])))
    required_axes = {"source_bank", "spatial", "angular", "energy"}
    missing_axes = sorted(required_axes - convergence_axes)
    if missing_axes:
        blockers.append(f"missing first-case convergence axes: {missing_axes}")

    artifacts = manifest.get("artifacts", {})
    checks: list[ArtifactCheck] = []
    for artifact_id, schema in REQUIRED_ARTIFACTS.items():
        declaration = artifacts.get(artifact_id)
        if not isinstance(declaration, dict):
            blockers.append(f"required artifact declaration is missing: {artifact_id}")
            continue
        check = check_artifact(
            artifact_id,
            declaration,
            schema,
            repository_root=repository_root,
            require_external_large_artifacts=True,
        )
        checks.append(check)
        blockers.extend(f"{artifact_id}: {reason}" for reason in check.reasons)

    for artifact_id, schema in OPTIONAL_ARTIFACTS.items():
        declaration = artifacts.get(artifact_id)
        if declaration is None:
            continue
        if not isinstance(declaration, dict):
            blockers.append(f"optional artifact declaration is malformed: {artifact_id}")
            continue
        check = check_artifact(
            artifact_id,
            declaration,
            schema,
            repository_root=repository_root,
            require_external_large_artifacts=True,
        )
        checks.append(check)
        blockers.extend(f"{artifact_id}: {reason}" for reason in check.reasons)

    material_family = str(manifest.get("material_family", ""))
    if material_family == "GdBCO" and "gdbco_cascade_data" not in artifacts:
        blockers.append("GdBCO replay requires complete evaluated Gd cascade data")
    if material_family == "YBCO" and "gdbco_cascade_data" in artifacts:
        warnings.append("Gd cascade data were supplied for a YBCO replay and will not be used")

    normalization = manifest.get("normalization", {})
    required_normalization = {
        "source_rate_per_s",
        "symmetry_factor",
        "source_basis",
        "normalization_hash",
    }
    missing_normalization = sorted(required_normalization - set(normalization))
    if missing_normalization:
        blockers.append(f"missing normalization fields: {missing_normalization}")
    if float(normalization.get("source_rate_per_s", 0.0) or 0.0) <= 0.0:
        blockers.append("source_rate_per_s must be positive")
    if float(normalization.get("symmetry_factor", 0.0) or 0.0) <= 0.0:
        blockers.append("symmetry_factor must be positive")

    source_hashes = {
        str(item.get("source_hash", ""))
        for item in artifacts.values()
        if isinstance(item, dict) and item.get("source_hash")
    }
    if len(source_hashes) > 1:
        blockers.append("artifact declarations contain inconsistent source_hash values")
    geometry_hashes = {
        str(item.get("geometry_hash", ""))
        for item in artifacts.values()
        if isinstance(item, dict) and item.get("geometry_hash")
    }
    if len(geometry_hashes) > 1:
        blockers.append("artifact declarations contain inconsistent geometry_hash values")

    return {
        "schema": RECEIPT_SCHEMA,
        "case_id": case_id,
        "manifest_path": str(manifest_path.resolve()),
        "manifest_sha256": sha256_file(manifest_path),
        "status": "READY_TO_EXECUTE_PHYSICAL_REPLAY" if not blockers else "BLOCKED",
        "ready_to_execute": not blockers,
        "physical_science_pass": False,
        "physical_science_pass_reason": (
            "Preflight validates readiness only; transport and matched-response evidence have not "
            "been executed by this command."
        ),
        "artifact_checks": [check.as_mapping() for check in checks],
        "blockers": blockers,
        "warnings": warnings,
        "required_artifacts": REQUIRED_ARTIFACTS,
        "optional_artifacts": OPTIONAL_ARTIFACTS,
        "required_protected_responses": sorted(REQUIRED_PROTECTED_RESPONSES),
        "required_first_case_convergence_axes": sorted(required_axes),
        "no_geant4_requirement": True,
    }


def template_manifest() -> dict[str, Any]:
    artifacts = {
        key: {
            "path": f"ABSOLUTE_EXTERNAL_PATH/{key}.h5",
            "sha256": "0" * 64,
            "schema": schema,
            "classification": "physical-candidate",
            "source_hash": "BOUND_SOURCE_HASH",
            "geometry_hash": "BOUND_GEOMETRY_HASH",
        }
        for key, schema in REQUIRED_ARTIFACTS.items()
    }
    return {
        "schema": SCHEMA,
        "case_id": "straight-ybco-first-physical-replay",
        "material_family": "YBCO",
        "classification": "physical-candidate",
        "runtime_pins": {
            "radiant_git_sha": "0" * 40,
            "opensn_git_sha": "0" * 40,
            "openmc_version": "0.16.0",
            "julia_version": "1.10.12",
            "python_version": "3.12",
        },
        "normalization": {
            "source_rate_per_s": 1.0,
            "symmetry_factor": 1.0,
            "source_basis": "per_source_history",
            "normalization_hash": "BOUND_NORMALIZATION_HASH",
        },
        "ownership": {
            "one_additive_owner_per_population_response": True,
            "comparison_responses_excluded_from_production_sum": True,
            "prompt_delayed_separation": True,
        },
        "protected_responses": sorted(REQUIRED_PROTECTED_RESPONSES),
        "required_convergence_axes": ["source_bank", "spatial", "angular", "energy"],
        "artifacts": artifacts,
    }


def _self_test() -> None:
    with tempfile.TemporaryDirectory() as directory:
        root = Path(directory)
        manifest_path = root / "manifest.json"
        manifest = template_manifest()
        manifest_path.write_text(json.dumps(manifest, indent=2), encoding="utf-8")
        result = validate_manifest(manifest, manifest_path=manifest_path)
        assert result["status"] == "BLOCKED"
        assert not result["ready_to_execute"]
        assert not result["physical_science_pass"]
        assert any("does not exist" in blocker for blocker in result["blockers"])
        assert result["no_geant4_requirement"]


def main(argv: Iterable[str] | None = None) -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("manifest", nargs="?", type=Path)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--repository-root", type=Path)
    parser.add_argument("--write-template", type=Path)
    parser.add_argument("--self-test", action="store_true")
    args = parser.parse_args(list(argv) if argv is not None else None)

    if args.self_test:
        _self_test()
        print("FIRST_YBCO_PHYSICAL_REPLAY_PREFLIGHT_SELF_TEST_PASS")
        return 0
    if args.write_template is not None:
        args.write_template.parent.mkdir(parents=True, exist_ok=True)
        args.write_template.write_text(
            json.dumps(template_manifest(), indent=2, sort_keys=True) + "\n",
            encoding="utf-8",
        )
        print(args.write_template)
        return 0
    if args.manifest is None:
        parser.error("manifest is required unless --self-test or --write-template is used")
    if not args.manifest.is_file():
        parser.error(f"manifest does not exist: {args.manifest}")

    manifest = json.loads(args.manifest.read_text(encoding="utf-8"))
    receipt = validate_manifest(
        manifest,
        manifest_path=args.manifest,
        repository_root=args.repository_root,
    )
    output = args.output or args.manifest.with_name(
        args.manifest.stem + ".preflight-receipt.json"
    )
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(receipt, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(output)
    return 0 if receipt["ready_to_execute"] else 2


if __name__ == "__main__":
    raise SystemExit(main())
