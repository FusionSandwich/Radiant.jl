#!/usr/bin/env python3
"""Create synthetic arrays and a v2 builder specification for cross-language tests."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

SPEC_SCHEMA = "radiant.hts.physical_reference_builder_spec/v2"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("output_directory", type=Path)
    args = parser.parse_args()
    output = args.output_directory.resolve()
    output.mkdir(parents=True, exist_ok=True)

    photoelectric_values = np.asarray([[1.0, 2.0], [3.0, 4.0]], dtype=np.float64)
    photoelectric_uncertainty = np.asarray(
        [[0.01, 0.02], [0.03, 0.04]], dtype=np.float64
    )
    compton_values = np.asarray([[0.5, 0.25], [0.125, 0.0625]], dtype=np.float64)
    compton_uncertainty = np.full((2, 2), 0.005, dtype=np.float64)

    np.save(output / "photoelectric-values.npy", photoelectric_values, allow_pickle=False)
    np.save(
        output / "photoelectric-uncertainty.npy",
        photoelectric_uncertainty,
        allow_pickle=False,
    )
    np.savetxt(output / "compton-values.csv", compton_values, delimiter=",")
    np.savetxt(
        output / "compton-uncertainty.csv", compton_uncertainty, delimiter="," 
    )

    defaults = {
        "classification": "synthetic",
        "particle_tag": "photon",
        "units": "MeV/g per transport source basis",
        "producer": "python-cross-language-fixture",
        "source_artifact_hash": "1" * 64,
        "geometry_hash": "2" * 64,
        "material_state_hash": "3" * 64,
        "normalization_hash": "4" * 64,
    }
    specification = {
        "schema": SPEC_SCHEMA,
        "require_unique_process_keys": True,
        "defaults": defaults,
        "metadata": {
            "classification": "software-verification",
            "synthetic_fixture_may_set_physical_pass": False,
            "expected_shape": "2,2",
        },
        "responses": [
            {
                "response_id": "photoelectric-layer-field",
                "process_key": "photoelectric",
                "scoring_semantics": "synthetic layer-by-sublayer deposition field",
                "values_path": "photoelectric-values.npy",
                "uncertainty_path": "photoelectric-uncertainty.npy",
                "result_artifact_hash": "5" * 64,
                "metadata": {"array_format": "npy"},
            },
            {
                "response_id": "compton-layer-field",
                "process_key": "compton",
                "scoring_semantics": "synthetic layer-by-sublayer deposition field",
                "values_path": "compton-values.csv",
                "uncertainty_path": "compton-uncertainty.csv",
                "result_artifact_hash": "6" * 64,
                "metadata": {"array_format": "csv"},
            },
        ],
    }
    spec_path = output / "physical-reference-spec.json"
    spec_path.write_text(
        json.dumps(specification, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(spec_path)


if __name__ == "__main__":
    main()
