#!/usr/bin/env python3
"""Download and checksum every public file in a Zenodo record using only the Python stdlib.

Large cascade repositories stay outside Git. The generated JSON manifest binds record metadata,
file names, sizes, and checksums for later Radiant qualification. Existing complete files are
reused only when their advertised checksum matches.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import shutil
import urllib.request


def _digest(path: Path, algorithm: str = "sha256") -> str:
    hasher = hashlib.new(algorithm)
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(8 * 1024 * 1024), b""):
            hasher.update(chunk)
    return hasher.hexdigest()


def _open_json(url: str) -> dict:
    request = urllib.request.Request(url, headers={"User-Agent": "RadiantHTS-data-adapter/1"})
    with urllib.request.urlopen(request, timeout=120) as response:
        return json.load(response)


def _download(url: str, destination: Path) -> None:
    partial = destination.with_suffix(destination.suffix + ".part")
    partial.parent.mkdir(parents=True, exist_ok=True)
    request = urllib.request.Request(url, headers={"User-Agent": "RadiantHTS-data-adapter/1"})
    with urllib.request.urlopen(request, timeout=120) as response, partial.open("wb") as output:
        shutil.copyfileobj(response, output, length=8 * 1024 * 1024)
    partial.replace(destination)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("record_id", type=int)
    parser.add_argument("output_directory", type=Path)
    parser.add_argument("--metadata-only", action="store_true")
    args = parser.parse_args()
    if args.record_id <= 0:
        parser.error("record_id must be positive")

    metadata_url = f"https://zenodo.org/api/records/{args.record_id}"
    metadata = _open_json(metadata_url)
    args.output_directory.mkdir(parents=True, exist_ok=True)
    records = []
    for file_entry in metadata.get("files", []):
        key = file_entry.get("key") or file_entry.get("filename")
        if not key:
            raise RuntimeError("Zenodo file entry is missing its key/name")
        links = file_entry.get("links", {})
        url = links.get("content") or links.get("self")
        if not url:
            raise RuntimeError(f"Zenodo file {key} has no content link")
        destination = args.output_directory / key
        advertised_checksum = str(file_entry.get("checksum", ""))
        expected_algorithm = "md5"
        expected_digest = advertised_checksum
        if ":" in advertised_checksum:
            expected_algorithm, expected_digest = advertised_checksum.split(":", 1)
        status = "metadata-only"
        if not args.metadata_only:
            reuse = destination.is_file()
            if reuse and expected_digest:
                reuse = _digest(destination, expected_algorithm) == expected_digest.lower()
            if not reuse:
                _download(url, destination)
            if expected_digest:
                calculated = _digest(destination, expected_algorithm)
                if calculated != expected_digest.lower():
                    destination.unlink(missing_ok=True)
                    raise RuntimeError(f"Checksum mismatch for {key}")
            status = "downloaded-and-verified"
        record = {
            "name": key,
            "size": int(file_entry.get("size", 0)),
            "advertised_checksum": advertised_checksum,
            "content_url": url,
            "status": status,
        }
        if destination.is_file():
            record["sha256"] = _digest(destination, "sha256")
            record["local_path"] = str(destination.resolve())
        records.append(record)

    manifest = {
        "schema": "radiant.external_data_download/v1",
        "record_id": args.record_id,
        "record_doi": metadata.get("doi"),
        "record_title": metadata.get("metadata", {}).get("title"),
        "record_updated": metadata.get("updated"),
        "record_url": metadata.get("links", {}).get("html"),
        "metadata_url": metadata_url,
        "files": records,
        "large_data_outside_git": True,
    }
    manifest_path = args.output_directory / "zenodo-download-manifest.json"
    manifest_path.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")
    print(manifest_path)


if __name__ == "__main__":
    main()
