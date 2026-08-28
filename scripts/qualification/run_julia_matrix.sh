#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
JULIA16_BIN="${JULIA16_BIN:-julia-1.6.7}"
JULIA110_BIN="${JULIA110_BIN:-julia-1.10.12}"

run_one() {
    local executable="$1"
    local expected="$2"
    if ! command -v "$executable" >/dev/null 2>&1 && [[ ! -x "$executable" ]]; then
        echo "Missing Julia ${expected} executable: ${executable}" >&2
        return 2
    fi
    local actual
    actual="$($executable --startup-file=no -e 'print(VERSION)')"
    if [[ "$actual" != "$expected" ]]; then
        echo "Expected Julia ${expected}, found ${actual} at ${executable}" >&2
        return 2
    fi
    echo "Running Radiant qualification with Julia ${actual}"
    "$executable" --startup-file=no --project="$ROOT" \
        "$ROOT/scripts/qualification/run_julia_qualification.jl"
}

run_one "$JULIA16_BIN" "1.6.7"
run_one "$JULIA110_BIN" "1.10.12"

echo "Pinned Julia qualification matrix passed."
