#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
REPO_ROOT="$(cd "$ROOT/.." && pwd)"

export PYTHON="${PYTHON:-$REPO_ROOT/lwa/bin/python}"
export JULIA_DEPOT_PATH="${JULIA_DEPOT_PATH:-/tmp/julia_ttcalsun:$HOME/.julia}"
export JULIA_NUM_THREADS="${JULIA_NUM_THREADS:-64}"
export OPENBLAS_NUM_THREADS="${OPENBLAS_NUM_THREADS:-64}"
export MKL_NUM_THREADS="${MKL_NUM_THREADS:-64}"
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-64}"

exec julia --project="$ROOT" "$ROOT/bin/ttcalsun.jl" "$@"
