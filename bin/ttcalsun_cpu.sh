#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

source "$ROOT/install.sh"

exec julia --project="$ROOT" "$ROOT/bin/ttcalsun.jl" "$@"
