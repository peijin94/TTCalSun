#!/usr/bin/env bash

# Execute this file to initialize TTCalSun:
#   ./install.sh
#
# Source this file to configure the runtime environment without installing:
#   source ./install.sh
#
# Both modes keep Julia package state, precompile caches, scratch data, and
# temp files under /fast/rtpipe/env/juliadepot by default.

ttcalsun_configure() {
    local ttcalsun_root
    local rtpipe_root
    ttcalsun_root="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
    rtpipe_root="$(cd "${ttcalsun_root}/.." && pwd)"

    export TTCALSUN_ROOT="${TTCALSUN_ROOT:-${ttcalsun_root}}"
    export RTPIPE_ROOT="${RTPIPE_ROOT:-${rtpipe_root}}"

    if [[ -f "${RTPIPE_ROOT}/use_julia.sh" ]]; then
        source "${RTPIPE_ROOT}/use_julia.sh"
    else
        export JULIA_DEPOT_ROOT="${JULIA_DEPOT_ROOT:-${RTPIPE_ROOT}/env/juliadepot}"
        export JULIA_DEPOT_PATH="${JULIA_DEPOT_PATH:-${JULIA_DEPOT_ROOT}}"
        export JULIA_PKG_DEVDIR="${JULIA_PKG_DEVDIR:-${JULIA_DEPOT_ROOT}/dev}"
        export JULIA_HISTORY="${JULIA_HISTORY:-${JULIA_DEPOT_ROOT}/logs/repl_history.jl}"
        export TMPDIR="${JULIA_TMPDIR:-${JULIA_DEPOT_ROOT}/tmp}"
        mkdir -p \
            "${JULIA_DEPOT_ROOT}" \
            "${JULIA_DEPOT_ROOT}/dev" \
            "${JULIA_DEPOT_ROOT}/logs" \
            "${TMPDIR}"
    fi

    if [[ -x "${RTPIPE_ROOT}/env/lwa/bin/python" ]]; then
        export PYTHON="${PYTHON:-${RTPIPE_ROOT}/env/lwa/bin/python}"
    else
        export PYTHON="${PYTHON:-${RTPIPE_ROOT}/lwa/bin/python}"
    fi

    export JULIA_NUM_THREADS="${JULIA_NUM_THREADS:-64}"
    export OPENBLAS_NUM_THREADS="${OPENBLAS_NUM_THREADS:-64}"
    export MKL_NUM_THREADS="${MKL_NUM_THREADS:-64}"
    export OMP_NUM_THREADS="${OMP_NUM_THREADS:-64}"
}

ttcalsun_install() {
    ttcalsun_configure
    echo "Installing TTCalSun Julia dependencies into ${JULIA_DEPOT_PATH}"
    julia --project="${TTCALSUN_ROOT}" -e 'using Pkg; Pkg.instantiate(); Pkg.precompile()'
}

if [[ "${BASH_SOURCE[0]}" == "$0" ]]; then
    set -euo pipefail
    ttcalsun_install "$@"
else
    ttcalsun_configure "$@"
fi
