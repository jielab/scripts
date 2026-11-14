#!/usr/bin/env bash
set -euo pipefail
args=()
if [[ -n ${LD_ZSTD_PREFIX:-} ]]; then
    args+=("-I$LD_ZSTD_PREFIX/include" "-L$LD_ZSTD_PREFIX/lib" "-Wl,-rpath,$LD_ZSTD_PREFIX/lib")
fi
exec "${LD_CXX:-g++}" "${args[@]}" "$@" -pthread -fopenmp
