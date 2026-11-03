#!/usr/bin/env bash
set -euo pipefail
prefix="${LD_ZSTD_PREFIX:-${CONDA_PREFIX:?Activate the LD environment first}}"
args=("-I$prefix/include")
for arg in "$@"; do
    if [[ "$arg" == -lzstd ]]; then
        args+=("-L$prefix/lib" "-Wl,-rpath,$prefix/lib" -lzstd -pthread)
    else
        args+=("$arg")
    fi
done
exec "${LD_CXX:-${CONDA_PREFIX}/bin/x86_64-conda-linux-gnu-c++}" "${args[@]}"
