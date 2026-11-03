#!/usr/bin/env bash
set -euo pipefail
prefix="${LD_ZSTD_PREFIX:-/share/apps/anaconda3/2020.7}"
args=("-I$prefix/include")
for arg in "$@"; do
    if [[ "$arg" == -lzstd ]]; then
        args+=("$prefix/lib/libzstd.a" -pthread)
    else
        args+=("$arg")
    fi
done
exec /share/base/gcc/11.2.0/bin/g++ "${args[@]}"
