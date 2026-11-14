#!/usr/bin/env bash
# Probe executable SIMD variants: some virtualized CPUs advertise AVX512 while
# their OS cannot run the bundled AVX512 build. Every variant is BWA-MEM2 2.2.1.
set -euo pipefail
HERE="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
for variant in avx512bw avx2 avx sse42 sse41; do
    exe="$HERE/bwa-mem2.$variant"
    if [[ -x "$exe" ]] && "$exe" version >/dev/null 2>&1; then
        exec "$exe" "$@"
    fi
done
echo 'No compatible BWA-MEM2 SIMD executable on this Linux x86_64 host' >&2
exit 2
