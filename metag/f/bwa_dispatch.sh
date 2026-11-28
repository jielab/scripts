#!/usr/bin/env bash

# Standalone help must return before runtime setup or dispatch.
if [[ ${BASH_SOURCE[0]} == "$0" && ( ${1:-} == -h || ${1:-} == --help ) ]]; then
    cat <<'QUICKSTART_HELP'
bash ./bwa_dispatch.sh -h  # 内部工具/单独阶段；调用要求见下文
# 从所属项目主入口选择模块；本文件不代表整个项目
# 必填输入沿用现有接口，不自动补造数据或调度参数

This is an internal helper. Use the owning project entry for analysis and scheduling.
QUICKSTART_HELP
    exit 0
fi
# Probe executable SIMD variants: some virtualized CPUs advertise AVX512 while
# their OS cannot run the bundled AVX512 build. Every variant is BWA-MEM2 2.2.1.
set -euo pipefail
HERE="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)"
for variant in avx512bw avx2 avx sse42 sse41; do
    exe="$HERE/bwa-mem2.$variant"
    if [[ -x "$exe" ]] && "$exe" version >/dev/null 2>&1; then
        exec "$exe" "$@"
    fi
done
echo 'No compatible BWA-MEM2 SIMD executable on this Linux x86_64 host' >&2
exit 2

