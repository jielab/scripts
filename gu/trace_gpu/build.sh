#!/usr/bin/env bash


# Direct CLI help must return before setup, downloads or analysis.
if [[ "${BASH_SOURCE[0]}" == "$0" ]]; then
	case "${1:-}" in
		-h | --help)
			cat <<'HELP'
source ./build.sh  # 内部函数库，由所属入口加载
# 无独立分析模块；请通过所属项目的主 .sh 选择模块
# 参数和环境变量沿用调用方；不要作为全流程提交

Modules:
    Default  Build and check the TRACE GPU library

Main options:
    TRACE_GPU_HOME  Installation directory
    GU_NVCC / GU_PYTHON  Compiler and runtime
    -h, --help  Show help

Example:
    ./build.sh
HELP
			exit 0
			;;
	esac
fi
set -euo pipefail
here=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
out=${TRACE_GPU_HOME:-${GU_SOFT:-$HOME/.local/share/gu}/trace-gpu}
nvcc=${GU_NVCC:-$(command -v nvcc || true)}
python=${GU_PYTHON:-$(command -v python3 || true)}
[[ -n $nvcc && -x $nvcc && -n $python && -x $python ]] || { echo "ERROR: CUDA nvcc and Python are required; set GU_NVCC and GU_PYTHON" >&2; exit 2; }
mkdir -p "$out"
"$nvcc" -O3 -std=c++17 --fmad=false -arch=native --shared -Xcompiler=-fPIC \
  "$here/kernels.cu" -o "$out/libtrace_gpu.so.part"
mv "$out/libtrace_gpu.so.part" "$out/libtrace_gpu.so"
"$python" - "$here" "$out" <<'PY'
import hashlib,json,pathlib,sys
import tracehmm.trace,trace_utils
here,out=map(pathlib.Path,sys.argv[1:])
files=[pathlib.Path(tracehmm.trace.__file__),pathlib.Path(trace_utils.__file__),out/'libtrace_gpu.so',here/'kernels.cu',here/'backend.py',here.parent/'f/trace.py']
(out/'build-manifest.json').write_text(json.dumps(dict(backend='cuda-emission-transition-v1',sha256={str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in files}),indent=2)+'\n')
PY
TRACE_BACKEND=cuda "$python" "$here/../f/trace.py" infer --gpu-check
