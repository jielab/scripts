#!/usr/bin/env bash

# Public directory options use --dir-XXX; legacy spellings remain compatible.
if [[ "${BASH_SOURCE[0]}" == "$0" ]]; then
    source /mnt/d/scripts/0f/cli_dirs.sh
    script_directory_args "$(realpath -- "${BASH_SOURCE[0]}")" "$@"
    set -- "${SCRIPT_DIRECTORY_ARGS[@]}"
    unset SCRIPT_DIRECTORY_ARGS
fi

# Direct CLI help must return before setup, downloads or analysis.
if [[ "${BASH_SOURCE[0]}" == "$0" ]]; then
	case "${1:-}" in
		-h | --help)
			cat <<'HELP'
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
out=${TRACE_GPU_HOME:-/mnt/d/software/gu/trace-gpu}
nvcc=${GU_NVCC:-/home/huangj/anaconda3/envs/as3_mamba/bin/nvcc}
python=${GU_PYTHON:-/home/huangj/anaconda3/envs/gu/bin/python}
mkdir -p "$out"
"$nvcc" -O3 -std=c++17 --fmad=false -arch=native --shared -Xcompiler=-fPIC \
  "$here/kernels.cu" -o "$out/libtrace_gpu.so.part"
mv "$out/libtrace_gpu.so.part" "$out/libtrace_gpu.so"
"$python" - "$here" "$out" <<'PY'
import hashlib,json,pathlib,sys
import tracehmm.trace,trace_utils
here,out=map(pathlib.Path,sys.argv[1:])
files=[pathlib.Path(tracehmm.trace.__file__),pathlib.Path(trace_utils.__file__),out/'libtrace_gpu.so',here/'kernels.cu',here/'backend.py',here.parent/'trace.py']
(out/'build-manifest.json').write_text(json.dumps(dict(backend='cuda-emission-transition-v1',sha256={str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in files}),indent=2)+'\n')
PY
TRACE_BACKEND=cuda "$python" "$here/../trace.py" infer --gpu-check
