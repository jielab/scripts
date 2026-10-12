#!/usr/bin/env bash
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
