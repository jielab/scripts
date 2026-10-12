#!/usr/bin/env bash
set -euo pipefail
here=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
upstream=${IBDMIX_RUNTIME:-/mnt/d/software/gu/ibdmix}
nvcc=${IBDMIX_NVCC:-/home/huangj/anaconda3/envs/as3_mamba/bin/nvcc}
out=${IBDMIX_GPU_BIN:-/mnt/d/software/gu/ibdmix-gpu/ibdmix-gpu}
[[ -x $nvcc ]] || { echo "ERROR: set IBDMIX_NVCC to a CUDA compiler" >&2; exit 2; }
mkdir -p "$(dirname -- "$out")"
tmp=$out.build.$$
trap 'rm -f "$tmp"' EXIT
"$nvcc" -O3 -std=c++17 --fmad=false -arch=native \
  -I"$upstream/include" -I"$upstream/build/_deps/cli11-src/include" \
  "$here/caller.cu" "$upstream/src/lod_calculator.cc" \
  "$upstream/src/Mask_Reader.cc" "$upstream/src/Sample_Mapper.cc" -o "$tmp"
"$tmp" --gpu-check
mv -f "$tmp" "$out"
python3 - "$here" "$upstream" "$nvcc" "$out" <<'PY'
import hashlib,json,subprocess,sys
from pathlib import Path
source,upstream,compiler,binary=map(Path,sys.argv[1:])
files=[source/'caller.cu',source/'build.sh',upstream/'src/lod_calculator.cc',upstream/'src/Mask_Reader.cc',upstream/'src/Sample_Mapper.cc']
manifest={'backend':'cuda-segment-v1',
 'upstream_commit':subprocess.check_output(['git','-C',str(upstream),'rev-parse','HEAD'],text=True).strip(),
 'compiler':subprocess.check_output([str(compiler),'--version'],text=True).strip(),
 'flags':['-O3','-std=c++17','--fmad=false','-arch=native'],
 'files':{str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in files},
 'binary_sha256':hashlib.sha256(binary.read_bytes()).hexdigest()}
(binary.parent/'build-manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
PY
printf 'Built %s\n' "$out"
