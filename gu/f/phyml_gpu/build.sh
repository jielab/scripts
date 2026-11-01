#!/usr/bin/env bash


# Direct CLI help must return before setup, downloads or analysis.
if [[ "${BASH_SOURCE[0]}" == "$0" ]]; then
	case "${1:-}" in
		-h | --help)
			cat <<'HELP'
Modules:
    Default  Build and check the PhyML GPU executable

Main options:
    PHYML_GPU_SOURCE / PHYML_GPU_HOME  Source/install paths
    GU_NVCC  CUDA compiler
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
soft=${GU_SOFT:-$HOME/.local/share/gu}
src=${PHYML_GPU_SOURCE:-$soft/phyml-gpu-src}
out=${PHYML_GPU_HOME:-$soft/phyml-gpu}
nvcc=$(command -v "${GU_NVCC:-nvcc}" || true)
[[ -n $nvcc ]] || { echo 'ERROR: nvcc is required; activate a CUDA compiler environment or set GU_NVCC' >&2; exit 2; }
nvcc=$(readlink -f -- "$nvcc")
commit=42bb4b52f479938a3427dfe4357b861f0e306686
mkdir -p "$out"
exec 9>"$out/.build.lock"
flock -x 9
if [[ ! -d $src/.git ]]; then
  git clone https://github.com/stephaneguindon/phyml.git "$src"
  git -C "$src" checkout "$commit"
fi
[[ $(git -C "$src" rev-parse HEAD) == "$commit" ]] || { echo "ERROR: PhyML source must be pinned to $commit; use a fresh PHYML_GPU_SOURCE" >&2; exit 2; }
python3 "$here/patch_source.py" "$src"
mkdir -p "$out"
cp "$here/gu_cuda.h" "$src/src/gu_cuda.h"
"$nvcc" -O3 -std=c++17 --fmad=false -arch=native -Xcompiler=-fPIC -c "$here/partials.cu" -o "$out/partials.o"
gcc -O3 -std=gnu99 -I"$src/src" -I"$src" -c "$here/adapter.c" -o "$out/adapter.o"
cuda_root=$(cd -- "$(dirname -- "$nvcc")/.." && pwd)
cuda_static=''
for candidate in "$cuda_root/lib/libcudart_static.a" "$cuda_root/lib64/libcudart_static.a" "$cuda_root/targets/x86_64-linux/lib/libcudart_static.a"; do
  if [[ -f $candidate ]]; then cuda_static=$candidate; break; fi
done
[[ -n $cuda_static ]] || { echo "ERROR: libcudart_static.a missing under $cuda_root" >&2; exit 2; }
( cd "$src"
  autoreconf --force --install >"$out/configure.log" 2>&1
  ./configure --enable-phyml --disable-native --disable-docs >>"$out/configure.log" 2>&1
  make -C src clean >"$out/build.log" 2>&1
  make -C src -j4 CFLAGS='-O3 -std=gnu99 -DDISABLE_NATIVE -DGU_CUDA -ffp-contract=off' \
    LIBS="$out/adapter.o $out/partials.o $cuda_static -lstdc++ -ldl -lrt -lpthread -lm" >>"$out/build.log" 2>&1
)
cp "$src/src/phyml" "$out/phyml-gpu.part"
mv "$out/phyml-gpu.part" "$out/phyml-gpu"
"$out/phyml-gpu" --gpu-check
python3 - "$here" "$src" "$out" <<'PY'
import hashlib,json,pathlib,subprocess,sys
here,src,out=map(pathlib.Path,sys.argv[1:])
files=[*here.glob('*.c'),*here.glob('*.cu'),here/'gu_cuda.h',here/'patch_source.py',out/'phyml-gpu']
(out/'build-manifest.json').write_text(json.dumps(dict(backend='cuda-pruning-v1',upstream_commit=subprocess.check_output(['git','-C',str(src),'rev-parse','HEAD'],text=True).strip(),sha256={str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in files}),indent=2)+'\n')
PY
printf 'Built %s\n' "$out/phyml-gpu"
