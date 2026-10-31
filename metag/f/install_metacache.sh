#!/usr/bin/env bash
set -euo pipefail

HERE="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
source "$HERE/cuda_env.sh"
metag_use_cuda || true
TOOLS="${METAG_TOOLS:-$HOME/.cache/metag/tools}"
MC="$TOOLS/metacache"
COMMIT="d7646eca4c4dc131262b16d2910923fce3f5d4fc"
THREADS="${THREADS:-16}"

mkdir -p "$TOOLS"
command -v git >/dev/null || { echo "git missing" >&2; exit 1; }
command -v make >/dev/null || { echo "make missing" >&2; exit 1; }
command -v g++ >/dev/null || { echo "g++ missing" >&2; exit 1; }

if [[ ! -d "$MC/.git" ]]; then
  git clone https://github.com/muellan/metacache.git "$MC"
fi

cd "$MC"
git fetch origin
git checkout --detach "$COMMIT"
git submodule update --init --recursive

# Preserve per-read reference evidence alongside upstream species classification.
if ! git apply --reverse --check "$HERE/metacache-reference-hits.patch" 2>/dev/null; then
  git apply --check "$HERE/metacache-reference-hits.patch"
  git apply "$HERE/metacache-reference-hits.patch"
fi

echo "[MetaCache CPU]"
make clean
make -j"$THREADS" cpu

GPU_OK=0
if command -v nvidia-smi >/dev/null 2>&1 && command -v nvcc >/dev/null 2>&1; then
  CUDA_VER="$(nvcc --version | sed -n 's/.*release \([0-9][0-9]*\.[0-9][0-9]*\).*/\1/p' | tail -1)"
  if python3 - "$CUDA_VER" <<'PY'
import sys
v=tuple(map(int,sys.argv[1].split('.')))
raise SystemExit(0 if v >= (12,8) else 1)
PY
  then
    CC="$(nvidia-smi --query-gpu=compute_cap --format=csv,noheader,nounits | head -1 | tr -d ' ')"
    ARCH="sm_${CC/./}"
    echo "[MetaCache GPU] CUDA=$CUDA_VER CC=$CC ARCH=$ARCH"
    # Do not destroy the already built CPU binary if GPU compile fails.
    if make -j"$THREADS" gpu CUDA_ARCH="$ARCH"; then
      GPU_OK=1
    else
      echo "WARNING: MetaCache GPU compile failed; CPU binary remains usable." >&2
    fi
  else
    echo "WARNING: CUDA toolkit $CUDA_VER is too old for Blackwell; need >=12.8. CPU MetaCache remains usable." >&2
  fi
else
  echo "WARNING: nvidia-smi or nvcc unavailable; CPU MetaCache only." >&2
fi

{
  echo "commit=$(git rev-parse HEAD)"
  echo "cpu_binary=$MC/metacache"
  echo "gpu_binary=$MC/metacache_gpu"
  echo "gpu_built=$GPU_OK"
  g++ --version | head -1
  nvcc --version 2>/dev/null | tail -1 || true
} > "$TOOLS/metacache.install.txt"

echo "PASS: MetaCache installed at $MC"
