#!/usr/bin/env bash
# Shared by installation and runtime commands; sourcing this file installs nothing.
metag_use_cuda() {
  local candidate output major minor bin_dir cuda_root
  METAG_CUDA_VERSION=""

  # WSL exposes the Windows driver tools here, sometimes outside PATH.
  if [[ -x /usr/lib/wsl/lib/nvidia-smi && ":$PATH:" != *:/usr/lib/wsl/lib:* ]]; then
    export PATH="/usr/lib/wsl/lib:$PATH"
  fi

  for candidate in "${CUDA_HOME:-/usr/local/cuda}/bin/nvcc" \
      "$(command -v nvcc 2>/dev/null || true)" \
      /usr/local/cuda/bin/nvcc /usr/local/cuda-*/bin/nvcc; do
    [[ -x "$candidate" ]] || continue
    output="$("$candidate" --version 2>/dev/null)" || continue
    [[ "$output" =~ release[[:space:]]+([0-9]+)\.([0-9]+) ]] || continue
    major="${BASH_REMATCH[1]}"
    minor="${BASH_REMATCH[2]}"
    ((major > 12 || (major == 12 && minor >= 8))) || continue

    bin_dir="$(cd -- "$(dirname -- "$candidate")" && pwd -P)" || continue
    cuda_root="$(dirname -- "$bin_dir")"
    # MetaCache's host compiler uses CUDA_HOME/include for CUDA headers.
    export CUDA_HOME="$cuda_root"
    if [[ "$PATH" != "$bin_dir" && "$PATH" != "$bin_dir:"* ]]; then
      export PATH="$bin_dir:$PATH"
    fi
    if [[ -d "$cuda_root/lib64" && ":${LD_LIBRARY_PATH:-}:" != *:"$cuda_root/lib64":* ]]; then
      export LD_LIBRARY_PATH="$cuda_root/lib64${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"
    fi
    # shellcheck disable=SC2034 # Read by install.sh after this function succeeds.
    METAG_CUDA_VERSION="$major.$minor"
    return 0
  done
  return 1
}
