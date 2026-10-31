#!/usr/bin/env bash
set -euo pipefail
HERE="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
SYSTEM=0
[[ "${1:-}" == "--system" ]] && SYSTEM=1

source "$HERE/f/cuda_env.sh"

# Apply network settings only to this installer; keep the APT download cache so
# rerunning after an interrupted download can reuse completed packages.
APT_OPTIONS=(
  -o Acquire::Retries=10
  -o Acquire::http::Timeout=120
  -o Acquire::https::Timeout=120
)
CUDA_CURL_OPTIONS=()
if [[ "${METAG_CUDA_DIRECT:-0}" == 1 ]]; then
  # NVIDIA can redirect downloads to its regional .cn endpoint.
  for cuda_host in developer.download.nvidia.com developer.download.nvidia.cn; do
    APT_OPTIONS+=(
      -o "Acquire::http::Proxy::$cuda_host=DIRECT"
      -o "Acquire::https::Proxy::$cuda_host=DIRECT"
    )
  done
  CUDA_CURL_OPTIONS+=(--noproxy "developer.download.nvidia.com,developer.download.nvidia.cn")
fi

system_apt() {
  local status
  if sudo apt-get "${APT_OPTIONS[@]}" "$@"; then
    return 0
  else
    status=$?
  fi
  echo "APT failed (exit $status); see the error above. Downloaded packages are kept for retry." >&2
  echo "After resolving it, rerun: ./install.sh --system" >&2
  if [[ "${METAG_CUDA_DIRECT:-0}" != 1 ]]; then
    echo "For NVIDIA proxy/download errors, try direct access: METAG_CUDA_DIRECT=1 ./install.sh --system" >&2
  fi
  return "$status"
}

install_cuda_toolkit() (
  # Keep OS variables and the temporary download local to this subshell.
  source /etc/os-release
  if [[ "${ID:-}" != ubuntu || "$(uname -m)" != x86_64 ]]; then
    echo "Automatic CUDA installation requires Ubuntu x86_64 (native or WSL2)." >&2
    echo "Install CUDA Toolkit >=12.8 for this platform, then rerun ./install.sh." >&2
    exit 1
  fi
  case "${VERSION_ID:-}" in
    20.04|22.04|24.04) ;;
    *) echo "Automatic CUDA 12.8 installation does not support Ubuntu $VERSION_ID." >&2; exit 1;;
  esac

  cuda_repo="ubuntu${VERSION_ID//./}"
  if [[ "$(uname -r)" == *[Mm]icrosoft* ]]; then
    cuda_repo="wsl-ubuntu"
  fi
  cuda_tmp="$(mktemp -d)"
  trap 'rm -rf -- "$cuda_tmp"' EXIT
  curl --fail --location --retry 3 --connect-timeout 30 \
    "${CUDA_CURL_OPTIONS[@]}" \
    --output "$cuda_tmp/cuda-keyring.deb" \
    "https://developer.download.nvidia.com/compute/cuda/repos/$cuda_repo/x86_64/cuda-keyring_1.1-1_all.deb"
  sudo dpkg -i "$cuda_tmp/cuda-keyring.deb"
  system_apt update
  # Toolkit only: WSL uses the Windows driver, never a Linux driver package.
  system_apt install -y --no-install-recommends cuda-toolkit-12-8
)

if ((SYSTEM)); then
  if [[ "${METAG_CUDA_DIRECT:-0}" == 1 ]]; then
    echo "Using direct connections for NVIDIA downloads."
  fi
  system_apt update
  system_apt install -y build-essential git zlib1g-dev curl wget aria2 ca-certificates
  if metag_use_cuda; then
    echo "Using existing CUDA Toolkit $METAG_CUDA_VERSION: $(command -v nvcc)"
  else
    echo "Installing CUDA Toolkit 12.8 from NVIDIA (required for Blackwell)."
    install_cuda_toolkit
    if ! metag_use_cuda; then
      echo "CUDA installation finished, but no nvcc >=12.8 was found." >&2
      exit 1
    fi
    echo "Using CUDA Toolkit $METAG_CUDA_VERSION: $(command -v nvcc)"
  fi
fi

CONDA=""
for c in "$HOME/miniforge3/bin/conda" "$HOME/miniconda3/bin/conda" "$(command -v conda 2>/dev/null || true)"; do
  if [[ -n "$c" && -x "$c" ]]; then CONDA="$c";break;fi
done
if [[ -z "$CONDA" ]]; then
  echo "Conda/Miniforge not found. Install Miniforge first." >&2
  exit 1
fi

ENV_PREFIX="$(dirname "$(dirname "$CONDA")")/envs/metag"
install_conda_env() (
  # Limit overrides to this operation and honor explicit caller settings.
  export CONDA_REMOTE_MAX_RETRIES="${CONDA_REMOTE_MAX_RETRIES:-5}"
  export CONDA_REMOTE_CONNECT_TIMEOUT_SECS="${CONDA_REMOTE_CONNECT_TIMEOUT_SECS:-30}"
  export CONDA_REMOTE_READ_TIMEOUT_SECS="${CONDA_REMOTE_READ_TIMEOUT_SECS:-120}"
  if [[ "${METAG_CONDA_DIRECT:-0}" == 1 ]]; then
    # requests uses no_proxy even when http(s)_proxy is set. Keep other bypasses.
    export no_proxy="${no_proxy:-${NO_PROXY:-}}"
    export no_proxy="${no_proxy:+$no_proxy,}conda.anaconda.org,repo.anaconda.com"
    export NO_PROXY="$no_proxy"
    echo "Using direct connections for Conda downloads."
  fi
  if [[ -f "$ENV_PREFIX/conda-meta/history" ]]; then
    "$CONDA" env update -p "$ENV_PREFIX" -f "$HERE/environment.yml" --prune
  else
    "$CONDA" env create -p "$ENV_PREFIX" -f "$HERE/environment.yml"
  fi
)
install_conda_env

# Use conda environment commands for the rest.
# Third-party Conda activation/deactivation hooks may read unset variables.
set +u
source "$(dirname "$(dirname "$CONDA")")/etc/profile.d/conda.sh"
conda activate "$ENV_PREFIX"
set -u

# Build pinned MetaCache. This never touches reference FASTA.
bash "$HERE/f/install_metacache.sh"

echo
echo "PASS software installation."
echo "After reference FASTA and NCBI metadata downloads finish (see 0.README.md):"
echo "  ./1.library.sh"
echo "Then classify and split samples:"
echo "  ./2.split.sh"
echo "Software, GPU and reference checks run automatically at startup."
