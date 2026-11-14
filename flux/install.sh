#!/usr/bin/env bash

# Display standalone help before loading runtimes or scheduling work.
if [[ ${BASH_SOURCE[0]} == "$0" && ( ${1:-} == -h || ${1:-} == --help ) ]]; then
    cat <<'PUBLIC_ENTRY_HELP'
FLUX runtime installation

Main options:
    --machine auto|local|hpc  Machine selection
    --device auto|cpu|gpu    PyTorch processor selection
    --check                  Check without installing

Installs into FLUX_ENV (default DIR_WORK/software/flux). Requires Python 3.10-3.12.
On HPC login, select --device gpu explicitly for the GPU wheel.

Examples:
    ./install.sh --machine hpc --device gpu --check
    ./install.sh --machine hpc --device gpu
PUBLIC_ENTRY_HELP
    exit 0
fi
# Explicit installer; flux.sh never downloads dependencies.
set -euo pipefail
FLUX_ROOT=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)
export FLUX_ROOT
device=auto
machine=auto
check=false
while (($#)); do
    case "$1" in
        --device) device=${2:?}; shift 2 ;;
        --machine) machine=${2:?}; shift 2 ;;
        --check) check=true; shift ;;
        -h|--help)
            bash "${BASH_SOURCE[0]}" --help
            exit 0 ;;
        *) echo "ERROR: Unknown option: $1" >&2; exit 2 ;;
    esac
done
source "$FLUX_ROOT/f/runtime.sh"
flux_runtime --machine "$machine"
case "$device" in auto|cpu|gpu) ;; *) echo 'Invalid --device' >&2; exit 2 ;; esac
if [[ $check == true ]]; then
    "$FLUX_PYTHON" "$FLUX_ROOT/f/check_install.py"
    exit
fi
if [[ $device == auto ]]; then
    if [[ $FLUX_MACHINE == hpc ]] || command -v nvidia-smi >/dev/null 2>&1; then device=gpu; else device=cpu; fi
fi
base_python=${FLUX_BASE_PYTHON:-${SHARED_PYTHON:-python3}}
"$base_python" -c 'import sys; assert (3,10) <= sys.version_info[:2] <= (3,12), "Use Python 3.10-3.12 for the pinned environment"'
if [[ ! -x "$FLUX_ENV/bin/python" ]]; then "$base_python" -m venv "$FLUX_ENV"; fi
py="$FLUX_ENV/bin/python"
"$py" -m pip install --upgrade pip
"$py" -m pip install -r "$FLUX_ROOT/f/require.txt"
if [[ $device == gpu ]]; then index=https://download.pytorch.org/whl/cu118; else index=https://download.pytorch.org/whl/cpu; fi
# CPU and CUDA wheels share the base version. Force selection from this index
# so a previous CPU install cannot satisfy a later explicit GPU installation.
"$py" -m pip install --force-reinstall torch==2.6.0 --index-url "${FLUX_TORCH_INDEX:-$index}"
"$py" "$FLUX_ROOT/f/check_install.py"
printf 'Installed: %s\n' "$FLUX_ENV"
