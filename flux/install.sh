#!/usr/bin/env bash
# Explicit installer: analysis commands never install packages.
set -euo pipefail
FLUX_ROOT=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)
export FLUX_ROOT PYTHONNOUSERSITE=1
usage() {
    cat <<'HELP'
FLUX: reuse the existing ai environment; install only missing/incompatible dependencies.
Usage: ./install.sh [run] [options]
  --machine auto|local|hpc   Resolve shared filesystem roots.
  --device auto|cpu|gpu      auto: GPU on HPC; detected NVIDIA GPU locally.
  --cuda auto|cu118|cu128    auto: cu128 for Blackwell; cu118 for HPC/older GPUs.
  --check                   Check the environment without installing.
  --plan-only               Print environment and wheel selection only.
  -h, --help                Show help.
Default: discover existing ai in active/Anaconda/Miniforge/site environments.
FLUX_ENV overrides the prefix. No environment or Python is created/replaced.
HPC login nodes may install/check GPU wheels without an allocated GPU.
Examples:
  ./install.sh run
  ./install.sh run --check
  ./install.sh run --machine hpc --device gpu
HELP
}
[[ ${1:-} != run ]] || shift
device=auto; machine=auto; cuda=auto; check=false; plan=false
while (($#)); do
    case "$1" in
        --device) device=${2:?Missing device}; shift 2 ;;
        --machine) machine=${2:?Missing machine}; shift 2 ;;
        --cuda) cuda=${2:?Missing CUDA wheel}; shift 2 ;;
        --check) check=true; shift ;;
        --plan-only) plan=true; shift ;;
        -h|--help|help) usage; exit 0 ;;
        *) echo "ERROR: Unknown option: $1" >&2; exit 2 ;;
    esac
done
case "$device" in auto|cpu|gpu) ;; *) echo 'Invalid --device' >&2; exit 2 ;; esac
case "$cuda" in auto|cu118|cu128) ;; *) echo 'Invalid --cuda' >&2; exit 2 ;; esac
case "$machine" in auto|local|hpc) ;; *) echo 'Invalid --machine' >&2; exit 2 ;; esac
source "$FLUX_ROOT/f/runtime.sh"
flux_paths --machine "$machine"
if [[ $device == auto ]]; then
    if [[ $FLUX_MACHINE == hpc ]] || nvidia-smi -L >/dev/null 2>&1; then device=gpu; else device=cpu; fi
fi
if [[ $device == cpu ]]; then
    wheel=cpu
elif [[ $cuda != auto ]]; then
    wheel=$cuda
else
    wheel=cu118
    # A login node may expose no GPU. HPC defaults remain compatible with V100.
    caps=$(nvidia-smi --query-gpu=compute_cap --format=csv,noheader 2>/dev/null || true)
    while IFS= read -r cap; do
        major=${cap%%.*}
        if [[ $major =~ ^[0-9]+$ ]] && ((major >= 10)); then wheel=cu128; fi
    done <<< "$caps"
fi
printf '[flux install] machine=%s prefix=%s fallback-wheel=%s\n' "$FLUX_MACHINE" "$FLUX_ENV" "$wheel"
py=$FLUX_ENV/bin/python
[[ -x $py && -d $FLUX_ENV/conda-meta ]] || {
    echo "ERROR: Existing ai environment not found. Set FLUX_ENV to its absolute path." >&2
    exit 2
}
export LD_LIBRARY_PATH="$FLUX_ENV/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"
"$py" -I -c 'import sys,pathlib; assert pathlib.Path(sys.prefix).resolve()==pathlib.Path(sys.argv[1]).resolve(), "Python prefix mismatch"' "$FLUX_ENV"
args=(--wheel "$wheel")
if [[ $check != true && $plan != true ]]; then args+=(--install); fi
"$py" -I "$FLUX_ROOT/f/dependencies.py" "${args[@]}"
[[ $plan != true ]] || exit 0
unset FLUX_PYTHON
flux_runtime --machine "$machine"
args=(--prefix "$FLUX_ENV")
# An explicit CUDA request is a validation requirement, never a silent replacement.
if [[ $device == gpu && $cuda != auto ]]; then args+=(--expected-wheel "$cuda"); fi
if [[ $device == gpu && $FLUX_MACHINE != hpc ]]; then args+=(--require-gpu); fi
"$FLUX_PYTHON" -I "$FLUX_ROOT/f/check_install.py" "${args[@]}"
if ! "$py" -I -m pip check; then
    echo 'WARNING: Shared environment has dependency conflicts; review pip check output. Unrelated packages were not repaired.' >&2
fi
printf 'Ready: %s\nRun ./flux.sh run (no conda activation required).\n' "$FLUX_ENV"
