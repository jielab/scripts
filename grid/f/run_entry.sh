#!/usr/bin/env bash
# Shared runtime setup; the caller chooses its independent dispatcher.
set -euo pipefail
umask 077
ROOT=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd -P)
entry=${1:?Missing dispatcher name}; shift
DRIVER="$ROOT/f/${entry}_driver.py"
[[ -f $DRIVER ]] || { echo "Unknown dispatcher: $entry" >&2; exit 2; }
for token in "$@"; do
    if [[ $token == --help || $token == -h ]]; then
        source "$ROOT/../0f/0os.sh"
        exec "${GRID_PLANNER_PYTHON:-$SHARED_PYTHON}" "$DRIVER" "$@"
    fi
done
requested_mode=${submit:-${SUBMIT:-auto}}
user_python=${GRID_PYTHON:-}
[[ $user_python != "${_OS_DEFAULT_GRID_PYTHON:-__unset__}" ]] || user_python=
plan_only=false
arguments=("$@")
core_options=()
for ((i=0; i<${#arguments[@]}; i++)); do
    if [[ ${arguments[i]} == --submit ]]; then requested_mode=${arguments[i+1]:?Missing --submit value}; fi
    if [[ ${arguments[i]} == --submit=* ]]; then requested_mode=${arguments[i]#--submit=}; fi
    if [[ ${arguments[i]} == --plan-only ]]; then plan_only=true; fi
    case "${arguments[i]}" in
        --cores|--job-cores) core_options+=(--cores "${arguments[i+1]:?Missing --cores value}") ;;
        --cores=*|--job-cores=*) core_options+=(--cores "${arguments[i]#*=}") ;;
    esac
done
export submit=$requested_mode SUBMIT=$requested_mode
[[ -f $ROOT/../0f/0os.sh ]] || { echo 'Install grid/ beside the repository shared 0f/ directory.' >&2; exit 2; }
source "$ROOT/../0f/0os.sh"
# Let the shared parser own CPU defaults, validation and explicit overrides.
# Other GRID flags (notably --wait) remain with the existing Python dispatcher.
workflow_args --submit "$requested_mode" "${core_options[@]}"
if [[ -n $user_python ]]; then
    PYTHON=$user_python
elif [[ -x ${DIR_WORK}/.conda/envs/grid/bin/python ]]; then
    PYTHON=${DIR_WORK}/.conda/envs/grid/bin/python
elif [[ -x $HOME/miniforge3/envs/grid/bin/python ]]; then
    PYTHON=$HOME/miniforge3/envs/grid/bin/python
elif [[ -x ${GRID_CONDA_ENV:-$HOME/miniforge3/envs/grid}/bin/python ]]; then
    PYTHON=${GRID_CONDA_ENV:-$HOME/miniforge3/envs/grid}/bin/python
else
    if [[ $submit == hpc ]]; then PYTHON=${DIR_WORK}/.conda/envs/grid/bin/python
    else PYTHON=python3
    fi
fi
export GRID_WORKER_PYTHON=$PYTHON
if command -v -- "$PYTHON" >/dev/null 2>&1; then
    PYTHON=$(command -v -- "$PYTHON")
elif [[ $plan_only == true ]]; then
    # A local HPC preview retains the intended worker interpreter in its plan.
    PYTHON=$(command -v python3)
else
    echo "Missing worker Python: $PYTHON. Activate grid or set GRID_PYTHON." >&2
    exit 2
fi
export GRID_PYTHON=$PYTHON
export PATH="$(dirname -- "$PYTHON"):$PATH"
# HPC system libstdc++ is older than the environment LightGBM binary.
python_lib="$(dirname -- "$(dirname -- "$PYTHON")")/lib"
if [[ -f $python_lib/libstdc++.so.6 ]]; then
    export LD_LIBRARY_PATH="$python_lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"
fi
export PYTHONDONTWRITEBYTECODE=1 PYTHONUNBUFFERED=1
export OPENBLAS_NUM_THREADS=${OPENBLAS_NUM_THREADS:-1} OMP_NUM_THREADS=${OMP_NUM_THREADS:-1} MKL_NUM_THREADS=${MKL_NUM_THREADS:-1}
exec "$PYTHON" "$DRIVER" "$@"
