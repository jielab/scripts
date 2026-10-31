#!/usr/bin/env bash
# GRID novel-method entry.
#
# Scientific modules:
#   blocknet -> LD-block representation + predictors
#   twins    -> fixed-radius genetic reference support / n_twins
#   assess   -> predictability + incremental-benefit assessment and evaluation
#   final    -> read-only integration with 1.pca/2.csx/3.disco/4.prsformer outputs
#
# Published benchmark methods live outside this entry:
#   1.pca.sh, 2.csx.sh, 3.disco.sh, 4.prsformer.sh
set -euo pipefail
umask 077

usage() {
cat <<'HELP'
GRID modules:
  blocknet   Prepare/verify CSx-derived LD blocks and fit csx-linear/csx-blocknet
             (Clinical and csx-moe remain benchmark rows).
  twins      Build fixed-radius genetic reference support and count ALL eligible twins.
  assess     Fit count-only predictability/benefit heads; freeze selections; evaluate/report.
  final      READ-ONLY article summary across GRID and 1-4 upstream scripts; no model fitting.
  architecture  Outcome-free regional dimensionality from CSx posterior weights.
  all        blocknet -> twins -> assess -> architecture -> final.
  shiny      Result viewer.
  cleanup    Preview/release verified disposable caches.
  self-check Run code checks.
  abm        Backward-compatible alias of assess.

Examples:
  ./grid.sh --traits height,ldl,t2dm  # same as all
  ./grid.sh blocknet --traits height,ldl,t2dm --device auto
  ./grid.sh twins    --traits height,ldl,t2dm --device auto
  ./grid.sh assess   --traits height,ldl,t2dm --device auto
  ./grid.sh final    --traits height,ldl,t2dm
  ./grid.sh all      --traits height,ldl,t2dm --device auto
  ./grid.sh shiny --analysis-root /mnt/d/analysis/grid --reindex

Options:
  --traits CSV          height,ldl,t2dm (in the supplied order)
  --device auto|cpu|cuda Default: auto
  --analysis-root DIR   Default: /mnt/d/analysis/grid
  --config FILE_OR_DIR  Default: config/abm.json; custom roots need matching block paths
  --threads N           Default: 4
  --dry-run             Preview without preparing inputs, training or deleting caches
  <module> --help       Module-specific options

Important:
  GRID never trains PRSformer/PRSformer2.  If the shared preparation cache is
  absent, run:
      ./4.prsformer.sh prepare --traits height,ldl,t2dm
  An active local/kinship/ACTIVE.json is verified and applied automatically.
  It keeps all test IDs and excludes their development relatives.
  The cache is infrastructure only; PRSformer model predictions are NOT inputs
  to the primary GRID blocknet/twin-count analysis.
HELP
}

case "${1:-all}" in
  -h|--help|help) usage; exit 0 ;;
esac

ROOT=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)

export SCRIPT_FOREGROUND=1
export SCRIPT_LOG_DIR=${SCRIPT_LOG_DIR:-/mnt/d/analysis/grid/logs}
source "$ROOT/../0f/console.sh"

runtime=${GRID_CONDA_ENV:-$HOME/miniforge3/envs/grid}
if [[ -x "$runtime/bin/python3" ]]; then
  default_python="$runtime/bin/python3"
else
  default_python=python3
fi
PYTHON=${GRID_PYTHON:-$default_python}
PYTHON_PATH=$(command -v -- "$PYTHON")
export GRID_PYTHON="$PYTHON_PATH"
export PATH="$(dirname -- "$PYTHON_PATH"):$PATH"
export PYTHONDONTWRITEBYTECODE=1 PYTHONUNBUFFERED=1
export PYTHONPATH="$ROOT/f${PYTHONPATH:+:$PYTHONPATH}"

module=all
if [[ $# -gt 0 && $1 != -* ]]; then module=$1; shift; fi
args=("$@")
if [[ $module == all ]]; then
  for token in "${args[@]}"; do
    case "$token" in -h|--help|help) usage; exit 0 ;; esac
  done
fi

for token in "${args[@]}"; do
  [[ $token != --through && $token != --through=* ]] || {
    echo "ERROR: --through is internal to the blocknet/twins/assess module split." >&2
    exit 2
  }
done

case "$module" in
  blocknet|twins|assess|abm|all)
    for ((j=0; j<${#args[@]}; j++)); do
      token=${args[$j]}
      case "$token" in
        --predictor|--representation)
          (( j+1 < ${#args[@]} )) || { echo "ERROR: missing value for $token" >&2; exit 2; }
          value=${args[$((j+1))]}
          ;;
        --predictor=*|--representation=*) value=${token#*=} ;;
        *) continue ;;
      esac
      if [[ $token == --predictor* && $value != blocknet || $token == --representation* && $value != csx_ldblocks_fixed ]]; then
        echo "The publication GRID workflow is frozen to csx-blocknet + csx_ldblocks_fixed." >&2
        exit 2
      fi
    done
    ;;
esac

# final accepts only its own safe/read-only subset.
architecture_config="$ROOT/config/abm.blocks.json"
final_args=()
i=0
while (( i < ${#args[@]} )); do
  a=${args[$i]}
  case "$a" in
    --config)
      (( i + 1 < ${#args[@]} )) || { echo "ERROR: missing value for $a" >&2; exit 2; }
      architecture_config=${args[$((i+1))]}
      if [[ -d $architecture_config ]]; then architecture_config="$architecture_config/abm.blocks.json"; else architecture_config="$(dirname -- "$architecture_config")/abm.blocks.json"; fi
      i=$((i+2))
      ;;
    --config=*)
      architecture_config=${a#*=}
      if [[ -d $architecture_config ]]; then architecture_config="$architecture_config/abm.blocks.json"; else architecture_config="$(dirname -- "$architecture_config")/abm.blocks.json"; fi
      i=$((i+1))
      ;;
    --traits|--analysis-root)
      (( i + 1 < ${#args[@]} )) || { echo "ERROR: missing value for $a" >&2; exit 2; }
      final_args+=("$a" "${args[$((i+1))]}")
      i=$((i+2))
      ;;
    --traits=*|--analysis-root=*|--dry-run|-h|--help)
      final_args+=("$a")
      i=$((i+1))
      ;;
    *)
      i=$((i+1))
      ;;
  esac
done

# One public workflow per analysis root. Children inherit the held descriptor.
analysis_root=/mnt/d/analysis/grid
preview=false
for ((j=0; j<${#args[@]}; j++)); do
  case "${args[$j]}" in
    --analysis-root) analysis_root=${args[$((j+1))]} ;;
    --analysis-root=*) analysis_root=${args[$j]#*=} ;;
    --dry-run|-h|--help) preview=true ;;
  esac
done
case "$module" in
  all|blocknet|twins|assess|abm)
    if [[ $preview == false && ${GRID_PUBLIC_WORKFLOW_ROOT:-} != "$analysis_root" ]]; then
      mkdir -p "$analysis_root/cache/abm"
      exec {grid_workflow_lock}>>"$analysis_root/cache/abm/.publication.lock"
      flock -n "$grid_workflow_lock" || { echo "Another GRID publication workflow is already using $analysis_root" >&2; exit 2; }
      export GRID_PUBLIC_WORKFLOW_ROOT=$analysis_root
    fi
    ;;
esac

run_abm_through() {
  local phase=$1
  shift
  exec "$PYTHON_PATH" "$ROOT/f/twin.py" abm "$@" --through "$phase"
}

run_final() {
  # The dispatcher owns compare -> paper summary -> browser, exactly once.
  exec "$PYTHON_PATH" "$ROOT/f/final_dispatch.py" "$@"
}

case "$module" in
  blocknet)
    run_abm_through fit_predictors "${args[@]}"
    ;;
  twins)
    run_abm_through build_support "${args[@]}"
    ;;
  assess|abm)
    run_abm_through report "${args[@]}"
    ;;
  architecture)
    exec "$PYTHON_PATH" "$ROOT/f/regional_dimensionality.py" "${args[@]}"
    ;;
  final)
    run_final "${final_args[@]}"
    ;;
  all)
    "$0" blocknet "${args[@]}"
    "$0" twins "${args[@]}"
    "$0" assess "${args[@]}"
    "$0" architecture "${final_args[@]}" --config "${architecture_config}"
    "$0" final "${final_args[@]}"
    ;;
  shiny|cleanup|self-check)
    exec "$PYTHON_PATH" "$ROOT/f/twin.py" "$module" "${args[@]}"
    ;;
  *)
    echo "ERROR: unknown GRID module: $module" >&2
    usage >&2
    exit 2
    ;;
esac
