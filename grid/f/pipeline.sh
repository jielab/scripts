#!/usr/bin/env bash
# Shared launcher for PCA, PRS-CSx and official DiscoDivas.
set -euo pipefail
ROOT=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd -P)
method=$1; shift
source "$ROOT/f/environment.sh"
source "$ROOT/f/common.sh" "$@"
GRID_CHECK=FALSE
GRID_STAGE=all
GRID_SCORE_DIR=${GRID_SCORE_DIR:-/mnt/d/data/ukb/pgs}
GRID_DISCO_A=1,1,1,1
GRID_REGRESS_PCA=TRUE
GRID_CSX_MODELS=populations,auto,meta
args=()
while (($#)); do
  case "$1" in
    --check) GRID_CHECK=TRUE; shift;;
    --models) _grid_need_value "$@"; GRID_CSX_MODELS=$2; shift 2;;
    --stage) _grid_need_value "$@"; GRID_STAGE=$2; shift 2;;
    --score-dir) _grid_need_value "$@"; GRID_SCORE_DIR=${2%/}; shift 2;;
    --a-list) _grid_need_value "$@"; GRID_DISCO_A=$2; shift 2;;
    --regress-pca) _grid_need_value "$@"; GRID_REGRESS_PCA=$(grid_require_bool "$1" "$2"); shift 2;;
    --trait) _grid_need_value "$@"; GRID_TRAITS=${2,,}; shift 2;;
    --traits) _grid_need_value "$@"; GRID_TRAITS=${2,,}; shift 2;;
    *) args+=("$1"); shift;;
  esac
done
grid_parse_args "${args[@]}"
if [[ $method == csx ]]; then
  [[ -n $GRID_CSX_MODELS ]] || _grid_die '--models cannot be empty'
  for model in $(grid_csv_words "$GRID_CSX_MODELS"); do
    [[ $model == populations || $model == auto || $model == meta ]] || _grid_die '--models accepts populations,auto,meta'
    [[ $model != meta || $GRID_PHI != auto ]] || _grid_die 'Use --models populations,auto with --phi auto; meta requires a fixed phi'
  done
fi
[[ $GRID_STAGE == all || $GRID_STAGE == weights || $GRID_STAGE == score ]] || _grid_die '--stage must be all, weights or score'
[[ $method == csx || $GRID_STAGE == all ]] || _grid_die '--stage applies only to CSx'
# Validate before creating any output or launching computation.
[[ $GRID_POPS == AFR,EAS,EUR,SAS ]] || _grid_die 'This workflow requires --pops AFR,EAS,EUR,SAS (order matches the supplied reference centers)'
read -r -a CHRS <<< "$(grid_expand_chrs "$GRID_CHRS")"
((${#CHRS[@]})) || _grid_die 'Empty chromosome list'
for c in "${CHRS[@]}"; do [[ $c =~ ^([1-9]|1[0-9]|2[0-2])$ ]] || _grid_die "Invalid autosome: $c"; done
[[ $method != pca || ${CHRS[*]} == "$(seq -s ' ' 1 22)" ]] || _grid_die 'PCA requires all 22 autosomes'
POPS=(AFR EAS EUR SAS)
for t in $(grid_csv_words "$GRID_TRAITS"); do [[ $t == height || $t == ldl || $t == t2dm ]] || _grid_die "Unknown trait: $t"; done
[[ -n $GRID_TRAITS ]] || _grid_die 'Empty trait list'
work="$GRID_OUTPUT_ROOT/$method"
mkdir -p "$work/cmd" "$work/log"
run_id=$(date +%Y%m%dT%H%M%S).$$
GRID_COMMAND_FILE="$work/cmd/$method.$run_id.sh"
GRID_RUN_LOG="$work/log/$method.$run_id.log"
printf '#!/usr/bin/env bash\nset -euo pipefail\nsource %q\n' "$ROOT/f/environment.sh" > "$GRID_COMMAND_FILE"
exec > >(tee -a "$GRID_RUN_LOG") 2>&1
trap 'rc=$?; if ((rc)); then echo "ERROR: exit=$rc; details: $GRID_RUN_LOG"; fi' EXIT
source "$ROOT/f/pipeline_lib.sh"
echo "[$method] work directory: $work"
echo "[$method] command file: $GRID_COMMAND_FILE"
echo "[$method] log file: $GRID_RUN_LOG"
if [[ $method == pca ]]; then
  source "$ROOT/f/pca.sh"
else
  for trait in $(grid_csv_words "$GRID_TRAITS"); do
    if [[ $method == csx ]]; then
      for csx_model in $(grid_csv_words "$GRID_CSX_MODELS"); do
        case "$csx_model" in
          populations) ( GRID_TRAIT=$trait; source "$ROOT/f/csx.sh" );;
          auto|meta) ( GRID_TRAIT=$trait; source "$ROOT/f/csx_combined.sh" );;
          *) _grid_die '--models accepts populations,auto,meta';;
        esac
      done
    else
      ( GRID_TRAIT=$trait; source "$ROOT/f/$method.sh" )
    fi
  done
fi
