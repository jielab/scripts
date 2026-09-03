#!/bin/bash


# 🚩 Proteome black-hole workflow
set -euo pipefail
export LC_ALL=C

SCRIPT_PATH="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" >/dev/null 2>&1 && pwd)/$(basename -- "${BASH_SOURCE[0]}")"
SCRIPT_DIR=$(dirname "$SCRIPT_PATH")
F_DIR="$SCRIPT_DIR/f"

usage() {
cat <<'USAGE'
Usage:
  phole.sh [step] [options]

Steps:
  prepare       read <category>/<protein>/magma/<protein>.genes.out into a graph matrix
  fit           fit cis + alpha*beta + sparse-gamma decomposition
  validate      target-position permutation test on held-out proteins
  plot          make prediction and black-hole diagnostic figures
  all           prepare,fit,validate,plot [default]
  selftest      generate MAGMA-shaped synthetic data and run the full pipeline

Core options:
  --data-root PATH       data root [/mnt/d]
  --gwas-root PATH       PROT project root [data.BIG/gwas/prot]
  --category NAME        gwas_post output category [common]
  --protein-bed FILE     CHR START END PROTEIN [GENE] BED [files/ppp_3k.38.bed]
  --protein-map FILE     optional PROTEIN GENE table; overrides BED gene labels
  --gwas LIST            comma-separated protein/GWAS names [all discovered]
  --out PATH             PHOLE output root [data.BIG/gwas/prot/phole]
  --python PATH          Python executable [auto-detect a compatible environment]
  --replace TRUE/FALSE   recreate existing step outputs [FALSE]

Matrix options:
  --score zstat|logp     MAGMA graph score [zstat]
  --scale robust|none    within-protein score scaling [robust]
  --score-cap FLOAT      winsorization cap [12]
  --protein-build auto|37|38  protein BED build [auto from filename]
  --cis-flank INT        cis window around the child's coding gene [100000]
  --protein-filter MODE  none|any-bonf|cis-bonf|coding-bonf|target-bonf [target-bonf]
  --protein-alpha FLOAT  per-protein family-wise alpha for the MAGMA gate [0.05]
  --require-magma-done TRUE/FALSE  require magma.done beside genes.out [TRUE]

Model options:
  --train-frac FLOAT     fraction of retained proteins used for training [0.80]
  --n-train INT          optional exact training count; overrides --train-frac
  --split random|chromosome  cold-start split [random]
  --seed INT             random seed [20260828]
  --rank auto|INT        alpha*beta rank [auto]
  --rank-grid LIST       ranks for internal cold-start selection [2,4,8,12,16,24]
  --inner-frac FLOAT     training fraction reserved for rank selection [0.20]
  --bandwidths LIST      coding-position kernel scales in bp [1e5,1e6,1e7,1e8]
  --global-weight FLOAT  cross-chromosome/background kernel weight [0.02]
  --gamma-z FLOAT        sparse direct-edge residual threshold [4]
  --gamma-max-edges INT  maximum gamma edges written [1000000]
  --metric-pairs INT     sampled held-out trans edges for metrics [1000000]
  --perm INT             target-position permutations [100]

Examples:
  /mnt/d/scripts/phole/phole.sh all --gwas-root /mnt/d/data.BIG/gwas/prot \
    --protein-bed /mnt/d/files/ppp_3k.38.bed --train-frac 0.80

  # Treat cvd_cad exactly like one PROT protein child:
  /mnt/d/scripts/phole/phole.sh prepare \
    --gwas-root /mnt/d/data.BIG/gwas/main --gwas cvd_cad \
    --protein-bed /path/to/protein.b38.bed --out /mnt/d/data.BIG/gwas/main/phole

  /mnt/d/scripts/phole/phole.sh selftest --out /tmp/phole.selftest --replace TRUE

Scientific contract:
  MAGMA ZSTAT is a gene-level association score, not a signed abundance effect.
  PHOLE therefore calls the MAGMA analysis an association-strength graph and
  never interprets its sign as protein up/down regulation or as causal mediation.
USAGE
}


# 🚩 Defaults
dir0=/mnt/d
gwas_root=""
category=common
protein_bed=""
protein_map=""
gwas=""
out=""
python=""
step=all
step_set=0
replace=FALSE

score=zstat
scale=robust
score_cap=12
protein_build=auto
cis_flank=100000
protein_filter=target-bonf
protein_alpha=0.05
require_magma_done=TRUE

train_frac=0.80
n_train=""
split=random
seed=20260828
rank=auto
rank_grid=2,4,8,12,16,24
inner_frac=0.20
bandwidths=100000,1000000,10000000,100000000
global_weight=0.02
gamma_z=4
gamma_max_edges=1000000
metric_pairs=1000000
perm=100


# 🚩 Arguments
need_arg_value() {
  local opt="$1" val="${2-}"
  if [[ -z "$val" || "$val" == --* ]]; then
    echo "ERROR: missing value for $opt" >&2
    usage >&2
    exit 2
  fi
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    prepare|fit|validate|plot|all|selftest|*,*)
      (( step_set == 0 )) || { echo "ERROR: multiple step modules: $step and $1" >&2; exit 2; }
      step="$1"; step_set=1; shift;;
    --data-root) need_arg_value "$1" "${2-}"; dir0="$2"; shift 2;;
    --gwas-root) need_arg_value "$1" "${2-}"; gwas_root="$2"; shift 2;;
    --category) need_arg_value "$1" "${2-}"; category="$2"; shift 2;;
    --protein-bed) need_arg_value "$1" "${2-}"; protein_bed="$2"; shift 2;;
    --protein-map) need_arg_value "$1" "${2-}"; protein_map="$2"; shift 2;;
    --gwas) need_arg_value "$1" "${2-}"; gwas="$2"; shift 2;;
    --out) need_arg_value "$1" "${2-}"; out="$2"; shift 2;;
    --python) need_arg_value "$1" "${2-}"; python="$2"; shift 2;;
    --replace) need_arg_value "$1" "${2-}"; replace="$2"; shift 2;;
    --score) need_arg_value "$1" "${2-}"; score="$2"; shift 2;;
    --scale) need_arg_value "$1" "${2-}"; scale="$2"; shift 2;;
    --score-cap) need_arg_value "$1" "${2-}"; score_cap="$2"; shift 2;;
    --protein-build) need_arg_value "$1" "${2-}"; protein_build="$2"; shift 2;;
    --cis-flank) need_arg_value "$1" "${2-}"; cis_flank="$2"; shift 2;;
    --protein-filter) need_arg_value "$1" "${2-}"; protein_filter="$2"; shift 2;;
    --protein-alpha) need_arg_value "$1" "${2-}"; protein_alpha="$2"; shift 2;;
    --require-magma-done) need_arg_value "$1" "${2-}"; require_magma_done="$2"; shift 2;;
    --train-frac) need_arg_value "$1" "${2-}"; train_frac="$2"; shift 2;;
    --n-train) need_arg_value "$1" "${2-}"; n_train="$2"; shift 2;;
    --split) need_arg_value "$1" "${2-}"; split="$2"; shift 2;;
    --seed) need_arg_value "$1" "${2-}"; seed="$2"; shift 2;;
    --rank) need_arg_value "$1" "${2-}"; rank="$2"; shift 2;;
    --rank-grid) need_arg_value "$1" "${2-}"; rank_grid="$2"; shift 2;;
    --inner-frac) need_arg_value "$1" "${2-}"; inner_frac="$2"; shift 2;;
    --bandwidths) need_arg_value "$1" "${2-}"; bandwidths="$2"; shift 2;;
    --global-weight) need_arg_value "$1" "${2-}"; global_weight="$2"; shift 2;;
    --gamma-z) need_arg_value "$1" "${2-}"; gamma_z="$2"; shift 2;;
    --gamma-max-edges) need_arg_value "$1" "${2-}"; gamma_max_edges="$2"; shift 2;;
    --metric-pairs) need_arg_value "$1" "${2-}"; metric_pairs="$2"; shift 2;;
    --perm) need_arg_value "$1" "${2-}"; perm="$2"; shift 2;;
    -v|--version) cat "$SCRIPT_DIR/VERSION"; exit 0;;
    -h|--help) usage; exit 0;;
    *) echo "ERROR: unknown step/module or option: $1" >&2; usage >&2; exit 2;;
  esac
done


# 🚩 Validation and runtime
upper(){ echo "$1" | tr '[:lower:]' '[:upper:]'; }
replace=$(upper "$replace")
require_magma_done=$(upper "$require_magma_done")
[[ "$replace" == TRUE || "$replace" == FALSE ]] || { echo "ERROR: --replace must be TRUE or FALSE" >&2; exit 2; }
[[ "$require_magma_done" == TRUE || "$require_magma_done" == FALSE ]] || { echo "ERROR: --require-magma-done must be TRUE or FALSE" >&2; exit 2; }
[[ "$score" == zstat || "$score" == logp ]] || { echo "ERROR: --score must be zstat or logp" >&2; exit 2; }
[[ "$scale" == robust || "$scale" == none ]] || { echo "ERROR: --scale must be robust or none" >&2; exit 2; }
[[ "$protein_build" == auto || "$protein_build" == 37 || "$protein_build" == 38 ]] || { echo "ERROR: --protein-build must be auto, 37, or 38" >&2; exit 2; }
[[ "$protein_filter" == none || "$protein_filter" == any-bonf || "$protein_filter" == cis-bonf || "$protein_filter" == coding-bonf || "$protein_filter" == target-bonf ]] || { echo "ERROR: --protein-filter must be none, any-bonf, cis-bonf, coding-bonf, or target-bonf" >&2; exit 2; }
[[ "$split" == random || "$split" == chromosome ]] || { echo "ERROR: --split must be random or chromosome" >&2; exit 2; }
[[ -z "$n_train" || "$n_train" =~ ^[1-9][0-9]*$ ]] || { echo "ERROR: --n-train must be a positive integer" >&2; exit 2; }
[[ "$cis_flank" =~ ^[0-9]+$ ]] || { echo "ERROR: --cis-flank must be a non-negative integer" >&2; exit 2; }
[[ "$perm" =~ ^[1-9][0-9]*$ ]] || { echo "ERROR: --perm must be a positive integer" >&2; exit 2; }
awk -v x="$train_frac" 'BEGIN{exit !(x+0>0 && x+0<1)}' || { echo "ERROR: --train-frac must be in (0,1)" >&2; exit 2; }
awk -v x="$protein_alpha" 'BEGIN{exit !(x+0>0 && x+0<1)}' || { echo "ERROR: --protein-alpha must be in (0,1)" >&2; exit 2; }

[[ -z "$gwas_root" ]] && gwas_root="$dir0/data.BIG/gwas/prot"
[[ -z "$protein_bed" ]] && protein_bed="$dir0/files/ppp_3k.38.bed"
[[ -z "$out" ]] && out="$gwas_root/phole"

log(){ echo "[$(date '+%F %T')] $*" >&2; }
need_file(){ [[ -s "$1" ]] || { echo "ERROR: missing or empty file: $1" >&2; exit 1; }; }
has_step(){ [[ ",$step," == *",$1,"* || "$step" == all ]]; }

for requested in ${step//,/ }; do
  case "$requested" in prepare|fit|validate|plot|all|selftest) ;; *) echo "ERROR: invalid step: $requested" >&2; exit 2;; esac
done
need_file "$F_DIR/common.py"
required_packages=(numpy pandas)
if has_step plot || [[ "$step" == selftest ]]; then
  required_packages+=(matplotlib)
fi
if [[ "$step" == selftest ]]; then
  required_packages+=(scipy)
fi
python_has_packages() {
  local candidate="$1"
  command -v "$candidate" >/dev/null 2>&1 || return 1
  "$candidate" - "${required_packages[@]}" <<'PY'
import importlib.util
import sys
missing = [name for name in sys.argv[1:] if importlib.util.find_spec(name) is None]
sys.exit(1 if missing else 0)
PY
}
if [[ -n "$python" ]]; then
  command -v "$python" >/dev/null 2>&1 || { echo "ERROR: Python not found: $python" >&2; exit 1; }
  python_has_packages "$python" || {
    echo "ERROR: Python lacks required package(s) for step '$step': ${required_packages[*]} ($python)" >&2
    exit 1
  }
else
  python_candidates=(python3)
  [[ -z "${PHOLE_PYTHON:-}" ]] || python_candidates=("$PHOLE_PYTHON" "${python_candidates[@]}")
  [[ -x "${HOME}/anaconda3/bin/python" ]] && python_candidates+=("${HOME}/anaconda3/bin/python")
  [[ -x "${HOME}/miniconda3/bin/python" ]] && python_candidates+=("${HOME}/miniconda3/bin/python")
  [[ -x "${HOME}/miniforge3/bin/python" ]] && python_candidates+=("${HOME}/miniforge3/bin/python")
  for candidate in "${python_candidates[@]}"; do
    if python_has_packages "$candidate"; then
      python="$candidate"
      break
    fi
  done
  [[ -n "$python" ]] || {
    echo "ERROR: no compatible Python found; required packages for step '$step': ${required_packages[*]}" >&2
    echo "Use --python PATH or set PHOLE_PYTHON." >&2
    exit 1
  }
fi
log "Python=$python packages=${required_packages[*]}"


# 🚩 Pipeline commands
matrix_dir="$out/matrix"
model_dir="$out/model"
validation_dir="$out/validation"
plot_dir="$out/plot"

prepare_args=(
  --gwas-root "$gwas_root" --category "$category" --protein-bed "$protein_bed" --out "$matrix_dir"
  --score "$score" --scale "$scale" --score-cap "$score_cap"
  --cis-flank "$cis_flank" --protein-build "$protein_build"
  --protein-filter "$protein_filter" --protein-alpha "$protein_alpha"
  --require-magma-done "$require_magma_done" --replace "$replace"
)
[[ -z "$protein_map" ]] || prepare_args+=(--protein-map "$protein_map")
[[ -z "$gwas" ]] || prepare_args+=(--gwas "$gwas")

fit_args=(
  --matrix-dir "$matrix_dir" --out "$model_dir" --train-frac "$train_frac"
  --split "$split" --seed "$seed" --rank "$rank" --rank-grid "$rank_grid"
  --inner-frac "$inner_frac" --cis-flank "$cis_flank" --gamma-z "$gamma_z"
  --gamma-max-edges "$gamma_max_edges" --bandwidths "$bandwidths"
  --global-weight "$global_weight" --metric-pairs "$metric_pairs" --replace "$replace"
)
[[ -z "$n_train" ]] || fit_args+=(--n-train "$n_train")

validate_args=(
  --matrix-dir "$matrix_dir" --model-dir "$model_dir" --out "$validation_dir"
  --perm "$perm" --seed "$seed" --replace "$replace"
)

plot_args=(
  --matrix-dir "$matrix_dir" --model-dir "$model_dir" --validation-dir "$validation_dir"
  --out "$plot_dir" --seed "$seed" --replace "$replace"
)

if [[ "$step" == selftest ]]; then
  log "Create a MAGMA-shaped synthetic proteome"
  "$python" "$F_DIR/simulate.py" --out "$out/synthetic" --n-source 600 --n-protein 72 --rank 3 --seed "$seed"
  gwas_root="$out/synthetic/gwas"
  protein_bed="$out/synthetic/protein.b38.bed"
  matrix_dir="$out/matrix"; model_dir="$out/model"; validation_dir="$out/validation"; plot_dir="$out/plot"
  "$python" "$F_DIR/prepare.py" --gwas-root "$gwas_root" --category common --protein-bed "$protein_bed" --out "$matrix_dir" \
    --score zstat --scale robust --score-cap 12 --cis-flank 100000 --protein-build 38 --replace "$replace"
  "$python" "$F_DIR/fit.py" --matrix-dir "$matrix_dir" --out "$model_dir" --train-frac 0.80 \
    --split random --seed "$seed" --rank auto --rank-grid 2,3,4,6 --inner-frac 0.25 \
    --cis-flank 100000 --gamma-z 4 --gamma-max-edges 100000 --bandwidths "$bandwidths" \
    --global-weight "$global_weight" --metric-pairs 100000 --replace "$replace"
  "$python" "$F_DIR/validate.py" --matrix-dir "$matrix_dir" --model-dir "$model_dir" \
    --out "$validation_dir" --perm 20 --seed "$seed" --replace "$replace"
  "$python" "$F_DIR/plot.py" --matrix-dir "$matrix_dir" --model-dir "$model_dir" \
    --validation-dir "$validation_dir" --out "$plot_dir" --seed "$seed" --replace "$replace"
  log "SELFTEST DONE: $out"
  exit 0
fi

mkdir -p "$out"
training_spec="fraction=$train_frac"
[[ -z "$n_train" ]] || training_spec="n-train=$n_train (override)"
log "PHOLE step=$step gwas-root=$gwas_root category=$category out=$out score=$score protein-filter=$protein_filter training=$training_spec split=$split rank=$rank"
if has_step prepare; then
  need_file "$protein_bed"
  [[ -z "$protein_map" ]] || need_file "$protein_map"
  "$python" "$F_DIR/prepare.py" "${prepare_args[@]}"
fi
if has_step fit; then
  "$python" "$F_DIR/fit.py" "${fit_args[@]}"
fi
if has_step validate; then
  "$python" "$F_DIR/validate.py" "${validate_args[@]}"
fi
if has_step plot; then
  "$python" "$F_DIR/plot.py" "${plot_args[@]}"
fi
log "DONE"
