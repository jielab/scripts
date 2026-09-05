#!/usr/bin/env bash
# Shared configuration and option parsing for GRID modules.
set -euo pipefail
export LC_ALL=C
export PYTHONDONTWRITEBYTECODE=1
export MKL_NUM_THREADS=${MKL_NUM_THREADS:-1}
export OPENBLAS_NUM_THREADS=${OPENBLAS_NUM_THREADS:-1}
export NUMEXPR_NUM_THREADS=${NUMEXPR_NUM_THREADS:-1}
export OMP_NUM_THREADS=${OMP_NUM_THREADS:-1}

GRID_ROOT=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd -P)
GRID_F="$GRID_ROOT/f"
GRID_ORIGINAL_ARGS=("$@")

# Core defaults.
GRID_DATA_ROOT=${GRID_DATA_ROOT:-/mnt/d}
GRID_TRAIT=${GRID_TRAIT:-height}
GRID_TRAITS=${GRID_TRAITS:-height,ldl,t2dm}
GRID_POPS=${GRID_POPS:-AFR,EAS,EUR,SAS}
GRID_CHRS=${GRID_CHRS:-1-22}
GRID_JOBS=${GRID_JOBS:-4}
GRID_THREADS=${GRID_THREADS:-8}
GRID_REPLACE=${GRID_REPLACE:-FALSE}
GRID_DRY_RUN=${GRID_DRY_RUN:-FALSE}
GRID_GWAS_DIR=${GRID_GWAS_DIR:-/mnt/d/data.BIG/gwas/4grid}
if [[ -z ${GRID_TARGET_DIR+x} ]]; then
  GRID_TARGET_DIR=/mnt/d/data/ukb/gen/typ
  if [[ ! -s $GRID_TARGET_DIR/chr1.pgen && ! -s $GRID_TARGET_DIR/chr1.bed && -s /mnt/h/ukbGen/37/hap/chr1.pgen ]]; then
    GRID_TARGET_DIR=/mnt/h/ukbGen/37/hap
  fi
fi
GRID_IMP_DIR=${GRID_IMP_DIR:-/mnt/h/ukbGen/37/imp}
GRID_PHE_FILE=${GRID_PHE_FILE:-/mnt/d/data/ukb/phe/Rdata/phe.rds}
GRID_OUTPUT_ROOT=${GRID_OUTPUT_ROOT:-/mnt/d/analysis/grid}
GRID_KEEP=${GRID_KEEP:-}
GRID_REMOVE=${GRID_REMOVE:-}
GRID_N_GWAS=${GRID_N_GWAS:-}

# PRS-CSx.
GRID_CSX_REF_DIR=${GRID_CSX_REF_DIR:-/mnt/d/data.BIG/refLD/csx}
GRID_CSX_SNPINFO=${GRID_CSX_SNPINFO:-$GRID_CSX_REF_DIR/snpinfo_mult_1kg_hm3}
GRID_CSX_BIM_PREFIX=${GRID_CSX_BIM_PREFIX:-$GRID_TARGET_DIR/ukb_array}
GRID_PHI=${GRID_PHI:-1e-2}
GRID_MCMC_ITER=${GRID_MCMC_ITER:-4000}
GRID_MCMC_BURNIN=${GRID_MCMC_BURNIN:-2000}
GRID_MCMC_THIN=${GRID_MCMC_THIN:-5}
GRID_SEED=${GRID_SEED:-20260904}
GRID_SUMSTATS_CHUNK=${GRID_SUMSTATS_CHUNK:-500000}

# PCA/ancestry/DiscoDivas.
GRID_MED_FILE=${GRID_MED_FILE:-/mnt/d/files/DiscoDivas/med.g1000.4pop.tsv}
GRID_PCA_WEIGHT=${GRID_PCA_WEIGHT:-/mnt/d/files/DiscoDivas/g1k_hm3_maf5_woamb_wolr.pca.weight}
GRID_COV_PCS=${GRID_COV_PCS:-20}
GRID_DISTANCE_PCS=${GRID_DISTANCE_PCS:-10}
GRID_PCA_FILE=${GRID_PCA_FILE:-}
GRID_ANCESTRY_FILE=${GRID_ANCESTRY_FILE:-}
GRID_ANCESTRY_PROB_MIN=${GRID_ANCESTRY_PROB_MIN:-0.90}
GRID_ANCHOR_PROB_MIN=${GRID_ANCHOR_PROB_MIN:-0.999}
GRID_ANCHOR_MAX_PER_GROUP=${GRID_ANCHOR_MAX_PER_GROUP:-10000}

# ARG-Needle. The phased BGENs are UKB Field 22438, GRCh37.
GRID_ARG_ACTION=${GRID_ARG_ACTION:-check}
GRID_ARG_HAP_DIR=${GRID_ARG_HAP_DIR:-/mnt/h/ukbGen/37/hap}
GRID_ARG_MAP_DIR=${GRID_ARG_MAP_DIR:-}
GRID_ARG_MAP_PATTERN=${GRID_ARG_MAP_PATTERN:-}
if [[ -z ${GRID_ARG_OUT+x} ]]; then
  GRID_ARG_OUT="$(dirname -- "$GRID_ARG_HAP_DIR")/arg"
  GRID_ARG_OUT_AUTO=TRUE
else
  GRID_ARG_OUT_AUTO=FALSE
fi
GRID_ARG_SCRATCH=${GRID_ARG_SCRATCH:-}
GRID_ARG_FULL=${GRID_ARG_FULL:-FALSE}
GRID_ARG_MAX_INDIVIDUALS=${GRID_ARG_MAX_INDIVIDUALS:-20000}
GRID_ARG_SEED_HAPLOTYPES=${GRID_ARG_SEED_HAPLOTYPES:-4000}
GRID_ARG_THREADS=${GRID_ARG_THREADS:-8}
GRID_ARG_JOBS=${GRID_ARG_JOBS:-1}
GRID_ARG_ANCHORS_PER_POP=${GRID_ARG_ANCHORS_PER_POP:-1000}
GRID_ARG_WINDOW_BP=${GRID_ARG_WINDOW_BP:-1000000}
GRID_ARG_NEEDLE_HOME=${GRID_ARG_NEEDLE_HOME:-${ARG_NEEDLE_HOME:-}}
GRID_ARG_NORMALIZE=${GRID_ARG_NORMALIZE:-TRUE}
GRID_ARG_MAF_MIN=${GRID_ARG_MAF_MIN:-0.001}
GRID_ARG_GENO_MAX=${GRID_ARG_GENO_MAX:-0.05}
GRID_ARG_KEEP=${GRID_ARG_KEEP:-}
GRID_ARG_TREES_DIR=${GRID_ARG_TREES_DIR:-$GRID_ARG_OUT/trees}
GRID_ARG_AFFINITY=${GRID_ARG_AFFINITY:-FALSE}

# Evolutionary transport model and GRID scoring.
GRID_LDSCORE_DIR=${GRID_LDSCORE_DIR:-}
GRID_REQUIRE_LD=${GRID_REQUIRE_LD:-TRUE}
GRID_EXTERNAL_AGE=${GRID_EXTERNAL_AGE:-}
GRID_MAX_SNPS_PER_CHR=${GRID_MAX_SNPS_PER_CHR:-0}
GRID_RIDGE_ALPHA=${GRID_RIDGE_ALPHA:-10}
GRID_CONSERVATION_MIN=${GRID_CONSERVATION_MIN:-0.05}
GRID_CONSERVATION_MAX=${GRID_CONSERVATION_MAX:-0.995}
GRID_LOCAL=${GRID_LOCAL:-FALSE}

# Evaluation.
GRID_EVAL_REPEATS=${GRID_EVAL_REPEATS:-100}
GRID_EVAL_FOLDS=${GRID_EVAL_FOLDS:-5}
GRID_MIN_GROUP=${GRID_MIN_GROUP:-100}
GRID_EVAL_TUNED=${GRID_EVAL_TUNED:-TRUE}
GRID_EVAL_ZERO_SHOT=${GRID_EVAL_ZERO_SHOT:-TRUE}
GRID_WRITE_PREDICTIONS=${GRID_WRITE_PREDICTIONS:-FALSE}
GRID_ETHNICITY_COL=${GRID_ETHNICITY_COL:-auto}
GRID_PHENOTYPE_COL=${GRID_PHENOTYPE_COL:-auto}
GRID_EVAL_COVARIATES=${GRID_EVAL_COVARIATES:-auto}

_grid_die(){ echo "ERROR: $*" >&2; exit 2; }
_grid_need_value(){ [[ $# -ge 2 && -n ${2:-} && ${2:-} != --* ]] || _grid_die "$1 requires a value"; }
grid_bool(){ case "${1^^}" in TRUE|T|1|YES|Y|ON) echo TRUE;; FALSE|F|0|NO|N|OFF) echo FALSE;; *) return 2;; esac; }
grid_is_true(){ [[ $(grid_bool "$1") == TRUE ]]; }
grid_require_bool(){ local x; x=$(grid_bool "$2") || _grid_die "$1 must be TRUE/FALSE"; printf '%s\n' "$x"; }

grid_parse_args(){
  while (($#)); do
    case "$1" in
      --trait) _grid_need_value "$@"; GRID_TRAIT=${2,,}; shift 2;;
      --traits) _grid_need_value "$@"; GRID_TRAITS=${2,,}; shift 2;;
      --pops) _grid_need_value "$@"; GRID_POPS=${2^^}; shift 2;;
      --chrs|--chr|--arg-chrs) _grid_need_value "$@"; GRID_CHRS=$2; shift 2;;
      --jobs) _grid_need_value "$@"; GRID_JOBS=$2; shift 2;;
      --threads) _grid_need_value "$@"; GRID_THREADS=$2; shift 2;;
      --data-root) _grid_need_value "$@"; GRID_DATA_ROOT=${2%/}; shift 2;;
      --dir-gwas) _grid_need_value "$@"; GRID_GWAS_DIR=${2%/}; shift 2;;
      --dir-gen) _grid_need_value "$@"; GRID_TARGET_DIR=${2%/}; shift 2;;
      --dir-imp) _grid_need_value "$@"; GRID_IMP_DIR=${2%/}; shift 2;;
      --phe-file|--phe) _grid_need_value "$@"; GRID_PHE_FILE=$2; shift 2;;
      --output-root|--output-dir) _grid_need_value "$@"; GRID_OUTPUT_ROOT=${2%/}; shift 2;;
      --keep|--target-ids) _grid_need_value "$@"; GRID_KEEP=$2; shift 2;;
      --remove|--exclude-ids) _grid_need_value "$@"; GRID_REMOVE=$2; shift 2;;
      --n-gwas) _grid_need_value "$@"; GRID_N_GWAS=$2; shift 2;;
      --replace) _grid_need_value "$@"; GRID_REPLACE=$(grid_require_bool "$1" "$2"); shift 2;;
      --dry-run) _grid_need_value "$@"; GRID_DRY_RUN=$(grid_require_bool "$1" "$2"); shift 2;;

      --csx-ref-dir|--ref-dir) _grid_need_value "$@"; GRID_CSX_REF_DIR=${2%/}; shift 2;;
      --csx-snpinfo) _grid_need_value "$@"; GRID_CSX_SNPINFO=$2; shift 2;;
      --csx-bim-prefix|--bim-prefix) _grid_need_value "$@"; GRID_CSX_BIM_PREFIX=$2; shift 2;;
      --phi) _grid_need_value "$@"; GRID_PHI=$2; shift 2;;
      --mcmc-iter) _grid_need_value "$@"; GRID_MCMC_ITER=$2; shift 2;;
      --mcmc-burnin) _grid_need_value "$@"; GRID_MCMC_BURNIN=$2; shift 2;;
      --mcmc-thin) _grid_need_value "$@"; GRID_MCMC_THIN=$2; shift 2;;
      --seed) _grid_need_value "$@"; GRID_SEED=$2; shift 2;;
      --sumstats-chunk) _grid_need_value "$@"; GRID_SUMSTATS_CHUNK=$2; shift 2;;

      --med-file) _grid_need_value "$@"; GRID_MED_FILE=$2; shift 2;;
      --pca-file) _grid_need_value "$@"; GRID_PCA_FILE=$2; shift 2;;
      --pca-weight) _grid_need_value "$@"; GRID_PCA_WEIGHT=$2; shift 2;;
      --cov-pcs) _grid_need_value "$@"; GRID_COV_PCS=$2; shift 2;;
      --distance-pcs) _grid_need_value "$@"; GRID_DISTANCE_PCS=$2; shift 2;;
      --ancestry-file|--auto-ancestry-file) _grid_need_value "$@"; GRID_ANCESTRY_FILE=$2; shift 2;;
      --ancestry-prob-min|--ancestry-prob-threshold) _grid_need_value "$@"; GRID_ANCESTRY_PROB_MIN=$2; shift 2;;
      --anchor-prob-min) _grid_need_value "$@"; GRID_ANCHOR_PROB_MIN=$2; shift 2;;
      --anchor-max-per-group) _grid_need_value "$@"; GRID_ANCHOR_MAX_PER_GROUP=$2; shift 2;;
      --ethnicity-col) _grid_need_value "$@"; GRID_ETHNICITY_COL=$2; shift 2;;
      --phenotype-col) _grid_need_value "$@"; GRID_PHENOTYPE_COL=$2; shift 2;;
      --eval-covariates) _grid_need_value "$@"; GRID_EVAL_COVARIATES=$2; shift 2;;

      --arg-action) _grid_need_value "$@"; GRID_ARG_ACTION=${2,,}; shift 2;;
      --arg-hap-dir|--arg-hap-root|--ukb-hap-root)
        _grid_need_value "$@"; GRID_ARG_HAP_DIR=${2%/}
        if [[ $GRID_ARG_OUT_AUTO == TRUE ]]; then GRID_ARG_OUT="$(dirname -- "$GRID_ARG_HAP_DIR")/arg"; GRID_ARG_TREES_DIR="$GRID_ARG_OUT/trees"; fi
        shift 2;;
      --arg-map-dir) _grid_need_value "$@"; GRID_ARG_MAP_DIR=${2%/}; shift 2;;
      --arg-map-pattern) _grid_need_value "$@"; GRID_ARG_MAP_PATTERN=$2; shift 2;;
      --arg-out|--arg-dir) _grid_need_value "$@"; GRID_ARG_OUT=${2%/}; GRID_ARG_OUT_AUTO=FALSE; GRID_ARG_TREES_DIR="$GRID_ARG_OUT/trees"; shift 2;;
      --arg-trees-dir) _grid_need_value "$@"; GRID_ARG_TREES_DIR=${2%/}; shift 2;;
      --arg-scratch) _grid_need_value "$@"; GRID_ARG_SCRATCH=${2%/}; shift 2;;
      --arg-full) _grid_need_value "$@"; GRID_ARG_FULL=$(grid_require_bool "$1" "$2"); shift 2;;
      --arg-max-individuals|--arg-max-samples) _grid_need_value "$@"; GRID_ARG_MAX_INDIVIDUALS=$2; shift 2;;
      --arg-seed-haplotypes) _grid_need_value "$@"; GRID_ARG_SEED_HAPLOTYPES=$2; shift 2;;
      --arg-threads) _grid_need_value "$@"; GRID_ARG_THREADS=$2; shift 2;;
      --arg-jobs) _grid_need_value "$@"; GRID_ARG_JOBS=$2; shift 2;;
      --arg-anchors-per-pop) _grid_need_value "$@"; GRID_ARG_ANCHORS_PER_POP=$2; shift 2;;
      --arg-window-bp) _grid_need_value "$@"; GRID_ARG_WINDOW_BP=$2; shift 2;;
      --arg-needle-home) _grid_need_value "$@"; GRID_ARG_NEEDLE_HOME=${2%/}; shift 2;;
      --arg-normalize) _grid_need_value "$@"; GRID_ARG_NORMALIZE=$(grid_require_bool "$1" "$2"); shift 2;;
      --arg-maf-min) _grid_need_value "$@"; GRID_ARG_MAF_MIN=$2; shift 2;;
      --arg-geno-max) _grid_need_value "$@"; GRID_ARG_GENO_MAX=$2; shift 2;;
      --arg-keep) _grid_need_value "$@"; GRID_ARG_KEEP=$2; shift 2;;
      --arg-affinity) _grid_need_value "$@"; GRID_ARG_AFFINITY=$(grid_require_bool "$1" "$2"); shift 2;;

      --ldscore-dir) _grid_need_value "$@"; GRID_LDSCORE_DIR=${2%/}; shift 2;;
      --require-ld) _grid_need_value "$@"; GRID_REQUIRE_LD=$(grid_require_bool "$1" "$2"); shift 2;;
      --external-age) _grid_need_value "$@"; GRID_EXTERNAL_AGE=$2; shift 2;;
      --grid-max-snps-per-chr) _grid_need_value "$@"; GRID_MAX_SNPS_PER_CHR=$2; shift 2;;
      --grid-ridge-alpha) _grid_need_value "$@"; GRID_RIDGE_ALPHA=$2; shift 2;;
      --conservation-min) _grid_need_value "$@"; GRID_CONSERVATION_MIN=$2; shift 2;;
      --conservation-max) _grid_need_value "$@"; GRID_CONSERVATION_MAX=$2; shift 2;;
      --grid-local) _grid_need_value "$@"; GRID_LOCAL=$(grid_require_bool "$1" "$2"); shift 2;;

      --eval-repeats) _grid_need_value "$@"; GRID_EVAL_REPEATS=$2; shift 2;;
      --eval-folds) _grid_need_value "$@"; GRID_EVAL_FOLDS=$2; shift 2;;
      --min-group) _grid_need_value "$@"; GRID_MIN_GROUP=$2; shift 2;;
      --eval-tuned) _grid_need_value "$@"; GRID_EVAL_TUNED=$(grid_require_bool "$1" "$2"); shift 2;;
      --eval-zero-shot) _grid_need_value "$@"; GRID_EVAL_ZERO_SHOT=$(grid_require_bool "$1" "$2"); shift 2;;
      --write-predictions) _grid_need_value "$@"; GRID_WRITE_PREDICTIONS=$(grid_require_bool "$1" "$2"); shift 2;;

      -h|--help) exec bash "$GRID_F/help.sh";;
      *) _grid_die "unknown option '$1'";;
    esac
  done
  GRID_OUTPUT_ROOT=${GRID_OUTPUT_ROOT%/}
  [[ -n $GRID_LDSCORE_DIR ]] || GRID_LDSCORE_DIR="$GRID_OUTPUT_ROOT/reference/ldscore"
  [[ -n $GRID_PCA_FILE ]] || GRID_PCA_FILE="$GRID_DATA_ROOT/data/ukb/phe/pca/ukb.discodivas.pca.tsv.gz"
  [[ -n $GRID_ANCESTRY_FILE ]] || GRID_ANCESTRY_FILE="$(dirname -- "$GRID_PCA_FILE")/ukb.ancestry.auto.tsv.gz"
  [[ $GRID_JOBS =~ ^[1-9][0-9]*$ ]] || _grid_die "--jobs must be a positive integer"
  [[ $GRID_THREADS =~ ^[1-9][0-9]*$ ]] || _grid_die "--threads must be a positive integer"
  [[ $GRID_ARG_THREADS =~ ^[1-9][0-9]*$ ]] || _grid_die "--arg-threads must be a positive integer"
  [[ $GRID_ARG_JOBS =~ ^[1-9][0-9]*$ ]] || _grid_die "--arg-jobs must be a positive integer"
}

grid_csv_words(){ printf '%s' "$1" | tr ',;' '  '; }
grid_expand_chrs(){
  python3 - "$1" <<'PY'
import re,sys
s=sys.argv[1].replace(',', ' ').replace(';',' ')
out=[]
for tok in s.split():
    m=re.fullmatch(r'(\d+)-(\d+)',tok)
    if m:
        a,b=map(int,m.groups()); out.extend(map(str,range(a,b+1)))
    else: out.append(tok.replace('chr','').replace('CHR',''))
seen=[]
for x in out:
    if x not in seen: seen.append(x)
print(' '.join(seen))
PY
}

grid_find_gwas(){
  local trait=${1,,} pop=${2^^} f
  f=$(find -L "$GRID_GWAS_DIR" -maxdepth 1 -type f -iname "${trait}.${pop}.gz" -print -quit 2>/dev/null || true)
  [[ -n $f ]] || return 1
  printf '%s\n' "$f"
}

grid_bgen_for(){ find "$GRID_ARG_HAP_DIR" -maxdepth 1 -type f -name "ukb22438_c${1}_b0_v2.bgen" -print -quit; }
grid_pfile_for(){
  local p="$GRID_ARG_HAP_DIR/chr${1}"
  [[ -s $p.pgen && ( -s $p.pvar || -s $p.pvar.zst ) && -s $p.psam ]] && printf '%s\n' "$p"
}
grid_sample_for(){
  local p
  p=$(find "$GRID_ARG_HAP_DIR" -maxdepth 1 -type f -name "ukb22438_c${1}_b0_v2_s*.sample" -print -quit)
  if [[ -n $p ]]; then printf '%s\n' "$p"
  elif [[ -s $GRID_ARG_HAP_DIR/chr${1}.psam ]]; then printf '%s\n' "$GRID_ARG_HAP_DIR/chr${1}.psam"
  fi
}

grid_map_for(){
  local c=$1 p
  if [[ -n $GRID_ARG_MAP_PATTERN ]]; then
    p=${GRID_ARG_MAP_PATTERN//\{chr\}/$c}; p=${p//%CHR%/$c}; [[ -s $p ]] && { printf '%s\n' "$p"; return; }
  fi
  [[ -n $GRID_ARG_MAP_DIR ]] || return 1
  for p in \
    "$GRID_ARG_MAP_DIR/chr${c}.map" "$GRID_ARG_MAP_DIR/chr${c}.txt" "$GRID_ARG_MAP_DIR/chr${c}.map.gz" \
    "$GRID_ARG_MAP_DIR/chr${c}.b37.gmap.gz" \
    "$GRID_ARG_MAP_DIR/genetic_map_GRCh37_chr${c}.txt" "$GRID_ARG_MAP_DIR/genetic_map_GRCh37_chr${c}.txt.gz" \
    "$GRID_ARG_MAP_DIR/genetic_map_chr${c}_combined_b37.txt" "$GRID_ARG_MAP_DIR/genetic_map_chr${c}_combined_b37.txt.gz" \
    "$GRID_ARG_MAP_DIR/genetic_map_chr${c}_combined_b37.txt.map"; do
    [[ -s $p ]] && { printf '%s\n' "$p"; return; }
  done
  return 1
}

grid_target_mode(){
  local c=$1 d=$GRID_TARGET_DIR
  if [[ -s $d/chr$c.pgen && ( -s $d/chr$c.pvar || -s $d/chr$c.pvar.zst ) && -s $d/chr$c.psam ]]; then echo pfile
  elif [[ -s $d/chr$c.bed && -s $d/chr$c.bim && -s $d/chr$c.fam ]]; then echo bfile
  else return 1; fi
}

grid_run(){
  printf '[GRID]'; printf ' %q' "$@"; printf '\n' >&2
  grid_is_true "$GRID_DRY_RUN" || "$@"
}

grid_run_logged(){
  local log=$1; shift
  mkdir -p "$(dirname "$log")"
  printf '[GRID]'; printf ' %q' "$@"; printf '\n' | tee -a "$log" >&2
  grid_is_true "$GRID_DRY_RUN" || "$@" >>"$log" 2>&1
}

grid_checksum_record(){
  local out=$1; shift
  mkdir -p "$(dirname "$out")"
  { date -Is; for f in "$@"; do [[ -e $f ]] && stat -c '%n\t%s\t%Y' "$f"; done; } > "$out"
}
