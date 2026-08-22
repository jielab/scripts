#!/usr/bin/env bash
set -euo pipefail

usage() {
    cat <<'USAGE'
Usage:
  ./birth.sh [order|season|all] [options]

Options:
  --phenome-source SOURCE  auto | phewas | fod       [default: auto]
  --min-cases N            minimum cases per phenotype [default: 500]
  --sib-cap N              sibling-count QC cap        [default: 10]
  --jobs N                 parallel workers            [default: 1; PheWAS is memory-heavy]
  --primary-model MODEL    m1 | m2 | m3                [default: m2]
  --monthly-top N          season: monthly ORs for top N phenotypes [default: 30]
  --strata-top N           stratified analyses for top N phenotypes [default: 20]
  --within MODE            auto | 1 | 0                [default: auto]
  -h, --help               show this help

Main model definitions:
  m1: sex + age + age^2 + 5-y birth cohort + ethnicity + follow-up duration
  m2: m1 + birth-TDI + birthplace easting/northing spatial terms   [primary]
  m3: m2 + current TDI + assessment centre                        [sensitivity]

Phenome source:
  auto   : use Rdata/PheWAS.rds if present; otherwise build phenotypes from
           all.rds fod_{icd10,icd9,srd,gp,ref}_* first-occurrence-date columns.
  phewas : require Rdata/PheWAS.rds.
  fod    : use the FOD phenotypes already merged into all.rds.

Useful environment overrides:
  DIR0=/mnt/d
  PHEDIR=$DIR0/data/ukb/phe
  GEN_DIR=$DIR0/data/ukb/gen
  UKB_OUT=$DIR0/analysis/birth
  BIRTH_R=./f/birth.R
  MIN_CASES_ORDER=500      # article-like between-family case threshold
  MIN_DISCORDANT_WITHIN=100 # article-like within-family discordant-pair threshold
  PHENO_MAP=$PHEDIR/common/birth.phenome.lst

Examples:
  ./birth.sh order
  ./birth.sh season --jobs 12
  ./birth.sh order --min-cases 500
  PHENOME_SOURCE=fod ./birth.sh all
USAGE
}

step="all"
step_set=0
phenome_source=${PHENOME_SOURCE:-auto}
min_cases=${MIN_CASES:-500}
sib_cap=${SIB_CAP:-10}
jobs=${N_JOBS:-1}
primary_model=${PRIMARY_MODEL:-m2}
monthly_top=${MONTHLY_TOP:-30}
strata_top=${STRATA_TOP:-20}
within_mode=${RUN_WITHIN:-auto}

while [[ $# -gt 0 ]]; do
    case "$1" in
        --phenome-source) phenome_source="$2"; shift 2 ;;
        --min-cases) min_cases="$2"; shift 2 ;;
        --sib-cap) sib_cap="$2"; shift 2 ;;
        --jobs) jobs="$2"; shift 2 ;;
        --primary-model) primary_model="$2"; shift 2 ;;
        --monthly-top) monthly_top="$2"; shift 2 ;;
        --strata-top) strata_top="$2"; shift 2 ;;
        --within) within_mode="$2"; shift 2 ;;
        -h|--help) usage; exit 0 ;;
        --*) echo "ERROR: unknown option: $1" >&2; usage >&2; exit 1 ;;
        order|season|all) (( step_set == 0 )) || { echo "ERROR: multiple modules specified: $step and $1" >&2; exit 1; }; step="$1"; step_set=1; shift ;;
        *) echo "ERROR: unknown argument: $1" >&2; usage >&2; exit 1 ;;
    esac
done

case "$step" in order|season|all) ;; *) echo "ERROR: module must be order, season, or all" >&2; exit 1 ;; esac
case "$phenome_source" in auto|phewas|fod) ;; *) echo "ERROR: --phenome-source must be auto, phewas, or fod" >&2; exit 1 ;; esac
case "$primary_model" in m1|m2|m3) ;; *) echo "ERROR: --primary-model must be m1, m2, or m3" >&2; exit 1 ;; esac
case "$within_mode" in auto|0|1) ;; *) echo "ERROR: --within must be auto, 0, or 1" >&2; exit 1 ;; esac

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 🚩 Paths and shared inputs (same convention as phe.sh)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
dir0=${DIR0:-/mnt/d}
phedir=${PHEDIR:-$dir0/data/ukb/phe}
script_dir=${SCRIPT_DIR:-$dir0/scripts}
gen_dir=${GEN_DIR:-$dir0/data/ukb/gen}
this_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
project_name="$(basename "$this_dir")"
task_name="$(basename "${BASH_SOURCE[0]}" .sh)"
outdir=${UKB_OUT:-$dir0/analysis/$project_name}
birth_R=${BIRTH_R:-$this_dir/f/birth.R}

[[ -s "$birth_R" ]] || { echo "ERROR: cannot find $birth_R" >&2; exit 1; }
[[ -s "$phedir/Rdata/all.rds" ]] || { echo "ERROR: cannot find $phedir/Rdata/all.rds" >&2; exit 1; }
mkdir -p "$outdir"

if [[ "$step" == "all" ]]; then
    run_steps="order,season"
else
    run_steps="$step"
fi

echo "========== BIRTH PIPELINE =========="
echo "STEP=$step"
echo "PHEDIR=$phedir"
echo "GEN_DIR=$gen_dir"
echo "OUT=$outdir"
echo "PHENOME_SOURCE=$phenome_source MIN_CASES=$min_cases SIB_CAP=$sib_cap N_JOBS=$jobs PRIMARY_MODEL=$primary_model"
echo

DIR0="$dir0" PHEDIR="$phedir" SCRIPT_DIR="$script_dir" GEN_DIR="$gen_dir" UKB_OUT="$outdir" \
RUN_STEPS="$run_steps" PHENOME_SOURCE="$phenome_source" MIN_CASES="$min_cases" \
MIN_CASES_ORDER="${MIN_CASES_ORDER:-$min_cases}" MIN_CASES_SEASON="${MIN_CASES_SEASON:-$min_cases}" \
MIN_DISCORDANT_WITHIN="${MIN_DISCORDANT_WITHIN:-100}" SIB_CAP="$sib_cap" N_JOBS="$jobs" \
PRIMARY_MODEL="$primary_model" MONTHLY_TOP="$monthly_top" STRATA_TOP="$strata_top" RUN_WITHIN="$within_mode" \
PHENO_MAP="${PHENO_MAP:-}" Rscript "$birth_R"
