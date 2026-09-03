#!/usr/bin/env bash
set -euo pipefail

usage() {
    cat <<'EOF'
Usage:
  ./assoc.sh                         Run forest, phewas, then medwas
  ./assoc.sh [OPTIONS]               Run all three modules with the options
  ./assoc.sh METHOD [ACTION] [OPTIONS]

Methods and actions:
  forest    [save_rds|plot]
                          save_rds: run assoc_reg() and save the result
                          plot: read the saved result and draw forest plots
                          No action: run save_rds, then plot
  phewas    [save_rds|plot]
                          save_rds: run PheWAS and save results/checkpoints
                          plot: read the saved result and draw PheWAS plots
                          No action: use existing RDS, or run save_rds, then plot
  medwas    [save_rds|plot]
                          save_rds: regress each selected X against img/met/prot/bbc traits
                          plot: draw four volcano panels per X
                          No action: run save_rds, then plot

Options:
  --Xs LIST               Comma-separated exposures/gen variables.
                          Default: all columns in gen.rds except the ID column.
  --varX LIST             Comma-separated covariates.
                          Default: c(vars.basic, "le8.sco").
  --mode MODE             Genotype model(s): add, dom, res, or a comma-separated
                          list such as res,dom,add. Requires explicit
                          --Xs genotype variables. Dosage is retained for add
                          and rounded to hard calls for dom/res.
  --Ys, --Y LIST          Comma-separated outcomes for forest.
                          Default: all fod_icd10_* outcomes in all.rds.
  --type TYPE             assoc_reg model type: bt, qt, t2e, t2e.tdc, ordinal.
                          Default: t2e.
  --scale-X, --scale_X BOOL
                          Scale numeric exposures: TRUE or FALSE. Default: TRUE.
  --out-prefix NAME       Result prefix. Default: gen.
  --no-resume             Recompute completed PheWAS checkpoints.
  --data-root DIR         Data/analysis root. Default: /mnt/d.
  --phenotype-dir DIR     UKB phenotype directory.
  --output-dir DIR        Association output directory.
  --script-dir DIR        Shared scripts root.
  --r-script FILE         R workflow entry point.
  -h, --help              Show this help message.

Usage examples:
  ./assoc.sh forest --Xs abo.AB.rs8176746_G,apoe.rs7412_C --Ys ihd,stroke --type t2e
  ./assoc.sh phewas save_rds --Xs alc.aldh2.rs671_G,fto.rs9939609_T --varX age,sex,PC1,PC2
  ./assoc.sh phewas plot --Xs alc.aldh2.rs671_G,fto.rs9939609_T
  ./assoc.sh medwas --Xs alc.aldh2.rs671_G,fto.rs9939609_T --varX age,sex,PC1,PC2
  ./assoc.sh phewas --Xs smell.rs1953558_C --mode res
  ./assoc.sh medwas --Xs smell.rs1953558_C --mode res --varX age,sex,tdi,PC1,PC2
  ./assoc.sh --Xs smell.rs1953558_C --mode res,dom,add --varX age,sex,tdi,PC1,PC2
EOF
}

this_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
dir0=/mnt/d
phedir=""
outdir=""
script_dir=""
assoc_R=""
filtered_args=()

while (( $# )); do
    case "$1" in
        --data-root) dir0=${2:?ERROR: --data-root requires DIR}; shift 2 ;;
        --phenotype-dir) phedir=${2:?ERROR: --phenotype-dir requires DIR}; shift 2 ;;
        --output-dir) outdir=${2:?ERROR: --output-dir requires DIR}; shift 2 ;;
        --script-dir) script_dir=${2:?ERROR: --script-dir requires DIR}; shift 2 ;;
        --r-script) assoc_R=${2:?ERROR: --r-script requires FILE}; shift 2 ;;
        *) filtered_args+=("$1"); shift ;;
    esac
done
set -- "${filtered_args[@]}"
phedir=${phedir:-$dir0/data/ukb/phe}
outdir=${outdir:-$dir0/analysis/ukb/assoc}
script_dir=${script_dir:-$dir0/scripts}
assoc_R=${assoc_R:-$this_dir/f/assoc.R}
path_args=(--data-root "$dir0" --phenotype-dir "$phedir" --output-dir "$outdir" --script-dir "$script_dir" --r-script "$assoc_R")

if [[ $# -eq 0 ]]; then
    for module in forest phewas medwas; do
        "$0" "$module" "${path_args[@]}"
    done
    exit 0
fi
case "$1" in -h|--help) usage; exit 0;; esac

# When METHOD is omitted, apply the same options to forest, phewas, and medwas.
if [[ "$1" == --* ]]; then
    all_options=("$@")
    for module in forest phewas medwas; do
        "$0" "$module" "${all_options[@]}" "${path_args[@]}"
    done
    exit 0
fi

method="$1"
shift
case "$method" in
    forest|phewas|medwas) default_action=all ;;
    *) echo "ERROR: unknown METHOD: $method" >&2; usage >&2; exit 1 ;;
esac

action="$default_action"
if [[ $# -gt 0 && "$1" != --* ]]; then action="$1"; shift; fi
case "$method:$action" in
    forest:all|forest:save_rds|forest:plot|phewas:all|phewas:save_rds|phewas:plot|medwas:all|medwas:save_rds|medwas:plot) ;;
    *) echo "ERROR: unsupported METHOD/ACTION: $method $action" >&2; usage >&2; exit 1 ;;
esac

# Expand a comma-separated --mode list into independent runs. Keeping each mode
# in its own R process also keeps checkpoints and result files isolated.
mode_value=""
args=("$@")
for ((i=0; i<${#args[@]}; i++)); do
    if [[ "${args[i]}" == "--mode" ]]; then
        if (( i + 1 >= ${#args[@]} )); then
            echo "ERROR: --mode requires a value" >&2
            exit 1
        fi
        mode_value="${args[i+1]}"
        break
    fi
done
if [[ "$mode_value" == *,* ]]; then
    IFS=',' read -ra modes <<< "$mode_value"
    for mode_one in "${modes[@]}"; do
        mode_one="${mode_one//[[:space:]]/}"
        [[ -n "$mode_one" ]] || { echo "ERROR: --mode contains an empty entry: $mode_value" >&2; exit 1; }
        run_args=("${args[@]}")
        for ((j=0; j<${#run_args[@]}; j++)); do
            if [[ "${run_args[j]}" == "--mode" ]]; then run_args[j+1]="$mode_one"; break; fi
        done
        "$0" "$method" "$action" "${run_args[@]}" "${path_args[@]}"
    done
    exit 0
fi

command -v Rscript >/dev/null 2>&1 || { echo "ERROR: Rscript not found" >&2; exit 1; }
[[ -f "$assoc_R" ]] || { echo "ERROR: f/assoc.R not found: $assoc_R" >&2; exit 1; }

echo "========== RUN ASSOC: METHOD=$method ACTION=$action =========="
DIR0="$dir0" PHEDIR="$phedir" UKB_OUT="$outdir" SCRIPT_DIR="$script_dir" \
    Rscript "$assoc_R" --method "$method" --action "$action" "$@"
echo "========== DONE ASSOC: METHOD=$method ACTION=$action =========="
