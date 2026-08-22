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
  -h, --help              Show this help message.

Environment:
  DIR0=/mnt/d
  PHEDIR=$DIR0/data/ukb/phe
  UKB_OUT=$DIR0/analysis/ukb/assoc
  SCRIPT_DIR=$DIR0/scripts
  PHEWAS_CORES=1

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

if [[ $# -eq 0 ]]; then
    for module in forest phewas medwas; do
        "$0" "$module"
    done
    exit 0
fi
case "$1" in -h|--help) usage; exit 0;; esac

# When METHOD is omitted, apply the same options to forest, phewas, and medwas.
if [[ "$1" == --* ]]; then
    all_options=("$@")
    for module in forest phewas medwas; do
        "$0" "$module" "${all_options[@]}"
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
        "$0" "$method" "$action" "${run_args[@]}"
    done
    exit 0
fi

this_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
dir0=${DIR0:-/mnt/d}
project_name="$(basename "$this_dir")"
task_name="$(basename "${BASH_SOURCE[0]}" .sh)"
outdir=${UKB_OUT:-$dir0/analysis/$project_name/$task_name}
assoc_R=${ASSOC_R:-$this_dir/f/assoc.R}
command -v Rscript >/dev/null 2>&1 || { echo "ERROR: Rscript not found" >&2; exit 1; }
[[ -f "$assoc_R" ]] || { echo "ERROR: f/assoc.R not found: $assoc_R" >&2; exit 1; }

echo "========== RUN ASSOC: METHOD=$method ACTION=$action =========="
DIR0="$dir0" UKB_OUT="$outdir" Rscript "$assoc_R" --method "$method" --action "$action" "$@"
echo "========== DONE ASSOC: METHOD=$method ACTION=$action =========="
