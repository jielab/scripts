#!/usr/bin/env bash

# Display standalone help before loading runtimes or scheduling work.
if [[ ${BASH_SOURCE[0]} == "$0" && ( ${1:-} == -h || ${1:-} == --help ) ]]; then
    cat <<'PUBLIC_ENTRY_HELP'
UKB age/date-of-disease analysis

Modules:
    date-check  Date consistency, recorded burden and clustering [default]

Main options:
    --Y CSV                    Traits from common/*.lst
    --phenotype-dir DIR         UKB Rdata/rap/common root [DIR_UKB]
    --output-dir DIR            Results [ANALYSIS_ROOT/ukb/aod/Y]
    --mode history|baseline     Feature window [history]
    --lookback-years N          Baseline lookback [5]
    --no-linked-records         Use saved AOD and summary evidence only
    --coverage FILE             Source observation intervals
    --input-dir DIR             Reuse an extracted bundle
    --threads N --slots N --memory-gb N --queue NAME --walltime HHH:MM
    --write-cmd                 Write commands without running
    --wait TRUE|FALSE

Run after phe.sh. Reads existing data without modifying Rdata/.
Report: date-check.html includes per-variable input min/max dates.
aod.tsv.gz contains participant identifiers.
Hospital summary AOD is first-per-code, not a complete admission history.

Examples:
    ./aod.sh date-check --Y cvd_cad --write-cmd
    ./aod.sh date-check --Y cvd_cad
    ./aod.sh --Y cvd_cad,cvd_afib,cvd_hfail
    ./aod.sh --Y cvd_cad --mode baseline --lookback-years 5
PUBLIC_ENTRY_HELP
    exit 0
fi
# Disease-date consistency, recorded-event burden and trajectory clustering.
# Phenotype generation remains in phe.sh; this entry only reads its outputs.
this_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)"
source "$this_dir/f/jobs.sh"
ukb_job_init "$@" || exit $?
set -- "${UKB_ARGS[@]}"
source "$shared_dir/cli_dirs.sh"
script_directory_args "$this_dir/aod.sh" "$@"
set -- "${SCRIPT_DIRECTORY_ARGS[@]}"
# Optional explicit spelling of this complete workflow, retained in help so the
# former date-check entry is easy to find. Workers use the canonical arguments.
[[ ${1:-} != date-check ]] || shift
worker_args=("$@")
set -euo pipefail

usage() {
    bash "${BASH_SOURCE[0]}" --help
}

Y=""; phedir=""; common_dir=""; outdir=""; output_given=0
input_dir=""; input_rds=""; coverage=""; as_of=""; official_maps=""
official_maps_given=0
gp_records=auto; hesin=auto; hesin_diag=auto; hesin_oper=auto
mode=history; lookback=5; landmark=5; episode_gap=30; k=auto
seed=20261111; min_cluster_n=20; cluster_sample=5000; stability_repeats=8
demo=0; keep_work=0; rscript=$UKB_RSCRIPT

while (($#)); do
    case "$1" in
        --Y) Y=${2:?ERROR: --Y requires NAME}; shift 2 ;;
        --phenotype-dir|--dir-phenotype) phedir=${2:?ERROR: --phenotype-dir requires DIR}; shift 2 ;;
        --common-dir|--dir-common) common_dir=${2:?ERROR: --common-dir requires DIR}; shift 2 ;;
        --output-dir|--dir-output) outdir=${2:?ERROR: --output-dir requires DIR}; output_given=1; shift 2 ;;
        --input-dir|--dir-input) input_dir=${2:?ERROR: --input-dir requires DIR}; shift 2 ;;
        --input-rds) input_rds=${2:?ERROR: --input-rds requires FILE}; shift 2 ;;
        --gp-records) gp_records=${2:?ERROR: --gp-records requires auto, none, or FILE}; shift 2 ;;
        --hesin) hesin=${2:?ERROR: --hesin requires auto, none, or FILE}; shift 2 ;;
        --hesin-diag) hesin_diag=${2:?ERROR: --hesin-diag requires auto, none, or FILE}; shift 2 ;;
        --hesin-oper) hesin_oper=${2:?ERROR: --hesin-oper requires auto, none, or FILE}; shift 2 ;;
        --no-linked-records) gp_records=none; hesin=none; hesin_diag=none; hesin_oper=none; shift ;;
        --official-maps) official_maps=${2:?ERROR: --official-maps requires FILE or none}; official_maps_given=1; shift 2 ;;
        --coverage) coverage=${2:?ERROR: --coverage requires FILE}; shift 2 ;;
        --as-of) as_of=${2:?ERROR: --as-of requires YYYY-MM-DD}; shift 2 ;;
        --mode) mode=${2:?ERROR: --mode requires history or baseline}; shift 2 ;;
        --lookback-years) lookback=${2:?ERROR: --lookback-years requires N}; shift 2 ;;
        --landmark-years) landmark=${2:?ERROR: --landmark-years requires N}; shift 2 ;;
        --episode-gap) episode_gap=${2:?ERROR: --episode-gap requires DAYS}; shift 2 ;;
        --k) k=${2:?ERROR: --k requires auto or 2..8}; shift 2 ;;
        --seed) seed=${2:?ERROR: --seed requires N}; shift 2 ;;
        --min-cluster-n) min_cluster_n=${2:?ERROR: --min-cluster-n requires N}; shift 2 ;;
        --cluster-sample) cluster_sample=${2:?ERROR: --cluster-sample requires N}; shift 2 ;;
        --stability-repeats) stability_repeats=${2:?ERROR: --stability-repeats requires N}; shift 2 ;;
        --rscript) rscript=${2:?ERROR: --rscript requires FILE}; shift 2 ;;
        --python) UKB_PYTHON=${2:?ERROR: --python requires FILE}; shift 2 ;;
        --demo) demo=1; shift ;;
        --keep-work) keep_work=1; shift ;;
        -h|--help) usage; exit 0 ;;
        *) echo "ERROR: unknown aod.sh option: $1" >&2; usage >&2; exit 2 ;;
    esac
done

[[ $Y =~ ^[A-Za-z0-9_]+(,[A-Za-z0-9_]+)*$ ]] || { echo 'ERROR: --Y NAME[,NAME...] is required' >&2; exit 2; }
[[ $mode == history || $mode == baseline ]] || { echo 'ERROR: --mode must be history or baseline' >&2; exit 2; }
[[ $k == auto || $k =~ ^[2-8]$ ]] || { echo 'ERROR: --k must be auto or 2..8' >&2; exit 2; }
for value in "$lookback" "$landmark"; do
    [[ $value =~ ^([1-9][0-9]*([.][0-9]+)?|0[.][0-9]*[1-9][0-9]*)$ ]] || { echo 'ERROR: year windows must be positive numbers' >&2; exit 2; }
done
for value in "$episode_gap" "$seed"; do
    [[ $value =~ ^(0|[1-9][0-9]*)$ ]] || { echo 'ERROR: episode gap and seed must be non-negative integers' >&2; exit 2; }
done
for value in "$min_cluster_n" "$cluster_sample" "$stability_repeats"; do
    [[ $value =~ ^[1-9][0-9]*$ ]] || { echo 'ERROR: cluster sizes and repeats must be positive integers' >&2; exit 2; }
done
[[ -z $as_of || $as_of =~ ^[0-9]{4}-[0-9]{2}-[0-9]{2}$ ]] || { echo 'ERROR: --as-of must be YYYY-MM-DD' >&2; exit 2; }
if [[ -z $as_of && -z $input_dir ]] && (( demo == 0 )); then
    # Freeze fresh-input runs. Extracted bundles retain their own export cutoff;
    # the backend validates it and rejects expansion beyond available evidence.
    # Neither boundary is a clinical release cutoff or a coverage assumption.
    as_of=$(date +%F)
    worker_args+=(--as-of "$as_of")
fi
[[ -z $as_of || $(date -d "$as_of" +%F 2>/dev/null) == "$as_of" ]] || { echo 'ERROR: --as-of is not a valid calendar date' >&2; exit 2; }
(( demo == 0 )) || [[ -z $input_dir && -z $input_rds && -z $coverage ]] || { echo 'ERROR: --demo cannot be combined with real input/coverage overrides' >&2; exit 2; }
[[ -z $input_dir || -z $input_rds ]] || { echo 'ERROR: use either --input-dir or --input-rds' >&2; exit 2; }
[[ -z $input_dir || $Y != *,* ]] || { echo 'ERROR: --input-dir accepts one trait only' >&2; exit 2; }
phedir=${phedir:-$DIR_UKB}; common_dir=${common_dir:-$phedir/common}
official_maps=${official_maps:-$common_dir/first_occurrence_maps.tsv}
if [[ $official_maps != none && ! -f $official_maps ]]; then
    if (( official_maps_given )); then
        echo "ERROR: explicit official mapping file does not exist: $official_maps" >&2; exit 2
    fi
    official_maps=none
fi
outdir=${outdir:-$ANALYSIS_ROOT/ukb/aod}
export UKB_PYTHON
ukb_dispatch aod "$outdir" "$this_dir/aod.sh" "${worker_args[@]}"
ukb_worker_runtime

run_rscript() {
    local executable
    executable=$(command -v "$rscript") || { echo "ERROR: Rscript not found: $rscript" >&2; return 2; }
    if [[ $(readlink -f "$executable") == /usr/bin/Rscript ]]; then
        env -u PROJ_LIB -u PROJ_DATA -u GDAL_DATA "$rscript" "$@"
    else
        "$rscript" "$@"
    fi
}
IFS=',' read -r -a traits <<< "$Y"
root=$outdir
declare -A seen=()
for Y in "${traits[@]}"; do
    [[ -z ${seen[$Y]:-} ]] || continue
    seen[$Y]=1
    (
        outdir=$root
        if (( output_given == 0 || ${#traits[@]} > 1 )); then outdir=$outdir/$Y; fi
        mkdir -p -- "$outdir"
        outdir="$(cd -- "$outdir" && pwd -P)"
        # A private staging directory prevents failed analysis from replacing a report.
        work=$(mktemp -d "$outdir/.aod.XXXXXX")
        cleanup() {
            local status=$?
            trap - EXIT
            if (( keep_work )); then
                printf 'AOD staging retained: %s\n' "$work"
            elif [[ $work == "$outdir"/.aod.* ]]; then
                rm -rf -- "$work"
            fi
            exit "$status"
        }
        trap cleanup EXIT
        trap 'exit 130' INT
        trap 'exit 143' TERM
        printf 'AOD: Y=%s; mode=%s; output=%s\n' "$Y" "$mode" "$outdir"
        analysis_args=(--input-dir "$work" --output-dir "$work" --Y "$Y"
            --mode "$mode" --lookback-years "$lookback" --landmark-years "$landmark"
            --episode-gap "$episode_gap" --k "$k" --seed "$seed"
            --min-cluster-n "$min_cluster_n" --cluster-sample "$cluster_sample"
            --stability-repeats "$stability_repeats")
        [[ -z $coverage ]] || analysis_args+=(--coverage "$coverage")
        [[ -z $as_of ]] || analysis_args+=(--as-of "$as_of")
        if (( demo )); then
            analysis_args+=(--demo)
        elif [[ -n $input_dir ]]; then
            for item in patients.tsv.gz events.tsv.gz metadata.tsv; do
                [[ -s $input_dir/$item ]] || { echo "ERROR: missing extracted input: $input_dir/$item" >&2; exit 2; }
                cp -- "$input_dir/$item" "$work/$item"
            done
        else
            [[ -d $phedir/Rdata || -n $input_rds ]] || { echo "ERROR: missing $phedir/Rdata" >&2; exit 2; }
            command -v "$rscript" >/dev/null || { echo "ERROR: Rscript not found: $rscript" >&2; exit 2; }
            linked_args=(linked --Y "$Y" --phenotype-dir "$phedir" --common-dir "$common_dir"
                --output-dir "$work" --gp-records "$gp_records" --hesin "$hesin"
                --hesin-diag "$hesin_diag" --hesin-oper "$hesin_oper" --official-maps "$official_maps")
            [[ -z $input_rds ]] || linked_args+=(--allow-undefined-trait)
            "$UKB_PYTHON" "$this_dir/f/date_check_io.py" "${linked_args[@]}"
            r_args=(--Y "$Y" --phenotype-dir "$phedir" --output-dir "$work")
            [[ -z $as_of ]] || r_args+=(--as-of "$as_of")
            [[ -z $input_rds ]] || r_args+=(--input-rds "$input_rds")
            run_rscript "$this_dir/f/aod.R" "${r_args[@]}"
        fi
        "$UKB_PYTHON" "$this_dir/f/aod.py" "${analysis_args[@]}"
        "$UKB_PYTHON" "$this_dir/f/aod_report.py" --input "$work/summary.json" --output "$work/date-check.html"
        # Publish only after every requested result exists; publish the report last.
        for item in aod.tsv.gz clusters.tsv model.json date-check.html; do
            [[ -s $work/$item ]] || { echo "ERROR: analysis did not produce $item" >&2; exit 1; }
        done
        for item in aod.tsv.gz clusters.tsv model.json date-check.html; do
            mv -f -- "$work/$item" "$outdir/$item"
        done
        printf 'AOD report: %s/date-check.html\n' "$outdir"
    )
done
