#!/usr/bin/env bash

# Display standalone help before loading runtimes or scheduling work.
if [[ ${BASH_SOURCE[0]} == "$0" && ( ${1:-} == -h || ${1:-} == --help ) ]]; then
    cat <<'PUBLIC_ENTRY_HELP'
./aod.sh  # 默认 date-check 全流程
./aod.sh date-check  # 单一分析模块
./aod.sh --mode history --lookback-years 5 --landmark-years 2

LE8 age/date-of-disease analysis

Modules:
    date-check  Recorded course, clusters, ADuLT and date consistency [default]

Main options:
    --Y CSV                    Traits from common/*.lst
    --phenotype-dir DIR         YGB Rdata/rap/common root [DIR_YGB]
    --output-dir DIR            Results [ANALYSIS_ROOT/le8/Y/aod]
    --mode history|baseline     Feature window [history]
    --lookback-years N          Baseline lookback [5]
    --landmark-years N          Onset-to-landmark window [2]
    --trajectory-source SOURCE  hospital|gp|combined [hospital]
    --apply-model FILE          Assign clusters with a frozen model
    --validate                 Test post-landmark mortality prediction (SciPy)
    --no-linked-records         Use saved AOD and summary evidence only
    --coverage FILE             Source observation intervals
    --input-dir DIR             Reuse an extracted bundle
    --threads N --slots N --queue NAME
    --write-cmd                 Write commands without running
    --wait TRUE|FALSE

Run after phe.sh. Reads existing data without modifying Rdata/.
Clinical follow-up: 0f/phenotype.R::date_follow_end; input audit dates are separate.
Report: main/reports/date-check.html leads with findings, course profiles and validation.
tmp/aod.tsv.gz contains participant identifiers.
Hospital summary AOD is first-per-code, not a complete admission history.

Examples:
    ./aod.sh date-check --Y cvd_cad --write-cmd
    ./aod.sh date-check --Y cvd_cad
    ./aod.sh --Y cvd_cad,cvd_afib,cvd_hfail
    ./aod.sh --Y cvd_afib --validate
    ./aod.sh --Y cvd_afib --landmark-years 5 --output-dir /tmp/afib-5y
    ./aod.sh --Y cvd_cad --mode baseline --lookback-years 5
PUBLIC_ENTRY_HELP
    exit 0
fi
# Disease-date consistency, recorded-event burden and trajectory clustering.
# Phenotype generation remains in phe.sh; this entry only reads its outputs.
entry_file="$(realpath -- "${BASH_SOURCE[0]}")"
this_dir="$(cd -- "$(dirname -- "$entry_file")/.." && pwd -P)"
source "$this_dir/../0f/phenotype_jobs.sh"
ygb_job_init "$@" || exit $?
set -- "${YGB_ARGS[@]}"
source "$shared_dir/cli_dirs.sh"
script_directory_args "$entry_file" "$@"
set -- "${SCRIPT_DIRECTORY_ARGS[@]}"
# Explicit analysis-module name. Workers use the canonical arguments.
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
mode=history; lookback=5; landmark=2; episode_gap=30; k=auto
trajectory_source=hospital; apply_model=""; validate=0
seed=20261111; min_cluster_n=20; cluster_sample=5000; stability_repeats=8
demo=0; keep_work=0; rscript=$YGB_RSCRIPT

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
        --trajectory-source) trajectory_source=${2:?ERROR: --trajectory-source requires hospital, gp, or combined}; shift 2 ;;
        --apply-model) apply_model=${2:?ERROR: --apply-model requires FILE}; shift 2 ;;
        --validate) validate=1; shift ;;
        --episode-gap) episode_gap=${2:?ERROR: --episode-gap requires DAYS}; shift 2 ;;
        --k) k=${2:?ERROR: --k requires auto or 2..8}; shift 2 ;;
        --seed) seed=${2:?ERROR: --seed requires N}; shift 2 ;;
        --min-cluster-n) min_cluster_n=${2:?ERROR: --min-cluster-n requires N}; shift 2 ;;
        --cluster-sample) cluster_sample=${2:?ERROR: --cluster-sample requires N}; shift 2 ;;
        --stability-repeats) stability_repeats=${2:?ERROR: --stability-repeats requires N}; shift 2 ;;
        --rscript) rscript=${2:?ERROR: --rscript requires FILE}; shift 2 ;;
        --python) YGB_PYTHON=${2:?ERROR: --python requires FILE}; shift 2 ;;
        --demo) demo=1; shift ;;
        --keep-work) keep_work=1; shift ;;
        -h|--help) usage; exit 0 ;;
        *) echo "ERROR: unknown aod.sh option: $1" >&2; usage >&2; exit 2 ;;
    esac
done

[[ $Y =~ ^[A-Za-z0-9_]+(,[A-Za-z0-9_]+)*$ ]] || { echo 'ERROR: --Y NAME[,NAME...] is required' >&2; exit 2; }
[[ $mode == history || $mode == baseline ]] || { echo 'ERROR: --mode must be history or baseline' >&2; exit 2; }
[[ $trajectory_source == hospital || $trajectory_source == gp || $trajectory_source == combined ]] || { echo 'ERROR: --trajectory-source must be hospital, gp, or combined' >&2; exit 2; }
[[ -z $apply_model || -s $apply_model ]] || { echo 'ERROR: --apply-model file is missing or empty' >&2; exit 2; }
[[ -z $apply_model || $Y != *,* ]] || { echo 'ERROR: --apply-model accepts one disease at a time' >&2; exit 2; }
(( validate == 0 )) || [[ -z $apply_model ]] || { echo 'ERROR: validate a train/test fit; --apply-model produces assignments only' >&2; exit 2; }
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
    # Freeze only the raw-input audit snapshot; clinical observation is capped
    # separately by the shared date_follow_end inside aod.R/aod.py.
    as_of=$(date +%F)
    worker_args+=(--as-of "$as_of")
fi
[[ -z $as_of || $(date -d "$as_of" +%F 2>/dev/null) == "$as_of" ]] || { echo 'ERROR: --as-of is not a valid calendar date' >&2; exit 2; }
(( demo == 0 )) || [[ -z $input_dir && -z $input_rds && -z $coverage ]] || { echo 'ERROR: --demo cannot be combined with real input/coverage overrides' >&2; exit 2; }
[[ -z $input_dir || -z $input_rds ]] || { echo 'ERROR: use either --input-dir or --input-rds' >&2; exit 2; }
[[ -z $input_dir || $Y != *,* ]] || { echo 'ERROR: --input-dir accepts one trait only' >&2; exit 2; }
phedir=${phedir:-$DIR_YGB}; common_dir=${common_dir:-$phedir/common}
official_maps=${official_maps:-$common_dir/first_occurrence_maps.tsv}
if [[ $official_maps != none && ! -f $official_maps ]]; then
    if (( official_maps_given )); then
        echo "ERROR: explicit official mapping file does not exist: $official_maps" >&2; exit 2
    fi
    official_maps=none
fi
outdir=${outdir:-$ANALYSIS_ROOT/le8}
export YGB_PYTHON
if (( YGB_EXECUTE == 0 )); then
    exec "$YGB_PYTHON" "$(dirname -- "$entry_file")/le8_scheduler.py" aod \
        "$outdir" "$Y" "$YGB_QUEUE" "$YGB_SLOTS" "$YGB_WRITE_ONLY" "$entry_file" -- \
        "${worker_args[@]}" --threads "$YGB_THREADS" --rscript "$rscript" --python "$YGB_PYTHON"
fi
ygb_worker_runtime

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
        if (( output_given == 0 )); then outdir=$root/$Y/aod; elif (( ${#traits[@]} > 1 )); then outdir=$root/$Y; fi
        mkdir -p -- "$outdir"
        outdir="$(cd -- "$outdir" && pwd -P)"
        # A private staging directory prevents failed analysis from replacing a report.
        mkdir -p "$outdir/tmp"
        work=$(mktemp -d "$outdir/tmp/.aod.XXXXXX")
        cleanup() {
            local status=$?
            trap - EXIT
            if (( keep_work )); then
                printf 'AOD staging retained: %s\n' "$work"
            elif [[ $work == "$outdir"/tmp/.aod.* ]]; then
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
            --trajectory-source "$trajectory_source"
            --episode-gap "$episode_gap" --k "$k" --seed "$seed"
            --min-cluster-n "$min_cluster_n" --cluster-sample "$cluster_sample"
            --stability-repeats "$stability_repeats")
        [[ -z $coverage ]] || analysis_args+=(--coverage "$coverage")
        [[ -z $apply_model ]] || analysis_args+=(--apply-model "$apply_model")
        [[ -z $as_of ]] || analysis_args+=(--as-of "$as_of")
        if (( demo )); then
            analysis_args+=(--demo)
        elif [[ -n $input_dir ]]; then
            for item in patients.tsv.gz events.tsv.gz metadata.tsv; do
                [[ -s $input_dir/$item ]] || { echo "ERROR: missing extracted input: $input_dir/$item" >&2; exit 2; }
                cp -- "$input_dir/$item" "$work/$item"
            done
            if [[ -z $coverage && -s $input_dir/coverage.tsv ]]; then
                cp -- "$input_dir/coverage.tsv" "$work/coverage.tsv"
                analysis_args+=(--coverage "$work/coverage.tsv")
            fi
        else
            [[ -d $phedir/Rdata || -n $input_rds ]] || { echo "ERROR: missing $phedir/Rdata" >&2; exit 2; }
            command -v "$rscript" >/dev/null || { echo "ERROR: Rscript not found: $rscript" >&2; exit 2; }
            linked_args=(linked --Y "$Y" --phenotype-dir "$phedir" --common-dir "$common_dir"
                --output-dir "$work" --gp-records "$gp_records" --hesin "$hesin"
                --hesin-diag "$hesin_diag" --hesin-oper "$hesin_oper" --official-maps "$official_maps")
            [[ -z $input_rds ]] || linked_args+=(--allow-undefined-trait)
            "$YGB_PYTHON" "$shared_dir/phenotype_date_io.py" "${linked_args[@]}"
            r_args=(--Y "$Y" --phenotype-dir "$phedir" --output-dir "$work")
            [[ -z $as_of ]] || r_args+=(--as-of "$as_of")
            [[ -z $input_rds ]] || r_args+=(--input-rds "$input_rds")
            run_rscript "$this_dir/f/aod.R" "${r_args[@]}"
        fi
        "$YGB_PYTHON" "$this_dir/f/aod.py" "${analysis_args[@]}"
        report_args=(--input "$work/summary.json" --output "$work/date-check.html")
        if (( validate )); then
            validation_args=(--input "$work/aod.tsv.gz" --Y "$Y" --output-dir "$work")
            (( demo == 0 )) || validation_args+=(--synthetic)
            "$YGB_PYTHON" "$this_dir/f/aod_validate.py" "${validation_args[@]}"
            report_args+=(--validation "$work/validation.json")
        fi
        "$YGB_PYTHON" "$this_dir/f/aod_report.py" "${report_args[@]}"
        # Publish only after every requested result exists; publish the report last.
        for item in aod.tsv.gz clusters.tsv model.json date-check.html; do
            [[ -s $work/$item ]] || { echo "ERROR: analysis did not produce $item" >&2; exit 1; }
        done
        if [[ -n ${LE8_DX_OUTPUT:-} ]]; then
            mkdir -p "$outdir/logs"
            exec {publish_lock}>"$outdir/logs/publish.lock"
            flock "$publish_lock"
            run_rscript "$(dirname -- "$entry_file")/dx_export.R" "$phedir/Rdata/all.rds" "$Y" "$work" "$shared_dir"
            "$YGB_PYTHON" "$(dirname -- "$entry_file")/dx_dates.py" "$work" "$outdir" "$LE8_DX_OUTPUT" "$Y"
        fi
        for item in aod.tsv.gz clusters.tsv model.json date-check.html; do
            case "$item" in
                date-check.html) bucket=main/reports ;;
                aod.tsv.gz) bucket=tmp ;;
                *) bucket=main/results ;;
            esac
            mkdir -p "$outdir/$bucket"
            mv -f -- "$work/$item" "$outdir/$bucket/$item"
        done
        if (( validate )); then
            for item in validation.json validation.tsv; do
                [[ ! -s $work/$item ]] || mv -f -- "$work/$item" "$outdir/main/results/$item"
            done
        fi
        printf 'AOD report: %s/main/reports/date-check.html\n' "$outdir"
    )
done

