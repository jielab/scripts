#!/usr/bin/env bash
set -euo pipefail

# Internal entry used by systemd-run before each analysis process.
run_memory_scope() {
    local scope ram swap cgroup base
    fail() { printf 'ERROR: DOSE memory scope verification failed: %s\n' "$*" >&2; exit 2; }
    (( $# >= 4 )) || fail "expected scope, RAM max, swap max, and command"
    scope=$1 ram=$2 swap=$3
    shift 3
    [[ $scope == dose-*.scope && $scope != */* ]] || fail "invalid scope name"
    [[ $ram =~ ^[1-9][0-9]*$ && $swap =~ ^[0-9]+$ ]] || fail "invalid limits"
    cgroup=$(awk -F: '$1 == "0" {print $3}' /proc/self/cgroup)
    [[ $cgroup == */"$scope" ]] || fail "analysis is outside the requested scope"
    base="/sys/fs/cgroup$cgroup"
    check_limit() {
        [[ -r $base/$1 && $(cat "$base/$1") == "$2" ]] || fail "$1 does not match $2"
    }
    check_limit memory.max "$ram"
    check_limit memory.swap.max "$swap"
    check_limit memory.oom.group 1
    printf 'Memory guard verified: scope=%s RAM=%s bytes swap=%s bytes (whole process tree)\n' "$scope" "$ram" "$swap" >&2
    exec "$@"
}
if [[ ${1:-} == --internal-memory-scope ]]; then
    shift
    run_memory_scope "$@"
    exit 0
fi

SCRIPT_DIR=$(dirname -- "$(realpath -- "${BASH_SOURCE[0]}")")

usage() {
    cat <<'HELP'
DOSE pipeline: ./dose.sh [MODULE] [OPTIONS]

Modules (default: all):
    foodwas     s1: dietary patterns, food associations and dose-response plots
                Starts with the two s0 basic preprocessing steps
    ml          s2: machine learning, feature selection and validation
    score       s3: DOSE score, stroke and other outcomes
    mechanism   s4: omics mediation, MRI and structural equation models
    signature   s5: proteomic/metabolomic LASSO signatures and Cox models
    nhanes      NHANES external validation
    all         foodwas -> ml -> score -> mechanism -> signature
                -> s4 summary figure -> nhanes -> UKB/NHANES stroke comparison

Options:
    --steps CSV        Run selected step names within ONE module, in pipeline order
    --list-steps       List steps and their execution order; no analysis
    --dry-run          Print commands and directories; no analysis or output writes
    --replace BOOL     TRUE rebuilds/overwrites results; FALSE checks and skips [FALSE]
    --skip-basic       Skip s0 in foodwas/all (requires prepared data)
    --jobs N           Python model/BLAS threads [4]
    --memory-limit-gb N Whole-step RAM limit, including all workers [16 GiB]
    --memory-swap-gb N  Whole-step swap limit; 0 disables swap [0 GiB]
    --dir-output DIR   Results [/mnt/d/analysis/dose]
    --dir-ukb DIR      UKB input root [/mnt/d/data/ukb]
    --dir-nhanes DIR   NHANES input root [/mnt/d/files/diet; contains cdc/]
    --dir-files DIR    Dictionaries [/mnt/d/files/diet]
    --dir-helper DIR   Shared R helpers [/mnt/d/scripts/0f]
    --rscript FILE     R interpreter [Rscript]
    --python FILE      Python interpreter [DOSE runtime venv, or python3]
    -- ARGS...        Pass script arguments to a single explicitly selected step
    -h, --help         Show help

Aliases: s1=foodwas, s2=ml, s3=score, s4=mechanism, s5=signature.
Standalone modules use existing upstream results; they do not rerun other modules.
The mechanism summary runs at the end of signature, once its inputs are available.
The UKB/NHANES stroke comparison runs after nhanes has generated its results.
Shared functions are in f/0.common.R and f/0.common.py.
Each analysis step runs in a separate systemd user scope with a shared process-tree
memory budget. An out-of-memory event at the hard limit stops that step and the
pipeline. Requires systemd and cgroup v2;
help, listing and dry runs do not require systemd. Concurrent runs have separate budgets.
Run ./dose.sh with no arguments for the entire workflow below.
Completed steps are skipped only when their inputs, code and output files match
the saved checkpoint. Interrupted steps with matching inputs can resume. Changed
or untracked existing results require --replace TRUE; this also rebuilds caches.
Keep .cache (NHANES scoring intermediates) and .checkpoints (verified results and
resume state) after success so later steps and reruns can reuse completed work.

Examples (run these six commands in order for the same complete workflow):
    ./dose.sh foodwas
    ./dose.sh ml
    ./dose.sh score
    ./dose.sh mechanism
    ./dose.sh signature
    ./dose.sh nhanes
HELP
}

die() { printf 'ERROR: %s\n' "$*" >&2; exit 2; }
value_required() { [[ $# -ge 2 && -n $2 && $2 != --* ]] || die "$1 requires a value"; }
MODULE=all
MODULE_SET=0
STEPS=""
DRY_RUN=0
LIST_STEPS=0
SKIP_BASIC=0
REPLACE=FALSE
OUTDIR=${DOSE_OUTDIR:-${DOSE_PROJECT:-/mnt/d/analysis/dose}}
UKBDIR=${DOSE_UKB:-/mnt/d/data/ukb}
NHANESDIR=${NHANES_DATA:-}
FILESDIR=${DOSE_FILES_DIR:-/mnt/d/files/diet}
HELPERDIR=${DOSE_HELPER_DIR:-/mnt/d/scripts/0f}
JOBS=${DOSE_THREADS:-4}
MEMORY_LIMIT_GB=${DOSE_MEMORY_LIMIT_GB:-16}
MEMORY_SWAP_GB=${DOSE_MEMORY_SWAP_GB:-0}
RSCRIPT_BIN=${RSCRIPT_BIN:-Rscript}
RUNTIME_DIR=${DOSE_RUNTIME_DIR:-$HOME/.local/share/dose}
DEFAULT_PYTHON=python3
if [[ -x $RUNTIME_DIR/venv/bin/python ]]; then DEFAULT_PYTHON=$RUNTIME_DIR/venv/bin/python; fi
PYTHON_BIN=${DOSE_PYTHON:-${RETICULATE_PYTHON:-$DEFAULT_PYTHON}}
SCRIPT_ARGS=()

while (($#)); do
    case "$1" in
        --steps) value_required "$@"; STEPS=$2; shift 2 ;;
        --jobs) value_required "$@"; JOBS=$2; shift 2 ;;
        --memory-limit-gb) value_required "$@"; MEMORY_LIMIT_GB=$2; shift 2 ;;
        --memory-swap-gb) value_required "$@"; MEMORY_SWAP_GB=$2; shift 2 ;;
        --dir-output|--output-dir) value_required "$@"; OUTDIR=$2; shift 2 ;;
        --dir-ukb) value_required "$@"; UKBDIR=$2; shift 2 ;;
        --dir-nhanes) value_required "$@"; NHANESDIR=$2; shift 2 ;;
        --dir-files) value_required "$@"; FILESDIR=$2; shift 2 ;;
        --dir-helper) value_required "$@"; HELPERDIR=$2; shift 2 ;;
        --rscript) value_required "$@"; RSCRIPT_BIN=$2; shift 2 ;;
        --python) value_required "$@"; PYTHON_BIN=$2; shift 2 ;;
        --dry-run) DRY_RUN=1; shift ;;
        --replace) value_required "$@"; REPLACE=${2^^}; shift 2 ;;
        --list-steps) LIST_STEPS=1; shift ;;
        --skip-basic) SKIP_BASIC=1; shift ;;
        -h|--help) usage; exit 0 ;;
        --) shift; SCRIPT_ARGS=("$@"); break ;;
        all|foodwas|ml|score|mechanism|signature|nhanes|s[1-5])
            ((MODULE_SET == 0)) || die "Specify only one module"
            MODULE=$1; MODULE_SET=1; shift ;;
        *) die "Unknown argument: $1 (see --help)" ;;
    esac
done
case "$MODULE" in
    s1) MODULE=foodwas ;; s2) MODULE=ml ;;
    s3) MODULE=score ;; s4) MODULE=mechanism ;; s5) MODULE=signature ;;
esac
[[ $JOBS =~ ^[1-9][0-9]*$ ]] || die "--jobs must be a positive integer"
[[ $REPLACE == TRUE || $REPLACE == FALSE ]] || die "--replace must be TRUE or FALSE"
[[ $MEMORY_LIMIT_GB =~ ^[1-9][0-9]{0,4}$ ]] || die "--memory-limit-gb must be a positive integer (GiB)"
[[ $MEMORY_SWAP_GB =~ ^(0|[1-9][0-9]{0,4})$ ]] || die "--memory-swap-gb must be a non-negative integer (GiB)"
MEMORY_MAX=$((MEMORY_LIMIT_GB * 1073741824))
MEMORY_SWAP_MAX=$((MEMORY_SWAP_GB * 1073741824))
[[ -z $STEPS || $MODULE != all ]] || die "--steps requires a single module"
[[ $SKIP_BASIC == 0 || $MODULE == all || $MODULE == foodwas ]] || die "--skip-basic is only for foodwas/all"

declare -a PLAN_MODULES=() PLAN_STEPS=() PLAN_KINDS=() PLAN_FILES=()
declare -A REQUESTED=() FOUND=()
if [[ -n $STEPS ]]; then
    [[ $STEPS != ,* && $STEPS != *, && $STEPS != *,,* ]] || die "Empty step in --steps"
    IFS=, read -r -a requested_steps <<< "$STEPS"
    for step in "${requested_steps[@]}"; do REQUESTED[$step]=1; done
fi
if ((${#SCRIPT_ARGS[@]})); then
    [[ -n $STEPS && ${#REQUESTED[@]} == 1 ]] || die "Script arguments require exactly one --steps entry"
fi

# The summary uses s5 outputs: place it between signature and the final NHANES module.
while IFS=$'\t' read -r module step kind file mode; do
    [[ -n $module && $module != \#* ]] || continue
    [[ $MODULE == all || $module == "$MODULE" ||
        ( $MODULE == signature && $mode == final && -z $STEPS ) ||
        ( $MODULE == nhanes && $mode == after_nhanes && -z $STEPS ) ]] || continue
    if [[ -n $STEPS ]]; then
        [[ -n ${REQUESTED[$step]:-} ]] || continue
        FOUND[$step]=1
    elif ((LIST_STEPS == 0)); then
        [[ $mode != optional ]] || continue
        [[ $mode != final || $MODULE == all || $MODULE == signature ]] || continue
        [[ $mode != after_nhanes || $MODULE == all || $MODULE == nhanes ]] || continue
        [[ $mode != basic || $SKIP_BASIC == 0 ]] || continue
    fi
    [[ -f $SCRIPT_DIR/f/$file ]] || die "Missing script: $SCRIPT_DIR/f/$file"
    PLAN_MODULES+=("$module"); PLAN_STEPS+=("$step")
    PLAN_KINDS+=("$kind"); PLAN_FILES+=("$file")
done <<'DOSE_STEPS'
# module	step	interpreter	file	mode
foodwas	preprocess	R	foodwas.R	basic
foodwas	cohort	R	foodwas.R	basic
foodwas	patterns	R	foodwas.R	default
foodwas	items	R	foodwas.R	default
foodwas	groups	R	foodwas.R	default
foodwas	trend	R	foodwas.R	default
foodwas	rcs	R	foodwas.R	default
foodwas	plot	R	foodwas.R	default
foodwas	plot_myth	R	foodwas.R	default
ml	compare	python	ml.py	default
ml	importance	python	ml.py	default
ml	select	python	ml.py	default
ml	cox	R	ml.R	default
ml	validate	python	ml.py	default
ml	bootstrap	python	ml.py	default
ml	plot	R	ml.R	default
score	build	R	score.R	default
score	distribution	R	score.R	default
score	stroke	R	score.R	default
score	sensitivity	R	score.R	default
score	plot	R	score.R	default
score	outcomes	R	score.R	default
mechanism	protein	R	mechanism.R	default
mechanism	pathways	R	mechanism.R	default
mechanism	network	R	mechanism.R	default
mechanism	pathway_plot	R	mechanism.R	default
mechanism	metabolite	R	mechanism.R	default
mechanism	metabolite_plot	R	mechanism.R	default
mechanism	inflammation	R	mechanism.R	default
mechanism	inflammation_plot	R	mechanism.R	default
mechanism	mri_prepare	python	mechanism.py	default
mechanism	mri	R	mechanism.R	default
mechanism	mri_refit	R	mechanism.R	default
mechanism	mri_plot	python	mechanism.py	default
mechanism	mri_check	python	mechanism.py	default
mechanism	mri_wmh	R	mechanism.R	optional
mechanism	sem	R	mechanism.R	default
mechanism	sem_plot	R	mechanism.R	default
mechanism	sem_sensitivity	R	mechanism.R	default
signature	prepare	R	signature.R	default
signature	protein_lasso	python	signature.py	default
signature	protein_cox	R	signature.R	default
signature	metabolite_lasso	python	signature.py	default
signature	metabolite_cox	R	signature.R	default
signature	metabolite_plot	R	signature.R	default
mechanism	summary	R	mechanism.R	final
nhanes	score	R	nhanes.R	default
nhanes	stroke	R	nhanes.R	default
nhanes	plot	R	nhanes.R	default
nhanes	flow	R	nhanes.R	default
score	stroke_comparison	R	score.R	after_nhanes
DOSE_STEPS
for step in "${!REQUESTED[@]}"; do
    [[ -n ${FOUND[$step]:-} ]] || die "Unknown step '$step' for $MODULE; use --list-steps"
done
((${#PLAN_STEPS[@]})) || die "No steps selected"
if ((LIST_STEPS)); then
    for i in "${!PLAN_STEPS[@]}"; do
        printf '%-11s %-34s f/%s\n' "${PLAN_MODULES[$i]}" "${PLAN_STEPS[$i]}" "${PLAN_FILES[$i]}"
    done
    exit 0
fi

# Resolve relative overrides before changing the working directory.
OUTDIR=$(realpath -m -- "$OUTDIR")
UKBDIR=$(realpath -m -- "$UKBDIR")
FILESDIR=$(realpath -m -- "$FILESDIR")
NHANESDIR=$(realpath -m -- "${NHANESDIR:-$FILESDIR}")
HELPERDIR=$(realpath -m -- "$HELPERDIR")
export DOSE_OUTDIR="$OUTDIR" DOSE_PROJECT="$OUTDIR" DOSE_SCRIPT_DIR="$SCRIPT_DIR/f"
export DOSE_UKB="$UKBDIR" NHANES_DATA="$NHANESDIR" DOSE_FILES_DIR="$FILESDIR"
export DOSE_HELPER_DIR="$HELPERDIR"
export DOSE_RUNTIME_DIR="$RUNTIME_DIR"
if [[ -d $RUNTIME_DIR/R-library ]]; then
    export R_LIBS="$RUNTIME_DIR/R-library${R_LIBS:+:$R_LIBS}"
fi
export DOSE_THREADS="$JOBS" OMP_NUM_THREADS="$JOBS" OPENBLAS_NUM_THREADS="$JOBS" MKL_NUM_THREADS="$JOBS"
export DOSE_REPLACE="$REPLACE"
export PYTHONDONTWRITEBYTECODE=1

if ((DRY_RUN == 0)); then
    # Fail before creating outputs if process-tree memory control is unavailable.
    if ! command -v systemd-run >/dev/null || ! command -v systemctl >/dev/null; then
        die "systemd-run and systemctl are required; refusing to run without memory limits"
    fi
    [[ -r /sys/fs/cgroup/cgroup.controllers ]] ||
        die "cgroup v2 is required for memory limits"
    export XDG_RUNTIME_DIR=${XDG_RUNTIME_DIR:-/run/user/$(id -u)}
    export DBUS_SESSION_BUS_ADDRESS=${DBUS_SESSION_BUS_ADDRESS:-unix:path=$XDG_RUNTIME_DIR/bus}
    systemctl --user show-environment >/dev/null 2>&1 ||
        die "Cannot reach the systemd user manager; refusing to run without memory limits"
    CACHE_PYTHON=$(command -v python3) || die "python3 is required for checkpoint validation"
    # Resolve only the interpreters actually used by this selection, before any analysis.
    for kind in "${PLAN_KINDS[@]}"; do
        case "$kind" in
            R) RSCRIPT_BIN=$(command -v -- "$RSCRIPT_BIN") || die "Rscript unavailable; use --rscript" ;;
            python) PYTHON_BIN=$(command -v -- "$PYTHON_BIN") || die "Python unavailable; use --python" ;;
            *) die "Unsupported interpreter: $kind" ;;
        esac
    done
    mkdir -p "$OUTDIR/logs"
fi
export DOSE_PYTHON="$PYTHON_BIN" RSCRIPT_BIN
printf 'DOSE module: %s\nResults: %s\nUKB: %s\nNHANES: %s\nDictionaries: %s\n' "$MODULE" "$OUTDIR" "$UKBDIR" "$NHANESDIR" "$FILESDIR"
printf 'Per-step memory budget: %s GiB RAM, %s GiB swap (all workers combined); separate systemd scopes\n' "$MEMORY_LIMIT_GB" "$MEMORY_SWAP_GB"
printf 'Replace existing results: %s (verified checkpoints are reused when FALSE)\n' "$REPLACE"

run_step() {
    local module=$1 step=$2 kind=$3 file=$4 status rc scope cache_action
    local -a command=()
    case "$kind" in
        R) command=("$RSCRIPT_BIN" --vanilla "$SCRIPT_DIR/f/$file") ;;
        python) command=("$PYTHON_BIN" "$SCRIPT_DIR/f/$file") ;;
    esac
    command+=("${SCRIPT_ARGS[@]}")
    printf '\n[%s / %s] DOSE_STEP=%q' "$module" "$step" "$step"
    printf ' %q' "${command[@]}"; printf '\n'
    if ((DRY_RUN)); then return 0; fi
    mkdir -p "$OUTDIR/logs/$module"
    status="$OUTDIR/logs/$module/$step.status"
    cache_action=$("$CACHE_PYTHON" "$SCRIPT_DIR/f/checkpoint.py" begin "$module" "$step" -- "${command[@]}") || return $?
    if [[ $cache_action == SKIP ]]; then
        printf 'SKIPPED: %s/%s (inputs, code and outputs verified)\n' "$module" "$step"
        printf 'SKIPPED\nchecked=%s\n' "$(date -Is)" > "$status"
        if [[ $step == mri_prepare ]]; then
            export DOSE_MRI_RAW="$OUTDIR/ukb/05mechanism/04mri/inputs/mri_selected.tsv"
        fi
        return 0
    fi
    scope="dose-${module}-${step}-$$-${RANDOM}.scope"
    printf 'RUNNING\nstarted=%s\nscope=%s\n' "$(date -Is)" "$scope" > "$status"
    # Keep this supervisor and tee outside the scope so OOM still records FAILED.
    # The checker verifies kernel limits before exec; no environment flag bypasses it.
    if (cd "$OUTDIR" && DOSE_STEP="$step" systemd-run --user --scope --quiet --collect \
        --expand-environment=no --unit="$scope" \
        -p MemoryAccounting=yes \
        -p "MemoryMax=$MEMORY_MAX" -p "MemorySwapMax=$MEMORY_SWAP_MAX" \
        -p OOMPolicy=kill -p KillMode=control-group -p TimeoutStopSec=10s \
        -- bash "$SCRIPT_DIR/dose.sh" --internal-memory-scope "$scope" \
        "$MEMORY_MAX" "$MEMORY_SWAP_MAX" "${command[@]}") \
        2>&1 | tee "$OUTDIR/logs/$module/$step.log"; then
        if ! "$CACHE_PYTHON" "$SCRIPT_DIR/f/checkpoint.py" finish "$module" "$step" -- "${command[@]}"; then
            printf 'FAILED\nreason=checkpoint_validation\nfinished=%s\n' "$(date -Is)" > "$status"
            return 2
        fi
        printf 'COMPLETE\nfinished=%s\nscope=%s\n' "$(date -Is)" "$scope" > "$status"
    else
        rc=$?
        printf 'FAILED\nexit_code=%s\nfinished=%s\nscope=%s\n' "$rc" "$(date -Is)" "$scope" > "$status"
        printf 'FAILED: %s (see %s/logs/%s/%s.log)\n' "$step" "$OUTDIR" "$module" "$step" >&2
        printf 'Memory budget: %s GiB RAM + %s GiB swap. Scope: %s\n' "$MEMORY_LIMIT_GB" "$MEMORY_SWAP_GB" "$scope" >&2
        printf 'Scope diagnostics: journalctl --user -u %s --no-pager\n' "$scope" >&2
        return "$rc"
    fi
    # Use extracted MRI data only after mri_prepare succeeds and verifies provenance.
    if [[ $step == mri_prepare ]]; then
        export DOSE_MRI_RAW="$OUTDIR/ukb/05mechanism/04mri/inputs/mri_selected.tsv"
    fi
}

for i in "${!PLAN_STEPS[@]}"; do
    run_step "${PLAN_MODULES[$i]}" "${PLAN_STEPS[$i]}" "${PLAN_KINDS[$i]}" "${PLAN_FILES[$i]}"
done
if ((DRY_RUN)); then
    printf '\nDry run complete; no analysis was started.\n'
else
    printf '\nDOSE complete: %s\n' "$OUTDIR"
fi
