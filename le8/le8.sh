#!/usr/bin/env bash

# Display standalone help before loading runtimes or scheduling work.
if [[ ${BASH_SOURCE[0]} == "$0" && ( ${1:-} == -h || ${1:-} == --help ) ]]; then
    cat <<'PUBLIC_ENTRY_HELP'
LE8 analysis

Modules:
    c1_correlate   Associations and selected prediction methods
    c1_reliability Evaluate frozen glmnet reliability
    c2_cause       MR and genetic decomposition (CPU)
    c3_coloc       Colocalization (CPU)
    c4_connect    Connections and panel validation (CPU)
    c5_cellulation Cell annotation (CPU)
    pgs_focus     Matched-PGS analysis
    final         Aggregate reports
    share         Portable viewer export
    shiny         Local viewer
    install       Install/check the runtime separately
    Default: c1_correlate,c2_cause,c3_coloc,c4_connect,c5_cellulation

Main options:
    --submit hpc|local|auto    Machine selection [auto]
    --Y CSV                   Outcomes [cvd_cad,ra]
    --biom CSV                Molecular layers [prot,met]
    --methods CSV             glmnet,lightgbm,abm,abm_deep [default]
    --abm-deep-run cpu|gpu|auto Device for selected abm_deep [cpu]; auto means cpu
    --replace TRUE|FALSE      Overwrite requested results [TRUE]; FALSE reuses complete outputs
    --analysis-root DIR       Analysis output directory
    --cores N                 CPU budget [8]
    --memory-limit-gb N        Host RAM budget [32 GiB]
    --queue NAME              CPU queue
    --gpu-queue NAME           Queue used only for explicitly selected deep GPU work
    --gpu-cores N              CPU slots in a deep GPU job [8]
    --gpu-memory-limit-gb N    GPU-node host RAM cap; not GPU VRAM
    --hours N                 Allocation time limit
    --jobs-N N                Concurrent local jobs / active submissions
    --wait TRUE|FALSE         Wait for the complete submitted DAG [FALSE]
    --plan-only               Write reviewable commands without executing/submitting
    --config FILE             Method parameters [config/c1_methods.json]
    --association run|skip    C1 association stage [run]
    --pgs TRUE|FALSE          Molecular PGS associations [FALSE]

ABM uses frozen molecular neighbours and always runs on CPU. Deep is a separate
Transformer model with neighbour support diagnostics. GPU training is followed
by CPU evaluation on the same frozen cohort. Scheduling is automatic.

Examples:
    ./le8.sh --Y ad,amr,asthma,copd,cvd_cad,cvd_htn,cvd_pad,depression,epilepsy,ra,sle,stroke --biom prot,met --submit hpc --queue 33 --cores 40
PUBLIC_ENTRY_HELP
    exit 0
fi

# le8
# LE8 FINAL interface. Default sequence: integrated multiroute C1, C2, C3, C4, C5.
set -euo pipefail
export PYTHONDONTWRITEBYTECODE=1
LE8_ROOT=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
le8_gpu_queue_override=${LE8_GPU_QUEUE:-}
if [[ -n ${_OS_DEFAULT_LE8_GPU_QUEUE:-} && $le8_gpu_queue_override == "$_OS_DEFAULT_LE8_GPU_QUEUE" ]]; then
    le8_gpu_queue_override=""
fi
source "$LE8_ROOT/../0f/workflow_env.sh"
workflow_init le8 "$@"
set -- "${WORKFLOW_ARGS[@]}"
source "$LE8_ROOT/f/runtime_env.sh"
unset le8_gpu_queue_override
# Public runs produce inspectable .cmd files per trait/biom/resource phase. Planning uses only Python
# stdlib; no R, cohort preparation or fitting can execute on an HPC login node.
if [[ ${WORKFLOW_WORKER:-0} != 1 ]]; then
    le8_plan=TRUE
    for le8_arg in "$@"; do
        case "$le8_arg" in -h|--help|help) le8_plan=FALSE ;; esac
    done
    case "${1:-}" in evidence|viewer) le8_plan=FALSE ;; esac
    if [[ $le8_plan == TRUE ]]; then
        exec "${LE8_PLANNER_PYTHON:-python3}" "$LE8_ROOT/f/submit.py" "$@"
    fi
    if [[ $submit == hpc && ${1:-} =~ ^(evidence|viewer)$ && -z ${LSB_JOBID:-} ]]; then
        echo 'Standalone utilities on HPC require an allocated LSF job. Use the scope workflow for scheduled analysis.' >&2
        exit 2
    fi
fi
# One logger covers the public dispatcher and all of its children.
export SCRIPT_FOREGROUND=1
export SCRIPT_LOG_DIR=${SCRIPT_LOG_DIR:-${DIR_RUNTIME_TMP}/le8-logs/runs}
if [[ $# == 0 && ${SCRIPT_CONSOLE_ACTIVE:-0} != 1 && ${SCRIPT_VERBOSE:-0} != 1 ]]; then
	exec python3 "$LE8_ROOT/../0f/console_run.py" --script "$LE8_ROOT/le8.sh" --
fi
source "$LE8_ROOT/../0f/console.sh"

le8_usage() {
	bash "${BASH_SOURCE[0]}" --help
}

# Standalone model/index utilities share the public launcher and interpreter.
le8_utility() {
	local mode=$1
	shift
	local python_bin=${PYTHON:-${LE8_REPORT_PYTHON:-${PYTHON_BIN:-${HOME}/anaconda3/envs/le8/bin/python3}}}
	command -v "$python_bin" >/dev/null 2>&1 || python_bin=python3
	export PYTHONPATH="$LE8_ROOT/f${PYTHONPATH:+:$PYTHONPATH}"
	case "$mode" in
		evidence) "$python_bin" "$LE8_ROOT/f/final.py" evidence "$@" ;;
		viewer) "${R_BIN:-Rscript}" "$LE8_ROOT/shiny/app.R" "$@" ;;
	esac
}


le8_main() {
	case "${1:-}" in
		evidence|viewer) le8_utility "$@"; return $? ;;
	esac
	# Shared defaults, deliberately visible at the public entry point.
	local modules=c1_correlate,c2_cause,c3_coloc,c4_connect,c5_cellulation trait_csv=${Y:-cvd_cad,ra} biom_csv=${BIOM:-prot,met}
	local analysis_root=${LE8_ANALYSIS_ROOT:-$ANALYSIS_ROOT/le8}
	local seed=${SEED:-2026} replace=TRUE r_bin=${R_BIN:-Rscript}
	local port=${LE8_SHINY_PORT:-3839} host=${LE8_SHINY_HOST:-127.0.0.1}
	# C2 retains the R analysis heap while external MR workers run. A 24 GiB
	# phase can leave less than one 4 GiB worker after the 4 GiB reserve.
	local cores=8 ukb_phe="" memory_limit=${LE8_MEMORY_LIMIT_GB:-32} memory_swap="" dry_run=FALSE
	local -a extra=() args=()
	# Analysis input paths and principal scientific/resource settings.
	# These are inherited by 0.engine.sh; explicit native CLI options override them.
	export UKB_PHE=${UKB_PHE:-${DIR_UKB:?Initialize paths with 0f/0os.sh}}
	export LE8_UKB_REL_FILE=${LE8_UKB_REL_FILE:-${DIR_DATA:-${DIR_BIG}}/gen/ukb/37/typ/ukb_rel.dat}
	export LE8_GWAS_DIR=${LE8_GWAS_DIR:-${DIR_DATA:-${DIR_BIG}}/gwas/main}
	export LE8_PQTL_IV_DIR=${LE8_PQTL_IV_DIR:-${DIR_DATA:-${DIR_BIG}}/gwas/prot}
	export LE8_MQTL_IV_DIR=${LE8_MQTL_IV_DIR:-${DIR_DATA:-${DIR_BIG}}/gwas/met}
	export LE8_REFGEN_ROOT=${LE8_REFGEN_ROOT:-${DIR_DATA:-${DIR_BIG}}/gen/1kg}
	export LE8_GRCH=${LE8_GRCH:-auto}
	# Pool capacity coordinates concurrent runs; each phase requests only its own budget.
	export LE8_TOTAL_CORES=${LE8_TOTAL_CORES:-16}
	export LE8_RESOURCE_MEMORY_GIB=${LE8_RESOURCE_MEMORY_GIB:-64}
	export RUN_MRlink2=${RUN_MRlink2:-Top} RUN_Dandelion=${RUN_Dandelion:-Top}
	export MRLINK2_EIGH_BACKEND=cpu
	export FINAL_NESTED_CV=${FINAL_NESTED_CV:-FALSE}
	export C2_FEATURE_SCOPE=${C2_FEATURE_SCOPE:-all_qtl} C2_MAX_FEATURES=${C2_MAX_FEATURES:-0}
	export C4_FOCUS_BUDGETS=${C4_FOCUS_BUDGETS:-5,10,50}
	export C4_MODULE_BOOT=${C4_MODULE_BOOT:-100} C4_MODULE_STABILITY=${C4_MODULE_STABILITY:-0.70}
	# Author-published CIGMA inference remains labelled as external evidence.
	if [[ -z ${C5_CIGMA_RESULTS:-} && -s $DIR_ANNOT/le8/cigma/chen2026/cigma.results.csv ]]; then
		export C5_CIGMA_RESULTS=$DIR_ANNOT/le8/cigma/chen2026/cigma.results.csv
	fi
	if [[ $# -gt 0 && $1 != -* ]]; then
		modules=$1
		shift
	fi
	while (($#)); do
		case "$1" in
		help | -h | --help)
			le8_usage
			return 0
			;;
		--pwas-workers | --pgs-workers | --mrlink2-workers | --mrlink2-inner-threads | --cpu-task-memory-gib)
			extra+=("$1" "${2:?option requires a value}")
			shift 2
			;;
		--dry-run)
			dry_run=TRUE
			extra+=("$1")
			shift
			;;
		--Y)
			trait_csv=${2:?--Y requires CSV}
			shift 2
			;;
		--biom)
			biom_csv=${2:?--biom requires CSV}
			shift 2
			;;
		--analysis-root)
			analysis_root=${2:?--analysis-root requires DIR}
			shift 2
			;;
		--seed)
			seed=${2:?--seed requires N}
			shift 2
			;;
		--replace)
			replace=${2:?--replace requires TRUE or FALSE}
			shift 2
			;;
		--r-bin)
			r_bin=${2:?--r-bin requires FILE}
			shift 2
			;;
		--port)
			port=${2:?--port requires N}
			shift 2
			;;
		--host)
			host=${2:?--host requires ADDRESS}
			shift 2
			;;
		--cores)
			cores=${2:?--cores requires N}
			shift 2
			;;
		--ukb-phe)
			ukb_phe=${2:?--ukb-phe requires DIR}
			shift 2
			;;
		--ukb-rel-file)
			extra+=("$1" "${2:?--ukb-rel-file requires FILE}")
			shift 2
			;;
		--memory-limit-gb)
			memory_limit=${2:?--memory-limit-gb requires N}
			shift 2
			;;
		--memory-swap-gb)
			memory_swap=${2:?--memory-swap-gb requires N}
			shift 2
			;;
		*)
			extra+=("$1")
			shift
			;;
		esac
	done
	[[ $modules != help ]] || {
		le8_usage
		return 0
	}
	local python_bin=${LE8_REPORT_PYTHON:-${PYTHON_BIN:-}}
	if [[ -z $python_bin ]]; then
		if [[ -x ${HOME}/anaconda3/envs/le8/bin/python3 ]]; then
			python_bin=${HOME}/anaconda3/envs/le8/bin/python3
		else python_bin=python3; fi
	fi
	command -v "$python_bin" >/dev/null || {
		echo "ERROR: Python not found: $python_bin" >&2
		return 2
	}
	export PYTHON_BIN=$python_bin LE8_ANALYSIS_ROOT=$analysis_root DIRSCRIPT=$LE8_ROOT LE8_FDIR=$LE8_ROOT/f
	args=("$modules" --Y "$trait_csv" --biom "$biom_csv" --analysis-root "$analysis_root"
		--seed "$seed" --replace "$replace" --r-bin "$r_bin" --port "$port" --host "$host")
	[[ -z $cores ]] || args+=(--cores "$cores")
	[[ -z $ukb_phe ]] || args+=(--ukb-phe "$ukb_phe")
	[[ -z $memory_limit ]] || args+=(--memory-limit-gb "$memory_limit")
	[[ -z $memory_swap ]] || args+=(--memory-swap-gb "$memory_swap")
	"$python_bin" "$LE8_ROOT/f/0.common.py" pipeline "${args[@]}" "${extra[@]}"
}
# Parse invocation and exit together before any long-running child starts.
{
	le8_main "$@"
	exit $?
}
