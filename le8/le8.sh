#!/usr/bin/env bash

# Public directory options use --dir-XXX; legacy spellings remain compatible.
if [[ "${BASH_SOURCE[0]}" == "$0" ]]; then
    source "$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/../0f" && pwd)/cli_dirs.sh"
    script_directory_args "$(realpath -- "${BASH_SOURCE[0]}")" "$@"
    set -- "${SCRIPT_DIRECTORY_ARGS[@]}"
    unset SCRIPT_DIRECTORY_ARGS
fi
# 🚩 le8
# LE8 FINAL interface. Default sequence: integrated multiroute C1, C2, C3, C4, C5.
set -euo pipefail
export PYTHONDONTWRITEBYTECODE=1
export TMPDIR=/tmp TMP=/tmp TEMP=/tmp PYTHONPYCACHEPREFIX=/tmp/python-cache
LE8_ROOT=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
# One logger covers the public dispatcher and all of its children.
export SCRIPT_FOREGROUND=1
export SCRIPT_LOG_DIR=${SCRIPT_LOG_DIR:-/tmp/le8-logs/runs}
if [[ $# == 0 && ${SCRIPT_CONSOLE_ACTIVE:-0} != 1 && ${SCRIPT_VERBOSE:-0} != 1 ]]; then
	exec python3 "$LE8_ROOT/../0f/console_run.py" --script "$LE8_ROOT/le8.sh" --
fi
source "$LE8_ROOT/../0f/console.sh"

le8_usage() {
	cat <<'HELP'
Modules:
    c1_correlate   Associations and all prediction methods (glmnet, LightGBM, ABM)
    c2_cause       MR and genetic decomposition
    c3_coloc       Colocalization
    c4_connect     Connections, panel validation and reconstruction
    c5_cellulation Cell annotation and CIGMA
    pgs_focus      Matched-PGS analysis
    final / share / shiny  Reports / portable viewer / local viewer
    Default: c1_correlate -> C2 -> C3 -> C4 -> C5

Main options:
    --Y CSV, --trait CSV  Outcomes [cvd_cad,ra]
    --biom CSV            Omic layers [prot,met]
    --dir-analysis-root DIR   Results [/mnt/d/analysis/le8]
    --dir-ukb-phe DIR         Phenotype/omics inputs [/mnt/d/data/ukb]
    --cores N             CPU budget [8]
    --memory-limit-gb N   Process-tree RAM cap [32]
    --replace TRUE|FALSE  Refit selected stages [FALSE]
    --methods CSV         glmnet,lightgbm,abm_case_preserving,abm_profile,abm_prevalent
                          abm_legacy,molecular_deviation,abm_deep (all eight by default)
    --device NAME         [auto], cpu or cuda; independent of scientific method
    --config FILE         C1 parameter file [config/c1_methods.json]
    --history-config FILE Historical-case assumptions [config/c1_history.json]
    --association MODE    [run], auto or skip the association analysis
    --skip-abm            Run only the conventional C1 prediction methods
    --prepared-source DIR Reuse a prepared cohort for C1 fitting
    --preflight / --dry-run  Check inputs / print commands
    --no-reindex          Read the existing Shiny index
    --port N / --host IP  Viewer address [3839 / 127.0.0.1]
    -h, --help            Help; full options and notes: 0.README.md

Example:
    cd /mnt/d/scripts/le8
    ./le8.sh --Y cvd_cad --biom prot,met --preflight
    ./le8.sh --Y cvd_cad --biom prot,met
    ./le8.sh c1_correlate --Y ra,sle2,ad,epilepsy,mi,cvd_afib,cvd_abaneu,cvd_htn,copd,asthma,masld,amr --biom prot,met --replace TRUE
    ./le8.sh final --Y cvd_cad --biom prot,met
    ./le8.sh shiny --Y cvd_cad --biom prot,met --no-reindex

    The C1 line explicitly refits; remove --replace TRUE when resuming existing results.
    The next line checks/skips completed C1 and continues C2-C5.
    All commands run in the foreground; Ctrl+C stops Shiny before the next project.
    # Optional reconstruction bootstrap and exploratory deep method:
    C4_EXPLAIN_BOOT=200 ./le8.sh c4_connect --Y cvd_cad,ra --biom prot,met
    ./le8.sh c1_correlate --methods abm_deep --device auto --Y cvd_cad --biom prot,met

    # Utility syntax; choose one subcommand (the | below means alternatives):
    # ./le8.sh abm prepare|fit|project|evaluate|demo [options]
    # ./le8.sh deviation prepare|fit|project|evaluate|merge-pgs|demo [options]
    # ./le8.sh pgs-merge --manifest FILE --output FILE
    # ./le8.sh evidence --analysis-root DIR [options]
    # ./le8.sh viewer --db FILE [--port 3839]

    # Utility help (choose a subcommand and supply its required input files):
    ./le8.sh abm --help
    ./le8.sh deviation --help
    ./le8.sh pgs-merge --help
    ./le8.sh evidence --help
    ./le8.sh viewer --help
    bash f/c2.cause.sh index-cis --help
    python f/c5.cellulation.py cigma --help
    python f/final.py report --help
    python f/final.py tables --help

    # Sensitivity methods: abm_legacy, abm_lenient, abm_deep, molecular_deviation.
    # c1_abm/c1_abm_GPU remain deprecated aliases and execute unified C1 only once.
    ./le8.sh --Y cvd_cad --cores 4
    ./le8.sh  # Default workflow
HELP
}

# Standalone model/index utilities share the public launcher and interpreter.
le8_utility() {
	local mode=$1
	shift
	local python_bin=${PYTHON:-${LE8_REPORT_PYTHON:-${PYTHON_BIN:-${HOME}/anaconda3/envs/le8/bin/python3}}}
	command -v "$python_bin" >/dev/null 2>&1 || python_bin=python3
	export PYTHONPATH="$LE8_ROOT/f${PYTHONPATH:+:$PYTHONPATH}"
	case "$mode" in
		abm)
			if [[ ${1:-} == prepare ]]; then
				shift
				"$python_bin" "$LE8_ROOT/f/abm.py" prepare "$@"
			else "$python_bin" "$LE8_ROOT/f/abm.py" "$@"; fi ;;
		deviation) "$python_bin" "$LE8_ROOT/f/deviation.py" "$@" ;;
		pgs-merge) "$python_bin" "$LE8_ROOT/f/deviation.py" merge-pgs "$@" ;;
		evidence) "$python_bin" "$LE8_ROOT/f/final.py" evidence "$@" ;;
		viewer) "${R_BIN:-Rscript}" "$LE8_ROOT/shiny/app.R" "$@" ;;
	esac
}

le8_coloc_cuda_preflight() {
	local r_bin=$1 python_bin=$2 r_libraries
	# R on Debian prepends its system library directories to LD_LIBRARY_PATH.
	# On WSL these can contain a Linux libcuda that cannot access the Windows GPU.
	# Preserve the selected R runtime's paths, but put the WSL driver first for R
	# and its Python/Conda children. Keep this runtime fix outside f/: changing
	# analysis-source fingerprints would invalidate completed scientific caches.
	if [[ -r /usr/lib/wsl/lib/libcuda.so.1 ]]; then
		r_libraries=${R_LD_LIBRARY_PATH:-}
		if [[ -z $r_libraries ]]; then
			r_libraries=$("$r_bin" --vanilla -e 'cat(Sys.getenv("LD_LIBRARY_PATH"))') || return 2
		fi
		export R_LD_LIBRARY_PATH="/usr/lib/wsl/lib${r_libraries:+:$r_libraries}"
		export LD_LIBRARY_PATH="/usr/lib/wsl/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"
	fi
	echo '[LE8] START GPU-coloc CUDA preflight (R -> Python)'
	if ! "$r_bin" --vanilla -e '
		a <- commandArgs(trailingOnly = TRUE)
		status <- system2(a[1], shQuote(c(a[2], "check-cuda")))
		if (status != 0L) cat("GPU-coloc R child LD_LIBRARY_PATH=", Sys.getenv("LD_LIBRARY_PATH"), "\n", sep = "", file = stderr())
		quit(save = "no", status = status)
	' "$python_bin" "$LE8_ROOT/f/c3.coloc_GPU.py"; then
		echo 'ERROR: GPU-coloc CUDA preflight failed before analysis.' >&2
		return 2
	fi
	echo '[LE8] DONE GPU-coloc CUDA preflight (R -> Python)'
}

le8_main() {
	case "${1:-}" in
		abm|deviation|pgs-merge|evidence|viewer) le8_utility "$@"; return $? ;;
	esac
	# Shared defaults, deliberately visible at the public entry point.
	local modules=c1_correlate,c2_cause,c3_coloc,c4_connect,c5_cellulation trait_csv=${Y:-cvd_cad,ra} biom_csv=${BIOM:-prot,met}
	local analysis_root=${LE8_ANALYSIS_ROOT:-/mnt/d/analysis/le8}
	local seed=${SEED:-2026} replace=FALSE r_bin=${R_BIN:-Rscript}
	local port=${LE8_SHINY_PORT:-3839} host=${LE8_SHINY_HOST:-127.0.0.1}
	# C2 retains the R analysis heap while external MR workers run. A 24 GiB
	# phase can leave less than one 4 GiB worker after the 4 GiB reserve.
	local cores=8 ukb_phe="" memory_limit=${LE8_MEMORY_LIMIT_GB:-32} memory_swap="" dry_run=FALSE
	local -a extra=() args=()
	# Analysis input paths and principal scientific/resource settings.
	# These are inherited by 0.engine.sh; explicit native CLI options override them.
	export UKB_PHE=${UKB_PHE:-/mnt/d/data/ukb}
	export LE8_UKB_REL_FILE=${LE8_UKB_REL_FILE:-/mnt/f/gen/ukb/37/typ/ukb_rel.dat}
	export LE8_GWAS_DIR=${LE8_GWAS_DIR:-/mnt/f/gwas/main}
	export LE8_PQTL_IV_DIR=${LE8_PQTL_IV_DIR:-/mnt/f/gwas/prot}
	export LE8_MQTL_IV_DIR=${LE8_MQTL_IV_DIR:-/mnt/f/gwas/met}
	export LE8_REFGEN_ROOT=${LE8_REFGEN_ROOT:-/mnt/f/gen/1kg}
	export LE8_GRCH=${LE8_GRCH:-auto}
	# Pool capacity coordinates concurrent runs; each phase requests only its own budget.
	export LE8_TOTAL_CORES=${LE8_TOTAL_CORES:-16}
	export LE8_RESOURCE_MEMORY_GIB=${LE8_RESOURCE_MEMORY_GIB:-32}
	export RUN_MRlink2=${RUN_MRlink2:-Top} RUN_Dandelion=${RUN_Dandelion:-Top}
	export MRLINK2_EIGH_BACKEND=${MRLINK2_EIGH_BACKEND:-auto}
	export RUN_GPU_COLOC=${RUN_GPU_COLOC:-TRUE} FINAL_NESTED_CV=${FINAL_NESTED_CV:-FALSE}
	export C2_FEATURE_SCOPE=${C2_FEATURE_SCOPE:-all_qtl} C2_MAX_FEATURES=${C2_MAX_FEATURES:-0}
	export C4_FOCUS_BUDGETS=${C4_FOCUS_BUDGETS:-5,10,50}
	export C4_MODULE_BOOT=${C4_MODULE_BOOT:-100} C4_MODULE_STABILITY=${C4_MODULE_STABILITY:-0.70}
	# Author-published CIGMA inference remains labelled as external evidence.
	if [[ -z ${C5_CIGMA_RESULTS:-} && -s /mnt/f/annot/cigma/chen2026/cigma.results.csv ]]; then
		export C5_CIGMA_RESULTS=/mnt/f/annot/cigma/chen2026/cigma.results.csv
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
		--pwas-workers | --pgs-workers | --mrlink2-workers | --mrlink2-inner-threads | --cpu-task-memory-gib | --gpu-jobs)
			extra+=("$1" "${2:?option requires a value}")
			shift 2
			;;
		--eigh-backend)
			case "${2:-}" in auto|cpu|cuda) export MRLINK2_EIGH_BACKEND=$2 ;; *) echo 'ERROR: --eigh-backend requires auto|cpu|cuda' >&2; return 2 ;; esac
			shift 2
			;;
		--run-gpu-coloc)
			local gpu_coloc_arg=${2:?--run-gpu-coloc requires TRUE or FALSE}
			case "${gpu_coloc_arg,,}" in
				true | 1 | yes | y) export RUN_GPU_COLOC=TRUE ;;
				false | 0 | no | n) export RUN_GPU_COLOC=FALSE ;;
				*) echo 'ERROR: --run-gpu-coloc requires TRUE or FALSE' >&2; return 2 ;;
			esac
			extra+=("$1" "$RUN_GPU_COLOC")
			shift 2
			;;
		--dry-run)
			dry_run=TRUE
			extra+=("$1")
			shift
			;;
		--Y | --trait | -Y)
			trait_csv=${2:?--Y requires CSV}
			shift 2
			;;
		--biom | -b)
			biom_csv=${2:?--biom requires CSV}
			shift 2
			;;
		--analysis-root | --dir-analysis-root)
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
		--ukb-phe | --dir-ukb-phe)
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
	if [[ ",$modules," == *,c3_coloc,* && $RUN_GPU_COLOC == TRUE && $dry_run != TRUE ]]; then
		le8_coloc_cuda_preflight "$r_bin" "$python_bin" || return 2
	fi
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
