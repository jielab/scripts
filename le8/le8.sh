#!/usr/bin/env bash
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
LE8 FINAL — Usage: ./le8.sh [module[,module...]] [options]

No module list: c1_correlate -> c1_abm (Twin/deviation) -> C2 -> C3 -> C4 (connect/panel validation/reconstruction) -> C5.
Run final,shiny separately to combine completed outcomes and omic layers.
Each analysis module writes cumulative performance tables/figures in its performance/ directory.
C1 refreshes comparisons with existing ABM fits; unavailable prediction models are marked not evaluated.
Completed results are checked and skipped, preserving original methods/scope.
C1 validation receipts in /tmp avoid table expansion while result files remain unchanged.
Input timestamps or code changes do not refit completed modules; use --replace TRUE.
ABM molecular-QC repair requires an explicit new fit; existing fits are not upgraded in place.
ABM uses Rdata/prot.rds and Rdata/met.rds, with no automatic plate/well/center adjustment.
Named RDS values keep their supplied scale; input/adjustment changes require a new ABM fit.
For an integrated refit use c1_abm --replace TRUE; successful replacements retain _history snapshots.
Keep the original --outer-roster. Previously inspected tests remain exploratory.
Family data are optional, including standalone prepare; absent: individual ID groups with kinship unverified.
Objective coronary endpoints are optional for reference calibration; missing endpoints disable diagnostic evaluation.
C1 and ABM publish each outcome/omic layer separately so completed layers survive later failures.
Result workbooks and reusable raw/model objects restore working inputs in /tmp.
Saved fits regenerate PNGs with same-name result workbooks; no refitting.
Native tables use XLSX and reusable objects use RDS. ABM uses aggregate XLSX/TSV and private LOCAL_ONLY arrays/models.

Examples (copy directly; run/resume C1-C5 for both omic layers):
  cd /mnt/d/scripts/le8
  ./le8.sh --Y cvd_cad --biom prot,met --preflight && \
  ./le8.sh --Y cvd_cad --biom prot,met

Run/resume C1 associations and the new ABM, then build the aggregate report:
  ./le8.sh c1_correlate,c1_abm --Y cvd_cad --biom prot,met
  ./le8.sh final --Y cvd_cad --biom prot,met

After that analysis succeeds, build the report and open Shiny:
  ./le8.sh final --Y cvd_cad --biom prot,met
  ./le8.sh shiny --Y cvd_cad --biom prot,met --no-reindex

Resume an interrupted run with valid checkpoints (no forced refit):
  ./le8.sh --Y cvd_cad --biom prot,met
Continue from C3; completed C1/C2 dependencies are checked and skipped:
  ./le8.sh c3_coloc,c4_connect,c5_cellulation --Y cvd_cad --biom prot

Modules:
  c1_correlate   C1 associations + basic clinical/biomarker prediction benchmark
  c1_abm         Twin confidence + genetic deviation + comparison with C1 benchmark
  c2_cause       MR and genetic decomposition
  c3_coloc       Colocalization and locus evidence
  c4_connect     LE8 connections, interactions, nonlinearity, panel validation
                 and paired reconstruction bootstrap (one entry point)
  c5_cellulation Cell annotation and explicit CIGMA jobs
  pgs_focus      Dedicated matched-PGS analysis
  final         Native reports + twin/deviation summaries + evidence browser; no fitting
  share         Export aggregate database, catalogue and a standalone Shiny viewer
  shiny         Rebuild the aggregate index, then serve the local Shiny app

Main options (defaults are set below in this script):
  --Y CSV, --trait CSV    cvd_cad,ra
  --biom CSV             prot,met
  --analysis-root DIR    /mnt/d/analysis/le8
  --ukb-phe DIR          External phenotype/omics input root (analysis only)
  --cores N              Phase CPU budget; default 8
  --pwas-workers N|auto Dynamic blocks of 4 features; auto up to 8, RAM may reduce to serial
  --pgs-workers N|auto  Six-model checkpoints; auto up to 4, RAM may reduce to serial
  --mrlink2-workers N   Concurrent MR tasks; default 4
  --mrlink2-inner-threads N  Threads per MR task; default 1
  --cpu-task-memory-gib N  Per-worker admission estimate (not a hard RAM cap)
  --gpu-jobs N          Exclusive GPU phases; this release requires 1
  --eigh-backend MODE   MR-link-2 eigendecomposition: auto (default, GPU preferred), cpu, cuda
  --seed N               2026
  --replace TRUE|FALSE   FALSE; check/skip completed outputs; TRUE refits selected stages
  --r-bin FILE           Rscript
  --port N               3839
  --host ADDRESS         127.0.0.1
  --memory-limit-gb N    Process-tree RAM cap; default 32 GiB; 0 disables it
  --memory-swap-gb N     Explicit swap cap
  --run-abm              Include ABM with an explicit module list
  --skip-abm             Skip ABM in the default C1-C5 sequence
  --fit-reference, --fit-joint, --fit-genetic  Explicit Final fitting
  --details       Regenerate detailed R plots from saved results
  --index-only           Build index without report figures or R
  --prepare-only         Prepare without launching the Shiny server
  --no-reindex           Read the existing Shiny index
  --preflight            Check dependencies / inputs (CUDA checked before native scans)
  --dry-run              Print planned commands
  --shiny-review         Validate prepared sources and exit
  --out DIR              Must be <analysis-root>/final
  --group-file FILE      Shared eid,group mapping (LE8_GROUP_FILE)
  --group-col COLUMN     Shared phenotype family column (LE8_GROUP_COLUMN)
  --end-date YYYY-MM-DD  Administrative cutoff (DATE_FOLLOW_END)
  --Y-date COLUMN       Outcome diagnosis date, shared with ABM
  --outer-roster FILE    Shared eid,role table (training/test), for paired comparisons
  --atlas, --universe, --panels, --contrasts FILE
  --matched-draws N      C5 matched draws: default 500 (C5_MATCHED_DRAWS); 0 or >=100
  --cigma-manifest, --cigma-results, --cigma-cells FILE
                         Downloaded Chen 2026 summaries are used when installed;
                         C5_CIGMA_RESULTS overrides the published-summary path.
  --allow-untested-cigma C5 explicit version override
  --no-plots            C5 only: skip annotation plots
  --max-table-mb N        128
  --abm-config FILE       Twin model configuration (default config/abm.json)
  --abm-host-args STRING  Explicit phenotype/omics input overrides; parsed without shell evaluation
  --deviation-config FILE  Deviation configuration; defaults per omic layer
                           Explicit known_cad_state / coronary_calcium tasks use separate result directories
  --objective-file FILE    Optional eid-keyed objective coronary endpoint table
  --pgs-root PATH     PGS file/directory; use {biom} for layer-specific paths
  --pgs-mapping FILE  Explicit PGS assay mapping
  --remerge-pgs       Merge retained native scores at full precision before C1
  --share-out ZIP        Portable viewer package (default: ../le8-share.zip)
  --strict               Fail on blocking report audit findings
  --help                 Show this help

Native data options: --gwas-dir, --pqtl-dir, --mqtl-dir, --refgen-root, --grch.
Native analysis options also pass through: --Y-date, --vars.adj,
--white-only, --run-mrlink2, --run-dandelion, --run-gpu-coloc, --nested-cv.
Outputs: <analysis-root>/<Y>/<biom>/<module>/; optional Final fits: final/<Y>/;
summary report: final/; aggregate database: evidence/; execution logs and caches: /tmp/.
Question-led C1 → C2 → C4 overview and figures: final/; Shiny opens the aggregate evidence browser.
LE8 reconstruction uncertainty: C4_EXPLAIN_BOOT=200 ./le8.sh c4_connect --Y cvd_cad,ra --biom prot,met
Final report files overwrite the fixed final/ destination.
Per-run defaults: 8 CPU cores and 32 GiB RAM.
The RAM cap includes the waiting R process, MR workers and scheduler headroom.
Shared pool across LE8 runs: LE8_TOTAL_CORES=16 and LE8_RESOURCE_MEMORY_GIB=32 (GiB).
All concurrent runs must use the same pool settings; per-run requests may be smaller.
Screen: major stages and one failure summary; full logs: /tmp/le8-logs/runs/.
Runs stay in the terminal foreground: Ctrl+C stops the job; Ctrl+Z suspends the whole job.
SCRIPT_VERBOSE=1 streams all diagnostics to the terminal.
PYTHON_BIN / LE8_REPORT_PYTHON and ABM_PYTHON select existing environments.
ABM_PYTHON defaults to the report interpreter; select a separate backend environment if needed.
GPU-coloc checks CUDA through R before analysis; WSL driver libraries take priority.
Twin confidence and molecular deviation run on CPU.

Additional module entry points:
  ./le8.sh abm prepare|fit|project|evaluate|demo [options]
  ./le8.sh deviation prepare|fit|project|evaluate|merge-pgs|demo [options]
  ./le8.sh pgs-merge --manifest FILE --output FILE
  ./le8.sh evidence --analysis-root DIR [options]
  ./le8.sh viewer --db FILE [--port 3839]
  bash f/c2.cause.sh index-cis --help
  python f/c5.cellulation.py cigma --help
  python f/final.py report --help
  python f/final.py tables --help
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
	local modules=c1_correlate,c1_abm,c2_cause,c3_coloc,c4_connect,c5_cellulation trait_csv=${Y:-cvd_cad,ra} biom_csv=${BIOM:-prot,met}
	local analysis_root=${LE8_ANALYSIS_ROOT:-/mnt/d/analysis/le8}
	local seed=${SEED:-2026} replace=FALSE r_bin=${R_BIN:-Rscript}
	local port=${LE8_SHINY_PORT:-3839} host=${LE8_SHINY_HOST:-127.0.0.1}
	# C2 retains the R analysis heap while external MR workers run. A 24 GiB
	# phase can leave less than one 4 GiB worker after the 4 GiB reserve.
	local cores=8 ukb_phe="" memory_limit=${LE8_MEMORY_LIMIT_GB:-32} memory_swap="" dry_run=FALSE
	local -a extra=() args=()
	# Analysis input paths and principal scientific/resource settings.
	# These are inherited by 0.engine.sh; explicit native CLI options override them.
	export UKB_PHE=${UKB_PHE:-/mnt/d/data/ukb/phe}
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
	if [[ -z ${C5_CIGMA_RESULTS:-} && -s /mnt/d/data/reference/cigma/chen2026/cigma.results.csv ]]; then
		export C5_CIGMA_RESULTS=/mnt/d/data/reference/cigma/chen2026/cigma.results.csv
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
