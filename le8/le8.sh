#!/usr/bin/env bash


# 🚩 le8
# Public LE8 interface. Default sequence: C1 (including reference ABM), C2, C3, C4, C5.
set -euo pipefail
export PYTHONDONTWRITEBYTECODE=1
export TMPDIR=/tmp TMP=/tmp TEMP=/tmp PYTHONPYCACHEPREFIX=/tmp/python-cache
LE8_ROOT=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)

le8_usage() {
	cat <<'HELP'
Usage: ./le8.sh [module[,module...]] [options]

No module list: C1 + reference ABM -> C2 -> C3 -> C4 (connect/panel validation) -> C5.
Run final,shiny separately to combine completed outcomes and omic layers.
Completed results keep their original methods/scope; --replace TRUE requests new fits.
Result workbooks and named participant RDS files restore numerical inputs in /tmp.
Saved fits regenerate PNGs with same-name result workbooks; no refitting.
Participant tables use descriptive names such as test_individuals.rds; no workbook export.

Examples:
  cd /mnt/d/scripts/le8
  ./install.sh                         # lightweight report + Shiny dependencies
  ./install.sh --abm                   # reference ABM dependencies for default analysis
  ./le8.sh final,shiny --Y cvd_cad,ra --biom prot,met
  ./le8.sh final --index-only --Y cvd_cad,ra --biom prot,met
  ./le8.sh shiny --no-reindex --port 3839
  ./le8.sh c1_correlate --Y cvd_cad --biom prot --preflight
  ./le8.sh c2_cause,c3_coloc --Y cvd_cad --biom prot --dry-run
  ./le8.sh c4_connect,c4_panel_validation --Y cvd_cad --biom prot
  ./le8.sh c4_explain --Y cvd_cad --biom prot
  ./le8.sh c5_cellulation --Y cvd_cad --biom prot
  ./le8.sh c1_abm --Y cvd_cad --biom prot
  ./le8.sh final --fit-joint --Y cvd_cad --biom prot,met --replace TRUE

Modules:
  c1_correlate   Measured/PGS associations, temporal analyses and enrichment
  c1_abm         Agent-based modeling (ABM): selective reference and TabICLv2
  c2_cause       MR and genetic decomposition
  c3_coloc       Colocalization and locus evidence
  c4_connect     LE8 connections, proxies, mediation, interactions and nonlinearity
  c4_panel_validation  Frozen matched-budget prediction and reconstruction
  c4_explain     Optional paired reconstruction bootstrap
  c5_cellulation Cell annotation and explicit CIGMA jobs
  pgs_focus      Dedicated matched-PGS analysis
  final         Existing aggregate reports; aggregate reporting, no fitting
  share         Export aggregate workbooks, PNGs and a standalone Shiny viewer
  shiny         Rebuild the aggregate index, then serve the local Shiny app

Main options (defaults are set below in this script):
  --Y CSV, --trait CSV    cvd_cad,ra
  --biom CSV             prot,met
  --analysis-root DIR    /mnt/d/analysis/le8
  --ukb-phe DIR          External phenotype/omics input root (analysis only)
  --cores N              Explicit worker count; native default 1, ABM 16
  --seed N               2026
  --replace TRUE|FALSE   FALSE; validate/reuse completed outputs; TRUE refits selected stages
  --r-bin FILE           Rscript
  --port N               3839
  --host ADDRESS         127.0.0.1
  --memory-limit-gb N    Explicit process-tree RAM cap; 0 disables it
  --memory-swap-gb N     Explicit swap cap
  --abm-backend METHOD   reference | tabicl | both (default reference)
  --abm-args STRING      Extra backend arguments, parsed without shell evaluation
  --run-abm              Include ABM with an explicit module list
  --skip-abm             Skip ABM in the default C1-C5 sequence
  --fit-reference, --fit-joint, --fit-genetic  Explicit Final fitting
  --details       Regenerate detailed R plots from saved results
  --index-only           Build index without report figures or R
  --prepare-only         Prepare without launching the Shiny server
  --no-reindex           Read the existing Shiny index
  --preflight            Check dependencies / inputs (ABM checked before native scans)
  --dry-run              Print planned commands
  --shiny-review         Validate prepared sources and exit
  --out DIR              Must be <analysis-root>/final
  --group-file FILE      Shared eid,group mapping (LE8_GROUP_FILE)
  --group-col COLUMN     Shared phenotype family column (LE8_GROUP_COLUMN)
  --end-date YYYY-MM-DD  Administrative cutoff (DATE_FOLLOW_END)
  --Y-date COLUMN       Outcome diagnosis date, shared with ABM
  --outer-roster FILE    Shared eid,role table (training/test), for paired comparisons
  --atlas, --universe, --panels, --contrasts FILE
  --matched-draws N      C5 matched draws: 0 or >=100
  --cigma-manifest, --cigma-results, --cigma-cells FILE
  --allow-untested-cigma C5 explicit version override
  --no-plots            C5 only: skip annotation plots
  --abm-root, --abm-tf-root DIR  Explicit external import sources
  --max-table-mb N        128
  --share-out ZIP        Portable viewer package (default: ../le8-share.zip)
  --strict               Fail on blocking report audit findings
  --help                 Show this help

Native data options: --gwas-dir, --pqtl-dir, --mqtl-dir, --refgen-root, --grch.
Native analysis options also pass through: --Y-date, --vars.adj,
--white-only, --run-mrlink2, --run-dandelion, --run-gpu-coloc, --nested-cv.
Outputs: <analysis-root>/<Y>/<biom>/<module>/; optional Final fits: final/<Y>/;
summary report: final/; index: shiny/; execution logs and caches: /tmp/.
Question-led overview and Fig6–8: final/; Shiny opens the research-question view.
LE8 reconstruction uncertainty: C4_EXPLAIN_BOOT=200 ./le8.sh c4_explain --Y cvd_cad,ra --biom prot,met
Final report files overwrite the fixed final/ destination.
PYTHON_BIN / LE8_REPORT_PYTHON and ABM_PYTHON select existing environments.
ABM_PYTHON defaults to the report interpreter; select a separate backend environment if needed.

Additional module entry points:
  python f/c1.abm.py annotations --help
  python f/c1.abm.py figures RUN_DIR
  bash f/c2.cause.sh index-cis --help
  python f/c5.cellulation.py cigma --help
  python f/final.py report --help
  python f/final.py tables --help
HELP
}

le8_main() {
	# Shared defaults, deliberately visible at the public entry point.
	local modules=c1_correlate,c1_abm,c2_cause,c3_coloc,c4_connect,c4_panel_validation,c5_cellulation trait_csv=${Y:-cvd_cad,ra} biom_csv=${BIOM:-prot,met}
	local analysis_root=${LE8_ANALYSIS_ROOT:-/mnt/d/analysis/le8}
	local seed=${SEED:-2026} replace=FALSE r_bin=${R_BIN:-Rscript}
	local port=${LE8_SHINY_PORT:-3839} host=${LE8_SHINY_HOST:-127.0.0.1}
	local backend=reference cores="" ukb_phe="" memory_limit="" memory_swap=""
	local -a extra=() args=()
	# Analysis input paths and principal scientific/resource settings.
	# These are inherited by 0.engine.sh; explicit native CLI options override them.
	export UKB_PHE=${UKB_PHE:-/mnt/d/data/ukb/phe}
	export LE8_GWAS_DIR=${LE8_GWAS_DIR:-/mnt/f/gwas/main}
	export LE8_PQTL_IV_DIR=${LE8_PQTL_IV_DIR:-/mnt/f/gwas/prot}
	export LE8_MQTL_IV_DIR=${LE8_MQTL_IV_DIR:-/mnt/f/gwas/met}
	export LE8_REFGEN_ROOT=${LE8_REFGEN_ROOT:-/mnt/f/gen/1kg}
	export LE8_GRCH=${LE8_GRCH:-auto}
	export RUN_MRlink2=${RUN_MRlink2:-Top} RUN_Dandelion=${RUN_Dandelion:-Top}
	export RUN_GPU_COLOC=${RUN_GPU_COLOC:-TRUE} FINAL_NESTED_CV=${FINAL_NESTED_CV:-FALSE}
	export C2_FEATURE_SCOPE=${C2_FEATURE_SCOPE:-all_qtl} C2_MAX_FEATURES=${C2_MAX_FEATURES:-0}
	export C4_FOCUS_BUDGETS=${C4_FOCUS_BUDGETS:-5,10,50}
	export C4_MODULE_BOOT=${C4_MODULE_BOOT:-100} C4_MODULE_STABILITY=${C4_MODULE_STABILITY:-0.70}
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
		--abm-backend | --abm-engine)
			backend=${2:?--abm-backend requires METHOD}
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
		--seed "$seed" --replace "$replace" --r-bin "$r_bin" --port "$port" --host "$host"
		--abm-backend "$backend")
	[[ -z $cores ]] || args+=(--cores "$cores")
	[[ -z $ukb_phe ]] || args+=(--ukb-phe "$ukb_phe")
	[[ -z $memory_limit ]] || args+=(--memory-limit-gb "$memory_limit")
	[[ -z $memory_swap ]] || args+=(--memory-swap-gb "$memory_swap")
	"$python_bin" "$LE8_ROOT/f/0.common.py" dispatch "${args[@]}" "${extra[@]}"
}
# Parse invocation and exit together before any long-running child starts.
{
	le8_main "$@"
	exit $?
}
