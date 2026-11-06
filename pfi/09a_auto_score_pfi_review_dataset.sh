#!/usr/bin/env bash
# BEGIN QIMING RESOURCE EXAMPLES (Admin, 2026-11-05)
# Queue 33: 40 cores / 192 GiB per node. Reserve 10/20/30/40 cores by
# rounding max(required_cores/40, required_memory_GiB/192) UP to 25% steps.
# Example for a 24-core, 160-GiB CPU worker:
#   bsub -q 33 -n 40 -W 24:00 -M 160GB -R "span[hosts=1] rusage[mem=160GB/host]" ...
# A small <=10-core, <=48-GiB job uses -q 33 -n 10; request its actual RAM.
# For exact CPU counts, use Admin-suggested v4-3t or v3-6t after queueinfo.
# The examples are allocation guidance, not measured peak RAM for every stage.
# This entry initializes LOCAL paths and is not a supported --submit hpc pipeline.
# The bsub example above applies only after its worker paths/runtime are ported.
# Details and queue limits: ../0f/HPC_SUBMISSION_POLICY.md
# END QIMING RESOURCE EXAMPLES

# Resolve shared roots from this checkout, independent of the current directory.
_paths_directory=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)
while [[ ! -f "$_paths_directory/0f/0os.sh" ]]; do
    [[ $_paths_directory != / ]] || { echo "Cannot locate 0f/0os.sh" >&2; return 2 2>/dev/null || exit 2; }
    _paths_directory=$(dirname -- "$_paths_directory")
done
submit=${submit:-${SUBMIT:-local}}
source "$_paths_directory/0f/0os.sh"
unset _paths_directory


# Resolve --dir-dat/--dir-big/--dir-tmp through the shared path table.
source "$SCRIPTS_ROOT/0f/local_paths.sh"
local_workflow_init "$@" || exit $?
set -- "${WORKFLOW_ARGS[@]}"

# Public directory options use --dir-XXX; legacy spellings remain compatible.
if [[ "${BASH_SOURCE[0]}" == "$0" ]]; then
    source "${SCRIPTS_ROOT}/0f/cli_dirs.sh"
    script_directory_args "$(realpath -- "${BASH_SOURCE[0]}")" "$@"
    set -- "${SCRIPT_DIRECTORY_ARGS[@]}"
    unset SCRIPT_DIRECTORY_ARGS
fi

# Direct CLI help must return before setup, downloads or analysis.
if [[ "${BASH_SOURCE[0]}" == "$0" ]]; then
	case "${1:-}" in
		-h | --help)
			cat <<'HELP'
Modules:
    Default  Create weak-supervision review scores

Main options:
    INPUT_CSV OUTPUT_CSV  Optional positional score paths
    ANALYSIS_ROOT / RAW_ROOT  Results/raw data roots
    PFI_PYTHON_BIN / PFI_RSCRIPT_BIN  Interpreters
    Defaults: conf/pfi_paths.sh; workflow order: ./pfi.sh -h
    -h, --help  Show help

Example:
    ./09a_auto_score_pfi_review_dataset.sh
HELP
			exit 0
			;;
	esac
fi


# 🚩 09a_auto_score_pfi_review_dataset
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PFI_ROOT="$SCRIPT_DIR"
cd "$PFI_ROOT"
source "$PFI_ROOT/conf/pfi_paths.sh"
mkdir -p "$DAT_DIR" "$LOG_DIR"

: "${PFI_ARTICLE_METRICS_CSV:=}"
: "${PFI_HIGH_ANCHOR_CSV:=}"
: "${PFI_LOW_ANCHOR_CSV:=}"
: "${PFI_AUTO_SCORE_MODE:=weak_supervision}"

in_csv="${1:-$PFI_PAPER_SCORE_REVIEW_CSV}"
out_csv="${2:-$PFI_PAPER_SCORE_REVIEW_PSEUDO}"

if [[ ! -s "$in_csv" ]]; then
	echo "ERROR: blinded review CSV is missing: $in_csv" >&2
	echo "Run: ./09a_make_pfi_review_inputs.sh" >&2
	exit 1
fi

extra_args=()
[[ -n "$PFI_ARTICLE_METRICS_CSV" && -s "$PFI_ARTICLE_METRICS_CSV" ]] && extra_args+=(--article_metrics "$PFI_ARTICLE_METRICS_CSV")
[[ -n "$PFI_HIGH_ANCHOR_CSV" && -s "$PFI_HIGH_ANCHOR_CSV" ]] && extra_args+=(--high_anchor_csv "$PFI_HIGH_ANCHOR_CSV")
[[ -n "$PFI_LOW_ANCHOR_CSV" && -s "$PFI_LOW_ANCHOR_CSV" ]] && extra_args+=(--low_anchor_csv "$PFI_LOW_ANCHOR_CSV")

"$PFI_PYTHON_BIN" -u f/auto_score_pfi_review_dataset.py \
	--review_csv "$in_csv" \
	--review_key "$PFI_PAPER_SCORE_REVIEW_KEY" \
	--features "$PFI_ARTICLE_FEATURES" \
	--out_csv "$out_csv" \
	--audit_csv "$PFI_PAPER_SCORE_REVIEW_PSEUDO_AUDIT" \
	--mode "$PFI_AUTO_SCORE_MODE" \
	"${extra_args[@]}" \
	2>&1 | tee "$LOG_DIR/09a_auto_score_pfi_review_dataset.log"

cat >&2 <<EOF
OK: pseudo-scored review file written to:
  $out_csv
Use it for training if manual expert labels are unavailable:
  ./pfi.sh prepare-training
EOF
