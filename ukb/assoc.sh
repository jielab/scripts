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
source "${SCRIPTS_ROOT}/0f/console.sh"
set -euo pipefail


# 🚩 Command-line help
usage() {
	cat <<'HELP'
Modules:
    forest / phewas / medwas  Association plots / PheWAS / mediation

Main options:
    --Xs CSV --Ys CSV  Exposures and outcomes
    --varX CSV  Covariates
    --type TYPE  Outcome type
    -h, --help  Show options and examples below

Example:
    cd ${SCRIPTS_ROOT}/0LOCAL/ukb

    ./assoc.sh forest --Xs abo.AB.rs8176746_G,apoe.rs7412_C --Ys ihd,stroke --type t2e
    ./assoc.sh phewas save_rds --Xs alc.aldh2.rs671_G,fto.rs9939609_T --varX age,sex,PC1,PC2
    ./assoc.sh phewas plot --Xs alc.aldh2.rs671_G,fto.rs9939609_T
    ./assoc.sh medwas --Xs alc.aldh2.rs671_G,fto.rs9939609_T --varX age,sex,PC1,PC2

    ./assoc.sh forest --Xs abo.AB.rs8176746_G --Ys ihd
    ./assoc.sh forest --Xs abo.AB.rs8176746_G --Ys ihd,stroke --type t2e
HELP
}

this_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
dir0=$DIR_DAT
phedir=""
outdir=""
script_dir=""
assoc_R=""
filtered_args=()


# 🚩 Command-line arguments
while (($#)); do
	case "$1" in
		--data-root)
			dir0=${2:?ERROR: --data-root requires DIR}
			shift 2
			;;
		--phenotype-dir)
			phedir=${2:?ERROR: --phenotype-dir requires DIR}
			shift 2
			;;
		--output-dir)
			outdir=${2:?ERROR: --output-dir requires DIR}
			shift 2
			;;
		--script-dir)
			script_dir=${2:?ERROR: --script-dir requires DIR}
			shift 2
			;;
		--r-script)
			assoc_R=${2:?ERROR: --r-script requires FILE}
			shift 2
			;;
		*)
			filtered_args+=("$1")
			shift
			;;
	esac
done
set -- "${filtered_args[@]}"
phedir=${phedir:-$DIR_UKB}
outdir=${outdir:-${ANALYSIS_ROOT:-${ANALYSIS_ROOT}}/ukb/assoc}
script_dir=${script_dir:-$SCRIPTS_LOCAL_ROOT}
assoc_R=${assoc_R:-$this_dir/f/assoc.R}
path_args=(--data-root "$dir0" --phenotype-dir "$phedir" --output-dir "$outdir" --script-dir "$script_dir" --r-script "$assoc_R")

if [[ $# -eq 0 ]]; then
	for module in forest phewas medwas; do
		"$0" "$module" "${path_args[@]}"
	done
	exit 0
fi
case "$1" in -h | --help)
	usage
	exit 0
	;;
esac

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
	forest | phewas | medwas) default_action=all ;;
	*)
		echo "ERROR: unknown METHOD: $method" >&2
		usage >&2
		exit 1
		;;
esac

action="$default_action"
if [[ $# -gt 0 && "$1" != --* ]]; then
	action="$1"
	shift
fi
case "$method:$action" in
	forest:all | forest:save_rds | forest:plot | phewas:all | phewas:save_rds | phewas:plot | medwas:all | medwas:save_rds | medwas:plot) ;;
	*)
		echo "ERROR: unsupported METHOD/ACTION: $method $action" >&2
		usage >&2
		exit 1
		;;
esac

# Expand a comma-separated --mode list into independent runs. Keeping each mode
# in its own R process also keeps checkpoints and result files isolated.
mode_value=""
args=("$@")
for ((i = 0; i < ${#args[@]}; i++)); do
	if [[ "${args[i]}" == "--mode" ]]; then
		if ((i + 1 >= ${#args[@]})); then
			echo "ERROR: --mode requires a value" >&2
			exit 1
		fi
		mode_value="${args[i + 1]}"
		break
	fi
done
if [[ "$mode_value" == *,* ]]; then
	IFS=',' read -ra modes <<<"$mode_value"
	for mode_one in "${modes[@]}"; do
		mode_one="${mode_one//[[:space:]]/}"
		[[ -n "$mode_one" ]] || {
			echo "ERROR: --mode contains an empty entry: $mode_value" >&2
			exit 1
		}
		run_args=("${args[@]}")
		for ((j = 0; j < ${#run_args[@]}; j++)); do
			if [[ "${run_args[j]}" == "--mode" ]]; then
				run_args[j + 1]="$mode_one"
				break
			fi
		done
		"$0" "$method" "$action" "${run_args[@]}" "${path_args[@]}"
	done
	exit 0
fi

command -v Rscript >/dev/null 2>&1 || {
	echo "ERROR: Rscript not found" >&2
	exit 1
}
[[ -f "$assoc_R" ]] || {
	echo "ERROR: f/assoc.R not found: $assoc_R" >&2
	exit 1
}

echo "========== RUN ASSOC: METHOD=$method ACTION=$action =========="
DIR0="$dir0" PHEDIR="$phedir" UKB_OUT="$outdir" SCRIPT_DIR="$script_dir" \
	Rscript "$assoc_R" --method "$method" --action "$action" "$@"
echo "========== DONE ASSOC: METHOD=$method ACTION=$action =========="
