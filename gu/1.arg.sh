#!/usr/bin/env bash
source "$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/../0f" && pwd)/output_dirs.sh"
script_output_directory_args "${BASH_SOURCE[0]}" "$@" || { return 2 2>/dev/null || exit 2; }
set -- "${SCRIPT_OUTPUT_DIRECTORY_ARGS[@]}"

# Public help returns before runtime setup, input access or job submission.
if [[ ${BASH_SOURCE[0]} == "$0" ]]; then
    case "${1:-}" in
        ""|-h|--help|help)
            cat <<'PUBLIC_ENTRY_HELP'
GU: build and validate ancestral recombination graphs.

Usage:
  ./1.arg.sh run|MODULE [options]

Modules:
  run | build                      Build configured methods and chromosomes.
  prep_gen | check                 Prepare genotype inputs or validate ARG outputs.
  submit | status | retry RUN_DIR  Manage an existing immutable job plan.

Permanent output directories (legacy options remain accepted):
  --dir-out-arg DIR
                                  Native ARG products [${GU_REF_ROOT}/TARGET/GRCH/arg.singer; HPC: analysis/gu/arg].
  --dir-out-gen4arg DIR
                                  Prepared ARG inputs [VCF_DIR.4arg.singer; HPC: analysis/gu/arg].

Options:
  --target NAME --grch 37|38       Dataset/build [1kg / 37]
  --dir-target DIR                 Default: OS data root/gen/<target>/<grch>
  --method LIST --format NAME      Methods: tsinfer,singer [tsinfer / native]
  --cores N                        Shared budget: 0f/0os.sh [HPC 40 / local 16]
  --threads N --jobs N             Default to the shared budget; SINGER also caps workers by memory
  --chr 1-22,X                     Chromosome selection [1-22,X]
  --replace TRUE|FALSE             Reuse verified outputs/checkpoints [FALSE]
  --dir-arg DIR --dir-gen4arg DIR  Isolated ARG and preparation outputs
  --resource-guard auto|systemd|scheduler|ulimit Resource control [auto]
  --memory-limit-gb N --memory-swap-gb N Local limits [32 / 4]
  --queue NAME --slots N --walltime HHH:MM LSF resources (slots >= threads)
  --run-dir DIR                    Optional explicit plan directory; normally managed automatically
  --plan-only                      Write job files without submitting; also accepted by submit/retry
  -h, --help                       Show this help.

Examples:
  ./1.arg.sh run --method tsinfer --target 1kg --grch 37 --chr 20 --plan-only
  ./1.arg.sh run --method tsinfer,singer --chr 1,20,X --format trace --queue 2t50c --memory-limit-gb 256
  ./1.arg.sh status RUN_DIR
PUBLIC_ENTRY_HELP
            exit 0 ;;
        all) printf 'Use run for the complete workflow; see -h.\n' >&2; exit 2 ;;
    esac
fi
if [[ ${1:-} == run ]]; then shift; set -- build "$@"; fi

set -euo pipefail
ROOT=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)
source "$ROOT/f/scratch.sh"
source "$ROOT/f/0.common.sh"
gu_load_config
arg_help=0
for arg_option in "$@"; do
	case "$arg_option" in -h|--help|--help-all) arg_help=1 ;; esac
done
if [[ $arg_help == 0 && ${GU_ARG_BATCH_WORKER:-0} != 1 ]]; then
    case "${1:-}" in
        status|submit|retry)
            exec "${GU_SUBMIT_PYTHON:-$SHARED_PYTHON}" "$ROOT/f/arg_batch.py" "$@" ;;
        build|prep_gen)
            arg_plan_only=0
            for arg_option in "$@"; do [[ $arg_option != --plan-only ]] || arg_plan_only=1; done
            if [[ $arg_plan_only == 1 || ( $submit == hpc && -z ${LSB_JOBID:-} ) ]]; then
                exec "${GU_SUBMIT_PYTHON:-$SHARED_PYTHON}" "$ROOT/f/arg_batch.py" "$@"
            fi
            unset arg_plan_only
            ;;
    esac
fi
unset arg_help arg_option
source "$ROOT/f/console.sh"


# Command-line help
# Keep public ARG examples together.
usage() {
	bash "${BASH_SOURCE[0]}" --help
}

command=${1:-help}
case "$command" in
	-h | --help | help)
		usage
		exit 0
		;;
	--help-all)
		usage
		exit 0
		;;
esac
case "$command" in
	prep_gen | build | check) shift ;;
	*)
		echo "ERROR: unknown command: $command (use ./1.arg.sh -h)" >&2
		exit 2
		;;
esac
for option in "$@"; do
	case "$option" in
		-h | --help)
			usage
			exit 0
			;;
		--help-all)
			usage
			exit 0
			;;
	esac
done
export REFGEN_GEN4ARG_ONLY=0
source "$ROOT/f/arg.sh"
arg_run_guarded "$command" "$@"
