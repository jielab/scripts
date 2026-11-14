#!/usr/bin/env bash

# Display standalone help before loading runtimes or scheduling work.
if [[ ${BASH_SOURCE[0]} == "$0" && ( ${1:-} == -h || ${1:-} == --help ) ]]; then
    cat <<'PUBLIC_ENTRY_HELP'
Modules:
    build / check  Build and validate ARG outputs
    submit RUN_DIR  Submit unsubmitted jobs from a saved plan
    status RUN_DIR  Show each method/chromosome job, state and elapsed seconds
    retry RUN_DIR   Resubmit only confirmed failed jobs; reuse checkpoints

Main options:
    --target NAME --grch 37|38  Dataset/build [1kg / 37]
    --dir-target DIR  Default: OS data root/gen/<target>/<grch>
                     --dir-gen remains a compatibility alias
    --method LIST --format NAME  Methods: threads,needle,tsinfer,singer [needle / native]
                                HPC login: one independent LSF job per method x chromosome
    --threads N --jobs N  Resources [8 / 1]; Threads processes one chromosome at a time
                         Needle/tsinfer: jobs chromosomes, threads per chromosome
                         SINGER: jobs windows, threads total CPU budget
    --chr 1-22,X  Chromosome selection [1-22,X]
    --replace TRUE|FALSE  Reuse verified outputs/checkpoints [FALSE]
    --dir-arg DIR --dir-gen4arg DIR  Isolated ARG and preparation outputs
    --resource-guard auto|systemd|scheduler|ulimit  Resource control [auto]
    --memory-limit-gb N --memory-swap-gb N  Local limits [32 / 4]
    --queue NAME --slots N --walltime HHH:MM  LSF resources (slots >= threads)
    --run-dir DIR  New batch directory; defaults under OS analysis root
    --plan-only    Write job files without submitting; also accepted by submit/retry
    -h, --help  Show options and examples below

HPC build/prep_gen on a login node submits method x chromosome jobs automatically.
Each task has isolated output/preparation directories, logs and attempt timings.
--dir-arg/--dir-gen4arg are batch bases; tasks use BASE/<method>/chr<chromosome>.
--jobs controls work INSIDE each job (SINGER windows), not the number submitted.
Local/allocated runs retain direct execution. --plan-only always creates an LSF plan.
Saved RUN_DIR/plan.json records exact output paths, commands and resources.
For TRACE, each task's exported ARG is OUTPUT/trace/<method>; pass that directory
with the matching --chr to gu.sh trace. No multi-chromosome merge job is submitted.

Examples:
    ./install.sh --arg
    ./arg.sh build --grch 37 --target 1kg --chr 1,20,X --method threads,needle,tsinfer,singer --format trace --threads 8 --jobs 1
    ./arg.sh status RUN_DIR
    ./arg.sh retry RUN_DIR
    # Preview the build above by appending --plan-only; then ./arg.sh submit RUN_DIR
PUBLIC_ENTRY_HELP
    exit 0
fi

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
		echo "ERROR: unknown command: $command (use ./arg.sh -h)" >&2
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
