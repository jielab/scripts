#!/usr/bin/env bash

# Omitted action uses this entry's configured workflow; explicit actions are unchanged.
if [[ ${BASH_SOURCE[0]} == "$0" ]]; then
    if (($# == 0)); then
        set -- build
    elif [[ $1 == --* && $1 != --help && $1 != --help-full && $1 != --help-all && $1 != --execute && $1 != --worker ]]; then
        set -- build "$@"
    fi
fi

# Display standalone help before loading runtimes or scheduling work.
if [[ ${BASH_SOURCE[0]} == "$0" && ( ${1:-} == -h || ${1:-} == --help ) ]]; then
    cat <<'PUBLIC_ENTRY_HELP'
./arg.sh  # 默认 build；构建配置的方法及染色体
./arg.sh MODULE  # build|prep_gen|check|submit|status|retry
./arg.sh --method tsinfer --target 1kg --grch 37

Modules:
    build / check  Build and validate ARG outputs
    submit RUN_DIR  Submit unsubmitted jobs from a saved plan
    status RUN_DIR  Show each method/chromosome job, state and elapsed seconds
    retry RUN_DIR   Resubmit only confirmed failed jobs; reuse checkpoints

Main options:
    --target NAME --grch 37|38  Dataset/build [1kg / 37]
    --dir-target DIR  Default: OS data root/gen/<target>/<grch>
                     --dir-gen remains a compatibility alias
    --method LIST --format NAME  Methods: tsinfer,singer [tsinfer / native]
                                HPC login: one independent LSF job per method x chromosome
    --cores N            Shared budget: 0f/0os.sh [HPC 40 / local 16]
    --threads N --jobs N  Default to the shared budget; SINGER also caps workers by memory
                         tsinfer: threads per chromosome
                         SINGER: jobs windows, threads total CPU budget
    --chr 1-22,X  Chromosome selection [1-22,X]
    --replace TRUE|FALSE  Reuse verified outputs/checkpoints [FALSE]
    --dir-arg DIR --dir-gen4arg DIR  Isolated ARG and preparation outputs
    --resource-guard auto|systemd|scheduler|ulimit  Resource control [auto]
    --memory-limit-gb N --memory-swap-gb N  Local limits [32 / 4]
    --queue NAME --slots N --walltime HHH:MM  LSF resources (slots >= threads)
    --run-dir DIR  Optional explicit plan directory; normally managed automatically
    --plan-only    Write job files without submitting; also accepted by submit/retry
    -h, --help  Show options and examples below

Current HPC 1kg GRCh37 restart (copy; reuse completed outputs and checkpoints):
    cd /work/sph-huangj/scripts/gu
    bash arg.sh build --target 1kg --grch 37 --method tsinfer,singer --chr 1,20,X --format trace --queue 2t50c --memory-limit-gb 256
    # 256 GiB retains this analysis's original memory budget; CPU defaults are shared.

HPC build/prep_gen on a login node submits method x chromosome jobs automatically.
Each task has isolated output/preparation directories, logs and attempt timings.
Default results: analysis/gu/arg/<target>/<build>/<method>/chr<chromosome>/.
Prepared data and job metadata use .prepared/ and .jobs/ inside the target/build directory.
Repeat build to verify/reuse completed outputs and retry failed tasks; active jobs are skipped.
Use --replace TRUE to recompute outputs. No manual submit/retry command is normally needed.
tsinfer checkpoints live under each output directory's .work/ and survive failures.
--dir-arg/--dir-gen4arg are batch bases; tasks use BASE/<method>/chr<chromosome>.
--jobs controls work INSIDE each job (SINGER windows), not the number submitted.
Local/allocated runs retain direct execution. --plan-only always creates an LSF plan.
Saved job plans record exact output paths, commands and resources.
For TRACE, each task's exported ARG is OUTPUT/trace/<method>; pass that directory
with the matching --chr to gu.sh trace. No multi-chromosome merge job is submitted.

Examples:
    ./install.sh --arg
    ./arg.sh build --grch 37 --target 1kg --chr 1,20,X --method tsinfer,singer --format trace --threads 8 --jobs 8 --memory-limit-gb 256 --queue 2t50c
    # SINGER: up to 8 windows/job; 256 GiB covers 8 x 24 GiB plus headroom.
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
