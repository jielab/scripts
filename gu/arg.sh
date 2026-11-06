#!/usr/bin/env bash
# BEGIN QIMING RESOURCE EXAMPLES (Admin, 2026-11-05)
# Queue 33: 40 cores / 192 GiB per node. Reserve 10/20/30/40 cores by
# rounding max(required_cores/40, required_memory_GiB/192) UP to 25% steps.
# Example for a 24-core, 160-GiB CPU worker:
#   bsub -q 33 -n 40 -W 24:00 -M 160GB -R "span[hosts=1] rusage[mem=160GB/host]" ...
# A small <=10-core, <=48-GiB job uses -q 33 -n 10; request its actual RAM.
# For exact CPU counts, use Admin-suggested v4-3t or v3-6t after queueinfo.
# The examples are allocation guidance, not measured peak RAM for every stage.
# Use gu.sh hpc plan (or gu-hpc.sh plan), then review and submit its RUN_DIR.
# CPU plan resource flags: --queue 33 --threads 24 --memory-gb 160
# The worker is allocated 40 slots; threads is the scientific thread count.
# Details and queue limits: ../0f/HPC_SUBMISSION_POLICY.md
# END QIMING RESOURCE EXAMPLES


set -euo pipefail
ROOT=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)
source "$ROOT/f/scratch.sh"
source "$ROOT/f/0.common.sh"
gu_load_config
arg_help=0
for arg_option in "$@"; do
	case "$arg_option" in -h|--help|--help-all) arg_help=1 ;; esac
done
if [[ $submit == hpc && -z ${LSB_JOBID:-} && $arg_help == 0 ]]; then
	case "${1:-}" in
		build|prep_gen)
			echo 'ERROR: ARG construction in HPC mode requires an allocated LSF job. Use bash gu.sh hpc plan arg ... followed by bash gu.sh hpc submit RUN_DIR, or run this entry inside an allocated job.' >&2
			exit 2 ;;
	esac
fi
unset arg_help arg_option
source "$ROOT/f/console.sh"


# 🚩 Command-line help
# Keep public ARG examples together.
usage() {
	cat <<'HELP'
Modules:
    build / check  Build and validate ARG outputs

Main options:
    --date-scratch YYYY-MM-DD  Scratch analysis date [2026-11-03]; pass the original date to resume
    --dir-gen DIR  Reference genotype root [$GU_REF_ROOT/1kg/GRCH]
    --method NAME --format NAME  Inference/output mode [needle / native]
    --threads N --jobs N  Resources [8 / 1]
    --chr 1-22,X --grch 37|38
    --dir-arg DIR --dir-gen4arg DIR  Isolated ARG and preparation outputs
    --resource-guard auto|systemd|scheduler|ulimit  Resource control [auto]
    --memory-limit-gb N --memory-swap-gb N  Local limits [32 / 4]
    -h, --help  Show options and examples below

Example:
    cd gu

    ./install.sh --arg
    ./arg.sh build --dir-gen $GU_REF_ROOT/1kg/37 --method threads --format trace --threads 8 --jobs 1
    ./arg.sh check --dir-gen $GU_REF_ROOT/1kg/37 --method threads --format trace

    Threads keeps one chromosome job. Each build writes log/arg.timing.chrN.json
    with prepare/infer/convert/check times, cache hits, and sampled process-tree RSS.
    Scratch defaults to ${DIR_RUNTIME_TMP}/refgen_arg_scratch-<uid>/ and must remain on native Linux storage.
    No generic ARG GPU backend is used.

    ./arg.sh build --method threads --format trace --threads 4
    ./arg.sh build
HELP
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
