#!/usr/bin/env bash

set -euo pipefail
ROOT=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)
source "$ROOT/f/scratch.sh"
source "$ROOT/f/0.common.sh"
gu_load_config
source "$ROOT/f/console.sh"


# 🚩 Command-line help
# Keep public ARG examples together.
usage() {
	cat <<'HELP'
Modules:
    build / check  Build and validate ARG outputs

Main options:
    --date-scratch YYYY-MM-DD  Scratch analysis date [2026-10-31]; pass the original date to resume
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
    Scratch defaults to /tmp/refgen_arg_scratch-<uid>/ and must remain on native Linux storage.
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
