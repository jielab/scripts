#!/usr/bin/env bash

# Public directory options use --dir-XXX; legacy spellings remain compatible.
if [[ "${BASH_SOURCE[0]}" == "$0" ]]; then
    source /mnt/d/scripts/0f/cli_dirs.sh
    script_directory_args "$(realpath -- "${BASH_SOURCE[0]}")" "$@"
    set -- "${SCRIPT_DIRECTORY_ARGS[@]}"
    unset SCRIPT_DIRECTORY_ARGS
fi
source /mnt/d/scripts/0f/console.sh
set -euo pipefail
ROOT=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)


# 🚩 Command-line help
# Keep public ARG examples together.
usage() {
	cat <<'HELP'
Modules:
    build / check  Build and validate ARG outputs

Main options:
    --dir-gen DIR  Reference genotype root [/mnt/f/gen/1kg/GRCH]
    --method NAME --format NAME  Inference/output mode [needle / native]
    --threads N --jobs N  Resources [8 / 1]
    -h, --help  Help; details: 0.README.md

Example:
    cd /mnt/d/scripts/gu

    ./install.sh --arg
    ./arg.sh build --dir-gen /mnt/f/gen/1kg/37 --method threads --format trace --threads 8 --jobs 1
    ./arg.sh check --dir-gen /mnt/f/gen/1kg/37 --method threads --format trace

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
