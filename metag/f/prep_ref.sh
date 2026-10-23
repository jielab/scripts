#!/usr/bin/env bash

# Public directory options use --dir-XXX; legacy spellings remain compatible.
if [[ "${BASH_SOURCE[0]}" == "$0" ]]; then
    source /mnt/d/scripts/0f/cli_dirs.sh
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
    Default  Shared reference preparation and audit

Main options:
    --dir-metag DIR --dir-annot DIR  Dataset/reference paths
    --threads N --download-threads N --max-db-gb N  Resources
    --mask-jobs N --masker-threads N  Parallel genomes / masking threads [2 / 1]
    --prepare-memory-gb N  Preparation-only RAM ceiling [12 GiB]
    --dry-run  Show plan
    -h, --help  Show help; workflow details: ../0.README.md

Example:
    ./metag.sh prep_ref --dry-run
    ./metag.sh prep_ref
HELP
			exit 0
			;;
	esac
fi
# Same reference preparation implementation as the automatic workflow.
set -euo pipefail
HERE="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
exec "${METAG_PYTHON:-python3}" "$HERE/workflow.py" --reference-only "$@"
