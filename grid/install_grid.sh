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
    Default  Compatibility installer; delegates to install.sh

Main options:
    --dry-run  Print installation commands
    -h, --help  Show help; dependency details: 0.README.md

Example:
    ./install_grid.sh --dry-run
    ./install_grid.sh
HELP
			exit 0
			;;
	esac
fi
# Compatibility entry; canonical installer: ./install.sh
set -euo pipefail
ROOT=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)
exec bash "$ROOT/install.sh" "$@"
