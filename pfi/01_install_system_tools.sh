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
    Default  Install system tools with apt (sudo required)

Main options:
    No analysis options; installs the packages listed below in this script.
    -h, --help  Show help

Example:
    ./01_install_system_tools.sh
HELP
			exit 0
			;;
	esac
fi


# 🚩 01_install_system_tools
set -euo pipefail

sudo apt-get update
sudo apt-get install -y wget curl aria2 unzip gzip tar libxml2-utils build-essential r-base r-base-dev
