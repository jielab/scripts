#!/usr/bin/env bash

# Display standalone help before loading runtimes or scheduling work.
if [[ ${BASH_SOURCE[0]} == "$0" && ( ${1:-} == -h || ${1:-} == --help ) ]]; then
    cat <<'PUBLIC_ENTRY_HELP'
Modules:
    Library  Shared functions/configuration; source from the owning entry script

Main options:
    No standalone analysis options. Read the variables/functions below.
    -h, --help  Show this library help

Examples:
    source ./console.sh
PUBLIC_ENTRY_HELP
    exit 0
fi

# Direct CLI help must return before setup, downloads or analysis.
if [[ "${BASH_SOURCE[0]}" == "$0" ]]; then
	case "${1:-}" in
		-h | --help)
			bash "${BASH_SOURCE[0]}" --help
			exit 0
			;;
	esac
fi


# Load the central runtime directory when sourced without an entry wrapper.
if [[ -z ${DIR_RUNTIME_TMP:-} ]]; then
    _console_paths=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)
    while [[ ! -f "$_console_paths/0f/0os.sh" ]]; do
        [[ $_console_paths != / ]] || { echo "Cannot locate 0f/0os.sh" >&2; return 2; }
        _console_paths=$(dirname -- "$_console_paths")
    done
    source "$_console_paths/0f/0os.sh"
    unset _console_paths
fi
# 🚩 Runtime temporary files
export TMPDIR="${TMPDIR:-$DIR_RUNTIME_TMP}" TMP="${TMP:-$DIR_RUNTIME_TMP}" TEMP="${TEMP:-$DIR_RUNTIME_TMP}"
export PYTHONDONTWRITEBYTECODE=1 PYTHONPYCACHEPREFIX=${DIR_RUNTIME_TMP}/python-cache


# 🚩 Console logging
# Shared console policy for user-facing analysis entry points.
# Screen output is event-driven: major stages and diagnostics, never timer reports.
# Internal workers inherit the marker and keep their existing output contracts.
if [[ ${SCRIPT_CONSOLE_ACTIVE:-0} != 1 ]]; then
	_console_bypass=0
	for _console_arg in "$@"; do
		case "$_console_arg" in
		-h | --help | --help-all | help | --dry-run | shiny)
			_console_bypass=1
			;;
		esac
	done
	[[ $# -gt 0 ]] || _console_bypass=1
	# Only an explicit background request may escape the launcher lifetime.
	_console_detach=0
	case " $* " in
		*' --foreground FALSE '* | *' --foreground false '* | *' --foreground=FALSE '* | *' --foreground=false '*) _console_detach=1 ;;
	esac
	if [[ $_console_detach == 0 ]]; then
		export SCRIPT_FOREGROUND=1
		if [[ $_console_bypass == 1 || ${SCRIPT_VERBOSE:-0} == 1 || ${SCRIPT_CONSOLE_PASSTHROUGH:-0} == 1 ]]; then
			exec python3 "$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)/console_run.py" --passthrough --script "${BASH_SOURCE[1]}" -- "$@"
		fi
		exec python3 "$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)/console_run.py" --script "${BASH_SOURCE[1]}" -- "$@"
	fi
	unset _console_bypass _console_arg _console_detach
fi
