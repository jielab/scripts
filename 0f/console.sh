#!/usr/bin/env bash

# Direct CLI help must return before setup, downloads or analysis.
if [[ "${BASH_SOURCE[0]}" == "$0" ]]; then
	case "${1:-}" in
		-h | --help)
			cat <<'HELP'
Modules:
    Library  Shared functions/configuration; source from the owning entry script

Main options:
    No standalone analysis options. Read the variables/functions below.
    -h, --help  Show this library help

Example:
    source ./console.sh
HELP
			exit 0
			;;
	esac
fi


# 🚩 Runtime temporary files
export TMPDIR=/tmp TMP=/tmp TEMP=/tmp
export PYTHONDONTWRITEBYTECODE=1 PYTHONPYCACHEPREFIX=/tmp/python-cache


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
			exec python3 /mnt/d/scripts/0f/console_run.py --passthrough --script "${BASH_SOURCE[1]}" -- "$@"
		fi
		exec python3 /mnt/d/scripts/0f/console_run.py --script "${BASH_SOURCE[1]}" -- "$@"
	fi
	unset _console_bypass _console_arg _console_detach
fi
