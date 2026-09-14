# Shared console policy for user-facing analysis entry points.
# Screen output is event-driven: major stages and diagnostics, never timer reports.
# Internal workers inherit the marker and keep their existing output contracts.
if [[ ${SCRIPT_CONSOLE_ACTIVE:-0} != 1 ]]; then
  _console_bypass=0
  for _console_arg in "$@"; do
    case "$_console_arg" in
      -h|--help|--help-all|help|--dry-run|shiny)
        _console_bypass=1 ;;
    esac
  done
  [[ $# -gt 0 ]] || _console_bypass=1
  # An explicitly detached job must not leave this foreground logger waiting.
  case " $* " in *' --foreground FALSE '*|*' --foreground false '*) _console_bypass=1;; esac
  if [[ $_console_bypass == 0 && ${SCRIPT_VERBOSE:-0} != 1 ]]; then
    exec python3 /mnt/d/scripts/0f/console_run.py --script "${BASH_SOURCE[1]}" -- "$@"
  fi
  unset _console_bypass _console_arg
fi
