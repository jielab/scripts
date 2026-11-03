#!/usr/bin/env bash
set -euo pipefail
HERE="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
if [[ ${1:-} == --software-only ]]; then
    shift
    source "$HERE/f/env.sh"
    workflow_args "$@"
    set -- "${WORKFLOW_ARGS[@]}"
    exec "$METAG_PYTHON" "$HERE/f/install_tools.py" "$@"
fi
case "${1:-}" in
    -h|--help|help)
        cat <<'HELP'
Usage: bash install.sh [setup] [options]
Install the pinned CPU tools through the local / LSF setup workflow.

  --software-only     Install/verify tools without reference or cohort data (put first)
  --config FILE       Runtime paths/resources (defaults: hpc.yml)
  --submit FALSE      Generate setup commands without running or submitting them
  --submit TRUE       Run/submit the generated setup commands
  submit=local|hpc     Environment variable selecting execution mode
  jobs_N=N            Environment variable controlling local concurrency

Examples:
  bash install.sh --config hpc.yml --submit FALSE
  submit=local bash install.sh --config hpc.yml --submit TRUE
  submit=hpc bash install.sh --config hpc.yml --submit TRUE

Configuration:
  software.yml  Source versions, archives and SHA256 checksums
  hpc.yml              Tool/output paths and job resources
  f/env.sh             Interpreter and site compiler environment
Requires Python >= 3.10 with PyYAML; gcc, g++, make, Perl and Bash on workers.
Bundled Kraken sources are compiled with project native extensions; optional
Kun-Peng uses its pinned release binary. Existing installations are verified.
Licenses: Kraken GPL-3.0 in its source archive; Kun-Peng in LICENSE.kunpeng.
README files are optional documentation and are never installation inputs.
HELP
        exit 0 ;;
esac
if [[ $# == 0 || ${1:0:1} == - ]]; then
    set -- setup "$@"
fi
exec bash "$HERE/f/entry.sh" "$@"
