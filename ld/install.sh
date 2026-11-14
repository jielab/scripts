#!/usr/bin/env bash

# Display standalone help before loading runtimes or scheduling work.
if [[ ${BASH_SOURCE[0]} == "$0" && ( ${1:-} == -h || ${1:-} == --help ) ]]; then
    cat <<'PUBLIC_ENTRY_HELP'
LD dependencies

Main options:
    --check    Check installed tools and dependencies
    --install  Install the configured environment

Default: --check. --install creates the Conda environment from env.yml.

Examples:
    ./install.sh --check
    ./install.sh --install
PUBLIC_ENTRY_HELP
    exit 0
fi
set -euo pipefail
ROOT=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)
case ${1:---check} in
    --install)
        conda env create --prefix "${LD_ENV_PREFIX:-$HOME/.conda/envs/ld}" --file "$ROOT/env.yml"
        echo 'Tool defaults are in ld.sh; environment variables can override them. Run ./install.sh --check.' ;;
    --check)
        source "$ROOT/ld.sh"
        source "$ROOT/f/env.sh"
        "$PYTHON" -c 'import sys; assert sys.version_info >= (3,11)'
        "$PLINK2" --version
        "$LD_CXX" --version
        if [[ -d $ROOT/tests ]]; then
            exec "$PYTHON" -m unittest discover -s "$ROOT/tests" -v
        fi
        echo 'Tool checks passed; development tests are not installed.' ;;
    *) echo 'Usage: ./install.sh --install|--check' >&2; exit 2 ;;
esac
