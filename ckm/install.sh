#!/usr/bin/env bash

# Public help returns before runtime setup, input access or job submission.
if [[ ${BASH_SOURCE[0]} == "$0" ]]; then
    case "${1:-}" in
        ""|-h|--help|help)
            cat <<'PUBLIC_ENTRY_HELP'
CKM: check or install analysis dependencies.

Usage:
  ./install.sh run [options]

Modules:
  run                              Check the configured dependencies.

Options:
  --check                          Check installed R/Python dependencies.
  --r                              Install missing R dependencies.
  --conda PREFIX                   Create or update a Conda environment.
  -h, --help                       Show this help.

Examples:
  ./install.sh run
  ./install.sh --check
PUBLIC_ENTRY_HELP
            exit 0 ;;
        all) printf 'Use run for the complete workflow; see -h.\n' >&2; exit 2 ;;
    esac
fi
if [[ ${1:-} == run ]]; then shift; set -- --check "$@"; fi
# Install separately. Analysis entries never install or update dependencies.
set -euo pipefail
ckm_install_root=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)
source "$ckm_install_root/../0f/0os.sh"
mode=${1:---check}
case "$mode" in
    --check)
        rscript=${CKM_RSCRIPT:-Rscript}
        "$rscript" -e 'for(p in c("data.table","survival","jsonlite")){if(!requireNamespace(p,quietly=TRUE))stop("Missing package: ",p);cat(p,as.character(packageVersion(p)),"\n")};cat(R.version.string,"\n")'
        "${CKM_PYTHON:-python3}" --version
        ;;
    --r)
        "${CKM_RSCRIPT:-Rscript}" -e 'p<-c("data.table","survival","jsonlite");m<-p[!vapply(p,requireNamespace,logical(1),quietly=TRUE)];if(length(m))install.packages(m,repos="https://cloud.r-project.org");for(x in p)stopifnot(requireNamespace(x,quietly=TRUE))'
        ;;
    --conda)
        prefix=${2:?Usage: install.sh --conda /absolute/path/to/envs/ckm}
        [[ $prefix == /* ]] || { echo 'Use an absolute environment prefix.' >&2; exit 2; }
        if command -v mamba >/dev/null; then manager=mamba
        elif command -v conda >/dev/null; then manager=conda
        else echo 'conda or mamba is required for --conda; use --r with an existing R instead.' >&2; exit 2; fi
        if [[ -d $prefix/conda-meta ]]; then
            "$manager" env update --prefix "$prefix" --file "$ckm_install_root/env.yml"
        else
            "$manager" env create --prefix "$prefix" --file "$ckm_install_root/env.yml"
        fi
        CKM_RSCRIPT=$prefix/bin/Rscript CKM_PYTHON=$prefix/bin/python3 bash "$ckm_install_root/install.sh" --check
        ;;
    *)
        echo 'Usage: install.sh --check | --r | --conda /absolute/env/prefix' >&2
        exit 2 ;;
esac
