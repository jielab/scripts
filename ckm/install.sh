#!/usr/bin/env bash
# Install separately. Analysis entries never install or update dependencies.
set -euo pipefail
ckm_install_root=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)
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
