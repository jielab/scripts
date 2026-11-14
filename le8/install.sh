#!/usr/bin/env bash

# Display standalone help before loading runtimes or scheduling work.
if [[ ${BASH_SOURCE[0]} == "$0" && ( ${1:-} == -h || ${1:-} == --help ) ]]; then
    cat <<'PUBLIC_ENTRY_HELP'
Modules:
    --profile report|analysis|deep|validation  Dependency profile [analysis]
    --env le8|cigma  Create/update a complete Conda environment

Main options:
    --submit auto|local|hpc --jobs-N N --wait TRUE|FALSE  Scheduled installation
    --plan-only --plan-dir DIR --queue NAME --hours N  Reviewable .cmd plan
    --python FILE --r-bin FILE  Interpreters
    --no-r  Skip R packages
    --print-environment le8|cigma  Print environment recipe
    --check  Check installed Python/R dependencies without data or downloads
    -h, --help  Show options and examples below

Examples:
    ./install.sh --profile analysis --check
    ./install.sh --profile analysis
    # Alternative: create/update the full Conda environment.
    ./install.sh --env le8
PUBLIC_ENTRY_HELP
    exit 0
fi

# Dependency entry point for LE8; install is a separate scheduled module.
# CIGMA uses a separate environment because its pandas/scipy requirements differ.
set -euo pipefail
ROOT=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
LE8_ROOT=$ROOT
source "$ROOT/../0f/workflow_env.sh"
workflow_init le8 "$@"
set -- "${WORKFLOW_ARGS[@]}"
source "$ROOT/f/runtime_env.sh"
if [[ -z ${CONDA_EXE:-} && -x $OS_GU_CONDA/bin/conda ]]; then
    export CONDA_EXE=$OS_GU_CONDA/bin/conda
fi
if [[ ${WORKFLOW_WORKER:-0} != 1 ]]; then
    inspect_only=FALSE
    for option in "$@"; do
        case "$option" in -h|--help|--print-environment) inspect_only=TRUE ;; esac
    done
    if [[ $inspect_only != TRUE ]]; then
        # The environment being created need not contain Python yet.
        export WORKFLOW_PYTHON=${LE8_PLANNER_PYTHON:-python3}
        exec "${LE8_PLANNER_PYTHON:-python3}" "$ROOT/f/submit.py" install "$@"
    fi
fi

# Standard Conda YAML, printed on demand; pip/cran comments select light installs.
environment_recipe() {
    case "$1" in
    le8) cat "$ROOT/env.yml"
        ;;
    cigma) cat <<'CIGMA_ENV'
# Full upstream CIGMA environment for new installations. The prepared-matrix
# adapter has also been numerically tested in the local fit-only runtime.
name: le8-cigma
channels:
  - conda-forge
dependencies:
  - python=3.12
  - pip
  - r-base
  - r-optparse
  - r-numderiv
  - r-matrix
  - pip:
      - cigma==1.1.0
CIGMA_ENV
        ;;
    *) echo "ERROR: environment must be le8 or cigma" >&2; return 2 ;;
    esac
}

profile=analysis
mode=pip
environment_name=le8
python_bin=${LE8_REPORT_PYTHON:-${PYTHON_BIN:-${HOME}/anaconda3/envs/le8/bin/python3}}
r_bin=${R_BIN:-Rscript}
install_r=TRUE
while (($#)); do
    case "$1" in
    --check) mode=check; shift ;;
    --profile) profile=${2:?--profile requires report, analysis, deep or validation}; shift 2 ;;
    --env) mode=conda; environment_name=${2:?--env requires le8 or cigma}; shift 2 ;;
    --print-environment) mode=print; environment_name=${2:?--print-environment requires le8 or cigma}; shift 2 ;;
    --python) python_bin=${2:?--python requires FILE}; shift 2 ;;
    --r-bin) r_bin=${2:?--r-bin requires FILE}; shift 2 ;;
    --no-r) install_r=FALSE; shift ;;
    -h|--help)
	bash "${BASH_SOURCE[0]}" --help
        exit 0 ;;
    *) echo "ERROR: Unknown installation option: $1" >&2; exit 2 ;;
    esac
done
case "$profile" in report|analysis|deep|validation) ;; *) echo 'ERROR: --profile requires report, analysis, deep or validation' >&2; exit 2 ;; esac
if [[ $mode == print ]]; then environment_recipe "$environment_name"; exit 0; fi
requirement_file=$(mktemp)
trap 'rm -f -- "$requirement_file"' EXIT
run() {
    "$@"
}
if [[ $mode == check ]]; then
    "$python_bin" - "$profile" <<'PY_CHECK'
import importlib, sys
modules = ["numpy", "scipy", "pandas", "matplotlib", "openpyxl"]
if sys.argv[1] != "report":
    modules += ["sklearn", "joblib", "threadpoolctl", "lightgbm", "pyreadr", "pyarrow", "bitarray", "duckdb"]
if sys.argv[1] == "deep":
    modules.append("torch")
if sys.argv[1] == "validation":
    modules.append("pytest")
for name in modules:
    importlib.import_module(name)
print("LE8 Python dependencies PASS: " + sys.argv[1])
PY_CHECK
    "$python_bin" -m pip check
    if [[ $install_r == TRUE ]]; then
        "$r_bin" --vanilla -e 'p<-c("data.table","ggplot2","openxlsx","brglm2","doParallel","callr","glmnet","coloc","susieR","clusterProfiler","org.Hs.eg.db","GO.db","DANDELION","shiny","jsonlite","digest"); m<-p[!vapply(p,requireNamespace,logical(1),quietly=TRUE)]; if(length(m)) stop("missing R packages: ",paste(m,collapse=", ")); cat("LE8 R dependencies PASS\n")'
    fi
    exit 0
fi
if [[ $mode == conda ]]; then
    environment_recipe "$environment_name" > "$requirement_file"
    [[ $environment_name != cigma ]] || environment_name=le8-cigma
    conda_bin=${CONDA_EXE:-conda}
    run "$conda_bin" env update --name "$environment_name" --file "$requirement_file"
    if [[ $environment_name == le8 ]]; then
        run "$conda_bin" run -n le8 Rscript --vanilla -e 'if(!requireNamespace("DANDELION",quietly=TRUE)) install.packages("https://cloud.r-project.org/src/contrib/DANDELION_0.1.0.tar.gz",repos=NULL,type="source")'
    fi
    exit 0
fi
if ! command -v "$python_bin" >/dev/null 2>&1; then
    echo "ERROR: selected Python does not exist: $python_bin. Create its environment or supply --python." >&2
    exit 2
fi
environment_recipe le8 | awk -v profile="$profile" '
  / # pip:/ {
    split($0, fields, " # pip:"); group=fields[2]; package=fields[1]
    sub(/^[[:space:]]*-[[:space:]]*/, "", package)
    if (package ~ /^[A-Za-z0-9_.-]+=[^=]/) sub(/=/, "==", package)
    if (group=="report" || (profile!="report" && group=="abm") ||
        (profile=="deep" && group=="deep") ||
        (profile=="validation" && group=="validation")) print package
  }
' > "$requirement_file"
if [[ $profile == deep ]]; then
    run "$python_bin" -m pip install --extra-index-url https://download.pytorch.org/whl/cu126 -r "$requirement_file"
else
    run "$python_bin" -m pip install -r "$requirement_file"
fi
if [[ $install_r == TRUE ]]; then
    mapfile -t r_packages < <(environment_recipe le8 | sed -n 's/.* # cran://p')
    run "$r_bin" - "${r_packages[@]}" <<'R_INSTALL'
# Explicit dependency installation; app.R never installs software.
p <- commandArgs(trailingOnly=TRUE)
m <- p[!vapply(p, requireNamespace, logical(1), quietly=TRUE)]
if (length(m)) install.packages(m, repos="https://cloud.r-project.org")
if(!requireNamespace("DANDELION",quietly=TRUE)) install.packages("https://cloud.r-project.org/src/contrib/DANDELION_0.1.0.tar.gz",repos=NULL,type="source")
R_INSTALL
    # Bioconductor annotation packages are not available from CRAN. Install
    # into the same library selected by runtime_env.sh on local/HPC machines.
    mapfile -t bioc_packages < <(environment_recipe le8 | sed -n 's/.* # bioc://p')
    run "$r_bin" --vanilla - "${bioc_packages[@]}" <<'BIOC_INSTALL'
p <- commandArgs(trailingOnly=TRUE)
m <- p[!vapply(p, requireNamespace, logical(1), quietly=TRUE)]
if (length(m)) {
    if (!requireNamespace("BiocManager", quietly=TRUE))
        install.packages("BiocManager", repos="https://cloud.r-project.org")
    BiocManager::install(m, ask=FALSE, update=FALSE)
}
m <- p[!vapply(p, requireNamespace, logical(1), quietly=TRUE)]
if (length(m)) stop("Missing Bioconductor packages after installation: ", paste(m, collapse=", "))
cat("LE8 Bioconductor dependencies PASS\n")
BIOC_INSTALL
fi
