#!/usr/bin/env bash
# BEGIN QIMING RESOURCE EXAMPLES (Admin, 2026-11-05)
# Queue 33: 40 cores / 192 GiB per node. Reserve 10/20/30/40 cores by
# rounding max(required_cores/40, required_memory_GiB/192) UP to 25% steps.
# Example for a 24-core, 160-GiB CPU worker:
#   bsub -q 33 -n 40 -W 24:00 -M 160GB -R "span[hosts=1] rusage[mem=160GB/host]" ...
# A small <=10-core, <=48-GiB job uses -q 33 -n 10; request its actual RAM.
# For exact CPU counts, use Admin-suggested v4-3t or v3-6t after queueinfo.
# The examples are allocation guidance, not measured peak RAM for every stage.
# CPU planning flags (also supported by the install dispatcher):
#   --submit hpc --queue 33 --cores 40 --memory-limit-gb 160 --hours 24 --plan-only
# New CPU plans auto-round slots. GPU stages keep their separate --gpu-queue resources.
# Details and queue limits: ../0f/HPC_SUBMISSION_POLICY.md
# END QIMING RESOURCE EXAMPLES


# Public directory options use --dir-XXX; legacy spellings remain compatible.
if [[ "${BASH_SOURCE[0]}" == "$0" ]]; then
    source "$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/../0f" && pwd)/cli_dirs.sh"
    script_directory_args "$(realpath -- "${BASH_SOURCE[0]}")" "$@"
    set -- "${SCRIPT_DIRECTORY_ARGS[@]}"
    unset SCRIPT_DIRECTORY_ARGS
fi
# Dependency entry point for LE8. This installer is never called by le8.sh.
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
        exec "${LE8_PLANNER_PYTHON:-python3}" "$ROOT/f/hpc_workflow.py" install "$@"
    fi
fi

# Standard Conda YAML, printed on demand; pip/cran comments select light installs.
environment_recipe() {
    case "$1" in
    le8) cat "$ROOT/environment.yml"
        ;;
    cigma|le8-cigma) cat <<'CIGMA_ENV'
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

profile=report
mode=pip
environment_name=le8
dry_run=FALSE
python_bin=${LE8_REPORT_PYTHON:-${PYTHON_BIN:-${HOME}/anaconda3/envs/le8/bin/python3}}
r_bin=${R_BIN:-Rscript}
install_r=TRUE
while (($#)); do
    case "$1" in
    --check) mode=check; shift ;;
    --report) profile=report; shift ;;
    --abm) profile=abm; shift ;;
    --all) profile=abm; shift ;;
    --validation) profile=validation; shift ;;
    --env) mode=conda; environment_name=${2:?--env requires le8 or cigma}; shift 2 ;;
    --print-environment) mode=print; environment_name=${2:?--print-environment requires le8 or cigma}; shift 2 ;;
    --python) python_bin=${2:?--python requires FILE}; shift 2 ;;
    --r-bin) r_bin=${2:?--r-bin requires FILE}; shift 2 ;;
    --no-r) install_r=FALSE; shift ;;
    --dry-run) dry_run=TRUE; shift ;;
    -h|--help)
	cat <<'HELP'
Modules:
    --report / --abm / --all / --validation  Dependency profiles
    --env le8|cigma  Create/update a complete Conda environment

Main options:
    --submit local|hpc --jobs-N N --wait TRUE|FALSE  Scheduled installation
    --plan-only --plan-dir DIR --queue NAME --hours N  Reviewable .cmd plan
    --python FILE --r-bin FILE  Interpreters
    --no-r  Skip R packages
    --print-environment le8|cigma  Print environment recipe
    --check  Check installed Python/R dependencies without data or downloads
    --dry-run  Print commands without installing
    -h, --help  Show options and examples below

Example:
    source "${SCRIPTS_ROOT}/0f/0os.sh"
    # Dependency groups (choose the group needed for this installation):
    ./install.sh --report
    ./install.sh --abm
    ./install.sh --all
    ./install.sh --validation
    ./install.sh --print-environment le8
    ./install.sh --print-environment cigma

    ./install.sh --env le8
    ./install.sh --env cigma --dry-run
HELP
        exit 0 ;;
    *) echo "ERROR: Unknown installation option: $1" >&2; exit 2 ;;
    esac
done
if [[ $mode == print ]]; then environment_recipe "$environment_name"; exit 0; fi
requirement_file=$(mktemp)
trap 'rm -f -- "$requirement_file"' EXIT
run() {
    if [[ $dry_run == TRUE ]]; then printf '%q ' "$@"; printf '\n'; else "$@"; fi
}
if [[ $mode == check ]]; then
    "$python_bin" -c 'import numpy,scipy,pandas,sklearn,lightgbm,pyreadr,pyarrow,bitarray,duckdb,torch; print("LE8 Python dependencies PASS")'
    "$python_bin" -m pip check
    "$r_bin" --vanilla -e 'p<-c("data.table","ggplot2","openxlsx","brglm2","doParallel","callr","glmnet","coloc","susieR","clusterProfiler","org.Hs.eg.db","GO.db","DANDELION","shiny","jsonlite","digest"); m<-p[!vapply(p,requireNamespace,logical(1),quietly=TRUE)]; if(length(m)) stop("missing R packages: ",paste(m,collapse=", ")); cat("LE8 R dependencies PASS\n")'
    exit 0
fi
if [[ $mode == conda ]]; then
    environment_recipe "$environment_name" > "$requirement_file"
    [[ $environment_name != cigma ]] || environment_name=le8-cigma
    if [[ $dry_run == TRUE ]]; then cat "$requirement_file"; fi
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
        (profile=="validation" && group=="validation")) print package
  }
' > "$requirement_file"
if [[ $dry_run == TRUE ]]; then cat "$requirement_file"; fi
run "$python_bin" -m pip install --extra-index-url https://download.pytorch.org/whl/cu126 -r "$requirement_file"
if [[ $install_r == TRUE ]]; then
    mapfile -t r_packages < <(environment_recipe le8 | sed -n 's/.* # cran://p')
    run "$r_bin" - "${r_packages[@]}" <<'R_INSTALL'
# Explicit dependency installation; app.R never installs software.
p <- commandArgs(trailingOnly=TRUE)
m <- p[!vapply(p, requireNamespace, logical(1), quietly=TRUE)]
if (length(m)) install.packages(m, repos="https://cloud.r-project.org")
if(!requireNamespace("DANDELION",quietly=TRUE)) install.packages("https://cloud.r-project.org/src/contrib/DANDELION_0.1.0.tar.gz",repos=NULL,type="source")
R_INSTALL
fi
