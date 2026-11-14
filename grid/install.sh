#!/usr/bin/env bash

# Display standalone help before loading runtimes or scheduling work.
if [[ ${BASH_SOURCE[0]} == "$0" && ( ${1:-} == -h || ${1:-} == --help ) ]]; then
    cat <<'PUBLIC_ENTRY_HELP'
Install/check the shared GRID and PRSformer environment

Main options:
    --environment-file FILE       Conda recipe [env.yml]
    --runtime-check-run cpu|gpu|auto  Runtime check processor [cpu]
    --check                       Check the existing environment
    --submit hpc|local|auto        Machine selection [auto]
    --plan-only                    Preview allocation and commands

GPU checks require an allocated GPU. Installation does not train models.

Examples:
    ./install.sh --check
    ./install.sh --check --runtime-check-run gpu
    # Install or update the environment.
    ./install.sh
PUBLIC_ENTRY_HELP
    exit 0
fi

# Single installer for the shared GRID / PRSformer environment and runtime checks.
set -euo pipefail
umask 077

usage() {
	bash "${BASH_SOURCE[0]}" --help
}

ROOT=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)
source "$ROOT/f/workflow_entry.sh"
grid_workflow_entry install "$@"
set -- "${GRID_ENTRY_ARGS[@]}"
DRY_RUN=FALSE
recipe=${GRID_ENV_RECIPE:-$ROOT/env.yml}
while (($#)); do
	case $1 in
		-h|--help|help) usage; exit 0 ;;
		--check) export GRID_INSTALL_ENV_READY=1; shift ;;
		--dry-run) DRY_RUN=TRUE; shift ;;
        --environment-file) recipe=${2:?--environment-file requires FILE}; shift 2 ;;
		*) printf 'ERROR: unknown argument %s; see --help\n' "$1" >&2; exit 2 ;;
	esac
done

unset R_HOME R_LIBS R_LIBS_USER R_LIBS_SITE R_ENVIRON_USER
ENV_NAME=${GRID_ENV_NAME:-${ENV_NAME:-grid}}
UPSTREAM=${PRSFORMER_UPSTREAM_DIR:-${GRID_PRSFORMER_UPSTREAM:-$DIR_SOFTWARE/PRSformer}}
REVISION=7dea4e1bb27975885c937f2be82f9243bc82bda7
[[ -x "$HOME/miniforge3/bin/conda" ]] && export PATH="$HOME/miniforge3/bin:$PATH"
[[ ! -x "$OS_GU_CONDA/bin/conda" ]] || export PATH="$OS_GU_CONDA/bin:$PATH"
if command -v mamba >/dev/null 2>&1; then solver=mamba
elif command -v conda >/dev/null 2>&1; then solver=conda
else
	printf 'ERROR: install Miniforge/Conda first.\n' >&2
	exit 2
fi

run() {
	printf 'RUN '; printf '%q ' "$@"; printf '\n'
	[[ $DRY_RUN == FALSE ]] || return 0
	"$@"
}

# Prepare the environment first: compute nodes may provide Git only through Conda.
if [[ -n ${GRID_CONDA_ENV:-} ]]; then
    env_target=(-p "$GRID_CONDA_ENV")
    if [[ -d "$GRID_CONDA_ENV/conda-meta" ]]; then env_action=update; else env_action=create; fi
else
    env_target=(-n "$ENV_NAME")
    if "$solver" env list | awk -v name="$ENV_NAME" '$1 == name {found=1} END {exit !found}'; then env_action=update; else env_action=create; fi
fi
[[ -s $recipe ]] || { echo "ERROR: Missing Conda recipe: $recipe" >&2; exit 2; }
# Online preparation may create the environment before an offline LSF worker.
# Reuse is explicit and still runs all dependency/device checks below.
if [[ ${GRID_INSTALL_ENV_READY:-0} == 1 ]]; then
    [[ -n ${GRID_CONDA_ENV:-} && -d $GRID_CONDA_ENV/conda-meta ]] || {
        echo 'ERROR: GRID_INSTALL_ENV_READY requires an existing GRID_CONDA_ENV' >&2; exit 2;
    }
else
    run "$solver" env "$env_action" "${env_target[@]}" -f "$recipe"
fi

# Check the pinned source through the selected environment. Never reset a checkout.
if [[ -e $UPSTREAM && $DRY_RUN == FALSE ]]; then
    [[ $("$solver" run "${env_target[@]}" git -C "$UPSTREAM" rev-parse HEAD) == "$REVISION" ]] || {
        printf 'ERROR: existing PRSformer source must be at %s.\n' "$REVISION" >&2; exit 2;
    }
    [[ -z $("$solver" run "${env_target[@]}" git -C "$UPSTREAM" status --porcelain) ]] || {
        printf 'ERROR: PRSformer checkout has local changes; inspect them first.\n' >&2; exit 2;
    }
fi


if [[ ! -e $UPSTREAM ]]; then
	run mkdir -p -- "$(dirname -- "$UPSTREAM")"
	run "$solver" run "${env_target[@]}" git clone https://github.com/23andMe/PRSformer.git "$UPSTREAM"
	run "$solver" run "${env_target[@]}" git -C "$UPSTREAM" checkout --detach "$REVISION"
fi

run "$solver" run "${env_target[@]}" python -m pip check
run "$solver" run "${env_target[@]}" python -c 'import yaml; print("YAML configuration runtime checked")'
run "$solver" run "${env_target[@]}" python -c 'import numpy,pandas,scipy,sklearn,matplotlib,openpyxl,rdata,joblib,torch,natten,pgenlib,pyreadr'
run "$solver" run "${env_target[@]}" Rscript --vanilla -e 'p <- c("shiny","DT","DBI","RSQLite","ggplot2","jsonlite","testthat"); stopifnot(all(vapply(p,requireNamespace,logical(1),quietly=TRUE)))'
run "$solver" run "${env_target[@]}" plink2 --version
check_device=${WORKFLOW_DEVICE:-auto}
if [[ $check_device == auto && $DRY_RUN == FALSE ]]; then
    if "$solver" run "${env_target[@]}" python -c 'import torch; raise SystemExit(0 if torch.cuda.is_available() else 1)'; then
        check_device=cuda
    else check_device=cpu; fi
fi
if [[ $check_device == cuda ]]; then
    run "$solver" run "${env_target[@]}" python "$ROOT/f/main.py" --check-device
    # GPU checks execute on an allocated GPU node, including a forward/backward.
    run "$solver" run "${env_target[@]}" bash "$ROOT/f/prsformer.sh" \
        --python python --upstream-dir "$UPSTREAM" --device cuda --check-runtime
else
    run "$solver" run "${env_target[@]}" python -c 'import torch; x=torch.eye(4); assert torch.equal(x@x,x); print("CPU runtime checked; CUDA/NATTEN device execution has not been validated")'
fi

if [[ $DRY_RUN == TRUE ]]; then
	printf 'PLAN ONLY: no environment, source checkout or data was changed.\n'
else
	printf 'Shared GRID / PRSformer environment %s runtime checks completed.\n' "${GRID_CONDA_ENV:-$ENV_NAME}"
fi
