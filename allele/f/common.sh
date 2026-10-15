#!/usr/bin/env bash
set -euo pipefail

allele_load_config() {
  local root=${1:?root required}
  ALLELE_ROOT=$root
  local config=${ALLELE_CONFIG:-$root/allele.env}
  if [[ -f "$config" ]]; then
    # shellcheck disable=SC1090,SC1091
    source "$config"
  fi
  export ALLELE_ENV_NAME=${ALLELE_ENV_NAME:-allele}
  export ALLELE_OUT=${ALLELE_OUT:-/mnt/d/analysis/allele}
  # Keep an explicitly configured legacy ALLELE_GNOMAD input working.
  export ALLELE_VCF_DIR=${ALLELE_VCF_DIR:-${ALLELE_GNOMAD:-/mnt/f/annot/gnomad}}
  export ALLELE_GNOMAD=${ALLELE_GNOMAD:-/mnt/f/annot/gnomad}
  export ALLELE_FASTA=${ALLELE_FASTA:-/mnt/f/gen/1hgp/38/GRCH38.fasta.gz}
  export ALLELE_ATLAS_SCORERS=${ALLELE_ATLAS_SCORERS:-ATAC,DNASE,CHIP_TF,CHIP_HISTONE,CAGE,PROCAP}
  export ALLELE_ATLAS_CHUNK=${ALLELE_ATLAS_CHUNK:-128}
  export ALLELE_BIALLELIC_RATIO=${ALLELE_BIALLELIC_RATIO:-2}
  export ALLELE_SEED=${ALLELE_SEED:-20261015}
  export ALLELE_API_KEY_FILE=${ALLELE_API_KEY_FILE:-/mnt/d/files/alphaG.api}
  if [[ -z ${ALPHA_GENOME_API_KEY:-} && -r "$ALLELE_API_KEY_FILE" ]]; then
    # Read one plain-text line; never source/execute a credential file.
    IFS= read -r ALPHA_GENOME_API_KEY < "$ALLELE_API_KEY_FILE" || true
    ALPHA_GENOME_API_KEY=${ALPHA_GENOME_API_KEY%$'\r'}
    export ALPHA_GENOME_API_KEY
  fi
}

allele_die() { echo "ERROR: $*" >&2; exit 2; }
allele_note() { printf '[allele] %s\n' "$*" >&2; }

allele_require_file() { [[ -s "$1" ]] || allele_die "missing file: $1"; }
allele_require_cmd() { command -v "$1" >/dev/null 2>&1 || allele_die "missing command: $1"; }

allele_find_conda() {
  local candidate
  for candidate in "${ALLELE_CONDA_EXE:-}" "${CONDA_EXE:-}" \
      "$(command -v conda || true)" "$HOME/anaconda3/bin/conda" \
      "$HOME/miniforge3/bin/conda" "$HOME/miniconda3/bin/conda"; do
    if [[ -n "$candidate" && -x "$candidate" ]]; then
      ALLELE_CONDA=$candidate
      return
    fi
  done
  allele_die 'conda not found; install Miniforge/Miniconda or set ALLELE_CONDA_EXE'
}

allele_env_exists() {
  "$ALLELE_CONDA" env list | awk -v name="$ALLELE_ENV_NAME" \
    '$1 == name { found=1 } END { exit !found }'
}

allele_require_env() {
  allele_find_conda
  allele_env_exists || allele_die "conda environment '$ALLELE_ENV_NAME' is missing. Run: $ALLELE_ROOT/install.sh"
}

allele_python() {
  # Do not exec here: `all` must continue after each Python stage returns.
  "$ALLELE_CONDA" run --no-capture-output -n "$ALLELE_ENV_NAME" python "$@"
}
