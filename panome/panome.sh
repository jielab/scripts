#!/usr/bin/env bash
set -euo pipefail
panome_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"

# Default analysis settings. Explicit command-line options take precedence.
panome_defaults=(--device cuda --cores 16 --quality-teacher all
                 --run-name v5_attention_allteachers)

if [[ -n "${PANOME_PYTHON:-}" ]]; then
  panome_python="$PANOME_PYTHON"
elif [[ -x "$HOME/venvs/panome-v3/bin/python" ]]; then
  # Reuse the existing external environment; never create a local .venv.
  panome_python="$HOME/venvs/panome-v3/bin/python"
else
  panome_python="python3"
fi

# Each outcome/layer pair gets its own Python process, defaults and output path.
panome_biom=prot
panome_y=cvd_cad
panome_args=()
panome_shared_path=false
panome_omics_path=false
panome_diagnosis_col=false
while [[ $# -gt 0 ]]; do
  case "$1" in
    --Y)
      panome_y=${2:?missing value for --Y}; shift 2 ;;
    --Y=*)
      panome_y=${1#*=}; shift ;;
    --biom)
      panome_biom=${2:?missing value for --biom}; shift 2 ;;
    --biom=*)
      panome_biom=${1#*=}; shift ;;
    --run-dir|--run-dir=*|--output|--output=*)
      panome_shared_path=true
      panome_args+=("$1"); shift ;;
    --omics-file|--omics-file=*)
      panome_omics_path=true
      panome_args+=("$1"); shift ;;
    --diagnosis-col|--diagnosis-col=*)
      panome_diagnosis_col=true
      panome_args+=("$1"); shift ;;
    *)
      panome_args+=("$1"); shift ;;
  esac
done

if [[ ! "$panome_y" =~ ^[^,[:space:]]+(,[^,[:space:]]+)*$ ]]; then
  echo "Error: --Y expects a nonempty outcome or comma-separated outcomes, e.g. cvd_cad,ra." >&2
  exit 2
fi
IFS=',' read -r -a panome_outcomes <<< "$panome_y"
declare -A panome_seen_outcomes=()
for panome_outcome in "${panome_outcomes[@]}"; do
  if [[ -n "${panome_seen_outcomes[$panome_outcome]:-}" ]]; then
    echo "Error: --Y must not repeat an outcome." >&2
    exit 2
  fi
  panome_seen_outcomes[$panome_outcome]=1
done

if [[ ! "$panome_biom" =~ ^(prot|met)(,(prot|met))*$ ]]; then
  echo "Error: --biom expects prot, met, or prot,met." >&2
  exit 2
fi
IFS=',' read -r -a panome_layers <<< "$panome_biom"
if [[ ${#panome_layers[@]} -gt 1 ]]; then
  if [[ ${#panome_layers[@]} -ne 2 || "${panome_layers[0]}" == "${panome_layers[1]}" ]]; then
    echo "Error: --biom must not repeat a molecular layer." >&2
    exit 2
  fi
  if [[ "$panome_omics_path" == true ]]; then
    echo "Error: --omics-file requires a single --biom." >&2
    exit 2
  fi
fi
if [[ ${#panome_outcomes[@]} -gt 1 && "$panome_diagnosis_col" == true ]]; then
  echo "Error: --diagnosis-col requires a single --Y; multiple outcomes use fod_icd10_<Y>." >&2
  exit 2
fi
if [[ ${#panome_outcomes[@]} -gt 1 || ${#panome_layers[@]} -gt 1 ]]; then
  if [[ "$panome_shared_path" == true ]]; then
    echo "Error: --run-dir and --output require a single --Y and --biom. Use --analysis-root and --run-name for separate runs." >&2
    exit 2
  fi
  for panome_outcome in "${panome_outcomes[@]}"; do
    for panome_layer in "${panome_layers[@]}"; do
      echo "[PANOME] START Y=$panome_outcome biom=$panome_layer" >&2
      "$panome_python" "$panome_dir/f/panome.py" "${panome_defaults[@]}" \
        "${panome_args[@]}" --Y "$panome_outcome" --biom "$panome_layer"
      echo "[PANOME] DONE Y=$panome_outcome biom=$panome_layer" >&2
    done
  done
else
  exec "$panome_python" "$panome_dir/f/panome.py" "${panome_defaults[@]}" \
    "${panome_args[@]}" --Y "${panome_outcomes[0]}" --biom "${panome_layers[0]}"
fi
