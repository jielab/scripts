#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<'HELP'
Usage: ./nstrain.sh [options]

  --data-root DIR       Data root [/mnt/d]
  --input-dir DIR       Virome input root [<data-root>/data/viome/raw]
  --output-dir DIR      Nextstrain output [<data-root>/analysis/viome]
  --conda-env NAME      Nextstrain Conda environment [nextstrain]
  --min-length INT      Minimum sequence length [11000]
  --view TRUE|FALSE     Start Auspice after export [FALSE]
  --gisaid-dir DIR      GISAID source directory [<input-dir>/gisaid]
  --fasta FILE          GISAID FASTA
  --patient-status FILE Patient metadata TSV
  --seq-tech FILE       Sequencing-technology TSV
  --centroids FILE      Country centroids TSV
  --ref-fasta FILE      Reference FASTA
  --ref-genbank FILE    Reference GenBank file
  --python FILE         Python interpreter [python3]
  -h, --help            Show this message
HELP
}

die() { echo "ERROR: $*" >&2; exit 2; }
bool01() {
  case "${1,,}" in
    true|1|yes|y) printf '1\n' ;;
    false|0|no|n) printf '0\n' ;;
    *) return 2 ;;
  esac
}

data_root=/mnt/d
input_dir=""
output_dir=""
conda_env=nextstrain
min_length=11000
view=0
gisaid_dir=""
raw_fasta=""
patient_status=""
seq_tech=""
centroids=""
ref_fasta=""
ref_gb=""
python_bin=python3

while (( $# )); do
  case "$1" in
    --data-root) data_root=${2:?}; shift 2 ;;
    --input-dir) input_dir=${2:?}; shift 2 ;;
    --output-dir) output_dir=${2:?}; shift 2 ;;
    --conda-env) conda_env=${2:?}; shift 2 ;;
    --min-length) min_length=${2:?}; shift 2 ;;
    --view) view=$(bool01 "${2:?}") || die "--view requires TRUE or FALSE"; shift 2 ;;
    --gisaid-dir) gisaid_dir=${2:?}; shift 2 ;;
    --fasta) raw_fasta=${2:?}; shift 2 ;;
    --patient-status) patient_status=${2:?}; shift 2 ;;
    --seq-tech) seq_tech=${2:?}; shift 2 ;;
    --centroids) centroids=${2:?}; shift 2 ;;
    --ref-fasta) ref_fasta=${2:?}; shift 2 ;;
    --ref-genbank) ref_gb=${2:?}; shift 2 ;;
    --python) python_bin=${2:?}; shift 2 ;;
    -h|--help) usage; exit 0 ;;
    *) die "unknown option: $1" ;;
  esac
done

[[ $min_length =~ ^[1-9][0-9]*$ ]] || die "--min-length must be a positive integer"
[[ -n $input_dir ]] || input_dir=$data_root/data/viome/raw
[[ -n $output_dir ]] || output_dir=$data_root/analysis/viome
[[ -n $gisaid_dir ]] || gisaid_dir=$input_dir/gisaid
[[ -n $raw_fasta ]] || raw_fasta=$gisaid_dir/gisaid.fasta
[[ -n $patient_status ]] || patient_status=$gisaid_dir/patient_status.tsv
[[ -n $seq_tech ]] || seq_tech=$gisaid_dir/seq_tech.tsv
[[ -n $centroids ]] || centroids=$input_dir/countries.centroids.tsv
[[ -n $ref_fasta ]] || ref_fasta=$input_dir/hchikv.ref.fasta
[[ -n $ref_gb ]] || ref_gb=$input_dir/hchikv.ref.gb


# 🚩 hChikV inputs and output layout
script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

mkdir -p "$output_dir" "$output_dir/auspice"
cd "$output_dir"

require_file() {
  local f="$1"
  local label="$2"
  if [[ ! -s "$f" ]]; then
    echo "ERROR: missing $label: $f" >&2
    exit 1
  fi
}

require_file "$raw_fasta" "GISAID FASTA"
require_file "$patient_status" "GISAID patient_status.tsv"
require_file "$seq_tech" "GISAID seq_tech.tsv"
require_file "$centroids" "countries.centroids.tsv"
require_file "$ref_fasta" "hChikV reference FASTA"


# 🚩 1. Activate Nextstrain environment
if [[ -n "$conda_env" ]]; then
  if command -v conda >/dev/null 2>&1; then
    # shellcheck disable=SC1091
    source "$(conda info --base)/etc/profile.d/conda.sh"
    conda activate "$conda_env"
  else
    echo "WARNING: conda not found; continuing with current environment." >&2
  fi
fi


# 🚩 2. Prepare FASTA and metadata
$python_bin "$script_dir/f/nstrain.py" prepare \
  --fasta "$raw_fasta" \
  --patient-status "$patient_status" \
  --seq-tech "$seq_tech" \
  --outdir "$output_dir"


# 🚩 3. Build sequence index
augur index \
  --sequences sequences.fasta \
  --output sequence_index.tsv


# 🚩 4. Filter low-quality / incomplete records
#    unknown.exclude removes records with unknown date, region, or country.
augur filter \
  --sequences sequences.fasta \
  --metadata metadata.tsv \
  --exclude unknown.exclude \
  --min-length "$min_length" \
  --output-sequences filtered.fasta \
  --output-metadata filtered_metadata.tsv


# 🚩 5. Align to reference genome
augur align \
  --sequences filtered.fasta \
  --reference-sequence "$ref_fasta" \
  --output alignment.fasta \
  --fill-gaps


# 🚩 6. Build phylogenetic tree
augur tree \
  --alignment alignment.fasta \
  --tree-builder-args '-ninit 10 -n 4' \
  --output tree_raw.nwk


# 🚩 7. Time-resolved tree refinement
augur refine \
  --tree tree_raw.nwk \
  --alignment alignment.fasta \
  --metadata filtered_metadata.tsv \
  --timetree \
  --output-tree tree.nwk \
  --output-node-data branch_lengths.json


# 🚩 8. Add country coordinates and write lat_longs.tsv
#    This reads your existing countries.centroids.tsv.
$python_bin "$script_dir/f/nstrain.py" add-coords \
  --metadata filtered_metadata.tsv \
  --centroids "$centroids" \
  --output filtered_metadata_with_coords.tsv \
  --lat-longs lat_longs.tsv \
  --missing-output missing_country_coords.tsv


# 🚩 9. Infer geographic traits when enough states are available
node_data=(branch_lengths.json)

n_country=$($python_bin "$script_dir/f/nstrain.py" count-values --metadata filtered_metadata_with_coords.tsv --column country --ignore-unknown)
if (( n_country > 1 )); then
  augur traits \
    --tree tree.nwk \
    --metadata filtered_metadata_with_coords.tsv \
    --columns country \
    --output-node-data traits_country.json
  node_data+=(traits_country.json)
else
  echo "Skipping country traits: fewer than two country states."
fi

n_division=$($python_bin "$script_dir/f/nstrain.py" count-values --metadata filtered_metadata_with_coords.tsv --column division --ignore-unknown)
if (( n_division > 1 )); then
  augur traits \
    --tree tree.nwk \
    --metadata filtered_metadata_with_coords.tsv \
    --columns division \
    --output-node-data traits_division.json
  node_data+=(traits_division.json)
else
  echo "Skipping division traits: fewer than two division states."
fi


# 🚩 10. Reconstruct ancestral nucleotide mutations
augur ancestral \
  --tree tree.nwk \
  --alignment alignment.fasta \
  --output-node-data nt_muts.json \
  --inference joint \
  --root-sequence "$ref_fasta"
node_data+=(nt_muts.json)


# 🚩 11. Translate amino-acid mutations
#     Use GenBank reference if available because it contains gene annotations.
if [[ -s "$ref_gb" ]]; then
  augur translate \
    --tree tree.nwk \
    --ancestral-sequences nt_muts.json \
    --reference-sequence "$ref_gb" \
    --output-node-data aa_muts.json
  node_data+=(aa_muts.json)
else
  echo "WARNING: $ref_gb not found; skipping augur translate / aa_muts.json." >&2
fi


# 🚩 12. Write Auspice config once
$python_bin "$script_dir/f/nstrain.py" write-config \
  --output auspice_config.json \
  --title "hChikV genomic analysis" \
  --maintainer "Jie Huang"


# 🚩 13. Export Auspice JSON
#     Map resolution is country because countries.centroids.tsv provides country coordinates.
export_args=(
  augur export v2
  --tree tree.nwk
  --metadata filtered_metadata_with_coords.tsv
  --metadata-id-columns strain
  --node-data "${node_data[@]}"
  --auspice-config auspice_config.json
  --output auspice/chikv_with_map.json
  --validation-mode warn
)

if [[ -s lat_longs.tsv ]]; then
  export_args+=(--geo-resolutions country --lat-longs lat_longs.tsv)
fi

if [[ -s colors.tsv ]]; then
  export_args+=(--colors colors.tsv)
fi

"${export_args[@]}"


# 🚩 14. View locally, optional
echo "Done: $output_dir/auspice/chikv_with_map.json"
echo "To view: cd $output_dir && auspice view --datasetDir auspice"

if [[ "$view" == "1" ]]; then
  auspice view --datasetDir auspice
fi
