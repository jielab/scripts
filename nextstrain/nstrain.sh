#!/usr/bin/env bash
set -euo pipefail




#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 🚩 Shell options, paths, and shared inputs
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 🚩 hChikV GISAID -> Nextstrain/Augur pipeline
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 🚩 0. Paths and parameters
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
dir0=${dir0:-/mnt/d}
dirdat=${dirdat:-$dir0/data/viome/raw}
dirout=${dirout:-$dir0/analysis/viome}
CONDA_ENV=${CONDA_ENV:-nextstrain}
MIN_LENGTH=${MIN_LENGTH:-11000}
VIEW=${VIEW:-0}

raw_gisaid_dir=${raw_gisaid_dir:-$dirdat/gisaid}
raw_fasta=${raw_fasta:-$raw_gisaid_dir/gisaid.fasta}
patient_status=${patient_status:-$raw_gisaid_dir/patient_status.tsv}
seq_tech=${seq_tech:-$raw_gisaid_dir/seq_tech.tsv}
centroids=${centroids:-$dirdat/countries.centroids.tsv}
ref_fasta=${ref_fasta:-$dirdat/hchikv.ref.fasta}
ref_gb=${ref_gb:-$dirdat/hchikv.ref.gb}

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
py=${py:-python3}

mkdir -p "$dirout" "$dirout/auspice"
cd "$dirout"

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


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 🚩 1. Activate Nextstrain environment
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
if [[ -n "$CONDA_ENV" ]]; then
  if command -v conda >/dev/null 2>&1; then
    # shellcheck disable=SC1091
    source "$(conda info --base)/etc/profile.d/conda.sh"
    conda activate "$CONDA_ENV"
  else
    echo "WARNING: conda not found; continuing with current environment." >&2
  fi
fi


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 🚩 2. Prepare FASTA and metadata
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#    - FASTA header is rewritten to EPI_ISL_xxx
#    - metadata.tsv uses strain = EPI_ISL_xxx
#    - country coordinates are NOT hard-coded here
$py "$script_dir/f/nstrain.py" prepare \
  --fasta "$raw_fasta" \
  --patient-status "$patient_status" \
  --seq-tech "$seq_tech" \
  --outdir "$dirout"


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 🚩 3. Build sequence index
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
augur index \
  --sequences sequences.fasta \
  --output sequence_index.tsv


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 🚩 4. Filter low-quality / incomplete records
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#    unknown.exclude removes records with unknown date, region, or country.
augur filter \
  --sequences sequences.fasta \
  --metadata metadata.tsv \
  --exclude unknown.exclude \
  --min-length "$MIN_LENGTH" \
  --output-sequences filtered.fasta \
  --output-metadata filtered_metadata.tsv


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 🚩 5. Align to reference genome
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
augur align \
  --sequences filtered.fasta \
  --reference-sequence "$ref_fasta" \
  --output alignment.fasta \
  --fill-gaps


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 🚩 6. Build phylogenetic tree
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
augur tree \
  --alignment alignment.fasta \
  --tree-builder-args '-ninit 10 -n 4' \
  --output tree_raw.nwk


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 🚩 7. Time-resolved tree refinement
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
augur refine \
  --tree tree_raw.nwk \
  --alignment alignment.fasta \
  --metadata filtered_metadata.tsv \
  --timetree \
  --output-tree tree.nwk \
  --output-node-data branch_lengths.json


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 🚩 8. Add country coordinates and write lat_longs.tsv
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#    This reads your existing countries.centroids.tsv.
$py "$script_dir/f/nstrain.py" add-coords \
  --metadata filtered_metadata.tsv \
  --centroids "$centroids" \
  --output filtered_metadata_with_coords.tsv \
  --lat-longs lat_longs.tsv \
  --missing-output missing_country_coords.tsv


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 🚩 9. Infer geographic traits when enough states are available
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
node_data=(branch_lengths.json)

n_country=$($py "$script_dir/f/nstrain.py" count-values --metadata filtered_metadata_with_coords.tsv --column country --ignore-unknown)
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

n_division=$($py "$script_dir/f/nstrain.py" count-values --metadata filtered_metadata_with_coords.tsv --column division --ignore-unknown)
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


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 🚩 10. Reconstruct ancestral nucleotide mutations
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
augur ancestral \
  --tree tree.nwk \
  --alignment alignment.fasta \
  --output-node-data nt_muts.json \
  --inference joint \
  --root-sequence "$ref_fasta"
node_data+=(nt_muts.json)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 🚩 11. Translate amino-acid mutations
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
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


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 🚩 12. Write Auspice config once
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
$py "$script_dir/f/nstrain.py" write-config \
  --output auspice_config.json \
  --title "hChikV genomic analysis" \
  --maintainer "Jie Huang"


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 🚩 13. Export Auspice JSON
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
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


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 🚩 14. View locally, optional
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
echo "Done: $dirout/auspice/chikv_with_map.json"
echo "To view: cd $dirout && auspice view --datasetDir auspice"

if [[ "$VIEW" == "1" ]]; then
  auspice view --datasetDir auspice
fi
