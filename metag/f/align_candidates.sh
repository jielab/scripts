#!/usr/bin/env bash

# Public directory options use --dir-XXX; legacy spellings remain compatible.
if [[ "${BASH_SOURCE[0]}" == "$0" ]]; then
    source /mnt/d/scripts/0f/cli_dirs.sh
    script_directory_args "$(realpath -- "${BASH_SOURCE[0]}")" "$@"
    set -- "${SCRIPT_DIRECTORY_ARGS[@]}"
    unset SCRIPT_DIRECTORY_ARGS
fi
# Optional independent ALIGNMENT EVIDENCE, NOT automatic species confirmation.
# The index should contain the candidate AND relevant close relatives/decoys.
set -euo pipefail
INDEX=""; R1=""; R2=""; OUT=""; THREADS=12
usage() {
	cat <<'HELP'
Modules:
    Default  Independent candidate alignment evidence

Main options:
    --index PREFIX --r1 FILE --r2 FILE  Index and reads; R2 optional
    --dir-out NEW_DIR --threads N  New output directory and threads [12]
    -h, --help  Show help

Example:
    # Paired-end form; replace the index/read paths and use a new output directory.
    bash metag.sh align --index /path/to/index --r1 /path/to/reads_R1.fastq.gz \
      --r2 /path/to/reads_R2.fastq.gz --out /path/to/new_output --threads 12

    ./metag.sh align --index /path/to/index --r1 /path/to/reads.fastq.gz --dir-out /path/to/new_output
HELP
}
while (($#)); do
  case "$1" in
    --index) INDEX="$2"; shift 2;;
    --r1) R1="$2"; shift 2;;
    --r2) R2="$2"; shift 2;;
    --out) OUT="$2"; shift 2;;
    --threads) THREADS="$2"; shift 2;;
    -h|--help) usage; exit 0;;
    *) usage >&2; exit 2;;
  esac
done
[[ -n "$INDEX" && -f "$R1" && -n "$OUT" ]] || exit 2
[[ -f "${INDEX}.1.bt2" || -f "${INDEX}.1.bt2l" ]] || { echo 'Missing Bowtie2 index' >&2; exit 1; }
[[ ! -e "$OUT" ]] || { echo 'Output directory already exists; use a new one.' >&2; exit 1; }
mkdir -p "$OUT/tmp"
CMD=(bowtie2 --very-sensitive -p "$THREADS" -x "$INDEX")
if [[ -n "$R2" ]]; then [[ -f "$R2" ]] || exit 2; CMD+=(-1 "$R1" -2 "$R2"); else CMD+=(-U "$R1"); fi
printf '%q ' "${CMD[@]}" > "$OUT/command.txt"; printf '\n' >> "$OUT/command.txt"
"${CMD[@]}" 2> "$OUT/bowtie2.log" | samtools sort -@ 2 -m 512M -T "$OUT/tmp/sort" -o "$OUT/all.bam" -
samtools index "$OUT/all.bam"
# 2308 = unmapped + secondary + supplementary; do not count every alignment as a molecule.
samtools view -b -q 30 -F 2308 -o "$OUT/primary.q30.bam" "$OUT/all.bam"
samtools index "$OUT/primary.q30.bam"
samtools coverage "$OUT/primary.q30.bam" > "$OUT/coverage.by_contig.tsv"
samtools flagstat "$OUT/all.bam" > "$OUT/flagstat.txt"
cat > "$OUT/INTERPRETATION.txt" <<'EOF'
This is alignment evidence only. MAPQ>=30 is meaningful relative to the references
present in this index; an isolated human/viral genome cannot exclude close relatives.
Inspect coverage breadth, multiple independent loci, mismatches, duplication,
negative controls and the competing-reference panel. This script does not infer
infection, living organisms, person identity, or confirmed species presence.
Output rows are contigs, NOT automatically aggregated species.
EOF
