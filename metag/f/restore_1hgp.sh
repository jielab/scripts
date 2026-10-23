#!/usr/bin/env bash
# Restore the two reference distributions used by allele3/gwas/gu/ukb.
# Existing build directories are never overwritten. No uncompressed whole
# genome is written: decompression is piped directly into BGZF compression.
set -euo pipefail

usage() {
    cat <<'EOF'
Usage: bash restore_1hgp.sh [--root DIR] [--build 37|38|both] [--threads N] [--dry-run]

Default destination: /mnt/f/annot/1hgp
Default builds: both; compression threads: 4
Requires curl, gzip, bgzip, samtools, md5sum, sha256sum, python3.
Publishes BUILD/GRCHBUILD.fasta.gz plus freshly built .fai and .gzi.
Existing BUILD directories cause an error, even if empty. Use a new root.
Failed staging directories are preserved and their location is printed.
EOF
}

root=/mnt/f/annot/1hgp
build_choice=both
threads=4
dry_run=0
while (($#)); do
    case "$1" in
        --root|--build|--threads)
            (($# >= 2)) || { usage >&2; exit 2; }
            case "$1" in
                --root) root=$2 ;;
                --build) build_choice=$2 ;;
                --threads) threads=$2 ;;
            esac
            shift 2 ;;
        --dry-run) dry_run=1; shift ;;
        -h|--help) usage; exit 0 ;;
        *) usage >&2; exit 2 ;;
    esac
done
[[ -n "$root" ]] || { echo 'Empty root is not allowed' >&2; exit 2; }
[[ "$threads" =~ ^[1-9][0-9]*$ ]] || { echo 'Invalid threads' >&2; exit 2; }
case "$build_choice" in
    both) builds=(37 38) ;;
    37|38) builds=("$build_choice") ;;
    *) echo 'Build must be 37, 38 or both' >&2; exit 2 ;;
esac

source_for() {
    if [[ "$1" == 37 ]]; then
        url=https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/human_g1k_v37.fasta.gz
        # Pinned original-download digest in the existing 2026-10-03 cleanup.tsv.
        expected_md5=45f81df94f0408d082363e34a081ed81
        expected_sha=2f9cd9e853a9284c53884e6a551b1c7284795dd053f255d630aeeb114d1fa81f
        count=84
        length=3101804739
        probe=1:10001-10040
        checksum_source='/mnt/f/gen/1hgp/cleanup.tsv (local historical receipt)'
    else
        base=https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids
        url=$base/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz
        # Verified against the NCBI md5checksums.txt on 2026-10-23.
        expected_md5=a08035b6a6e31780e96a34008ff21bd6
        expected_sha=9cce8b926416dd96b152deea85188495b75f7ac8d634cc723a017067be8702b7
        count=195
        length=3099922541
        probe=chr1:10001-10040
        checksum_source=$base/md5checksums.txt
    fi
}

for build in "${builds[@]}"; do
    source_for "$build"
    printf 'GRCh%s\nSource: %s\nDestination: %s/%s\nMD5 (original download): %s\n' \
        "$build" "$url" "$root" "$build" "$expected_md5"
done
((dry_run == 0)) || exit 0
for tool in curl gzip bgzip samtools md5sum sha256sum python3; do
    command -v "$tool" >/dev/null || { echo "Missing tool: $tool" >&2; exit 1; }
done
for build in "${builds[@]}"; do
    [[ ! -e "$root/$build" && ! -L "$root/$build" ]] || {
        echo "Destination exists; nothing downloaded: $root/$build" >&2
        exit 1
    }
done
mkdir -p -- "$root"
stage=''
trap 'rc=$?; if ((rc != 0)) && [[ -n "$stage" ]]; then printf "Failed staging directory retained: %s\n" "$stage" >&2; fi' EXIT
for build in "${builds[@]}"; do
    source_for "$build"
    stage=$(mktemp -d "$root/.restore-$build.XXXXXXXX")
    printf 'Staging: %s\n' "$stage"
    mkdir -- "$stage/ready"
    curl --fail --location --retry 5 --retry-all-errors --connect-timeout 30 \
        --output "$stage/original.fa.gz" "$url"
    printf '%s  %s\n' "$expected_md5" "$stage/original.fa.gz" | md5sum --check -
    # The sequence-text SHA256 also detects any corruption during transcoding.
    fa="$stage/ready/GRCH$build.fasta.gz"
    gzip -dc -- "$stage/original.fa.gz" | bgzip -@ "$threads" -c > "$fa"
    actual_sha=$(gzip -dc -- "$fa" | sha256sum)
    [[ "${actual_sha%% *}" == "$expected_sha" ]] || {
        echo 'Uncompressed FASTA SHA256 mismatch' >&2; exit 1;
    }
    samtools faidx "$fa"
    [[ -s "$fa.fai" && -s "$fa.gzi" ]]
    python3 - "$fa.fai" "$count" "$length" <<'PY'
import sys
from pathlib import Path
rows = [line.split('\t') for line in Path(sys.argv[1]).read_text().splitlines()]
if len(rows) != int(sys.argv[2]) or sum(int(r[1]) for r in rows) != int(sys.argv[3]):
    raise SystemExit('Reference sequence count or total length mismatch')
PY
    samtools faidx "$fa" "$probe" > "$stage/ready/index-check.fa"
    printf 'source_url\t%s\nsource_md5\t%s\nsource_md5_provenance\t%s\nuncompressed_sha256\t%s\n' \
        "$url" "$expected_md5" "$checksum_source" "$expected_sha" > "$stage/ready/source.tsv"
    (cd "$stage/ready" && sha256sum "GRCH$build.fasta.gz" "GRCH$build.fasta.gz.fai" \
        "GRCH$build.fasta.gz.gzi" > SHA256SUMS)
    # -T avoids nesting under an unexpected destination; -n preserves existing data.
    mv -T -n -- "$stage/ready" "$root/$build"
    [[ ! -d "$stage/ready" ]] || { echo 'Destination appeared; not overwritten' >&2; exit 1; }
    # Remove only this script's verified redundant download and its empty stage.
    rm -- "$stage/original.fa.gz"
    rmdir -- "$stage"
    stage=''
    printf 'Ready: %s/%s/GRCH%s.fasta.gz\n' "$root" "$build" "$build"
done
