#!/usr/bin/env bash
set -euo pipefail
case ${1:-help} in
    -h|--help|help) echo 'Usage: ./0.mask.sh run   Install the bundled UCSC hg38 genomicSuperDups BED under DIR_FILES/allele3.'; exit 0 ;;
    run) [[ $# == 1 ]] || exit 2 ;;
    *) echo 'Usage: ./0.mask.sh run' >&2; exit 2 ;;
esac
root=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
source "${ALLELE_OS_CONFIG:-$root/../0f/0os.sh}"
dest=${ALLELE_REPEAT_MASK_BED:-$DIR_FILES/allele3/hg38.segmental_duplications.bed.gz}
source_bed="$root/config/hg38.segmental_duplications.bed.gz"
# Source: https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/genomicSuperDups.txt.gz
# Gzip SHA-256: 3f4388ed3eb7462034e89c1f45314160bdb14558502092326063105c095061df
# Columns 2,3,4; canonical chr1-22,X,Y.
# Source and BED are zero-based half-open; exact duplicate triples are removed.
checksum=5654a6543ddabe472253f530f0ce9d0efdd988680a4e6901fb56bd0f8dde8196
valid() { [[ -f $1 ]] && [[ $(sha256sum < "$1") == "$checksum  -" ]]; }
valid "$source_bed" || { echo "Invalid bundled mask: $source_bed" >&2; exit 1; }
if valid "$dest"; then echo "[mask ready] $dest"; exit 0; fi
mkdir -p -- "$(dirname -- "$dest")"
tmp=$(mktemp "$(dirname -- "$dest")/.segdup.XXXXXX")
trap 'rm -f -- "$tmp"' EXIT
cp -- "$source_bed" "$tmp"
gzip -t "$tmp"
valid "$tmp"
mv -fT -- "$tmp" "$dest"
echo "[mask ready] $dest"
