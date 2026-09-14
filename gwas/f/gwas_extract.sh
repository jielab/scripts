#!/usr/bin/env bash
set -euo pipefail
here=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
if [[ ${1:-} == --intersect ]]; then
    requested=$2 shared=$3 output=$4
else
    counts=$1 output=$2 maf=$3
fi
tmp=$(mktemp "${output}.tmp.XXXXXX")
trap 'rm -f -- "$tmp" "$tmp.counts.tsv"' EXIT
if [[ ${1:-} == --intersect ]]; then
    awk 'FILENAME == ARGV[1] { if (NF) wanted[$1]=1; next } $1 in wanted { print $1 }' "$requested" "$shared" > "$tmp"
else
    LC_ALL=C awk -v maf="$maf" -v summary="$tmp.counts.tsv" -f "$here/gwas_extract.awk" "$counts" > "$tmp"
fi
if [[ ! -s "$tmp" ]]; then
    echo "ERROR: No variants pass for $output; stopping before association" >&2
    exit 2
fi
[[ ! -f "$tmp.counts.tsv" ]] || mv -f -- "$tmp.counts.tsv" "$output.counts.tsv"
mv -f -- "$tmp" "$output"
