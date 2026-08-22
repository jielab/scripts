#!/usr/bin/env bash
set -euo pipefail
ROOT=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd); F=$ROOT/f
ACTION=${1:-build}; shift || true
BACKEND=${GU_ARG_BACKEND:-tsinfer}
IN_DIR=${GU_TARGET_VCF_DIR:?set GU_TARGET_VCF_DIR}
PREP_DIR=${GU_ARG_VCF_DIR:-${GU_TARGET_ROOT:-.}/vcf.4arg}
OUT_DIR=${GU_ARG_DIR:-${GU_TARGET_ROOT:-.}/arg}
TMP=${GU_ARG_TMP_ROOT:-${XDG_CACHE_HOME:-$HOME/.cache}/gu/arg}
THREADS=${GU_ARG_THREADS:-8}
SAMPLE_MATCH_MIN_WORK=${GU_ARG_SAMPLE_MATCH_MIN_WORK:-50000000}
CHRS=${GU_CHRS:-"1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X"}
MISS=${GU_ARG_MISSING_MAX:-0.05}; MAC=${GU_ARG_MAC_MIN:-2}
mkdir -p "$PREP_DIR" "$OUT_DIR" "$TMP"

usage(){ cat <<EOF2
Usage: ./gu.sh arg [prepare|build|check]
GU_ARG_BACKEND=tsinfer|external

Recommended:
  1KG exploratory/engineering: GU_ARG_BACKEND=tsinfer
  TRACE sensitivity validation: build posterior ARGs with SINGER, convert to .trees/.tsz,
                                then use GU_ARG_BACKEND=external and GU_ARG_DIR=...

GU v2 does NOT automatically split chr1/chrX into coordinate chunks for TRACE.
Coordinate chunks are different genomic regions, not posterior ARG replicates.
The default chromosome set is 1..22 and X. An explicit GU_CHRS value is used unchanged.
EOF2
}
[[ $ACTION == help || $ACTION == -h || $ACTION == --help ]] && { usage; exit 0; }

prep_one(){
  local c=$1 in=$IN_DIR/chr${c}.vcf.gz out=$PREP_DIR/chr${c}.vcf.gz
  [[ -s $in ]] || { echo "ERROR missing $in" >&2; return 1; }
  [[ -s $out && -s $out.tbi ]] && return 0
  bash "$F/arg_vcf_prep.sh" --in "$in" --out "$out" --threads "$THREADS" --missing-max "$MISS" --mac-min "$MAC" --tmp-root "$TMP"
}

build_one(){
  local c=$1 in=$PREP_DIR/chr${c}.vcf.gz out=$OUT_DIR/chr${c}.trees
  [[ -s $out ]] && return 0
  local extra=()
  [[ ${GU_ARG_ONE_BIT:-0} == 1 ]] && extra+=(--one-bit)
  python3 "$F/arg_tsinfer.py" --vcf "$in" --out "$out" --sample-map "$OUT_DIR/chr${c}.sample_map.tsv" \
    --qc-json "$OUT_DIR/chr${c}.arg_qc.json" --threads "$THREADS" --tmp-root "$TMP" \
    --sample-match-min-work "$SAMPLE_MATCH_MIN_WORK" \
    --recombination-rate "${GU_RECOMB_RATE:-1e-8}" --mutation-rate "${GU_MUT_RATE:-1.25e-8}" \
    "${extra[@]}"
}

check_external(){
  local c n bad=0
  for c in $CHRS; do
    n=$(find "$OUT_DIR" -type f \( -name '*.trees' -o -name '*.tsz' \) -print | awk -v c="$c" 'BEGIN{IGNORECASE=1}{b=$0;gsub(/.*\//,"",b);if(b ~ ("(^|[^A-Za-z0-9])chr" c "([^A-Za-z0-9]|$)")) n++}END{print n+0}')
    ((n>0)) || { echo "ERROR: external ARG missing chr$c under $OUT_DIR" >&2; bad=1; }
  done
  ((bad==0))
}

case "$BACKEND:$ACTION" in
  tsinfer:prepare) for c in $CHRS; do prep_one "$c"; done ;;
  tsinfer:build) for c in $CHRS; do prep_one "$c"; build_one "$c"; done ;;
  tsinfer:check) for c in $CHRS; do [[ -s $OUT_DIR/chr${c}.trees ]] || { echo "missing chr$c" >&2; exit 1; }; done ;;
  external:check|external:build|external:prepare) check_external ;;
  *) echo "ERROR unsupported GU_ARG_BACKEND/ACTION: $BACKEND/$ACTION" >&2; usage; exit 2 ;;
esac
