#!/usr/bin/env bash
# Transfer exact-match ancestral alleles from 1KG.
set -euo pipefail

usage(){ cat <<'USAGE'
Usage:
  add_aa_from_1kg.sh --in target.vcf.gz --ref 1kg.chrN.vcf.gz --out out.vcf.gz [--threads N]

The reference VCF must contain INFO/AA. INFO/AA values such as A||| are normalized to
one base. A target site is retained only when CHROM, POS, REF and ALT match the 1KG
record exactly and AA is one of REF/ALT. No strand complement, REF/ALT guessing, or
liftover is performed.
USAGE
}

IN="" REF="" OUT="" THREADS=4
while (($#)); do
  case "$1" in
    --in) IN=$2; shift 2 ;;
    --ref) REF=$2; shift 2 ;;
    --out) OUT=$2; shift 2 ;;
    --threads) THREADS=$2; shift 2 ;;
    -h|--help) usage; exit 0 ;;
    *) echo "ERROR: unknown option $1" >&2; exit 2 ;;
  esac
done
[[ -s ${IN:-} && -s ${REF:-} && -n ${OUT:-} ]] || { usage >&2; exit 1; }
for x in bcftools bgzip tabix awk; do command -v "$x" >/dev/null || { echo "ERROR: missing $x" >&2; exit 1; }; done
bcftools view -h "$REF" | grep -q '^##INFO=<ID=AA,' || { echo "ERROR: reference VCF lacks INFO/AA: $REF" >&2; exit 1; }
mkdir -p "$(dirname "$OUT")"
TMP=$(mktemp -d "${TMPDIR:-/tmp}/gu-aa.XXXXXX")
trap 'rm -rf "$TMP"' EXIT

AA=$TMP/aa.tsv
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/AA\n' "$REF" | \
awk 'BEGIN{FS=OFS="\t"}
{
  aa=toupper($5); sub(/\|.*/,"",aa); sub(/,.*/,"",aa)
  ref=toupper($3); alt=toupper($4)
  if(aa~/^[ACGT]$/ && (aa==ref || aa==alt)) print $1,$2,$3,$4,aa
}' > "$AA"
[[ -s $AA ]] || { echo "ERROR: no usable INFO/AA records in $REF" >&2; exit 1; }
bgzip -c "$AA" > "$AA.gz"
tabix -f -s 1 -b 2 -e 2 "$AA.gz"
printf '##INFO=<ID=AA,Number=1,Type=String,Description="Ancestral allele transferred from build-matched 1000 Genomes INFO/AA by exact CHROM,POS,REF,ALT match">\n' > "$TMP/aa.hdr"

input_n=$(bcftools index -n "$IN" 2>/dev/null || echo 0)
PART=${OUT%.vcf.gz}.part.$$.vcf.gz
rm -f "$PART" "$PART.tbi"

# Remove any existing AA before annotating, so provenance is unambiguous.
bcftools annotate --threads "$THREADS" -x INFO/AA "$IN" -Ou 2>/dev/null | \
  bcftools annotate --threads "$THREADS" -a "$AA.gz" -h "$TMP/aa.hdr" -c CHROM,POS,REF,ALT,INFO/AA -Ou | \
  bcftools view --threads "$THREADS" -i 'INFO/AA!="." && (INFO/AA=REF || INFO/AA=ALT)' -Oz -o "$PART"
tabix -f -p vcf "$PART"
N=$(bcftools index -n "$PART")
if [[ $N -le 0 ]]; then
  echo "ERROR: no variants retained after exact AA transfer." >&2
  echo "Check genome build, chromosome naming (e.g. 22 vs chr22), REF orientation, and the 1KG source VCF." >&2
  exit 1
fi

mv -f "$PART" "$OUT"
mv -f "$PART.tbi" "$OUT.tbi"
rate=$(awk -v n="$N" -v d="$input_n" 'BEGIN{if(d>0)printf "%.6f",n/d;else print "NA"}')
printf 'metric\tvalue\ninput_target_variants\t%s\nretained_exact_AA_variants\t%s\nretention_fraction\t%s\nsource_1kg_vcf\t%s\n' \
  "$input_n" "$N" "$rate" "$REF" > "${OUT%.vcf.gz}.aa_qc.tsv"

echo "AA transfer complete: $OUT retained=$N/$input_n fraction=$rate"
awk -v r="$rate" 'BEGIN{if(r!="NA" && r+0<0.25) exit 1; exit 0}' || \
  echo "WARNING: <25% of target variants received exact 1KG AA annotation; verify build/REF/chromosome conventions." >&2
