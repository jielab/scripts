#!/usr/bin/env bash
# Prepare phased VCF input for ARG inference.
set -euo pipefail

usage() {
	cat <<'USAGE'
Usage:
  arg.prep_vcf.sh --in chrN.vcf.gz --out chrN.vcf.gz [options]

Options:
  --threads N           bcftools threads (default: 4)
  --missing-max X       require F_MISSING < X (default: 0.05)
  --mac-min N           require minor allele count >= N (default: 2)
  --tmp-root DIR        temporary work root (default: /tmp/gu-arg-vcf-$USER)
  --ref-fasta FILE      optional reference FASTA for REF validation/normalization
  --replace TRUE|FALSE  replace an existing output (default: FALSE)

The input must contain INFO/AA and phased GT data. INFO/AA values such as
"A|||" are normalized to "A". Sites whose AA is unknown or does not equal
REF or ALT are removed.
USAGE
}

die() { echo "ERROR: $*" >&2; exit 1; }
need() { command -v "$1" >/dev/null 2>&1 || die "required command not found: $1"; }
ensure_fill_tags() {
	if bcftools plugin -l 2>/dev/null | grep -Fx fill-tags >/dev/null; then return 0; fi
	if env -u BCFTOOLS_PLUGINS bcftools plugin -l 2>/dev/null | grep -Fx fill-tags >/dev/null; then
		unset BCFTOOLS_PLUGINS
		return 0
	fi
	return 1
}

in="" out="" threads=4 missing_max=0.05 mac_min=2
tmp_root=${ARG_TMP_ROOT:-/tmp/gu-arg-vcf-${USER:-user}}
ref_fasta="" replace=0 chrom="" build=37 x_psam=""
while (( $# )); do
	case "$1" in
		--in) [[ $# -ge 2 ]] || die "--in requires a file"; in=$2; shift 2 ;;
		--out) [[ $# -ge 2 ]] || die "--out requires a file"; out=$2; shift 2 ;;
		--threads) [[ $# -ge 2 && $2 =~ ^[1-9][0-9]*$ ]] || die "--threads requires a positive integer"; threads=$2; shift 2 ;;
		--missing-max) [[ $# -ge 2 ]] || die "--missing-max requires a number"; missing_max=$2; shift 2 ;;
		--mac-min) [[ $# -ge 2 && $2 =~ ^[1-9][0-9]*$ ]] || die "--mac-min requires a positive integer"; mac_min=$2; shift 2 ;;
		--tmp-root) [[ $# -ge 2 ]] || die "--tmp-root requires a directory"; tmp_root=$2; shift 2 ;;
		--ref-fasta) [[ $# -ge 2 ]] || die "--ref-fasta requires a file"; ref_fasta=$2; shift 2 ;;
		--chr) chrom=$2; shift 2 ;;
		--build) build=$2; shift 2 ;;
		--x-psam) x_psam=$2; shift 2 ;;
		--replace)
			[[ $# -ge 2 ]] || die "--replace requires TRUE or FALSE"
			case "${2,,}" in true|1|yes) replace=1 ;; false|0|no) replace=0 ;; *) die "--replace requires TRUE or FALSE" ;; esac
			shift 2 ;;
		-h|--help) usage; exit 0 ;;
		*) die "unknown option: $1" ;;
	esac
done

[[ -n $in && -n $out ]] || { usage >&2; exit 1; }
[[ $out == *.vcf.gz ]] || die "output must end in .vcf.gz: $out"
if [[ -s $out && -s $out.tbi && $replace -eq 0 ]]; then
	echo "SKIP: prepared ARG VCF already exists: $out"
	exit 0
fi
[[ -s $in ]] || die "input VCF is missing or empty: $in"
[[ -z $ref_fasta || -s $ref_fasta ]] || die "reference FASTA is missing: $ref_fasta"
for cmd in bcftools bgzip tabix awk sed stat; do need "$cmd"; done
ensure_fill_tags || die "bcftools fill-tags plugin is unavailable"
ARG_BCFTOOLS_BIN=$(type -P bcftools)
bcftools() {
	"$ARG_BCFTOOLS_BIN" "$@" 2> >(sed '/^Warning: cannot determine contig names given the \.csi index alone$/d' >&2)
}
bcftools view -h "$in" | grep '^##INFO=<ID=AA,' >/dev/null || die "input VCF has no INFO/AA header: $in"
python3 - "$missing_max" <<'PY' || die "--missing-max must be in [0,1)"
import sys
x=float(sys.argv[1]); assert 0 <= x < 1
PY

mkdir -p "$(dirname "$out")" "$tmp_root"
meta=$out.input.meta.tsv
completion_record() {
	printf 'input\t%s\nmissing_max\t%s\nmac_min\t%s\nref_fasta\t%s\n' \
		"$(stat -c '%n:%s:%Y' "$in")" "$missing_max" "$mac_min" \
		"$([[ -n $ref_fasta ]] && stat -c '%n:%s:%Y' "$ref_fasta" || printf none)"
}
lock=$out.lock
exec 9>"$lock"
if command -v flock >/dev/null 2>&1; then flock -n 9 || die "another process is preparing $out"; fi
work=$(mktemp -d "$tmp_root/argvcf.XXXXXX")
cleanup() { rm -rf -- "$work"; }
trap cleanup EXIT
trap 'exit 130' INT
trap 'exit 143' TERM

raw_bcf=$work/pass.bcf
tagged_bcf=$work/tagged.bcf
filtered_bcf=$work/filtered.bcf
aa_tsv=$work/aa.tsv
aa_tsv_gz=$work/aa.tsv.gz
aa_hdr=$work/aa.hdr
info_only_bcf=$work/info-only.bcf
# Keep the standard .vcf.gz suffix so htslib/bcftools can discover the
# temporary file's .tbi index during validation.
out_part=${out%.vcf.gz}.part.$$.vcf.gz
rm -f "$out_part" "$out_part.tbi"

echo "[ARG-VCF] PASS + biallelic SNP selection: $in"
source_vcf=$in
view_args=(view --threads "$threads" -f PASS -m2 -M2 -v snps -Ob -o "$raw_bcf")
if [[ $chrom == X ]]; then
	# Keep filtering in bcftools; validate biological ploidy after variant QC,
	# avoiding a Python pass over billions of raw, mostly rare genotypes.
	x_contig=$(bcftools index -s "$in" | awk '$1=="X" || $1=="chrX" || $1=="23" {print $1}')
	[[ -n $x_contig && $x_contig != *$'\n'* ]] || die 'X input must have exactly one X contig'
	x_lo=2699521; x_hi=154931043
	if [[ $build == 38 ]]; then x_lo=2781480; x_hi=155701382; fi
	view_args+=(-t "$x_contig:$x_lo-$x_hi")
fi
if [[ -n $ref_fasta ]]; then
	# SNP-only data need no left alignment, but this excludes REF mismatches and
	# removes duplicate records at the same position before downstream filtering.
	bcftools "${view_args[@]}" "$source_vcf"
	bcftools index -f "$raw_bcf"
	bcftools norm --threads "$threads" -f "$ref_fasta" -c x -d all -Ob -o "$work/norm.bcf" "$raw_bcf"
	mv "$work/norm.bcf" "$raw_bcf"
	bcftools index -f "$raw_bcf"
else
	bcftools "${view_args[@]}" "$source_vcf"
	bcftools index -f "$raw_bcf"
	bcftools norm --threads "$threads" -d all -Ob -o "$work/dedup.bcf" "$raw_bcf"
	mv "$work/dedup.bcf" "$raw_bcf"
	bcftools index -f "$raw_bcf"
fi

echo "[ARG-VCF] Recalculate allele counts and missingness from GT"
# Some bcftools fill-tags versions do not provide MAC. For the biallelic sites
# selected above, MAC >= N is equivalent to AC >= N and AN - AC >= N.
bcftools +fill-tags "$raw_bcf" --threads "$threads" -Ob -o "$tagged_bcf" -- -t AC,AN,F_MISSING
bcftools index -f "$tagged_bcf"
bcftools view --threads "$threads" \
	-i "INFO/F_MISSING<=${missing_max} && INFO/AC>=${mac_min} && INFO/AN-INFO/AC>=${mac_min} && INFO/AA!='.'" \
	-Ob -o "$filtered_bcf" "$tagged_bcf"
bcftools index -f "$filtered_bcf"

echo "[ARG-VCF] Normalize INFO/AA and retain only AA matching REF or ALT"
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/AA\n' "$filtered_bcf" | \
awk 'BEGIN{FS=OFS="\t"}
{
	aa=toupper($5); sub(/\|.*/,"",aa); sub(/,.*/,"",aa)
	ref=toupper($3); alt=toupper($4)
	if (aa ~ /^[ACGT]$/ && (aa==ref || aa==alt)) print $1,$2,$3,$4,aa
}' > "$aa_tsv"
[[ -s $aa_tsv ]] || die "no variants retained after INFO/AA normalization: $in"
bgzip -c "$aa_tsv" > "$aa_tsv_gz"
tabix -f -s 1 -b 2 -e 2 "$aa_tsv_gz"
printf '##INFO=<ID=AA,Number=1,Type=String,Description="Normalized ancestral allele parsed from original INFO/AA">\n' > "$aa_hdr"

# Apply one documented complement operation per annotation family. This avoids
# relying on ambiguous parsing of multiple caret-prefixed families in one list.
bcftools annotate --threads "$threads" -x INFO/AA "$filtered_bcf" -Ou | \
	bcftools annotate --threads "$threads" -a "$aa_tsv_gz" -h "$aa_hdr" \
		-c CHROM,POS,REF,ALT,INFO/AA -Ou | \
	bcftools view --threads "$threads" -i 'INFO/AA!="."' -Ou | \
	bcftools annotate --threads "$threads" -x '^INFO/AA' -Ob -o "$info_only_bcf"
format_remove=$(bcftools view -h "$info_only_bcf" | sed -n 's/^##FORMAT=<ID=\([^,>]*\).*/FORMAT\/\1/p' | grep -Fvx FORMAT/GT | paste -sd, - || true)
if [[ -n $format_remove ]]; then
	bcftools annotate --threads "$threads" -x "$format_remove" -Oz -o "$out_part" "$info_only_bcf"
else
	bcftools view --threads "$threads" -Oz -o "$out_part" "$info_only_bcf"
fi
tabix -f -p vcf "$out_part"

[[ $(bcftools query -l "$out_part" | wc -l) -gt 0 ]] || die "output has no samples: $out_part"
final_n=$(bcftools index -n "$out_part")
[[ $final_n -gt 0 ]] || die "output has no variants: $out_part"
unexpected_info=$(bcftools view -h "$out_part" | sed -n 's/^##INFO=<ID=\([^,>]*\).*/\1/p' | grep -Fvx AA || true)
unexpected_format=$(bcftools view -h "$out_part" | sed -n 's/^##FORMAT=<ID=\([^,>]*\).*/\1/p' | grep -Fvx GT || true)
[[ -z $unexpected_info ]] || die "output retained unexpected INFO annotations: $(tr '\n' ',' <<< "$unexpected_info")"
[[ -z $unexpected_format ]] || die "output retained unexpected FORMAT annotations: $(tr '\n' ',' <<< "$unexpected_format")"
bad_aa=$(bcftools query -f '%REF\t%ALT\t%INFO/AA\n' "$out_part" | \
	awk 'BEGIN{n=0} $3!=$1 && $3!=$2{n++} END{print n+0}')
[[ $bad_aa -eq 0 ]] || die "output contains $bad_aa sites with AA not matching REF/ALT"

validation_args=(--chr "$chrom" --build "$build")
[[ -z $x_psam ]] || validation_args+=(--x-psam "$x_psam")
python3 "$(dirname -- "${BASH_SOURCE[0]}")/arg.validate_vcf.py" "$out_part" "${validation_args[@]}"
mv -f "$out_part" "$out"
rm -f "$out.csi"
mv -f "$out_part.tbi" "$out.tbi"
completion_record > "$meta"
raw_n=$(bcftools index -n "$raw_bcf")
filtered_n=$(bcftools index -n "$filtered_bcf")
sample_n=$(bcftools query -l "$out" | wc -l)
printf 'metric\tvalue\ninput_pass_biallelic_snp\t%s\npost_mac_missing\t%s\npost_valid_AA\t%s\nsamples\t%s\nmissing_max\t%s\nmac_min\t%s\n' \
	"$raw_n" "$filtered_n" "$final_n" "$sample_n" "$missing_max" "$mac_min" > "${out%.vcf.gz}.qc.tsv"
echo "[ARG-VCF] DONE: $out variants=$final_n samples=$sample_n"
