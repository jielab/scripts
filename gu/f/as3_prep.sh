#!/usr/bin/env bash
# ArchaicSeeker3 preprocessing with explicit target/reference separation,
# all-target support, cross-fitted target/reference overlap handling,
# and optional experimental chrX PAR/non-PAR units.
set -euo pipefail

need() { command -v "$1" >/dev/null 2>&1 || { echo "FATAL ERROR: missing command: $1" >&2; exit 1; }; }
ensure_fill_tags() {
    if bcftools plugin -l 2>/dev/null | grep -Fx fill-tags >/dev/null; then return 0; fi
    if env -u BCFTOOLS_PLUGINS bcftools plugin -l 2>/dev/null | grep -Fx fill-tags >/dev/null; then
        unset BCFTOOLS_PLUGINS
        return 0
    fi
    return 1
}
for cmd in bcftools bgzip tabix awk sort comm md5sum python3; do need "$cmd"; done
ensure_fill_tags || { echo "FATAL ERROR: bcftools fill-tags plugin is unavailable" >&2; exit 1; }

script_dir=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
dir0=${dir0:-/mnt/d}
dir_ref=${dir_ref:-$dir0/data.BIG/refGen}
dirfunc=${dirfunc:-$dir0/scripts/0f}
phe_sh=${phe_sh:-$dirfunc/0phe.f.sh}
[[ -s $phe_sh ]] || { echo "FATAL ERROR: shared helper not found: $phe_sh" >&2; exit 1; }
# shellcheck source=/dev/null
source "$phe_sh"

# Parse preparation options passed by gu.sh. Environment variables remain valid
# overrides for advanced/custom layouts.
dir_as3=${dir_archaic:-/mnt/d/data.BIG/refGen/archaic/GRCH${GRCH:-38}}
while (( $# )); do
    case "$1" in
        --dir-as3) [[ $# -ge 2 ]] || { echo "FATAL ERROR: --dir-as3 requires a directory" >&2; exit 1; }; dir_as3=$2; shift 2 ;;
        --out) [[ $# -ge 2 ]] || { echo "FATAL ERROR: --out requires a directory" >&2; exit 1; }; export AS3_DATA_OUT=$2; shift 2 ;;
        --chr|--chrs) [[ $# -ge 2 ]] || { echo "FATAL ERROR: $1 requires a chromosome list" >&2; exit 1; }; export AS3_CHRS=${2//,/ }; shift 2 ;;
        --threads) [[ $# -ge 2 && $2 =~ ^[1-9][0-9]*$ ]] || { echo "FATAL ERROR: --threads requires a positive integer" >&2; exit 1; }; export AS3_THREADS_PER_JOB=$2; shift 2 ;;
        --max-procs) [[ $# -ge 2 && $2 =~ ^[1-9][0-9]*$ ]] || { echo "FATAL ERROR: --max-procs requires a positive integer" >&2; exit 1; }; export AS3_MAX_PROCS=$2; shift 2 ;;
        --tmp-root) [[ $# -ge 2 ]] || { echo "FATAL ERROR: --tmp-root requires a directory" >&2; exit 1; }; export AS3_TMP_ROOT=$2; shift 2 ;;
        --target-vcf-dir) [[ $# -ge 2 ]] || { echo "FATAL ERROR: $1 requires a directory" >&2; exit 1; }; export AS3_TARGET_VCF_DIR=$2; shift 2 ;;
        --target-vcf-template) [[ $# -ge 2 ]] || { echo "FATAL ERROR: $1 requires a template" >&2; exit 1; }; export AS3_TARGET_VCF_TEMPLATE=$2; shift 2 ;;
        --afr-ref-vcf-dir) [[ $# -ge 2 ]] || { echo "FATAL ERROR: $1 requires a directory" >&2; exit 1; }; export AS3_AFR_REFERENCE_VCF_DIR=$2; shift 2 ;;
        --afr-ref-vcf-template) [[ $# -ge 2 ]] || { echo "FATAL ERROR: $1 requires a template" >&2; exit 1; }; export AS3_AFR_REFERENCE_VCF_TEMPLATE=$2; shift 2 ;;
        --target-samples) [[ $# -ge 2 ]] || { echo "FATAL ERROR: $1 requires FILE or ALL" >&2; exit 1; }; export AS3_TARGET_SAMPLES_FILE=$2; shift 2 ;;
        --afr-ref-samples) [[ $# -ge 2 ]] || { echo "FATAL ERROR: $1 requires a file" >&2; exit 1; }; export AS3_AFR_REFERENCE_SAMPLES_FILE=$2; shift 2 ;;
        --target-panel) [[ $# -ge 2 ]] || { echo "FATAL ERROR: $1 requires a file" >&2; exit 1; }; export AS3_TARGET_PANEL_FILE=$2; shift 2 ;;
        --afr-ref-panel) [[ $# -ge 2 ]] || { echo "FATAL ERROR: $1 requires a file" >&2; exit 1; }; export AS3_AFR_REFERENCE_PANEL_FILE=$2; shift 2 ;;
        --target-sample-col) [[ $# -ge 2 ]] || { echo "FATAL ERROR: $1 requires a column" >&2; exit 1; }; export AS3_TARGET_SAMPLE_COL=$2; shift 2 ;;
        --target-sex-col) [[ $# -ge 2 ]] || { echo "FATAL ERROR: $1 requires a column" >&2; exit 1; }; export AS3_TARGET_SEX_COL=$2; shift 2 ;;
        --afr-ref-sample-col) [[ $# -ge 2 ]] || { echo "FATAL ERROR: $1 requires a column" >&2; exit 1; }; export AS3_AFR_REFERENCE_SAMPLE_COL=$2; shift 2 ;;
        --afr-ref-sex-col) [[ $# -ge 2 ]] || { echo "FATAL ERROR: $1 requires a column" >&2; exit 1; }; export AS3_AFR_REFERENCE_SEX_COL=$2; shift 2 ;;
        --genome-build) [[ $# -ge 2 ]] || { echo "FATAL ERROR: $1 requires b37 or b38" >&2; exit 1; }; export AS3_GENOME_BUILD=$2; shift 2 ;;
        --allow-chrx) export AS3_ALLOW_CHRX=1; shift ;;
        --overlap-policy) [[ $# -ge 2 ]] || { echo "FATAL ERROR: $1 requires error or crossfit" >&2; exit 1; }; export AS3_OVERLAP_POLICY=$2; shift 2 ;;
        --overlap-folds) [[ $# -ge 2 && $2 =~ ^[0-9]+$ && $2 -ge 2 ]] || { echo "FATAL ERROR: $1 requires an integer >=2" >&2; exit 1; }; export AS3_OVERLAP_FOLDS=$2; shift 2 ;;
        --check-only) export AS3_PREP_CHECK_ONLY=1; shift ;;
        *) echo "FATAL ERROR: unknown AS3 preparation option: $1" >&2; exit 1 ;;
    esac
done

default_1kg_panel=${dir_1kg:-/mnt/d/data.BIG/refGen/1kg/GRCH${GRCH:-38}}/vcf/samples_v3.ALL.panel
export AS3_TARGET_PANEL_FILE=${AS3_TARGET_PANEL_FILE:-$default_1kg_panel}
export AS3_AFR_REFERENCE_PANEL_FILE=${AS3_AFR_REFERENCE_PANEL_FILE:-$AS3_TARGET_PANEL_FILE}
export AS3_TARGET_SAMPLE_COL=${AS3_TARGET_SAMPLE_COL:-sample}
export AS3_TARGET_SEX_COL=${AS3_TARGET_SEX_COL:-gender}
export AS3_AFR_REFERENCE_SAMPLE_COL=${AS3_AFR_REFERENCE_SAMPLE_COL:-sample}
export AS3_AFR_REFERENCE_SEX_COL=${AS3_AFR_REFERENCE_SEX_COL:-gender}
if [[ " ${AS3_CHRS:-} " =~ [[:space:]](X|x|23)([[:space:]]|$) ]]; then export AS3_ALLOW_CHRX=1; fi
export AS3_DATA_OUT=${AS3_DATA_OUT:-$dir_as3/as3/preprocessed}

MAX_PROCS=${AS3_MAX_PROCS:-1}
THREADS=${AS3_THREADS_PER_JOB:-6}
TMP_ROOT=${AS3_TMP_ROOT:-/tmp/gu-as3-${USER:-user}}
OUT=${AS3_DATA_OUT:-$dir_as3/as3/preprocessed}
REF_FASTA=${AS3_REF_GENOME:-$dir_ref/fasta/GRCH${GRCH:-38}.fasta}
TARGET_DIR=${AS3_TARGET_VCF_DIR:-${AS3_MODERN_HUMAN_VCF_DIR:-${AS3_MODERN_VCF_DIR:-${dir_1kg:-$dir_ref/1kg/GRCH${GRCH:-38}}/vcf}}}
AFR_DIR=${AS3_AFR_REFERENCE_VCF_DIR:-$TARGET_DIR}
DEN_DIR=${AS3_DENISOVAN_VCF_DIR:-$dir_as3/Denisova}
NEA_DIR=${AS3_NEANDERTHAL_VCF_DIR:-$dir_as3/Altai}
MASK_DIR=${AS3_MASK_BED_DIR:-$dir_as3/high_quality_region_GRCh38}
TARGET_TEMPLATE=${AS3_TARGET_VCF_TEMPLATE:-}
[[ -n $TARGET_TEMPLATE ]] || TARGET_TEMPLATE='chr{CHR}.vcf.gz'
AFR_TEMPLATE=${AS3_AFR_REFERENCE_VCF_TEMPLATE:-}
[[ -n $AFR_TEMPLATE ]] || AFR_TEMPLATE=$TARGET_TEMPLATE
DEN_TEMPLATE=${AS3_DENISOVAN_VCF_TEMPLATE:-}
[[ -n $DEN_TEMPLATE ]] || DEN_TEMPLATE='output_se.hg38.chr{CHR}.vcf.gz'
NEA_TEMPLATE=${AS3_NEANDERTHAL_VCF_TEMPLATE:-}
[[ -n $NEA_TEMPLATE ]] || NEA_TEMPLATE='output_se.hg38.chr{CHR}.vcf.gz'
TARGET_SPEC=${AS3_TARGET_SAMPLES_FILE:-ALL}
AFR_SPEC=${AS3_AFR_REFERENCE_SAMPLES_FILE:-${AS3_YRI_SAMPLES_FILE:-}}
SAMPLE_LISTS_DIR=${AS3_SAMPLE_LISTS_DIR:-}
TARGET_PANEL=${AS3_TARGET_PANEL_FILE:-}
AFR_PANEL=${AS3_AFR_REFERENCE_PANEL_FILE:-$TARGET_PANEL}
TARGET_SAMPLE_COL=${AS3_TARGET_SAMPLE_COL:-}
TARGET_SEX_COL=${AS3_TARGET_SEX_COL:-}
AFR_SAMPLE_COL=${AS3_AFR_REFERENCE_SAMPLE_COL:-}
AFR_SEX_COL=${AS3_AFR_REFERENCE_SEX_COL:-}
GENOME_BUILD=${AS3_GENOME_BUILD:-b${GRCH:-38}}
MAF_MIN=${AS3_TARGET_MAF_MIN:-0.001}
OVERLAP_POLICY=${AS3_OVERLAP_POLICY:-crossfit}
OVERLAP_FOLDS=${AS3_OVERLAP_FOLDS:-5}
ALLOW_X=${AS3_ALLOW_CHRX:-0}
CHECK_ONLY=${AS3_PREP_CHECK_ONLY:-0}
MASK_PREFIXES_STR=${AS3_MASK_PREFIXES:-"N D"}
read -r -a MASK_PREFIXES <<< "$MASK_PREFIXES_STR"
read -r -a REQUESTED <<< "${AS3_CHRS:-$(seq 1 22)}"

[[ $MAX_PROCS =~ ^[1-9][0-9]*$ ]] || { echo "FATAL ERROR: AS3_MAX_PROCS must be a positive integer" >&2; exit 1; }
[[ $THREADS =~ ^[1-9][0-9]*$ ]] || { echo "FATAL ERROR: AS3_THREADS_PER_JOB must be a positive integer" >&2; exit 1; }
[[ $OVERLAP_FOLDS =~ ^[2-9][0-9]*$|^[2-9]$ ]] || { echo "FATAL ERROR: AS3_OVERLAP_FOLDS must be >=2" >&2; exit 1; }
[[ $GENOME_BUILD == b37 || $GENOME_BUILD == b38 ]] || { echo "FATAL ERROR: AS3_GENOME_BUILD must be b37 or b38" >&2; exit 1; }
[[ -z ${GRCH:-} || $GENOME_BUILD == b$GRCH ]] || { echo "FATAL ERROR: GRCH=$GRCH conflicts with AS3_GENOME_BUILD=$GENOME_BUILD" >&2; exit 1; }
[[ $GENOME_BUILD == b38 ]] || echo "WARNING: the public AS3 preprocessing examples use GRCh38; b37 is an explicit nonstandard analysis and requires build-matched archaic VCFs, FASTA, and masks" >&2
[[ $OVERLAP_POLICY == error || $OVERLAP_POLICY == crossfit ]] || { echo "FATAL ERROR: AS3_OVERLAP_POLICY must be error or crossfit" >&2; exit 1; }
[[ $ALLOW_X == 0 || $ALLOW_X == 1 ]] || { echo "FATAL ERROR: AS3_ALLOW_CHRX must be 0 or 1" >&2; exit 1; }
[[ $CHECK_ONLY == 0 || $CHECK_ONLY == 1 ]] || { echo "FATAL ERROR: AS3_PREP_CHECK_ONLY must be 0 or 1" >&2; exit 1; }
[[ -s $REF_FASTA ]] || { echo "FATAL ERROR: AS3_REF_GENOME missing: $REF_FASTA" >&2; exit 1; }
[[ -s ${REF_FASTA}.fai ]] || { echo "FATAL ERROR: FASTA index missing: ${REF_FASTA}.fai" >&2; exit 1; }
for d in "$TARGET_DIR" "$AFR_DIR" "$DEN_DIR" "$NEA_DIR" "$MASK_DIR"; do
    [[ -d $d ]] || { echo "FATAL ERROR: required directory missing: $d" >&2; exit 1; }
done
[[ -n $AFR_SPEC ]] || { echo "FATAL ERROR: set AS3_AFR_REFERENCE_SAMPLES_FILE" >&2; exit 1; }

resolve_list_spec() {
    local spec=$1
    if [[ $spec == ALL ]]; then printf 'ALL\n'; return 0; fi
    if [[ -s $spec ]]; then readlink -f "$spec"; return 0; fi
    if [[ -n $SAMPLE_LISTS_DIR && -s $SAMPLE_LISTS_DIR/$spec ]]; then readlink -f "$SAMPLE_LISTS_DIR/$spec"; return 0; fi
    echo "FATAL ERROR: sample list not found: $spec" >&2
    return 1
}
TARGET_SPEC=$(resolve_list_spec "$TARGET_SPEC")
AFR_SPEC=$(resolve_list_spec "$AFR_SPEC")

mkdir -p "$OUT"/{Final_Target_VCFs,Final_Ref_VCFs,Reference_Maps,Master_Target_Sites,samples/plans,qc,log} "$TMP_ROOT"
[[ $TMP_ROOT == /mnt/* ]] && echo "WARNING: use a Linux filesystem such as /tmp for AS3_TMP_ROOT" >&2

list_ids() { awk 'NF && $1 !~ /^#/ {print $1}' "$1"; }
count_lines() { awk 'NF{n++} END{print n+0}' "$1"; }
file_stamp() { stat -c '%n:%s:%Y' "$1"; }
valid_vcf() {
    local f=$1
    [[ -s $f && -s $f.tbi ]] || return 1
    bcftools view -h "$f" >/dev/null 2>&1 || return 1
    [[ $(bcftools index -n "$f" 2>/dev/null || echo 0) -gt 0 ]]
}
render() {
    local dir=$1 tpl=$2 chr=$3 x
    x=${tpl//\{CHR\}/$chr}
    x=${x//\{chr\}/$chr}
    [[ $x == /* ]] && printf '%s\n' "$x" || printf '%s/%s\n' "$dir" "$x"
}
resolve_vcf() {
    local role=$1 dir=$2 tpl=$3 chr=$4 f
    local -a candidates=()
    f=$(render "$dir" "$tpl" "$chr")
    [[ -s $f ]] && { printf '%s\n' "$f"; return 0; }
    case "$role" in
        target|African-reference)
            candidates=(
                "$dir/chr${chr}.vcf.gz"
                "$dir/KGP.GRCh38.PASS.snps.biallelic.shapeit4.${chr}.part.vcf.gz"
                "$dir/KGP.GRCh37.PASS.snps.biallelic.shapeit4.${chr}.part.vcf.gz"
            )
            ;;
        Denisovan)
            candidates=(
                "$dir/Den.hg38.${chr}.part.vcf.gz"
                "$dir/output_se.hg38.chr${chr}.vcf.gz"
                "$dir/Den.hg19.${chr}.part.vcf.gz"
                "$dir/chr${chr}_mq25_mapab100.vcf.gz"
                "$dir/chr${chr}.Den25.L35MQ25.B30.map35_100.vcf.gz"
                "$dir/chr${chr}.vcf.gz"
            )
            ;;
        Neanderthal)
            candidates=(
                "$dir/Nean.hg38.${chr}.part.vcf.gz"
                "$dir/output_se.hg38.chr${chr}.vcf.gz"
                "$dir/Nean.hg19.${chr}.part.vcf.gz"
                "$dir/chr${chr}.noRB.vcf.gz"
                "$dir/chr${chr}_mq25_mapab100.vcf.gz"
                "$dir/chr${chr}.vcf.gz"
            )
            ;;
        *)
            echo "FATAL ERROR: unknown VCF role: $role" >&2
            return 1
            ;;
    esac
    for f in "${candidates[@]}"; do
        [[ -s $f ]] && { printf '%s\n' "$f"; return 0; }
    done
    echo "FATAL ERROR: cannot resolve $role VCF for chr$chr under $dir using template $tpl" >&2
    return 1
}
source_chr() { case "$1" in XPAR|XNONPAR_F|XNONPAR_M) echo X ;; *) echo "$1" ;; esac; }
unit_region() {
    local var
    case "$1" in
        XPAR) var=X_PAR_${GENOME_BUILD}; printf '%s\n' "${!var}" ;;
        XNONPAR_F|XNONPAR_M) var=X_NONPAR_${GENOME_BUILD}; printf '%s\n' "${!var}" ;;
        *) printf '\n' ;;
    esac
}
first_contig() {
    python3 - "$1" <<'PYVCF'
import gzip
import sys
path = sys.argv[1]
with open(path, "rb") as raw:
    magic = raw.read(2)
opener = gzip.open if magic == b"\x1f\x8b" else open
with opener(path, "rt", encoding="utf-8", errors="replace") as handle:
    for line in handle:
        if line and not line.startswith("#"):
            print(line.split("\t", 1)[0].strip())
            break
    else:
        raise SystemExit(f"no VCF records found: {path}")
PYVCF
}
adapt_region() {
    local vcf=$1 region=$2 contig
    [[ -n $region ]] || { printf '\n'; return 0; }
    contig=$(first_contig "$vcf")
    [[ -n $contig ]] || { echo "FATAL ERROR: cannot determine contig from $vcf" >&2; return 1; }
    printf '%s\n' "${region//X:/$contig:}"
}
check_members() {
    local list=$1 vcf=$2 label=$3 work=$4
    list_ids "$list" | sort -u > "$work/wanted.txt"
    bcftools query -l "$vcf" | sort -u > "$work/available.txt"
    comm -23 "$work/wanted.txt" "$work/available.txt" > "$work/missing.txt"
    if [[ -s $work/missing.txt ]]; then
        echo "ERROR [$label]: sample IDs absent from $vcf" >&2
        head -20 "$work/missing.txt" >&2
        return 1
    fi
}
make_base_list() {
    local spec=$1 vcf=$2 out=$3 label=$4 work=$5
    if [[ $spec == ALL ]]; then bcftools query -l "$vcf" > "$out"; else list_ids "$spec" > "$out"; fi
    [[ -s $out ]] || { echo "FATAL ERROR: empty $label sample list" >&2; return 1; }
    [[ $(sort "$out" | uniq -d | wc -l) -eq 0 ]] || { echo "FATAL ERROR: duplicate IDs in $label sample list" >&2; return 1; }
    check_members "$out" "$vcf" "$label" "$work"
}
normalize_sex_panel() {
    local panel=$1 out=$2 listdir=$3 sample_col=$4 sex_col=$5 label=$6
    [[ -s $panel ]] || { echo "FATAL ERROR: $label panel missing: $panel" >&2; return 1; }
    local -a args=(--input "$panel" --output "$out" --list-dir "$listdir" --allow-missing-pop --require-sex)
    [[ -n $sample_col ]] && args+=(--sample-col "$sample_col")
    [[ -n $sex_col ]] && args+=(--sex-col "$sex_col")
    python3 "$script_dir/sample_panel.py" "${args[@]}"
}
filter_by_sex() {
    local base=$1 panel=$2 sex=$3 out=$4
    awk -F'\t' -v sex="$sex" 'NR==FNR{want[$1]=1; next} NR>1 && $4==sex && want[$1]{print $1}' "$base" "$panel" > "$out"
    [[ -s $out ]] || { echo "FATAL ERROR: no $sex samples remain after intersecting $base with $panel" >&2; return 1; }
}
check_panel_covers_list() {
    local list=$1 panel=$2 label=$3 work=$4
    list_ids "$list" | sort -u > "$work/wanted.txt"
    awk -F'\t' 'NR>1{print $1}' "$panel" | sort -u > "$work/panel.txt"
    comm -23 "$work/wanted.txt" "$work/panel.txt" > "$work/missing.txt"
    if [[ -s $work/missing.txt ]]; then
        echo "FATAL ERROR: $label sex panel omits samples; examples:" >&2
        head -20 "$work/missing.txt" >&2
        return 1
    fi
}

# Expand requested chromosomes into analysis units and remove duplicates.
units=()
for chr in "${REQUESTED[@]}"; do
    chr=${chr#chr}
    case "$chr" in
        X|23)
            [[ $ALLOW_X == 1 ]] || { echo "FATAL ERROR: chrX requested but AS3_ALLOW_CHRX=0; use --allow-chrx for an explicitly experimental analysis" >&2; exit 1; }
            units+=(XPAR XNONPAR_F XNONPAR_M)
            ;;
        [1-9]|1[0-9]|2[0-2]) units+=("$chr") ;;
        *) echo "FATAL ERROR: unsupported AS3 chromosome: $chr" >&2; exit 1 ;;
    esac
done
mapfile -t units < <(printf '%s\n' "${units[@]}" | awk '!seen[$0]++')
(( ${#units[@]} > 0 )) || { echo "FATAL ERROR: no AS3 analysis units requested" >&2; exit 1; }
if printf '%s\n' "${units[@]}" | grep -q '^X'; then
    [[ -s $TARGET_PANEL ]] || { echo "FATAL ERROR: --target-panel is required for chrX" >&2; exit 1; }
    [[ -s $AFR_PANEL ]] || { echo "FATAL ERROR: --afr-ref-panel is required for chrX" >&2; exit 1; }
    echo "WARNING: AS3 chrX input support is experimental; the public pretrained model is not documented as chrX-calibrated." >&2
    normalize_sex_panel "$TARGET_PANEL" "$OUT/samples/target.panel.normalized.tsv" "$OUT/samples/target.panel.lists" "$TARGET_SAMPLE_COL" "$TARGET_SEX_COL" target
    normalize_sex_panel "$AFR_PANEL" "$OUT/samples/afr.panel.normalized.tsv" "$OUT/samples/afr.panel.lists" "$AFR_SAMPLE_COL" "$AFR_SEX_COL" African-reference
fi

# Build exact target/reference universes and a non-overlapping task plan for each unit.
printf 'panel_id\tunit\ttarget_list\tafr_reference_list\n' > "$OUT/tasks.tsv"
printf 'unit\tsource_chr\tregion\tn_target_universe\tn_afr_reference_universe\tn_panels\n' > "$OUT/samples/unit_summary.tsv"
for unit in "${units[@]}"; do
    chr=$(source_chr "$unit")
    target_vcf=$(resolve_vcf target "$TARGET_DIR" "$TARGET_TEMPLATE" "$chr")
    afr_vcf=$(resolve_vcf African-reference "$AFR_DIR" "$AFR_TEMPLATE" "$chr")
    work=$(mktemp -d "$TMP_ROOT/as3_plan_${unit}.XXXXXX")
    trap 'rm -rf -- "$work"' EXIT
    target_base=$OUT/samples/${unit}.target.base.txt
    afr_base=$OUT/samples/${unit}.afr.base.txt
    make_base_list "$TARGET_SPEC" "$target_vcf" "$target_base" "$unit target" "$work"
    make_base_list "$AFR_SPEC" "$afr_vcf" "$afr_base" "$unit African reference" "$work"
    if [[ $unit == XPAR || $unit == XNONPAR_F || $unit == XNONPAR_M ]]; then
        check_panel_covers_list "$target_base" "$OUT/samples/target.panel.normalized.tsv" target "$work"
        check_panel_covers_list "$afr_base" "$OUT/samples/afr.panel.normalized.tsv" African-reference "$work"
    fi
    target_unit=$OUT/samples/${unit}.target.txt
    afr_unit=$OUT/samples/${unit}.afr.txt
    case "$unit" in
        XNONPAR_F)
            filter_by_sex "$target_base" "$OUT/samples/target.panel.normalized.tsv" female "$target_unit"
            filter_by_sex "$afr_base" "$OUT/samples/afr.panel.normalized.tsv" female "$afr_unit"
            ;;
        XNONPAR_M)
            filter_by_sex "$target_base" "$OUT/samples/target.panel.normalized.tsv" male "$target_unit"
            # Female African references retain two observed non-PAR X haplotypes and avoid
            # double-weighting duplicated male reference haplotypes.
            filter_by_sex "$afr_base" "$OUT/samples/afr.panel.normalized.tsv" female "$afr_unit"
            ;;
        *) cp "$target_base" "$target_unit"; cp "$afr_base" "$afr_unit" ;;
    esac
    plan_dir=$OUT/samples/plans/$unit
    rm -rf "$plan_dir"
    python3 "$script_dir/as3_plan.py" --target "$target_unit" --reference "$afr_unit" --outdir "$plan_dir" --policy "$OVERLAP_POLICY" --folds "$OVERLAP_FOLDS"
    nplan=0
    while IFS=$'\t' read -r panel tlist rlist nt nr; do
        [[ $panel == panel_id ]] && continue
        printf '%s\t%s\t%s\t%s\n' "$panel" "$unit" "$tlist" "$rlist" >> "$OUT/tasks.tsv"
        nplan=$((nplan + 1))
    done < "$plan_dir/panels.tsv"
    (( nplan > 0 )) || { echo "FATAL ERROR: no AS3 panels created for $unit" >&2; exit 1; }

    # Every target in the unit must occur exactly once across the panel target lists.
    tail -n +2 "$plan_dir/panels.tsv" | cut -f2 | while IFS= read -r f; do cat "$f"; done | sort > "$work/planned.txt"
    sort "$target_unit" > "$work/expected.txt"
    cmp -s "$work/planned.txt" "$work/expected.txt" || { echo "FATAL ERROR: target coverage mismatch in AS3 plan for $unit" >&2; exit 1; }
    [[ $(uniq -d "$work/planned.txt" | wc -l) -eq 0 ]] || { echo "FATAL ERROR: duplicated target assignment in AS3 plan for $unit" >&2; exit 1; }
    printf '%s\t%s\t%s\t%s\t%s\t%s\n' "$unit" "$chr" "$(unit_region "$unit")" "$(count_lines "$target_unit")" "$(count_lines "$afr_unit")" "$nplan" >> "$OUT/samples/unit_summary.tsv"
    rm -rf "$work"
    trap - EXIT
done

# Resolve all remaining inputs before any expensive preprocessing.
for unit in "${units[@]}"; do
    chr=$(source_chr "$unit")
    resolve_vcf Denisovan "$DEN_DIR" "$DEN_TEMPLATE" "$chr" >/dev/null
    resolve_vcf Neanderthal "$NEA_DIR" "$NEA_TEMPLATE" "$chr" >/dev/null
    for prefix in "${MASK_PREFIXES[@]}"; do
        [[ -n $prefix ]] || continue
        [[ -s $MASK_DIR/${prefix}_chr${chr}.bed ]] || { echo "FATAL ERROR: mask missing: $MASK_DIR/${prefix}_chr${chr}.bed" >&2; exit 1; }
    done
done
if [[ $CHECK_ONLY == 1 ]]; then
    echo "CHECK PASSED: AS3 preprocessing plan=$OUT/tasks.tsv units=${units[*]} target_policy=$TARGET_SPEC overlap_policy=$OVERLAP_POLICY"
    exit 0
fi

input_signature() {
    local unit=$1 source=$2 list=$3
    {
        printf 'unit\t%s\nsource\t%s\nlist_md5\t%s\nregion\t%s\ngenome_build\t%s\nmaf_min\t%s\nref_fasta\t%s\n' \
            "$unit" "$(file_stamp "$source")" "$(md5sum "$list" | awk '{print $1}')" "$(unit_region "$unit")" "$GENOME_BUILD" "$MAF_MIN" "$(file_stamp "$REF_FASTA")"
    }
}
prepare_master() (
    set -euo pipefail
    local unit=$1 chr target_src target_list region adapted out sig work tagged
    chr=$(source_chr "$unit")
    target_src=$(resolve_vcf target "$TARGET_DIR" "$TARGET_TEMPLATE" "$chr")
    target_list=$OUT/samples/${unit}.target.base.txt
    out=$OUT/Master_Target_Sites/target_sites.chr${unit}.vcf.gz
    sig=$OUT/Master_Target_Sites/target_sites.chr${unit}.input.tsv
    if valid_vcf "$out" && [[ -s $sig ]] && cmp -s <(input_signature "$unit" "$target_src" "$target_list") "$sig"; then
        echo "[AS3 master chr${unit}] SKIP"
        exit 0
    fi
    work=$(mktemp -d "$TMP_ROOT/as3_master_${unit}.XXXXXX")
    trap 'rm -rf -- "$work"' EXIT
    region=$(adapt_region "$target_src" "$(unit_region "$unit")")
    local -a view_args=(view --threads "$THREADS" -f PASS -m2 -M2 -v snps -S "$target_list")
    [[ -n $region ]] && view_args+=(-r "$region")
    view_args+=("$target_src")
    bcftools "${view_args[@]}" -Ou | \
        bcftools norm --threads "$THREADS" -f "$REF_FASTA" -c s -d all -Ob -o "$work/target.norm.bcf"
    bcftools index -f "$work/target.norm.bcf"
    tagged=$work/target.tagged.bcf
    bcftools +fill-tags "$work/target.norm.bcf" --threads "$THREADS" -Ob -o "$tagged" -- -t AC,AN,MAF
    bcftools index -f "$tagged"
    bcftools view --threads "$THREADS" -i "INFO/MAF>=${MAF_MIN}" -Oz -o "$work/target.sites.vcf.gz" "$tagged"
    tabix -f -p vcf "$work/target.sites.vcf.gz"
    valid_vcf "$work/target.sites.vcf.gz" || { echo "ERROR: no valid target sites retained for $unit" >&2; exit 1; }
    mv -f "$work/target.sites.vcf.gz" "$out"
    mv -f "$work/target.sites.vcf.gz.tbi" "$out.tbi"
    input_signature "$unit" "$target_src" "$target_list" > "$sig"
    echo "[AS3 master chr${unit}] DONE variants=$(bcftools index -n "$out") samples=$(bcftools query -l "$out" | wc -l)"
)

status=0
running=0
for unit in "${units[@]}"; do
    prepare_master "$unit" &
    running=$((running + 1))
    if (( running >= MAX_PROCS )); then wait -n || status=1; running=$((running - 1)); fi
done
while (( running > 0 )); do wait -n || status=1; running=$((running - 1)); done
(( status == 0 )) || { echo "FATAL ERROR: master target preparation failed" >&2; exit 1; }

prepare_external_panel() {
    local source=$1 sample_list=$2 target_master=$3 region=$4 target_contig=$5 out=$6 mode=$7 log_file=$8
    local source_contig rename_map
    source_contig=$(first_contig "$source")
    [[ -n $source_contig ]] || { echo "ERROR: cannot determine contig from $source" >&2; return 1; }
    rename_map=$(mktemp "$TMP_ROOT/as3.rename.XXXXXX")
    printf '%s\t%s\n' "$source_contig" "$target_contig" > "$rename_map"
    local -a view_args=(view --threads "$THREADS" -m2 -M2 -v snps)
    [[ -n $sample_list ]] && view_args+=(-S "$sample_list")
    [[ -n $region ]] && view_args+=(-r "$region")
    view_args+=("$source")
    local -a fix=(--chrom "$target_contig")
    case "$mode" in
        modern_diploid) fix+=(--require-diploid --require-phased) ;;
        modern_haploid) fix+=(--duplicate-haploid pipe --require-phased) ;;
        archaic) fix+=(--duplicate-haploid pipe) ;;
        *) rm -f "$rename_map"; echo "ERROR: unknown panel mode $mode" >&2; return 1 ;;
    esac
    bcftools "${view_args[@]}" -Ou | \
        bcftools norm --threads "$THREADS" -f "$REF_FASTA" -c s -d all -Ou | \
        bcftools annotate --rename-chrs "$rename_map" -Ou | \
        bcftools view --threads "$THREADS" -T "$target_master" -Ou | \
        bcftools annotate --threads "$THREADS" -x INFO,^FORMAT/GT -Ov | \
        python3 "$script_dir/vcf_gt_fix.py" "${fix[@]}" 2> "$log_file" | \
        bgzip -@ "$THREADS" -c > "$out"
    rm -f "$rename_map"
    tabix -f -p vcf "$out"
    valid_vcf "$out"
}

process_task() (
    set -euo pipefail
    local panel=$1 unit=$2 target_list=$3 afr_list=$4 chr key target_master target_src afr_src den_src nea_src
    local target_contig region target_region afr_region den_region nea_region work target_out ref_out map_out qc_out sig_out
    chr=$(source_chr "$unit")
    key=${panel}.chr${unit}
    target_master=$OUT/Master_Target_Sites/target_sites.chr${unit}.vcf.gz
    target_src=$(resolve_vcf target "$TARGET_DIR" "$TARGET_TEMPLATE" "$chr")
    afr_src=$(resolve_vcf African-reference "$AFR_DIR" "$AFR_TEMPLATE" "$chr")
    den_src=$(resolve_vcf Denisovan "$DEN_DIR" "$DEN_TEMPLATE" "$chr")
    nea_src=$(resolve_vcf Neanderthal "$NEA_DIR" "$NEA_TEMPLATE" "$chr")
    target_out=$OUT/Final_Target_VCFs/target_panel.${key}.vcf.gz
    ref_out=$OUT/Final_Ref_VCFs/ref_panel.${key}.vcf.gz
    map_out=$OUT/Reference_Maps/reference.${key}.map
    qc_out=$OUT/qc/${key}.qc.tsv
    sig_out=$OUT/qc/${key}.input.tsv
    {
        printf 'panel\t%s\nunit\t%s\nsource_chr\t%s\nregion\t%s\ngenome_build\t%s\ntarget_master\t%s\ntarget_list_md5\t%s\nafr_list_md5\t%s\ntarget_source\t%s\nafr_source\t%s\nden_source\t%s\nnea_source\t%s\nref_fasta\t%s\nmask_prefixes\t%s\n' \
            "$panel" "$unit" "$chr" "$(unit_region "$unit")" "$GENOME_BUILD" "$(file_stamp "$target_master")" \
            "$(md5sum "$target_list" | awk '{print $1}')" "$(md5sum "$afr_list" | awk '{print $1}')" \
            "$(file_stamp "$target_src")" "$(file_stamp "$afr_src")" "$(file_stamp "$den_src")" "$(file_stamp "$nea_src")" \
            "$(file_stamp "$REF_FASTA")" "$MASK_PREFIXES_STR"
        for prefix in "${MASK_PREFIXES[@]}"; do
            [[ -n $prefix ]] || continue
            printf 'mask_%s\t%s\n' "$prefix" "$(file_stamp "$MASK_DIR/${prefix}_chr${chr}.bed")"
        done
    } > "$TMP_ROOT/${key}.signature.$$"
    if valid_vcf "$target_out" && valid_vcf "$ref_out" && [[ -s $map_out && -s $qc_out && -s $sig_out ]] && cmp -s "$TMP_ROOT/${key}.signature.$$" "$sig_out"; then
        rm -f "$TMP_ROOT/${key}.signature.$$"
        echo "[AS3 $key] SKIP"
        exit 0
    fi
    work=$(mktemp -d "$TMP_ROOT/as3_${key}.XXXXXX")
    trap 'rm -rf -- "$work" "$TMP_ROOT/${key}.signature.$$"' EXIT
    check_members "$target_list" "$target_master" "$key target" "$work"
    check_members "$afr_list" "$afr_src" "$key African-reference" "$work"
    comm -12 <(sort -u "$target_list") <(sort -u "$afr_list") > "$work/overlap.txt"
    [[ ! -s $work/overlap.txt ]] || { echo "ERROR: target/reference leakage in $key" >&2; exit 1; }

    target_contig=$(first_contig "$target_master")
    region=$(unit_region "$unit")
    target_region=$(adapt_region "$target_src" "$region")
    afr_region=$(adapt_region "$afr_src" "$region")
    den_region=$(adapt_region "$den_src" "$region")
    nea_region=$(adapt_region "$nea_src" "$region")

    local target_mode=modern_diploid
    [[ $unit == XNONPAR_M ]] && target_mode=modern_haploid
    # Extract target calls from the already normalized master sites and validate phase/ploidy.
    local -a target_fix=()
    [[ $target_mode == modern_haploid ]] && target_fix=(--duplicate-haploid pipe --require-phased) || target_fix=(--require-diploid --require-phased)
    bcftools view --threads "$THREADS" -S "$target_list" "$target_master" -Ou | \
        bcftools annotate --threads "$THREADS" -x INFO,^FORMAT/GT -Ov | \
        python3 "$script_dir/vcf_gt_fix.py" "${target_fix[@]}" 2> "$OUT/log/${key}.target.gt.log" | \
        bgzip -@ "$THREADS" -c > "$work/target.vcf.gz"
    tabix -f -p vcf "$work/target.vcf.gz"
    valid_vcf "$work/target.vcf.gz" || { echo "ERROR: invalid target VCF for $key" >&2; exit 1; }

    prepare_external_panel "$afr_src" "$afr_list" "$target_master" "$afr_region" "$target_contig" "$work/afr.vcf.gz" modern_diploid "$OUT/log/${key}.afr.gt.log" || { echo "ERROR: African reference preparation failed for $key" >&2; exit 1; }
    prepare_external_panel "$den_src" "" "$target_master" "$den_region" "$target_contig" "$work/den.vcf.gz" archaic "$OUT/log/${key}.den.gt.log" || { echo "ERROR: Denisovan preparation failed for $key" >&2; exit 1; }
    prepare_external_panel "$nea_src" "" "$target_master" "$nea_region" "$target_contig" "$work/nea.vcf.gz" archaic "$OUT/log/${key}.nea.gt.log" || { echo "ERROR: Neanderthal preparation failed for $key" >&2; exit 1; }

    bcftools merge --threads "$THREADS" --missing-to-ref -0 -Oz -o "$work/refs.vcf.gz" "$work/den.vcf.gz" "$work/nea.vcf.gz" "$work/afr.vcf.gz"
    tabix -f -p vcf "$work/refs.vcf.gz"
    mkdir "$work/isec"
    bcftools isec -p "$work/isec" -n=2 -c none --threads "$THREADS" -Oz "$work/target.vcf.gz" "$work/refs.vcf.gz"
    tabix -f -p vcf "$work/isec/0000.vcf.gz"
    tabix -f -p vcf "$work/isec/0001.vcf.gz"
    bcftools merge --threads "$THREADS" --force-samples -Oz -o "$work/merged.vcf.gz" "$work/isec/0000.vcf.gz" "$work/isec/0001.vcf.gz"
    bcftools view --threads "$THREADS" -i 'N_ALT>0' -Oz -o "$work/filtered.vcf.gz" "$work/merged.vcf.gz"
    tabix -f -p vcf "$work/filtered.vcf.gz"

    masked=$work/filtered.vcf.gz
    for prefix in "${MASK_PREFIXES[@]}"; do
        [[ -n $prefix ]] || continue
        mask=$MASK_DIR/${prefix}_chr${chr}.bed
        [[ -s $mask ]] || { echo "ERROR: mask missing: $mask" >&2; exit 1; }
        normalized_mask=$work/${prefix}.normalized.bed
        awk -v chr="$target_contig" 'BEGIN{OFS="\t"} NF>=3 && $1 !~ /^#/ {$1=chr; print}' "$mask" > "$normalized_mask"
        [[ -s $normalized_mask ]] || { echo "ERROR: empty normalized mask: $mask" >&2; exit 1; }
        next=$work/masked.${prefix}.vcf.gz
        bcftools view --threads "$THREADS" -T "$normalized_mask" -Oz -o "$next" "$masked"
        tabix -f -p vcf "$next"
        masked=$next
    done

    { list_ids "$afr_list"; bcftools query -l "$work/den.vcf.gz"; bcftools query -l "$work/nea.vcf.gz"; } | awk 'NF&&!seen[$1]++' > "$work/ref.list"
    { list_ids "$afr_list" | awk '{print $1"\tAFR"}'; bcftools query -l "$work/den.vcf.gz" | awk '{print $1"\tDEN"}'; bcftools query -l "$work/nea.vcf.gz" | awk '{print $1"\tNEAN"}'; } | awk 'NF==2&&!seen[$1]++' > "$work/reference.map"
    [[ -s $work/reference.map ]] || { echo "ERROR: empty reference map for $key" >&2; exit 1; }

    bcftools view --threads "$THREADS" -S "$target_list" -Oz -o "$work/target.final.vcf.gz" "$masked"
    tabix -f -p vcf "$work/target.final.vcf.gz"
    bcftools view --threads "$THREADS" -S "$work/ref.list" -Oz -o "$work/ref.final.vcf.gz" "$masked"
    tabix -f -p vcf "$work/ref.final.vcf.gz"
    valid_vcf "$work/target.final.vcf.gz" || { echo "ERROR: invalid final target VCF for $key" >&2; exit 1; }
    valid_vcf "$work/ref.final.vcf.gz" || { echo "ERROR: invalid final reference VCF for $key" >&2; exit 1; }

    # Final target and reference sample sets must remain disjoint and exact.
    bcftools query -l "$work/target.final.vcf.gz" | sort -u > "$work/final.target.samples"
    bcftools query -l "$work/ref.final.vcf.gz" | sort -u > "$work/final.ref.samples"
    comm -12 "$work/final.target.samples" "$work/final.ref.samples" > "$work/final.overlap"
    [[ ! -s $work/final.overlap ]] || { echo "ERROR: final target/reference overlap in $key" >&2; exit 1; }
    cmp -s <(sort -u "$target_list") "$work/final.target.samples" || { echo "ERROR: final target sample mismatch in $key" >&2; exit 1; }
    cmp -s <(sort -u "$work/ref.list") "$work/final.ref.samples" || { echo "ERROR: final reference sample mismatch in $key" >&2; exit 1; }

    target_part=${target_out%.vcf.gz}.part.$$.vcf.gz
    ref_part=${ref_out%.vcf.gz}.part.$$.vcf.gz
    map_part=${map_out}.part.$$
    qc_part=${qc_out}.part.$$
    rm -f "$target_part" "$target_part.tbi" "$ref_part" "$ref_part.tbi" "$map_part" "$qc_part"
    mv "$work/target.final.vcf.gz" "$target_part"; mv "$work/target.final.vcf.gz.tbi" "$target_part.tbi"
    mv "$work/ref.final.vcf.gz" "$ref_part"; mv "$work/ref.final.vcf.gz.tbi" "$ref_part.tbi"
    cp "$work/reference.map" "$map_part"
    printf 'panel_id\t%s\nunit\t%s\nsource_chr\t%s\nregion\t%s\ntarget_samples\t%s\nafr_reference_samples\t%s\nreference_samples_total\t%s\nvariants\t%s\nexperimental_chrX\t%s\n' \
        "$panel" "$unit" "$chr" "$region" "$(count_lines "$target_list")" "$(count_lines "$afr_list")" \
        "$(bcftools query -l "$ref_part" | wc -l)" "$(bcftools index -n "$target_part")" "$([[ $unit == X* ]] && echo 1 || echo 0)" > "$qc_part"
    mv -f "$target_part" "$target_out"; mv -f "$target_part.tbi" "$target_out.tbi"
    mv -f "$ref_part" "$ref_out"; mv -f "$ref_part.tbi" "$ref_out.tbi"
    mv -f "$map_part" "$map_out"; mv -f "$qc_part" "$qc_out"
    mv -f "$TMP_ROOT/${key}.signature.$$" "$sig_out"
    trap - EXIT
    rm -rf "$work"
    echo "[AS3 $key] DONE variants=$(bcftools index -n "$target_out") target_samples=$(bcftools query -l "$target_out" | wc -l)"
)

status=0
running=0
while IFS=$'\t' read -r panel unit target_list afr_list; do
    [[ $panel == panel_id ]] && continue
    process_task "$panel" "$unit" "$target_list" "$afr_list" &
    running=$((running + 1))
    if (( running >= MAX_PROCS )); then wait -n || status=1; running=$((running - 1)); fi
done < "$OUT/tasks.tsv"
while (( running > 0 )); do wait -n || status=1; running=$((running - 1)); done
(( status == 0 )) || { echo "FATAL ERROR: one or more AS3 preprocessing tasks failed" >&2; exit 1; }

printf 'panel_id\tunit\ttarget_vcf\treference_vcf\treference_map\ttarget_list\tafr_reference_list\texperimental_chrX\n' > "$OUT/manifest.tsv"
while IFS=$'\t' read -r panel unit target_list afr_list; do
    [[ $panel == panel_id ]] && continue
    key=${panel}.chr${unit}
    printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
        "$panel" "$unit" \
        "$OUT/Final_Target_VCFs/target_panel.${key}.vcf.gz" \
        "$OUT/Final_Ref_VCFs/ref_panel.${key}.vcf.gz" \
        "$OUT/Reference_Maps/reference.${key}.map" \
        "$target_list" "$afr_list" "$([[ $unit == X* ]] && echo 1 || echo 0)"
done < "$OUT/tasks.tsv" >> "$OUT/manifest.tsv"

printf 'analysis_units\t%s\noverlap_policy\t%s\noverlap_folds\t%s\ngenome_build\t%s\ntarget_maf_min\t%s\n' \
    "${units[*]}" "$OVERLAP_POLICY" "$OVERLAP_FOLDS" "$GENOME_BUILD" "$MAF_MIN" > "$OUT/preparation.summary.tsv"
echo "AS3 preprocessing complete: $OUT manifest=$OUT/manifest.tsv"
