#!/usr/bin/env bash
# HF × diet fizzy drinks liking：eQTLGen / GTEx Whole Blood cis-SMR。
set -euo pipefail

# 1. Paths
PROJECT="/mnt/d/projects/flt"
OUT_PATH="$PROJECT/04flt_smr/G732-mdd_2025"
OUTCOME_FILE="/mnt/d/data/gwas/brain/common/mdd_2025/gwas/mdd_2025.ma"
FLT_FILE="/mnt/d/data/gwas/flt/common/G732/gwas/G732.ma"

SMR_PATH="/mnt/d/software/linux/smr/smr-1.3.1-linux-x86_64/smr"
BFILE_PATH="/mnt/d/data/ldref/opengwas/1kg.v3/EUR"
EQTLGEN_PATH="/mnt/d/data/gwas/eqtl/eqtlgen_v1/cis-eQTLs-full_eQTLGen_AF_incl_nr_formatted_20191212.new.txt_besd-dense"
GTEX_PATH="/mnt/d/data/gwas/eqtl/gtex_v8_whole_blood/Whole_Blood/Whole_Blood"
SCRIPT_PATH="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
THREADS=10

# 2. Check inputs（bash s7_2_flt_smr.sh --check 仅检查，不运行分析）
[[ $# -eq 0 || ( $# -eq 1 && "$1" == "--check" ) ]] || { echo "用法：bash s7_2_flt_smr.sh [--check]" >&2; exit 1; }
[[ -x "$SMR_PATH" ]] || { echo "找不到 SMR：$SMR_PATH" >&2; exit 1; }
command -v python3 >/dev/null
for FILE in "$OUTCOME_FILE" "$FLT_FILE" "$SCRIPT_PATH/s7_2_collect_smr.py"     "$BFILE_PATH".{bed,bim,fam} "$EQTLGEN_PATH".{besd,epi,esi} "$GTEX_PATH".{besd,epi,esi}; do
    [[ -r "$FILE" && -s "$FILE" ]] || { echo "缺少输入：$FILE" >&2; exit 1; }
done
for FILE in "$OUTCOME_FILE" "$FLT_FILE"; do
    awk 'NR==1 {sub(/\r$/, ""); exit !(tolower($0) ~ /^snp[ \t]+a1[ \t]+a2[ \t]+freq[ \t]+b[ \t]+se[ \t]+p[ \t]+n$/)}' "$FILE"         || { echo "GWAS 表头应为 SNP A1 A2 freq b se p n：$FILE" >&2; exit 1; }
done
printf 'HF: %s\nFLT: %s\n输入检查通过。\n' "$OUTCOME_FILE" "$FLT_FILE"
[[ "${1:-}" != "--check" ]] || exit 0

# 每次运行独立目录，避免覆盖已有结果。
mkdir -p "$OUT_PATH"
OUT_PATH=$(mktemp -d "$OUT_PATH/run_$(date +%Y%m%d_%H%M%S)_XXXXXX")
echo "本次输出：$OUT_PATH"

# 3. Cis-SMR：保留原脚本参数；不额外添加 trans 分析。
# diff-freq=0.20、diff-freq-prop=0.05，未照搬参考脚本的 0.999。

## eQTLGen: hf
BEQTL_PATH="$EQTLGEN_PATH"
FILE="$OUTCOME_FILE"
OUTFILE="$OUT_PATH/eqtlgen_v1__hf"
echo "SMR: eqtlgen_v1 / hf"
"$SMR_PATH" --bfile "$BFILE_PATH" --gwas-summary "$FILE" --beqtl-summary "$BEQTL_PATH" \
    --out "$OUTFILE" --cis-wind 2000 --peqtl-smr 5e-8 --heidi-max-m 20 \
    --diff-freq 0.20 --diff-freq-prop 0.05 --thread-num "$THREADS" > "$OUTFILE.log" 2>&1
[[ -s "$OUTFILE.smr" ]] || { echo "没有生成结果，请查看 $OUTFILE.log" >&2; exit 1; }

## eQTLGen: diet_fizzy_drinks
BEQTL_PATH="$EQTLGEN_PATH"
FILE="$FLT_FILE"
OUTFILE="$OUT_PATH/eqtlgen_v1__diet_fizzy_drinks"
echo "SMR: eqtlgen_v1 / diet_fizzy_drinks"
"$SMR_PATH" --bfile "$BFILE_PATH" --gwas-summary "$FILE" --beqtl-summary "$BEQTL_PATH" \
    --out "$OUTFILE" --cis-wind 2000 --peqtl-smr 5e-8 --heidi-max-m 20 \
    --diff-freq 0.20 --diff-freq-prop 0.05 --thread-num "$THREADS" > "$OUTFILE.log" 2>&1
[[ -s "$OUTFILE.smr" ]] || { echo "没有生成结果，请查看 $OUTFILE.log" >&2; exit 1; }

## GTEx v8 Whole Blood: hf
BEQTL_PATH="$GTEX_PATH"
FILE="$OUTCOME_FILE"
OUTFILE="$OUT_PATH/gtex_v8_whole_blood__hf"
echo "SMR: gtex_v8_whole_blood / hf"
"$SMR_PATH" --bfile "$BFILE_PATH" --gwas-summary "$FILE" --beqtl-summary "$BEQTL_PATH" \
    --out "$OUTFILE" --cis-wind 2000 --peqtl-smr 5e-8 --heidi-max-m 20 \
    --diff-freq 0.20 --diff-freq-prop 0.05 --thread-num "$THREADS" > "$OUTFILE.log" 2>&1
[[ -s "$OUTFILE.smr" ]] || { echo "没有生成结果，请查看 $OUTFILE.log" >&2; exit 1; }

## GTEx v8 Whole Blood: diet_fizzy_drinks
BEQTL_PATH="$GTEX_PATH"
FILE="$FLT_FILE"
OUTFILE="$OUT_PATH/gtex_v8_whole_blood__diet_fizzy_drinks"
echo "SMR: gtex_v8_whole_blood / diet_fizzy_drinks"
"$SMR_PATH" --bfile "$BFILE_PATH" --gwas-summary "$FILE" --beqtl-summary "$BEQTL_PATH" \
    --out "$OUTFILE" --cis-wind 2000 --peqtl-smr 5e-8 --heidi-max-m 20 \
    --diff-freq 0.20 --diff-freq-prop 0.05 --thread-num "$THREADS" > "$OUTFILE.log" 2>&1
[[ -s "$OUTFILE.smr" ]] || { echo "没有生成结果，请查看 $OUTFILE.log" >&2; exit 1; }

# 4. 两边均满足 P_SMR<0.05、P_HEIDI>0.05 后，按 probeID 取交集。
python3 "$SCRIPT_PATH/s7_2_collect_smr.py" "$OUT_PATH"
