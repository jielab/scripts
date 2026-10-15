#!/usr/bin/env bash
set -euo pipefail
export PYTHONDONTWRITEBYTECODE=1
ROOT=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
F="$ROOT/f"
# shellcheck source=f/common.sh
source "$F/common.sh"
allele_load_config "$ROOT"

usage() {
cat <<'HELP'
用法：./allele.sh COMMAND [options]

命令
  check-input 检查待分析 VCF 是否齐全、下载是否完成（不扫描全文件）
  inspect    检查环境、依赖和数据路径
  test       离线回归测试，不下载数据、不调用 API
  prepare    从群体 VCF 的频率/计数和 GRCh38 FASTA 构建候选变异
  select     选择多等位位点及匹配的双等位对照
  atlas      获取 AlphaGenome 效应向量，需要 API key
  geometry   计算效应向量的几何特征
  stage1     拟合突变机会模型
  stage2     拟合保留关联模型
  final      汇总图表
  all        从 prepare 到 final 顺序运行，需要数据和 API key

正式运行示例（gnomAD v4.1 genomes；按顺序逐条复制，某一步报错先停止）
  以下命令在 /mnt/d/scripts/allele 目录运行。
  当前机器已安装环境，不需重跑 install.sh，也不需 conda activate。

  # 1. 等下载完成后检查 chr1–22 的 VCF 和索引；成功后才继续
  ./allele.sh check-input --vcf-dir /mnt/f/annot/gnomad --chr autosomes

  # 2. 构建候选表；完整读取 chr1–22，不限制区间。首次使用一个新输出目录
  ./allele.sh prepare --vcf-dir /mnt/f/annot/gnomad --chr autosomes \
    --fasta /mnt/f/gen/1hgp/38/GRCH38.fasta.gz --out /mnt/d/analysis/allele_gnomad_v4.1

  # 3. 选择全部多等位位点和匹配对照；打印需要查询的 Atlas 规模
  ./allele.sh select --out /mnt/d/analysis/allele_gnomad_v4.1

  # 4. 突变机会模型；这一步不需要 API key
  ./allele.sh stage1 --out /mnt/d/analysis/allele_gnomad_v4.1

  # 5. 自动从 /mnt/d/files/alphaG.api 读取 key；不需要再输入或 export
  ./allele.sh atlas --out /mnt/d/analysis/allele_gnomad_v4.1

  # 6. Atlas 完成后计算几何、保留关联并出图
  ./allele.sh geometry --out /mnt/d/analysis/allele_gnomad_v4.1
  ./allele.sh stage2 --out /mnt/d/analysis/allele_gnomad_v4.1
  ./allele.sh final --out /mnt/d/analysis/allele_gnomad_v4.1

  # 结果目录：/mnt/d/analysis/allele_gnomad_v4.1/final/
  # atlas 中断后，复制同一条 atlas 命令即可续跑

可选：1KG 小区间测试（与正式 gnomAD 输出分开，不是上面步骤的前置要求）
  ./allele.sh prepare --vcf /mnt/f/gen/1kg/38/vcf/chr22.vcf.gz \
    --region chr22:20000001-21000000 --out /mnt/d/analysis/allele_1kg_try
  ./allele.sh select --out /mnt/d/analysis/allele_1kg_try
  ./allele.sh stage1 --out /mnt/d/analysis/allele_1kg_try
  ./allele.sh test

常用选项
  --vcf FILE / --vcf-dir DIR            单个 VCF / 目录内所有 VCF，二选一
  --region CHR:START-END                只读取指定区间（1-based，需 VCF 索引）
  --chr autosomes / all / 1,2,3         目录输入的染色体；gnomAD 默认 1–22
  --fasta FILE --out DIR                覆盖参考序列/输出路径
  --biallelic-ratio N                   双等位对照比例，默认 2
  --scorers A,B,... --chunk N           Atlas 模态及批量大小，默认 128

数据和 key 已齐备时，可用 all 加相同输入参数从头运行。
默认正式输入匹配 gnomad.genomes.v4.1.sites.chr*.vcf.bgz；检查失败会停止。
X/Y 不在常染色体分析中，显式 --chr all 才加入；需另行考虑倍性。
prepare 已有结果时不会覆盖，请使用新的 --out。后续步骤必须使用同一个 --out。
全基因组 Atlas 的时间和空间取决于 select 的位点数，离线测试不验证 API 授权。
更多说明见 README.md；自定义路径可复制 config.example.env 为 allele.env。
HELP
printf '\n当前路径\n  VCF 目录: %s\n  FASTA:    %s\n  输出:     %s\n' \
  "$ALLELE_VCF_DIR" "$ALLELE_FASTA" "$ALLELE_OUT"
}

OUT=$ALLELE_OUT
VCF_DIR=$ALLELE_VCF_DIR
VCF_FILE=
VCF_DIR_SET=
REGION=
CHROMS=
FASTA=$ALLELE_FASTA
SCORERS=$ALLELE_ATLAS_SCORERS
CHUNK=$ALLELE_ATLAS_CHUNK
BIRATIO=$ALLELE_BIALLELIC_RATIO

METHOD=${1:-help}; shift || true
case "$METHOD" in
  help|-h|--help) usage; exit 0;;
  inspect|check-input|prepare|select|atlas|geometry|stage1|stage2|final|all|test) ;;
  *) allele_die "unknown command: $METHOD";;
esac
while (($#)); do
  case "$1" in
    --out|--vcf|--vcf-dir|--region|--chr|--fasta|--scorers|--chunk|--biallelic-ratio)
      [[ $# -ge 2 && -n "$2" && "$2" != --* ]] || allele_die "missing value for $1";;
  esac
  case "$1" in
    --out) OUT=$2; shift 2;;
    --vcf) VCF_FILE=$2; shift 2;;
    --vcf-dir) VCF_DIR=$2; VCF_DIR_SET=1; shift 2;;
    --region) REGION=$2; shift 2;;
    --chr) CHROMS=$2; shift 2;;
    --fasta) FASTA=$2; shift 2;;
    --scorers) SCORERS=$2; shift 2;;
    --chunk) CHUNK=$2; shift 2;;
    --biallelic-ratio) BIRATIO=$2; shift 2;;
    -h|--help) usage; exit 0;;
    *) allele_die "unknown argument: $1";;
  esac
done
[[ "$CHUNK" =~ ^[1-9][0-9]*$ ]] || allele_die '--chunk must be a positive integer'
[[ "$BIRATIO" =~ ^[0-9]+([.][0-9]+)?$ ]] || allele_die '--biallelic-ratio must be nonnegative'
[[ -z "$VCF_FILE" || -z "$VCF_DIR_SET" ]] || allele_die 'use --vcf or --vcf-dir, not both'
[[ -z "$REGION" || -n "$VCF_FILE" ]] || allele_die '--region requires --vcf'
[[ -z "$CHROMS" || -z "$VCF_FILE" ]] || allele_die '--chr is for directory input; use --region with --vcf'
export ALLELE_OUT="$OUT" ALLELE_VCF_DIR="$VCF_DIR" ALLELE_FASTA="$FASTA"
allele_require_env
if [[ "$METHOD" != inspect && "$METHOD" != test && "$METHOD" != check-input ]]; then
  mkdir -p "$OUT"/{work,atlas,results,final,logs}
fi

run_prepare() {
  allele_require_file "$FASTA"
  local -a input=(--vcf-dir "$VCF_DIR")
  if [[ -n "$VCF_FILE" ]]; then
    allele_require_file "$VCF_FILE"
    input=(--vcf "$VCF_FILE")
  fi
  [[ -z "$REGION" ]] || input+=(--region "$REGION")
  [[ -z "$CHROMS" ]] || input+=(--chr "$CHROMS")
  [[ "$METHOD" != check-input ]] || input+=(--check-only)
  allele_note "prepare: VCF=${VCF_FILE:-$VCF_DIR} REGION=${REGION:-all} FASTA=$FASTA"
  allele_python "$F/extract_sites.py" "${input[@]}" --fasta "$FASTA" \
    --out "$OUT/work/candidates.parquet" --sites-out "$OUT/work/sites.parquet"
}
run_select() {
  allele_require_file "$OUT/work/sites.parquet"
  allele_python "$F/select_atlas_sites.py" --sites "$OUT/work/sites.parquet" \
    --ratio "$BIRATIO" --seed "$ALLELE_SEED" --out "$OUT/work/atlas_sites.parquet"
}
run_atlas() {
  allele_require_file "$OUT/work/atlas_sites.parquet"
  [[ -n ${ALPHA_GENOME_API_KEY:-} ]] || allele_die "API key missing: save one plain-text line to $ALLELE_API_KEY_FILE, or set ALPHA_GENOME_API_KEY"
  allele_python "$F/atlas_query.py" --sites "$OUT/work/atlas_sites.parquet" \
    --out-dir "$OUT/atlas" --scorers "$SCORERS" --chunk "$CHUNK"
}
run_geometry() {
  allele_require_file "$OUT/atlas/rows.parquet"
  allele_python "$F/geometry.py" --candidates "$OUT/work/candidates.parquet" \
    --atlas-sites "$OUT/work/atlas_sites.parquet" --atlas-dir "$OUT/atlas" \
    --out-dir "$OUT/results"
}
run_stage1() {
  allele_require_file "$OUT/work/candidates.parquet"
  allele_python "$F/stage1_mutation.py" --candidates "$OUT/work/candidates.parquet" \
    --out-dir "$OUT/results" --seed "$ALLELE_SEED"
}
run_stage2() {
  allele_require_file "$OUT/results/geometry_candidates.parquet"
  allele_require_file "$OUT/results/stage1_candidates.parquet"
  allele_python "$F/stage2_retention.py" --geometry "$OUT/results/geometry_candidates.parquet" \
    --sites "$OUT/results/geometry_sites.parquet" --mutation "$OUT/results/stage1_candidates.parquet" \
    --out-dir "$OUT/results" --seed "$ALLELE_SEED"
}
run_final() {
  allele_python "$F/final.py" --work "$OUT/work" --results "$OUT/results" --out-dir "$OUT/final"
}

case "$METHOD" in
  help|-h|--help) usage;;
  inspect) allele_python "$F/validate.py" --root "$ROOT" --out "$OUT";;
  test) allele_python -m unittest discover -s "$ROOT/tests" -v;;
  check-input) run_prepare;;
  prepare) run_prepare;;
  select) run_select;;
  atlas) run_atlas;;
  geometry) run_geometry;;
  stage1) run_stage1;;
  stage2) run_stage2;;
  final) run_final;;
  all) run_prepare; run_select; run_stage1; run_atlas; run_geometry; run_stage2; run_final;;
  *) allele_die "unknown command: $METHOD";;
esac
