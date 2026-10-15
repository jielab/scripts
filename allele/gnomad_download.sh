#!/usr/bin/env bash
# Standalone gnomAD v4.1 genome sites downloader. Requires Bash and aria2c.
set -euo pipefail

usage() {
  cat <<'HELP'
用法：./gnomad_download.sh --chr 1,2,3,4,5 [options]

  --chr LIST      必填：1..22,X,Y 的逗号列表，或 all
  --proxy URL     HTTP 代理，默认 http://127.0.0.1:7897；none 表示直连
  --out DIR       保存目录，默认 /mnt/f/annot/gnomad/4.1/genomes
  --jobs N        同时下载的文件数，默认 2
  -h, --help      显示帮助

多台电脑分别运行（各自只需本脚本和 aria2c）：
  ./gnomad_download.sh --proxy http://127.0.0.1:7897 --chr 1,2,3,4,5
  ./gnomad_download.sh --proxy http://127.0.0.1:7897 --chr 6,7,8,9,10 --out /data/gnomad

单个染色体续传 / 全量下载：
  ./gnomad_download.sh --chr 22
  ./gnomad_download.sh --chr all --proxy none

下载 gnomAD v4.1 genomes sites summary 和 .tbi 索引，无个体 genotype。
内置官方 MD5，使用 Google Cloud/AWS 镜像；每文件最多 16 连接。
Ctrl+C 后重跑同一命令即可续传；保留旁边的 .aria2 控制文件。
不同电脑分配不同染色体；完成后汇总 .vcf.bgz 和 .tbi 文件。
不读取 allele.env，不需要 conda；环境变量 GNOMAD_OUT/GNOMAD_PROXY 可设置默认值。
HELP
}

die() { printf 'ERROR: %s\n' "$*" >&2; exit 2; }
need_value() { [[ -n "${2:-}" && "$2" != --* ]] || die "missing value for $1"; }

chromosomes=
proxy=${GNOMAD_PROXY:-http://127.0.0.1:7897}
out=${GNOMAD_OUT:-/mnt/f/annot/gnomad/4.1/genomes}
jobs=2
if (($# == 0)); then usage; exit 0; fi
while (($#)); do
  case "$1" in
    --chr) need_value "$@"; chromosomes=$2; shift 2;;
    --proxy) need_value "$@"; proxy=$2; shift 2;;
    --out) need_value "$@"; out=$2; shift 2;;
    --jobs) need_value "$@"; jobs=$2; shift 2;;
    -h|--help) usage; exit 0;;
    *) die "unknown argument: $1";;
  esac
done
[[ -n "$chromosomes" ]] || die '--chr is required; use e.g. --chr 1,2,3 or --chr all'
if [[ "$chromosomes" == all ]]; then
  chrs=({1..22} X Y)
  chromosomes=$(IFS=,; printf '%s' "${chrs[*]}")
fi
[[ "$chromosomes" != *, && "$chromosomes" != ,* && "$chromosomes" != *,,* ]] || die 'invalid --chr list'
IFS=',' read -r -a chrs <<< "$chromosomes"
for chrom in "${chrs[@]}"; do
  [[ "$chrom" =~ ^([1-9]|1[0-9]|2[0-2]|X|Y)$ ]] || die "invalid chromosome: $chrom (use 1..22,X,Y)"
done
[[ "$jobs" =~ ^[1-9][0-9]*$ ]] || die '--jobs must be a positive integer'

if [[ "$proxy" == none ]]; then
  proxy_args=(--http-proxy= --https-proxy= --all-proxy= --no-proxy=)
else
  [[ "$proxy" =~ ^http://([A-Za-z0-9._-]+|\[[0-9A-Fa-f:]+\]):([0-9]{1,5})/?$ ]] || \
    die '--proxy must be http://HOST:PORT (without credentials), or none'
  proxy_port=$((10#${BASH_REMATCH[2]}))
  ((proxy_port >= 1 && proxy_port <= 65535)) || die '--proxy port must be 1..65535'
  proxy_args=(--http-proxy="$proxy" --https-proxy="$proxy" --all-proxy="$proxy" --no-proxy=)
fi
command -v aria2c >/dev/null 2>&1 || die 'aria2c not found; install aria2 (Ubuntu/WSL: sudo apt install aria2)'

# Official object MD5 hashes retrieved from the public GCS JSON API, 2026-10-15.
# Columns: chromosome, VCF MD5, tabix index MD5. Keep this pinned to release 4.1.
catalog() {
  cat <<'DATA'
22 a7bcf712a6b8d29e690468bf1dd8913d 44f5ef4f99bfed0ec833686c8ba3abb0
1 328b4578212afec2cde394a1b02d544f 02c67bb77109f1638650c0cfa389efd3
2 518ca01e6757a68bc0abe76f85af644d 495c2e6b578719dc77c6fdca2e799875
3 716c181431a3c11a2eb18c5f50a3542d 7c1561073864f92ce76682f5a8049d17
4 d9b913f3e30c8f410f9ce7dee5dee6d4 2ae3b1436d71057eb513a233c73b6278
5 a2a3b9014af5c8f9bbaaf743968e48f8 4cb33aebe76d6daecf97f61e52da80dc
6 e65c2aa321c5e272548fb3d901bae382 59132364c566098d6c9640fda7374b9c
7 58ee22cf3dcc8cb8b493d218e19432cc 7a73677be3fe16da3bdd051602c1f3de
8 9854d9df22977cf5bac0f9fc05f4e8f5 dd27f62d7160200e4b30b939f45b69c0
9 6adfc9c47000cf66d1305392051b391d c8db53ed8db04a0b5a95ff485e602e14
10 c2cd760130d2339f7135fc70700db1e1 8e1b974d66d41a2105c8aeb6c2b30bab
11 d1e7a4dcf3ff62eeffca57afefc5a33a cfcc3e35ebba2aba1fa01d5511f1677c
12 644bdbc5c53d9112edbacad401a28d1a bfac8b89bfebb82c155875ca28bff2e8
13 84b12f299210d2a7e390c56a234b7b68 fba3906439434b68c0e48c14d8a2d548
14 0e524b414faf5a51b74d153939b3ddcd 12f370988f353e643133c516bd7f00c4
15 00d1386eadbfcb653eb810a6c08ed250 c7c41dcc9e053954af7562ea0200a040
16 2b106c12ca8cca9bd58e98d2c248ef4d 94c1a5a500f088eef4670ebff98b0f41
17 9d58b459e75312b487c660666ea540c4 1cd1d1ba0e320cb609f033b7b6e5da47
18 b95d30a6f5a45242d834fc7351afd760 565b78fa172f1e20617ec9c66bdabc45
19 73fa1c7c09072e8ca5e26a1443f0af2a b8af5ad748507cc4291165f705b4056b
20 6c6f326c67b288ec99932905da72c1e6 c6d3ce76c55d9bd41d849b7dd5ccad05
21 9b6584cfe62f6e8bd0c3cea90a4ce56a 83cf6b974c7e13c12279c0f195d4ce48
X 8b91766906865b0795c653af51cb73b8 4895e179e13fe52ff411e52b7fd695ca
Y 1ffb9c683674f41ff7cf524e5bb56bb8 3e156c7a8d7c3c8c4e9e14cacd7e6556
DATA
}

gcp=https://storage.googleapis.com/gcp-public-data--gnomad/release/4.1/vcf/genomes
aws=https://gnomad-public-us-east-1.s3.amazonaws.com/release/4.1/vcf/genomes
download_manifest=$(mktemp)
trap 'rm -f -- "$download_manifest"' EXIT
while read -r chrom vcf_md5 index_md5; do
  [[ ",$chromosomes," == *",$chrom,"* ]] || continue
  file="gnomad.genomes.v4.1.sites.chr${chrom}.vcf.bgz"
  printf '%s/%s\t%s/%s\n  out=%s\n  checksum=md5=%s\n' \
    "$gcp" "$file" "$aws" "$file" "$file" "$vcf_md5"
  file+=.tbi
  printf '%s/%s\t%s/%s\n  out=%s\n  checksum=md5=%s\n' \
    "$gcp" "$file" "$aws" "$file" "$file" "$index_md5"
done < <(catalog) > "$download_manifest"
[[ -s "$download_manifest" ]] || die 'no files matched the requested chromosomes'
mkdir -p -- "$out"
printf '[gnomAD v4.1] chromosomes=%s\nOutput: %s\nProxy: %s\n' "$chromosomes" "$out" "$proxy" >&2
aria2c --no-conf=true -c "--max-concurrent-downloads=$jobs" -s16 -x8 -k4M \
  "${proxy_args[@]}" \
  --uri-selector=adaptive --file-allocation=none \
  --auto-file-renaming=false --allow-overwrite=false --always-resume=true \
  --check-integrity=true --connect-timeout=30 --timeout=60 \
  --max-tries=20 --retry-wait=10 --auto-save-interval=30 \
  --dir="$out" --input-file="$download_manifest"
