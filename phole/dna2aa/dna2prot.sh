#!/usr/bin/env bash

# 从 DNA + 注释模型，推出某些候选转录本对应的预测蛋白序列。
# 输出应该理解为 transcript-specific predicted proteins，
# 而不是“该个体真实体内唯一存在的蛋白序列”。

set -euo pipefail

dir0=/mnt/d
diref=$dir0/data/annot/ensembl
dirgen=$dir0/data/gen/1kg/high_coverage
dirwrk=$dir0/analysis/vcf2prot
mkdir -p "$dirwrk"


# bcftools query -l $dirgen/chr1.vcf.gz > $dir0/files/1kg.sam.txt
awk '$4=="ABO"' "$dir0/files/ppp.b38.bed" | while IFS=$'\t' read -r chr pos1 pos2 gene; do
    echo "=============================="
    echo "Processing gene: $gene  region: $chr:$pos1-$pos2"
    echo "=============================="

    outdir=$dirwrk/$gene
    mkdir -p "$outdir"/{vcf,dna,aa,alphaG,logs}

    loc=${chr}:${pos1}-${pos2}
    vcf_chr=$dirgen/chr${chr}.vcf.gz
    fa_chr=$diref/fa/Homo_sapiens.GRCh38.dna.chromosome.${chr}.fa
    gff3_chr=$diref/gff3/Homo_sapiens.GRCh38.114.chromosome.${chr}.gff3.gz
    gtf_all=$diref/Homo_sapiens.GRCh38.114.gtf.gz

    # 1) 提取基因区域 VCF
    # 输入 VCF 染色体名是 1,2,3...，无 chr 前缀，所以这里直接用 $loc
    bcftools view "$vcf_chr"  --regions "$loc" -Oz -o "$outdir/vcf/${gene}.vcf.gz"
    tabix -f "$outdir/vcf/${gene}.vcf.gz"

    # 2) 为区域 VCF 添加 consequence 注释，适合 phased VCF
    bcftools csq "$outdir/vcf/${gene}.vcf.gz" -f "$fa_chr" -g "$gff3_chr" -Oz -o "$outdir/vcf/${gene}.csq.vcf.gz"
    tabix -f "$outdir/vcf/${gene}.csq.vcf.gz"

    # 3) 生成每个样本两条 haplotype DNA FASTA
    # GT是phased，可直接用 -H 1 和 -H 2。若将来用phased/unphased 混合数据，可考虑 1pIu / 2pIu
    parallel -j16 --bar \
      "samtools faidx '$fa_chr' '$loc' | bcftools consensus -s {1} -H {2} '$outdir/vcf/${gene}.vcf.gz' > '$outdir/dna/{1}.H{2}.fa'" \
      :::: "$dir0/files/1kg.sam.txt" ::: 1 2

    # 4) 用你自己的 transcript 翻译脚本，从 haplotype DNA FASTA -> AA FASTA
    # 这一步不用再依赖 vcf2prot examples/reference_sequences.fasta，而是统一使用当前同版 GTF/FASTA 注释。
    parallel -j8 --bar \
      "python $dir0/scripts/main/dna2prot.py \
        --gene '$gene' --chr '$chr' --start '$pos1' --end '$pos2' \
        --hap-fa '$outdir/dna/{1}.H{2}.fa' --gtf '$gtf_all' \
        --out-fa '$outdir/aa/{1}.H{2}.aa.fa' \
        > '$outdir/logs/{1}.H{2}.dna2prot.log' 2>&1" \
      :::: "$dir0/files/1kg.sam.txt" ::: 1 2

    # 5) 用 AlphaGenome 🚩 做 RNA/splicing 辅助预测
    # AlphaGenome 接口/教程通常使用 chr 前缀坐标 chr${chr}
    parallel -j4 --bar \
      "python $dir0/scripts/main/alphaG_batch.py \
        --gene '$gene' --chr 'chr{3}' --start '$pos1' --end '$pos2' \
        --hap-fa '$outdir/dna/{1}.H{2}.fa' --gtf '$gtf_all' \
        --outdir '$outdir/alphaG/{1}.H{2}' \
        > '$outdir/logs/{1}.H{2}.alphaG.log' 2>&1" \
      :::: "$dir0/files/1kg.sam.txt" ::: 1 2 ::: "$chr"

done