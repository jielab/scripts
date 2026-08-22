#!/usr/bin/env bash


# 🚩 文件保存在/mnt/d/data.BIG/annot

# GENCODE
wget -c https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_19/gencode.v19.annotation.gtf.gz
wget -c https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_50/gencode.v50.basic.annotation.gtf.gz

zcat gencode.v50.basic.annotation.gtf.gz \
  | awk 'BEGIN{FS=OFS="\t"} $3=="gene" {
      gene_id=""; gene_name=""; gene_type="";
      if (match($9, /gene_id "([^"]+)"/, a)) gene_id=a[1];
      if (match($9, /gene_name "([^"]+)"/, b)) gene_name=b[1];
      if (match($9, /gene_type "([^"]+)"/, c)) gene_type=c[1];
      print $1, $4, $5, gene_id, gene_name, gene_type, $7
    }' \
  > gencode.v50.GRCh38.genes.tsv

zcat gencode.v19.annotation.gtf.gz \
  | perl -F'\t' -lane '
      next if /^#/;
      next unless $F[2] eq "gene";
      ($gene_id)   = $F[8] =~ /gene_id "([^"]+)"/;
      ($gene_name) = $F[8] =~ /gene_name "([^"]+)"/;
      ($gene_type) = $F[8] =~ /gene_type "([^"]+)"/;
      print join("\t", $F[0], $F[3], $F[4], $gene_id, $gene_name, $gene_type, $F[6]);
    ' \
  > gencode.v19.GRCh37.genes.tsv
  

# Ensembl
wget -c https://ftp.ensembl.org/pub/release-116/gtf/homo_sapiens/Homo_sapiens.GRCh38.116.gtf.gz
wget -c https://ftp.ensembl.org/pub/grch37/release-87/gtf/homo_sapiens/Homo_sapiens.GRCh37.87.gtf.gz

zcat Homo_sapiens.GRCh38.116.gtf.gz \
  | perl -F'\t' -lane '
      next if /^#/;
      next unless $F[2] eq "gene";
      ($gene_id)   = $F[8] =~ /gene_id "([^"]+)"/;
      ($gene_name) = $F[8] =~ /gene_name "([^"]+)"/;
      ($biotype)   = $F[8] =~ /gene_biotype "([^"]+)"/;
      print join("\t", $F[0], $F[3], $F[4], $gene_id, $gene_name, $biotype, $F[6]);
    ' \
  > ensembl.GRCh38.genes.tsv