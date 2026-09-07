#!/usr/bin/env bash
set -euo pipefail
here=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
usage() {
  cat <<'EOF'
GWAS comparison and LDSC modules:
  ./gwas_compare.sh compare [options]
  ./gwas_compare.sh ldsc [options]

compare options:
  --gwas-files CSV       Two or more standardized GWAS files, in comparison order.
  --output-dir DIR       Output directory [./gwas_compare].
  --labels CSV           Unique plot labels [input filenames].
  --grch 37|38           Optional declared build; all input builds must agree.
  --mplot TRUE|FALSE     Shared-axis Manhattan comparison [TRUE].
  --compare-beta TRUE|FALSE  Significant-variant BETA scatter plots [TRUE].
  --compare-EAF TRUE|FALSE   Significant-variant EAF scatter plots [TRUE].
  --p-threshold FLOAT    Significance threshold [5e-8].
  --significant first|either|both  Select significance in the first file,
                         either file, or both files [first].
  Each follower is compared with the FIRST file. Required columns: SNP CHR POS
  EA NEA P, plus BETA/EAF when requested. Alleles are aligned before comparison;
  swaps change BETA to -BETA and EAF to 1-EAF. Duplicate allele/position keys and
  palindromic SNPs are excluded from scatter plots. Inputs must use the same
  build; mismatched .grch metadata is rejected. Check legacy inputs without it.
  Outputs: manhattan.compare.png, 01_vs_NN.BETA/EAF.png, harmonized TSVs,
  input_qc.tsv, comparison_qc.tsv and compare.log. No original GWAS is changed.

ldsc options:
  --gwas-files CSV       Standardized GWAS files (SNP EA NEA BETA P and N).
  --output-dir DIR       Output directory [./ldsc].
  --merge-alleles FILE   Required HapMap3 w_hm3.snplist.
  --ref-ld-chr PREFIX    Required reference LD-score prefix, with trailing / or dot.
  --w-ld-chr PREFIX      Required regression-weight prefix.
  --N NUMBER            Explicit fallback only when an input has no N column.
  --ldsc-software-dir DIR  Software checkout [/mnt/d/software/ldsc].
  --conda PATH           Conda executable [auto].
  --conda-env NAME       LDSC conda environment [ldsc].
  --python PATH         LDSC-compatible Python interpreter; bypasses conda.
  --run-munge TRUE|FALSE Run munge_sumstats.py [TRUE].
  --run-h2 TRUE|FALSE    Run ldsc.py --h2 for each input [TRUE].
  --run-rg TRUE|FALSE    Run all pairwise ldsc.py --rg comparisons [TRUE].
  --run TRUE|FALSE       FALSE writes commands only [TRUE].
  Outputs: ldsc.cmd.sh, numbered *.sumstats.gz, *.h2.log and *.rg.log.
  Invalid P is removed and P=0 is capped at 1e-300 in a separate munging copy.

Examples:
  # Compare the old HPC bald result to the new REGENIE result:
  ./gwas_compare.sh compare \
    --gwas-files /mnt/d/data.BIG/gwas/main/common/bald/gwas/bald.gz,/mnt/d/data.BIG/gwas/self/common/bald/gwas/bald.regenie.gz \
    --mplot TRUE --compare-beta TRUE --compare-EAF TRUE \
    --output-dir /mnt/d/data.BIG/gwas/self/common/bald/qc/compare

  # Compare three methods: PLINK2 vs REGENIE, then PLINK2 vs SAIGE:
  ./gwas_compare.sh compare \
    --gwas-files /mnt/d/data.BIG/gwas/self/common/height/gwas/height.plink2.gz,/mnt/d/data.BIG/gwas/self/common/height/gwas/height.regenie.gz,/mnt/d/data.BIG/gwas/self/common/height/gwas/height.saige.gz \
    --labels PLINK2,REGENIE,SAIGE --significant either --p-threshold 5e-8 \
    --output-dir /mnt/d/data.BIG/gwas/self/common/height/qc/compare

  # LDSC: replace /path/... with your actual GWAS and reference files:
  ./gwas_compare.sh ldsc --gwas-files /path/A.gz,/path/B.gz \
    --output-dir /path/ldsc-results --merge-alleles /path/w_hm3.snplist \
    --ref-ld-chr /path/eur_w_ld_chr/ --w-ld-chr /path/eur_w_ld_chr/
  # Add --run FALSE to inspect commands first, or --run-rg FALSE for h2 only.

Data processing (format/magma/liftover/cis/lead/mplot/pgs):
  /mnt/d/scripts/gwas/format_gwas.sh --help
EOF
}
case "${1:-}" in
  compare) shift; exec Rscript "$here/f/gwas_compare.R" "$@" ;;
  ldsc) shift; exec python3 "$here/f/gwas_ldsc.py" "$@" ;;
  -h|--help|"")
    usage
    ;;
  *) echo "ERROR: unknown module: $1 (expected compare or ldsc; see --help)" >&2; exit 1 ;;
esac
