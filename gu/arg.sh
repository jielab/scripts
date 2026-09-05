#!/usr/bin/env bash
set -euo pipefail
ROOT=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)

# Help and every public example live HERE, not in a sourced helper.
usage(){ cat <<'HELP'
Run prep_gen, then build. Build includes output validation.
Both methods default to all supported chromosomes: 1-22 and non-PAR X.
Full command output and errors are saved under <arg-dir>/log/arg.<command>.*.log.

cd /mnt/d/scripts/gu

# 1KG: all QC-passing SNPs; prepared data in pfile.4arg.
./arg.sh prep_gen --dir-gen /mnt/d/data.BIG/refGen/1kg/37 \
  --method needle --format trace --threads 8 --jobs 4

./arg.sh build --dir-gen /mnt/d/data.BIG/refGen/1kg/37 \
  --method needle --format trace --threads 8 --jobs 4

# UKB: up to 20,000 individuals; prepared data in hap.4arg.
./arg.sh prep_gen --dir-gen /mnt/h/ukbGen/37 --dir-pfile /mnt/h/ukbGen/37/hap \
  --method needle --format trace --max-individuals 20000 --threads 8 --jobs 4

./arg.sh build --dir-gen /mnt/h/ukbGen/37 --dir-pfile /mnt/h/ukbGen/37/hap \
  --method needle --format trace --max-individuals 20000 --threads 8 --jobs 4

All options and alternatives (HapMap3 / tsinfer): ./arg.sh --help-all
HELP
}

usage_all(){ cat <<'HELP'
Usage: ./arg.sh prep_gen --dir-gen DIR [options]
       ./arg.sh build    --dir-gen DIR [options]

Commands:
  prep_gen  Prepare phased genetic inputs only; does not infer an ARG.
  build     Infer/convert ARGs, optionally export TRACE, then validate outputs.
  check     Optional read-only revalidation; build already runs this check.

  --method needle|tsinfer   ARG-Needle is the default; tsinfer is an alternate TRACE backend
  --format native|trace     native method output, or add GU TRACE files [native]
  --dir-gen DIR             dataset/build root, or its pfile/hap directory
  --arg-dir DIR             reusable ARG root [dataset/build root/arg]
  --dir-pfile DIR           phased PGEN directory used by ARG-Needle
  --dir-vcf DIR             indexed phased VCF directory used by tsinfer
  --map-dir DIR             GRCh genetic maps used by ARG-Needle
  --ancestry-file FILE      sample ancestry table used to balance the panel
  --scratch DIR             native-Linux ARG-Needle scratch directory
                            default includes dataset/input/output identity hash
  --gen4arg-dir DIR         prepared genetic input: <pfile-dir>.4arg for needle,
                            <vcf-dir>.4arg for tsinfer
  --keep-snv FILE          Optional Needle SNP-ID list (e.g. HapMap3); default: all QC-passing SNPs
                            --extract is a compatibility alias; --keep selects SAMPLES
  --needle-mode auto|array|sequence  [auto]: 1KG without --keep-snv => sequence;
                            other inputs / SNP lists => array. Explicit setting overrides auto.
  --maf-min FLOAT          Needle MAF filter [0.001]
  --missing-max FLOAT      missing-call fraction [0 needle; 0.05 tsinfer]
  --mac-min N              tsinfer minimum minor-allele count [2]
  --normalize TRUE|FALSE   Needle time normalization [TRUE]
  --normalize-demography FILE  optional haploid population-size demography
  --seed N                 panel selection and inference seed [20260904]
  --features TRUE|FALSE    GRID four-population features [FALSE for trace, TRUE for native]
  --chr LIST               default: 1-22,X for both methods (autosomes + non-PAR X)
  --grch 37|38 --replace TRUE|FALSE
  --chr X                  non-PAR X supported by both methods; 1-22,X includes X
  --x-odd-male error|drop-last  Needle only: ASMC requires even haplotype count [drop-last]
                            drop-last excludes the last selected male when needed;
                            logs the excluded IID in panel/X/selection.qc.json
  X uses combined chrX.pgen/psam or chrX.vcf.gz, not a merge of male/female files.
  X keeps one haplotype per male and two per female, with original sample IDs.
  Needle X normalization is disabled unless --normalize-demography is supplied.
  --action ACTION           needle: all|prepare|infer|convert|features|affinity|check
                            tsinfer: build|prepare|check
  --max-individuals N --seed-haplotypes N --anchors-per-pop N
  --full TRUE|FALSE --threads N --jobs N
  --threads controls preparation/tsinfer; installed Needle inference is single-threaded.
  --jobs controls concurrent chromosomes (memory grows accordingly).
  prep_gen only prepares input, regardless of --format; no ARG is constructed.

Output layout under --arg-dir:
  chrN.trees + sidecars       tsinfer/tsdate trees used by GU TRACE
  argn/chrN.argn              ARG-Needle serialized ARG
  trees/chrN.trees            ARG-Needle tskit conversion
  trees/chrN.variants.tsv.gz  GRID genealogy features
  trace/METHOD/chrN.trees     optional GU TRACE-format tree
  trace/METHOD/chrN.sample_map.tsv

ARG workflow examples (WSL bash; cd /mnt/d/scripts/gu first):
  # Each command below is standalone; no Bash arrays or prior variables required.
  # Run prep_gen, then build (includes validation). Inspect chr22 before starting all supported chromosomes.
  # A and B are the main workflows. C and D are optional alternatives.
  # SNPs = biallelic A/C/G/T SNPs passing MAF/missingness QC, not every raw VCF record.
  # Default Needle MAF >= 0.001, missingness = 0; no SNP-ID list is applied by default.
  # 1KG without --keep-snv automatically uses sequence; UKB/HapMap3 use array.
  # To override automatic mode selection: --needle-mode array or sequence on BOTH steps.
  # Use forward slashes in WSL paths: /mnt/d/... (not backslashes).

  # A. 1KG + Needle: all QC-passing SNPs (default); output pfile.4arg.
  ./arg.sh prep_gen \
    --dir-gen /mnt/d/data.BIG/refGen/1kg/37 --method needle --format trace \
    --chr 22 --threads 8 --jobs 4

  ./arg.sh build \
    --dir-gen /mnt/d/data.BIG/refGen/1kg/37 --method needle --format trace \
    --chr 22 --threads 8 --jobs 4


  # Continue here only after reviewing chr22 QC, memory use and runtime.
  ./arg.sh prep_gen \
    --dir-gen /mnt/d/data.BIG/refGen/1kg/37 --method needle --format trace \
    --threads 8 --jobs 4

  ./arg.sh build \
    --dir-gen /mnt/d/data.BIG/refGen/1kg/37 --method needle --format trace \
    --threads 8 --jobs 4



  # B. UKB + Needle: phased array SNPs; output hap.4arg (20,000 individuals).
  ./arg.sh prep_gen \
    --dir-gen /mnt/h/ukbGen/37 --dir-pfile /mnt/h/ukbGen/37/hap \
    --method needle --format trace --max-individuals 20000 \
    --chr 22 --threads 8 --jobs 4

  ./arg.sh build \
    --dir-gen /mnt/h/ukbGen/37 --dir-pfile /mnt/h/ukbGen/37/hap \
    --method needle --format trace --max-individuals 20000 \
    --chr 22 --threads 8 --jobs 4


  # Continue here only after reviewing chr22 QC, memory use and runtime.
  ./arg.sh prep_gen \
    --dir-gen /mnt/h/ukbGen/37 --dir-pfile /mnt/h/ukbGen/37/hap \
    --method needle --format trace --max-individuals 20000 \
    --threads 8 --jobs 4

  ./arg.sh build \
    --dir-gen /mnt/h/ukbGen/37 --dir-pfile /mnt/h/ukbGen/37/hap \
    --method needle --format trace --max-individuals 20000 \
    --threads 8 --jobs 4



  # C. OPTIONAL alternative: 1KG + tsinfer/tsdate; output vcf.4arg using INFO/AA.
  ./arg.sh prep_gen \
    --dir-gen /mnt/d/data.BIG/refGen/1kg/37 --method tsinfer --format trace \
    --chr 22 --threads 8 --jobs 4

  ./arg.sh build \
    --dir-gen /mnt/d/data.BIG/refGen/1kg/37 --method tsinfer --format trace \
    --chr 22 --threads 8 --jobs 4


  # Continue here only after reviewing chr22 QC, memory use and runtime.
  ./arg.sh prep_gen \
    --dir-gen /mnt/d/data.BIG/refGen/1kg/37 --method tsinfer --format trace \
    --threads 8 --jobs 4

  ./arg.sh build \
    --dir-gen /mnt/d/data.BIG/refGen/1kg/37 --method tsinfer --format trace \
    --threads 8 --jobs 4



  # D. OPTIONAL alternative to A: 1KG HapMap3 + Needle (array mode).
  # Separate directories preserve the full-SNP results from A.
  ./arg.sh prep_gen \
    --dir-gen /mnt/d/data.BIG/refGen/1kg/37 --method needle --format trace \
    --keep-snv /mnt/d/data.BIG/refGen/hm3/hapmap3_r3.snp \
    --gen4arg-dir /mnt/d/data.BIG/refGen/1kg/37/pfile.4arg.hm3 \
    --arg-dir /mnt/d/data.BIG/refGen/1kg/37/arg.hm3 \
    --chr 22 --threads 8 --jobs 4

  ./arg.sh build \
    --dir-gen /mnt/d/data.BIG/refGen/1kg/37 --method needle --format trace \
    --keep-snv /mnt/d/data.BIG/refGen/hm3/hapmap3_r3.snp \
    --gen4arg-dir /mnt/d/data.BIG/refGen/1kg/37/pfile.4arg.hm3 \
    --arg-dir /mnt/d/data.BIG/refGen/1kg/37/arg.hm3 \
    --chr 22 --threads 8 --jobs 4


  # Continue here only after reviewing chr22 QC, memory use and runtime.
  ./arg.sh prep_gen \
    --dir-gen /mnt/d/data.BIG/refGen/1kg/37 --method needle --format trace \
    --keep-snv /mnt/d/data.BIG/refGen/hm3/hapmap3_r3.snp \
    --gen4arg-dir /mnt/d/data.BIG/refGen/1kg/37/pfile.4arg.hm3 \
    --arg-dir /mnt/d/data.BIG/refGen/1kg/37/arg.hm3 \
    --threads 8 --jobs 4

  ./arg.sh build \
    --dir-gen /mnt/d/data.BIG/refGen/1kg/37 --method needle --format trace \
    --keep-snv /mnt/d/data.BIG/refGen/hm3/hapmap3_r3.snp \
    --gen4arg-dir /mnt/d/data.BIG/refGen/1kg/37/pfile.4arg.hm3 \
    --arg-dir /mnt/d/data.BIG/refGen/1kg/37/arg.hm3 \
    --threads 8 --jobs 4


  # --keep-snv must be identical on preparation and inference commands.
  # Omitting it selects the full QC-passing SNP set; it does not inherit a prior list.
  # --keep selects individuals; --keep-snv selects variants by ID (one ID per line).
  # Full-density sequence inference retains more information but can be MUCH slower.
  # Installed Needle is single-threaded per chromosome; --threads does not change that.
  # --jobs 2 runs two chromosomes concurrently and increases memory requirements.
  # Do not reuse the old ~/refgen_arg_scratch/37.b37 directory.
  # Default CEU normalization is a modelling assumption for mixed populations.
  # HAPS in pfile.4arg/hap.4arg are not pfiles; do not pass them as --dir-pfile.

  # Non-PAR X. Needle's ASMC backend needs an even number of real haplotypes.
  # Default drop-last excludes one male if the selected male count is odd.
  # Use --x-odd-male error to require an already-even panel without exclusions.
  # Use the same --x-odd-male value for prep_gen and build.
  ./arg.sh build --dir-gen /mnt/d/data.BIG/refGen/1kg/37 \
    --method needle --format trace --chr X --x-odd-male drop-last --threads 8 --jobs 1

  # tsinfer retains all males and females, including an odd total haplotype count.
  ./arg.sh build --dir-gen /mnt/d/data.BIG/refGen/1kg/37 \
    --method tsinfer --format trace --chr X --threads 8 --jobs 1


HELP
}

command=${1:-help}
case "$command" in
  -h|--help|help) usage; exit 0;;
  --help-all) usage_all; exit 0;;
esac
case "$command" in
  prep_gen|build|check) shift;;
  *) echo "ERROR: unknown command: $command (use ./arg.sh -h)" >&2; exit 2;;
esac
for option in "$@"; do
  case "$option" in
    -h|--help) usage; exit 0;;
    --help-all) usage_all; exit 0;;
  esac
done
export REFGEN_GEN4ARG_ONLY=0
case "$command" in
  prep_gen) exec bash "$ROOT/f/arg.run_logged.sh" prep_gen "$ROOT/f/arg.prep_gen.sh" "$@";;
  build) exec bash "$ROOT/f/arg.run_logged.sh" build "$ROOT/f/arg.cli.sh" "$@";;
  check) exec bash "$ROOT/f/arg.run_logged.sh" check "$ROOT/f/arg.cli.sh" "$@" --action check;;
esac
