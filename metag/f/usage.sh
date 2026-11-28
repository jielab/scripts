#!/usr/bin/env bash

# Standalone help must return before runtime setup or dispatch.
if [[ ${BASH_SOURCE[0]} == "$0" && ( ${1:-} == -h || ${1:-} == --help ) ]]; then
    cat <<'QUICKSTART_HELP'
bash ./usage.sh -h  # 内部工具/单独阶段；调用要求见下文
# 从所属项目主入口选择模块；本文件不代表整个项目
# 必填输入沿用现有接口，不自动补造数据或调度参数

Usage: ./metag.sh split [OPTIONS]
  --config FILE                 Runtime paths and CPU/memory resources
  --dir-metag DIR               Raw root containing dataset subdirectories
  --dat-list CSV                Dataset names, e.g. Wxiamen,england
  --sample NAME                 Restrict sample or dataset/sample; repeatable
  --batch-size N                Samples per job [10]; one resident index load per batch
  --batch LIST                  One-based batches per dataset: 1-2, 1,3-5, 3-, all
  --species-min-fragments N     Separate species file threshold [20]; counts unchanged
  --species-dir DIR             Permanent flat sample outputs
  --dir-out DIR                 Runtime jobs/logs/tmp/qc root
  --replace TRUE|FALSE          Default FALSE: reuse complete outputs; rebuild incomplete ones
  --submit hpc|local|auto        Execution machine [auto]
  --plan-only                   Write jobs without running/submitting
  --memory-mode mapped|ram       Demand paging or whole-index preload [mapped]
  --mem-gib N                    Classification allocation [128 GiB on HPC]
  --queue NAME                  Override the classification queue
  --with-pilot                  Optional prefix pilot (not a sensitivity validation)
Incomplete input mate pairs are skipped with a warning; malformed selected reads fail.

./metag.sh split --dat-list Wxiamen --batch-size 10 --batch 1-2 --submit hpc --plan-only
./metag.sh split --dat-list Wxiamen --batch-size 10 --batch 1-2 --submit hpc
./metag.sh split --dat-list Ndesert,Wxiamen,changj,england --batch-size 10
# Default --replace FALSE skips complete samples; use --replace TRUE to rerun.
The final partial batch is allowed. Larger batches reduce index reloads but increase
per-job wall time and the retry scope if input classification fails. Default: CPU only.
QUICKSTART_HELP
    exit 0
fi
set -euo pipefail
case "${1:-metag}" in
metag|metag_advanced) bash "$(dirname -- "${BASH_SOURCE[0]}")/../metag.sh" --help
;;
split) cat <<'HELP'
Usage: ./metag.sh split [OPTIONS]
  --config FILE                 Runtime paths and CPU/memory resources
  --dir-metag DIR               Raw root containing dataset subdirectories
  --dat-list CSV                Dataset names, e.g. Wxiamen,england
  --sample NAME                 Restrict sample or dataset/sample; repeatable
  --batch-size N                Samples per job [10]; one resident index load per batch
  --batch LIST                  One-based batches per dataset: 1-2, 1,3-5, 3-, all
  --species-min-fragments N     Separate species file threshold [20]; counts unchanged
  --species-dir DIR             Permanent flat sample outputs
  --dir-out DIR                 Runtime jobs/logs/tmp/qc root
  --replace TRUE|FALSE          Default FALSE: reuse complete outputs; rebuild incomplete ones
  --submit hpc|local|auto        Execution machine [auto]
  --plan-only                   Write jobs without running/submitting
  --memory-mode mapped|ram       Demand paging or whole-index preload [mapped]
  --mem-gib N                    Classification allocation [128 GiB on HPC]
  --queue NAME                  Override the classification queue
  --with-pilot                  Optional prefix pilot (not a sensitivity validation)
Incomplete input mate pairs are skipped with a warning; malformed selected reads fail.

./metag.sh split --dat-list Wxiamen --batch-size 10 --batch 1-2 --submit hpc --plan-only
./metag.sh split --dat-list Wxiamen --batch-size 10 --batch 1-2 --submit hpc
./metag.sh split --dat-list Ndesert,Wxiamen,changj,england --batch-size 10
# Default --replace FALSE skips complete samples; use --replace TRUE to rerun.
The final partial batch is allowed. Larger batches reduce index reloads but increase
per-job wall time and the retry scope if input classification fails. Default: CPU only.
HELP
;;
report|split_summary) cat <<'HELP'
Usage: bash f/analysis_entry.sh report [--sample DATASET/SAMPLE] [OPTIONS]
  --config FILE --input-root ROOT --manifest jobs.json --out DIR
  --sample may repeat; --manifest includes missing/failed requested samples.
  --compare-root ROOT compares completed results for the same input/database.
  --submit hpc|local|auto [auto] --plan-only --replace TRUE|FALSE [TRUE]
  --queue NAME --cores N --mem-gib N
Checks completion, file sizes and count conservation. No FASTQ/index reread.
Read report.html; exact tables include samples.tsv, species-by-sample.tsv and groups.tsv.
All fractions use all fragments (one paired fragment = two reads). Classification
fractions are not cell counts or confirmed organism abundances. Missing/invalid
requested results cause nonzero exit after reports are written.
HELP
;;
validate|human_validate) cat <<'HELP'
Usage: ./metag.sh validate --sample DATASET/SAMPLE [OPTIONS]
  --all runs the complete validation suite (resumes completed checks by default).
  --all --without-thresholds avoids opening the large Kraken indexes entirely.
  --all --candidates-only validates only retained human/plant candidate pairs.
  --thresholds compares c=0.1/g=2, c=0/g=2, c=0/g=1 on the same random sample.
  Default without these selectors: all saved human candidates.
  --scope full: all original sample reads (human compatibility screen).
  --plants selects bounded plant competitive validation; --reads selects read diagnostics.
  --plants --group human selects human/chimp/gorilla competitive validation.
  --all resources: 32 cores / 256 GiB / 120 hours, steps execute serially.
  --all scans all input for sequence metrics; --qc-fragments N requests a prefix.
  --output-dir SAMPLE_DIR or --target-fastq R1 [R2] may select explicit inputs.
  --reference FASTA --bwa-mem2 PATH override the prepared independent human index.
  --min-mapq 20 --min-aligned-fraction 0.80 --min-identity 0.95
  --pair-policy any_mate|both_mates [any_mate] --threads 24
  --config FILE --out DIR --replace TRUE|FALSE [TRUE]
  --submit hpc|local|auto [auto] --plan-only
  --queue 38 --cores 32 --mem-gib 128 --hours 24

Prepare the human index once: ./1.lib.sh human --submit hpc
./metag.sh validate --sample england/ERR14374050 --submit hpc
./metag.sh validate --all --sample england/ERR14374050 --submit hpc
./metag.sh validate --sample england/ERR14374050 --scope full --submit hpc
Candidate validation measures support among selected candidates, not missed reads.
Human-only alignment means human-compatible, not species-specific truth. Viral
decoys chrEBV/NC_007605.1 are excluded. Original Kraken labels remain intact.
Read human-validation.tsv and the detailed JSON. Full/candidate outputs use stable
separate directories; replacement is explicit in --replace (default TRUE).
HELP
;;
validation_suite) cat <<'HELP'
Usage: ./metag.sh validate --all --sample DATASET/SAMPLE [OPTIONS]
  --config FILE --submit hpc|local|auto --plan-only
  --replace TRUE|FALSE [FALSE]   Resume completed checks for the same request
  --threads 24 --cores 32 --mem-gib 256 --hours 120
  --candidates-only             Saved human/plant candidates; no raw scans
  --without-thresholds          All independent checks; do not query Kraken
  --qc-fragments N [0]          0 checks every pair; positive N is a prefix
  --threshold-fragments 100000 --seed 20261114
  --max-species 60 --max-reference-gib 24  Per competitive panel
  External inputs must be complete before submission; otherwise log and exit.
Prepare minimap2 once with ./install.sh --qc-tools and GRCh38 with ./1.lib.sh human.
Reports go to work_dir/qc/validate/DATASET/SAMPLE/report.html. Every step is
attempted unless its prerequisite failed. Partial reference coverage is explicit.
Classification labels and existing species FASTQs are retained.
HELP
;;
threshold_diagnose) cat <<'HELP'
Usage: ./metag.sh validate --thresholds --sample DATASET/SAMPLE [OPTIONS]
  --config FILE --fragments 100000 --seed 20261114 --threads 22
  --submit hpc|local|auto --plan-only --replace TRUE|FALSE [TRUE]
  --keep-sample --out DIR
Uses the same seeded random sample across all input ordinals for three serial
full-union queries (c=0.1/g=2, c=0/g=2, c=0/g=1). Always uses mmap. This diagnoses
threshold sensitivity; extra classified reads are not automatically correct.
HELP
;;
plant_validate) cat <<'HELP'
Usage: ./metag.sh validate --plants --sample DATASET/SAMPLE [OPTIONS]
  --taxid TAXID (repeatable) --include-broad --representatives FILE
  --competitor-taxid TAXID --max-species 60 --max-reference-gib 24
  --max-fragments 1000000 --minimap2 PATH --threads 32 --out DIR
  --group plant|human [plant]; serial bounded panels cover all selected taxa.
  --keep-alignment-files retains temporary per-panel FASTA/SAM files.
  --submit hpc|local|auto [auto] --plan-only --replace TRUE|FALSE [TRUE]
More than 60 targets are split into serial panels. Missing or oversized references
are explicitly untested; panel MAPQ is not evidence against omitted competitors.
No count floor is imposed on candidate validation. Shared species_other files are
filtered by exact per-read assignments. Pairs remain together. MAPQ/support refer
to the declared competitive panel; omitted competitors remain a limitation.
HELP
;;
sample_diagnose) cat <<'HELP'
Usage: ./metag.sh validate --reads --sample DATASET/SAMPLE [OPTIONS]
  --limit N [2000000]; 0 scans the entire original input and verifies its digest.
  --config FILE --output-dir SAMPLE_DIR --r1 FASTQ --r2 FASTQ --out DIR
  --submit hpc|local|auto [auto] --plan-only --replace TRUE|FALSE [TRUE]
Writes sample-qc.tsv and sample-diagnostic.json. Prefix QC is not a random sample.
Original raw reads are required; no new reads are downloaded.
HELP
;;
refine|species_refine) cat <<'HELP'
Usage: ./metag.sh refine [OPTIONS]
  --representatives FILE --counts counts.tsv --assignments read-assignments.tsv.gz
  --r1 FASTQ --r2 FASTQ --out DIR --minimap2 PATH --threads 32
  --candidate-min-fragments N --max-species N --max-reference-gib N
  --rescue-status unclassified|higher_rank (repeatable)
  --materialize --materialize-min-fragments 20
  --sylph-candidates FILE --include-human (human excluded by default)
  --submit hpc|local|auto [auto] --plan-only --replace TRUE|FALSE [TRUE]
Existing species assignments are retained. Higher-rank rescues must match lineage;
conflicting positive mates are not rescued. More assigned reads is not proof of
higher accuracy. --profile-only selects separate Sylph candidate discovery.
HELP
;;
sylph_profile) cat <<'HELP'
Usage: ./metag.sh refine --profile-only [OPTIONS]
  --representatives FILE --db FILE (repeatable) --r1 FASTQ --r2 FASTQ --out DIR
  --sylph PATH --threads 32 --min-ani 95
  --submit hpc|local|auto [auto] --plan-only --replace TRUE|FALSE [TRUE]
Sylph profiles genome candidates; it does not classify individual reads.
Prepare this optional library once with: ./1.lib.sh sylph --buckets prok
HELP
;;
group_summary) cat <<'HELP'
Usage: bash f/analysis_entry.sh report --corrected-groups --sample DATASET/SAMPLE --out DIR
Repeat --sample as needed. Recomputes display groups using the configured taxonomy;
original exact taxids/FASTQs are retained. --replace TRUE|FALSE [TRUE]
--submit hpc|local|auto [auto] --plan-only
HELP
;;
build_lib)
    bash "$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd -P)/1.lib.sh" --help
;;
build_kraken) cat <<'HELP'
Usage: ./1.lib.sh kraken [--config FILE] [OPTIONS]
  --submit hpc|local|auto [auto] --plan-only
  --with-split --dat-list CSV --sample NAME --batch-size N --replace TRUE|FALSE
Per-category resource budgets live in config.json. Publication depends on every
category build and performs the cheap final check. Run setup first on a fresh machine.
After publication succeeds, run build --index-use optimized.
Existing completed library
components are reused. Do not rebuild a complete library to classify new samples.
HELP
;;
build_human) cat <<'HELP'
Usage: ./1.lib.sh human [--config FILE] [--root DIR]
  --submit hpc|local|auto [auto] --plan-only
Prepares NCBI GRCh38 no-alt and pinned BWA-MEM2 2.2.1. Complete indexes are reused.
--plan-only performs no downloads/builds; execute its generated command when ready.
HUMAN_DIR, BWA_DIR, GRCH38_URL, HUMAN_QUEUE, HUMAN_CORES, HUMAN_MEM_GIB,
HUMAN_HOURS and METAG_CONFIG may override site settings.
BUILD.COMPLETE plus all expected sidecars identify a completed reference.
Human validation does not download sample reads or rebuild the Kraken library.
HELP
;;
build_sylph) cat <<'HELP'
Usage: ./1.lib.sh sylph [--config FILE] [--buckets prok|virus|euk|LIST]
  --submit hpc|local|auto [auto] --plan-only
Optional compact representative databases; default prok, using existing FASTAs.
--plan-only performs no installs/downloads/builds. Complete buckets are reused.
SYLPH_DIR, SYLPH_TOOL_DIR, SYLPH_BUCKETS, SYLPH_QUEUE, SYLPH_THREADS,
SYLPH_MEM_GIB, SYLPH_HOURS and AUDIT_DIR may override site settings.
HELP
;;
check) cat <<'HELP'
Usage: ./1.lib.sh check [--index-use optimized|full] [--full]
  ./1.lib.sh check --index-use optimized
  ./1.lib.sh check --index-use full --full  # submit a full union SHA256 scan
  ./1.lib.sh status                        # inspect availability; never builds
  --config FILE --root DIR --submit hpc|local|auto --plan-only
Quick checks inspect catalogs, headers and recorded file state. --full additionally
hashes the full union and validates current publication receipts.
See ./1.lib.sh --help for setup -> kraken -> build ordering and optional modules.
HELP
;;
compare) cat <<'HELP'
Usage: ./1.lib.sh check --benchmark [OPTIONS]
  --config FILE --dat-list CSV --sample NAME --all-samples
  --methods CSV --profiles FILE --repeats N --fragments N
  --submit hpc|local|auto [auto] --plan-only
Requires separately completed compatible native classifier databases. Uses matched
inputs and reports agreement/timing; agreement is not biological ground truth.
HELP
;;
reference_audit) cat <<'HELP'
Usage: ./1.lib.sh check --references [OPTIONS]
Inspect configured reference metadata and write representative/audit tables.
This optional operation does not rebuild the classification library.
Use --help-full for reference-audit worker options.
HELP
;;
*) echo "Unknown help topic: $1" >&2; exit 2 ;;
esac

