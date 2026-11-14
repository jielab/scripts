#!/usr/bin/env bash

# Display standalone help before loading runtimes or scheduling work.
if [[ ${BASH_SOURCE[0]} == "$0" && ( ${1:-} == -h || ${1:-} == --help ) ]]; then
    cat <<'PUBLIC_ENTRY_HELP'
Metagenomics analysis

Modules:
    split     Classify paired reads, preserving each pair in one bin
    report    Summarize committed counts
    validate  Independently validate candidates or inspect original reads
    refine    Optional competitive alignment or Sylph profiling

Main options:
    --config FILE              Runtime paths and resources [config.json]
    --dat-list CSV             Dataset selection
    --sample DATASET/SAMPLE    Select one sample
    --replace TRUE|FALSE       Split default FALSE (reuse); other modules default TRUE
    --submit hpc|local|auto     Machine selection [auto]
    --plan-only                Preview commands without running
    MODULE -h / --help-full    Module options / full worker options

Reuse existing libraries; build separately when needed with ./1.lib.sh.

Examples:
    # All batches; completed samples are skipped by default.
    ./metag.sh split --dat-list Ndesert,Wxiamen,changj,england --batch-size 10
    ./metag.sh split --config config.json --dat-list england --plan-only
    ./metag.sh split --config config.json --dat-list england
    # After classification completes, summarize or validate one sample.
    ./metag.sh report --config config.json --sample england/ERR14374050
    ./metag.sh validate --config config.json --sample england/ERR14374050
PUBLIC_ENTRY_HELP
    exit 0
fi
set -euo pipefail

usage() {
case "${1:-metag}" in
metag|metag_advanced) bash "${BASH_SOURCE[0]}" --help
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
  --queue NAME                  Classification queue; CPU/memory unchanged
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
Usage: ./metag.sh report [--sample DATASET/SAMPLE] [OPTIONS]
  --config FILE --input-root ROOT --manifest jobs.json --out DIR
  --sample may repeat; --manifest includes missing/failed requested samples.
  --compare-root ROOT compares completed results for the same input/database.
  --submit hpc|local|auto [auto] --plan-only --replace TRUE|FALSE [TRUE]
  --wait-for JOBID --wait-mode done|ended --queue NAME --cores N --mem-gib N
Checks completion, file sizes and count conservation. No FASTQ/index reread.
Read report.html; exact tables include samples.tsv, species-by-sample.tsv and groups.tsv.
All fractions use all fragments (one paired fragment = two reads). Classification
fractions are not cell counts or confirmed organism abundances. Missing/invalid
requested results cause nonzero exit after reports are written.
HELP
;;
validate|human_validate) cat <<'HELP'
Usage: ./metag.sh validate --sample DATASET/SAMPLE [OPTIONS]
  Default: all saved human candidates. --scope full: all original sample reads.
  --plants selects plant competitive validation; --reads selects read diagnostics.
  --output-dir SAMPLE_DIR or --target-fastq R1 [R2] may select explicit inputs.
  --reference FASTA --bwa-mem2 PATH override the prepared independent human index.
  --min-mapq 20 --min-aligned-fraction 0.80 --min-identity 0.95
  --pair-policy any_mate|both_mates [any_mate] --threads 24
  --config FILE --out DIR --replace TRUE|FALSE [TRUE]
  --submit hpc|local|auto [auto] --plan-only --wait-for JOBID
  --queue 38 --cores 32 --mem-gib 128 --hours 24

Prepare the human index once: ./1.lib.sh human --submit hpc
./metag.sh validate --sample england/ERR14374050 --submit hpc
./metag.sh validate --sample england/ERR14374050 --scope full --submit hpc
Candidate validation measures support among selected candidates, not missed reads.
Human-only alignment means human-compatible, not species-specific truth. Viral
decoys chrEBV/NC_007605.1 are excluded. Original Kraken labels remain intact.
Read human-validation.tsv and the detailed JSON. Full/candidate outputs use stable
separate directories; replacement is explicit in --replace (default TRUE).
HELP
;;
plant_validate) cat <<'HELP'
Usage: ./metag.sh validate --plants --sample DATASET/SAMPLE [OPTIONS]
  --taxid TAXID (repeatable) --include-broad --representatives FILE
  --competitor-taxid TAXID --max-species 60 --max-reference-gib 24
  --max-fragments 1000000 --minimap2 PATH --threads 32 --out DIR
  --submit hpc|local|auto [auto] --plan-only --replace TRUE|FALSE [TRUE]
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
Usage: ./metag.sh report --corrected-groups --sample DATASET/SAMPLE --out DIR
Repeat --sample as needed. Recomputes display groups using the configured taxonomy;
original exact taxids/FASTQs are retained. --replace TRUE|FALSE [TRUE]
--submit hpc|local|auto [auto] --plan-only
HELP
;;
build_lib) cat <<'HELP'
Usage: ./1.lib.sh check|kraken|human|sylph [OPTIONS]
  check       Check the existing library; never builds it or submits split.
  kraken      Build partitioned indexes when needed.
  human       Prepare the independent human reference.
  sylph       Prepare the optional Sylph library.

Checks:
  ./1.lib.sh check                  Quick catalog/files/header check
  ./1.lib.sh check --full           Submit one full SHA256 check on HPC
  ./1.lib.sh check --full --plan-only  Preview without submission
  ./1.lib.sh check --full -o /work/sph-huangj/analysis/metag/lib
  ./1.lib.sh check --status         Read saved build status
  --config FILE                    Defaults to this project's config.json
  -o, --out-dir DIR                 Full-check scheduler log directory
Default logs: ${ANALYSIS_ROOT}/metag/lib/check.%J.out and check.%J.err.
HPC/local is detected by 0f/0os.sh. HPC full check: metag.check, queue 33,
40 cores, 160 GiB, 120 hours. Inside an allocation, check runs directly.
Full check is independent and need only be run once for this library.
It never submits metag.sh. Submit sample analysis yourself when ready:
  ./metag.sh split --dat-list Ndesert,Wxiamen,changj,england --batch-size 10 --batch 1
Keep config.json: it supplies reference paths, resources and analysis settings.

HELP
;;
build_kraken) cat <<'HELP'
Usage: ./1.lib.sh kraken [--config FILE] [OPTIONS]
  --submit hpc|local|auto [auto] --plan-only
  --with-split --dat-list CSV --sample NAME --batch-size N --replace TRUE|FALSE
Per-category resource budgets live in config.json. Publication depends on every
category build and performs the cheap final check. Existing completed library
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
Usage: ./1.lib.sh check|kraken|human|sylph [OPTIONS]
  check       Check the existing library; never builds it or submits split.
  kraken      Build partitioned indexes when needed.
  human       Prepare the independent human reference.
  sylph       Prepare the optional Sylph library.

Checks:
  ./1.lib.sh check                  Quick catalog/files/header check
  ./1.lib.sh check --full           Submit one full SHA256 check on HPC
  ./1.lib.sh check --full --plan-only  Preview without submission
  ./1.lib.sh check --full -o /work/sph-huangj/analysis/metag/lib
  ./1.lib.sh check --status         Read saved build status
  --config FILE                    Defaults to this project's config.json
  -o, --out-dir DIR                 Full-check scheduler log directory
Default logs: ${ANALYSIS_ROOT}/metag/lib/check.%J.out and check.%J.err.
HPC/local is detected by 0f/0os.sh. HPC full check: metag.check, queue 33,
40 cores, 160 GiB, 120 hours. Inside an allocation, check runs directly.
Full check is independent and need only be run once for this library.
It never submits metag.sh. Submit sample analysis yourself when ready:
  ./metag.sh split --dat-list Ndesert,Wxiamen,changj,england --batch-size 10 --batch 1
Keep config.json: it supplies reference paths, resources and analysis settings.

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

}
HERE="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
module=${1:---help}
if (($#)); then shift; fi
case "$module" in
    -h|--help) usage "metag"; exit 0 ;;
    --help-full) usage "metag_advanced"; exit 0 ;;
    report|validate|refine) exec bash "$HERE/f/analysis_entry.sh" "$module" "$@" ;;
    audit-reference)
        if [[ ${1:-} == --help-full ]]; then
            shift; source "$HERE/f/env.sh"
            exec "$METAG_PYTHON" "$HERE/f/reference_audit.py" --help "$@"
        fi
        for arg in "$@"; do
            if [[ $arg == -h || $arg == --help ]]; then usage "reference_audit"; exit 0; fi
        done
        source "$HERE/f/env.sh"
        exec "$METAG_PYTHON" "$HERE/f/reference_audit.py" "$@" ;;
    split|compare|status|audit|doctor|submit|selftest) ;;
    *) echo "Unknown analysis module: $module. Use 'bash metag.sh -h'; database builds use 'bash 1.lib.sh -h'." >&2; exit 2 ;;
esac
if [[ $module == split || $module == compare ]]; then
    for arg in "$@"; do
        if [[ $arg == -h || $arg == --help ]]; then usage "$module"; exit 0; fi
    done
fi
if [[ ${1:-} == --help-full ]]; then shift; set -- --help "$@"; fi
set -- "$module" "$@"
METAG_CODE_ROOT="$HERE"
source "$METAG_CODE_ROOT/f/env.sh"
workflow_args "$@"
set -- "${WORKFLOW_ARGS[@]}"
"$METAG_PYTHON" -c 'import sys; sys.exit(0 if sys.version_info >= (3, 10) else "metag HPC needs Python >= 3.10; load a site Python module or set METAG_PYTHON in env.sh")'
"$METAG_PYTHON" -c 'import yaml' || { echo 'PyYAML is required in METAG_PYTHON (python -m pip install PyYAML).' >&2; exit 1; }
exec "$METAG_PYTHON" "$METAG_CODE_ROOT/f/parallel_cli.py" "$@"

