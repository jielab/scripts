#!/usr/bin/env bash
# BEGIN QIMING RESOURCE EXAMPLES (Admin, 2026-11-05)
# Queue 33: 40 cores / 192 GiB per node. Reserve 10/20/30/40 cores by
# rounding max(required_cores/40, required_memory_GiB/192) UP to 25% steps.
# Example for a 24-core, 160-GiB CPU worker:
#   bsub -q 33 -n 40 -W 24:00 -M 160GB -R "span[hosts=1] rusage[mem=160GB/host]" ...
# A small <=10-core, <=48-GiB job uses -q 33 -n 10; request its actual RAM.
# For exact CPU counts, use Admin-suggested v4-3t or v3-6t after queueinfo.
# The examples are allocation guidance, not measured peak RAM for every stage.
# Add to the analysis command below for a 160-GiB job:
#   --submit hpc --queue 33 --job-cores 40 --job-memory-gb 160 --job-hours 24 --write-cmd
# New HPC plans auto-round slots; remove --write-cmd to submit the reviewed scope.
# Details and queue limits: ../0f/HPC_SUBMISSION_POLICY.md
# END QIMING RESOURCE EXAMPLES

# HERV analyses; the permanent completed scan is a read-only input.
set -euo pipefail
ROOT=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
WORKFLOW_NAME=allele3
# shellcheck source=f/0.run.sh
source "$ROOT/f/0.run.sh"

usage() {
    cat <<'HELP'
Usage: ./allele3.sh [MODULE] [OPTIONS]

Modules:
  (omitted)        Run ALL main analysis modules in order (default)
  primary          Compute primary enrichment from current scan counters
  archaic          Optional GU/IBDmix exploration; requires --ibdmix TRUE
  summary / spatial / windows   Enrichment and spatial analyses
  alleles / time / selection     Allele catalog and ancient analyses
  candidates / genotypes / ld    Candidates, genotype validation and LD
  batch MODULE       allofus / atlas / mfi / herv (optional analysis branches)
  research MODULE    data / allofus / allofus-cache / mfi / mapped-targets / panel / atlas / atlas-allofus / balance / space / domains / cohort / archaic-alleles
  report           Assemble the results

Default order (no MODULE argument):
  primary -> summary -> spatial -> windows -> alleles -> time
  -> selection -> candidates -> genotypes -> ld -> report
  One all.cmd and one HPC job. Modules run inside the same allocation.
  Each module starts only after the preceding module succeeds;
  a failure blocks later modules. Valid matching results are reused.
  Memory is the peak needed by any module; --job-hours covers the whole run.
  IBDmix is disabled by default, including cached archaic labels and evidence.
  With --ibdmix TRUE, archaic runs after primary as an exploratory module.
  First complete 1.pca3D.sh, 2.herv.sh and 3.scan.sh.

Options:
  --out DIR              Output [${ANALYSIS_ROOT}/allele3]
  --scan-dir DIR         Completed scan [${ANALYSIS_ROOT}/allele3/stores/gnomad_scan]
  --ibdmix TRUE|FALSE    Include exploratory IBDmix evidence [FALSE]
  --archaic-root DIR     GU cohorts [${DIR_GEN}; 1kg/37/ibdmix, ukb/37/cohortNNN/ibdmix]
  --archaic-datasets LIST  Names/globs [1kg]; auto selects all complete cohorts
  --archaic-references S  Altai,Chagyr,Denisova,Denisova25,Vindija by default
  --archaic-chromosomes S  auto includes available X/Y separately; or autosomes,X,Y / X,Y
  --archaic-cache DIR    Persistent cohort products [${ANALYSIS_ROOT}/allele3/stores/archaic_overlap]
  --stages LIST          Comma-separated stages; only requested stages run
  --metadata TSV         Sample metadata [allele3.env: ${DIR_GEN}/ancient/allele3/pca3d/metadata/sample_metadata.tsv]
  --coordinates FILE    Shared PCA coordinates [allele3.env: ${DIR_GEN}/ancient/allele3/pca3d/pca/all.sscore]
  --preflight, --dry-run Check paths and dependencies without running analysis
  --screen-log-interval S  Progress interval in seconds [300]
  --verbose-screen       Show full output on screen (always saved in logs)
  --help-full            Show detailed options
  -h, --help             Show this help

Execution (shared with ld/metag 0f/0os.sh):
  --submit local|hpc    Run locally or submit one .cmd/bsub job for main modules
  --date-scratch DATE   Dated HPC scratch root; specify the active allocation date
  --write-cmd           Generate reviewable .cmd files without submitting
  --max-running N       HPC concurrency ceiling [32]; independent of --jobs-N
  --jobs-N N            Local parallel jobs [1; maximum 8]
  --wait TRUE|FALSE     Wait for submitted jobs [FALSE]
  --queue NAME --job-cores N --job-memory-gb N --job-hours N
  --job-dir DIR         Stable submission receipts; new plan/input -> new directory
  --retry-failed        Explicit EXIT retry with pending dependency repair

Qiming HPC resource policy (administrator guidance, 2026-11-05):
  Queue 33: one node has 40 CPU cores and 192 GB memory.
  Use the larger of CPU_needed/40 and memory_needed_GB/192, then round UP:
    Larger resource fraction     --job-cores / bsub -n
       0 < fraction <= 25%        10
      25% < fraction <= 50%       20
      50% < fraction <= 75%       30
      75% < fraction <= 100%      40
  Equivalently: cores = 10 * max(ceil(CPU_needed/10), ceil(memory_GB/48)).
  Examples: 4 cores + 32 GB -> 10 cores; 4 + 64 -> 20; 24 + 160 -> 40.
  Keep the requested memory at the amount needed; round the CPU allocation.
  More than 40 cores or 192 GB does not fit one queue-33 node.
  Check queueinfo / queueinfo -l QUEUE before choosing a different queue.
  The administrator suggests v4-3t or v3-6t for high-memory work when keeping
  fewer CPU slots; their node sizes/limits differ from the queue-33 table.
  New HPC plans automatically round CPU slots up using this policy; explicit
  resource flags below make the intended allocation visible before submission.

  For a 24-core / 160-GB workload on queue 33, use:
    --queue 33 --job-cores 40 --job-memory-gb 160 --job-hours 24
  These generate the following LSF directives (40 slots on one 40-core node):
    #BSUB -q 33
    #BSUB -n 40
    #BSUB -W 24:00
    #BSUB -R "span[hosts=1] rusage[mem=160GB/host]"
  Use allele3.sh --submit hpc (or its receipt-aware submit.sh) to submit once.

Examples:
  cd /work/sph-huangj/scripts/allele3

  # Recommended: omit MODULE; one all.cmd runs all 11 main modules in one job.
  # IBDmix is disabled; no GU results are needed for this main analysis.
  # 64 GB requires 20 CPU slots on queue 33 under the administrator's policy.
  bash allele3.sh --submit hpc \
    --out /work/sph-huangj/analysis/allele3/runs/herv_20261104 \
    --ibdmix FALSE \
    --job-dir "/work/sph-huangj/analysis/allele3/jobs/all_single_$(date +%Y%m%d_%H%M%S)" \
    --queue 33 --job-cores 20 --job-memory-gb 64 --job-hours 72 \
    --wait FALSE

  # Individual-module examples below are for targeted reruns, not the normal sequence.
  # Candidate selection: 32 GB fits the first quarter of a queue-33 node.
  bash allele3.sh candidates --submit hpc \
    --out /work/sph-huangj/analysis/allele3/runs/herv_20261104 \
    --low-mappability-bed /work/sph-huangj/files/allele3/umap/hg38.k100.lt1.bed.gz \
    --job-dir "/work/sph-huangj/analysis/allele3/jobs/candidates_$(date +%Y%m%d_%H%M%S)" \
    --queue 33 --job-cores 10 --job-memory-gb 32 --job-hours 24 \
    --wait FALSE

  # After candidates: validate the selected sites with 1KG genotypes.
  # Submitted 2026-11-05: job 4387432 DONE; 4,952 targets, 2,340 samples.
  # The timestamp creates a fresh submission directory when copied next time.
  bash allele3.sh genotypes --submit hpc \
    --out /work/sph-huangj/analysis/allele3/runs/herv_20261104 \
    --targets /work/sph-huangj/analysis/allele3/runs/herv_20261104/candidates/targets.tsv \
    --kg-vcf-dir /work/sph-huangj/data/gen/1kg/38/vcf \
    --kg-sample-panel /work/sph-huangj/data/gen/1kg/38/samples.txt \
    --job-dir "/work/sph-huangj/analysis/allele3/jobs/genotypes_$(date +%Y%m%d_%H%M%S)" \
    --queue 33 --job-cores 10 --job-memory-gb 32 --job-hours 24 \
    --wait FALSE

  # High-memory example: if this workload needs 24 cores / 160 GB, reserve 40.
  # Opt-in exploration only; use a separate output from the main analysis.
  bash allele3.sh archaic --ibdmix TRUE --submit hpc --archaic-datasets 1kg \
    --archaic-root /work/sph-huangj/data/gen/1kg/37/ibdmix \
    --out /work/sph-huangj/analysis/allele3/runs/ibdmix_exploratory \
    --job-dir "/work/sph-huangj/analysis/allele3/jobs/archaic_$(date +%Y%m%d_%H%M%S)" \
    --queue 33 --job-cores 40 --job-memory-gb 160 --job-hours 24 \
    --wait FALSE
  # To retain 24 cores for the same 160 GB on an administrator-suggested queue,
  # replace the resource flags with --queue v4-3t --job-cores 24 --job-memory-gb 160
  # after checking queueinfo -l v4-3t (v3-6t is another suggested option).

  # Atlas lookup: small memory request, but queue 33 still needs 10 CPU slots.
  # Submitted after quota recovery on 2026-11-05: job 4387445 DONE.
  # All 4,952 targets annotated; observed for over 10 minutes after submission.
  bash allele3.sh batch atlas --submit hpc \
    --targets /work/sph-huangj/analysis/allele3/runs/herv_20261104/candidates/targets.tsv \
    --atlas-resources /work/sph-huangj/scripts/allele3/config/alphagenome_resources.uploaded.tsv \
    --atlas-source-version uploaded-20261104-release-unverified \
    --out /work/sph-huangj/analysis/allele3/runs/herv_20261104/atlas \
    --job-dir "/work/sph-huangj/analysis/allele3/jobs/atlas_$(date +%Y%m%d_%H%M%S)" \
    --queue 33 --job-cores 10 --job-memory-gb 8 --job-hours 6 \
    --wait FALSE

  # LD/HERV associations using the existing GRCh38 LD maps for AFR/EAS/EUR/SAS.
  # 2026-11-05: job 4387858 DONE; requested resubmission 4387881 reused valid results.
  # First calculation peaked at 2.7 GB; 8 GB RAM with 10 CPU slots fits queue 33.
  # Without --ld-covariates, these are descriptive associations.
  bash allele3.sh ld --submit hpc \
    --out /work/sph-huangj/analysis/allele3/runs/herv_20261104 \
    --ld-block-dir /work/sph-huangj/data/LD/block \
    --erv-joined /work/sph-huangj/data/annot/allele3/repeats/hg38/rmskJoinedCurrent.txt.gz \
    --chrom-sizes /work/sph-huangj/data/gen/1hgp/38/GRCH38.fasta.gz.fai \
    --job-dir "/work/sph-huangj/analysis/allele3/jobs/ld_$(date +%Y%m%d_%H%M%S)" \
    --queue 33 --job-cores 10 --job-memory-gb 8 --job-hours 24 \
    --wait FALSE
  # Completed jobs disappear from plain bjobs; this does not indicate failure.
  # Check bjobs -a JOBID or bhist -l JOBID, plus the worker/output receipts.

  # Local workstation: load its shared path defaults before expanding variables.
  source "${SCRIPTS_ROOT}/0f/0os.sh"
  ./allele3.sh batch atlas --submit local --jobs-N 1 \
    --targets ${ANALYSIS_ROOT}/allele3/final/herv_focus_20261028/candidates/targets.tsv \
    --out ${ANALYSIS_ROOT}/allele3/research/atlas \
    --job-dir ${ANALYSIS_ROOT}/allele3/jobs/atlas

  Atlas requires extracted alphagenome_variant_impact_score_snvs.tsv.gz
  and a matching .tbi or .csi under DIR_DATA/annot/allele3/alphagenome/hg38/.
  Missing required inputs report full paths and exit nonzero before submission.
  Add --write-cmd to generate scripts only; generate HPC scripts on the HPC.
HELP
}

research_main() {
    local module=${1:-help} arg python=${ALLELE_PYTHON:-}
    (($# == 0)) || shift
    local memory=${ALLELE_MEMORY_LIMIT_GB:-16} interval=${ALLELE_SCREEN_LOG_INTERVAL:-300}
    local out=${ALLELE_ANALYSIS_ROOT:-${ANALYSIS_ROOT}/allele3}/research/$module
    if [[ $module == cohort && ${1:-} =~ ^(plan|summarize)$ ]]; then out=$out/$1; fi
    local dry_run=0 native=0 verbose=0 out_explicit=0
    local -a args=("$module")
    case "$module" in
        help|-h|--help) native=1 ;;
        data|allofus|allofus-cache|mfi|mapped-targets|panel|atlas|atlas-allofus|balance|space|domains|cohort|archaic-alleles) ;;
        *) fail "Unknown research module: $module (use research --help)" ;;
    esac
    while (($#)); do
        if [[ $1 == --*=* ]]; then arg=$1; shift; set -- "${arg%%=*}" "${arg#*=}" "$@"; fi
        case "$1" in
            --python) need_value "$@"; python=$2; shift 2 ;;
            --memory-limit-gb) need_value "$@"; memory=$2; shift 2 ;;
            --screen-log-interval) need_value "$@"; interval=$2; shift 2 ;;
            --verbose-screen) verbose=1; shift ;;
            --dry-run) dry_run=1; shift ;;
            --out) need_value "$@"; out=$2; out_explicit=1; args+=("$1" "$2"); shift 2 ;;
            -h|--help) native=1; args+=("$1"); shift ;;
            *) args+=("$1"); shift ;;
        esac
    done
    if (( ! native && ! out_explicit )); then args+=(--out "$out"); fi
    positive_integer --memory-limit-gb "$memory"
    positive_integer --screen-log-interval "$interval"
    workflow_memory_check "$memory"
    python=$(workflow_python "$python" ${DIR_HOME}/anaconda3/envs/allele3-evolution/bin/python3 \
        "$HOME/anaconda3/envs/allele3-evolution/bin/python3")
    WORKFLOW_CONSOLE="$ROOT/runtime/pipeline_console.py"
    # Derived modules atomically replace their result directory. Keep a running
    # log outside that directory so its published path remains valid.
    WORKFLOW_LOG_DIR="${out%/}.logs"
    launch_workflow "research-$module" "$python" "$ROOT/f/research.py" "$out" "$memory" "$interval" \
        "$verbose" "$dry_run" "$((1-native))" "$native" "${args[@]}"
}

main() {
    if [[ ${1:-} == research ]]; then shift; research_main "$@"; return; fi
    case "${1:-}" in
        -h|--help|help) usage; return 0 ;;
        --version) printf '6.0.1\n'; return 0 ;;
        pipeline|primary|archaic|summary|spatial|alleles|candidates|genotypes|windows|time|selection|ld|report|"") ;;
        -*) ;;
        *) fail "Unknown analysis stage: $1. PCA, individual HERV and scan have separate 1.pca3D.sh / 2.herv.sh / 3.scan.sh entry points." ;;
    esac
    local command=pipeline arg python=${ALLELE_PYTHON:-}
    local memory=${ALLELE_MEMORY_LIMIT_GB:-16} interval=${ALLELE_SCREEN_LOG_INTERVAL:-300}
    local analysis_root=${ALLELE_ANALYSIS_ROOT:-${ANALYSIS_ROOT}/allele3}
    local out=${ALLELE_HERV_FOCUS_OUT:-} verbose=0 native=0
    local -a args=()
    if (($#)) && [[ $1 != -* ]]; then command=$1; shift; fi
    args+=("$command")
    [[ -z ${ALLELE_METADATA:-} ]] || args+=(--metadata "$ALLELE_METADATA")
    [[ -z ${ALLELE_COORDINATES:-} ]] || args+=(--coordinates "$ALLELE_COORDINATES")
    while (($#)); do
        if [[ $1 == --*=* ]]; then arg=$1; shift; set -- "${arg%%=*}" "${arg#*=}" "$@"; fi
        case "$1" in
            -h|--help) usage; return 0 ;;
            --help-full) args=(--help); native=1; shift ;;
            --version) printf '6.0.1\n'; return 0 ;;
            --python) need_value "$@"; python=$2; shift 2 ;;
            --memory-limit-gb) need_value "$@"; memory=$2; shift 2 ;;
            --screen-log-interval) need_value "$@"; interval=$2; shift 2 ;;
            --verbose-screen) verbose=1; shift ;;
            --preflight|--dry-run) args+=("$1"); native=1; shift ;;
            --dir-analysis-root) need_value "$@"; analysis_root=$2; args+=("$1" "$2"); shift 2 ;;
            --dir-out|--out) need_value "$@"; out=$2; shift 2 ;;
            *) args+=("$1"); shift ;;
        esac
    done
    positive_integer --memory-limit-gb "$memory"
    positive_integer --screen-log-interval "$interval"
    workflow_memory_check "$memory"
    python=$(workflow_python "$python" ${DIR_HOME}/anaconda3/envs/allele3-evolution/bin/python3 \
        "$HOME/anaconda3/envs/allele3-evolution/bin/python3")
    out=$("$python" -c 'import sys; sys.path.insert(0, sys.argv[1]); from herv_pipeline import resolve_output_dir; print(resolve_output_dir(sys.argv[2], sys.argv[3]))' \
        "$ROOT/f" "$analysis_root" "$out")
    args+=(--dir-out "$out")
    WORKFLOW_CONSOLE="$ROOT/runtime/pipeline_console.py"
    launch_workflow "$command" "$python" "$ROOT/f/herv_pipeline.py" "$out" "$memory" "$interval" \
        "$verbose" 0 "$((1-native))" "$native" "${args[@]}"
}

main "$@"

