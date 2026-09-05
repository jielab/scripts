#!/usr/bin/env bash

set -euo pipefail

# Non-interactive WSL shells do not always execute the conda initialization
# block in ~/.bashrc. Discover the retained per-user Anaconda installation so
# preflight checks and `conda run -n le8` work from scripts and schedulers too.
if ! command -v conda >/dev/null 2>&1 && [[ -x "${HOME}/anaconda3/bin/conda" ]]; then
  export PATH="${HOME}/anaconda3/condabin:${PATH}"
fi

# Keep the original command line so the script can relaunch itself inside a
# cgroup after parsing the requested resource limits.
original_args=("$@")


# 🚩 Life's Essential 8 supervised 5C omics pipeline
usage() {
  cat <<'HELP'
Usage:
  ./le8.sh [options]

Core options:
  -Y, --trait TRAIT       Outcome under <main-project>/common/<TRAIT>/gwas/<TRAIT>.gz.
                          Default: cvd_cad. Multiple: cvd_cad,cvd_stroke_i
  -b, --biom LAYERS       prot, met, or prot,met. Default: prot,met
  -s, --steps JOBS        Comma-separated jobs.
      --from JOB           Run from JOB to the end.
      --to JOB             Stop after JOB.
      --replace TRUE|FALSE Ignore module caches and recompute. Default: FALSE.
      --run-mrlink2 Top|All|None Execute MR-link-2. Default: Top.
      --run-dandelion Top|All|None Execute protein DANDELION. Default: Top.
      --dandelion-gene-p FILE    Gene-level disease P values (WES burden preferred).
      --dandelion-allow-magma TRUE|FALSE  Permit MAGMA as sensitivity analysis. Default: TRUE.
      --dandelion-max-target-fraction FLOAT  Broad-selection audit threshold. Default: 0.25.
      --dandelion-gene-annotation FILE  Genome-wide gene coordinates for cis exclusion.
      --dandelion-snp-file FILE  Independent disease SNPs. Default: outcome *.jma.cojo.
      --dandelion-snp-gene-map FILE  Optional SNP-to-cis-gene mapping.
      --dandelion-fdr FLOAT      Target pair FDR. Default: 0.10.
      --run-gpu-coloc TRUE|FALSE Execute GPU-coloc. Default: TRUE.
      --nested-cv TRUE|FALSE     Run the expensive C5 nested CV. Default: FALSE.
      --direction-anchors CSV   C1/C2/C5 anchors. Default: PCSK9,LPA,GDF15,NTPROBNP,MMP12.
      --c4-module-boot N        C4 supervised-module bootstraps. Default: 100.
      --c4-module-stability X   C4 stable-core threshold. Default: 0.70.
      --c5-topn-max N           Largest top-N score benchmark. Default: 30.
      --c5-topn-boot N          Bootstrap replicates for score stability. Default: 150.
      --c5-lead-min-cases N     Minimum cases at headline lead time. Default: 100.
      --c5-mechanism-n N        N for score-trajectory mechanism panel. Default: 10.
      --inc_prot FILE      Restrict C5 protein prediction to biomarkers in FILE.
      --inc_met FILE       Restrict C5 metabolite prediction to biomarkers in FILE.
      --match-memory-mb N       Per-process cap for match_GRCH/match_SNP. Default: 4096.
      --match-sort-memory-mb N  GNU sort buffer per invocation. Default: 512.
      --match-tmp-dir DIR       Disk workspace for external sort. Default: /mnt/d/tmp.
      --analysis-root DIR  Analysis output root. Default: /mnt/d/analysis/le8.
      --gwas-dir DIR       Outcome GWAS project root.
      --pqtl-dir DIR       pQTL project root.
      --mqtl-dir DIR       mQTL project root.
      --prot-bed FILE       Protein annotation BED.
      --refgen-root DIR    1KG root containing <37|38>/pfile.
      --grch BUILD         auto, 37, or 38. Default: auto.
      --ukb-phe DIR        UKB phenotype directory.
      --outcome-gwas FILE  Explicit outcome GWAS for a single trait.
      --mrlink2-ref-bed FILE
      --mrlink2-ref-pfile-dir DIR
      --mrlink2-ref-pop POP       1KG super-population. Default: EUR.
      --mrlink2-ref-id-dir DIR    Persistent <POP>.id.2col directory.
      --mrlink2-ref-samples FILE  Sample table; default: sibling samples.txt.
      --shared-shell FILE  Shared 0phe.f.sh library.
      --r-bin FILE         Rscript executable.
      --cores N            R worker count. Default: 4.
      --memory-limit-gb N  Hard RAM cap for the complete LE8 process tree.
                          Default: 32. Use 0 to disable the guard.
      --memory-swap-gb N   Additional swap cap under cgroup v2. Default: 4.
      --seed N             Random seed. Default: 2026.
      --gpu-coloc-bin FILE GPU-coloc executable. Default: gpu-coloc.
      --preflight          Check scripts, paths and software, then exit.
      --dry-run            Print the execution plan, then exit.
  -h, --help               Show this help.

Job order:
  c1_correlate       Observational ProtWAS/MWAS and pre-diagnostic patterns
  c2_cause           cis/local + trans/distal MR, MR-link-2, DANDELION disease drivers
  c3_coloc           Locus-specific colocalization and fine-mapping evidence
  c4_connect         LE8-supervised proxy discovery, clustering and mediation
  c5_consolidate     Prediction panels; optional nested-CV consolidation
  s1_interact        Sex-specific omics and pairwise LE8 interactions
  s2_nonlin          Non-linear associations and CV penalty algorithms

gwas_post layout consumed by LE8:
  <project>/common/<trait>/gwas/<trait>.gz
  <project>/common/<trait>/gwas/<trait>.{cis.gz,jma.cojo}
  <main-project>/common/<trait>/magma/<trait>.genes.out

For a single-trait diagnostic run, --outcome-gwas may point directly to an
outcome file. Project-root options remain the preferred API.

Examples:
  ./le8.sh -Y cvd_cad --biom prot,met --run-mrlink2 Top
  ./le8.sh -Y cvd_cad --biom prot --steps c1_correlate,c2_cause --run-dandelion Top
  ./le8.sh -Y cvd_cad --biom prot --steps c2_cause --dandelion-gene-p /path/cvd_cad.WES_burden.tsv --dandelion-gene-annotation /path/genes.GRCh37.tsv
  ./le8.sh -Y cvd_cad,masld --biom prot,met --steps c1_correlate,c5_consolidate
  ./le8.sh -Y cvd_cad --biom met --from c4_connect
HELP
}


# 🚩 Argument normalization
bool_word() {
  case "${1,,}" in
    true|1|yes|y) printf 'TRUE\n' ;;
    false|0|no|n) printf 'FALSE\n' ;;
    *) return 2 ;;
  esac
}

run_mode() {
  # TRUE/FALSE are accepted for old command lines, but the exported contract
  # always uses the exact Top/All/None spelling.
  case "${1,,}" in
    top) printf 'Top\n' ;;
    all|true|1|yes|y) printf 'All\n' ;;
    none|false|0|no|n) printf 'None\n' ;;
    *) return 2 ;;
  esac
}


# 🚩 Defaults
trait_csv=cvd_cad
biom_csv=prot,met
step_csv=""
from_job=""
to_job=""
dry_run=FALSE
preflight_only=FALSE
replace=FALSE
RUN_MRlink2="${RUN_MRlink2:-Top}"
RUN_Dandelion="${RUN_Dandelion:-Top}"
C2_DANDELION_GENE_P_FILE=""
C2_DANDELION_GENE_ANNOTATION=""
C2_DANDELION_SNP_FILE=""
C2_DANDELION_SNP_GENE_MAP=""
C2_DANDELION_FDR=0.10
C2_DANDELION_ALLOW_MAGMA=TRUE
C2_DANDELION_MAX_TARGET_FRACTION=0.25
RUN_GPU_COLOC=TRUE
C5_NESTED_CV=FALSE
C1_DIRECTION_ANCHORS="${C1_DIRECTION_ANCHORS:-PCSK9,LPA,GDF15,NTPROBNP,MMP12}"
C1_LE4_COVARS="${C1_LE4_COVARS:-diet.pts,pa.pts,smoke.pts,sleep.pts}"
C1_FULL_LE8_SENSITIVITY="${C1_FULL_LE8_SENSITIVITY:-TRUE}"
# drug.lipid is the baseline binary indicator derived in ukb/f/phe.R from
# category 1 in UKB fields 6153 and 6177 across all available assessment
# instances.  It is preferable here to the raw pipe-delimited drug.big3*
# fields because the model expects an analysis-ready treatment covariate.
C1_TREATMENT_VARS="${C1_TREATMENT_VARS:-drug.lipid}"
C2_HERITABILITY_FILE="${C2_HERITABILITY_FILE:-}"
C2_PROT_HERITABILITY_FILE="${C2_PROT_HERITABILITY_FILE:-/mnt/d/files/ppp.h2.csv}"
C2_MET_HERITABILITY_FILE="${C2_MET_HERITABILITY_FILE:-}"
C4_MODULE_BOOT="${C4_MODULE_BOOT:-100}"
C4_MODULE_STABILITY="${C4_MODULE_STABILITY:-0.70}"
C5_TOPN_MAX="${C5_TOPN_MAX:-30}"
C5_TOPN_BOOT="${C5_TOPN_BOOT:-150}"
C5_LEAD_MIN_CASES="${C5_LEAD_MIN_CASES:-100}"
C5_MECHANISM_N="${C5_MECHANISM_N:-10}"
C5_INC_PROT=""
C5_INC_MET=""
MATCH_SNP_MEMORY_MB="${MATCH_SNP_MEMORY_MB:-4096}"
MATCH_SNP_SORT_MEMORY_MB="${MATCH_SNP_SORT_MEMORY_MB:-512}"
MATCH_SNP_TMP_DIR="${MATCH_SNP_TMP_DIR:-/mnt/d/tmp}"
LE8_ANALYSIS_ROOT=/mnt/d/analysis/le8
LE8_GWAS_DIR=/mnt/d/data.BIG/gwas/main
LE8_PQTL_IV_DIR=/mnt/d/data.BIG/gwas/prot
LE8_MQTL_IV_DIR=/mnt/d/data.BIG/gwas/met
LE8_PROT_BED=""
LE8_REFGEN_ROOT=/mnt/d/data.BIG/refGen/1kg
LE8_GRCH=auto
UKB_PHE=/mnt/d/data/ukb/phe
Y_GWAS=""
MRLINK2_REF_BED=""
MRLINK2_REF_PFILE_DIR=""
MRLINK2_REF_POP="${MRLINK2_REF_POP:-EUR}"
MRLINK2_REF_ID_DIR="${MRLINK2_REF_ID_DIR:-}"
MRLINK2_REF_SAMPLES="${MRLINK2_REF_SAMPLES:-}"
PHE_F=/mnt/d/scripts/0f/0phe.f.sh
R_BIN=Rscript
N_CORES=4
LE8_MEMORY_LIMIT_GB="${LE8_MEMORY_LIMIT_GB:-32}"
LE8_MEMORY_SWAP_GB="${LE8_MEMORY_SWAP_GB:-4}"
SEED=2026
GPU_COLOC_BIN=gpu-coloc

while [[ $# -gt 0 ]]; do
  case "$1" in
    -Y|--trait) trait_csv="$2"; shift 2 ;;
    -b|--biom) biom_csv="$2"; shift 2 ;;
    -s|--steps) step_csv="$2"; shift 2 ;;
    --from) from_job="$2"; shift 2 ;;
    --to) to_job="$2"; shift 2 ;;
    --replace) replace=$(bool_word "${2:?ERROR: --replace requires TRUE or FALSE}") || { echo "ERROR: --replace requires TRUE or FALSE" >&2; exit 2; }; shift 2 ;;
    --run-mrlink2) RUN_MRlink2=$(run_mode "${2:?ERROR: --run-mrlink2 requires Top, All, or None}") || { echo "ERROR: --run-mrlink2 requires Top, All, or None" >&2; exit 2; }; shift 2 ;;
    --run-dandelion) RUN_Dandelion=$(run_mode "${2:?ERROR: --run-dandelion requires Top, All, or None}") || { echo "ERROR: --run-dandelion requires Top, All, or None" >&2; exit 2; }; shift 2 ;;
    --dandelion-gene-p) C2_DANDELION_GENE_P_FILE="${2:?ERROR: --dandelion-gene-p requires FILE}"; shift 2 ;;
    --dandelion-gene-annotation) C2_DANDELION_GENE_ANNOTATION="${2:?ERROR: --dandelion-gene-annotation requires FILE}"; shift 2 ;;
    --dandelion-snp-file) C2_DANDELION_SNP_FILE="${2:?ERROR: --dandelion-snp-file requires FILE}"; shift 2 ;;
    --dandelion-snp-gene-map) C2_DANDELION_SNP_GENE_MAP="${2:?ERROR: --dandelion-snp-gene-map requires FILE}"; shift 2 ;;
    --dandelion-fdr) C2_DANDELION_FDR="${2:?ERROR: --dandelion-fdr requires FLOAT}"; shift 2 ;;
    --dandelion-allow-magma) C2_DANDELION_ALLOW_MAGMA=$(bool_word "${2:?ERROR: --dandelion-allow-magma requires TRUE or FALSE}") || { echo "ERROR: --dandelion-allow-magma requires TRUE or FALSE" >&2; exit 2; }; shift 2 ;;
    --dandelion-max-target-fraction) C2_DANDELION_MAX_TARGET_FRACTION="${2:?ERROR: --dandelion-max-target-fraction requires FLOAT}"; shift 2 ;;
    --run-gpu-coloc) RUN_GPU_COLOC=$(bool_word "${2:?ERROR: --run-gpu-coloc requires TRUE or FALSE}") || { echo "ERROR: --run-gpu-coloc requires TRUE or FALSE" >&2; exit 2; }; shift 2 ;;
    --nested-cv) C5_NESTED_CV=$(bool_word "${2:?ERROR: --nested-cv requires TRUE or FALSE}") || { echo "ERROR: --nested-cv requires TRUE or FALSE" >&2; exit 2; }; shift 2 ;;
    --direction-anchors) C1_DIRECTION_ANCHORS="${2:?ERROR: --direction-anchors requires CSV}"; shift 2 ;;
    --c4-module-boot) C4_MODULE_BOOT="${2:?ERROR: --c4-module-boot requires N}"; shift 2 ;;
    --c4-module-stability) C4_MODULE_STABILITY="${2:?ERROR: --c4-module-stability requires FLOAT}"; shift 2 ;;
    --c5-topn-max) C5_TOPN_MAX="${2:?ERROR: --c5-topn-max requires N}"; shift 2 ;;
    --c5-topn-boot) C5_TOPN_BOOT="${2:?ERROR: --c5-topn-boot requires N}"; shift 2 ;;
    --c5-lead-min-cases) C5_LEAD_MIN_CASES="${2:?ERROR: --c5-lead-min-cases requires N}"; shift 2 ;;
    --c5-mechanism-n) C5_MECHANISM_N="${2:?ERROR: --c5-mechanism-n requires N}"; shift 2 ;;
    --inc_prot) C5_INC_PROT="$2"; shift 2 ;;
    --inc_met) C5_INC_MET="$2"; shift 2 ;;
    --match-memory-mb) MATCH_SNP_MEMORY_MB="${2:?ERROR: --match-memory-mb requires N}"; shift 2 ;;
    --match-sort-memory-mb) MATCH_SNP_SORT_MEMORY_MB="${2:?ERROR: --match-sort-memory-mb requires N}"; shift 2 ;;
    --match-tmp-dir) MATCH_SNP_TMP_DIR="${2:?ERROR: --match-tmp-dir requires DIR}"; shift 2 ;;
    --analysis-root) LE8_ANALYSIS_ROOT="$2"; shift 2 ;;
    --gwas-dir) LE8_GWAS_DIR="$2"; shift 2 ;;
    --pqtl-dir) LE8_PQTL_IV_DIR="$2"; shift 2 ;;
    --mqtl-dir) LE8_MQTL_IV_DIR="$2"; shift 2 ;;
    --prot-bed) LE8_PROT_BED="$2"; shift 2 ;;
    --refgen-root) LE8_REFGEN_ROOT="$2"; shift 2 ;;
    --grch) LE8_GRCH="$2"; shift 2 ;;
    --ukb-phe) UKB_PHE="$2"; shift 2 ;;
    --outcome-gwas) Y_GWAS="$2"; shift 2 ;;
    --mrlink2-ref-bed) MRLINK2_REF_BED="$2"; shift 2 ;;
    --mrlink2-ref-pfile-dir) MRLINK2_REF_PFILE_DIR="$2"; shift 2 ;;
    --mrlink2-ref-pop) MRLINK2_REF_POP="$2"; shift 2 ;;
    --mrlink2-ref-id-dir) MRLINK2_REF_ID_DIR="$2"; shift 2 ;;
    --mrlink2-ref-samples) MRLINK2_REF_SAMPLES="$2"; shift 2 ;;
    --shared-shell) PHE_F="$2"; shift 2 ;;
    --r-bin) R_BIN="$2"; shift 2 ;;
    --cores) N_CORES="$2"; shift 2 ;;
    --memory-limit-gb) LE8_MEMORY_LIMIT_GB="${2:?ERROR: --memory-limit-gb requires N}"; shift 2 ;;
    --memory-swap-gb) LE8_MEMORY_SWAP_GB="${2:?ERROR: --memory-swap-gb requires N}"; shift 2 ;;
    --seed) SEED="$2"; shift 2 ;;
    --gpu-coloc-bin) GPU_COLOC_BIN="$2"; shift 2 ;;
    --preflight) preflight_only=TRUE; shift ;;
    --dry-run) dry_run=TRUE; shift ;;
    -h|--help) usage; exit 0 ;;
    *) echo "ERROR: unknown option: $1" >&2; usage >&2; exit 2 ;;
  esac
done

script_path="$(realpath -- "${BASH_SOURCE[0]}")"
script_dir="$(dirname -- "$script_path")"
fdir="$script_dir/f"
export DIRSCRIPT="$script_dir"
export LE8_FDIR="$fdir"


# 🚩 Runtime resources
[[ "$LE8_MEMORY_LIMIT_GB" =~ ^[0-9]+$ ]] || {
  echo "ERROR: --memory-limit-gb must be a non-negative integer; got '$LE8_MEMORY_LIMIT_GB'." >&2
  exit 2
}
[[ "$LE8_MEMORY_SWAP_GB" =~ ^[0-9]+$ ]] || {
  echo "ERROR: --memory-swap-gb must be a non-negative integer; got '$LE8_MEMORY_SWAP_GB'." >&2
  exit 2
}

# A cgroup cap accounts for the R parent, forked workers, and external tools as
# one unit.  This is the important distinction from ulimit, which can only cap
# each process independently.  systemd-run --scope keeps the terminal itself
# outside the capped cgroup and propagates the worker's exit status.
memory_limit_mode=disabled
if (( LE8_MEMORY_LIMIT_GB > 0 )); then
  if [[ "${LE8_MEMORY_SCOPE_ACTIVE:-0}" == 1 ]]; then
    memory_limit_mode=cgroup-v2
  else
    user_systemd_state=""
    if command -v systemd-run >/dev/null 2>&1 && command -v systemctl >/dev/null 2>&1; then
      user_systemd_state="$(systemctl --user is-system-running 2>/dev/null || true)"
    fi
    if [[ "$user_systemd_state" == running || "$user_systemd_state" == degraded ]]; then
      echo "LE8 memory guard: ${LE8_MEMORY_LIMIT_GB} GiB RAM + ${LE8_MEMORY_SWAP_GB} GiB swap (cgroup process-tree hard cap)." >&2
      if systemd-run --user --scope --quiet --collect \
          -p "MemoryMax=${LE8_MEMORY_LIMIT_GB}G" \
          -p "MemorySwapMax=${LE8_MEMORY_SWAP_GB}G" \
          -p OOMPolicy=kill \
          -- env LE8_MEMORY_SCOPE_ACTIVE=1 \
            LE8_MEMORY_LIMIT_GB="$LE8_MEMORY_LIMIT_GB" \
            LE8_MEMORY_SWAP_GB="$LE8_MEMORY_SWAP_GB" \
            bash "$script_path" "${original_args[@]}"; then
        exit 0
      else
        memory_scope_status=$?
        if (( memory_scope_status == 137 )); then
          echo "ERROR: LE8 was killed; it likely reached the ${LE8_MEMORY_LIMIT_GB} GiB cgroup memory limit." >&2
          echo "       Reduce --cores or raise --memory-limit-gb before retrying." >&2
        fi
        exit "$memory_scope_status"
      fi
    fi

    # Portable fallback for WSL distributions without a running user systemd
    # manager.  It is weaker than a cgroup because the cap applies per process,
    # but it still prevents a single R process from exhausting the whole VM.
    memory_limit_kib=$((LE8_MEMORY_LIMIT_GB * 1024 * 1024))
    if ulimit -v "$memory_limit_kib" 2>/dev/null; then
      memory_limit_mode=per-process-ulimit
      echo "WARNING: cgroup memory control is unavailable; using a ${LE8_MEMORY_LIMIT_GB} GiB per-process ulimit." >&2
    else
      echo "ERROR: unable to enforce the requested LE8 memory limit; refusing to run without a guard." >&2
      exit 2
    fi
  fi
fi
export LE8_MEMORY_LIMIT_GB LE8_MEMORY_SWAP_GB

# Forked R workers can each materialize large model frames. Four workers keep
# the 440k-participant metabolomics scans stable when other WSL jobs are active;
# callers can still override this deliberately with N_CORES.
export N_CORES
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export SEED
export LE8_REPLACE="$replace"
RUN_MRlink2=$(run_mode "$RUN_MRlink2") || { echo "ERROR: RUN_MRlink2 requires Top, All, or None" >&2; exit 2; }
RUN_Dandelion=$(run_mode "$RUN_Dandelion") || { echo "ERROR: RUN_Dandelion requires Top, All, or None" >&2; exit 2; }
export RUN_MRlink2 RUN_Dandelion RUN_GPU_COLOC C5_NESTED_CV C5_INC_PROT C5_INC_MET
export C1_DIRECTION_ANCHORS C4_MODULE_BOOT C4_MODULE_STABILITY C5_TOPN_MAX C5_TOPN_BOOT C5_LEAD_MIN_CASES C5_MECHANISM_N
export C1_LE4_COVARS C1_FULL_LE8_SENSITIVITY C1_TREATMENT_VARS
export C2_HERITABILITY_FILE C2_PROT_HERITABILITY_FILE C2_MET_HERITABILITY_FILE
export MATCH_SNP_MEMORY_MB MATCH_SNP_SORT_MEMORY_MB MATCH_SNP_TMP_DIR
export C2_DANDELION_GENE_P_FILE C2_DANDELION_GENE_ANNOTATION C2_DANDELION_SNP_FILE C2_DANDELION_SNP_GENE_MAP C2_DANDELION_FDR
export C2_DANDELION_ALLOW_MAGMA C2_DANDELION_MAX_TARGET_FRACTION
export LE8_ANALYSIS_ROOT LE8_GWAS_DIR LE8_PQTL_IV_DIR LE8_MQTL_IV_DIR LE8_REFGEN_ROOT UKB_PHE
if [[ -z "$LE8_PROT_BED" ]]; then
  for prot_bed_candidate in "$LE8_PQTL_IV_DIR/ppp_3k.b38.bed" /mnt/d/files/ppp_3k.38.bed; do
    if [[ -s "$prot_bed_candidate" ]]; then LE8_PROT_BED="$prot_bed_candidate"; break; fi
  done
  LE8_PROT_BED=${LE8_PROT_BED:-$LE8_PQTL_IV_DIR/ppp_3k.b38.bed}
fi
export LE8_PROT_BED Y_GWAS MRLINK2_REF_BED MRLINK2_REF_PFILE_DIR MRLINK2_REF_POP MRLINK2_REF_ID_DIR MRLINK2_REF_SAMPLES GPU_COLOC_BIN
case "${LE8_GRCH,,}" in
  auto|"") LE8_GRCH=auto ;;
  37|b37|grch37|hg19) LE8_GRCH=37 ;;
  38|b38|grch38|hg38) LE8_GRCH=38 ;;
  *) echo "ERROR: --grch must be auto, 37, or 38; got '$LE8_GRCH'." >&2; exit 2 ;;
esac
MRLINK2_REF_POP="${MRLINK2_REF_POP^^}"
if [[ -z "$MRLINK2_REF_POP" || "$MRLINK2_REF_POP" == *[!A-Z0-9_-]* ]]; then
  echo "ERROR: --mrlink2-ref-pop must be a simple population label; got '$MRLINK2_REF_POP'." >&2
  exit 2
fi
export LE8_GRCH
awk -v x="$C2_DANDELION_FDR" 'BEGIN{exit !(x ~ /^[0-9]*[.]?[0-9]+([eE][-+]?[0-9]+)?$/ && x+0>0 && x+0<1)}' || {
  echo "ERROR: --dandelion-fdr must be a number in (0,1); got '$C2_DANDELION_FDR'." >&2
  exit 2
}
awk -v x="$C2_DANDELION_MAX_TARGET_FRACTION" 'BEGIN{exit !(x ~ /^[0-9]*[.]?[0-9]+([eE][-+]?[0-9]+)?$/ && x+0>0 && x+0<=1)}' || {
  echo "ERROR: --dandelion-max-target-fraction must be a number in (0,1]; got '$C2_DANDELION_MAX_TARGET_FRACTION'." >&2
  exit 2
}
for pair in "--c4-module-boot:$C4_MODULE_BOOT" "--c5-topn-max:$C5_TOPN_MAX" "--c5-topn-boot:$C5_TOPN_BOOT" "--c5-lead-min-cases:$C5_LEAD_MIN_CASES" "--c5-mechanism-n:$C5_MECHANISM_N" "--match-memory-mb:$MATCH_SNP_MEMORY_MB" "--match-sort-memory-mb:$MATCH_SNP_SORT_MEMORY_MB"; do
  name=${pair%%:*}; value=${pair#*:}
  [[ "$value" =~ ^[1-9][0-9]*$ ]] || { echo "ERROR: $name must be a positive integer; got '$value'." >&2; exit 2; }
done
awk -v x="$C4_MODULE_STABILITY" 'BEGIN{exit !(x ~ /^[0-9]*[.]?[0-9]+([eE][-+]?[0-9]+)?$/ && x+0>=0 && x+0<=1)}' || {
  echo "ERROR: --c4-module-stability must be in [0,1]; got '$C4_MODULE_STABILITY'." >&2
  exit 2
}
declare -A trait_grch=()
declare -A trait_gwas=()
export PHE_F
[[ -s "$PHE_F" ]] || { echo "ERROR: missing $PHE_F" >&2; exit 2; }
# shellcheck source=/mnt/d/scripts/0f/0phe.f.sh
source "$PHE_F"
for phe_fn in check_GRCH match_SNP match_GRCH; do
  declare -F "$phe_fn" >/dev/null 2>&1 || { echo "ERROR: $phe_fn is missing from $PHE_F" >&2; exit 2; }
done

# The analysis tree may be deliberately removed between clean runs. Recreate
# its root explicitly before preflight/log setup so all later mkdir calls have
# a stable parent directory.
analysis_root="$LE8_ANALYSIS_ROOT"
mkdir -p "$analysis_root"


# 🚩 Biomarker layers
biom_csv="$(echo "$biom_csv" | tr '[:upper:]' '[:lower:]' | tr -d ' ')"
IFS=',' read -ra biom_raw <<< "$biom_csv"
biom_layers=()
for x in "${biom_raw[@]}"; do
  [[ -z "$x" ]] && continue
  case "$x" in
    prot|protein) x="prot" ;;
    met|metabolite|metabolites) x="met" ;;
    *) echo "ERROR: --biom accepts only prot, met, or prot,met; got '$x'." >&2; exit 2 ;;
  esac
  [[ " ${biom_layers[*]} " == *" $x "* ]] || biom_layers+=("$x")
done
(( ${#biom_layers[@]} > 0 )) || { echo "ERROR: no valid biomarker layer selected." >&2; exit 2; }
export PROT_DO=FALSE
export MET_DO=FALSE
for x in "${biom_layers[@]}"; do
  [[ "$x" == prot ]] && export PROT_DO=TRUE
  [[ "$x" == met ]] && export MET_DO=TRUE
done
export BIOM="$(IFS=,; echo "${biom_layers[*]}")"

jobs=(c1_correlate c2_cause c3_coloc c4_connect c5_consolidate s1_interact s2_nonlin)
files=(c1_correlate.R c2_cause.R c3_coloc.R c4_connect.R c5_consolidate.R s1_interact.R s2_nonlin.R)

job_index() {
  local q="$1" i
  for i in "${!jobs[@]}"; do
    [[ "${jobs[$i]}" == "$q" ]] && { echo "$i"; return 0; }
  done
  echo "ERROR: unknown job '$q'. Allowed: ${jobs[*]}" >&2
  return 1
}

selected=()
if [[ -n "$step_csv" ]]; then
  IFS=',' read -ra raw_selected <<< "$step_csv"
  for j in "${raw_selected[@]}"; do
    j="$(echo "$j" | xargs)"; [[ -z "$j" ]] && continue
    job_index "$j" >/dev/null || exit 2
    selected+=("$j")
  done
else
  start=0; end=$((${#jobs[@]} - 1))
  [[ -n "$from_job" ]] && start="$(job_index "$from_job")"
  [[ -n "$to_job" ]] && end="$(job_index "$to_job")"
  (( start <= end )) || { echo "ERROR: --from comes after --to." >&2; exit 2; }
  for ((i=start; i<=end; i++)); do selected+=("${jobs[$i]}"); done
fi
(( ${#selected[@]} > 0 )) || { echo "ERROR: no jobs selected." >&2; exit 2; }


# 🚩 Step selection
requested=("${selected[@]}")
if [[ " ${selected[*]} " == *" c2_cause "* || " ${selected[*]} " == *" c3_coloc "* || " ${selected[*]} " == *" c4_connect "* || " ${selected[*]} " == *" s1_interact "* ]]; then
  [[ " ${selected[*]} " == *" c1_correlate "* ]] || selected+=("c1_correlate")
fi
if [[ " ${selected[*]} " == *" c5_consolidate "* ]]; then
  for dep in c1_correlate c2_cause c3_coloc c4_connect; do
    [[ " ${selected[*]} " == *" $dep "* ]] || selected+=("$dep")
  done
fi
ordered=()
for j in "${jobs[@]}"; do
  [[ " ${selected[*]} " == *" $j "* ]] && ordered+=("$j")
done
selected=("${ordered[@]}")

outcome_gwas_file() {
  local trait="$1" explicit="$Y_GWAS"
  if [[ -n "$explicit" ]]; then
    printf '%s\n' "$explicit"
  else
    printf '%s/common/%s/gwas/%s.gz\n' "${LE8_GWAS_DIR%/}" "$trait" "$trait"
  fi
}

bed_prefix_complete() {
  local prefix="$1"
  [[ -s "${prefix}.bed" && -s "${prefix}.bim" && -s "${prefix}.fam" ]]
}

pfile_prefix_complete() {
  local prefix="$1"
  [[ -s "${prefix}.pgen" && -s "${prefix}.psam" ]] &&
    [[ -s "${prefix}.pvar" || -s "${prefix}.pvar.zst" ]]
}

plink_prefix_probe() {
  local prefix="$1" ext
  for ext in .bim .pvar .pvar.zst; do
    [[ ! -s "${prefix}${ext}" ]] || { printf '%s\n' "${prefix}${ext}"; return 0; }
  done
  return 1
}


# 🚩 Preflight
preflight() {
  local bad=0 f need_c1=FALSE need_c5=FALSE need_mrlink2=FALSE need_dandelion=FALSE need_gpu_coloc=FALSE
  local need_qtl=FALSE need_prot_annotation=FALSE explicit_outcome="$Y_GWAS"
  local dan_gene_p dan_leads common0 trait_base
  [[ " ${selected[*]} " == *" c1_correlate "* ]] && need_c1=TRUE
  [[ " ${selected[*]} " == *" c5_consolidate "* ]] && need_c5=TRUE
  [[ " ${selected[*]} " == *" c2_cause "* && "$RUN_MRlink2" != None ]] && need_mrlink2=TRUE
  [[ " ${selected[*]} " == *" c2_cause "* && "$RUN_Dandelion" != None && "$PROT_DO" == TRUE ]] && need_dandelion=TRUE
  [[ " ${selected[*]} " == *" c3_coloc "* && "$RUN_GPU_COLOC" == TRUE ]] && need_gpu_coloc=TRUE
  [[ " ${selected[*]} " == *" c2_cause "* || " ${selected[*]} " == *" c3_coloc "* ]] && need_qtl=TRUE
  [[ "$PROT_DO" == TRUE && ( " ${selected[*]} " == *" c2_cause "* || " ${selected[*]} " == *" c3_coloc "* || " ${selected[*]} " == *" c4_connect "* ) ]] && need_prot_annotation=TRUE
  echo "LE8 5C preflight"
  echo "  root: $script_dir"
  echo "  layers: $BIOM (PROT_DO=$PROT_DO, MET_DO=$MET_DO)"
  echo "  traits: $trait_csv"
  echo "  jobs: ${selected[*]}"
  echo "  parallel workers: $N_CORES (OMP/BLAS threads per worker: $OMP_NUM_THREADS)"
  if [[ "$memory_limit_mode" == cgroup-v2 ]]; then
    echo "  memory guard: ${LE8_MEMORY_LIMIT_GB} GiB RAM + ${LE8_MEMORY_SWAP_GB} GiB swap (cgroup process tree)"
  elif [[ "$memory_limit_mode" == per-process-ulimit ]]; then
    echo "  memory guard: ${LE8_MEMORY_LIMIT_GB} GiB virtual memory per process (ulimit fallback)"
  else
    echo "  memory guard: disabled"
  fi
  echo "  MR-link-2: $RUN_MRlink2; DANDELION: $RUN_Dandelion; GPU-coloc: $RUN_GPU_COLOC"
  echo "  C5 nested CV: $C5_NESTED_CV"
  echo "  directionality anchors: $C1_DIRECTION_ANCHORS"
  echo "  C1 behavioral LE4 covariates: $C1_LE4_COVARS"
  echo "  C1 full-LE8 sensitivity: $C1_FULL_LE8_SENSITIVITY"
  echo "  C1/C2 treatment covariates: $C1_TREATMENT_VARS"
  [[ -z "$C2_PROT_HERITABILITY_FILE" ]] || echo "  C2 protein SNP-heritability file: $C2_PROT_HERITABILITY_FILE"
  [[ -z "$C2_MET_HERITABILITY_FILE" ]] || echo "  C2 metabolite SNP-heritability file: $C2_MET_HERITABILITY_FILE"
  [[ -z "$C2_HERITABILITY_FILE" ]] || echo "  C2 all-layer SNP-heritability override: $C2_HERITABILITY_FILE"
  echo "  C4 module bootstrap/stability: $C4_MODULE_BOOT / $C4_MODULE_STABILITY"
  echo "  C5 top-N max/bootstrap: $C5_TOPN_MAX / $C5_TOPN_BOOT"
  echo "  SNP matching resources: memory=${MATCH_SNP_MEMORY_MB} MiB; sort=${MATCH_SNP_SORT_MEMORY_MB} MiB; tmp=$MATCH_SNP_TMP_DIR"
  echo "  GWAS projects: outcome=$LE8_GWAS_DIR; pQTL=$LE8_PQTL_IV_DIR; mQTL=$LE8_MQTL_IV_DIR"
  echo "  1KG reference root: $LE8_REFGEN_ROOT"
  [[ -z "$C5_INC_PROT" ]] || echo "  C5 protein include list: $C5_INC_PROT"
  [[ -z "$C5_INC_MET" ]] || echo "  C5 metabolite include list: $C5_INC_MET"
  [[ -z "$C5_INC_PROT" || -s "$C5_INC_PROT" ]] || { echo "  MISSING: --inc_prot file ($C5_INC_PROT)"; bad=1; }
  [[ -z "$C5_INC_MET" || -s "$C5_INC_MET" ]] || { echo "  MISSING: --inc_met file ($C5_INC_MET)"; bad=1; }
  [[ "$PROT_DO" != TRUE || -z "$C2_PROT_HERITABILITY_FILE" || -s "$C2_PROT_HERITABILITY_FILE" ]] || { echo "  MISSING: C2 protein SNP-heritability file ($C2_PROT_HERITABILITY_FILE)"; bad=1; }
  [[ "$MET_DO" != TRUE || -z "$C2_MET_HERITABILITY_FILE" || -s "$C2_MET_HERITABILITY_FILE" ]] || { echo "  MISSING: C2 metabolite SNP-heritability file ($C2_MET_HERITABILITY_FILE)"; bad=1; }
  for f in "$C2_DANDELION_GENE_P_FILE" "$C2_DANDELION_GENE_ANNOTATION" "$C2_DANDELION_SNP_FILE" "$C2_DANDELION_SNP_GENE_MAP"; do
    [[ -z "$f" || -s "$f" ]] || { echo "  MISSING: DANDELION input $f"; bad=1; }
  done
  command -v "$R_BIN" >/dev/null 2>&1 || { echo "  MISSING: Rscript ($R_BIN)"; bad=1; }
  if [[ "$need_c5" == TRUE ]] && command -v "$R_BIN" >/dev/null 2>&1 && ! "$R_BIN" -e 'quit(status=if(requireNamespace("lightgbm",quietly=TRUE))0 else 1)' >/dev/null 2>&1; then
    echo "  MISSING: R package lightgbm (required by C5 Yu-fair; update environment.yml)."
    bad=1
  fi
  if [[ "$need_dandelion" == TRUE ]] && command -v "$R_BIN" >/dev/null 2>&1 && ! "$R_BIN" -e 'quit(status=if(requireNamespace("DANDELION",quietly=TRUE))0 else 1)' >/dev/null 2>&1; then
    echo "  MISSING: R package DANDELION."
    bad=1
  fi
  [[ "$need_c1" != TRUE ]] || command -v python3 >/dev/null 2>&1 || { echo "  MISSING: python3 (required by C1 enrichment fallback)"; bad=1; }
  if [[ "$need_mrlink2" == TRUE || "$need_gpu_coloc" == TRUE ]]; then
    command -v conda >/dev/null 2>&1 || { echo "  MISSING: conda (required by selected MR-link-2/GPU-coloc steps)"; bad=1; }
  fi
  if [[ "$need_mrlink2" == TRUE ]] && ! conda run -n le8 python3 -c 'import bitarray, duckdb, numpy, pandas, pyarrow, scipy' >/dev/null 2>&1; then
    echo "  MISSING: MR-link-2 Python dependencies (including bitarray). Run: conda env update -n le8 -f $script_dir/environment.yml"
    bad=1
  fi
  for f in "${files[@]}" comm.f.R c5_pred.R c1_enrich.py c2_mr_link2.py c2_mr_link2.sh c3_coloc_GPU.py c3_coloc_GPU.sh; do
    [[ -s "$fdir/$f" ]] || { echo "  MISSING: $fdir/$f"; bad=1; }
  done
  for p in "$UKB_PHE" "$LE8_GWAS_DIR"; do
    [[ -e "$p" ]] || echo "  WARNING: path not visible now: $p"
  done
  if [[ "$need_qtl" == TRUE && "$PROT_DO" == TRUE && ! -d "$LE8_PQTL_IV_DIR" ]]; then
    echo "  MISSING: pQTL gwas_post project $LE8_PQTL_IV_DIR"
    bad=1
  fi
  if [[ "$need_qtl" == TRUE && "$MET_DO" == TRUE && ! -d "$LE8_MQTL_IV_DIR" ]]; then
    echo "  MISSING: mQTL gwas_post project $LE8_MQTL_IV_DIR"
    bad=1
  fi
  if [[ "$need_prot_annotation" == TRUE && ! -s "$LE8_PROT_BED" ]]; then
    echo "  MISSING: protein BED annotation $LE8_PROT_BED"
    bad=1
  fi
  local trait outcome expected_grch n_traits=0 ref_dir ref_probe ref_bed missing_chr chr
  IFS=',' read -ra _check_traits <<< "$trait_csv"
  for trait in "${_check_traits[@]}"; do
    trait="$(echo "$trait" | xargs)"
    [[ -z "$trait" ]] || n_traits=$((n_traits + 1))
  done
  if [[ $n_traits -eq 0 ]]; then
    echo "  ERROR: no non-empty --trait value was provided."
    return 1
  fi
  if [[ -n "$explicit_outcome" && $n_traits -ne 1 ]]; then
    echo "  ERROR: --outcome-gwas is a single file and cannot be combined with multiple --trait values."
    return 1
  fi
  if [[ $n_traits -ne 1 && ( -n "$C2_DANDELION_GENE_P_FILE" || -n "$C2_DANDELION_SNP_FILE" ) ]]; then
    echo "  ERROR: explicit --dandelion-gene-p/--dandelion-snp-file can only be used with one --trait."
    return 1
  fi
  for trait in "${_check_traits[@]}"; do
    trait="$(echo "$trait" | xargs)"
    [[ -n "$trait" ]] || continue
    outcome="$(outcome_gwas_file "$trait")"
    [[ -s "$outcome" ]] || { echo "  MISSING: outcome GWAS $outcome"; bad=1; continue; }
    outcome="$(realpath -- "$outcome")"
    trait_gwas["$trait"]="$outcome"
    if [[ "$need_dandelion" == TRUE ]]; then
      common0="$(dirname -- "$(dirname -- "$outcome")")"
      trait_base="$(basename -- "$outcome" .gz)"
      dan_gene_p="$C2_DANDELION_GENE_P_FILE"
      if [[ -z "$dan_gene_p" && "$C2_DANDELION_ALLOW_MAGMA" == TRUE ]]; then
        dan_gene_p="$common0/magma/${trait_base}.genes.out"
      fi
      dan_leads="${C2_DANDELION_SNP_FILE:-$(dirname -- "$outcome")/${trait_base}.jma.cojo}"
      if [[ -n "$dan_gene_p" && -s "$dan_gene_p" ]]; then
        if [[ -z "$C2_DANDELION_GENE_P_FILE" ]]; then
          echo "  DANDELION gene-level disease evidence: $dan_gene_p (MAGMA sensitivity; excluded from primary C5 evidence)"
        else
          echo "  DANDELION gene-level disease evidence: $dan_gene_p"
        fi
      else
        echo "  MISSING: eligible DANDELION gene-level disease evidence${dan_gene_p:+ $dan_gene_p}"
        echo "           Pass --dandelion-gene-p (WES burden preferred), or enable --dandelion-allow-magma TRUE for sensitivity analysis."
        bad=1
      fi
      if [[ -s "$dan_leads" ]]; then
        echo "  DANDELION independent disease loci: $dan_leads"
      else
        echo "  MISSING: DANDELION independent disease loci $dan_leads"
        echo "           Run gwas_post.sh lead or pass --dandelion-snp-file."
        bad=1
      fi
      [[ -n "$C2_DANDELION_GENE_ANNOTATION" ]] || echo "  WARNING: DANDELION uses assayed-protein coordinates only; pass --dandelion-gene-annotation for genome-wide cis exclusion and locus labels."
    fi
    expected_grch=""
    [[ "$LE8_GRCH" == auto ]] || expected_grch="$LE8_GRCH"
    if check_GRCH "$outcome" "$expected_grch"; then
      trait_grch["$trait"]="$CHECK_GRCH_RESULT"
      echo "  trait GWAS: $trait=$outcome (GRCh$CHECK_GRCH_RESULT)"
    else
      bad=1
      continue
    fi
    if [[ "$need_mrlink2" == TRUE ]]; then
      ref_bed="$MRLINK2_REF_BED"
      if [[ -n "$ref_bed" ]]; then
        ref_probe=""
        for f in "${ref_bed}.bed" "${ref_bed}.bim" "${ref_bed}.fam"; do
          [[ -s "$f" ]] || { echo "  MISSING: MR-link-2 reference file $f"; bad=1; }
        done
        if [[ -s "${ref_bed}.bim" ]]; then
          if check_GRCH "$ref_bed" "$CHECK_GRCH_RESULT"; then
            echo "  MR-link-2 reference: $ref_bed (GRCh$CHECK_GRCH_RESULT)"
          else
            bad=1
          fi
        fi
      else
        ref_pop="$MRLINK2_REF_POP"
        ref_dir="${MRLINK2_REF_PFILE_DIR:-${LE8_REFGEN_ROOT%/}/$CHECK_GRCH_RESULT/pfile}"
        ref_id_dir="${MRLINK2_REF_ID_DIR:-${ref_dir%/}/../id}"
        ref_samples="${MRLINK2_REF_SAMPLES:-${ref_dir%/}/../samples.txt}"
        missing_chr=""
        population_keep_required=FALSE
        ref_probe=""
        for chr in {1..22}; do
          pop_prefix="$ref_dir/${ref_pop}.chr${chr}"
          global_prefix="$ref_dir/chr${chr}"
          pop_bfile_prefix="${ref_dir%/}/../bfile/${ref_pop}/chr${chr}"
          if bed_prefix_complete "$pop_bfile_prefix"; then
            [[ -n "$ref_probe" ]] || ref_probe="${ref_dir%/}/../bfile/${ref_pop}"
            continue
          fi

          if [[ "$ref_pop" != ALL ]] &&
             { bed_prefix_complete "$pop_prefix" || pfile_prefix_complete "$pop_prefix"; }; then
            if [[ -z "$ref_probe" ]]; then
              ref_probe="$(plink_prefix_probe "$pop_prefix")"
            fi
            continue
          fi
          if bed_prefix_complete "$global_prefix" || pfile_prefix_complete "$global_prefix"; then
            [[ -n "$ref_probe" ]] || ref_probe="$ref_dir"
            [[ "$ref_pop" == ALL ]] || population_keep_required=TRUE
            continue
          fi
          missing_chr="${missing_chr}${missing_chr:+,}${chr}"
        done
        if [[ -n "$missing_chr" ]]; then
          if [[ "$preflight_only" == TRUE || -n "$MRLINK2_REF_PFILE_DIR" ]]; then
            echo "  MISSING: complete MR-link-2 1KG files (.pvar or .pvar.zst) for chromosome(s) $missing_chr under $ref_dir"
            bad=1
          else
            echo "  WARNING: incomplete default MR-link-2 reference under $ref_dir (missing chromosome(s) $missing_chr)."
            echo "           Skip MR-link-2 only; core C2 cis/trans MR and the remaining 5C jobs will continue."
            RUN_MRlink2=None
            need_mrlink2=FALSE
          fi
        else
          population_ready=TRUE
          sample_count=0
          if [[ "$population_keep_required" == TRUE ]]; then
            persistent_keep="${ref_id_dir%/}/${ref_pop}.id.2col"
            if [[ -s "$persistent_keep" ]]; then
              if ! sample_count="$(awk 'BEGIN{FS="[ \t]+"}
                  {a=$1;b=$2;gsub(/\r/,"",a);gsub(/\r/,"",b);if(NF!=2||a!=b)bad=1;n++}
                  END{if(!bad)print n+0;exit bad?2:0}' "$persistent_keep")"; then
                echo "  MISSING: invalid PLINK two-column keep file: $persistent_keep"
                population_ready=FALSE
              elif (( sample_count == 0 )); then
                echo "  MISSING: empty PLINK keep file: $persistent_keep"
                population_ready=FALSE
              else
                echo "  MR-link-2 population: $ref_pop ($sample_count samples; --keep $persistent_keep)"
              fi
            elif [[ ! -s "$ref_samples" ]]; then
              echo "  MISSING: 1KG sample table for $ref_pop --keep: $ref_samples"
              population_ready=FALSE
            elif ! sample_count="$(awk -v target="$ref_pop" 'BEGIN{FS="[ \t]+"}
                NR==1{h=toupper($NF);gsub(/\r/,"",h);if(h!="SUPER_POP"){bad=1;exit 2};next}
                {p=toupper($NF);gsub(/\r/,"",p);if(p==target)n++}
                END{if(!bad)print n+0}' "$ref_samples")"; then
              echo "  MISSING: invalid 1KG sample table (last column must be super_pop): $ref_samples"
              population_ready=FALSE
            elif (( sample_count == 0 )); then
              echo "  MISSING: no $ref_pop samples in $ref_samples"
              population_ready=FALSE
            else
              echo "  MR-link-2 population: $ref_pop ($sample_count samples; dynamic --keep fallback)"
            fi
          fi
          if [[ "$population_ready" != TRUE ]]; then
            if [[ "$preflight_only" == TRUE || -n "$MRLINK2_REF_PFILE_DIR" || -n "$MRLINK2_REF_ID_DIR" || -n "$MRLINK2_REF_SAMPLES" ]]; then
              bad=1
            else
              echo "           Skip MR-link-2 only; core C2 cis/trans MR and the remaining 5C jobs will continue."
              RUN_MRlink2=None
              need_mrlink2=FALSE
            fi
          elif check_GRCH "$ref_probe" "$CHECK_GRCH_RESULT"; then
            echo "  MR-link-2 reference: $ref_dir (GRCh$CHECK_GRCH_RESULT; prebuilt BED or .pvar[.zst])"
          else
            bad=1
          fi
        fi
      fi
    fi
  done
  [[ "$need_mrlink2" != TRUE ]] || conda run -n le8 bash -c 'command -v plink' >/dev/null 2>&1 || echo "  WARNING: MR-link-2 is selected but plink is not available in conda environment le8."
  [[ "$need_mrlink2" != TRUE || -n "$MRLINK2_REF_BED" ]] || conda run -n le8 bash -c 'command -v plink2' >/dev/null 2>&1 || echo "  WARNING: MR-link-2 pfile conversion requires plink2 in conda environment le8."
  [[ "$need_gpu_coloc" != TRUE ]] || conda run -n le8 bash -c 'command -v "$GPU_COLOC_BIN"' >/dev/null 2>&1 || echo "  WARNING: GPU-coloc is selected but $GPU_COLOC_BIN is not available in conda environment le8."
  return "$bad"
}

preflight
[[ "$preflight_only" == TRUE ]] && exit 0

IFS=',' read -ra traits <<< "$trait_csv"


# 🚩 Execution
run_one() {
  local trait="$1" job="$2" file="" i log_root log_file grch ygwas mrlink2_ref_pfile_dir mrlink2_ref_id_dir mrlink2_ref_samples
  grch="${trait_grch[$trait]:-}"
  ygwas="${trait_gwas[$trait]:-}"
  mrlink2_ref_pfile_dir="${MRLINK2_REF_PFILE_DIR:-${LE8_REFGEN_ROOT%/}/$grch/pfile}"
  mrlink2_ref_id_dir="${MRLINK2_REF_ID_DIR:-${mrlink2_ref_pfile_dir%/}/../id}"
  mrlink2_ref_samples="${MRLINK2_REF_SAMPLES:-${mrlink2_ref_pfile_dir%/}/../samples.txt}"
  [[ "$grch" == 37 || "$grch" == 38 ]] || { echo "ERROR: no detected GRCh build for trait '$trait'." >&2; exit 3; }
  [[ -n "$ygwas" && -s "$ygwas" ]] || { echo "ERROR: no resolved outcome GWAS for trait '$trait'." >&2; exit 3; }
  for i in "${!jobs[@]}"; do [[ "${jobs[$i]}" == "$job" ]] && { file="${files[$i]}"; break; }; done
  [[ -n "$file" && -s "$fdir/$file" ]] || { echo "ERROR: missing worker for $job" >&2; exit 3; }
  log_root="$analysis_root/$trait/logs"
  mkdir -p "$log_root"
  log_file="$log_root/${job}.$(date +%Y%m%d-%H%M%S).log"
  echo
  echo "=============================================================================="
  echo "JOB=$job  TRAIT=$trait  BIOM=$BIOM"
  echo "GRCH=$grch"
  echo "OUTCOME_GWAS=$ygwas"
  echo "SCRIPT_ROOT=$script_dir"
  echo "SCRIPT=$fdir/$file"
  echo "LOG=$log_file"
  echo "=============================================================================="
  [[ "$dry_run" == TRUE ]] && return 0
  Y="$trait" Y_GWAS="$ygwas" BIOM="$BIOM" LE8_GRCH="$grch" MRLINK2_REF_PFILE_DIR="$mrlink2_ref_pfile_dir" \
    MRLINK2_REF_ID_DIR="$mrlink2_ref_id_dir" \
    MRLINK2_REF_POP="$MRLINK2_REF_POP" MRLINK2_REF_SAMPLES="$mrlink2_ref_samples" \
    LE8_JOB="$job" LE8_SCRIPT="$fdir/$file" \
    "$R_BIN" -e 'Y <- Sys.getenv("Y"); BIOM <- Sys.getenv("BIOM"); LE8_JOB <- Sys.getenv("LE8_JOB"); source(Sys.getenv("LE8_SCRIPT"), chdir=TRUE)' \
    2>&1 | tee "$log_file"
}

for tr in "${traits[@]}"; do
  tr="$(echo "$tr" | xargs)"; [[ -z "$tr" ]] && continue
  for job in "${selected[@]}"; do run_one "$tr" "$job"; done
done

[[ "$dry_run" == TRUE ]] && { echo; echo "Dry-run complete."; exit 0; }

echo
echo "Done. Outputs are under $analysis_root/<trait>/ (including prot/ and met/)."
