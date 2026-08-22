#!/usr/bin/env bash

set -euo pipefail

# Non-interactive WSL shells do not always execute the conda initialization
# block in ~/.bashrc. Discover a standard per-user Miniforge installation so
# preflight checks and `conda run -n le8` work from scripts and schedulers too.
if ! command -v conda >/dev/null 2>&1 && [[ -x "${HOME}/miniforge3/bin/conda" ]]; then
  export PATH="${HOME}/miniforge3/condabin:${PATH}"
fi

# Life's Essential 8 supervised 5C omics pipeline
# Root script: /mnt/d/scripts/le8/le8.sh
# Worker scripts: /mnt/d/scripts/le8/f/

usage() {
  cat <<'HELP'
Usage:
  ./le8.sh [options]

Core options:
  -Y, --trait TRAIT       Outcome name under main/clean/<TRAIT>/<TRAIT>.gz.
                          Default: cvd_cad. Multiple: cvd_cad,cvd_stroke_i
  -b, --biom LAYERS       prot, met, or prot,met. Default: prot,met
  -s, --steps JOBS        Comma-separated jobs.
      --only JOBS         Alias of --steps.
      --from JOB           Run from JOB to the end.
      --to JOB             Stop after JOB.
      --force              Ignore module cache files and recompute.
      --resume             Reuse valid module cache files (default).
      --run-mrlink2        Execute MR-link-2 (default: enabled).
      --no-mrlink2         Skip MR-link-2.
      --run-dandelion      Execute DANDELION disease-driver prioritization (default: enabled).
      --no-dandelion       Skip DANDELION.
      --run-gpu-coloc      Execute GPU-coloc (default: enabled).
      --no-gpu-coloc       Skip GPU-coloc.
      --nested-cv          Run the expensive C5 nested CV (default: disabled).
      --no-nested-cv       Skip C5 nested CV and its performance figures.
      --inc_prot FILE      Restrict C5 protein prediction to biomarkers in FILE.
      --inc_met FILE       Restrict C5 metabolite prediction to biomarkers in FILE.
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

Main paths (overridable by environment variables):
  UKB_PHE=/mnt/d/data/ukb/phe
  UKB_COMMON=/mnt/d/data/ukb/phe/common
  LE8_ANALYSIS_ROOT=/mnt/d/analysis/le8       # default; override if desired
  LE8_GWAS_DIR=/mnt/d/data.BIG/gwas/main/clean
  LE8_PQTL_IV_DIR=/mnt/d/data.BIG/gwas/ppp/clean
  LE8_MQTL_IV_DIR=/mnt/d/data.BIG/gwas/met/clean
  LE8_PPP_BED=/mnt/d/data.BIG/gwas/ppp/ppp_3k.b38.bed
  UKB_BGEN_DIR=/mnt/d/data/ukb/gen/imp

Examples:
  ./le8.sh -Y cvd_cad --biom prot,met
  ./le8.sh -Y cvd_cad,masld --biom prot,met --steps c1_correlate,c5_consolidate
  ./le8.sh -Y cvd_cad --biom met --from c4_connect
HELP
}

trait_csv="${Y:-cvd_cad}"
biom_csv="${BIOM:-prot,met}"
step_csv=""
from_job=""
to_job=""
dry_run=FALSE
preflight_only=FALSE
force=FALSE
export RUN_MRLINK2="${RUN_MRLINK2:-TRUE}"
export RUN_DANDELION="${RUN_DANDELION:-TRUE}"
export RUN_GPU_COLOC="${RUN_GPU_COLOC:-TRUE}"
export C5_NESTED_CV="${C5_NESTED_CV:-FALSE}"
export C5_INC_PROT="${C5_INC_PROT:-}"
export C5_INC_MET="${C5_INC_MET:-}"

while [[ $# -gt 0 ]]; do
  case "$1" in
    -Y|--trait) trait_csv="$2"; shift 2 ;;
    -b|--biom) biom_csv="$2"; shift 2 ;;
    -s|--steps|--only) step_csv="$2"; shift 2 ;;
    --from|--from-step) from_job="$2"; shift 2 ;;
    --to|--to-step) to_job="$2"; shift 2 ;;
    --force) force=TRUE; shift ;;
    --resume) force=FALSE; shift ;;
    --run-mrlink2) export RUN_MRLINK2=TRUE; shift ;;
    --no-mrlink2) export RUN_MRLINK2=FALSE; shift ;;
    --run-dandelion) export RUN_DANDELION=TRUE; shift ;;
    --no-dandelion) export RUN_DANDELION=FALSE; shift ;;
    --run-gpu-coloc) export RUN_GPU_COLOC=TRUE; shift ;;
    --no-gpu-coloc) export RUN_GPU_COLOC=FALSE; shift ;;
    --nested-cv) export C5_NESTED_CV=TRUE; shift ;;
    --no-nested-cv) export C5_NESTED_CV=FALSE; shift ;;
    --inc_prot) export C5_INC_PROT="$2"; shift 2 ;;
    --inc_met) export C5_INC_MET="$2"; shift 2 ;;
    --preflight) preflight_only=TRUE; shift ;;
    --dry-run) dry_run=TRUE; shift ;;
    -h|--help) usage; exit 0 ;;
    *) echo "ERROR: unknown option: $1" >&2; usage >&2; exit 2 ;;
  esac
done

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
fdir="$script_dir/f"
R_BIN="${R_BIN:-Rscript}"
export DIRSCRIPT="$script_dir"
export LE8_FDIR="$fdir"
# Forked R workers can each materialize large model frames. Four workers keep
# the 440k-participant metabolomics scans stable when other WSL jobs are active;
# callers can still override this deliberately with N_CORES.
export N_CORES="${N_CORES:-4}"
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-1}"
export OPENBLAS_NUM_THREADS="${OPENBLAS_NUM_THREADS:-1}"
export MKL_NUM_THREADS="${MKL_NUM_THREADS:-1}"
export SEED="${SEED:-2026}"
export LE8_FORCE="$force"
export LE8_ANALYSIS_ROOT="${LE8_ANALYSIS_ROOT:-/mnt/d/analysis/le8}"

# The analysis tree may be deliberately removed between clean runs. Recreate
# its root explicitly before preflight/log setup so all later mkdir calls have
# a stable parent directory.
analysis_root="$LE8_ANALYSIS_ROOT"
mkdir -p "$analysis_root"

# Normalize and validate biomarker layers.
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

# Expand required upstream jobs, then restore canonical execution order. C5
# connection-guided prediction is not valid without final C1 and C4 results.
requested=("${selected[@]}")
if [[ " ${selected[*]} " == *" c5_consolidate "* ]]; then
  for dep in c1_correlate c4_connect; do
    [[ " ${selected[*]} " == *" $dep "* ]] || selected+=("$dep")
  done
fi
ordered=()
for j in "${jobs[@]}"; do
  [[ " ${selected[*]} " == *" $j "* ]] && ordered+=("$j")
done
selected=("${ordered[@]}")

preflight() {
  local bad=0 f need_c1=FALSE need_c5=FALSE need_mrlink2=FALSE need_dandelion=FALSE need_gpu_coloc=FALSE
  [[ " ${selected[*]} " == *" c1_correlate "* ]] && need_c1=TRUE
  [[ " ${selected[*]} " == *" c5_consolidate "* ]] && need_c5=TRUE
  [[ " ${selected[*]} " == *" c2_cause "* && "$RUN_MRLINK2" == TRUE ]] && need_mrlink2=TRUE
  [[ " ${selected[*]} " == *" c2_cause "* && "$RUN_DANDELION" == TRUE && "$PROT_DO" == TRUE ]] && need_dandelion=TRUE
  [[ " ${selected[*]} " == *" c3_coloc "* && "$RUN_GPU_COLOC" == TRUE ]] && need_gpu_coloc=TRUE
  echo "LE8 5C preflight"
  echo "  root: $script_dir"
  echo "  layers: $BIOM (PROT_DO=$PROT_DO, MET_DO=$MET_DO)"
  echo "  traits: $trait_csv"
  echo "  jobs: ${selected[*]}"
  echo "  parallel workers: $N_CORES (OMP/BLAS threads per worker: $OMP_NUM_THREADS)"
  echo "  MR-link-2: $RUN_MRLINK2; DANDELION: $RUN_DANDELION; GPU-coloc: $RUN_GPU_COLOC"
  echo "  C5 nested CV: $C5_NESTED_CV"
  [[ -z "$C5_INC_PROT" ]] || echo "  C5 protein include list: $C5_INC_PROT"
  [[ -z "$C5_INC_MET" ]] || echo "  C5 metabolite include list: $C5_INC_MET"
  [[ -z "$C5_INC_PROT" || -s "$C5_INC_PROT" ]] || { echo "  MISSING: --inc_prot file ($C5_INC_PROT)"; bad=1; }
  [[ -z "$C5_INC_MET" || -s "$C5_INC_MET" ]] || { echo "  MISSING: --inc_met file ($C5_INC_MET)"; bad=1; }
  command -v "$R_BIN" >/dev/null 2>&1 || { echo "  MISSING: Rscript ($R_BIN)"; bad=1; }
  if [[ "$need_c5" == TRUE ]] && command -v "$R_BIN" >/dev/null 2>&1 && ! "$R_BIN" -e 'quit(status=if(requireNamespace("lightgbm",quietly=TRUE))0 else 1)' >/dev/null 2>&1; then
    echo "  MISSING: R package lightgbm (required by C5 Yu-fair; update environment.yml)."
    bad=1
  fi
  if [[ "$need_dandelion" == TRUE ]] && command -v "$R_BIN" >/dev/null 2>&1 && ! "$R_BIN" -e 'quit(status=if(requireNamespace("DANDELION",quietly=TRUE))0 else 1)' >/dev/null 2>&1; then
    echo "  MISSING: R package DANDELION. Install once with: Rscript $script_dir/setup_dandelion.R"
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
  for p in "${UKB_PHE:-/mnt/d/data/ukb/phe}" "${LE8_GWAS_DIR:-/mnt/d/data.BIG/gwas/main/clean}"; do
    [[ -e "$p" ]] || echo "  WARNING: path not visible now: $p"
  done
  [[ "$need_mrlink2" != TRUE ]] || conda run -n le8 bash -c 'command -v plink' >/dev/null 2>&1 || echo "  WARNING: MR-link-2 is selected but plink is not available in conda environment le8."
  [[ "$need_gpu_coloc" != TRUE ]] || conda run -n le8 bash -c 'command -v "${GPU_COLOC_BIN:-gpu-coloc}"' >/dev/null 2>&1 || echo "  WARNING: GPU-coloc is selected but gpu-coloc is not available in conda environment le8."
  return "$bad"
}

preflight
[[ "$preflight_only" == TRUE ]] && exit 0

IFS=',' read -ra traits <<< "$trait_csv"
run_one() {
  local trait="$1" job="$2" file="" i log_root log_file
  for i in "${!jobs[@]}"; do [[ "${jobs[$i]}" == "$job" ]] && { file="${files[$i]}"; break; }; done
  [[ -n "$file" && -s "$fdir/$file" ]] || { echo "ERROR: missing worker for $job" >&2; exit 3; }
  log_root="$analysis_root/$trait/logs"
  mkdir -p "$log_root"
  log_file="$log_root/${job}.$(date +%Y%m%d-%H%M%S).log"
  echo
  echo "=============================================================================="
  echo "JOB=$job  TRAIT=$trait  BIOM=$BIOM"
  echo "SCRIPT=$fdir/$file"
  echo "LOG=$log_file"
  echo "=============================================================================="
  [[ "$dry_run" == TRUE ]] && return 0
  Y="$trait" BIOM="$BIOM" LE8_JOB="$job" LE8_SCRIPT="$fdir/$file" \
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
