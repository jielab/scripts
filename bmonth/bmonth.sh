#!/usr/bin/env bash
set -euo pipefail
ORIGINAL_ARGS=("$@")
# The shared logger bypasses an empty argument list. Normalize the default
# invocation so a delayed `./bmonth.sh` run also keeps a persistent log.
if [[ ${#ORIGINAL_ARGS[@]} -eq 0 ]]; then ORIGINAL_ARGS=(all); fi
usage(){ cat <<'HELP'
Birth-month / birth-location analysis in UK Biobank

cd /mnt/d/scripts/bmonth
./bmonth.sh all
./bmonth.sh prep
./bmonth.sh phewas
./bmonth.sh assoc
./bmonth.sh med
./bmonth.sh final
./bmonth.sh plot
./bmonth.sh check

Options:
  --out DIR              output directory [ /mnt/d/analysis/bmonth ]
  --in DIR               phenotype directory containing Rdata/ [ /mnt/d/data/ukb/phe ]
  --rscript FILE         Rscript binary
  --outcomes CSV         incident outcomes [cvd_cad,ra]
  --exposures CSV        scalar exposures for PheWAS/mediation
  --phewas-cores N|auto  PheWAS workers [auto: up to 8, limited by available RAM]
  --med-boot N           mediation bootstrap replicates [100]
  --follow-end DATE      override reconstructed-outcome end
  --replace              recompute matching PheWAS caches

All analysis stages share a 24 GiB RAM cap across the entire process tree,
with swap disabled. Requires systemd user scopes and cgroup v2.
HELP
}
STEP=all
OUTDIR=/mnt/d/analysis/bmonth
INDIR=/mnt/d/data/ukb/phe
RSCRIPT_BIN=""
OUTCOMES="cvd_cad,ra"
EXPOSURES="bmonth_cos,bmonth_sin,bmonth_winter_summer,bmonth_spring_autumn,birth_north_100km,birth_east_100km"
PHEWAS_CORES=auto
MED_BOOT=100
FOLLOW_END=""
REPLACE=FALSE
while [[ $# -gt 0 ]]; do
  case "$1" in
    prep|assoc|phewas|med|final|plot|all|check) STEP="$1"; shift ;;
    --out) OUTDIR="${2:?}"; shift 2 ;;
    --in) INDIR="${2:?}"; shift 2 ;;
    --rscript) RSCRIPT_BIN="${2:?}"; shift 2 ;;
    --outcomes) OUTCOMES="${2:?}"; shift 2 ;;
    --exposures) EXPOSURES="${2:?}"; shift 2 ;;
    --phewas-cores) PHEWAS_CORES="${2:?}"; shift 2 ;;
    --med-boot) MED_BOOT="${2:?}"; shift 2 ;;
    --follow-end) FOLLOW_END="${2:?}"; shift 2 ;;
    --replace) REPLACE=TRUE; shift ;;
    -h|--help) usage; exit 0 ;;
    *) echo "ERROR: unknown argument $1" >&2; usage >&2; exit 2 ;;
  esac
done
[[ "$PHEWAS_CORES" == auto || "$PHEWAS_CORES" =~ ^[1-9][0-9]*$ ]] || { echo 'ERROR: --phewas-cores must be auto or a positive integer' >&2; exit 2; }
[[ "$MED_BOOT" =~ ^(0|[1-9][0-9]*)$ ]] || { echo 'ERROR: --med-boot must be nonnegative integer' >&2; exit 2; }
RSCRIPT_BIN=${RSCRIPT_BIN:-$(command -v Rscript || true)}
[[ -n "$RSCRIPT_BIN" ]] || { echo "ERROR: Rscript not found" >&2; exit 1; }
SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
mkdir -p "$OUTDIR"
export BMONTH_OUTDIR="$OUTDIR" BMONTH_INDIR="$INDIR" BMONTH_OUTCOMES="$OUTCOMES" BMONTH_EXPOSURES="$EXPOSURES"
export BMONTH_PHEWAS_CORES="$PHEWAS_CORES" BMONTH_MED_BOOT="$MED_BOOT" BMONTH_FOLLOW_END="$FOLLOW_END"
export BMONTH_REPLACE="$REPLACE" BMONTH_FDIR="$SCRIPT_DIR/f"
export SCRIPT_LOG_DIR="${SCRIPT_LOG_DIR:-$OUTDIR/logs}"
source /mnt/d/scripts/0f/console.sh "${ORIGINAL_ARGS[@]}" 2>/dev/null || true
# Isolate R and every worker from the WSL terminal's init.scope. The shared
# guard verifies the kernel limits on re-entry and refuses to run uncapped.
source "$SCRIPT_DIR/../0f/memory_cap.sh"
memory_cap_enter bmonth 24 0 "$SCRIPT_DIR/bmonth.sh" "${ORIGINAL_ARGS[@]}"
run_r(){
  local s="$1"
  echo "[bmonth] START ${s##*/}"
  if [[ ${s##*/} == 3.phewas.R ]]; then
    # Parallelize independent models, not BLAS threads inside every worker.
    OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 BLIS_NUM_THREADS=1 \
      VECLIB_MAXIMUM_THREADS=1 "$RSCRIPT_BIN" --vanilla "$s"
  else
    "$RSCRIPT_BIN" --vanilla "$s"
  fi
  echo "[bmonth] DONE ${s##*/}"
}
case "$STEP" in
  check) run_r "$SCRIPT_DIR/f/00.check.R" ;;
  prep) run_r "$SCRIPT_DIR/f/1.prep.R" ;;
  assoc) run_r "$SCRIPT_DIR/f/2.assoc.R" ;;
  phewas) run_r "$SCRIPT_DIR/f/3.phewas.R" ;;
  med) run_r "$SCRIPT_DIR/f/4.med.R" ;;
  final) run_r "$SCRIPT_DIR/f/5.final.R" ;;
  plot) run_r "$SCRIPT_DIR/f/6.replot.R" ;;
  all)
    run_r "$SCRIPT_DIR/f/00.check.R"
    run_r "$SCRIPT_DIR/f/1.prep.R"
    run_r "$SCRIPT_DIR/f/3.phewas.R"
    run_r "$SCRIPT_DIR/f/2.assoc.R"
    run_r "$SCRIPT_DIR/f/4.med.R"
    run_r "$SCRIPT_DIR/f/5.final.R"
    ;;
esac
