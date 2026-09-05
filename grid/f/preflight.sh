#!/usr/bin/env bash
set -euo pipefail
ROOT=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd -P)
# --scope belongs only to preflight; remove it before the common parser.
scope=all
args=()
while (($#)); do
  case "$1" in
    --scope) [[ $# -ge 2 ]] || { echo 'ERROR: --scope requires core|arg|grid|all' >&2; exit 2; }; scope=${2,,}; shift 2;;
    *) args+=("$1"); shift;;
  esac
done
# shellcheck source=f/common.sh
source "$ROOT/f/common.sh" "${args[@]}"
grid_parse_args "${args[@]}"
case "$scope" in core|arg|grid|all) ;; *) echo "ERROR: bad --scope=$scope" >&2; exit 2;; esac

report="$GRID_OUTPUT_ROOT/preflight.${scope}.tsv"
mkdir -p "$(dirname "$report")"
printf 'scope\tstatus\titem\tdetail\n' > "$report"
fail=0
record(){ local st=$1 item=$2 detail=${3:-}; printf '%s\t%s\t%s\t%s\n' "$scope" "$st" "$item" "$detail" | tee -a "$report"; [[ $st != FAIL ]] || fail=1; }
check_cmd(){ if command -v "$1" >/dev/null 2>&1; then record PASS "command:$1" "$(command -v "$1")"; else record FAIL "command:$1" missing; fi; }
check_file(){ if [[ -s $1 ]]; then record PASS "$2" "$1"; else record FAIL "$2" "$1"; fi; }
check_dir(){ if [[ -d $1 ]]; then record PASS "$2" "$1"; else record FAIL "$2" "$1"; fi; }
check_py(){ local mod=$1; if python3 - "$mod" <<'PY' >/dev/null 2>&1
import importlib,sys
importlib.import_module(sys.argv[1])
PY
  then record PASS "python:$mod" importable; else record FAIL "python:$mod" not_importable; fi
}
check_r(){ local pkg=$1; if Rscript -e "quit(status=ifelse(requireNamespace('$pkg',quietly=TRUE),0,1))" >/dev/null 2>&1; then record PASS "R:$pkg" installed; else record FAIL "R:$pkg" missing; fi; }

if [[ $scope == core || $scope == all ]]; then
  for x in bash awk sed grep sort gzip python3 Rscript plink2; do check_cmd "$x"; done
  if find "$GRID_TARGET_DIR" -maxdepth 1 -name '*.pvar.zst' -print -quit 2>/dev/null | grep -q .; then check_cmd zstdcat; fi
  for m in numpy pandas scipy h5py; do check_py "$m"; done
  for p in data.table ggplot2 patchwork openxlsx; do check_r "$p"; done
  check_file "$GRID_PHE_FILE" phenotype_rds
  check_file "$GRID_PCA_WEIGHT" DiscoDivas_pca_weight
  check_file "$GRID_MED_FILE" DiscoDivas_reference_centers
  check_dir "$GRID_CSX_REF_DIR" PRSCSx_reference_root
  check_file "$GRID_CSX_SNPINFO" PRSCSx_snpinfo
  check_file "$GRID_CSX_BIM_PREFIX.bim" target_validation_bim
  prscx=$(find "$GRID_ROOT/f/csx" -maxdepth 3 -type f \( -iname 'PRScsx.py' -o -iname 'prscsx.py' \) -print -quit 2>/dev/null || true)
  [[ -n $prscx ]] && record PASS PRSCSx_program "$prscx" || record FAIL PRSCSx_program "$GRID_ROOT/f/csx"
  for c in $(grid_expand_chrs "$GRID_CHRS" | awk '{print $1}'); do
    if grid_target_mode "$c" >/dev/null 2>&1; then record PASS "target_genotype_chr$c" "$GRID_TARGET_DIR"; else record FAIL "target_genotype_chr$c" "$GRID_TARGET_DIR"; fi
  done
  for t in $(grid_csv_words "$GRID_TRAITS"); do
    for p in $(grid_csv_words "$GRID_POPS"); do
      if f=$(grid_find_gwas "$t" "$p"); then record PASS "gwas:${t}.${p}" "$f"; else record FAIL "gwas:${t}.${p}" "$GRID_GWAS_DIR/${t}.${p}.gz"; fi
    done
  done
fi

if [[ $scope == arg || $scope == all ]]; then
  for x in gzip python3; do check_cmd "$x"; done
  check_py tskit
  for c in $(grid_expand_chrs "$GRID_CHRS"); do
    check_file "$GRID_ARG_OUT/argn/chr$c.argn" "ARG_Needle_argn_chr$c"
    check_file "$GRID_ARG_TREES_DIR/chr$c.trees" "ARG_Needle_trees_chr$c"
    check_file "$GRID_ARG_TREES_DIR/chr$c.sample_map.tsv" "ARG_Needle_sample_map_chr$c"
    check_file "$GRID_ARG_TREES_DIR/chr$c.anchors.tsv" "ARG_Needle_anchors_chr$c"
    check_file "$GRID_ARG_TREES_DIR/chr$c.variants.tsv.gz" "ARG_Needle_features_chr$c"
  done
  if ((fail)); then
    record WARN ARG_build_command "bash /mnt/d/scripts/gu/arg.sh build --method needle --dir-gen $(dirname -- "$GRID_ARG_HAP_DIR") --chr $GRID_CHRS"
  fi
fi

if [[ $scope == grid || $scope == all ]]; then
  for m in numpy pandas scipy h5py tskit; do check_py "$m"; done
  check_file "$GRID_MED_FILE" global_PCA_centers
  # The LD score files are generated from the already downloaded PRS-CSx HDF5 panels.
  for p in $(grid_csv_words "$GRID_POPS"); do
    low=${p,,}; d=$(find "$GRID_CSX_REF_DIR" -maxdepth 1 -type d -iname "ldblk*${low}*" -print -quit 2>/dev/null || true)
    [[ -n $d ]] && record PASS "PRSCSx_HDF5_LD:$p" "$d" || record FAIL "PRSCSx_HDF5_LD:$p" "$GRID_CSX_REF_DIR"
  done
  for c in $(grid_expand_chrs "$GRID_CHRS"); do
    f="$GRID_ARG_TREES_DIR/chr$c.variants.tsv.gz"
    [[ -s $f ]] && record PASS "ARG_features_chr$c" "$f" || record WARN "ARG_features_chr$c" 'created by arg.sh build --method needle'
  done
fi

if ((fail)); then
  echo "PREFLIGHT FAILED: $report" >&2
  exit 1
fi
echo "PREFLIGHT PASS (warnings may remain): $report"
