#!/usr/bin/env bash
set -euo pipefail

# ARG construction belongs to reference-data preparation. This private helper
# lets TRACE validate refGen outputs before starting its analysis scheduler.
ROOT=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd -P)
REFGEN_SH=${GU_REFGEN_SH:-$ROOT/../0data/refGen.sh}
ACTION=${1:-check}
[[ $ACTION == check ]] || { echo "ERROR: internal ARG validation is check-only; build ARG data with /mnt/d/scripts/0data/refGen.sh make-arg" >&2; exit 2; }
[[ -s $REFGEN_SH ]] || { echo "ERROR: refGen.sh is missing: $REFGEN_SH" >&2; exit 2; }

backend=${GU_ARG_METHOD:-${GU_ARG_BACKEND:-tsinfer}}
case "${backend,,}" in
  tsinfer) method=tsinfer; ref_action=check;;
  external) method=tsinfer; ref_action=check;;
  needle) method=needle; ref_action=check;;
  *) echo "ERROR: GU ARG method must be needle or tsinfer" >&2; exit 2;;
esac
format=${GU_ARG_FORMAT:-native}
case "$format" in native|trace) ;; *) echo "ERROR: GU_ARG_FORMAT must be native or trace" >&2; exit 2;; esac
replace=FALSE
chrs=${GU_CHRS:-1-22,X}

args=(make-arg --method "$method" --format "$format" --action "$ref_action"
      --dir-gen "${GU_TARGET_ROOT:?}" --arg-dir "${GU_ARG_DIR:-$GU_TARGET_ROOT/arg}"
      --chr "$chrs" --grch "${GU_BUILD:-37}" --replace "$replace"
      --threads "${GU_ARG_THREADS:-8}" --jobs "${GU_UNIT_JOBS:-1}"
      --sample-match-min-work "${GU_ARG_SAMPLE_MATCH_MIN_WORK:-50000000}")
[[ -z ${GU_TARGET_VCF_DIR:-} ]] || args+=(--dir-vcf "$GU_TARGET_VCF_DIR")
[[ -z ${GU_SAMPLE_PANEL:-} ]] || args+=(--sample-file "$GU_SAMPLE_PANEL")
if [[ $method == needle && -n ${GU_TARGET_GEN_PREFIX:-} ]]; then
  args+=(--dir-pfile "$(dirname -- "$GU_TARGET_GEN_PREFIX")")
fi
export REFGEN_RECOMB_RATE=${GU_RECOMB_RATE:-1e-8}
export REFGEN_MUT_RATE=${GU_MUT_RATE:-1.25e-8}
export REFGEN_ARG_ONE_BIT=${GU_ARG_ONE_BIT:-0}

echo "[GU CHECK] internal read-only ARG validation through:"
printf '  bash %q' "$REFGEN_SH"; printf ' %q' "${args[@]}"; printf '\n'
exec bash "$REFGEN_SH" "${args[@]}"
