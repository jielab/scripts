#!/usr/bin/env bash
set -euo pipefail
export PYTHONDONTWRITEBYTECODE=1
ROOT=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
F=$ROOT/f
GU_ORIGINAL_ARGS=("$@")
# shellcheck source=f/regions.sh
source "$F/regions.sh"
# shellcheck source=f/target.sh
source "$F/target.sh"
[[ -s "$ROOT/gu.env" ]] && source "$ROOT/gu.env"


# 🚩 Driver and CLI
# Parse the complete driver before starting any long-running child process.
# Keeping both the invocation and the final exit in one parsed command group
# prevents a live process from reading newly appended/rearranged script bytes
# if gu.sh is updated in place while an analysis is still running.
gu_main(){
usage(){ cat <<'HELP'
GU archaic introgression workflow

Usage:
  ./gu.sh phyml [check|match] [--loci FILE | --chr LIST] [--grch 37|38] [options]
  ./gu.sh ibdmix [check] [--loci FILE | --chr LIST] [--grch 37|38] [options]
  ./gu.sh trace [check|extract|infer|segments] [--loci FILE | --chr LIST] [--grch 37|38] [options]
  ./gu.sh as3 [check] [--loci FILE | --chr LIST] [--grch 38] [options]
  ./gu.sh arg build [--chr LIST] [--grch 37|38] [options]
  ./gu.sh normalize
  ./gu.sh shiny
  ./gu.sh ukb inspect-hap|make-panel|batches|hap-vcf|hap-arg-vcf|inspect-typed

Region selection:
  --loci FILE   BED: chr, 0-based start, half-open end, optional ID in column 4.
  --loci-flank SIZE  Search-interval flank [100kb]; accepts bp/kb/mb.
  --chr LIST    Comma- or space-separated chromosomes, for example 1,2,3,X.
  --grch BUILD  Genome build: 37 or 38.
  --target NAME Stable target dataset name used in output paths (for example 1kg or ukb).
  --target-dir PREFIX
                Per-chromosome VCF or PLINK2 pfile prefix; use together with --target.
                Whole/pure non-PAR X resolves to <PREFIX>X.male. An X locus whose
                analysis interval is wholly inside PAR1/PAR2 resolves to <PREFIX>XY;
                missing chrXY, boundary-crossing, or mixed pure-X/PAR loci are errors.
                Compressed .pvar.zst is auto-detected and read with PLINK2 'vzs'.
  --auto-normalize TRUE|FALSE
                Rebuild Shiny/SQLite after analysis [FALSE].

Method options:
  phyml:  --plot-phy TRUE|FALSE --replace-phyml TRUE|FALSE
          --phyml-region-mode core|ld --phyml-ld-r2 FLOAT
          --phyml-window-bp INT --phyml-jobs INT
          Exact BED core is the default; ld enables anchor-SNP LD shrinkage.
  ibdmix: --replace-ibdmix TRUE|FALSE
  trace:  --replace-trace TRUE|FALSE --trace-loci-mode posthoc|extract
  as3:    --replace-as3 TRUE|FALSE --as3-target-chunk-size INT
          Target samples per in-memory AS3 chunk [256].
  arg:    --replace-arg TRUE|FALSE --arg-dir DIR
          --target-vcf-dir DIR
  ukb:    --ukb-hap-root DIR --ukb-ref-fasta FILE --sample-panel FILE
          --ukb-keep FILE --ukb-1kg-vcf-dir DIR
          --ukb-vcf-out DIR --ukb-arg-vcf-out DIR
          --ukb-batch-size INT --ukb-anchors-per-group INT

Command orchestration (phyml, ibdmix, trace, as3, arg):
  --run-cmd TRUE|FALSE
          Always write one permanent command per analysis unit. TRUE executes
          generated commands locally; FALSE only writes them [TRUE].
  --foreground TRUE|FALSE
          Keep the request scheduler attached to the terminal [FALSE]. By
          default it is launched with setsid/nohup and reports log/PID paths.
  --jobs INT
          Concurrent local command files when --run-cmd TRUE [4].

Examples:
  ./gu.sh phyml --loci /mnt/d/files/gu.37.bed --plot-phy TRUE --jobs 4
  ./gu.sh ibdmix --chr X --grch 37 --target 1kg --target-dir /mnt/d/data.BIG/refGen/1kg/37/pfile/chr --jobs 4
  ./gu.sh arg build --chr 22,X --grch 37 --target 1kg --target-dir /mnt/d/data.BIG/refGen/1kg/37/pfile/chr --jobs 4
  ./gu.sh trace --chr 22,X --grch 37 --target 1kg --target-dir /mnt/d/data.BIG/refGen/1kg/37/pfile/chr --jobs 4
  ./gu.sh as3 --chr 3,22 --grch 38 --target 1kg --target-dir /mnt/d/data.BIG/refGen/1kg/38/pfile/chr --jobs 4
  ./gu.sh normalize
  ./gu.sh shiny
HELP
}

case "${1:-help}" in help|-h|--help) usage; exit 0;; esac
METHOD=$1; shift
case "$METHOD" in phyml|ibdmix|trace|as3|arg|normalize|shiny|ukb) ;; *) echo "ERROR: unknown method: $METHOD" >&2; usage >&2; exit 2;; esac

ACTION=""
case "$METHOD" in
  phyml) valid_actions=' check match ';;
  ibdmix) valid_actions=' check ';;
  trace) valid_actions=' check extract infer segments ';;
  as3) valid_actions=' check ';;
  arg) valid_actions=' build ';;
  normalize) valid_actions=' ';;
  shiny) valid_actions=' ';;
  ukb) valid_actions=' inspect-hap make-panel batches hap-vcf hap-arg-vcf inspect-typed ';;
esac
if [[ $# -gt 0 && $1 != --* ]]; then
  [[ $valid_actions == *" $1 "* ]] || { echo "ERROR: invalid $METHOD action: $1" >&2; exit 2; }
  ACTION=$1; shift
fi

gu_bool() {
  case "${1,,}" in
    true|1|yes|y|on) printf 'TRUE\n' ;;
    false|0|no|n|off) printf 'FALSE\n' ;;
    *) return 2 ;;
  esac
}

LOCI_INPUT=""
LOCI_FLANK_INPUT=${GU_LOCI_FLANK:-100kb}
LOCI_FLANK_SET=0
CHR_INPUT=""
GRCH_INPUT=""
TARGET_INPUT=""
TARGET_DIR_INPUT=""
TARGET_INPUT_SET=0
TARGET_DIR_INPUT_SET=0
PLOT_PHY_INPUT=${PHYML_PLOT_PHY:-FALSE}
PLOT_PHY_SET=0
REPLACE_PHYML_INPUT=FALSE
REPLACE_PHYML_SET=0
REPLACE_IBDMIX_INPUT=FALSE
REPLACE_IBDMIX_SET=0
REPLACE_TRACE_INPUT=FALSE
REPLACE_TRACE_SET=0
REPLACE_AS3_INPUT=FALSE
REPLACE_AS3_SET=0
AS3_TARGET_CHUNK_SIZE_INPUT=${AS3_TARGET_CHUNK_SIZE:-64}
AS3_TARGET_CHUNK_SIZE_SET=0
REPLACE_ARG_INPUT=FALSE
REPLACE_ARG_SET=0
AUTO_NORMALIZE_INPUT=FALSE
TRACE_LOCI_MODE_INPUT=posthoc
TRACE_LOCI_MODE_SET=0
PHYML_WINDOW_BP_INPUT=500000
PHYML_WINDOW_BP_SET=0
PHYML_JOBS_INPUT=1
PHYML_JOBS_SET=0
PHYML_LD_R2_INPUT=${PHYML_LD_R2:-0.98}
PHYML_LD_R2_SET=0
PHYML_REGION_MODE_INPUT=${PHYML_REGION_MODE:-core}
PHYML_REGION_MODE_SET=0
RUN_CMD_INPUT=TRUE
RUN_CMD_SET=0
FOREGROUND_INPUT=FALSE
FOREGROUND_SET=0
UNIT_JOBS_INPUT=4
UNIT_JOBS_SET=0
ARG_DIR_INPUT=""
ARG_DIR_SET=0
TARGET_VCF_DIR_INPUT=""
TARGET_VCF_DIR_SET=0
SAMPLE_PANEL_INPUT=""
UKB_OPTION_SET=0
UKB_HAP_ROOT_INPUT=""
UKB_TYPED_ROOT_INPUT=""
UKB_WORK_INPUT=""
UKB_THREADS_INPUT=""
UKB_REF_FASTA_INPUT=""
UKB_KEEP_INPUT=""
UKB_1KG_VCF_DIR_INPUT=""
UKB_VCF_OUT_INPUT=""
UKB_ARG_VCF_OUT_INPUT=""
UKB_BATCH_SIZE_INPUT=""
UKB_ANCHORS_PER_GROUP_INPUT=""
UKB_ANCHOR_LIST_INPUT=""
UKB_PANEL_CHR_INPUT=""
EXTRA=()
while (( $# )); do
  case "$1" in
    --loci) [[ $# -ge 2 && -n ${2:-} ]] || { echo "ERROR: --loci requires FILE" >&2; exit 2; }; LOCI_INPUT=$2; shift 2;;
    --loci-flank) [[ $# -ge 2 && -n ${2:-} ]] || { echo "ERROR: --loci-flank requires SIZE" >&2; exit 2; }; LOCI_FLANK_INPUT=$2; LOCI_FLANK_SET=1; shift 2;;
    --chr) [[ $# -ge 2 && -n ${2:-} ]] || { echo "ERROR: --chr requires LIST" >&2; exit 2; }; CHR_INPUT=$2; shift 2;;
    --grch) [[ $# -ge 2 && -n ${2:-} ]] || { echo "ERROR: --grch requires 37 or 38" >&2; exit 2; }; GRCH_INPUT=$2; shift 2;;
    --target) [[ $# -ge 2 && -n ${2:-} ]] || { echo "ERROR: --target requires a dataset name" >&2; exit 2; }; TARGET_INPUT=$2; TARGET_INPUT_SET=1; shift 2;;
    --target-dir) [[ $# -ge 2 && -n ${2:-} ]] || { echo "ERROR: --target-dir requires a per-chromosome prefix" >&2; exit 2; }; TARGET_DIR_INPUT=$2; TARGET_DIR_INPUT_SET=1; shift 2;;
    --target-gen) echo "ERROR: --target-gen was replaced by --target NAME --target-dir PREFIX" >&2; exit 2;;
    --plot-phy) [[ $# -ge 2 && -n ${2:-} && $2 != --* ]] || { echo "ERROR: --plot-phy requires TRUE or FALSE" >&2; exit 2; }; PLOT_PHY_INPUT=$2; PLOT_PHY_SET=1; shift 2;;
    --replace-phyml) [[ $# -ge 2 && -n ${2:-} && $2 != --* ]] || { echo "ERROR: --replace-phyml requires TRUE or FALSE" >&2; exit 2; }; REPLACE_PHYML_INPUT=$2; REPLACE_PHYML_SET=1; shift 2;;
    --replace-ibdmix) [[ $# -ge 2 && -n ${2:-} && $2 != --* ]] || { echo "ERROR: --replace-ibdmix requires TRUE or FALSE" >&2; exit 2; }; REPLACE_IBDMIX_INPUT=$2; REPLACE_IBDMIX_SET=1; shift 2;;
    --replace-trace) [[ $# -ge 2 && -n ${2:-} && $2 != --* ]] || { echo "ERROR: --replace-trace requires TRUE or FALSE" >&2; exit 2; }; REPLACE_TRACE_INPUT=$2; REPLACE_TRACE_SET=1; shift 2;;
    --replace-as3) [[ $# -ge 2 && -n ${2:-} && $2 != --* ]] || { echo "ERROR: --replace-as3 requires TRUE or FALSE" >&2; exit 2; }; REPLACE_AS3_INPUT=$2; REPLACE_AS3_SET=1; shift 2;;
    --as3-target-chunk-size) [[ $# -ge 2 && -n ${2:-} ]] || { echo "ERROR: --as3-target-chunk-size requires a positive integer" >&2; exit 2; }; AS3_TARGET_CHUNK_SIZE_INPUT=$2; AS3_TARGET_CHUNK_SIZE_SET=1; shift 2;;
    --replace-arg) [[ $# -ge 2 && -n ${2:-} && $2 != --* ]] || { echo "ERROR: --replace-arg requires TRUE or FALSE" >&2; exit 2; }; REPLACE_ARG_INPUT=$2; REPLACE_ARG_SET=1; shift 2;;
    --auto-normalize) [[ $# -ge 2 && -n ${2:-} && $2 != --* ]] || { echo "ERROR: --auto-normalize requires TRUE or FALSE" >&2; exit 2; }; AUTO_NORMALIZE_INPUT=$2; shift 2;;
    --trace-loci-mode) [[ $# -ge 2 && -n ${2:-} ]] || { echo "ERROR: --trace-loci-mode requires posthoc or extract" >&2; exit 2; }; TRACE_LOCI_MODE_INPUT=$2; TRACE_LOCI_MODE_SET=1; shift 2;;
    --phyml-window-bp) [[ $# -ge 2 && -n ${2:-} ]] || { echo "ERROR: --phyml-window-bp requires a positive integer" >&2; exit 2; }; PHYML_WINDOW_BP_INPUT=$2; PHYML_WINDOW_BP_SET=1; shift 2;;
    --phyml-jobs) [[ $# -ge 2 && -n ${2:-} ]] || { echo "ERROR: --phyml-jobs requires a positive integer" >&2; exit 2; }; PHYML_JOBS_INPUT=$2; PHYML_JOBS_SET=1; shift 2;;
    --phyml-ld-r2) [[ $# -ge 2 && -n ${2:-} ]] || { echo "ERROR: --phyml-ld-r2 requires a value from 0 to 1" >&2; exit 2; }; PHYML_LD_R2_INPUT=$2; PHYML_LD_R2_SET=1; shift 2;;
    --phyml-region-mode) [[ $# -ge 2 && -n ${2:-} ]] || { echo "ERROR: --phyml-region-mode requires core or ld" >&2; exit 2; }; PHYML_REGION_MODE_INPUT=$2; PHYML_REGION_MODE_SET=1; shift 2;;
    --run-cmd) [[ $# -ge 2 && -n ${2:-} && $2 != --* ]] || { echo "ERROR: --run-cmd requires TRUE or FALSE" >&2; exit 2; }; RUN_CMD_INPUT=$2; RUN_CMD_SET=1; shift 2;;
    --foreground) [[ $# -ge 2 && -n ${2:-} && $2 != --* ]] || { echo "ERROR: --foreground requires TRUE or FALSE" >&2; exit 2; }; FOREGROUND_INPUT=$2; FOREGROUND_SET=1; shift 2;;
    --jobs) [[ $# -ge 2 && -n ${2:-} ]] || { echo "ERROR: --jobs requires a positive integer" >&2; exit 2; }; UNIT_JOBS_INPUT=$2; UNIT_JOBS_SET=1; shift 2;;
    --arg-dir) [[ $# -ge 2 && -n ${2:-} ]] || { echo "ERROR: --arg-dir requires DIR" >&2; exit 2; }; ARG_DIR_INPUT=$2; ARG_DIR_SET=1; shift 2;;
    --target-vcf-dir) [[ $# -ge 2 && -n ${2:-} ]] || { echo "ERROR: --target-vcf-dir requires DIR" >&2; exit 2; }; TARGET_VCF_DIR_INPUT=$2; TARGET_VCF_DIR_SET=1; shift 2;;
    --sample-panel) [[ $# -ge 2 && -n ${2:-} ]] || { echo "ERROR: --sample-panel requires FILE" >&2; exit 2; }; SAMPLE_PANEL_INPUT=$2; shift 2;;
    --ukb-hap-root) [[ $# -ge 2 && -n ${2:-} ]] || { echo "ERROR: --ukb-hap-root requires DIR" >&2; exit 2; }; UKB_HAP_ROOT_INPUT=$2; UKB_OPTION_SET=1; shift 2;;
    --ukb-typed-root) [[ $# -ge 2 && -n ${2:-} ]] || { echo "ERROR: --ukb-typed-root requires DIR" >&2; exit 2; }; UKB_TYPED_ROOT_INPUT=$2; UKB_OPTION_SET=1; shift 2;;
    --ukb-work) [[ $# -ge 2 && -n ${2:-} ]] || { echo "ERROR: --ukb-work requires DIR" >&2; exit 2; }; UKB_WORK_INPUT=$2; UKB_OPTION_SET=1; shift 2;;
    --ukb-threads) [[ $# -ge 2 && -n ${2:-} ]] || { echo "ERROR: --ukb-threads requires a positive integer" >&2; exit 2; }; UKB_THREADS_INPUT=$2; UKB_OPTION_SET=1; shift 2;;
    --ukb-ref-fasta) [[ $# -ge 2 && -n ${2:-} ]] || { echo "ERROR: --ukb-ref-fasta requires FILE" >&2; exit 2; }; UKB_REF_FASTA_INPUT=$2; UKB_OPTION_SET=1; shift 2;;
    --ukb-keep) [[ $# -ge 2 && -n ${2:-} ]] || { echo "ERROR: --ukb-keep requires FILE" >&2; exit 2; }; UKB_KEEP_INPUT=$2; UKB_OPTION_SET=1; shift 2;;
    --ukb-1kg-vcf-dir) [[ $# -ge 2 && -n ${2:-} ]] || { echo "ERROR: --ukb-1kg-vcf-dir requires DIR" >&2; exit 2; }; UKB_1KG_VCF_DIR_INPUT=$2; UKB_OPTION_SET=1; shift 2;;
    --ukb-vcf-out) [[ $# -ge 2 && -n ${2:-} ]] || { echo "ERROR: --ukb-vcf-out requires DIR" >&2; exit 2; }; UKB_VCF_OUT_INPUT=$2; UKB_OPTION_SET=1; shift 2;;
    --ukb-arg-vcf-out) [[ $# -ge 2 && -n ${2:-} ]] || { echo "ERROR: --ukb-arg-vcf-out requires DIR" >&2; exit 2; }; UKB_ARG_VCF_OUT_INPUT=$2; UKB_OPTION_SET=1; shift 2;;
    --ukb-batch-size) [[ $# -ge 2 && -n ${2:-} ]] || { echo "ERROR: --ukb-batch-size requires a positive integer" >&2; exit 2; }; UKB_BATCH_SIZE_INPUT=$2; UKB_OPTION_SET=1; shift 2;;
    --ukb-anchors-per-group) [[ $# -ge 2 && -n ${2:-} ]] || { echo "ERROR: --ukb-anchors-per-group requires a non-negative integer" >&2; exit 2; }; UKB_ANCHORS_PER_GROUP_INPUT=$2; UKB_OPTION_SET=1; shift 2;;
    --ukb-anchor-list) [[ $# -ge 2 && -n ${2:-} ]] || { echo "ERROR: --ukb-anchor-list requires FILE" >&2; exit 2; }; UKB_ANCHOR_LIST_INPUT=$2; UKB_OPTION_SET=1; shift 2;;
    --ukb-panel-chr) [[ $# -ge 2 && -n ${2:-} ]] || { echo "ERROR: --ukb-panel-chr requires CHR" >&2; exit 2; }; UKB_PANEL_CHR_INPUT=$2; UKB_OPTION_SET=1; shift 2;;
    *) EXTRA+=("$1"); shift;;
  esac
done
[[ -z $LOCI_INPUT || -z $CHR_INPUT ]] || { echo "ERROR: --loci and --chr are mutually exclusive" >&2; exit 2; }
[[ $LOCI_FLANK_SET == 0 || -n $LOCI_INPUT ]] || { echo "ERROR: --loci-flank requires --loci" >&2; exit 2; }
if (( PLOT_PHY_SET )) && [[ $METHOD != phyml || ${ACTION:-run} != run ]]; then
  echo "ERROR: --plot-phy is only valid with the default 'phyml' analysis" >&2
  exit 2
fi
if (( REPLACE_PHYML_SET )) && [[ $METHOD != phyml || ${ACTION:-run} != run ]]; then
  echo "ERROR: --replace-phyml is only valid with the default 'phyml' analysis" >&2
  exit 2
fi
if (( REPLACE_IBDMIX_SET )) && [[ $METHOD != ibdmix || ${ACTION:-run} != run ]]; then echo "ERROR: --replace-ibdmix is only valid with the default ibdmix analysis" >&2; exit 2; fi
if (( REPLACE_TRACE_SET )) && [[ $METHOD != trace ]]; then echo "ERROR: --replace-trace is only valid with trace" >&2; exit 2; fi
if (( REPLACE_AS3_SET )) && [[ $METHOD != as3 || ${ACTION:-run} != run ]]; then echo "ERROR: --replace-as3 is only valid with the default AS3 analysis" >&2; exit 2; fi
if (( AS3_TARGET_CHUNK_SIZE_SET )) && [[ $METHOD != as3 ]]; then echo "ERROR: --as3-target-chunk-size is only valid with as3" >&2; exit 2; fi
if (( REPLACE_ARG_SET )) && [[ $METHOD != arg ]]; then echo "ERROR: --replace-arg is only valid with arg" >&2; exit 2; fi
if (( ARG_DIR_SET )) && [[ $METHOD != arg && $METHOD != trace ]]; then echo "ERROR: --arg-dir is only valid with arg or trace" >&2; exit 2; fi
if (( TARGET_VCF_DIR_SET )) && [[ $METHOD != arg ]]; then echo "ERROR: --target-vcf-dir is only valid with arg" >&2; exit 2; fi
if (( UKB_OPTION_SET )) && [[ $METHOD != ukb ]]; then echo "ERROR: --ukb-* options are only valid with ukb" >&2; exit 2; fi
if (( TRACE_LOCI_MODE_SET )) && [[ $METHOD != trace ]]; then echo "ERROR: --trace-loci-mode is only valid with trace" >&2; exit 2; fi
if (( PHYML_WINDOW_BP_SET || PHYML_JOBS_SET || PHYML_LD_R2_SET || PHYML_REGION_MODE_SET )) && [[ $METHOD != phyml ]]; then echo "ERROR: --phyml-window-bp/--phyml-jobs/--phyml-ld-r2/--phyml-region-mode are only valid with phyml" >&2; exit 2; fi
if (( RUN_CMD_SET || FOREGROUND_SET || UNIT_JOBS_SET )); then
  command_action=${ACTION:-run}; [[ $METHOD != arg ]] || command_action=${ACTION:-build}
  case "$METHOD:$command_action" in
    phyml:run|ibdmix:run|trace:run|trace:extract|trace:infer|trace:segments|as3:run|arg:build) ;;
    *) echo "ERROR: --run-cmd/--foreground/--jobs are valid only with phyml, ibdmix, trace analysis actions, as3, or arg build" >&2; exit 2 ;;
  esac
fi
if (( TARGET_INPUT_SET != TARGET_DIR_INPUT_SET )); then
  echo "ERROR: --target and --target-dir must be supplied together" >&2
  exit 2
fi
if (( TARGET_INPUT_SET && TARGET_VCF_DIR_SET )); then
  echo "ERROR: use either --target/--target-dir or --target-vcf-dir for arg, not both" >&2
  exit 2
fi
if (( TARGET_INPUT_SET )); then
  case "$METHOD" in phyml|ibdmix|trace|as3|arg) ;; *) echo "ERROR: --target/--target-dir are only valid with phyml, ibdmix, trace, as3, or arg" >&2; exit 2;; esac
  [[ $TARGET_INPUT =~ ^[A-Za-z0-9][A-Za-z0-9._-]*$ && $TARGET_INPUT != "." && $TARGET_INPUT != ".." ]] || {
    echo "ERROR: --target must contain only letters, numbers, '.', '_' or '-', and must not be '.' or '..'" >&2
    exit 2
  }
fi
if [[ $METHOD == arg && -n $LOCI_INPUT ]]; then
  echo "ERROR: arg requires full chromosomes; --loci would create a partial ARG that TRACE cannot use. Use --chr LIST." >&2
  exit 2
fi

PLOT_PHY_INPUT=$(gu_bool "$PLOT_PHY_INPUT") || { echo "ERROR: --plot-phy must be TRUE or FALSE" >&2; exit 2; }
REPLACE_PHYML_INPUT=$(gu_bool "$REPLACE_PHYML_INPUT") || { echo "ERROR: --replace-phyml must be TRUE or FALSE" >&2; exit 2; }
REPLACE_IBDMIX_INPUT=$(gu_bool "$REPLACE_IBDMIX_INPUT") || { echo "ERROR: --replace-ibdmix must be TRUE or FALSE" >&2; exit 2; }
REPLACE_TRACE_INPUT=$(gu_bool "$REPLACE_TRACE_INPUT") || { echo "ERROR: --replace-trace must be TRUE or FALSE" >&2; exit 2; }
REPLACE_AS3_INPUT=$(gu_bool "$REPLACE_AS3_INPUT") || { echo "ERROR: --replace-as3 must be TRUE or FALSE" >&2; exit 2; }
REPLACE_ARG_INPUT=$(gu_bool "$REPLACE_ARG_INPUT") || { echo "ERROR: --replace-arg must be TRUE or FALSE" >&2; exit 2; }
AUTO_NORMALIZE_INPUT=$(gu_bool "$AUTO_NORMALIZE_INPUT") || { echo "ERROR: --auto-normalize must be TRUE or FALSE" >&2; exit 2; }
RUN_CMD_INPUT=$(gu_bool "$RUN_CMD_INPUT") || { echo "ERROR: --run-cmd must be TRUE or FALSE" >&2; exit 2; }
FOREGROUND_INPUT=$(gu_bool "$FOREGROUND_INPUT") || { echo "ERROR: --foreground must be TRUE or FALSE" >&2; exit 2; }
[[ $TRACE_LOCI_MODE_INPUT == posthoc || $TRACE_LOCI_MODE_INPUT == extract ]] || { echo "ERROR: --trace-loci-mode must be posthoc or extract" >&2; exit 2; }
[[ $PHYML_REGION_MODE_INPUT == core || $PHYML_REGION_MODE_INPUT == ld ]] || { echo "ERROR: --phyml-region-mode must be core or ld" >&2; exit 2; }
[[ $PHYML_WINDOW_BP_INPUT =~ ^[1-9][0-9]*$ ]] || { echo "ERROR: --phyml-window-bp must be a positive integer" >&2; exit 2; }
[[ $PHYML_JOBS_INPUT =~ ^[1-9][0-9]*$ ]] || { echo "ERROR: --phyml-jobs must be a positive integer" >&2; exit 2; }
[[ $UNIT_JOBS_INPUT =~ ^[1-9][0-9]*$ ]] || { echo "ERROR: --jobs must be a positive integer" >&2; exit 2; }
awk -v x="$PHYML_LD_R2_INPUT" 'BEGIN{exit !(x ~ /^([0-9]+([.][0-9]*)?|[.][0-9]+)$/ && x >= 0 && x <= 1)}' || { echo "ERROR: --phyml-ld-r2 must be between 0 and 1" >&2; exit 2; }
[[ $AS3_TARGET_CHUNK_SIZE_INPUT =~ ^[1-9][0-9]*$ ]] || { echo "ERROR: --as3-target-chunk-size must be a positive integer" >&2; exit 2; }
[[ -z $UKB_THREADS_INPUT || $UKB_THREADS_INPUT =~ ^[1-9][0-9]*$ ]] || { echo "ERROR: --ukb-threads must be a positive integer" >&2; exit 2; }
[[ -z $UKB_BATCH_SIZE_INPUT || $UKB_BATCH_SIZE_INPUT =~ ^[1-9][0-9]*$ ]] || { echo "ERROR: --ukb-batch-size must be a positive integer" >&2; exit 2; }
[[ -z $UKB_ANCHORS_PER_GROUP_INPUT || $UKB_ANCHORS_PER_GROUP_INPUT =~ ^[0-9]+$ ]] || { echo "ERROR: --ukb-anchors-per-group must be a non-negative integer" >&2; exit 2; }

GU_AUTO_NORMALIZE=0
[[ $AUTO_NORMALIZE_INPUT == FALSE ]] || GU_AUTO_NORMALIZE=1
GU_ARG_REPLACE=0; [[ $REPLACE_ARG_INPUT == FALSE ]] || GU_ARG_REPLACE=1
IBDMIX_REPLACE=0; [[ $REPLACE_IBDMIX_INPUT == FALSE ]] || IBDMIX_REPLACE=1
TRACE_REPLACE=0; [[ $REPLACE_TRACE_INPUT == FALSE ]] || TRACE_REPLACE=1
AS3_REPLACE=0; [[ $REPLACE_AS3_INPUT == FALSE ]] || AS3_REPLACE=1
AS3_TARGET_CHUNK_SIZE=$AS3_TARGET_CHUNK_SIZE_INPUT
TRACE_LOCI_MODE=$TRACE_LOCI_MODE_INPUT
PHYML_CHR_WINDOW_BP=$PHYML_WINDOW_BP_INPUT
PHYML_JOBS=$PHYML_JOBS_INPUT
PHYML_LD_R2=$PHYML_LD_R2_INPUT
PHYML_REGION_MODE=$PHYML_REGION_MODE_INPUT
GU_RUN_CMD=$RUN_CMD_INPUT
GU_FOREGROUND=$FOREGROUND_INPUT
GU_UNIT_JOBS=$UNIT_JOBS_INPUT
export GU_AUTO_NORMALIZE GU_ARG_REPLACE IBDMIX_REPLACE TRACE_REPLACE AS3_REPLACE AS3_TARGET_CHUNK_SIZE
export TRACE_LOCI_MODE PHYML_CHR_WINDOW_BP PHYML_JOBS PHYML_LD_R2 PHYML_REGION_MODE GU_RUN_CMD GU_FOREGROUND GU_UNIT_JOBS
if [[ $METHOD == phyml ]]; then
  PHYML_PLOT_PHY=$PLOT_PHY_INPUT
  PHYML_REPLACE=$REPLACE_PHYML_INPUT
  export PHYML_PLOT_PHY PHYML_REPLACE
fi
if (( ${#EXTRA[@]} )); then echo "ERROR: unknown $METHOD option: ${EXTRA[*]}" >&2; exit 2; fi


# 🚩 Runtime environment
# Re-enter the general GU environment before creating request-scoped temporary
# files.  Re-exec does not run EXIT traps, so doing this later would strand the
# first process's region workspace.
ENV_NAME=${GU_ENV_NAME:-gu}
if [[ ${CONDA_DEFAULT_ENV:-} != "$ENV_NAME" && ${GU_ENV_REEXEC:-0} != 1 ]]; then
  conda_owner=""
  for conda_exe in "${GU_CONDA_EXE:-}" "$HOME/anaconda3/bin/conda" "$HOME/miniconda3/bin/conda" /opt/conda/bin/conda "$(command -v conda 2>/dev/null || true)"; do
    [[ -n $conda_exe && -x $conda_exe ]] || continue
    if "$conda_exe" run -n "$ENV_NAME" true >/dev/null 2>&1; then conda_owner=$conda_exe; break; fi
  done
  [[ -n $conda_owner ]] || { echo "ERROR: Conda environment '$ENV_NAME' is not installed or discoverable; run ./install.sh" >&2; exit 1; }
  export GU_ENV_REEXEC=1
  exec "$conda_owner" run --no-capture-output -n "$ENV_NAME" bash "$ROOT/gu.sh" "${GU_ORIGINAL_ARGS[@]}"
fi

GU_DATA_ROOT=${GU_DATA_ROOT:-/mnt/d}
GU_REF_ROOT=${GU_REF_ROOT:-$GU_DATA_ROOT/data.BIG/refGen}
GU_ANALYSIS_ROOT=${GU_ANALYSIS_ROOT:-$GU_DATA_ROOT/analysis/gu}
GU_RUN_TMP_ROOT=${GU_RUN_TMP_ROOT:-$GU_DATA_ROOT/tmp/gu}
GU_SOFT=${GU_SOFT:-$GU_DATA_ROOT/software/gu}
default_build=37; [[ $METHOD == as3 ]] && default_build=38
GU_BUILD=${GRCH_INPUT:-${GU_BUILD:-$default_build}}
case "${GU_BUILD,,}" in
  37|b37|grch37|hg19) GU_BUILD=37;;
  38|b38|grch38|hg38) GU_BUILD=38;;
  *) echo "ERROR: --grch must be 37 or 38" >&2; exit 2;;
esac
AS3_EXPERIMENTAL_BUILD=0
if [[ $METHOD == as3 && $GU_BUILD != 38 ]]; then
  echo "WARNING: AS3 stopped before preprocessing: the public data support GRCh38 or CHM13, and this GU entry supports GRCh38 only (requested GRCh$GU_BUILD)." >&2
  exit 2
fi
GU_TARGET=${TARGET_INPUT:-${GU_TARGET:-1kg}}
[[ $GU_TARGET =~ ^[A-Za-z0-9][A-Za-z0-9._-]*$ && $GU_TARGET != "." && $GU_TARGET != ".." ]] || {
  echo "ERROR: target name must contain only letters, numbers, '.', '_' or '-', and must not be '.' or '..'" >&2
  exit 2
}
if [[ -n $TARGET_DIR_INPUT ]]; then
  GU_TARGET_GEN_PREFIX=$TARGET_DIR_INPUT
  GU_TARGET_ROOT=${GU_TARGET_ROOT:-$(dirname -- "$(dirname -- "$GU_TARGET_GEN_PREFIX")")}
else
  GU_TARGET_ROOT=${GU_TARGET_ROOT:-$GU_REF_ROOT/$GU_TARGET/$GU_BUILD}
  GU_TARGET_GEN_PREFIX=${GU_TARGET_DIR:-${GU_TARGET_GEN_PREFIX:-$GU_TARGET_ROOT/pfile/chr}}
fi


# 🚩 Target inputs
# PFILE is the primary target format. Methods which require VCF receive a
# single-run VCF under the analysis temporary directory.  AS3 is VCF-native, so it directly
# reuses indexed chrN.vcf.gz inputs instead of needlessly round-tripping them
# through PGEN.  The explicit chrX.male.vcf.gz sibling remains available as a
# validated VCF input when a chrX.male pfile is absent.
GU_TARGET_NATIVE_VCF_PREFIX=${GU_TARGET_NATIVE_VCF_PREFIX-$GU_TARGET_ROOT/vcf/chr}
GU_TARGET_DIR=$GU_TARGET_GEN_PREFIX
if [[ $METHOD == as3 ]]; then
  AS3_RUNTIME=$F/as3_upstream
  AS3_REFERENCE_PANEL_DIR=${AS3_REFERENCE_PANEL_DIR:-$GU_REF_ROOT/archaic/38/vcf}
  AS3_REFERENCE_MAP=${AS3_REFERENCE_MAP:-$AS3_REFERENCE_PANEL_DIR/Ref_Panel.map.txt}
  AS3_MASK_DIR=${AS3_MASK_DIR:-$GU_REF_ROOT/archaic/38/mask}
  AS3_MODEL_DIR=${AS3_MODEL_DIR:-$GU_REF_ROOT/archaic/38/models}
  GU_ARCHAIC_ROOT=${GU_ARCHAIC_ROOT:-$AS3_REFERENCE_PANEL_DIR}
else
  archaic_build_root_var=GU_ARCHAIC${GU_BUILD}_ROOT
  GU_ARCHAIC_ROOT=${GU_ARCHAIC_ROOT:-${!archaic_build_root_var:-$GU_REF_ROOT/archaic/$GU_BUILD/vcf}}
fi
AS3_PUBLISHED_CALLS_DIR=${AS3_PUBLISHED_CALLS_DIR:-$GU_REF_ROOT/archaic/38/database_introgression_call}
AS3_REFERENCE_CALLSET_CACHE=${AS3_REFERENCE_CALLSET_CACHE:-$GU_ANALYSIS_ROOT/normalize/reference-callsets.sqlite}
GU_SAMPLE_PANEL=${SAMPLE_PANEL_INPUT:-${GU_SAMPLE_PANEL:-}}
GU_ARG_DIR=${ARG_DIR_INPUT:-${GU_ARG_DIR:-$GU_TARGET_ROOT/arg}}
if [[ $METHOD == phyml && -z $LOCI_INPUT && -z $CHR_INPUT ]]; then
  LOCI_INPUT=${GU_DEFAULT_LOCI:-/mnt/d/files/gu.$GU_BUILD.bed}
fi

GU_LOCI_FILE=""; GU_LOCI_CORE_FILE=""; GU_LOCI_MAP_FILE=""; GU_LOCI_FLANK_BP=0
GU_REGION_TMP_ROOT=""
GU_TARGET_TMP_DIR=""
gu_cleanup_region_tmp(){
  [[ -n ${GU_REGION_TMP_ROOT:-} && -d $GU_REGION_TMP_ROOT ]] || return 0
  case "$GU_REGION_TMP_ROOT" in
    "${TMPDIR:-/tmp}"/gu-regions.*) rm -rf -- "$GU_REGION_TMP_ROOT" ;;
    *) echo "WARNING: refusing to remove unexpected GU region temp path: $GU_REGION_TMP_ROOT" >&2 ;;
  esac
}
gu_cleanup_target_tmp(){
  [[ -n ${GU_TARGET_TMP_DIR:-} && -d $GU_TARGET_TMP_DIR ]] || return 0
  local tmp_base=$GU_RUN_TMP_ROOT parent
  case "$GU_TARGET_TMP_DIR" in
    "$tmp_base"/*/*/run.*)
      parent=$(dirname -- "$GU_TARGET_TMP_DIR")
      rm -rf -- "$GU_TARGET_TMP_DIR"
      while [[ $parent == "$tmp_base"/* ]]; do
        rmdir -- "$parent" 2>/dev/null || break
        parent=$(dirname -- "$parent")
      done
      rmdir -- "$tmp_base" 2>/dev/null || true
      ;;
    *) echo "WARNING: refusing to remove unexpected GU target temp path: $GU_TARGET_TMP_DIR" >&2 ;;
  esac
}
gu_cleanup_all_tmp(){
  gu_cleanup_region_tmp
  gu_cleanup_target_tmp
}
gu_cleanup_on_signal(){
  local signal=${1:?} code=1
  case "$signal" in HUP) code=129 ;; INT) code=130 ;; TERM) code=143 ;; esac
  trap - EXIT HUP INT TERM
  gu_cleanup_all_tmp
  exit "$code"
}
trap gu_cleanup_all_tmp EXIT
trap 'gu_cleanup_on_signal HUP' HUP
trap 'gu_cleanup_on_signal INT' INT
trap 'gu_cleanup_on_signal TERM' TERM
if [[ -n $LOCI_INPUT ]]; then
  GU_REGION_TMP_ROOT=$(mktemp -d "${TMPDIR:-/tmp}/gu-regions.XXXXXX")
  tmp_core="$GU_REGION_TMP_ROOT/.core.bed"; tmp_loci="$GU_REGION_TMP_ROOT/.analysis.bed"; tmp_map="$GU_REGION_TMP_ROOT/.map.tsv"
  gu_normalize_loci "$LOCI_INPUT" "$tmp_core"
  GU_LOCI_FLANK_BP=$(python3 "$F/expand_loci.py" --input "$tmp_core" --output "$tmp_loci" --map "$tmp_map" --flank "$LOCI_FLANK_INPUT" --build "$GU_BUILD")
  scope_key=$(gu_scope_id "$tmp_loci" "" "$LOCI_INPUT")
  scope_label=$(gu_scope_label "$tmp_loci" "" "$LOCI_INPUT")
  GU_LOCI_FILE="$GU_REGION_TMP_ROOT/$scope_key.bed"
  GU_LOCI_CORE_FILE="$GU_REGION_TMP_ROOT/$scope_key.core.bed"
  GU_LOCI_MAP_FILE="$GU_REGION_TMP_ROOT/$scope_key.map.tsv"
  mv -f "$tmp_loci" "$GU_LOCI_FILE"
  mv -f "$tmp_core" "$GU_LOCI_CORE_FILE"
  mv -f "$tmp_map" "$GU_LOCI_MAP_FILE"
  GU_CHRS=$(gu_loci_chrs "$GU_LOCI_FILE")
elif [[ -n $CHR_INPUT ]]; then
  GU_CHRS=$(gu_normalize_chrs "$CHR_INPUT")
  [[ -n $GU_CHRS ]] || { echo "ERROR: --chr produced an empty chromosome list" >&2; exit 2; }
  scope_key=$(gu_scope_id "" "$GU_CHRS")
  scope_label=$(gu_scope_label "" "$GU_CHRS")
else
  GU_CHRS=${GU_CHRS:-}
  [[ -n $GU_CHRS ]] && GU_CHRS=$(gu_normalize_chrs "$GU_CHRS")
  scope_key=$(gu_scope_id "" "$GU_CHRS")
  scope_label=$(gu_scope_label "" "$GU_CHRS")
fi
if [[ $METHOD == as3 && " $GU_CHRS " == *" X "* ]]; then
  echo "ERROR: AS3 supports GRCh38 autosomes 1-22 only; chrX is unavailable" >&2
  exit 2
fi
if [[ $METHOD == trace && ( ${ACTION:-run} == check || ${GU_CMD_WORKER:-0} == 1 ) ]]; then
  trace_precheck_chrs=${GU_CHRS:-"1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X"}
  [[ -d $GU_ARG_DIR ]] || { echo "ERROR: TRACE ARG directory is missing: $GU_ARG_DIR; run ./gu.sh arg build first" >&2; exit 2; }
  for c in $trace_precheck_chrs; do
    if ! find "$GU_ARG_DIR" -type f \( -name '*.trees' -o -name '*.tsz' \) -print 2>/dev/null |
      awk -v c="$c" 'BEGIN{IGNORECASE=1;found=0}{b=$0;gsub(/.*\//,"",b);if(b~("(^|[^A-Za-z0-9])chr"c"([^A-Za-z0-9]|$)")||b~("(^|[^A-Za-z0-9])"c"([^A-Za-z0-9]|$)"))found=1}END{exit !found}'; then
      echo "ERROR: TRACE requires a full-chromosome ARG for chr$c under $GU_ARG_DIR; run ./gu.sh arg build --chr $c --grch $GU_BUILD first" >&2
      exit 2
    fi
  done
fi

GU_SCOPE_KEY=$scope_key
GU_SCOPE_ID=$scope_label
GU_SCOPE_LABEL=$scope_label
GU_CHRX_TARGET_MODE=none
GU_CHRXY_EXTRACT_BED=""
if [[ " $GU_CHRS " == *" X "* ]]; then
  if [[ -n $GU_LOCI_FILE ]]; then
    GU_CHRX_TARGET_MODE=$(gu_classify_chrx_loci "$GU_LOCI_FILE" "$GU_BUILD") || exit 2
  else
    GU_CHRX_TARGET_MODE=male
  fi
fi
GU_CHRX_MALE_ONLY=0; GU_CHRX_PAR_DIPLOID=0
[[ $GU_CHRX_TARGET_MODE != male ]] || GU_CHRX_MALE_ONLY=1
[[ $GU_CHRX_TARGET_MODE != par ]] || GU_CHRX_PAR_DIPLOID=1
if (( GU_CHRX_PAR_DIPLOID )); then
  [[ -n $GU_REGION_TMP_ROOT && -n $GU_LOCI_FILE ]] || { echo "ERROR: chrXY routing is valid only for a PAR-scoped --loci request" >&2; exit 2; }
  GU_CHRXY_EXTRACT_BED=$GU_REGION_TMP_ROOT/${scope_key}.chrxy-extract.bed
  gu_make_chrxy_extract_bed "$GU_LOCI_FILE" "$GU_CHRXY_EXTRACT_BED" "$GU_BUILD" || exit 2
fi
GU_TARGET_NAMESPACE=$GU_TARGET


# 🚩 Command orchestration
# Command orchestration deliberately sits outside every method implementation.
# Each worker re-enters gu.sh with exactly one native analysis unit, so the
# method scripts and their scientific/output contracts remain unchanged.
gu_command_action(){
  if [[ $METHOD == arg ]]; then printf '%s\n' "${ACTION:-build}"; else printf '%s\n' "${ACTION:-run}"; fi
}

gu_command_orchestration_enabled(){
  local action; action=$(gu_command_action)
  case "$METHOD:$action" in
    phyml:run|ibdmix:run|trace:run|trace:extract|trace:infer|trace:segments|as3:run|arg:build) return 0 ;;
    *) return 1 ;;
  esac
}

gu_command_output(){
  local unit_label=$1 request_units=$2 override=""
  if [[ $METHOD == arg ]]; then printf '%s\n' "$GU_ARG_DIR"; return; fi
  case "$METHOD" in
    phyml) override=${PHYML_OUT:-} ;;
    ibdmix) override=${IBDMIX_OUT:-} ;;
    trace) override=${TRACE_OUT:-} ;;
    as3) override=${AS3_OUT:-} ;;
  esac
  gu_chr_result_dir "$GU_ANALYSIS_ROOT" "$METHOD" "$unit_label" "$GU_TARGET_NAMESPACE" "$override" "$request_units"
}

gu_write_analysis_unit_cmd(){
  local request_units=$1 unit_label=$2 kind=$3 chr=$4 core_start=${5:-} core_end=${6:-} locus=${7:-}
  local out cmd cmd_tmp unit_bed bed_tmp worker_action output_var=""
  local -a cmd_args
  out=$(gu_command_output "$unit_label" "$request_units")
  mkdir -p "$out"

  worker_action=$(gu_command_action)
  cmd_args=("$ROOT/gu.sh" "$METHOD")
  if [[ -n $ACTION || $METHOD == arg ]]; then cmd_args+=("$worker_action"); fi
  case "$kind" in
    locus)
      unit_bed=$out/$unit_label.bed
      bed_tmp=$unit_bed.tmp.$$
      printf '%s\t%s\t%s\t%s\n' "$chr" "$core_start" "$core_end" "$locus" > "$bed_tmp"
      ;;
    loci_chr)
      unit_bed=$out/${scope_label}.bed
      bed_tmp=$unit_bed.tmp.$$
      awk -F'\t' -v c="$chr" 'BEGIN{OFS="\t"}$1==c{print}' "$GU_LOCI_CORE_FILE" > "$bed_tmp"
      [[ -s $bed_tmp ]] || { rm -f "$bed_tmp"; echo "ERROR: no loci available for chr$chr" >&2; return 2; }
      ;;
    chromosome) unit_bed="" ;;
    *) echo "ERROR: internal command unit kind is invalid: $kind" >&2; return 2 ;;
  esac
  if [[ -n $unit_bed ]]; then
    if [[ -s $unit_bed ]] && cmp -s "$bed_tmp" "$unit_bed"; then rm -f "$bed_tmp"; else mv -f "$bed_tmp" "$unit_bed"; fi
    cmd_args+=(--loci "$unit_bed" --loci-flank "$LOCI_FLANK_INPUT")
  else
    cmd_args+=(--chr "$chr")
  fi
  cmd_args+=(--grch "$GU_BUILD")
  if [[ $METHOD == arg && $TARGET_VCF_DIR_SET == 1 ]]; then
    cmd_args+=(--target-vcf-dir "$TARGET_VCF_DIR_INPUT")
  else
    cmd_args+=(--target "$GU_TARGET" --target-dir "$GU_TARGET_GEN_PREFIX")
  fi
  case "$METHOD" in
    phyml)
      cmd_args+=(--plot-phy "$PLOT_PHY_INPUT" --replace-phyml "$REPLACE_PHYML_INPUT"
                 --phyml-ld-r2 "$PHYML_LD_R2_INPUT" --phyml-window-bp "$PHYML_WINDOW_BP_INPUT"
                 --phyml-jobs "$PHYML_JOBS_INPUT" --phyml-region-mode "$PHYML_REGION_MODE_INPUT")
      output_var=PHYML_OUT
      ;;
    ibdmix) cmd_args+=(--replace-ibdmix "$REPLACE_IBDMIX_INPUT"); output_var=IBDMIX_OUT ;;
    trace)
      cmd_args+=(--replace-trace "$REPLACE_TRACE_INPUT" --trace-loci-mode "$TRACE_LOCI_MODE_INPUT" --arg-dir "$GU_ARG_DIR")
      output_var=TRACE_OUT
      ;;
    as3) cmd_args+=(--replace-as3 "$REPLACE_AS3_INPUT" --as3-target-chunk-size "$AS3_TARGET_CHUNK_SIZE_INPUT"); output_var=AS3_OUT ;;
    arg) cmd_args+=(--replace-arg "$REPLACE_ARG_INPUT" --arg-dir "$GU_ARG_DIR") ;;
  esac
  # A unit command is already the leaf worker. Keep it attached to the local
  # scheduler (or an HPC scheduler) instead of recursively launching itself.
  cmd_args+=(--foreground TRUE)
  [[ $METHOD == arg ]] || cmd_args+=(--auto-normalize FALSE)
  [[ -z ${GU_SAMPLE_PANEL:-} ]] || cmd_args+=(--sample-panel "$GU_SAMPLE_PANEL")

  cmd=$out/$unit_label.cmd
  cmd_tmp=$cmd.tmp.$$
  {
    printf '#!/usr/bin/env bash\nset -euo pipefail\n'
    printf 'export GU_CMD_WORKER=1\n'
    printf 'export GU_ANALYSIS_ROOT=%q\n' "$GU_ANALYSIS_ROOT"
    [[ -z $output_var ]] || printf 'export %s=%q\n' "$output_var" "$out"
    printf 'exec'
    printf ' %q' "${cmd_args[@]}"
    printf '\n'
  } > "$cmd_tmp"
  chmod +x "$cmd_tmp"
  if [[ -s $cmd ]] && cmp -s "$cmd_tmp" "$cmd"; then rm -f "$cmd_tmp"; else mv -f "$cmd_tmp" "$cmd"; fi
  printf '%s\n' "$cmd"
}

gu_completed_analysis_cmd_output(){
  local out=$1 final
  case "$METHOD" in
    ibdmix)
      (( IBDMIX_REPLACE == 0 )) || return 1
      final=$out/final/all_archaic_refs.lod${IBDMIX_LOD:-4}.len${IBDMIX_MIN_BP:-50000}.segments.tsv.gz
      [[ -e $out/.complete && -s $final ]] || return 1
      printf '%s\n' "$final"
      ;;
    *) return 1 ;;
  esac
}

gu_run_one_analysis_cmd(){
  local cmd=$1 out base log err rc completed_output=""
  out=$(dirname -- "$cmd")
  base=$(basename -- "$cmd" .cmd)
  log=$out/$base.log
  err=$out/$base.err
  rm -f "$err"
  if completed_output=$(gu_completed_analysis_cmd_output "$out"); then
    printf '[GU CMD] SKIP unit=%s reason=output_complete final=%s\n' "$base" "$completed_output" >&2
    return 0
  fi
  printf '[GU CMD] START unit=%s cmd=%s\n' "$base" "$cmd" >&2
  if bash "$cmd" > "$log" 2>&1; then
    printf '[GU CMD] DONE unit=%s log=%s\n' "$base" "$log" >&2
    return 0
  else
    rc=$?
  fi
  printf 'ERROR: %s unit %s failed with exit=%s\nERROR: log=%s\n' "$METHOD" "$base" "$rc" "$log" > "$err"
  printf '[GU CMD] FAIL unit=%s exit=%s detail=%s\n' "$base" "$rc" "$err" >&2
  return "$rc"
}

gu_run_analysis_cmds_local(){
  local list=$1 jobs=$2 cmd running=0 status=0
  while IFS= read -r cmd; do
    [[ -s $cmd ]] || continue
    gu_run_one_analysis_cmd "$cmd" &
    running=$((running+1))
    if (( running >= jobs )); then
      if ! wait -n; then status=1; fi
      running=$((running-1))
    fi
  done < "$list"
  while (( running )); do
    if ! wait -n; then status=1; fi
    running=$((running-1))
  done
  (( status == 0 )) || { echo "ERROR: one or more $METHOD command files failed; see per-unit .err/.log files" >&2; return 1; }
}

gu_orchestrate_analysis_cmds(){
  local request_units cmd_root list list_tmp chr core_start core_end analysis_start analysis_end locus flank unit_label unit_kind
  local command_chrs=$GU_CHRS effective_jobs=$GU_UNIT_JOBS
  if [[ -z $GU_LOCI_MAP_FILE && -z $command_chrs ]]; then
    case "$METHOD" in
      as3) command_chrs="1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22" ;;
      ibdmix|trace|arg) command_chrs="1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X" ;;
    esac
  fi
  if [[ -n $GU_LOCI_MAP_FILE ]]; then
    case "$METHOD" in
      phyml|ibdmix) request_units=$(awk 'END{print (NR > 0 ? NR - 1 : 0)}' "$GU_LOCI_MAP_FILE"); unit_kind=locus ;;
      trace|as3) request_units=$(awk -F'\t' 'NR>1&&!seen[$1]++{n++}END{print n+0}' "$GU_LOCI_MAP_FILE"); unit_kind=loci_chr ;;
    esac
  else
    request_units=$(awk '{print NF}' <<< "$command_chrs"); unit_kind=chromosome
  fi
  (( request_units > 0 )) || { echo "ERROR: no $METHOD analysis units selected" >&2; return 2; }
  if [[ $METHOD == arg ]]; then cmd_root=$GU_ARG_DIR; else cmd_root=$GU_ANALYSIS_ROOT/$METHOD/$GU_TARGET_NAMESPACE; fi
  mkdir -p "$cmd_root"
  list=$cmd_root/${scope_label}.cmd.list
  list_tmp=$list.tmp.$$
  : > "$list_tmp"
  if [[ $unit_kind == locus ]]; then
    while IFS=$'\t' read -r chr core_start core_end analysis_start analysis_end locus flank; do
      unit_label=$(gu_locus_unit_label "$GU_LOCI_MAP_FILE" "$chr" "$core_start" "$core_end" "$locus")
      gu_write_analysis_unit_cmd "$request_units" "$unit_label" locus "$chr" "$core_start" "$core_end" "$locus" >> "$list_tmp"
    done < <(awk 'NR>1' "$GU_LOCI_MAP_FILE")
  elif [[ $unit_kind == loci_chr ]]; then
    while IFS= read -r chr; do
      unit_label=chr${chr}_${scope_label}
      gu_write_analysis_unit_cmd "$request_units" "$unit_label" loci_chr "$chr" >> "$list_tmp"
    done < <(awk -F'\t' 'NR>1&&!seen[$1]++{print $1}' "$GU_LOCI_MAP_FILE")
  else
    for chr in $command_chrs; do
      unit_label=chr$chr
      gu_write_analysis_unit_cmd "$request_units" "$unit_label" chromosome "$chr" >> "$list_tmp"
    done
  fi
  mv -f "$list_tmp" "$list"
  echo "[GU CMD] created=$request_units list=$list"
  if [[ $GU_RUN_CMD == FALSE ]]; then
    echo "[GU CMD] run_cmd=FALSE; generated command files only"
    return 0
  fi
  if [[ $METHOD == as3 && $effective_jobs -gt 1 ]]; then
    effective_jobs=1
    echo "[GU CMD] scheduler_cap method=as3 requested_jobs=$GU_UNIT_JOBS effective_jobs=$effective_jobs reason=shared_GPU_assignment"
  fi
  echo "[GU CMD] run_cmd=TRUE mode=local jobs=$effective_jobs"
  gu_run_analysis_cmds_local "$list" "$effective_jobs"
  if (( GU_AUTO_NORMALIZE )) && [[ $METHOD != arg ]]; then bash "$ROOT/gu.sh" normalize; fi
}

gu_launch_analysis_request_background(){
  local log_root stamp base log pid_file status_file pid
  local runner
  if [[ $METHOD == arg ]]; then
    log_root=$GU_ARG_DIR/log
  else
    log_root=$GU_ANALYSIS_ROOT/$METHOD/$GU_TARGET_NAMESPACE/log
  fi
  mkdir -p "$log_root"
  stamp=$(date '+%Y%m%d-%H%M%S')
  # METHOD-specific roots allow different modules to run concurrently.  The
  # launcher PID suffix also prevents collisions between two identical requests
  # submitted during the same second.
  base=$log_root/${METHOD}.${scope_label}.${stamp}.$$.background
  log=$base.log
  pid_file=$base.pid
  status_file=$base.status
  runner='status_file=$1; shift
tmp_status=${status_file}.tmp.$$
printf "state\trunning\npid\t%s\nstarted_at\t%s\n" "$$" "$(date -Is)" > "$tmp_status"
mv -f -- "$tmp_status" "$status_file"
"$@"
rc=$?
tmp_status=${status_file}.tmp.$$
printf "state\tfinished\nexit_code\t%s\nfinished_at\t%s\n" "$rc" "$(date -Is)" > "$tmp_status"
mv -f -- "$tmp_status" "$status_file"
exit "$rc"'
  if command -v setsid >/dev/null 2>&1; then
    nohup setsid bash -c "$runner" gu-background "$status_file" \
      bash "$ROOT/gu.sh" "${GU_ORIGINAL_ARGS[@]}" --foreground TRUE \
      </dev/null >"$log" 2>&1 &
  else
    nohup bash -c "$runner" gu-background "$status_file" \
      bash "$ROOT/gu.sh" "${GU_ORIGINAL_ARGS[@]}" --foreground TRUE \
      </dev/null >"$log" 2>&1 &
  fi
  pid=$!
  printf '%s\n' "$pid" > "$pid_file"
  disown "$pid" 2>/dev/null || true
  echo "[GU BG] STARTED method=$METHOD scope=$scope_label pid=$pid"
  echo "[GU BG] log=$log"
  echo "[GU BG] pid_file=$pid_file"
  echo "[GU BG] status_file=$status_file"
}

if gu_command_orchestration_enabled && [[ ${GU_CMD_WORKER:-0} != 1 && $GU_RUN_CMD == TRUE && $GU_FOREGROUND == FALSE ]]; then
  gu_launch_analysis_request_background
  exit $?
fi

if gu_command_orchestration_enabled && [[ ${GU_CMD_WORKER:-0} != 1 ]]; then
  gu_orchestrate_analysis_cmds
  exit $?
fi

USES_TARGET=0
case "$METHOD" in
  phyml|ibdmix|trace|as3) USES_TARGET=1 ;;
  arg) (( TARGET_VCF_DIR_SET )) || USES_TARGET=1 ;;
esac

declare -A GU_TARGET_SOURCE=() GU_TARGET_INDEX=() GU_TARGET_PVAR=() GU_TARGET_PSAM=() GU_TARGET_FORMAT=()
GU_TARGET_GEN_FORMAT=""
GU_TARGET_BUILD_SOURCE=""
GU_TARGET_BUILD_FORMAT=""
GU_CHECK_LOG=""


# 🚩 Target validation
gu_check_log(){ printf '[GU CHECK] %s\n' "$*" | tee -a "$GU_CHECK_LOG"; }
gu_check_fail(){ gu_check_log "ERROR: $*" >&2; exit 2; }
gu_detect_target_chr(){
  local c=$1 base native_base native_provenance_ok=0 pvar=""
  base=$(gu_target_genotype_base "$GU_TARGET_GEN_PREFIX" "$c")
  if [[ $METHOD != as3 && -s $base.pgen && -s $base.psam ]]; then
    for pvar in "$base.pvar" "$base.pvar.zst"; do [[ -s $pvar ]] && break; done
    [[ -s $pvar ]] || return 4
    GU_DETECT_FORMAT=pfile; GU_DETECT_SOURCE=$base; GU_DETECT_INDEX=""
    GU_DETECT_PVAR=$pvar; GU_DETECT_PSAM=$base.psam
    return 0
  fi
  if [[ $c == X && ${GU_CHRX_TARGET_MODE:-male} == male && -n ${GU_TARGET_NATIVE_VCF_PREFIX:-} ]]; then
    native_base=${GU_TARGET_NATIVE_VCF_PREFIX}X.male
    if [[ -s $native_base.vcf.gz ]]; then
      if [[ -s $native_base.vcf.gz.tbi ]]; then GU_DETECT_INDEX=$native_base.vcf.gz.tbi
      elif [[ -s $native_base.vcf.gz.csi ]]; then GU_DETECT_INDEX=$native_base.vcf.gz.csi
      else return 3
      fi
      GU_DETECT_FORMAT=vcf; GU_DETECT_SOURCE=$native_base.vcf.gz
      GU_DETECT_PVAR=""; GU_DETECT_PSAM=""
      return 0
    fi
  fi
  if [[ $c != X && -n ${GU_TARGET_NATIVE_VCF_PREFIX:-} ]]; then
    native_base=$GU_TARGET_NATIVE_VCF_PREFIX$c
    if [[ -s $native_base.vcf.gz ]]; then
      if [[ -s $native_base.vcf.gz.tbi ]]; then GU_DETECT_INDEX=$native_base.vcf.gz.tbi
      elif [[ -s $native_base.vcf.gz.csi ]]; then GU_DETECT_INDEX=$native_base.vcf.gz.csi
      else native_base=""
      fi
      # When a pfile is also present, require its PLINK log to identify this
      # exact VCF as the import source. This avoids silently substituting an
      # unrelated sibling VCF in a custom target layout.
      if [[ -n $native_base ]]; then
        if [[ ! -s $base.pgen ]]; then native_provenance_ok=1
        elif [[ -s $base.log ]] && grep -Fqx "  --vcf $native_base.vcf.gz" "$base.log"; then
          native_provenance_ok=1
          for pvar in "$base.pvar" "$base.pvar.zst"; do [[ -s $pvar ]] && break; done
          [[ -s $pvar ]] || pvar=""
        elif [[ ${GU_TRUST_NATIVE_TARGET_VCF:-0} == 1 ]]; then native_provenance_ok=1
        fi
      fi
      if (( native_provenance_ok )); then
        GU_DETECT_FORMAT=vcf; GU_DETECT_SOURCE=$native_base.vcf.gz; GU_DETECT_PVAR=$pvar
        [[ -z $pvar ]] && GU_DETECT_PSAM="" || GU_DETECT_PSAM=$base.psam
        return 0
      fi
    fi
  fi
  if [[ $c == X && ${GU_CHRX_TARGET_MODE:-male} == par && -s $base.pgen && -s $base.psam ]]; then
    for pvar in "$base.pvar" "$base.pvar.zst"; do [[ -s $pvar ]] && break; done
    [[ -s $pvar ]] || return 4
    GU_DETECT_FORMAT=pfile; GU_DETECT_SOURCE=$base; GU_DETECT_INDEX=""; GU_DETECT_PVAR=$pvar; GU_DETECT_PSAM=$base.psam; return 0
  fi
  if [[ -s $base.vcf.gz ]]; then
    if [[ -s $base.vcf.gz.tbi ]]; then GU_DETECT_INDEX=$base.vcf.gz.tbi
    elif [[ -s $base.vcf.gz.csi ]]; then GU_DETECT_INDEX=$base.vcf.gz.csi
    else return 3
    fi
    GU_DETECT_FORMAT=vcf; GU_DETECT_SOURCE=$base.vcf.gz; GU_DETECT_PVAR=""; GU_DETECT_PSAM=""; return 0
  fi
  if [[ -s $base.pgen && -s $base.psam ]]; then
    for pvar in "$base.pvar" "$base.pvar.zst"; do [[ -s $pvar ]] && break; done
    [[ -s $pvar ]] || return 4
    GU_DETECT_FORMAT=pfile; GU_DETECT_SOURCE=$base; GU_DETECT_INDEX=""; GU_DETECT_PVAR=$pvar; GU_DETECT_PSAM=$base.psam; return 0
  fi
  return 1
}

if (( USES_TARGET )); then
  target_tmp_parent=$GU_RUN_TMP_ROOT/$GU_TARGET_NAMESPACE/$GU_BUILD
  mkdir -p "$target_tmp_parent"
  GU_TARGET_TMP_DIR=$(mktemp -d "$target_tmp_parent/run.XXXXXX")
  GU_TARGET_VCF_DIR="$GU_TARGET_TMP_DIR/vcf"
  GU_METHOD_LOG_DIR=$GU_ANALYSIS_ROOT/$METHOD/$GU_TARGET_NAMESPACE/log
  mkdir -p "$GU_METHOD_LOG_DIR" "$GU_TARGET_VCF_DIR"
  GU_CHECK_LOG="$GU_METHOD_LOG_DIR/${METHOD}.${GU_SCOPE_ID}.log"
  : > "$GU_CHECK_LOG"
  gu_check_log "method=$METHOD action=${ACTION:-run}"
  gu_check_log "GRCH=$GU_BUILD"
  [[ -z $GU_LOCI_FILE ]] || gu_check_log "loci core=$GU_LOCI_CORE_FILE analysis=$GU_LOCI_FILE flank_bp=$GU_LOCI_FLANK_BP"
  gu_check_log "target=$GU_TARGET_NAMESPACE"
  gu_check_log "target genotype prefix=$GU_TARGET_GEN_PREFIX"
  [[ -z ${GU_TARGET_NATIVE_VCF_PREFIX:-} ]] || gu_check_log "native target VCF prefix=$GU_TARGET_NATIVE_VCF_PREFIX (autosomes chrN; male X chrX.male; run temp uses symlinks)"
  (( GU_CHRX_MALE_ONLY == 0 )) || gu_check_log "chrX default=male-only haploid non-PAR; PFILE preferred=$(gu_target_genotype_base "$GU_TARGET_GEN_PREFIX" X)"
  (( GU_CHRX_PAR_DIPLOID == 0 )) || gu_check_log "chrX locus contract=PAR diploid; source=$(gu_target_genotype_base "$GU_TARGET_GEN_PREFIX" X); PLINK PAR1/PAR2 export=internal chrX.vcf.gz"
  [[ $METHOD == trace || $METHOD == arg || -d $GU_ARCHAIC_ROOT ]] || gu_check_fail "archaic reference root is missing: $GU_ARCHAIC_ROOT"
  gu_check_log "archaic reference root=$GU_ARCHAIC_ROOT$([[ $METHOD == trace || $METHOD == arg ]] && printf ' (not used by %s)' "${METHOD^^}")"
  if [[ $METHOD == as3 ]]; then
    gu_check_log "AS3 GRCh38 reference panel=$AS3_REFERENCE_PANEL_DIR"
    gu_check_log "AS3 reference map=$AS3_REFERENCE_MAP"
    gu_check_log "AS3 3N1D mask release=$AS3_MASK_DIR (inventory-validated; not reapplied to preprocessed Ref1028)"
    gu_check_log "AS3 bundled code=$AS3_RUNTIME models=$AS3_MODEL_DIR"
    if [[ -d $AS3_PUBLISHED_CALLS_DIR ]] && find "$AS3_PUBLISHED_CALLS_DIR" -type f \( -name '*.bed' -o -name '*.bed.gz' \) -size +0c -print -quit | grep -q .; then
      gu_check_log "AS3 external callsets=$AS3_PUBLISHED_CALLS_DIR role=external_reference"
    else
      gu_check_log "AS3 external callsets=missing path=$AS3_PUBLISHED_CALLS_DIR"
    fi
    as3_preflight_chrs=${GU_CHRS:-"1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22"}
    [[ -s $AS3_REFERENCE_MAP ]] || gu_check_fail "AS3 Ref1028 map missing or empty: $AS3_REFERENCE_MAP"
    [[ -s $AS3_RUNTIME/ArchaicSeeker3.1-mamba && -d $AS3_RUNTIME/src ]] || gu_check_fail "bundled AS3 source is incomplete: $AS3_RUNTIME"
    [[ -s $AS3_RUNTIME/MODEL_SHA256SUMS ]] || gu_check_fail "bundled AS3 model checksum manifest is missing: $AS3_RUNTIME/MODEL_SHA256SUMS"
    (cd "$AS3_MODEL_DIR" && sha256sum -c "$AS3_RUNTIME/MODEL_SHA256SUMS" >/dev/null) || gu_check_fail "AS3 model files are missing or fail SHA-256 validation under $AS3_MODEL_DIR"
    read -r as3_afr as3_den as3_nean as3_other < <(awk 'BEGIN{FS="[[:space:]]+"}
      NF{total++; if(NF==2)n[$2]++}
      END{other=total-(n["AFR"]+n["DEN"]+n["NEAN"]);print n["AFR"]+0,n["DEN"]+0,n["NEAN"]+0,other+0}' "$AS3_REFERENCE_MAP")
    (( as3_afr == 146 && as3_den == 1 && as3_nean == 3 && as3_other == 0 )) || \
      gu_check_fail "AS3 Ref1028 map must contain exactly AFR=146, DEN=1, NEAN=3 in two columns: $AS3_REFERENCE_MAP"
    gu_check_log "AS3 reference map classes AFR=$as3_afr DEN=$as3_den NEAN=$as3_nean"
    for c in $as3_preflight_chrs; do
      c=${c#chr}; c=${c#CHR}
      [[ $c =~ ^([1-9]|1[0-9]|2[0-2])$ ]] || gu_check_fail "AS3 Ref1028 supports GRCh38 chr1-22 only; got chr$c"
      as3_vcf=$AS3_REFERENCE_PANEL_DIR/Ref_Panel.chr${c}.vcf.gz
      [[ -s $as3_vcf ]] || gu_check_fail "AS3 Ref1028 VCF missing or empty: $as3_vcf"
      [[ -s $as3_vcf.tbi || -s $as3_vcf.csi ]] || gu_check_fail "AS3 Ref1028 VCF index missing: $as3_vcf.tbi or $as3_vcf.csi"
      bcftools view -h "$as3_vcf" >/dev/null 2>&1 || gu_check_fail "AS3 Ref1028 VCF/header is unreadable: $as3_vcf"
      [[ $(bcftools index -n "$as3_vcf" 2>/dev/null || echo 0) -gt 0 ]] || gu_check_fail "AS3 Ref1028 VCF contains no indexed records: $as3_vcf"
      for as3_mask in "$AS3_MASK_DIR/chr${c}_mask.bed" "$AS3_MASK_DIR/chr${c}_mask.wchr.bed"; do
        [[ -s $as3_mask ]] || gu_check_fail "AS3 3N1D mask release is incomplete: $as3_mask"
      done
      awk -v chr="$c" 'NF{n++;if(NF<3||$1!~("^" chr "($|_)")||$2!~/^[0-9]+$/||$3!~/^[0-9]+$/||$2<0||$3<=$2)bad=1}END{exit !(n>0&&!bad)}' \
        "$AS3_MASK_DIR/chr${c}_mask.bed" || gu_check_fail "invalid numeric-contig 3N1D mask: $AS3_MASK_DIR/chr${c}_mask.bed"
      awk -v chr="$c" 'NF{n++;if(NF<3||$1!~("^chr" chr "($|_)")||$2!~/^[0-9]+$/||$3!~/^[0-9]+$/||$2<0||$3<=$2)bad=1}END{exit !(n>0&&!bad)}' \
        "$AS3_MASK_DIR/chr${c}_mask.wchr.bed" || gu_check_fail "invalid chr-prefixed 3N1D mask: $AS3_MASK_DIR/chr${c}_mask.wchr.bed"
      gu_check_log "AS3 Ref1028 chr$c=$as3_vcf masks=available"
    done
    gu_check_log "AS3 Ref1028/map/mask release preflight=PASS"
  fi

  if [[ -n $GU_CHRS ]]; then
    TARGET_CHECK_CHRS=$GU_CHRS
  elif [[ $METHOD == as3 ]]; then
    TARGET_CHECK_CHRS="1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22"
  else
    TARGET_CHECK_CHRS="1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X"
  fi
  mapfile -t target_check_array < <(tr ' ' '\n' <<< "$TARGET_CHECK_CHRS" | awk 'NF&&!seen[$1]++')
  (( ${#target_check_array[@]} )) || gu_check_fail "target chromosome list is empty"
  if [[ " $TARGET_CHECK_CHRS " == *" X "* && $GU_CHRX_TARGET_MODE == none ]]; then
    GU_CHRX_TARGET_MODE=male; GU_CHRX_MALE_ONLY=1
    gu_check_log "chrX default=male-only haploid non-PAR; PFILE preferred=$(gu_target_genotype_base "$GU_TARGET_GEN_PREFIX" X)"
  fi
  GU_BUILD_CHECK_CHR=${GU_BUILD_CHECK_CHR:-1}
  for c in "${target_check_array[@]}"; do
    if gu_detect_target_chr "$c"; then :; else
      rc=$?
      case "$rc" in
        3) gu_check_fail "VCF .tbi/.csi index missing for $(gu_target_genotype_base "$GU_TARGET_GEN_PREFIX" "$c").vcf.gz";;
        4) gu_check_fail "incomplete pfile trio for $(gu_target_genotype_base "$GU_TARGET_GEN_PREFIX" "$c")";;
        *) gu_check_fail "target chr$c not found; expected $(gu_target_genotype_base "$GU_TARGET_GEN_PREFIX" "$c").vcf.gz or $(gu_target_genotype_base "$GU_TARGET_GEN_PREFIX" "$c").{pgen,pvar,psam}";;
      esac
    fi
    GU_TARGET_SOURCE[$c]=$GU_DETECT_SOURCE; GU_TARGET_INDEX[$c]=${GU_DETECT_INDEX:-}; GU_TARGET_PVAR[$c]=$GU_DETECT_PVAR; GU_TARGET_PSAM[$c]=$GU_DETECT_PSAM; GU_TARGET_FORMAT[$c]=$GU_DETECT_FORMAT
    if [[ $c == X && $GU_DETECT_FORMAT == pfile ]]; then
      if (( GU_CHRX_PAR_DIPLOID )); then
        gu_validate_chrxy_pfile "$GU_DETECT_SOURCE" "$GU_BUILD" | tee -a "$GU_CHECK_LOG" || gu_check_fail "invalid chrXY target pfile: $GU_DETECT_SOURCE"
      else
        gu_validate_chrx_male_pfile "$GU_DETECT_SOURCE" "$GU_BUILD" | tee -a "$GU_CHECK_LOG" || gu_check_fail "invalid chrX.male target pfile: $GU_DETECT_SOURCE"
      fi
    fi
    if [[ $c == "$GU_BUILD_CHECK_CHR" ]]; then
      if [[ -n $GU_DETECT_PVAR ]]; then
        GU_TARGET_BUILD_SOURCE=$GU_DETECT_PVAR; GU_TARGET_BUILD_FORMAT=pfile
      else
        GU_TARGET_BUILD_SOURCE=$GU_DETECT_SOURCE; GU_TARGET_BUILD_FORMAT=$GU_DETECT_FORMAT
      fi
    fi
    gu_check_log "modern genome file chr$c=$GU_DETECT_SOURCE format=$GU_DETECT_FORMAT"
  done
  if [[ -z ${GU_TARGET_SOURCE[$GU_BUILD_CHECK_CHR]:-} ]]; then
    if gu_detect_target_chr "$GU_BUILD_CHECK_CHR"; then
      GU_TARGET_SOURCE[$GU_BUILD_CHECK_CHR]=$GU_DETECT_SOURCE
      GU_TARGET_INDEX[$GU_BUILD_CHECK_CHR]=${GU_DETECT_INDEX:-}
      GU_TARGET_PVAR[$GU_BUILD_CHECK_CHR]=$GU_DETECT_PVAR
      GU_TARGET_PSAM[$GU_BUILD_CHECK_CHR]=$GU_DETECT_PSAM
      GU_TARGET_FORMAT[$GU_BUILD_CHECK_CHR]=$GU_DETECT_FORMAT
      if [[ -n $GU_DETECT_PVAR ]]; then
        GU_TARGET_BUILD_SOURCE=$GU_DETECT_PVAR; GU_TARGET_BUILD_FORMAT=pfile
      else
        GU_TARGET_BUILD_SOURCE=$GU_DETECT_SOURCE; GU_TARGET_BUILD_FORMAT=$GU_DETECT_FORMAT
      fi
      gu_check_log "build-check genome file chr$GU_BUILD_CHECK_CHR=$GU_DETECT_SOURCE format=$GU_DETECT_FORMAT"
    else
      GU_BUILD_CHECK_CHR=${target_check_array[0]}
      if [[ -n ${GU_TARGET_PVAR[$GU_BUILD_CHECK_CHR]:-} ]]; then
        GU_TARGET_BUILD_SOURCE=${GU_TARGET_PVAR[$GU_BUILD_CHECK_CHR]}; GU_TARGET_BUILD_FORMAT=pfile
      else
        GU_TARGET_BUILD_SOURCE=${GU_TARGET_SOURCE[$GU_BUILD_CHECK_CHR]}; GU_TARGET_BUILD_FORMAT=${GU_TARGET_FORMAT[$GU_BUILD_CHECK_CHR]}
      fi
      gu_check_log "WARNING: chr1 is unavailable; build check will use requested chr$GU_BUILD_CHECK_CHR"
    fi
  fi
  target_formats=$(for c in "${target_check_array[@]}"; do printf '%s\n' "${GU_TARGET_FORMAT[$c]}"; done | sort -u | paste -sd+ -)
  GU_TARGET_GEN_FORMAT=$target_formats
  gu_check_log "target genotype format=${GU_TARGET_GEN_FORMAT^^}; build-check format=${GU_TARGET_BUILD_FORMAT^^}"

  if [[ -z $GU_SAMPLE_PANEL ]]; then
    target_parent=$(dirname -- "$GU_TARGET_GEN_PREFIX")
    for candidate in "$target_parent/samples.txt" "$(dirname -- "$target_parent")/samples.txt"; do
      if [[ -s $candidate ]]; then GU_SAMPLE_PANEL=$candidate; break; fi
    done
  fi
  GU_SAMPLE_PANEL_SOURCE_CHR=$GU_BUILD_CHECK_CHR
  if [[ ${GU_TARGET_FORMAT[$GU_SAMPLE_PANEL_SOURCE_CHR]:-} != pfile ]]; then
    for c in "${target_check_array[@]}"; do
      if [[ ${GU_TARGET_FORMAT[$c]} == pfile ]]; then GU_SAMPLE_PANEL_SOURCE_CHR=$c; break; fi
    done
  fi
  if [[ -z $GU_SAMPLE_PANEL && ${GU_TARGET_FORMAT[$GU_SAMPLE_PANEL_SOURCE_CHR]:-} == pfile ]]; then
    GU_SAMPLE_PANEL="$GU_TARGET_TMP_DIR/samples.txt"
    awk 'BEGIN{FS="[[:space:]]+";OFS="\t"}
      NR==1{for(i=1;i<=NF;i++){x=toupper($i);sub(/^#/,"",x);if(x=="IID"||x=="SAMPLE")id=i;if(x=="SEX")sex=i}if(!id||!sex){print "ERROR: pfile PSAM needs IID and SEX columns" > "/dev/stderr";exit 2};print "sample","sex";next}
      NF{print $id,$sex}' "${GU_TARGET_PSAM[$GU_SAMPLE_PANEL_SOURCE_CHR]}" > "$GU_SAMPLE_PANEL" || gu_check_fail "could not derive sample/sex metadata from ${GU_TARGET_PSAM[$GU_SAMPLE_PANEL_SOURCE_CHR]}"
  fi
  if [[ -n $GU_SAMPLE_PANEL ]]; then gu_check_log "sample metadata=$GU_SAMPLE_PANEL"
  elif [[ $METHOD == phyml ]]; then gu_check_log "WARNING: sample metadata was not found; phyml will infer chrX haploidy from GT width when possible"
  elif [[ $METHOD == trace || $METHOD == arg ]]; then gu_check_log "sample metadata was not found; ${METHOD^^} will use VCF/tree sample identities"
  elif [[ $METHOD == as3 ]]; then gu_check_log "WARNING: sample metadata was not found; AS3 will use target VCF sample IDs, but population-stratified normalization will be unavailable"
  else gu_check_fail "$METHOD requires sample metadata; use --sample-panel FILE or place samples.txt beside the target files"
  fi
  if (( GU_CHRX_MALE_ONLY )) && [[ ${GU_TARGET_FORMAT[X]:-} == vcf ]]; then
    [[ -n $GU_SAMPLE_PANEL && -s $GU_SAMPLE_PANEL ]] || gu_check_fail "native chrX.male.vcf.gz input requires --sample-panel with sample and sex columns"
    gu_validate_chrx_male_vcf "${GU_TARGET_SOURCE[X]}" "$GU_BUILD" "$GU_SAMPLE_PANEL" | tee -a "$GU_CHECK_LOG" || gu_check_fail "invalid chrX.male target VCF: ${GU_TARGET_SOURCE[X]}"
    gu_check_log "native chrX.male VCF contract=PASS: ${GU_TARGET_SOURCE[X]}"
  fi
  if (( GU_CHRX_PAR_DIPLOID )) && [[ ${GU_TARGET_FORMAT[X]:-} != pfile ]]; then
    gu_check_fail "PAR-scoped chrX loci require the PLINK2 chrXY pfile so PAR1/PAR2 and sex-aware ploidy are preserved"
  fi
  gu_check_log "preflight file check=PASS; log=$GU_CHECK_LOG"
else
  GU_TARGET_VCF_DIR=${TARGET_VCF_DIR_INPUT:-${GU_TARGET_VCF_DIR:-$GU_TARGET_ROOT/vcf}}
  if [[ $METHOD != ukb && -z $GU_SAMPLE_PANEL ]]; then
    GU_SAMPLE_PANEL=$GU_TARGET_ROOT/samples.txt
  fi
  if [[ $METHOD == arg ]]; then
    arg_check_chrs=${GU_CHRS:-"1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X"}
    if [[ " $arg_check_chrs " == *" X "* ]]; then
      GU_CHRX_MALE_ONLY=1
      arg_x_vcf=$GU_TARGET_VCF_DIR/chrX.vcf.gz
      [[ -s $arg_x_vcf && ( -s $arg_x_vcf.tbi || -s $arg_x_vcf.csi ) ]] || {
        echo "ERROR: ARG chrX requires an indexed pure-male haploid VCF at $arg_x_vcf" >&2
        exit 2
      }
      [[ -s $GU_SAMPLE_PANEL ]] || {
        echo "ERROR: ARG chrX validation requires --sample-panel with sample and sex columns" >&2
        exit 2
      }
      gu_validate_chrx_male_vcf "$arg_x_vcf" "$GU_BUILD" "$GU_SAMPLE_PANEL" || exit 2
    fi
  fi
fi

[[ -z $UKB_HAP_ROOT_INPUT ]] || export UKB_HAP_ROOT=$UKB_HAP_ROOT_INPUT
[[ -z $UKB_TYPED_ROOT_INPUT ]] || export UKB_TYPED_ROOT=$UKB_TYPED_ROOT_INPUT
[[ -z $UKB_WORK_INPUT ]] || export UKB_WORK=$UKB_WORK_INPUT
[[ -z $UKB_THREADS_INPUT ]] || export UKB_THREADS=$UKB_THREADS_INPUT
[[ -z $UKB_REF_FASTA_INPUT ]] || export UKB_REF_FASTA=$UKB_REF_FASTA_INPUT
[[ -z $UKB_KEEP_INPUT ]] || export UKB_KEEP=$UKB_KEEP_INPUT
[[ -z $UKB_1KG_VCF_DIR_INPUT ]] || export UKB_1KG_VCF_DIR=$UKB_1KG_VCF_DIR_INPUT
[[ -z $UKB_VCF_OUT_INPUT ]] || export UKB_VCF_OUT=$UKB_VCF_OUT_INPUT
[[ -z $UKB_ARG_VCF_OUT_INPUT ]] || export UKB_ARG_VCF_OUT=$UKB_ARG_VCF_OUT_INPUT
[[ -z $UKB_BATCH_SIZE_INPUT ]] || export GU_BATCH_SIZE=$UKB_BATCH_SIZE_INPUT
[[ -z $UKB_ANCHORS_PER_GROUP_INPUT ]] || export GU_ANCHORS_PER_GROUP=$UKB_ANCHORS_PER_GROUP_INPUT
[[ -z $UKB_ANCHOR_LIST_INPUT ]] || export GU_ANCHOR_LIST=$UKB_ANCHOR_LIST_INPUT
[[ -z $UKB_PANEL_CHR_INPUT ]] || export UKB_PANEL_CHR=$UKB_PANEL_CHR_INPUT

export GU_DATA_ROOT GU_REF_ROOT GU_ANALYSIS_ROOT GU_SOFT GU_BUILD GU_TARGET GU_TARGET_NAMESPACE GU_TARGET_DIR GU_TARGET_ROOT GU_TARGET_VCF_DIR GU_TARGET_GEN_PREFIX GU_TARGET_NATIVE_VCF_PREFIX GU_TARGET_GEN_FORMAT GU_CHRX_MALE_ONLY GU_CHRX_PAR_DIPLOID GU_CHRX_TARGET_MODE GU_CHRXY_EXTRACT_BED
export GU_TARGET_TMP_DIR
export GU_ARCHAIC_ROOT GU_SAMPLE_PANEL GU_ARG_DIR GU_LOCI_FILE GU_LOCI_CORE_FILE GU_LOCI_MAP_FILE GU_LOCI_FLANK_BP GU_CHRS GU_SCOPE_ID GU_SCOPE_KEY GU_SCOPE_LABEL GU_CHECK_LOG GU_BUILD_CHECK_CHR
export AS3_RUNTIME AS3_REFERENCE_PANEL_DIR AS3_REFERENCE_MAP AS3_MASK_DIR AS3_MODEL_DIR AS3_PUBLISHED_CALLS_DIR AS3_REFERENCE_CALLSET_CACHE
export dir0="$GU_DATA_ROOT" dir_ref="$GU_REF_ROOT" dir_archaic="$GU_ARCHAIC_ROOT" dir_1kg="$GU_REF_ROOT/1kg/$GU_BUILD" GRCH="$GU_BUILD"

unset R_LIBS R_LIBS_USER BCFTOOLS_PLUGINS
export R_ENVIRON_USER=/dev/null R_PROFILE_USER=/dev/null

PHE_F=${PHE_F:-$ROOT/../0f/0phe.f.sh}
[[ -s $PHE_F ]] || { echo "ERROR: missing shared helper: $PHE_F" >&2; exit 2; }
# shellcheck source=/mnt/d/scripts/0f/0phe.f.sh
source "$PHE_F"


# 🚩 Genome-build validation
# Use coordinate sentinels when pvars do not contain rsIDs.
gu_check_grch_pvar_positions(){
  local src=$1 expected=$2 gold=${CHECK_GRCH_SNP_LIST:-/mnt/d/data/ukb/phe/common/snp.lst}
  local n37 n38 inferred expected_n
  [[ -s $src && -s $gold ]] || { echo "ERROR: positional GRCh check input is missing: pvar=$src gold=$gold" >&2; return 2; }
  read -r n37 n38 < <(awk 'BEGIN{FS="[[:space:]]+"}
    NR==FNR{if(FNR>1){c=$3;sub(/^chr/,"",c);b37[c SUBSEP $4]=1;b38[c SUBSEP $5]=1}next}
    /^#/||NF<2{next}
    {c=$1;sub(/^chr/,"",c);k=c SUBSEP $2
     if((k in b37)&&!s37[k]++)n37++
     if((k in b38)&&!s38[k]++)n38++}
    END{print n37+0,n38+0}' "$gold" <(phe_zcat "$src"))
  if (( n37 > n38 )); then inferred=37; else inferred=38; fi
  if (( (n37 < 2 && n38 < 2) || n37 == n38 || (n37 > n38 ? n37-n38 : n38-n37) < 2 )); then
    echo "ERROR: positional GRCh check is inconclusive for $src (GRCh37 sentinels=$n37, GRCh38 sentinels=$n38)" >&2
    return 3
  fi
  expected_n=$(printf '%s' "$expected" | tr '[:upper:]' '[:lower:]')
  case "$expected_n" in 37|b37|grch37|hg19) expected_n=37;; 38|b38|grch38|hg38) expected_n=38;; *) echo "ERROR: invalid expected GRCh: $expected" >&2; return 2;; esac
  [[ $inferred == "$expected_n" ]] || { echo "ERROR: configured GRCH=$expected_n but positional check inferred GRCH=$inferred from $src" >&2; return 4; }
  echo "POSITION CHECK: GRCh37 sentinels=$n37, GRCh38 sentinels=$n38; inferred GRCh$inferred from coordinate-ID pvar"
}

gu_target_vcf_conversion_record(){
  local c=$1 pbase=$2 pvar=$3 psam=$4
  local x_mode=diploid
  if [[ $c == X ]]; then
    if (( GU_CHRX_PAR_DIPLOID )); then x_mode=par_diploid_chrxy
    else x_mode=male_haploid_nonpar
    fi
  fi
  printf 'schema\t4\nformat\tpfile\nchromosome\t%s\nx_mode\t%s\n' "$c" "$x_mode"
  stat -c 'pgen\t%n:%s:%Y' "${pbase}.pgen"
  stat -c 'pvar\t%n:%s:%Y' "$pvar"
  stat -c 'psam\t%n:%s:%Y' "$psam"
  if [[ -n $GU_LOCI_FILE && $METHOD != as3 ]]; then
    printf 'loci_content\t%s\n' "$(cksum < "$GU_LOCI_FILE")"
  else
    printf 'loci\twhole_chromosome\n'
  fi
}

gu_target_native_vcf_record(){
  local c=$1 source=$2 index=$3
  printf 'schema\t4\nformat\tvcf\nchromosome\t%s\n' "$c"
  stat -Lc 'vcf\t%n:%s:%Y' "$source"
  stat -Lc 'index\t%n:%s:%Y' "$index"
}

gu_prepare_target_vcfs(){
  local c chr_tmp source prepared prepared_index dest src_index dest_index part_prefix info pbase pvar psam info_prefix meta meta_part need_pfile=0 extract_bed
  command -v bcftools >/dev/null 2>&1 || gu_check_fail "bcftools is required"
  for c in "${target_check_array[@]}"; do [[ ${GU_TARGET_FORMAT[$c]} != pfile ]] || need_pfile=1; done
  if (( need_pfile )); then
    command -v plink2 >/dev/null 2>&1 || gu_check_fail "plink2 is required for pfile input"
    command -v tabix >/dev/null 2>&1 || gu_check_fail "tabix is required for pfile input"
  fi
  for c in "${target_check_array[@]}"; do
    chr_tmp=$GU_TARGET_TMP_DIR/chr$c
    mkdir -p "$chr_tmp"
    prepared=$chr_tmp/chr$c.vcf.gz
    dest=$GU_TARGET_VCF_DIR/chr$c.vcf.gz
    if [[ ${GU_TARGET_FORMAT[$c]} == vcf ]]; then
      source=${GU_TARGET_SOURCE[$c]}; src_index=${GU_TARGET_INDEX[$c]}
      [[ -s $src_index ]] || gu_check_fail "VCF index disappeared: $src_index"
      if [[ $src_index == *.tbi ]]; then prepared_index=$prepared.tbi; dest_index=$dest.tbi
      else prepared_index=$prepared.csi; dest_index=$dest.csi
      fi
      gu_check_log "link native target VCF chr$c=$source -> $prepared (no data copy)"
      gu_link_if_needed "$source" "$prepared" || gu_check_fail "could not link target VCF $source -> $prepared"
      if [[ $prepared_index == "$prepared.tbi" ]]; then rm -f "$prepared.csi" "$dest.csi"; else rm -f "$prepared.tbi" "$dest.tbi"; fi
      gu_link_if_needed "$src_index" "$prepared_index" || gu_check_fail "could not link target VCF index $src_index -> $prepared_index"
      gu_link_if_needed "$prepared" "$dest" || gu_check_fail "could not link prepared target VCF $prepared -> $dest"
      if [[ $dest_index == "$dest.tbi" ]]; then rm -f "$dest.csi"; else rm -f "$dest.tbi"; fi
      gu_link_if_needed "$prepared_index" "$dest_index" || gu_check_fail "could not link prepared target VCF index $prepared_index -> $dest_index"
      meta=$chr_tmp/source.tsv
      meta_part=$chr_tmp/.source.$$
      gu_target_native_vcf_record "$c" "$source" "$src_index" > "$meta_part"
      mv -f "$meta_part" "$meta"
    else
      pbase=${GU_TARGET_SOURCE[$c]}; pvar=${GU_TARGET_PVAR[$c]}; psam=${GU_TARGET_PSAM[$c]}
      local -a pfile_args=(--pfile "$pbase")
      [[ $pvar != *.zst ]] || pfile_args+=(vzs)
      if [[ $c == X && $GU_CHRX_MALE_ONLY == 1 ]]; then
        gu_check_log "chrX PFILE is already haploid; each male contributes one non-PAR X haplotype"
      elif [[ $METHOD == phyml || $METHOD == as3 || $METHOD == trace || $METHOD == arg ]]; then
        info_prefix=$chr_tmp/.pgen-info
        info=$(plink2 "${pfile_args[@]}" --pgen-info --out "$info_prefix" 2>&1) || { printf '%s\n' "$info" | tee -a "$GU_CHECK_LOG"; gu_check_fail "plink2 --pgen-info failed for $pbase"; }
        rm -f "$info_prefix.log"
        printf '%s\n' "$info" >> "$GU_CHECK_LOG"
        if printf '%s\n' "$info" | grep -Eqi 'no (hardcalls are explicitly phased|phased hardcalls)|phased hardcalls.*(no|absent)'; then
          gu_check_fail "$METHOD requires phased haplotypes, but $pbase has no phased hardcalls"
        fi
        printf '%s\n' "$info" | grep -Eqi 'explicitly phased hardcalls present|phased hardcalls.*present|phased hardcalls.*yes' || gu_check_fail "could not verify phased hardcalls in $pbase with plink2 --pgen-info"
      fi
      meta=$chr_tmp/source.tsv
      meta_part=$chr_tmp/.source.$$
      gu_target_vcf_conversion_record "$c" "$pbase" "$pvar" "$psam" > "$meta_part"
      part_prefix=$chr_tmp/.chr${c}.part.$$
      local -a extract_args=()
      if [[ -n $GU_LOCI_FILE && $METHOD != as3 ]]; then
        extract_bed=$GU_LOCI_FILE
        [[ $c != X || $GU_CHRX_PAR_DIPLOID != 1 ]] || extract_bed=$GU_CHRXY_EXTRACT_BED
        extract_args=(--extract bed0 "$extract_bed")
        gu_check_log "convert loci from pfile chr$c=$pbase -> $prepared (single-run temporary VCF)"
      else
        gu_check_log "convert whole chromosome pfile chr$c=$pbase -> $prepared (single-run temporary VCF)"
      fi
      plink2 "${pfile_args[@]}" "${extract_args[@]}" --export vcf bgz id-paste=iid --out "$part_prefix" >> "$GU_CHECK_LOG" 2>&1 || gu_check_fail "pfile-to-VCF conversion failed for $pbase"
      tabix -f -p vcf "$part_prefix.vcf.gz" >> "$GU_CHECK_LOG" 2>&1 || gu_check_fail "VCF indexing failed for $part_prefix.vcf.gz"
      mv -f "$part_prefix.vcf.gz" "$prepared"; mv -f "$part_prefix.vcf.gz.tbi" "$prepared.tbi"
      gu_link_if_needed "$prepared" "$dest" || gu_check_fail "could not link prepared target VCF $prepared -> $dest"
      gu_link_if_needed "$prepared.tbi" "$dest.tbi" || gu_check_fail "could not link prepared target VCF index $prepared.tbi -> $dest.tbi"
      mv -f "$meta_part" "$meta"
      rm -f "$part_prefix.log"
    fi
    bcftools view -h "$prepared" >/dev/null 2>&1 || gu_check_fail "unreadable prepared target VCF: $prepared"
    if [[ $c == X && $GU_CHRX_MALE_ONLY == 1 && ${GU_TARGET_FORMAT[$c]} == pfile ]]; then
      gu_validate_chrx_male_vcf "$prepared" "$GU_BUILD" "${GU_TARGET_PSAM[$c]}" | tee -a "$GU_CHECK_LOG" || \
        gu_check_fail "PLINK2 export did not preserve pure male haploid X: ${GU_TARGET_SOURCE[$c]}"
      gu_check_log "prepared chrX VCF contract=PASS: $prepared"
    fi
    gu_check_log "prepared modern genome chr$c=$prepared samples=$(bcftools query -l "$prepared" | awk 'NF{n++}END{print n+0}')"
  done
}

if (( USES_TARGET )); then
  gu_prepare_target_vcfs
  GU_BUILD_CHECK_INPUT=${GU_BUILD_CHECK_INPUT:-$GU_TARGET_BUILD_SOURCE}
  gu_check_log "verify modern genome build from $GU_BUILD_CHECK_INPUT"
  build_check_ok=0
  build_check_detail=$(mktemp "${TMPDIR:-/tmp}/gu-build-check.XXXXXX")
  if check_GRCH "$GU_BUILD_CHECK_INPUT" "$GU_BUILD" > "$build_check_detail" 2>&1; then
    cat "$build_check_detail" | tee -a "$GU_CHECK_LOG"
    build_check_ok=1
  elif [[ $GU_TARGET_BUILD_FORMAT == pfile ]] && gu_check_grch_pvar_positions "$GU_BUILD_CHECK_INPUT" "$GU_BUILD" >> "$build_check_detail" 2>&1; then
    gu_check_log "rsID build sentinels are absent; using strict coordinate-sentinel fallback"
    tail -n 1 "$build_check_detail" | tee -a "$GU_CHECK_LOG"
    build_check_ok=1
  else
    cat "$build_check_detail" | tee -a "$GU_CHECK_LOG" >&2
  fi
  rm -f "$build_check_detail"
  if (( ! build_check_ok )); then
    gu_check_fail "modern target build check failed; see $GU_CHECK_LOG"
  fi
  gu_check_log "GRCH check=PASS; configured=$GU_BUILD; source=${GU_TARGET_SOURCE[$GU_BUILD_CHECK_CHR]}"
  gu_check_log "GRCH $GU_BUILD; modern genome prefix=$GU_TARGET_GEN_PREFIX (${GU_TARGET_GEN_FORMAT^^}); prepared=$GU_TARGET_VCF_DIR/chr*.vcf.gz"
fi

if [[ -z ${IBDMIX_RUNTIME:-} ]]; then IBDMIX_RUNTIME=$GU_SOFT/ibdmix; fi
if [[ ! -d ${IBDMIX_RUNTIME:-} ]]; then IBDMIX_RUNTIME=$GU_SOFT/IBDmix; fi
GU_REQUEST_CHRS=${GU_CHRS:-}
if [[ -z $GU_REQUEST_CHRS ]]; then
  case "$METHOD" in
    ibdmix|trace) GU_REQUEST_CHRS="1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X" ;;
    as3) GU_REQUEST_CHRS="1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22" ;;
  esac
fi
GU_REQUEST_SCOPE_KEY=$GU_SCOPE_KEY
GU_REQUEST_SCOPE_LABEL=$scope_label
GU_REQUEST_LOCI_FILE=$GU_LOCI_FILE
GU_REQUEST_LOCI_CORE_FILE=$GU_LOCI_CORE_FILE
GU_REQUEST_LOCI_MAP_FILE=$GU_LOCI_MAP_FILE
if [[ -n $GU_REQUEST_LOCI_MAP_FILE ]]; then
  # Parentheses avoid awk parsing `print NR > ...` as output redirection.
  GU_REQUEST_UNIT_COUNT=$(awk 'END{print (NR > 0 ? NR - 1 : 0)}' "$GU_REQUEST_LOCI_MAP_FILE")
else
  GU_REQUEST_UNIT_COUNT=$(awk '{print NF}' <<< "$GU_REQUEST_CHRS")
fi
PHYML_OUT_OVERRIDE=${PHYML_OUT:-}
IBDMIX_OUT_OVERRIDE=${IBDMIX_OUT:-}
TRACE_OUT_OVERRIDE=${TRACE_OUT:-}
AS3_OUT_OVERRIDE=${AS3_OUT:-}

ensure_reference_callset_cache(){
  local -a cache_args=(--source "$AS3_PUBLISHED_CALLS_DIR" --cache "$AS3_REFERENCE_CALLSET_CACHE" --dataset-id "${AS3_REFERENCE_DATASET_ID:-AS3_1KG}")
  if python3 "$F/reference_callset_cache.py" check "${cache_args[@]}" --fast >/dev/null 2>&1; then
    echo "[GU NORMALIZE] reference callset cache=ready: $AS3_REFERENCE_CALLSET_CACHE"
    return 0
  fi
  echo "[GU NORMALIZE] reference callset cache missing or invalid; preparing automatically"
  python3 "$F/reference_callset_cache.py" prepare "${cache_args[@]}"
}

refresh_normalized_outputs(){
  local normalize_dir=${GU_NORMALIZE_DIR:-${GU_SHINY_DATA_DIR:-${GU_RSHINY_DIR:-$GU_ANALYSIS_ROOT/normalize}}}
  mkdir -p "$normalize_dir"
  local -a norm_args=(--analysis-root "$GU_ANALYSIS_ROOT" --output-dir "$normalize_dir" --database "${GU_SQLITE:-$GU_ANALYSIS_ROOT/gu.sqlite}" --build "GRCh$GU_BUILD" --reciprocal-overlap "${GU_CATALOG_RECIP_OVERLAP:-0.5}")
  if [[ -d ${AS3_PUBLISHED_CALLS_DIR:-} ]]; then
    ensure_reference_callset_cache
    norm_args+=(--reference-cache "$AS3_REFERENCE_CALLSET_CACHE")
  fi
  local pop_panel=${GU_NORMALIZE_SAMPLE_PANEL:-$GU_REF_ROOT/1kg/38/samples.txt}
  [[ ! -s $pop_panel ]] || norm_args+=(--sample-panel "$pop_panel")
  python3 "$F/normalize_results.py" "${norm_args[@]}"
}

maybe_refresh_normalized_outputs(){
  case "${GU_AUTO_NORMALIZE:-0}" in
    1) refresh_normalized_outputs ;;
    0) echo "[GU NORMALIZE] deferred; run ./gu.sh normalize after all chromosome/batch jobs (use --auto-normalize TRUE to rebuild after each run)" ;;
    *) echo "ERROR: internal auto-normalize state must be 0 or 1" >&2; return 2 ;;
  esac
}

stage_loci_metadata(){
  local out=$1 tmp
  [[ -n $GU_LOCI_FILE ]] || return 0
  mkdir -p "$out"
  tmp=$out/.request.loci.analysis.$$.bed; cp "$GU_LOCI_FILE" "$tmp"; gu_replace_if_changed "$tmp" "$out/request.loci.analysis.bed"
  tmp=$out/.request.loci.core.$$.bed; cp "$GU_LOCI_CORE_FILE" "$tmp"; gu_replace_if_changed "$tmp" "$out/request.loci.core.bed"
  tmp=$out/.request.loci.map.$$.tsv; cp "$GU_LOCI_MAP_FILE" "$tmp"; gu_replace_if_changed "$tmp" "$out/request.loci.map.tsv"
}

gu_replace_if_changed(){
  local tmp=$1 dest=$2
  if [[ ! -s $dest ]] || ! cmp -s "$tmp" "$dest"; then mv -f "$tmp" "$dest"; else rm -f "$tmp"; fi
}

use_staged_loci_metadata(){
  local out=$1
  [[ -n $GU_LOCI_FILE ]] || return 0
  GU_LOCI_FILE=$out/request.loci.analysis.bed
  GU_LOCI_CORE_FILE=$out/request.loci.core.bed
  GU_LOCI_MAP_FILE=$out/request.loci.map.tsv
  export GU_LOCI_FILE GU_LOCI_CORE_FILE GU_LOCI_MAP_FILE
}

gu_activate_chr(){
  local chr=$1
  GU_CHRS=$chr
  GU_SCOPE_KEY="${GU_REQUEST_SCOPE_KEY}.chr${chr}"
  GU_SCOPE_ID="chr${chr}"
  GU_UNIT_LABEL="chr${chr}"
  GU_OUTPUT_LAYOUT=per_chromosome
  GU_LOCUS_ID=""
  GU_LOCI_FILE=""; GU_LOCI_CORE_FILE=""; GU_LOCI_MAP_FILE=""
  export GU_CHRS GU_SCOPE_KEY GU_SCOPE_ID GU_UNIT_LABEL GU_OUTPUT_LAYOUT GU_LOCUS_ID GU_LOCI_FILE GU_LOCI_CORE_FILE GU_LOCI_MAP_FILE
}

gu_activate_locus(){
  local chr=$1 core_start=$2 core_end=$3 analysis_start=$4 analysis_end=$5 locus=$6 flank=$7 unit_label=$8
  local unit_key base tmp
  unit_key=$(printf '%s\n' "$chr|$core_start|$core_end|$locus" | cksum | awk '{print $1}')
  [[ -n ${GU_REGION_TMP_ROOT:-} ]] || { echo "ERROR: GU region temporary workspace is unavailable" >&2; return 2; }
  base=$GU_REGION_TMP_ROOT/${GU_REQUEST_SCOPE_KEY}.unit-${unit_key}
  mkdir -p "$(dirname -- "$base")"
  tmp=$base.analysis.tmp.$$
  printf '%s\t%s\t%s\t%s\n' "$chr" "$analysis_start" "$analysis_end" "$locus" > "$tmp"
  gu_replace_if_changed "$tmp" "$base.bed"
  tmp=$base.core.tmp.$$
  printf '%s\t%s\t%s\t%s\n' "$chr" "$core_start" "$core_end" "$locus" > "$tmp"
  gu_replace_if_changed "$tmp" "$base.core.bed"
  tmp=$base.map.tmp.$$
  printf 'chr\tcore_start\tcore_end\tanalysis_start\tanalysis_end\tlocus_id\tflank_bp\n%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
    "$chr" "$core_start" "$core_end" "$analysis_start" "$analysis_end" "$locus" "$flank" > "$tmp"
  gu_replace_if_changed "$tmp" "$base.map.tsv"
  GU_CHRS=$chr
  GU_SCOPE_KEY="${GU_REQUEST_SCOPE_KEY}.unit-${unit_key}"
  GU_SCOPE_ID=$unit_label
  GU_UNIT_LABEL=$unit_label
  GU_OUTPUT_LAYOUT=per_locus
  GU_LOCUS_ID=$locus
  GU_LOCI_FILE=$base.bed
  GU_LOCI_CORE_FILE=$base.core.bed
  GU_LOCI_MAP_FILE=$base.map.tsv
  export GU_CHRS GU_SCOPE_KEY GU_SCOPE_ID GU_UNIT_LABEL GU_OUTPUT_LAYOUT GU_LOCUS_ID GU_LOCI_FILE GU_LOCI_CORE_FILE GU_LOCI_MAP_FILE
}


# 🚩 Analysis units
# Group AS3 loci into one task per chromosome.
gu_activate_loci_chr(){
  local chr=$1 base tmp label
  [[ -n ${GU_REGION_TMP_ROOT:-} ]] || { echo "ERROR: GU region temporary workspace is unavailable" >&2; return 2; }
  base=$GU_REGION_TMP_ROOT/${GU_REQUEST_SCOPE_KEY}.group-chr${chr}
  mkdir -p "$(dirname -- "$base")"
  tmp=$base.analysis.tmp.$$
  awk -F'\t' -v c="$chr" 'BEGIN{OFS="\t"}$1==c{print}' "$GU_REQUEST_LOCI_FILE" > "$tmp"
  [[ -s $tmp ]] || { rm -f "$tmp"; echo "ERROR: no loci for chr$chr" >&2; return 2; }
  gu_replace_if_changed "$tmp" "$base.bed"
  tmp=$base.core.tmp.$$
  awk -F'\t' -v c="$chr" 'BEGIN{OFS="\t"}$1==c{print}' "$GU_REQUEST_LOCI_CORE_FILE" > "$tmp"
  gu_replace_if_changed "$tmp" "$base.core.bed"
  tmp=$base.map.tmp.$$
  awk -F'\t' -v c="$chr" 'NR==1||$1==c' "$GU_REQUEST_LOCI_MAP_FILE" > "$tmp"
  gu_replace_if_changed "$tmp" "$base.map.tsv"
  label="chr${chr}_${GU_REQUEST_SCOPE_LABEL}"
  GU_CHRS=$chr
  GU_SCOPE_KEY="${GU_REQUEST_SCOPE_KEY}.group-chr${chr}"
  GU_SCOPE_ID=$label
  GU_UNIT_LABEL=$label
  GU_OUTPUT_LAYOUT=per_chromosome_loci
  GU_LOCUS_ID=""
  GU_LOCI_FILE=$base.bed
  GU_LOCI_CORE_FILE=$base.core.bed
  GU_LOCI_MAP_FILE=$base.map.tsv
  export GU_CHRS GU_SCOPE_KEY GU_SCOPE_ID GU_UNIT_LABEL GU_OUTPUT_LAYOUT GU_LOCUS_ID GU_LOCI_FILE GU_LOCI_CORE_FILE GU_LOCI_MAP_FILE
}

gu_output_for_unit(){
  local method=$1 override=${2:-}
  gu_chr_result_dir "$GU_ANALYSIS_ROOT" "$method" "$GU_UNIT_LABEL" "$GU_TARGET_NAMESPACE" "$override" "$GU_REQUEST_UNIT_COUNT"
}

run_per_analysis_unit(){
  local runner=$1 chr core_start core_end analysis_start analysis_end locus flank unit_label status=0
  local -a failed=()
  (( GU_REQUEST_UNIT_COUNT > 0 )) || { echo "ERROR: no analysis units selected for $METHOD" >&2; return 2; }
  if [[ -n $GU_REQUEST_LOCI_MAP_FILE ]]; then
    while IFS=$'\t' read -r chr core_start core_end analysis_start analysis_end locus flank; do
      unit_label=$(gu_locus_unit_label "$GU_REQUEST_LOCI_MAP_FILE" "$chr" "$core_start" "$core_end" "$locus")
      if (
        gu_activate_locus "$chr" "$core_start" "$core_end" "$analysis_start" "$analysis_end" "$locus" "$flank" "$unit_label"
        echo "[GU RUN] method=$METHOD analysis_unit=$GU_UNIT_LABEL chromosome=chr$chr locus=$locus layout=per-locus"
        "$runner"
      ); then
        echo "[GU RUN] method=$METHOD analysis_unit=$unit_label status=complete"
      else
        status=1; failed+=("$unit_label")
        echo "[GU RUN] method=$METHOD analysis_unit=$unit_label status=failed" >&2
      fi
    done < <(awk 'NR>1' "$GU_REQUEST_LOCI_MAP_FILE")
  else
    [[ -n $GU_REQUEST_CHRS ]] || { echo "ERROR: no chromosomes selected for $METHOD" >&2; return 2; }
    for chr in $GU_REQUEST_CHRS; do
      if (
        gu_activate_chr "$chr"
        echo "[GU RUN] method=$METHOD analysis_unit=$GU_UNIT_LABEL chromosome=chr$chr layout=per-chromosome"
        "$runner"
      ); then
        echo "[GU RUN] method=$METHOD analysis_unit=chr$chr status=complete"
      else
        status=1; failed+=("chr$chr")
        echo "[GU RUN] method=$METHOD analysis_unit=chr$chr status=failed" >&2
      fi
    done
  fi
  if (( status )); then
    echo "ERROR: $METHOD failed analysis units: ${failed[*]}; completed unit directories were retained" >&2
    return 1
  fi
}

run_loci_grouped_by_chr(){
  local runner=$1 chr status=0
  local -a failed=()
  for chr in $GU_REQUEST_CHRS; do
    if (
      gu_activate_loci_chr "$chr"
      echo "[GU RUN] method=$METHOD analysis_unit=$GU_UNIT_LABEL chromosome=chr$chr layout=per-chromosome-loci"
      "$runner"
    ); then
      echo "[GU RUN] method=$METHOD analysis_unit=chr$chr status=complete"
    else
      status=1; failed+=("chr$chr")
    fi
  done
  if (( status )); then
    echo "ERROR: $METHOD failed chromosome groups: ${failed[*]}" >&2
    return 1
  fi
}

run_phyml_one(){
  export PHYML_OUT
  PHYML_OUT=$(gu_output_for_unit phyml "$PHYML_OUT_OVERRIDE")
  echo "[GU RUN] output=$PHYML_OUT"
  stage_loci_metadata "$PHYML_OUT"
  use_staged_loci_metadata "$PHYML_OUT"
  bash "$F/phyml.sh" "${ACTION:-run}"
}

run_ibdmix_one(){
  local a=${ACTION:-run} ga
  case "$a" in run)ga=ibdmix_run;;check)ga=ibdmix_check;;esac
  local IBDMIX_OUT
  IBDMIX_OUT=$(gu_output_for_unit ibdmix "$IBDMIX_OUT_OVERRIDE")
  echo "[GU RUN] output=$IBDMIX_OUT"
  stage_loci_metadata "$IBDMIX_OUT"
  use_staged_loci_metadata "$IBDMIX_OUT"
  env GU_ACTION="$ga" IBDMIX_LOCUS_FLANK_BP=0 dir0="$GU_DATA_ROOT" dir_ref="$GU_REF_ROOT" dir_archaic="$GU_ARCHAIC_ROOT" dirarch="$GU_ARCHAIC_ROOT" dirmod="$GU_TARGET_ROOT" sample_file="$GU_SAMPLE_PANEL" dirscript="$F" dirsoft="$IBDMIX_RUNTIME" GRCH="$GU_BUILD" genome_build="b$GU_BUILD" dirout="$IBDMIX_OUT" selected_loci="$GU_LOCI_FILE" chrs="$GU_CHRS" refs="${IBDMIX_REFS:-Altai Chagyr Vindija Denisova}" lod_cut="${IBDMIX_LOD:-4}" len_cut="${IBDMIX_MIN_BP:-50000}" job_of_chr="${IBDMIX_JOB_OF_CHR:-1}" job_in_chr="${IBDMIX_JOB_IN_CHR:-8}" bash "$F/ibdmix.sh"
}

run_trace_request(){
  GU_CHRS=$GU_REQUEST_CHRS
  GU_SCOPE_KEY=$GU_REQUEST_SCOPE_KEY
  GU_SCOPE_ID=$GU_REQUEST_SCOPE_LABEL
  GU_UNIT_LABEL=$GU_REQUEST_SCOPE_LABEL
  GU_OUTPUT_LAYOUT=request
  GU_LOCUS_ID=""
  GU_LOCI_FILE=$GU_REQUEST_LOCI_FILE
  GU_LOCI_CORE_FILE=$GU_REQUEST_LOCI_CORE_FILE
  GU_LOCI_MAP_FILE=$GU_REQUEST_LOCI_MAP_FILE
  export GU_CHRS GU_SCOPE_KEY GU_SCOPE_ID GU_UNIT_LABEL GU_OUTPUT_LAYOUT GU_LOCUS_ID
  export GU_LOCI_FILE GU_LOCI_CORE_FILE GU_LOCI_MAP_FILE
  export TRACE_CHRS=$GU_CHRS TRACE_OUT
  TRACE_OUT=$(gu_chr_result_dir "$GU_ANALYSIS_ROOT" trace "$GU_UNIT_LABEL" "$GU_TARGET_NAMESPACE" "$TRACE_OUT_OVERRIDE" 1)
  echo "[GU RUN] method=trace analysis_request=$GU_UNIT_LABEL chromosomes=$GU_CHRS layout=request output=$TRACE_OUT"
  stage_loci_metadata "$TRACE_OUT"
  use_staged_loci_metadata "$TRACE_OUT"
  export TRACE_LOCI_FILE=$GU_LOCI_FILE
  bash "$F/trace.sh" "${ACTION:-run}"
}

run_trace_one(){
  export TRACE_CHRS=$GU_CHRS TRACE_OUT
  TRACE_OUT=$(gu_output_for_unit trace "$TRACE_OUT_OVERRIDE")
  echo "[GU RUN] method=trace analysis_unit=$GU_UNIT_LABEL chromosomes=$GU_CHRS layout=$GU_OUTPUT_LAYOUT output=$TRACE_OUT"
  stage_loci_metadata "$TRACE_OUT"
  use_staged_loci_metadata "$TRACE_OUT"
  export TRACE_LOCI_FILE=$GU_LOCI_FILE
  bash "$F/trace.sh" "${ACTION:-run}"
}

run_as3_one(){
  local a=${ACTION:-run} experimental_x=0 as3_loci_tmp as3_tmp_label
  export AS3_CHRS=$GU_CHRS AS3_ALLOW_CHRX=0
  [[ " $AS3_CHRS " != *" X "* ]] || { echo "ERROR: AS3 chrX is unavailable because the official GRCh38 Ref1028 panel contains chr1-22 only" >&2; return 2; }
  # Preserve chromosome-scale model context. A loci scope is applied when calls
  # are normalized, rather than by concatenating small disjoint VCF slices.
  export AS3_LOCI_FILE=""
  export AS3_TARGET_VCF_DIR=${AS3_TARGET_VCF_DIR:-$GU_TARGET_VCF_DIR}
  export AS3_DATA_OUT AS3_OUT AS3_RUNTIME_WORK_DIR
  # AS3 loci use whole-chromosome model context. The validated manifest and
  # upstream working directory exist only for this top-level GU invocation.
  if [[ -n $GU_REQUEST_LOCI_FILE ]]; then
    as3_tmp_label=chr$GU_CHRS
  else
    as3_tmp_label=$GU_UNIT_LABEL
  fi
  AS3_DATA_OUT=$GU_TARGET_TMP_DIR/as3-input/$as3_tmp_label
  AS3_RUNTIME_WORK_DIR=$GU_TARGET_TMP_DIR/as3-runtime/$as3_tmp_label
  AS3_OUT=$(gu_output_for_unit as3 "$AS3_OUT_OVERRIDE")
  echo "[GU RUN] output=$AS3_OUT temporary_input=$AS3_DATA_OUT"
  mkdir -p "$AS3_OUT"
  printf 'key\tvalue\ntarget\t%s\ngenome_build\tGRCh%s\nscope_id\t%s\nchromosome\t%s\nlocus_id\t%s\noutput_layout\t%s\nexperimental_chrX\t%s\nexperimental_build\t%s\nmodel_context\twhole_chromosome\ntarget_chunk_size\t%s\ncgroup_memory_high\t%s\ncgroup_memory_max\t%s\ncgroup_swap_max\t%s\n' \
    "$GU_TARGET_NAMESPACE" "$GU_BUILD" "$GU_SCOPE_ID" "$GU_CHRS" "$GU_LOCUS_ID" "$GU_OUTPUT_LAYOUT" "$experimental_x" "$AS3_EXPERIMENTAL_BUILD" "$AS3_TARGET_CHUNK_SIZE" "${AS3_MEMORY_HIGH:-20G}" "${AS3_MEMORY_MAX:-24G}" "${AS3_MEMORY_SWAP_MAX:-8G}" > "$AS3_OUT/run.meta.tsv"
  if [[ -n $GU_LOCI_FILE ]]; then
    as3_loci_tmp=$AS3_OUT/.input.loci.$$.bed; cp "$GU_LOCI_FILE" "$as3_loci_tmp"; gu_replace_if_changed "$as3_loci_tmp" "$AS3_OUT/input.loci.bed"
    as3_loci_tmp=$AS3_OUT/.input.loci.core.$$.bed; cp "$GU_LOCI_CORE_FILE" "$as3_loci_tmp"; gu_replace_if_changed "$as3_loci_tmp" "$AS3_OUT/input.loci.core.bed"
    as3_loci_tmp=$AS3_OUT/.input.loci.map.$$.tsv; cp "$GU_LOCI_MAP_FILE" "$as3_loci_tmp"; gu_replace_if_changed "$as3_loci_tmp" "$AS3_OUT/input.loci.map.tsv"
    AS3_CLIP_LOCI_FILE=$AS3_OUT/input.loci.bed
  else
    rm -f "$AS3_OUT/input.loci.bed" "$AS3_OUT/input.loci.core.bed" "$AS3_OUT/input.loci.map.tsv"
    AS3_CLIP_LOCI_FILE=""
  fi
  export AS3_CLIP_LOCI_FILE
  local -a prep_args=(--reference-panel-dir "$AS3_REFERENCE_PANEL_DIR" --reference-map "$AS3_REFERENCE_MAP" --out "$AS3_DATA_OUT" --chr "$AS3_CHRS" --genome-build "b$GU_BUILD")
  if [[ $a == check ]]; then
    bash "$F/as3_prep.sh" "${prep_args[@]}" --check-only || return
    (
      export GU_ACTION=as3_env_check AS3_DATA_IN="$AS3_DATA_OUT" AS3_RUNTIME
      bash "$F/as3.sh"
    ) || return
    return
  fi
  # Always re-enter preparation: it cheaply revalidates direct target/Ref1028
  # inputs and atomically refreshes the manifest when provenance changes.
  bash "$F/as3_prep.sh" "${prep_args[@]}" || return
  (
    export GU_ACTION="as3_$a" AS3_DATA_IN="$AS3_DATA_OUT" AS3_RUNTIME
    bash "$F/as3.sh"
  ) || return
}


# 🚩 Method dispatch
case "$METHOD" in
  phyml)
    run_per_analysis_unit run_phyml_one
    [[ ${ACTION:-run} == check ]] || maybe_refresh_normalized_outputs
    ;;
  ibdmix)
    run_per_analysis_unit run_ibdmix_one
    [[ ${ACTION:-run} == check ]] || maybe_refresh_normalized_outputs
    ;;
  trace)
    if [[ " $GU_REQUEST_CHRS " == *" X "* ]] && (( $(awk '{print NF}' <<< "$GU_REQUEST_CHRS") > 1 )); then
      echo "[GU RUN] TRACE mixed autosome/chrX request: running chromosome-scoped outputs because ARG sample-node ploidy differs"
      if [[ -n $GU_REQUEST_LOCI_FILE ]]; then
        run_loci_grouped_by_chr run_trace_one
      else
        run_per_analysis_unit run_trace_one
      fi
    else
      run_trace_request
    fi
    if [[ ${ACTION:-run} == run || ${ACTION:-run} == infer || ${ACTION:-run} == segments ]]; then maybe_refresh_normalized_outputs; fi
    ;;
  as3)
    if [[ -n $GU_REQUEST_LOCI_FILE ]]; then
      run_loci_grouped_by_chr run_as3_one
    else
      run_per_analysis_unit run_as3_one
    fi
    if [[ ${ACTION:-run} == run ]]; then maybe_refresh_normalized_outputs; fi
    ;;
  arg) bash "$F/arg.sh" "${ACTION:-build}" ;;
  normalize) refresh_normalized_outputs ;;
  shiny)
    export GU_SQLITE=${GU_SQLITE:-$GU_ANALYSIS_ROOT/gu.sqlite}
    exec Rscript -e "shiny::runApp('$ROOT/shiny', host=Sys.getenv('GU_SHINY_HOST','127.0.0.1'), port=as.integer(Sys.getenv('GU_SHINY_PORT','3838')), launch.browser=interactive())"
    ;;
  ukb) exec bash "$F/ukb.sh" "$ACTION" ;;
esac
}

{ gu_main "$@"; exit $?; }
