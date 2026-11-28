#!/usr/bin/env bash
set -euo pipefail

usage() {
    cat <<'HELP'
LE8: molecular association, prediction and causal follow-up.

Usage:
  ./le8.sh run|MODULE[,MODULE...] [options]

Modules:
  run                             C1 timing/twins and Any1 neighbor benchmarks; C2 cause -> C3 coloc -> C4 connect -> C5 cellulation.
  c1_correlate                    Associations, first-year ICD profiles and individual predictivity.
  c1_reliability                  Frozen common-selector reliability across models (opt-in).
  c2_cause                        MR and genetic decomposition (CPU).
  c3_coloc                        Colocalization (CPU).
  c4_connect                      Connections, panel validation and explanations (CPU).
  c5_cellulation                  Cell annotation (CPU).
  pgs_focus                       Matched-PGS analysis (opt-in).
  install                         Install/check runtime; run separately from analysis.
  evidence | viewer               Standalone evidence / viewer utilities; use MODULE --help.
  Publication                     Use ./le8.paper.sh for paper, share and shiny.

Results contract: only this code release is accepted (RELEASE.json).
Start a new release in an empty analysis directory; no migration, repair or old-result fallback.

Options:
  --Y CSV                         Outcomes [Y or cvd_cad,ra].
  --dx-date CSV                   Primary ICD with Any1 neighbor sensitivity for all traits. Optional sensitivity CSV.
  --biom CSV                      Molecular layers: prot,met [BIOM or prot,met].
  --analysis-root DIR             Results root [LE8_ANALYSIS_ROOT or ANALYSIS_ROOT/le8].
  --replace TRUE|FALSE            Recompute / reuse complete results [TRUE].
  --config FILE                   C1 method configuration [config/c1_methods.json].
  --python FILE                   Scientific Python [shared LE8 runtime].
  --r-bin FILE                    Rscript; HPC workers use shared 0f/0os.r.sh.
  --plan-dir DIR                  Plan root [analysis-root/jobs].
  --seed N                        Scientific random seed [2026].
  --submit auto|local|hpc         Select execution environment [auto].
  --queue NAME                    CPU queue.
  --cores N                       CPU budget from shared 0f/0os.sh.
  --jobs-N N                      Shared workflow concurrency setting.
  --wait TRUE|FALSE               Wait for submitted workflow [FALSE].
  --plan-only                     Write commands/plan without executing/submitting.
  --methods CSV                   C1 methods; requires c1_correlate. Config defaults:
                                  clinical,glmnet,lightgbm,abm,abm_clinical,abm_deep,abm_prevalent.
  --abm-deep-run cpu|gpu|auto     Device for selected abm_deep [cpu]; gpu requests CUDA.
  --gpu-queue NAME                GPU queue [LE8_GPU_QUEUE / HPC_GPU_QUEUE].
  --gpu-cores N                   CPU slots per deep GPU job [LE8_GPU_CORES or 8].
  --association run|skip          C1 association stage [run].
  --pgs TRUE|FALSE                Molecular PGS associations [FALSE].
  --ygb-phe DIR                   YGB phenotype root [shared YGB_PHE].
  --pwas-workers N                Forward PWAS worker budget to the scientific engine.
  --pgs-workers N                 Forward PGS worker budget to the scientific engine.
  --mrlink2-workers N             Forward MRlink2 worker budget.
  --mrlink2-inner-threads N       Native threads per MRlink2 worker.
  -h, --help                      Show this help.

Examples:
  ./le8.sh run --Y cvd_cad,cvd_stroke_i,cvd_pad,cvd_afib,cvd_hfail,mi --biom prot,met
  ./le8.sh run --Y cvd_cad,cvd_hfail --biom prot,met
  nohup ./le8.sh c1_correlate --Y cvd_cad,cvd_stroke_i --dx-date icd,ref,icd1,ref1,any1 --biom prot,met &
  ./le8.sh c1_correlate --Y cvd_cad,cvd_stroke_i --dx-date any1 --biom prot,met
  ./le8.sh c1_correlate --Y cvd_cad --dx-date any1 --biom prot --replace FALSE
  # Preview the full analysis before submission.
  ./le8.sh run --Y cvd_cad,cvd_htn,cvd_pad,stroke,copd --biom prot,met --queue 33 --cores 40 --plan-only
  ./le8.sh run --Y cvd_cad,cvd_htn,cvd_pad,stroke,copd --biom prot,met --submit hpc --queue 33 --cores 40
  # C1 with CPU methods only.
  ./le8.sh c1_correlate --Y cvd_cad --biom prot --methods clinical,glmnet,lightgbm,abm,abm_clinical --submit hpc
  ./le8.sh c1_correlate --Y cvd_cad --biom prot --methods abm_deep --abm-deep-run gpu --plan-only
  ./le8.sh c2_cause,c3_coloc --Y cvd_cad --biom prot --plan-only
  # After analysis jobs finish successfully:
  ./le8.paper.sh run --Y cvd_cad --biom prot --plan-only
HELP
}
case "${1:-}" in
    ""|-h|--help|help) usage; exit 0 ;;
esac
if [[ $1 == --plan-only ]]; then set -- run "$@"; fi
if [[ $1 == run ]]; then shift; set -- c1_correlate,c2_cause,c3_coloc,c4_connect,c5_cellulation "$@"; fi
if [[ $1 == *, ]]; then usage >&2; exit 2; fi
IFS=, read -ra selected_modules <<< "$1"
for module in "${selected_modules[@]}"; do
    case "$module" in
        install|evidence|viewer|c1_correlate|c1_reliability|c2_cause|c3_coloc|c4_connect|c5_cellulation|pgs_focus) ;;
        *) printf 'Unknown module: %s\n' "$module" >&2; usage >&2; exit 2 ;;
    esac
done
LE8_ROOT=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)
source "$LE8_ROOT/f/entry.sh"
set -- "${WORKFLOW_ARGS[@]}"
source "$LE8_ROOT/../0f/console.sh"
{
    le8_main "$@"
    exit $?
}
