#!/usr/bin/env bash
set -euo pipefail
usage() {
    cat <<'HELP'
HealthFlux v3 | UK Biobank independent reproduction

Usage: ./flux.sh run|MODULE [options]
  run         Fit the full model, run held-out evaluation, baselines, blackout,
              annual rollout and both architecture ablations (8 jobs on HPC).
  prepare | train | readout | evaluate | baseline | dynamics | blackout | rollout
  ablate      Run just the two ablations using an existing trained main model.
  demo        End-to-end synthetic smoke test, including both forward experiments.

  --dir-ygb DIR          Input root [0f/0os.sh: DIR_YGB]
  --dir-out DIR          Project output [ANALYSIS_ROOT/flux]; demo uses flux/demo.
  --config FILE          YAML overrides. Default global training batch 3072 on one GPU.
                         config/prepare-memory.yml: in-memory prepare on v3-6t (256 GB, 8 cores).
  --endpoint-file FILE   Prespecified one-column ICD-10 three-character codes and DEATH.
                         With no file, use the local training vocabulary.
  --device auto|gpu|cpu  Auto selects GPU on HPC and available hardware locally.
  --gpus N               GPUs for training [1]; accumulation adjusts to preserve global batch.
  --cores N              Per-job CPU budget [site default].
  --no-ablation          Skip the two additional fits (7 jobs on HPC).
                         By default no-ODE disables propagation between observations;
                         no-update disables recurrent state revision at observations.
  --resume               train only: resume this release's interrupted training.
  --replace              Refit individual stages.
  --plan-only            Write commands under jobs/ without submitting or running.
  main/results/          Model, risk arrays and endpoint tables.

Full YGB analysis (default)
  ./flux.sh run
  Start with a new output directory. On HPC, train releases readout, dynamics and ablate.
  After readout, evaluate-baseline (CPU), blackout and rollout can run independently.
After all analysis jobs have finished successfully, generate the reports separately
  ./flux.paper.sh run
Preview the analysis jobs
  ./flux.sh run --plan-only
Small synthetic example
  ./flux.sh demo --device cpu --cores 2
HELP
}
case "${1:-}" in ""|help|-h|--help) usage; exit 0;; esac
for arg in "$@"; do case "$arg" in -h|--help) usage; exit 0;; esac; done
case "$1" in report|paper) echo 'Run ./flux.paper.sh run to generate summaries and figures.' >&2; exit 2;; esac
FLUX_ROOT=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)
export FLUX_ROOT PYTHONDONTWRITEBYTECODE=1
source "$FLUX_ROOT/f/runtime.sh"
flux_runtime "$@"
exec "$FLUX_PYTHON" "$FLUX_ROOT/f/flux.py" "$@"
