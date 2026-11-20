#!/usr/bin/env bash
set -euo pipefail
ROOT=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)
for arg in "$@"; do
    case "$arg" in -h|--help)
        cat <<'HELP'
GRID: regional models, ancestry-aware prediction and evaluation.

Quick start (after CSx score/fit AND Disco have completed):
  ./grid.sh all --grch 38 --traits height,ldl,t2dm
  ./grid.sh all --grch 38 --traits height,ldl,t2dm --plan-only

Usage: ./grid.sh [module] [options]
Modules (default: all):
  all           regions -> architecture / abm / block_models -> twins / evaluate -> final
                Includes abm_deep on GPU by default; excludes robustness
  regions       Build genome-wide block score features from existing CSx weights
  block_models  Fit regional glmnet and LightGBM models (CPU)
  abm           Fit ancestry-aware neighbour model (CPU)
  twins         Export ABM neighbour diagnostics
  abm_deep      Fit block-token deep model (CPU or explicit GPU)
  evaluate      Evaluate completed predictions
  architecture  Summarize regional architecture
  robustness    Run geometry sensitivity analysis (opt-in)
  final         Aggregate completed reports
  check         Check required inputs; do not submit analysis jobs
Single modules require their upstream artifacts to exist; they do not schedule
prerequisite modules. Use all to schedule the internal dependency graph.

Main options:
  --grch 38               Target genome build [38 only]
  --traits CSV            height,ldl,t2dm [all three]
  --config FILE           YAML/JSON configuration [config/grid.yml]
  --analysis-root DIR     Output/input root [${ANALYSIS_ROOT}/grid/38]
  --genotype-dir DIR      Imputed chrN pfiles; alias: --dir-imp
  --consensus-dir DIR     Existing GRCh38 chr*.bed catalog; no LD rebuilding
  --geometry MODE         consensus [default], regional, global_prs, legacy_unit_block
  --bootstrap N           Evaluation bootstrap replicates [200]
  --with-deep             Include abm_deep [enabled in default config]
  --without-deep          Omit abm_deep from all
  --abm-deep-run MODE     cpu|gpu|auto; auto means GPU [config: cuda]
                         cpu keeps deep enabled but runs it on CPU
  --submit MODE           auto|hpc|local [auto]
  --queue NAME            CPU queue [auto]
  --cores N               CPU slots per job [shared 0f/0os.sh WORKFLOW_CORES]
  --gpu-queue NAME        Deep GPU queue [config: 2v100-32-e5]
  --gpu-cores N           CPU slots per deep GPU job [8]
  --gpu-memory-limit-gb N Deep host RAM cap [HPC_GPU_MEMORY_GIB or 160 GiB];
                         actual reservation is min(cap, configured job memory)
  --job-hours N           Maximum duration per job [120 h]
  --threads N             Scientific CPU budget [--cores]
  --jobs N                Concurrent local jobs [1]; not an HPC concurrency limit
  --replace TRUE|FALSE    Recompute stages / reuse matching completed outputs [FALSE]
  --wait TRUE|FALSE       Wait for submitted HPC workflow [FALSE]; local waits
  --plan-only             Write commands/plan without submitting; skips input validation
  -h, --help              Show help without loading Python or accessing inputs

Prerequisites and parallelism:
Run 2.csx.sh all to completion (chromosome inference -> merge -> score -> fit),
then 3.disco.sh all to completion, before submitting GRID for the same traits.
GRID consumes their existing cohort, weights, benchmark and Disco predictions;
it neither launches upstream workflows nor attaches dependencies to active CSx
or Disco jobs. Missing required artifacts stop submission. A plan-only preview
is allowed earlier and does not prove that upstream results are ready.
You may process one fully finished trait while other traits are still running.

GRID has NO --chr option. Chromosomes come from config/grid.yml [1-22]; use the
same chromosome set as upstream. Regional features and models combine these
chromosomes. LSF jobs are scheduled by trait and ready stage. Within each score/
regions job, independent chromosomes run in parallel and are merged in order.
One all command schedules this graph automatically; --jobs 22 is unnecessary.

CPU execution within a job:
  --cores 40 reserves and budgets 40 CPU cores; --threads can lower the
    scientific budget. Defaults come only from shared 0f/0os.sh, including
    environment overrides via WORKFLOW_CORES (currently HPC 40 / local 16).
  Score/regions: normally 10 chromosome processes x 4 PLINK threads with the
    default 8 GiB PLINK memory; concurrency decreases if the RAM budget requires.
  ABM consensus search: independent query partitions, up to 40 processes;
    each process uses one BLAS thread and read-only shared arrays.
  Baseline/Disco: independent models and family-held-out folds run in parallel.
  Evaluation: paired bootstrap replicates run in parallel with the same RNG draws.
  glmnet: independent alpha paths run in parallel (3 by default); each lambda
    path remains sequential. LightGBM already receives the full CPU budget.
  Native thread settings also reach libraries loaded later and subprocesses.
  Logs print effective CPU/process budgets. Parsing, output assembly, I/O and
    dependent optimizer iterations can still use fewer cores; 40 is a budget,
    not a promise of constant 100% utilization. Memory limits bound concurrency.
  --jobs controls concurrent local workflows, not these within-job processes.

Default GPU, CPU-only and custom GPU examples:
  ./grid.sh all --grch 38 --traits height,ldl,t2dm --submit hpc
  ./grid.sh all --grch 38 --traits height --without-deep
  ./grid.sh all --grch 38 --traits height --abm-deep-run cpu
  ./grid.sh all --grch 38 --traits height --gpu-queue 2v100-32-e5 --gpu-cores 8
  ./grid.sh abm_deep --grch 38 --traits height --abm-deep-run gpu
  ./grid.sh evaluate --grch 38 --traits height
  ./grid.sh final --grch 38 --traits height
The standalone deep example requires existing regions and upstream predictions;
wait for deep to finish before evaluate, then wait for evaluate before final.
Only abm_deep requests a GPU (one CUDA GPU per job). ABM, glmnet, LightGBM and
other GRID stages use CPU. GPU queue: --gpu-queue, runtime.gpu_queue in config,
HPC_GPU_QUEUE, then 4a100-40. One exclusive GPU per deep job, with separate CPU
slots and host RAM limits, as in LE8. Evaluation waits for deep and runs on CPU.
Host RAM is not
GPU VRAM. PRSformer is a separate optional GPU workflow via 4.prsformer.sh;
grid.sh all does not launch it. GPU execution requires CUDA-enabled PyTorch; unavailable CUDA fails without
silent CPU fallback. For a CPU-only machine, use --abm-deep-run cpu or --without-deep.
HELP
        exit 0
        ;;
    esac
done
exec bash "$ROOT/f/run_entry.sh" grid "$@"
