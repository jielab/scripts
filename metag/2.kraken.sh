#!/usr/bin/env bash
# BEGIN QIMING RESOURCE EXAMPLES (Admin, 2026-11-05)
# Queue 33: 40 cores / 192 GiB per node. Reserve 10/20/30/40 cores by
# rounding max(required_cores/40, required_memory_GiB/192) UP to 25% steps.
# Example for a 24-core, 160-GiB CPU worker:
#   bsub -q 33 -n 40 -W 24:00 -M 160GB -R "span[hosts=1] rusage[mem=160GB/host]" ...
# A small <=10-core, <=48-GiB job uses -q 33 -n 10; request its actual RAM.
# For exact CPU counts, use Admin-suggested v4-3t or v3-6t after queueinfo.
# The examples are allocation guidance, not measured peak RAM for every stage.
# Use this entry with submit=hpc and --config hpc.yml; it submits its own DAG.
# Do not wrap the planner in bsub. --submit FALSE writes commands for review.
# New .cmd files round CPU slots automatically. With the current hpc.yml:
# prepare: -q 38 -n 16 / 64 GiB; microbes: -q 38 -n 64 / 448 GiB;
# plants/invertebrates: -q 1t81c -n 120 / 896 GiB;
# viruses: -q 33 -n 10 / 48 GiB; vertebrates: -q v3-6t -n 64 / 1792 GiB.
# Details and queue limits: ../0f/HPC_SUBMISSION_POLICY.md
# END QIMING RESOURCE EXAMPLES

set -euo pipefail
exec bash "$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)/f/entry.sh" split "$@"
