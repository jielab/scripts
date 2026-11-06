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
HERE="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
if [[ ${1:-} == --software-only ]]; then
    shift
    source "$HERE/f/env.sh"
    workflow_args "$@"
    set -- "${WORKFLOW_ARGS[@]}"
    exec "$METAG_PYTHON" "$HERE/f/install_tools.py" "$@"
fi
case "${1:-}" in
    -h|--help|help)
        cat <<'HELP'
Usage: bash install.sh [setup] [options]
Install the pinned CPU tools through the local / LSF setup workflow.

  --software-only     Install/verify tools without reference or cohort data (put first)
  --config FILE       Runtime paths/resources (defaults: hpc.yml)
  --submit FALSE      Generate setup commands without running or submitting them
  --submit TRUE       Run/submit the generated setup commands
  submit=local|hpc     Environment variable selecting execution mode
  jobs_N=N            Environment variable controlling local concurrency

Examples:
  bash install.sh --config hpc.yml --submit FALSE
  submit=local bash install.sh --config hpc.yml --submit TRUE
  submit=hpc bash install.sh --config hpc.yml --submit TRUE

Configuration:
  software.yml  Source versions, archives and SHA256 checksums
  hpc.yml              Tool/output paths and job resources
  f/env.sh             Interpreter and site compiler environment
Requires Python >= 3.10 with PyYAML; gcc, g++, make, Perl and Bash on workers.
Bundled Kraken sources are compiled with project native extensions; optional
Kun-Peng uses its pinned release binary. Existing installations are verified.
Licenses: Kraken GPL-3.0 in its source archive; Kun-Peng in LICENSE.kunpeng.
README files are optional documentation and are never installation inputs.
HELP
        exit 0 ;;
esac
if [[ $# == 0 || ${1:0:1} == - ]]; then
    set -- setup "$@"
fi
exec bash "$HERE/f/entry.sh" "$@"
