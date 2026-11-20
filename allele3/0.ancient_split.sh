#!/usr/bin/env bash

# Public help returns before runtime setup, input access or job submission.
if [[ ${BASH_SOURCE[0]} == "$0" ]]; then
    case "${1:-}" in
        ""|-h|--help|help)
            cat <<'PUBLIC_ENTRY_HELP'
ALLELE3: split ancient pfiles and lift coordinates to GRCh38.

Usage:
  ./0.ancient_split.sh run [options]

Modules:
  run                              Split GRCh37 chromosomes and produce GRCh38 pfiles.

Options:
  --grch 38                        Target genome build.
  --input-prefix PREFIX            Source ancient pfile prefix.
  --out37 DIR --out38 DIR          Output directories.
  --runtime DIR                    Job and log directory.
  --scratch-root DIR               Temporary workspace [/tmp on the worker].
  --chain FILE --fasta FILE        LiftOver chain and reference FASTA.
  --chromosomes auto|CSV           Chromosomes to process [auto].
  --jobs N --threads N             Concurrent chromosomes / PLINK threads [4 / 4].
  --memory-mb N                    Memory per PLINK process [8192].
  --submit auto|local|hpc          Execution mode [auto].
  --queue NAME --cores N           Queue and CPU budget from shared 0f settings.
  --job-hours N                    Maximum job duration.
  --replace TRUE|FALSE             Recompute or reuse completed outputs [FALSE].
  --plan-only                      Write commands without running or submitting.
  -h, --help                       Show this help.

Examples:
  ./0.ancient_split.sh run --grch 38 --plan-only
  ./0.ancient_split.sh run --grch 38 --chromosomes 1,2
PUBLIC_ENTRY_HELP
            exit 0 ;;
        all) printf 'Use run for the complete workflow; see -h.\n' >&2; exit 2 ;;
    esac
fi
if [[ ${1:-} == run ]]; then shift; fi
export submit=hpc SUBMIT=hpc
source "${ALLELE_OS_CONFIG:-$ROOT/../0f/0os.sh}"
export PATH="${DIR_WORK}/.conda/envs/allele3-evolution/bin:${DIR_SOFTWARE}/bin:$PATH"
exec "${ALLELE_PYTHON:-${DIR_WORK}/.conda/envs/allele3-evolution/bin/python3}" "$ROOT/f/ancient_split.py" "$@"
