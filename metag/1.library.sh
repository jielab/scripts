#!/usr/bin/env bash
set -euo pipefail
HERE="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
case "${1:-}" in
  -h|--help|help)
    cat <<'EOF'
Build permanent reference libraries (includes reference/taxid checks).

Usage:
  ./1.library.sh                         Build sourmash + MetaCache
  ./1.library.sh --workers 16            Set sourmash CPU worker processes
  ./1.library.sh --prepare-workers 2     Parallel MetaCache input preparation (1-4)
  ./1.library.sh --part-gib 3             Set MetaCache partition target
  ./1.library.sh --memory-limit-gb 16     Whole-build RAM cap (1-16 GiB; swap off)
  ./1.library.sh --check-only [--force]   Check references without building
  ./1.library.sh --sourmash-only          Build only sourmash
  ./1.library.sh --metacache-only         Build only MetaCache
  ./1.library.sh --kraken [--max-db-gb 72] Optional comparator

Default permanent paths (f/metag.conf; outside the analysis output directory):
  READ / retain FASTA: /mnt/f/annot/metag/raw
  READ source/taxonomy metadata: /mnt/f/annot/metag/0.metadata
  WRITE inventory/check records: /mnt/f/annot/metag/lib/reference
  WRITE sourmash library: /mnt/f/annot/metag/lib/sourmash
  WRITE MetaCache databases/taxonomy: /mnt/f/annot/metag/lib/metacache
  Library root (also optional comparators): /mnt/f/annot/metag/lib
  Configure DIR_REF / DIR_META / DIR_LIB in f/metag.conf to relocate these paths.
  Changing DIR_OUT alone does not relocate or erase these permanent libraries.

FASTA: /mnt/f/annot/metag/raw
Libraries: /mnt/f/annot/metag/lib
Settings: f/metag.conf
MetaCache GPU builds require user systemd/cgroup v2 isolation.
GPU build failures stop; there is no automatic CPU fallback.
Default RAM cap: 16 GiB shared by the entire workflow, with no swap.
Input preparers and the single GPU builder inherit that same cgroup.
Limit verification failure refuses to start; a killed build stops without retry.
Partition GiB groups compressed sources, not RAM. Derived inputs split long
records into ~1 MiB pieces and <=512 MiB bases per cached input shard.
Up to four unbuilt input shards form one <=2 GiB GPU database.
Two preparation processes overlap input generation with GPU construction.
Source FASTA is retained; completed matching databases and input shards are reused.
Rerun this command to resume. Completed, matching outputs are reused.
EOF
    exit 0;;
esac
command=build-db
case "${1:-}" in
  --check-only) command=check-ref; shift;;
  --sourmash-only) command=build-sourmash; shift;;
  --metacache-only) command=build-metacache; shift;;
  --kraken) command=build-kraken; shift;;
esac
export METAG_PROG="./1.library.sh"
exec bash "$HERE/f/run_workflow.sh" "$command" "$@"
