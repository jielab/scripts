#!/usr/bin/env bash
# Push one configured tree; private paths are excluded at every depth.
set -euo pipefail
source "$(dirname -- "${BASH_SOURCE[0]}")/0os.sh"
folder=scripts
apply=0
while (($#)); do
    case "$1" in
        --folder) folder=${2:?directory category required}; shift 2 ;;
        --apply) apply=1; shift ;;
        --dry-run) apply=0; shift ;;
        -h|--help)
            printf '%s\n' 'sync.sh [--folder scripts|phe|files|software|analysis|gen|gwas|metag|LD|annot] [--apply]' \
                'Default: preview scripts LOCAL -> HPC. --apply copies; never deletes destination files.' \
                'All paths containing 🈲, caches and tests are excluded. SSH uses your agent/key or prompts for a password.'
            exit 0 ;;
        *) echo "Unknown argument: $1" >&2; exit 2 ;;
    esac
done
[[ $submit == local ]] || { echo 'Run sync.sh on LOCAL.' >&2; exit 2; }
case "$folder" in
    scripts) source_dir=$SCRIPTS_ROOT; destination=$HPC_ROOT/scripts ;;
    phe|files|software|analysis) source_dir=$DIR_DAT/$folder; destination=$HPC_ROOT/$folder ;;
    gen|gwas|metag|LD|annot) source_dir=$DIR_BIG/$folder; destination=$HPC_ROOT/data/$folder ;;
    *) echo "Unknown directory category: $folder" >&2; exit 2 ;;
esac
[[ $source_dir != *"$LOCAL_ONLY_MARK"* && -d $source_dir ]] || { echo 'Missing or private source directory.' >&2; exit 2; }
options=(-rt --omit-dir-times --itemize-changes --protect-args
    --exclude="*$LOCAL_ONLY_MARK*"
    --exclude='__pycache__/' --exclude='*.pyc' --exclude='*.pyo'
    --exclude='.pytest_cache/' --exclude='.mypy_cache/' --exclude='.ruff_cache/'
    --exclude='tests/' --exclude='nohup.out')
# Compare small source files by content, preserving mtimes with -t.
# Large data trees retain the efficient size/time comparison.
[[ $folder != scripts ]] || options+=(--checksum)
((apply)) || options+=(--dry-run)
# No symlink traversal and no --delete: running jobs and host-specific results survive.
exec rsync "${options[@]}" -e "ssh -p $HPC_PORT" "$source_dir/" "$HPC_USER@$HPC_HOST:$destination/"
