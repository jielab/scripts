#!/usr/bin/env bash
# Single OS-specific R launcher. Environment changes stay in this subprocess.
set -e
os_r_dir=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)
source "$os_r_dir/0os.sh"
if [[ $submit == hpc ]]; then
    os_r_setup
    exec "$OS_RSCRIPT_HPC" "$@"
fi
exec Rscript "$@"
