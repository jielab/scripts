#!/usr/bin/env bash
# Public entry-point resource guard. Kept separate from arg.sh helpers so
# resource-policy changes do not invalidate prepared-genotype fingerprints.
arg_run_guarded() {
  local action=$1; shift
  local ram=${ARG_MEMORY_LIMIT_GB:-32} swap=${ARG_MEMORY_SWAP_GB:-4}
  local value state status
  local -a args=()
  while (($#)); do
    case "$1" in
      --memory-limit-gb|--memory-swap-gb)
        [[ $# -ge 2 && $2 =~ ^[0-9]{1,6}$ ]] || {
          echo "ERROR: $1 requires an integer between 0 and 999999." >&2; return 2;
        }
        if [[ $1 == --memory-limit-gb ]]; then ram=$2; else swap=$2; fi
        shift 2;;
      --jobs|--arg-jobs)
        [[ $# -ge 2 && $2 =~ ^[0-9]{1,6}$ ]] || {
          echo "ERROR: $1 requires a positive integer." >&2; return 2;
        }
        value=$((10#$2))
        ((value>0)) || { echo "ERROR: $1 must be positive." >&2; return 2; }
        if [[ $action == build ]] && ((value>1)); then
          echo "ARG memory guard: reducing build --jobs $value to 1; full-chromosome inference can exceed 20 GiB per worker." >&2
          value=1
        fi
        args+=("$1" "$value"); shift 2;;
      *) args+=("$1"); shift;;
    esac
  done
  [[ $ram =~ ^[0-9]{1,6}$ && $swap =~ ^[0-9]{1,6}$ ]] || {
    echo "ERROR: ARG memory limits must be integers between 0 and 999999 GiB." >&2; return 2;
  }
  ram=$((10#$ram)); swap=$((10#$swap))
  ((ram>0)) || { echo "ERROR: ARG RAM cap must be positive; unguarded runs are disabled." >&2; return 2; }
  [[ -f /sys/fs/cgroup/cgroup.controllers ]] || {
    echo "ERROR: ARG requires cgroup v2 memory control; refusing an unguarded run." >&2; return 2;
  }
  state=$(systemctl --user is-system-running 2>/dev/null || true)
  if ! command -v systemd-run >/dev/null 2>&1 || [[ $state != running && $state != degraded ]]; then
    echo "ERROR: user systemd is unavailable; refusing an unguarded ARG run." >&2
    return 2
  fi
  # Set native-library threads before any Python, NumPy or external tool starts.
  export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
  [[ $action != check ]] || args+=(--action check)
  if [[ $action == build ]]; then
    echo "ARG memory guard: $ram GiB RAM + $swap GiB swap for the entire task tree; build concurrency is limited to 1." >&2
  else
    echo "ARG memory guard: $ram GiB RAM + $swap GiB swap for the entire task tree." >&2
  fi
  # Launch the internal driver directly: no environment flag can bypass the
  # guard, and the calling terminal stays outside this scope.
  if systemd-run --user --scope --quiet --collect \
      -p "MemoryMax=${ram}G" -p "MemorySwapMax=${swap}G" -p OOMPolicy=kill \
      bash "$ROOT/f/arg.sh" logged "$action" "${args[@]}"; then
    return 0
  else
    status=$?
    echo "ERROR: ARG exited with status $status (configured RAM cap: $ram GiB)." >&2
    if ((status==137)); then
      echo "SIGKILL: check kernel logs for a cgroup OOM, global OOM, or external kill; this exit code alone does not identify the cause." >&2
    fi
    echo "Inspect journalctl -k -b --no-pager; after a WSL restart use -b -1. Verified checkpoints remain reusable." >&2
    return "$status"
  fi
}
