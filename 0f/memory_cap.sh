#!/usr/bin/env bash
# Whole-process cgroup limits shared by GU and GRID; sourcing only defines functions.
# A scope preserves the caller's environment, working directory and terminal.
# Limit the entire process tree; ulimit -v is unsuitable for CUDA address spaces.


# 🚩 memory_cap_enter
memory_cap_enter() {
	local label=$1 ram=$2 swap=$3 script=$4
	shift 4
	local original=("$@") help=false
	MEMORY_CAP_ARGS=()
	while (($#)); do
		case "$1" in
			--memory-cap-gb | --swap-cap-gb)
				if (($# < 2)); then
					echo "ERROR: missing value for $1" >&2
					exit 2
				fi
				if [[ $1 == --memory-cap-gb ]]; then ram=$2; else swap=$2; fi
				shift 2
				;;
			--memory-cap-gb=*)
				ram=${1#*=}
				shift
				;;
			--swap-cap-gb=*)
				swap=${1#*=}
				shift
				;;
			-h | --help | help)
				help=true
				MEMORY_CAP_ARGS+=("$1")
				shift
				;;
			*)
				MEMORY_CAP_ARGS+=("$1")
				shift
				;;
		esac
	done
	if [[ ! $ram =~ ^[1-9][0-9]{0,4}$ || ! $swap =~ ^(0|[1-9][0-9]{0,4})$ ]]; then
		echo 'ERROR: --memory-cap-gb must be a positive integer; --swap-cap-gb must be a nonnegative integer (GiB).' >&2
		exit 2
	fi
	# Public outputs consumed by the sourcing launchers.
	# shellcheck disable=SC2034
	MEMORY_CAP_GB=$ram
	# shellcheck disable=SC2034
	MEMORY_SWAP_CAP_GB=$swap
	[[ $help == false ]] || return 0
	local ram_bytes=$((ram * 1073741824)) swap_bytes=$((swap * 1073741824))
	local cg unit
	cg=$(awk -F: '$1 == "0" {print $3}' /proc/self/cgroup)
	if [[ -n ${_PIPELINE_MEMORY_SCOPE:-} && $cg == */"$_PIPELINE_MEMORY_SCOPE" ]]; then
		# Verify kernel limits after re-entry, never trust the environment marker alone.
		local base="/sys/fs/cgroup$cg"
		if [[ ! -r $base/memory.max || ! -r $base/memory.swap.max || ! -r $base/memory.oom.group ]] ||
			[[ $(cat "$base/memory.max") != "$ram_bytes" ||
			$(cat "$base/memory.swap.max") != "$swap_bytes" ||
			$(cat "$base/memory.oom.group") != 1 ]]; then
			echo "ERROR: $label memory limits were not applied; refusing to run uncapped." >&2
			exit 2
		fi
		echo "[$label] MEMORY CAP: RAM=${ram} GiB, swap=${swap} GiB, whole process tree; PID=$$; scope=$_PIPELINE_MEMORY_SCOPE" >&2
		return 0
	fi
	if ! command -v systemd-run >/dev/null || [[ ! -r /sys/fs/cgroup/cgroup.controllers ]]; then
		echo "ERROR: $label requires systemd with cgroup v2 for the memory cap; refusing to run uncapped." >&2
		exit 2
	fi
	export XDG_RUNTIME_DIR=${XDG_RUNTIME_DIR:-/run/user/$(id -u)}
	export DBUS_SESSION_BUS_ADDRESS=${DBUS_SESSION_BUS_ADDRESS:-unix:path=$XDG_RUNTIME_DIR/bus}
	if ! systemctl --user show-environment >/dev/null 2>&1; then
		echo "ERROR: $label cannot reach the systemd user manager; refusing to run uncapped." >&2
		exit 2
	fi
	unit="$label-memory-$$-$RANDOM.scope"
	# Keep this tiny supervisor outside the capped scope so a killed task can
	# report its limits and journal location instead of only "Terminated".
	local status
	if systemd-run --user --scope --quiet --collect --unit="$unit" \
		-p MemoryAccounting=yes \
		-p "MemoryMax=$ram_bytes" -p "MemorySwapMax=$swap_bytes" -p OOMPolicy=kill \
		-- env "_PIPELINE_MEMORY_SCOPE=$unit" bash "$script" "${original[@]}"; then
		exit 0
	else
		status=$?
		echo "[$label] Task exited with status $status (RAM cap=${ram} GiB, swap cap=${swap} GiB)." >&2
		echo "[$label] Diagnose this scope: journalctl --user -u $unit --no-pager" >&2
		exit "$status"
	fi
}


# 🚩 resource_cap_main
resource_cap_main() (
	# A scope accounts for the whole process tree, including detached workers.
	set -euo pipefail

	memory_bytes() {
		local value=${1^^} number
		[[ $value =~ ^([1-9][0-9]*)(G|GB|GIB)$ ]] || {
			echo 'ERROR: --memory-cap requires a positive whole number of GiB (for example 32G).' >&2
			return 2
		}
		number=${BASH_REMATCH[1]}
		[[ ${#number} -le 6 ]] || return 2
		printf '%s\n' "$((10#$number * 1073741824))"
	}

	inside_cap() {
		local bytes=$1 current expected=${GU_MEMORY_CGROUP:-}
		current=$(awk -F: '$1=="0" {print $3}' /proc/self/cgroup)
		[[ -n $expected && $expected != / && ($current == "$expected" || $current == "$expected/"*) ]] || return 1
		[[ $(cat "/sys/fs/cgroup$expected/memory.max") == "$bytes" &&
		$(cat "/sys/fs/cgroup$expected/memory.swap.max") == 0 &&
		$(cat "/sys/fs/cgroup$expected/memory.oom.group") == 1 ]]
	}

	mode=${1:-}
	shift
	if [[ $mode == --check ]]; then
		bytes=$(memory_bytes "${1:?}")
		inside_cap "$bytes"
		exit
	fi
	if [[ $mode == --inside ]]; then
		size=${1:?}
		shift
		bytes=$(memory_bytes "$size")
		GU_MEMORY_CGROUP=$(awk -F: '$1=="0" {print $3}' /proc/self/cgroup)
		[[ $GU_MEMORY_CGROUP == */gu-memory-*.scope ]] || {
			echo 'ERROR: memory cap scope was not created.' >&2
			exit 1
		}
		# systemd 255 does not expose memory.oom.group as a scope property.
		# Kill the entire run on OOM instead of leaving partial workers running.
		echo 1 >"/sys/fs/cgroup$GU_MEMORY_CGROUP/memory.oom.group"
		export GU_MEMORY_CGROUP
		inside_cap "$bytes" || {
			echo 'ERROR: memory cap verification failed.' >&2
			exit 1
		}
		echo "[GU MEMORY] cap=$size (GiB units) swap=0 scope=$GU_MEMORY_CGROUP; shared by all workers; exceeding cap terminates this run" >&2
		exec "$@"
	fi
	size=$mode
	bytes=$(memory_bytes "$size")
	if inside_cap "$bytes"; then exec "$@"; fi
	[[ -f /sys/fs/cgroup/cgroup.controllers ]] && command -v systemd-run >/dev/null &&
		systemctl --user show-environment >/dev/null 2>&1 || {
		echo 'ERROR: cannot enforce --memory-cap: cgroup v2 and a running systemd user manager are required. Analysis was not started.' >&2
		exit 1
	}
	exec systemd-run --user --scope --quiet --unit="gu-memory-$$-$RANDOM" \
		-p "MemoryMax=$bytes" -p MemorySwapMax=0 -p KillMode=control-group -p TimeoutStopSec=10s \
		bash "${BASH_SOURCE[0]}" --inside "$size" "$@"
)

if [[ ${BASH_SOURCE[0]} == "$0" ]]; then
	resource_cap_main "$@"
fi
