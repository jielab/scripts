#!/usr/bin/env bash

# Display standalone help before loading runtimes or scheduling work.
if [[ ${BASH_SOURCE[0]} == "$0" && ( ${1:-} == -h || ${1:-} == --help ) ]]; then
    cat <<'PUBLIC_ENTRY_HELP'
source ./memory_cap.sh  # 内部函数库，由所属入口加载
# 无独立分析模块；请通过所属项目的主 .sh 选择模块
# 参数和环境变量沿用调用方；不要作为全流程提交

Modules:
    Library  Shared functions/configuration; source from the owning entry script

Main options:
    No standalone analysis options. Read the variables/functions below.
    -h, --help  Show this library help

Examples:
    source ./memory_cap.sh
PUBLIC_ENTRY_HELP
    exit 0
fi

# Direct CLI help must return before setup, downloads or analysis.
if [[ "${BASH_SOURCE[0]}" == "$0" ]]; then
	case "${1:-}" in
		-h | --help)
			bash "${BASH_SOURCE[0]}" --help
			exit 0
			;;
	esac
fi
# Whole-process cgroup limits shared by GU and GRID; sourcing only defines functions.
# A scope preserves the caller's environment, working directory and terminal.
# Compatibility arguments expose working-set budgets; no additional OS memory cap.


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
	export MEMORY_CAP_MODE=none
	return 0
}


# Compatibility wrapper: execute without an extra memory limit.
resource_cap_main() (
    set -euo pipefail
    local_mode=${1:?memory wrapper mode or size required}
    shift
    case "$local_mode" in
        --check) exit 0 ;;
        --inside) shift ;;
    esac
    exec "$@"
)

if [[ ${BASH_SOURCE[0]} == "$0" ]]; then
	resource_cap_main "$@"
fi
