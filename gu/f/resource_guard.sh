#!/usr/bin/env bash
# Resource policy for a complete GU invocation. Source from gu.sh.
# LSF owns limits on compute nodes; systemd owns whole-process-tree limits locally.
# The CPU fallback is explicitly a per-process virtual-memory ceiling.
gu_resource_guard() {
	local script=$1 mode=${GU_RESOURCE_GUARD:-auto} size=${GU_MEMORY_CAP:-32G}
	shift
	local number bytes kb current effective
	[[ ${size^^} =~ ^([1-9][0-9]{0,5})(G|GB|GIB)$ ]] || {
		echo 'ERROR: --memory-cap requires a positive integer GiB, e.g. 32G' >&2; return 2;
	}
	number=${BASH_REMATCH[1]}; bytes=$((10#$number * 1073741824)); kb=$((bytes / 1024))
	case "$mode" in auto|systemd|scheduler|ulimit) ;; *)
		echo "ERROR: unknown resource guard: $mode (auto|systemd|scheduler|ulimit)" >&2; return 2 ;;
	esac
	if [[ $mode == auto ]]; then
		if [[ -n ${LSB_JOBID:-} ]]; then
			mode=scheduler
		elif command -v systemd-run >/dev/null && systemctl --user show-environment >/dev/null 2>&1; then
			mode=systemd
		else
			mode=ulimit
		fi
	fi
	if [[ ${GU_RUN_CMD:-TRUE} == FALSE ]]; then
		echo "[GU RESOURCE] command generation only; no inference resource scope entered" >&2
		return 0
	fi
	case "$mode" in
		scheduler)
			[[ -n ${LSB_JOBID:-} ]] || {
				echo 'ERROR: scheduler guard requires an actual LSF job (LSB_JOBID). Use auto locally.' >&2; return 2;
			}
			echo "[GU RESOURCE] LSF job=$LSB_JOBID; budget=$size; CPU/memory enforcement belongs to the submitted LSF reservation" >&2
			;;
		systemd)
			if ! bash "$ROOT/f/memory_cap.sh" --check "$size"; then
				exec bash "$ROOT/f/memory_cap.sh" "$size" bash "$script" "$@"
			fi
			;;
		ulimit)
			if [[ ${IBDMIX_BACKEND:-cpu} == cuda || ${PHYML_BACKEND:-cpu} == cuda || ${TRACE_BACKEND:-cpu} == cuda ]]; then
				echo 'ERROR: CUDA requires the systemd guard locally; ulimit -v cannot safely bound CUDA virtual address space. Enable the systemd user manager and use --resource-guard systemd.' >&2
				return 2
			fi
			current=$(ulimit -S -v); effective=$kb
			[[ $current == unlimited || $current -ge $effective ]] || effective=$current
			ulimit -S -v "$effective" || return 2
			echo "[GU RESOURCE] per-process virtual-memory limit=${effective} KiB; this is not an aggregate RSS limit; jobs=${GU_UNIT_JOBS:-1}" >&2
			;;
	esac
	export GU_RESOURCE_GUARD=$mode
}
