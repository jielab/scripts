#!/usr/bin/env bash
# UKB workers share one cgroup. PLINK --memory only budgets its workspace;
# Python, metadata, compression and the downstream caller also need RAM.
gu_plan_ukb_resources() {
	local cap=${1^^} requested_jobs=$2 number cap_mb reserve_mb usable_mb max_plink_mb memory_mb jobs
	[[ $cap =~ ^([1-9][0-9]{0,5})(G|GB|GIB)$ ]] || {
		echo 'ERROR: --memory-cap requires a positive whole number of GiB.' >&2
		return 2
	}
	number=${BASH_REMATCH[1]}
	cap_mb=$((10#$number * 1024))
	# Reserve 25% of the group (at least 1 GiB), plus 2 GiB per worker
	# outside PLINK. This is scheduling headroom, not a per-process RSS cap.
	reserve_mb=$((cap_mb / 4))
	((reserve_mb >= 1024)) || reserve_mb=1024
	usable_mb=$((cap_mb - reserve_mb))
	max_plink_mb=$((usable_mb - 2048))
	if [[ -n ${GU_PREP_MEMORY_MB:-} ]]; then
		[[ $GU_PREP_MEMORY_MB =~ ^[0-9]{1,9}$ ]] || {
			echo 'ERROR: GU_PREP_MEMORY_MB must be an integer number of MiB.' >&2
			return 2
		}
		memory_mb=$((10#$GU_PREP_MEMORY_MB))
	else
		memory_mb=8192
		((memory_mb <= max_plink_mb)) || memory_mb=$max_plink_mb
	fi
	if ((memory_mb < 640 || memory_mb > max_plink_mb)); then
		echo "ERROR: UKB PLINK workspace must be at least 640 MiB and fit with process overhead under --memory-cap $cap; use a cap of at least 4G and unset or reduce GU_PREP_MEMORY_MB." >&2
		return 2
	fi
	jobs=$((usable_mb / (memory_mb + 2048)))
	((jobs <= requested_jobs)) || jobs=$requested_jobs
	GU_UKB_EFFECTIVE_JOBS=$jobs
	GU_PREP_MEMORY_MB=$memory_mb
	export GU_PREP_MEMORY_MB
	if [[ ${GU_CMD_WORKER:-0} != 1 ]]; then
		echo "[GU MEMORY] UKB cap=$cap requested_jobs=$requested_jobs effective_jobs=$jobs plink_memory_mb=$memory_mb reserve_mb=$reserve_mb worker_overhead_mb=2048"
	fi
}
