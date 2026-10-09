#!/usr/bin/env bash

# Share successful reference checks across workers. Fingerprints include the
# actual source, sentinel list and validator code; changed inputs are rechecked.
gu_build_check_key() {
	python3 - "$1" "$2" "$3" "${CHECK_GRCH_SNP_LIST:-/mnt/d/data/ukb/phe/common/snp.lst}" \
		"$(declare -f check_GRCH phe_zcat phe_check_grch_tabix_rows gu_check_grch_pvar_positions)" <<'PY'
import hashlib, json, sys
from pathlib import Path
source, build, kind, gold, code = sys.argv[1:]
def stamp(name):
    p = Path(name).resolve(strict=True)
    if not p.is_file():
        raise ValueError('build-check caching requires a regular file')
    s = p.stat()
    return [str(p), s.st_dev, s.st_ino, s.st_size, s.st_mtime_ns, s.st_ctime_ns]
request = [1, build, kind, stamp(source), stamp(gold), hashlib.sha256(Path(gold).read_bytes()).hexdigest(), code]
print(hashlib.sha256(json.dumps(request, sort_keys=True).encode()).hexdigest())
PY
}

gu_check_build_uncached() {
	local input=$1 expected=$2 format=$3
	if check_GRCH "$input" "$expected"; then
		return 0
	elif [[ $format == pfile ]] && gu_check_grch_pvar_positions "$input" "$expected"; then
		printf '[GU CHECK] rsID sentinels unavailable; strict coordinate-sentinel fallback passed\n'
		return 0
	fi
	return 1
}

gu_check_build_cached() (
	local input=$1 expected=$2 format=$3 key after detail cache fd
	# Directories/legacy prefixes retain the uncached validation path.
	if [[ ! -f $input ]]; then
		gu_check_build_uncached "$input" "$expected"; return $?
	fi
	cache=${GU_BUILD_CHECK_CACHE_DIR:-/tmp/gu-build-check-cache-$UID}
	mkdir -p "$cache" || return
	key=$(gu_build_check_key "$input" "$expected" "$format") || return
	exec {fd}>"$cache/$key.lock" || return
	flock -x "$fd" || return
	# The source may have changed while this worker was waiting for the lock.
	after=$(gu_build_check_key "$input" "$expected" "$format") || return
	if [[ $key != "$after" ]]; then
		flock -u "$fd"
		gu_check_build_cached "$input" "$expected" "$format"; return $?
	fi
	if [[ -s $cache/$key.pass ]]; then
		printf '[GU CHECK] GRCH cache=HIT build=%s source=%s\n' "$expected" "$input"
		return 0
	fi
	detail=$(mktemp "$cache/$key.XXXXXX") || return
	trap 'rm -f -- "$detail"' EXIT
	printf '[GU CHECK] GRCH cache=MISS build=%s source=%s\n' "$expected" "$input"
	if ! gu_check_build_uncached "$input" "$expected" "$format" >"$detail" 2>&1; then
		cat "$detail"; return 1
	fi
	cat "$detail"
	after=$(gu_build_check_key "$input" "$expected" "$format") || return
	[[ $key == "$after" ]] || { echo 'ERROR: reference changed during build validation' >&2; return 1; }
	printf 'PASS\n' >>"$detail"
	mv -f "$detail" "$cache/$key.pass"
)
