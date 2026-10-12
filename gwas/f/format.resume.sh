#!/usr/bin/env bash

# Resume means reuse existing output names. Do not open GWAS data, indexes or
# receipts, compare timestamps/settings, or invalidate an old rsID fingerprint.
# --replace TRUE explicitly requests recomputation instead.
declare -A gwas_resume_done=()

# Read each output directory once, then use in-memory lookups in the main loop.
# Variables below are supplied by format.sh.
# shellcheck disable=SC2154
gwas_format_resume_prepare() {
	local names_file="$1" completed_file trait module
	local -a modules
	gwas_resume_done=()
	[[ "$replace" != TRUE && "$liftOver" != TRUE && "$delete_raw" != TRUE ]] || return 0
	IFS=',' read -r -a modules <<<"$step"
	for module in "${modules[@]}"; do
		case "$module" in format | thin | magma | lead | mplot) ;; *) return 0 ;; esac
	done
	completed_file=$(mktemp "$dir_cmd/gwas_resume.XXXXXXXX") || return 1
	log "START [resume-check] batch file-name lookup"
	if ! python3 - "$names_file" "$dir_out" "$category" "$dir_clean_arg" "$dir_magma" \
		"$step" "$thin" "$rsid_mode" "$write_sig" "$add_panel" "$jobs" >"$completed_file" <<'PY'
import os
import sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

(names_file, project, category, clean, magma, step, thin, rsid,
 write_sig, panel, jobs) = sys.argv[1:]
root = Path(project)
modules = set(step.split(','))

def names(path):
    try:
        # listdir reads names only: no per-file stat, data reads or checksums.
        return set(os.listdir(path))
    except (FileNotFoundError, NotADirectoryError):
        return set()

plots = names(root / 'mplot') if 'mplot' in modules else set()
states = names(root / '.project' / category / 'mplot') if 'mplot' in modules else set()
flags = names(root / '.project' / category / 'mplot' / 'flag') if 'mplot' in modules else set()
shared_magma = names(magma) if magma else None

def complete(trait):
    gwas_dir = Path(clean) if clean else root / category / trait / 'gwas'
    files = names(gwas_dir)
    if trait + '.gz' not in files or not {trait + '.gz.tbi', trait + '.gz.csi'} & files:
        return False
    if 'format' in modules and rsid == 'TRUE':
        if trait + '.rsid.done' not in names(gwas_dir.parent / 'qc'):
            return False
    if thin == 'TRUE' and modules & {'format', 'thin', 'mplot'}:
        required = {trait + '.thin.gz', trait + '.thin.gz.tbi'}
        if modules & {'format', 'thin'}:
            required.add(trait + '.thin.gz.done')
        if not required <= files:
            return False
    if 'magma' in modules or ('mplot' in modules and panel == 'magma'):
        gene_files = shared_magma if shared_magma is not None else names(gwas_dir.parent / 'magma')
        if not {'magma.done', trait + '.genes.out', trait + '.genes.raw'} <= gene_files:
            return False
    if 'lead' in modules:
        # Both markers are written only after their phase succeeds, including
        # valid empty selections that have no .jma.cojo/.clumps result table.
        if not {trait + '.awk.snp', trait + '.clump.done', trait + '.cojo.done'} <= files:
            return False
        if {'clump', 'cojo'} & files:
            return False
    if 'mplot' in modules:
        if not {trait + '.png', '0flag.tsv'} <= plots:
            return False
        if trait + '.state' not in states or trait + '.tsv' not in flags:
            return False
        if write_sig == 'TRUE' and trait + '.sig.txt' not in plots:
            return False
    return True

with open(names_file) as source:
    traits = [line.rstrip('\n') for line in source if line.strip()]
with ThreadPoolExecutor(max_workers=min(int(jobs), 8)) as pool:
    for trait, done in zip(traits, pool.map(complete, traits)):
        if done:
            print(trait)
PY
	then
		rm -f -- "$completed_file"
		return 1
	fi
	while IFS= read -r trait; do
		[[ -z "$trait" ]] || gwas_resume_done["$trait"]=1
	done <"$completed_file"
	rm -f -- "$completed_file"
}

gwas_format_resume_complete() {
	[[ -n "${gwas_resume_done[$1]:-}" ]]
}
