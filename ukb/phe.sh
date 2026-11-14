#!/usr/bin/env bash

# Display standalone help before loading runtimes or scheduling work.
if [[ ${BASH_SOURCE[0]} == "$0" && ( ${1:-} == -h || ${1:-} == --help ) ]]; then
    cat <<'PUBLIC_ENTRY_HELP'
Modules:
    phe         Raw preparation -> linked data -> QC -> GWAS phenotypes

Phe options (phe.sh phe):
    --start-step NAME --end-step NAME  Stage range
    --step NAME[,NAME...]  Selected phe stages (--steps is also accepted)
    --run-r TRUE|FALSE  Run R processing

Shared execution options:
    --submit local|hpc|auto  Optional override; HPC/local is detected automatically
    --write-cmd  Write the .cmd only; --run-cmd FALSE is equivalent
    --queue NAME --slots N --memory-gb N --walltime HHH:MM  LSF resources
    --threads N  Analysis threads [4]; allocation rounded by shared QiMing policy
    --wait TRUE|FALSE  Wait for verified LSF completion [FALSE]
    --dir-jobs DIR  Command/log directory [OUTPUT_DIR/logs]
    UKB_RSCRIPT / UKB_PYTHON  Runtime overrides; --rscript takes precedence
    phe4gwas writes PHEDIR/ukb.phe (default: $DIR_UKB/ukb.phe)
    --reextract-fields CSV  Re-extract and replace only these multi fields, e.g. p20002,p20008
    -h, --help  Show options and examples below

Examples:
    ./phe.sh phe --step phe4gwas --memory-gb 32 --write-cmd
    ./phe.sh phe --step phe4gwas --memory-gb 32
    # Alternative: run the default data_qc -> phe4gwas workflow.
    ./phe.sh phe
PUBLIC_ENTRY_HELP
    exit 0
fi
# Resolve local/HPC paths and job options through the shared 0f/0os.sh table.
this_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)"
source "$this_dir/f/jobs.sh"
ukb_job_init "$@" || exit $?
set -- "${UKB_ARGS[@]}"
source "$shared_dir/cli_dirs.sh"
script_directory_args "$this_dir/phe.sh" "$@"
set -- "${SCRIPT_DIRECTORY_ARGS[@]}"
worker_args=("$@")
set -euo pipefail


# Command-line help
usage() {
	bash "${BASH_SOURCE[0]}" --help
}

bool01() {
	case "${1,,}" in
		true | 1 | yes | y) printf '1\n' ;;
		false | 0 | no | n) printf '0\n' ;;
		*) return 2 ;;
	esac
}

bool_word() {
	case "${1,,}" in
		true | 1 | yes | y) printf 'TRUE\n' ;;
		false | 0 | no | n) printf 'FALSE\n' ;;
		*) return 2 ;;
	esac
}

this_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
MODULE=phe
case ${1:-} in
    phe) shift ;;
    date-check)
        echo "NOTICE: forwarding to the full AOD workflow; open date-check.html (tab 01: date consistency)." >&2
        shift
        legacy_args=("$@")
        (( UKB_EXECUTE == 0 )) || legacy_args=(--execute "${legacy_args[@]}")
        (( UKB_WRITE_ONLY == 0 )) || legacy_args=(--write-cmd "${legacy_args[@]}")
        export UKB_THREADS UKB_SLOTS UKB_MEMORY_GIB UKB_QUEUE UKB_WALLTIME UKB_PYTHON UKB_RSCRIPT
        [[ ! -v UKB_JOB_DIR ]] || export UKB_JOB_DIR
        exec bash "$this_dir/aod.sh" "${legacy_args[@]}"
        ;;
    ''|-*) ;; # Existing option-only calls retain the default phe module.
    *) echo "ERROR: unknown module '$1'; use phe; date analysis is in aod.sh" >&2; exit 2 ;;
esac
# Generated workers always record the explicit module, including legacy calls.
worker_args=("$MODULE" "$@")
RSCRIPT=$UKB_RSCRIPT
run_rscript() {
	# System R uses system GDAL/PROJ. Conda activation can export paths to
	# an incompatible database; let the linked system libraries find theirs.
	if [[ $(readlink -f "$(command -v "$RSCRIPT")") == /usr/bin/Rscript ]]; then
		env -u PROJ_LIB -u PROJ_DATA -u GDAL_DATA "$RSCRIPT" "$@"
	else
		"$RSCRIPT" "$@"
	fi
}
START_STEP=data_qc
END_STEP=""
RUN_STEPS=""
RUN_R=1
REEXTRACT_FIELDS=""
DATS="vip bbc img met death job"
ICD_SOURCES=icd10,icd9,opcs4,srd
SAVE_AOD_TEXT=TRUE
SAVE_AOD_LIST=FALSE
MAKE_PHEWAS=FALSE
dir0=$DIR_DAT
phedir=""
script_dir=""
gen_dir=""
outdir=""
phe_R=""
PHENO_IN=pheno.tsv.gz


# Command-line arguments
while (($#)); do
	case "$1" in
		--rscript) RSCRIPT=${2:?ERROR: --rscript requires executable}; shift 2 ;;
		--reextract-fields) REEXTRACT_FIELDS=${2:?ERROR: --reextract-fields requires CSV}; shift 2 ;;
		--start-step)
			START_STEP=${2:?ERROR: --start-step requires STEP}
			shift 2
			;;
		--end-step)
			END_STEP=${2:?ERROR: --end-step requires STEP}
			shift 2
			;;
		--step | --steps)
			RUN_STEPS=${2:?ERROR: --step requires NAME or CSV}
			shift 2
			;;
		--run-r)
			RUN_R=$(bool01 "${2:?ERROR: --run-r requires TRUE or FALSE}") || {
				echo "ERROR: --run-r requires TRUE or FALSE" >&2
				exit 2
			}
			shift 2
			;;
		--datasets)
			DATS=${2:?ERROR: --datasets requires LIST}
			shift 2
			;;
		--icd-sources)
			ICD_SOURCES=${2:?ERROR: --icd-sources requires CSV}
			shift 2
			;;
		--save-aod-text)
			SAVE_AOD_TEXT=$(bool_word "${2:?ERROR: --save-aod-text requires TRUE or FALSE}") || {
				echo "ERROR: --save-aod-text requires TRUE or FALSE" >&2
				exit 2
			}
			shift 2
			;;
		--save-aod-list)
			SAVE_AOD_LIST=$(bool_word "${2:?ERROR: --save-aod-list requires TRUE or FALSE}") || {
				echo "ERROR: --save-aod-list requires TRUE or FALSE" >&2
				exit 2
			}
			shift 2
			;;
		--make-phewas)
			MAKE_PHEWAS=$(bool_word "${2:?ERROR: --make-phewas requires TRUE or FALSE}") || {
				echo "ERROR: --make-phewas requires TRUE or FALSE" >&2
				exit 2
			}
			shift 2
			;;
		--data-root | --dir-data-root)
			dir0=${2:?ERROR: --data-root requires DIR}
			shift 2
			;;
		--phenotype-dir | --dir-phenotype)
			phedir=${2:?ERROR: --phenotype-dir requires DIR}
			shift 2
			;;
		--genotype-dir | --dir-genotype)
			gen_dir=${2:?ERROR: --genotype-dir requires DIR}
			shift 2
			;;
		--output-dir | --dir-output)
			outdir=${2:?ERROR: --output-dir requires DIR}
			shift 2
			;;
		--script-dir | --dir-script)
			script_dir=${2:?ERROR: --script-dir requires DIR}
			shift 2
			;;
		--r-script)
			phe_R=${2:?ERROR: --r-script requires FILE}
			shift 2
			;;
		--phenotype-input)
			PHENO_IN=${2:?ERROR: --phenotype-input requires FILE}
			shift 2
			;;
		-h | --help)
			usage
			exit 0
			;;
		*)
			echo "Unknown option: $1" >&2
			usage >&2
			exit 2
			;;
	esac
done
DATS=${DATS//,/ }


# Shell options, paths, and shared inputs
phedir=${phedir:-$DIR_UKB}
script_dir=${script_dir:-$SCRIPTS_ROOT}
gen_dir=${gen_dir:-$DIR_GEN/ukb/37}
project_name="$(basename "$this_dir")"
task_name="$(basename "${BASH_SOURCE[0]}" .sh)"
outdir=${outdir:-$ANALYSIS_ROOT/$project_name/$task_name}

phe_R=${phe_R:-$this_dir/f/phe.R}
if command -v python3 >/dev/null 2>&1; then python_cmd=python3; else python_cmd=python; fi

shell_steps=(prep_raw fetch_fod extract_fod extract_pheno extract_multi process_linked)
r_steps=(data_qc vip_phe ses move death srd icd fod_ref gp_prep gp biom hla_pca pgs merge drug phe4gwas)

all_steps=(all "${shell_steps[@]}" "${r_steps[@]}")
contains() {
	local x="$1"
	shift
	for y in "$@"; do [[ "$x" == "$y" ]] && return 0; done
	return 1
}

if ! contains "$START_STEP" "${all_steps[@]}"; then
	echo "ERROR: unknown START_STEP: $START_STEP" >&2
	usage >&2
	exit 1
fi
if [[ -n "$END_STEP" ]] && ! contains "$END_STEP" "${all_steps[@]}"; then
	echo "ERROR: unknown END_STEP: $END_STEP" >&2
	usage >&2
	exit 1
fi
has_r_run_steps=0
r_run_steps=""
if [[ -n "$RUN_STEPS" ]]; then
	IFS=',' read -ra _run_steps_check <<<"$RUN_STEPS"
	for _st in "${_run_steps_check[@]}"; do
		_st="${_st// /}"
		[[ -z "$_st" ]] && continue
		if ! contains "$_st" "${all_steps[@]}"; then
			echo "ERROR: unknown RUN_STEPS entry: $_st" >&2
			usage >&2
			exit 1
		fi
		if [[ "$_st" == "all" ]]; then
			has_r_run_steps=1
			r_run_steps="all"
		elif contains "$_st" "${r_steps[@]}"; then
			has_r_run_steps=1
			[[ "$r_run_steps" == "all" ]] || r_run_steps="${r_run_steps:+$r_run_steps,}$_st"
		fi
	done
fi

# All stage names have been checked before a command file is submitted.
ukb_dispatch phe "$outdir" "$this_dir/phe.sh" "${worker_args[@]}"
ukb_worker_runtime

idx_of() {
	local x="$1"
	shift
	local i=0
	for y in "$@"; do
		[[ "$x" == "$y" ]] && {
			echo "$i"
			return 0
		}
		i=$((i + 1))
	done
	echo -1
}
in_run_steps() {
	[[ -z "$RUN_STEPS" ]] && return 1
	[[ ",${RUN_STEPS// /,}," == *",all,"* || ",${RUN_STEPS// /,}," == *",$1,"* ]]
}
should_run_shell() {
	local step="$1"
	if [[ -n "$RUN_STEPS" ]]; then
		in_run_steps "$step"
		return $?
	fi
	[[ "$START_STEP" == "all" ]] && return 0
	contains "$START_STEP" "${shell_steps[@]}" || return 1
	local i s e
	i=$(idx_of "$step" "${shell_steps[@]}")
	s=$(idx_of "$START_STEP" "${shell_steps[@]}")
	if [[ -n "$END_STEP" ]] && contains "$END_STEP" "${shell_steps[@]}"; then e=$(idx_of "$END_STEP" "${shell_steps[@]}"); else e=$((${#shell_steps[@]} - 1)); fi
	[[ $i -ge $s && $i -le $e ]]
}
run_shell_step() {
	local step="$1"
	shift
	if should_run_shell "$step"; then
		echo
		echo "========== RUN SHELL STEP: $step =========="
		"$@"
		echo "========== DONE SHELL STEP: $step =========="
	else echo "Skip shell step: $step"; fi
}
zcat_any() {
	local f="$1"
	if [[ "$f" == *.gz ]]; then zcat "$f"; else cat "$f"; fi
}

merge_events() {
	local id="$1" cols="$2"
	awk -F'\t' -v OFS=$'\t' -v id="$id" -v cols="$cols" '
        BEGIN {
            print id, "events"
            n_req = split(cols, req, /[[:space:]]+/)
        }

        NR == 1 {
            for (i = 1; i <= NF; i++) {
                name[i] = $i
                col_idx[$i] = i
            }

            id_col = (id in col_idx) ? col_idx[id] : 1
            for (i = 1; i <= n_req; i++) {
                if (req[i] != "" && req[i] in col_idx && col_idx[req[i]] != id_col) {
                    use[++n_use] = col_idx[req[i]]
                }
            }
            if (n_use == 0) {
                for (i = 1; i <= NF; i++) {
                    if (i != id_col) use[++n_use] = i
                }
            }
            next
        }

        {
            row_id = $id_col
            events = ""
            for (i = 1; i <= n_use; i++) {
                idx = use[i]
                value = $idx
                if (value == "" || value == "NA") continue
                event = name[idx] "=" value
                events = (events == "" ? event : events "; " event)
            }

            if (row_id != previous_id) {
                if (previous_id != "") print ""
                printf "%s%s%s", row_id, OFS, events
                previous_id = row_id
            } else if (events != "") {
                printf " | %s", events
            }
        }

        END {
            if (previous_id != "") print ""
        }
    '
}

cd "$phedir/rap"


# Get raw phenotype data.
prep_raw() {
	mkdir -p raw
	local in=${PHENO_IN:-pheno.tsv.gz}
	[[ -s "$in" ]] || [[ -s "raw/pheno.tsv.gz" ]] || {
		echo "ERROR: cannot find pheno.tsv.gz or raw/pheno.tsv.gz" >&2
		exit 1
	}
	[[ -s "$in" ]] || in="raw/pheno.tsv.gz"
	zcat_any "$in" | awk -F'\t' 'NR==1 {for (i=1; i<=NF; i++) colname[i]=$i; ncols=NF; next} {for (i=1; i<=ncols; i++) {if ($i == "") missing[i]++}} END {for (i=1; i<=ncols; i++) {print colname[i], missing[i]+0}}' >pheno.missing
	zcat_any "$in" | sed -e 's/^\t/NA\t/; s/\t\t/\tNA\t/g; s/\t\t/\tNA\t/g; s/\t$/\tNA/' | gzip -f >raw/pheno.tab.gz
	ln -sf raw/pheno.tab.gz pheno.tab.gz
	# Reading only the gzip header avoids a zcat SIGPIPE under pipefail when
	# head exits early on a large cohort export.
	"$python_cmd" - raw/pheno.tab.gz >pheno.id <<'PY'
import gzip
import sys
with gzip.open(sys.argv[1], "rt", encoding="utf-8-sig") as inp:
    sys.stdout.write(inp.readline().rstrip("\r\n").replace("\t", "\n") + "\n")
PY
	nl -ba pheno.id | awk 'BEGIN{OFS="\t"}{print $1,$2}' >pheno.id.2col
	sed 's/_.*//' pheno.id | sort | uniq >pheno.id.uniq
	wc -l pheno.id.uniq
}
run_shell_step prep_raw prep_raw


# Fetch FOD fields.
fetch_fod() {
	mkdir -p "$outdir/fetch_fod"
	: >fod.lst.tmp
	for page in {2401..2417}; do
		fn="https://biobank.ndph.ox.ac.uk/ukb/label.cgi?id=$page"
		if curl --output /dev/null --silent --head --fail "$fn"; then
			wget "$fn" -O - >>fod.lst.tmp
		else
			echo "Page $page not found"
		fi
	done
	cat fod.lst.tmp | grep -E '<tr class="row_odd"' | grep 'field.cgi' | sed -E 's/<[^>]*>//g' | sed -E 's/^([0-9]+)([A-Za-z]+.*)/\1 \2/' | cut -d " " -f 1,3 | sed 's/^/p/; s/ /\t/' >"$outdir/fetch_fod/fod.each.lst"
}
run_shell_step fetch_fod fetch_fod


# Extract FOD fields.
extract_fod() {
	local mapping="$outdir/fetch_fod/fod.each.lst"
	[[ -s "$mapping" ]] || mapping="$phedir/common/fod.each.lst"
	[[ -s "$mapping" ]] || { echo "ERROR: missing fod.each.lst in fetch_fod or common" >&2; exit 1; }
	"$python_cmd" "$this_dir/f/date_check_io.py" fod --input raw/pheno.tab.gz --mapping "$mapping" \
		--date-output fod.date.tab.gz --source-output fod.source.tab.gz
}
run_shell_step extract_fod extract_fod


# Extract phenotype data.
extract_pheno() {
	local dat
	for dat in $DATS; do
		[[ -s "$phedir/common/$dat.lst" ]] || {
			echo "Skip $dat: missing $phedir/common/$dat.lst"
			continue
		}
		echo RUN: "$dat"
		"$python_cmd" "$this_dir/f/date_check_io.py" fields --input raw/pheno.tab.gz \
			--mapping "$phedir/common/$dat.lst" --output "$dat.tab.gz"
	done
}
run_shell_step extract_pheno extract_pheno


# Extract built-in multi-column data.
extract_multi() {
	local -a Arr1=(p87 p20002 p20003 p20008 p22009 p22182 p40002 p41270 p41271 p41272 p41280 p41281 p41282)
	local -a Arr2=(srd.time srd med srd.year pca hla death_icd10 icd10.code icd9.code opcs4.code icd10.date icd9.date opcs4.date)
	local -a requested=()
	local i field name
	if [[ -n "$REEXTRACT_FIELDS" ]]; then
		IFS=',' read -ra requested <<<"${REEXTRACT_FIELDS// /}"
		((${#requested[@]})) || { echo "ERROR: --reextract-fields requires field IDs" >&2; return 2; }
		for field in "${requested[@]}"; do
			contains "$field" "${Arr1[@]}" || {
				echo "ERROR: --reextract-fields has unsupported multi field: $field" >&2
				return 2
			}
		done
	fi
	for i in "${!Arr1[@]}"; do
		field=${Arr1[$i]}
		name=${Arr2[$i]}
		if [[ -n "$REEXTRACT_FIELDS" ]]; then
			contains "$field" "${requested[@]}" || continue
		elif [[ -s $name.tab.gz ]]; then
			echo DONE: "$field" "$name"
			continue
		fi
		echo RUN: "$field" "$name"
		if [[ -z "$REEXTRACT_FIELDS" ]] && ! awk -v wanted="$field" '{sub(/_.*/, ""); if($0 == wanted) found=1} END {exit !found}' pheno.id; then
			echo "WARN: no columns selected for $field"
			continue
		fi
		# Explicit re-extraction reads the current raw header, not a potentially
		# stale pheno.id. Missing requested fields fail; existing outputs cannot
		# silently stand in for a successful extraction.
		"$python_cmd" "$this_dir/f/date_check_io.py" fields --input raw/pheno.tab.gz --fields "$field" --output "$name.tab.gz"
	done
}
run_shell_step extract_multi extract_multi


# Process GP and HES files.
process_linked() {
	local gp_clinical=""
	[[ -s raw/gp_clinical.tsv.gz ]] && gp_clinical=raw/gp_clinical.tsv.gz
	[[ -z "$gp_clinical" && -s raw/gp_clinical.tsv ]] && gp_clinical=raw/gp_clinical.tsv
	if [[ -n "$gp_clinical" ]]; then
		zcat_any "$gp_clinical" | LC_ALL=C sort -t$'\t' -k 1,1n -k 3,3 |
			sed '1s/data_provider/dp/; 1s/event_dt/dt/; 1s/read_/r/g; 1s/value/v/g' |
			awk -F'\t' -v OFS=$'\t' 'BEGIN { print "eid\tevents" } NR==1 {for (i=1; i<=NF; i++) name[i] = $i; next} { ev = ""; for (i=2; i<=NF; i++) {v = $i; if (v == "") continue; key = name[i]; p = key "=" v; ev = (ev=="" ? p : ev "; " p)} if ($1 != prev) { if (prev != "") print ""; printf "%s\t%s", $1, ev; prev = $1} else if (ev != "") {printf " | %s", ev }} END { if (prev != "") print "" }' |
			gzip -f >gp_clinical.NEW.tsv.gz
	else
		echo "Skip GP clinical: raw/gp_clinical.tsv[.gz] not found"
	fi

	local gp_script=""
	[[ -s raw/gp_scripts.tsv.gz ]] && gp_script=raw/gp_scripts.tsv.gz
	[[ -z "$gp_script" && -s raw/gp_scripts.tsv ]] && gp_script=raw/gp_scripts.tsv
	[[ -z "$gp_script" && -s raw/gp_script.tsv.gz ]] && gp_script=raw/gp_script.tsv.gz
	[[ -z "$gp_script" && -s raw/gp_script.tsv ]] && gp_script=raw/gp_script.tsv
	if [[ -n "$gp_script" ]]; then
		zcat_any "$gp_script" | LC_ALL=C sort -t$'\t' -k 1,1n -k 3,3 |
			sed '1s/data_provider/dp/; 1s/issue_date/dt/; 1s/read_/r/g; 1s/_code//g; 1s/_name//g; 1s/quantity/qt/' |
			merge_events eid "dp dt r2 bnf dmd drug qt" |
			gzip -f >gp_scripts.NEW.tsv.gz
	else
		echo "Skip GP scripts: raw/gp_script(s).tsv[.gz] not found"
	fi

	local hes=""
	[[ -s raw/hesin_diag.tsv.gz ]] && hes=raw/hesin_diag.tsv.gz
	[[ -z "$hes" && -s raw/hesin_diag.tsv ]] && hes=raw/hesin_diag.tsv
	if [[ -n "$hes" ]]; then
		zcat_any "$hes" | cut -f 3- | LC_ALL=C sort -t$'\t' -k 1,1n -k 2,3n |
			sed '1s/_index//g; 1s/diag_//g; 1s/level/lv/g' |
			merge_events eid "ins arr lv icd9 icd9_nb icd10 icd10_nb" |
			gzip -f >hesin_diag.NEW.tsv.gz
	else
		echo "Skip HES diagnosis: raw/hesin_diag.tsv[.gz] not found"
	fi
}
run_shell_step process_linked process_linked


# Call f/phe.R
if [[ "$RUN_R" == "1" ]]; then
	if [[ -n "$RUN_STEPS" && "$has_r_run_steps" != "1" ]]; then
		echo "Skip R: RUN_STEPS has no R steps"
		exit 0
	fi
	if contains "$START_STEP" "${r_steps[@]}" || [[ "$START_STEP" == "all" ]]; then
		r_start="$START_STEP"
	else
		r_start=data_qc
	fi
	r_end=${END_STEP:-}
	if [[ -n "$r_end" ]] && ! contains "$r_end" "${r_steps[@]}"; then r_end=""; fi
	r_run_env="$RUN_STEPS"
	if [[ -n "$RUN_STEPS" ]]; then
		r_run_env="$r_run_steps"
		[[ "$r_run_env" == "all" ]] && r_run_env=""
	fi
	echo
	echo "========== RUN R: f/phe.R START_STEP=$r_start END_STEP=${r_end:-<default>} RUN_STEPS=${r_run_env:-<none>} =========="
	echo "ICD_SOURCES=$ICD_SOURCES SAVE_AOD_TEXT=$SAVE_AOD_TEXT SAVE_AOD_LIST=$SAVE_AOD_LIST MAKE_PHEWAS=$MAKE_PHEWAS"
	DIR0="$dir0" PHEDIR="$phedir" UKB_OUT="$outdir" SCRIPT_DIR="$script_dir" GEN_DIR="$gen_dir" START_STEP="$r_start" END_STEP="$r_end" RUN_STEPS="$r_run_env" \
		ICD_SOURCES="$ICD_SOURCES" SAVE_AOD_TEXT="$SAVE_AOD_TEXT" SAVE_AOD_LIST="$SAVE_AOD_LIST" MAKE_PHEWAS="$MAKE_PHEWAS" \
		run_rscript "$phe_R"
fi

