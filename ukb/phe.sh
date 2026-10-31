#!/usr/bin/env bash

this_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# Public directory options use --dir-XXX; legacy spellings remain compatible.
if [[ "${BASH_SOURCE[0]}" == "$0" ]]; then
    source "$this_dir/../0f/cli_dirs.sh"
    script_directory_args "$(realpath -- "${BASH_SOURCE[0]}")" "$@"
    set -- "${SCRIPT_DIRECTORY_ARGS[@]}"
    unset SCRIPT_DIRECTORY_ARGS
fi
# The shared logger uses the production /mnt/d launcher. An isolated checkout
# can still run the read-only audit and --help without that launcher installed.
if [[ -f /mnt/d/scripts/0f/console_run.py ]]; then source "$this_dir/../0f/console.sh"; fi
set -euo pipefail


# 🚩 Command-line help
usage() {
	cat <<'HELP'
Modules:
    Pipeline  Raw preparation -> linked data -> QC -> GWAS phenotypes
    date-check --Y NAME[,NAME...]  Compare saved dates, raw records and UKB reference

Main options:
    --start-step NAME --end-step NAME  Stage range
    --steps CSV  Explicit stages
    --run-r TRUE|FALSE  Run R processing
    --reextract-fields CSV  Re-extract and replace only these multi fields, e.g. p20002,p20008
    -h, --help  Help; details: 0.README.md

Date-check options:
    --Y NAME[,NAME...]          Trait names in common/*.lst; run sequentially
    --phenotype-dir DIR          UKB data root containing Rdata/, rap/, common/
    --common-dir DIR             Definition lists (default: PHEDIR/common)
    --output-dir DIR             Single Y: exact folder; multiple Y: parent folder
                                Default: analysis/ukb/phe/Y
    --gp-records auto|none|FILE   Read every matching GP record; default auto
    --hesin auto|none|FILE       Hospital episode table; default auto
    --hesin-diag auto|none|FILE  Hospital diagnosis table; default auto
    --no-linked-records         Quick audit using existing RDS/AOD only
    --as-of YYYY-MM-DD           Flag future records after this date
    --official-maps FILE        Published code map (bundled by default)
    --rscript FILE              Rscript executable; default Rscript

Example:
    cd /mnt/d/scripts/ukb

    ./phe.sh date-check --Y cvd_afib
    ./phe.sh date-check --Y cvd_afib,cvd_cad
    ./phe.sh date-check --Y cvd_afib --no-linked-records

    ./phe.sh --start-step prep_raw --end-step process_linked --run-r FALSE
    ./phe.sh --start-step data_qc --end-step phe4gwas

    # Optional ADuLT phenotype: requires complete population CIP (trait sex age cip source).
    ADULT_ENABLE=TRUE ADULT_TRAITS=cvd_cad ADULT_CIP_FILE=/mnt/d/data/ukb/common/adu_cip.tsv \
      ./phe.sh --steps phe4gwas

    ./phe.sh --start-step data_qc --end-step phe4gwas
    ./phe.sh  # Default: start at data_qc
HELP
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
MODULE=pipeline
if [[ ${1:-} == date-check ]]; then MODULE=date-check; shift; fi
Y=""
common_dir=""
GP_RECORDS=auto
HESIN=auto
HESIN_DIAG=auto
AS_OF=""
OFFICIAL_MAPS="$this_dir/f/first_occurrence_maps.tsv"
RSCRIPT=${RSCRIPT:-Rscript}
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
dir0=/mnt/d
phedir=""
script_dir=""
gen_dir=""
outdir=""
OUTPUT_DIR_GIVEN=0
phe_R=""
PHENO_IN=pheno.tsv.gz


# 🚩 Command-line arguments
while (($#)); do
	case "$1" in
		--Y) Y=${2:?ERROR: --Y requires a trait name}; shift 2 ;;
		--common-dir | --dir-common) common_dir=${2:?ERROR: --common-dir requires DIR}; shift 2 ;;
		--gp-records) GP_RECORDS=${2:?ERROR: --gp-records requires auto, none, or FILE}; shift 2 ;;
		--hesin) HESIN=${2:?ERROR: --hesin requires auto, none, or FILE}; shift 2 ;;
		--hesin-diag) HESIN_DIAG=${2:?ERROR: --hesin-diag requires auto, none, or FILE}; shift 2 ;;
		--no-linked-records) GP_RECORDS=none; HESIN=none; HESIN_DIAG=none; shift ;;
		--as-of) AS_OF=${2:?ERROR: --as-of requires YYYY-MM-DD}; shift 2 ;;
		--official-maps) OFFICIAL_MAPS=${2:?ERROR: --official-maps requires FILE}; shift 2 ;;
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
		--steps)
			RUN_STEPS=${2:?ERROR: --steps requires CSV}
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
			OUTPUT_DIR_GIVEN=1
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


# 🚩 Shell options, paths, and shared inputs
phedir=${phedir:-$dir0/data/ukb}
script_dir=${script_dir:-$dir0/scripts}
gen_dir=${gen_dir:-/mnt/f/gen/ukb/37}
project_name="$(basename "$this_dir")"
task_name="$(basename "${BASH_SOURCE[0]}" .sh)"
outdir=${outdir:-$dir0/analysis/$project_name/$task_name}

if [[ "$MODULE" == date-check ]]; then
	[[ "$Y" =~ ^[A-Za-z0-9_]+(,[A-Za-z0-9_]+)*$ ]] || { echo 'ERROR: date-check requires --Y NAME[,NAME...] from common/*.lst (no empty entries)' >&2; exit 2; }
	IFS=',' read -r -a date_check_traits <<< "$Y"
	date_check_root=$outdir
	declare -A date_check_seen=()
	common_dir=${common_dir:-$phedir/common}
	command -v "$RSCRIPT" >/dev/null 2>&1 || { echo "ERROR: Rscript not found: $RSCRIPT" >&2; exit 2; }
	if command -v python3 >/dev/null 2>&1; then python_cmd=python3; else python_cmd=python; fi
	[[ -d "$phedir/Rdata" ]] || { echo "ERROR: missing $phedir/Rdata" >&2; exit 2; }
	for Y in "${date_check_traits[@]}"; do
		[[ -z "${date_check_seen[$Y]:-}" ]] || continue
		date_check_seen[$Y]=1
		# Each trait has its own output and cleanup trap. Stop on the first failure.
		(
	outdir=$date_check_root
	if [[ "$OUTPUT_DIR_GIVEN" == 0 ]] || ((${#date_check_traits[@]} > 1)); then outdir="$outdir/$Y"; fi
	mkdir -p -- "$outdir"
	outdir="$(cd "$outdir" && pwd)"
	date_check_work=$(mktemp -d "$outdir/.date-check.XXXXXX")
	cleanup_date_check() {
		local status=$?
		trap - EXIT
		if [[ -n "${date_check_work:-}" && "$date_check_work" == "$outdir"/.date-check.* ]]; then rm -rf -- "$date_check_work"; fi
		exit "$status"
	}
	trap cleanup_date_check EXIT
	echo "Date check: Y=$Y; data=$phedir; report=$outdir"
	"$python_cmd" "$this_dir/f/date_check_io.py" linked --Y "$Y" \
		--phenotype-dir "$phedir" --common-dir "$common_dir" --output-dir "$date_check_work" \
		--gp-records "$GP_RECORDS" --hesin "$HESIN" --hesin-diag "$HESIN_DIAG" --official-maps "$OFFICIAL_MAPS"
	date_r_args=(--Y "$Y" --phenotype-dir "$phedir" --output-dir "$date_check_work")
	[[ -z "$AS_OF" ]] || date_r_args+=(--as-of "$AS_OF")
	"$RSCRIPT" "$this_dir/f/date_check.R" "${date_r_args[@]}"
	"$python_cmd" "$this_dir/f/date_check_report.py" --input-dir "$date_check_work" \
		--output "$date_check_work/date-check.html" --Y "$Y"
	for date_check_file in patients.tsv.gz events.tsv.gz summary.json date-check.html; do
		mv -f -- "$date_check_work/$date_check_file" "$outdir/$date_check_file"
	done
	echo "Date-check report: $outdir/date-check.html"
	exit 0
		)
	done
	exit 0
fi

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


# 🚩 Get raw phenotype data.
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


# 🚩 Fetch FOD fields.
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


# 🚩 Extract FOD fields.
extract_fod() {
	local mapping="$outdir/fetch_fod/fod.each.lst"
	[[ -s "$mapping" ]] || mapping="$phedir/common/fod.each.lst"
	[[ -s "$mapping" ]] || { echo "ERROR: missing fod.each.lst in fetch_fod or common" >&2; exit 1; }
	"$python_cmd" "$this_dir/f/date_check_io.py" fod --input raw/pheno.tab.gz --mapping "$mapping" \
		--date-output fod.date.tab.gz --source-output fod.source.tab.gz
}
run_shell_step extract_fod extract_fod


# 🚩 Extract phenotype data.
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


# 🚩 Extract built-in multi-column data.
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


# 🚩 Process GP and HES files.
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


# 🚩 Call f/phe.R
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
		Rscript "$phe_R"
fi
