#!/usr/bin/env bash
set -euo pipefail



usage() {
    cat <<'EOF'
Usage:
  ./phe.sh [options]

Options:
  --start-step STEP       Start from STEP. Default: data_qc.
  --end-step STEP         Stop after STEP.
  --steps CSV             Run only the listed steps.
  --run-r TRUE|FALSE      Run the R stage. Default: TRUE.
  --datasets LIST         Shell extraction datasets. Comma- or space-separated.
  --icd-sources CSV       ICD sources. Default: icd10,icd9,opcs4,srd.
  --save-aod-text TRUE|FALSE  Save aodtxt_* columns. Default: TRUE.
  --save-aod-list TRUE|FALSE  Save list-column aod_* values. Default: FALSE.
  --make-phewas TRUE|FALSE    Build PheWAS.rds. Default: FALSE.
  --data-root DIR         Data/analysis root. Default: /mnt/d.
  --phenotype-dir DIR     UKB phenotype directory.
  --genotype-dir DIR      UKB genotype directory.
  --output-dir DIR        Pipeline output directory.
  --script-dir DIR        Shared scripts root.
  --r-script FILE         R workflow entry point.
  --phenotype-input FILE  Input for extract_pheno. Default: pheno.tsv.gz.
  -h, --h, --help        Show this help message.

Shell steps:
  prep_raw fetch_fod extract_fod extract_pheno extract_multi process_linked

R steps:
  data_qc vip_phe ses move death srd icd audit_i67 fod_ref gp_prep gp biom hla_pca pgs merge ckm audit_ckm audit_dates drug phe4gwas

Association analyses (forest and PheWAS) have moved to ./assoc.sh.

Examples:
  ./phe.sh --datasets vip --steps extract_pheno
  ./phe.sh --save-aod-list TRUE --steps icd,merge
  ./phe.sh --start-step prep_raw --end-step process_linked --run-r FALSE
EOF
}

bool01() {
    case "${1,,}" in
        true|1|yes|y) printf '1\n' ;;
        false|0|no|n) printf '0\n' ;;
        *) return 2 ;;
    esac
}

bool_word() {
    case "${1,,}" in
        true|1|yes|y) printf 'TRUE\n' ;;
        false|0|no|n) printf 'FALSE\n' ;;
        *) return 2 ;;
    esac
}

this_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
START_STEP=data_qc
END_STEP=""
RUN_STEPS=""
RUN_R=1
DATS="vip bbc img.ckm met death job"
ICD_SOURCES=icd10,icd9,opcs4,srd
SAVE_AOD_TEXT=TRUE
SAVE_AOD_LIST=FALSE
MAKE_PHEWAS=FALSE
dir0=/mnt/d
phedir=""
script_dir=""
gen_dir=""
outdir=""
phe_R=""
PHENO_IN=pheno.tsv.gz

while (( $# )); do
    case "$1" in
        --start-step) START_STEP=${2:?ERROR: --start-step requires STEP}; shift 2 ;;
        --end-step) END_STEP=${2:?ERROR: --end-step requires STEP}; shift 2 ;;
        --steps) RUN_STEPS=${2:?ERROR: --steps requires CSV}; shift 2 ;;
        --run-r) RUN_R=$(bool01 "${2:?ERROR: --run-r requires TRUE or FALSE}") || { echo "ERROR: --run-r requires TRUE or FALSE" >&2; exit 2; }; shift 2 ;;
        --datasets) DATS=${2:?ERROR: --datasets requires LIST}; shift 2 ;;
        --icd-sources) ICD_SOURCES=${2:?ERROR: --icd-sources requires CSV}; shift 2 ;;
        --save-aod-text) SAVE_AOD_TEXT=$(bool_word "${2:?ERROR: --save-aod-text requires TRUE or FALSE}") || { echo "ERROR: --save-aod-text requires TRUE or FALSE" >&2; exit 2; }; shift 2 ;;
        --save-aod-list) SAVE_AOD_LIST=$(bool_word "${2:?ERROR: --save-aod-list requires TRUE or FALSE}") || { echo "ERROR: --save-aod-list requires TRUE or FALSE" >&2; exit 2; }; shift 2 ;;
        --make-phewas) MAKE_PHEWAS=$(bool_word "${2:?ERROR: --make-phewas requires TRUE or FALSE}") || { echo "ERROR: --make-phewas requires TRUE or FALSE" >&2; exit 2; }; shift 2 ;;
        --data-root) dir0=${2:?ERROR: --data-root requires DIR}; shift 2 ;;
        --phenotype-dir) phedir=${2:?ERROR: --phenotype-dir requires DIR}; shift 2 ;;
        --genotype-dir) gen_dir=${2:?ERROR: --genotype-dir requires DIR}; shift 2 ;;
        --output-dir) outdir=${2:?ERROR: --output-dir requires DIR}; shift 2 ;;
        --script-dir) script_dir=${2:?ERROR: --script-dir requires DIR}; shift 2 ;;
        --r-script) phe_R=${2:?ERROR: --r-script requires FILE}; shift 2 ;;
        --phenotype-input) PHENO_IN=${2:?ERROR: --phenotype-input requires FILE}; shift 2 ;;
        -h|--help) usage; exit 0 ;;
        *) echo "Unknown option: $1" >&2; usage >&2; exit 2 ;;
    esac
done
DATS=${DATS//,/ }


# 🚩 Shell options, paths, and shared inputs
phedir=${phedir:-$dir0/data/ukb/phe}
script_dir=${script_dir:-$dir0/scripts}
gen_dir=${gen_dir:-$dir0/data/ukb/gen}
project_name="$(basename "$this_dir")"
task_name="$(basename "${BASH_SOURCE[0]}" .sh)"
outdir=${outdir:-$dir0/analysis/$project_name/$task_name}

privacy_cleanup_completed_run=FALSE
cleanup_ukb_individual_outputs() {
    local status=$? f
    trap - EXIT
    if (( status != 0 )) || [[ "$privacy_cleanup_completed_run" != TRUE ]]; then
        echo "Phenotype pipeline did not reach full completion (exit=$status); preserving participant-level audit data for resume." >&2
        exit "$status"
    fi
    case "$outdir" in ""|/|/mnt|/mnt/|/mnt/d|/mnt/d/|/mnt/d/analysis|/mnt/d/analysis/)
        echo "PRIVACY ERROR: unsafe UKB_OUT: '$outdir'" >&2; exit 70 ;;
    esac
    if [[ -d "$outdir" ]]; then
        echo "Privacy cleanup: removing phenotype audit rows containing UKB eids"
        while IFS= read -r -d '' f; do rm -f -- "$f"; done < <(find "$outdir" -type f \( \
            -name 'I67_trajectory_audit.csv' -o -name 'dates_audit.csv' -o \
            -name 'dates_audit.rds' \) -print0)
    fi
    exit "$status"
}
trap cleanup_ukb_individual_outputs EXIT
phe_R=${phe_R:-$this_dir/f/phe.R}
if command -v python3 >/dev/null 2>&1; then python_cmd=python3; else python_cmd=python; fi
privacy_cleanup_completed_run=TRUE

shell_steps=(prep_raw fetch_fod extract_fod extract_pheno extract_multi process_linked)
r_steps=(data_qc vip_phe ses move death srd icd audit_i67 fod_ref gp_prep gp biom hla_pca pgs merge ckm audit_ckm audit_dates drug phe4gwas)

all_steps=(all "${shell_steps[@]}" "${r_steps[@]}")
contains() { local x="$1"; shift; for y in "$@"; do [[ "$x" == "$y" ]] && return 0; done; return 1; }

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
    IFS=',' read -ra _run_steps_check <<< "$RUN_STEPS"
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

idx_of() { local x="$1"; shift; local i=0; for y in "$@"; do [[ "$x" == "$y" ]] && { echo "$i"; return 0; }; i=$((i+1)); done; echo -1; }
in_run_steps() { [[ -z "$RUN_STEPS" ]] && return 1; [[ ",${RUN_STEPS// /,}," == *",all,"* || ",${RUN_STEPS// /,}," == *",$1,"* ]]; }
should_run_shell() {
    local step="$1"
    if [[ -n "$RUN_STEPS" ]]; then in_run_steps "$step"; return $?; fi
    [[ "$START_STEP" == "all" ]] && return 0
    contains "$START_STEP" "${shell_steps[@]}" || return 1
    local i s e
    i=$(idx_of "$step" "${shell_steps[@]}"); s=$(idx_of "$START_STEP" "${shell_steps[@]}")
    if [[ -n "$END_STEP" ]] && contains "$END_STEP" "${shell_steps[@]}"; then e=$(idx_of "$END_STEP" "${shell_steps[@]}"); else e=$((${#shell_steps[@]}-1)); fi
    [[ $i -ge $s && $i -le $e ]]
}
run_shell_step() {
    local step="$1"; shift
    if should_run_shell "$step"; then echo; echo "========== RUN SHELL STEP: $step =========="; "$@"; echo "========== DONE SHELL STEP: $step =========="; else echo "Skip shell step: $step"; fi
}
zcat_any() { local f="$1"; if [[ "$f" == *.gz ]]; then zcat "$f"; else cat "$f"; fi; }

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
    [[ -s "$in" ]] || [[ -s "raw/pheno.tsv.gz" ]] || { echo "ERROR: cannot find pheno.tsv.gz or raw/pheno.tsv.gz" >&2; exit 1; }
    [[ -s "$in" ]] || in="raw/pheno.tsv.gz"
    zcat_any "$in" | awk -F'\t' 'NR==1 {for (i=1; i<=NF; i++) colname[i]=$i; ncols=NF; next} {for (i=1; i<=ncols; i++) {if ($i == "") missing[i]++}} END {for (i=1; i<=ncols; i++) {print colname[i], missing[i]+0}}' > pheno.missing
    zcat_any "$in" | sed -e 's/^\t/NA\t/; s/\t\t/\tNA\t/g; s/\t\t/\tNA\t/g; s/\t$/\tNA/' | gzip -f > raw/pheno.tab.gz
    ln -sf raw/pheno.tab.gz pheno.tab.gz
    zcat raw/pheno.tab.gz | head -1 | tr '\t' '\n' > pheno.id
    nl -ba pheno.id | awk 'BEGIN{OFS="\t"}{print $1,$2}' > pheno.id.2col
    sed 's/_.*//' pheno.id | sort | uniq > pheno.id.uniq; wc -l pheno.id.uniq
}
run_shell_step prep_raw prep_raw


# 🚩 Fetch FOD fields.
fetch_fod() {
    : > fod.lst.tmp
    for page in {2401..2417}; do
        fn="https://biobank.ndph.ox.ac.uk/ukb/label.cgi?id=$page"
        if curl --output /dev/null --silent --head --fail "$fn"; then
            wget "$fn" -O - >> fod.lst.tmp
        else
            echo "Page $page not found"
        fi
    done
    cat fod.lst.tmp | grep -E '<tr class="row_odd"' | grep 'field.cgi' | sed -E 's/<[^>]*>//g' | sed -E 's/^([0-9]+)([A-Za-z]+.*)/\1 \2/' | cut -d " " -f 1,3 | sed 's/^/p/; s/ /\t/' > "$phedir/common/fod.each.lst"
}
run_shell_step fetch_fod fetch_fod


# 🚩 Extract FOD fields.
extract_fod() {
    local dat=fod
    awk 'BEGIN{print "eid"}{print $1}' "$phedir/common/$dat.each.lst" > "$dat.id"
    cols=$(grep -nwf "$dat.id" pheno.id.2col | sed 's/:/ /' | awk '{print $1}' | tr '\n' ',' | sed 's/,$//')
    [[ -n "$cols" ]] || { echo "ERROR: no columns found for $dat.each.lst" >&2; exit 1; }
    zcat raw/pheno.tab.gz | cut -f "$cols" | gzip -f > "$dat.date.tab.gz"
}
run_shell_step extract_fod extract_fod


# 🚩 Extract phenotype data.
extract_pheno() {
    for dat in $DATS; do
        [[ -s "$phedir/common/$dat.lst" ]] || { echo "Skip $dat: missing $phedir/common/$dat.lst"; continue; }
        echo RUN: "$dat"
        awk 'BEGIN{print "eid"}{print $1}' "$phedir/common/$dat.lst" > "$dat.id"
        dos2unix "$dat.id" >/dev/null 2>&1 || true
        sort "$dat.id" | uniq -d || true
        grep -wf "$dat.id" pheno.id.uniq | wc -l
        grep -nwf "$dat.id" pheno.id.2col | sed 's/:/ /; s/ /\t/g' > "$dat.id.tmp1" || true
        "$python_cmd" "$script_dir/0f/0join_file.py" \
          --inputs "$dat.id.tmp1,TAB,1" "$phedir/common/$dat.lst,TAB,0" \
          --output "$dat.id.tmp2"
        cols=$(awk '{if ($6 == "all" || ($3 ~ /_i0$/ || $3 ~ /_i0_a0$/)) printf $1","}' "$dat.id.tmp2" | sed 's/,$//')
        [[ -n "$cols" ]] || { echo "WARN: no columns selected for $dat"; continue; }
        zcat raw/pheno.tab.gz | cut -f "$cols" | gzip -f > "$dat.tab.gz"
    done
}
run_shell_step extract_pheno extract_pheno


# 🚩 Extract built-in multi-column data.
extract_multi() {
    Arr1=(p87 p20002 p20003 p22009 p22182 p40002 p41270 p41271 p41272 p41280 p41281 p41282)
    Arr2=(srd.time srd med pca hla death_icd10 icd10.code icd9.code opcs4.code icd10.date icd9.date opcs4.date)
    for i in ${!Arr1[@]}; do
        field=${Arr1[$i]}; name=${Arr2[$i]}
        if [[ -e $name.tab.gz ]]; then echo DONE: "$field" "$name"; continue; else echo RUN: "$field" "$name"; fi
        cols=$(fgrep -nw "$field" pheno.id.2col | awk -F ":" 'BEGIN{print 1}{print $1}' | tr '\n' ',' | sed 's/,$//')
        [[ -n "$cols" ]] || { echo "WARN: no columns selected for $field"; continue; }
        zcat raw/pheno.tab.gz | cut -f "$cols" | gzip -f > "$name.tab.gz"
    done
}
run_shell_step extract_multi extract_multi


# 🚩 Process GP and HES files.
process_linked() {
    local gp_clinical=""
    [[ -s raw/gp_clinical.tsv.gz ]] && gp_clinical=raw/gp_clinical.tsv.gz
    [[ -z "$gp_clinical" && -s raw/gp_clinical.tsv ]] && gp_clinical=raw/gp_clinical.tsv
    if [[ -n "$gp_clinical" ]]; then
        zcat_any "$gp_clinical" | LC_ALL=C sort -t$'\t' -k 1,1n -k 3,3 | \
            sed '1s/data_provider/dp/; 1s/event_dt/dt/; 1s/read_/r/g; 1s/value/v/g' | \
            awk -F'\t' -v OFS=$'\t' 'BEGIN { print "eid\tevents" } NR==1 {for (i=1; i<=NF; i++) name[i] = $i; next} { ev = ""; for (i=2; i<=NF; i++) {v = $i; if (v == "") continue; key = name[i]; p = key "=" v; ev = (ev=="" ? p : ev "; " p)} if ($1 != prev) { if (prev != "") print ""; printf "%s\t%s", $1, ev; prev = $1} else if (ev != "") {printf " | %s", ev }} END { if (prev != "") print "" }' | \
            gzip -f > gp_clinical.NEW.tsv.gz
    else
        echo "Skip GP clinical: raw/gp_clinical.tsv[.gz] not found"
    fi

    local gp_script=""
    [[ -s raw/gp_scripts.tsv.gz ]] && gp_script=raw/gp_scripts.tsv.gz
    [[ -z "$gp_script" && -s raw/gp_scripts.tsv ]] && gp_script=raw/gp_scripts.tsv
    [[ -z "$gp_script" && -s raw/gp_script.tsv.gz ]] && gp_script=raw/gp_script.tsv.gz
    [[ -z "$gp_script" && -s raw/gp_script.tsv ]] && gp_script=raw/gp_script.tsv
    if [[ -n "$gp_script" ]]; then
        zcat_any "$gp_script" | LC_ALL=C sort -t$'\t' -k 1,1n -k 3,3 | \
            sed '1s/data_provider/dp/; 1s/issue_date/dt/; 1s/read_/r/g; 1s/_code//g; 1s/_name//g; 1s/quantity/qt/' | \
            merge_events eid "dp dt r2 bnf dmd drug qt" | \
            gzip -f > gp_scripts.NEW.tsv.gz
    else
        echo "Skip GP scripts: raw/gp_script(s).tsv[.gz] not found"
    fi

    local hes=""
    [[ -s raw/hesin_diag.tsv.gz ]] && hes=raw/hesin_diag.tsv.gz
    [[ -z "$hes" && -s raw/hesin_diag.tsv ]] && hes=raw/hesin_diag.tsv
    if [[ -n "$hes" ]]; then
        zcat_any "$hes" | cut -f 3- | LC_ALL=C sort -t$'\t' -k 1,1n -k 2,3n | \
            sed '1s/_index//g; 1s/diag_//g; 1s/level/lv/g' | \
            merge_events eid "ins arr lv icd9 icd9_nb icd10 icd10_nb" | \
            gzip -f > hesin_diag.NEW.tsv.gz
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
