#!/usr/bin/env bash
# Analysis and library-check routing. Workers and scheduling retain their APIs.
set -euo pipefail
HERE="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)"
domain=${1:?Missing analysis module}
shift
selection=""
full_help=FALSE
short_help=FALSE
has_scope=FALSE
forward=()

select_action() {
    local requested=$1
    if [[ -n $selection && $selection != "$requested" ]]; then
        echo "Choose only one $domain action: --$selection and --$requested cannot be combined." >&2
        exit 2
    fi
    selection=$requested
}

while (($#)); do
    case "$domain:$1" in
        validate:--human|validate:--reads|validate:--plants|report:--corrected-groups)
            select_action "${1#--}"; shift ;;
        refine:--profile-only)
            select_action profile-only; shift ;;
        check:--status|check:--verify|check:--references|check:--selftest|check:--environment|check:--benchmark)
            select_action "${1#--}"; shift ;;
        *:--help-full) full_help=TRUE; shift ;;
        *:-h|*:--help) short_help=TRUE; forward+=("$1"); shift ;;
        validate:--scope|validate:--scope=*) has_scope=TRUE; forward+=("$1"); shift ;;
        *) forward+=("$1"); shift ;;
    esac
done

if [[ $domain == validate && $selection == reads && $has_scope == TRUE ]]; then
    echo "--scope is a human-validation option; --reads uses --limit (0 scans all)." >&2
    exit 2
fi
if [[ -z $selection && $short_help == TRUE && $full_help != TRUE ]]; then
    exec bash "$HERE/f/usage.sh" "$domain"
fi

worker=""
backend=""
defaults=()
case "$domain:${selection:-default}" in
    report:corrected-groups) worker=group_summary ;;
    validate:plants) worker=plant_validate ;;
    report:default) worker=split_summary; defaults=(--require-complete) ;;
    validate:default|validate:human) worker=human_validate; defaults=(--scope candidates) ;;
    validate:reads) worker=sample_diagnose ;;
    refine:default) worker=species_refine ;;
    refine:profile-only) worker=sylph_profile; defaults=(profile) ;;
    check:status) backend=status ;;
    check:default|check:verify) backend=audit ;;
    check:references) backend=audit-reference ;;
    check:selftest) backend=selftest ;;
    check:environment) backend=doctor ;;
    check:benchmark) backend=compare ;;
    *) echo "Unknown analysis action: $domain/$selection" >&2; exit 2 ;;
esac

args=()
# qc_entry and public frontends expect --help-full first.
[[ $full_help == TRUE ]] && args+=(--help-full)
if [[ -n $worker ]]; then
    exec bash "$HERE/f/qc_entry.sh" "$worker" "${args[@]}" "${defaults[@]}" "${forward[@]}"
fi
exec bash "$HERE/metag.sh" "$backend" "${args[@]}" "${forward[@]}"
