#!/usr/bin/env bash
# Sourced by gwas_compare.sh; CLI validation only. Project and statistical
# functions remain in gwas_compare_project.py, gwas_compare.R and gwas_ldsc.py.

compare_die() { echo "ERROR: $*" >&2; exit 2; }

compare_need_arg_value() {
  local option="$1" value="${2-}"
  [[ -n "$value" && "$value" != --* ]] || compare_die "Missing value for $option"
}

gwas_compare_validate_options() {
  local key
  local -a flags=(check_only)
  if [[ -n "$project_dir" && -n "$gwas_files" ]] || [[ -z "$project_dir" && -z "$gwas_files" ]]; then
    compare_die 'Specify one of --dir-gwas/--project-dir or --gwas-files'
  fi
  for key in grch require_grch; do
    case "${!key}" in
      ''|37|38) ;;
      *) compare_die "--${key//_/-} must be 37 or 38" ;;
    esac
  done
  if [[ -n "$grch" && -n "$require_grch" && "$grch" != "$require_grch" ]]; then
    compare_die '--grch and --require-grch must agree'
  fi

  case "$module" in
    shiny)
      [[ "$ld_max_snps" =~ ^[0-9]+$ ]] && (( 10#$ld_max_snps >= 2 && 10#$ld_max_snps <= 1000 )) || compare_die "--ld-max-snps must be 2–1000"
      flags+=(launch_browser)
      [[ "$port" =~ ^[0-9]+$ ]] && (( 10#$port > 0 && 10#$port <= 65535 )) || compare_die '--port must be 1-65535'
      [[ "$max_points" =~ ^[0-9]+$ ]] && (( 10#$max_points >= 100 && 10#$max_points <= 100000 )) || compare_die '--max-points must be 100-100000'
      awk -v p="$p_threshold" 'BEGIN{exit !(p ~ /^[0-9]*[.]?[0-9]+([eE][-+]?[0-9]+)?$/ && p+0>0 && p+0<=1)}' ||
        compare_die '--p-threshold must be a number in (0,1]'
      ;;
    compare)
      flags+=(compare_beta compare_eaf)
      case "$significant" in first|either|both) ;; *) compare_die 'Invalid --significant: use first, either or both' ;; esac
      awk -v p="$p_threshold" 'BEGIN{exit !(p ~ /^[0-9]*[.]?[0-9]+([eE][-+]?[0-9]+)?$/ && p+0>0 && p+0<=1)}' ||
        compare_die '--p-threshold must be a number in (0,1]'
      ;;
    ldsc)
      flags+=(run_munge run_h2 run_rg run)
      case "$missing_n" in error|skip) ;; *) compare_die '--missing-n must be error or skip' ;; esac
      if [[ -n "$sample_size" ]]; then
        awk -v n="$sample_size" 'BEGIN{exit !(n ~ /^[0-9]*[.]?[0-9]+([eE][-+]?[0-9]+)?$/ && n+0>0)}' ||
          compare_die '--N must be a positive number'
      fi
      ;;
  esac
  for key in "${flags[@]}"; do
    case "${!key}" in
      [Tt][Rr][Uu][Ee]) printf -v "$key" '%s' TRUE ;;
      [Ff][Aa][Ll][Ss][Ee]) printf -v "$key" '%s' FALSE ;;
      *) compare_die "--${key//_/-} must be TRUE or FALSE" ;;
    esac
  done
}
