#!/usr/bin/env bash
# Shared thin completion checks for the coordinator and workers.
# Fingerprints use paths, sizes and nanosecond mtimes; no full GWAS scan.

# Plot-only runs consume an existing artifact, independently of the producer's
# path-sensitive cache marker. Never regenerate data as a plotting side effect.
gwas_thin_plot_only() {
  [[ ",$1," == *,mplot,* && "$1" != all && ! ",$1," =~ ,(format|thin|liftover), ]]
}

gwas_thin_plot_ready() {
  local input="$1" output="$2"
  [[ -s "$input" && -s "$output" && -s "$output.tbi" && ! "$input" -nt "$output" && ! "$output" -nt "$output.tbi" ]] || return 1
  tabix -l "$output" >/dev/null 2>&1
}

gwas_thin_signature() {
  local input="$1" output="$2" grch="$3" cap="$4" hm3="$5" pos="$6" script="$7" phe="$8" f
  local facts
  facts=$(
    printf 'thin_state_v1\t%s\t%s\n' "$grch" "$cap"
    for f in "$input" "$output" "$output.tbi" "$hm3" "$script" "${script%/*}/thin.awk" "$phe"; do
      [[ -s "$f" ]] || exit 1
      stat -Lc '%n|%s|%y' -- "$f" || exit 1
    done
    if [[ -n "$pos" ]]; then
      [[ -s "$pos" ]] || exit 1
      stat -Lc '%n|%s|%y' -- "$pos" || exit 1
    fi
  ) || return 1
  printf '%s\n' "$facts" | sha256sum | cut -d ' ' -f 1
}

gwas_thin_complete() {
  local input="$1" output="$2" grch="$3" cap="$4" hm3="$5" pos="$6" script="$7" phe="$8" legacy="${9:-}"
  local signature saved f
  [[ -s "$input" && -s "$output" && -s "$output.tbi" && "$output" -nt "$input" && ! "$output" -nt "$output.tbi" ]] || return 1
  if [[ -e "$output.done" ]]; then
    signature=$(gwas_thin_signature "$input" "$output" "$grch" "$cap" "$hm3" "$pos" "$script" "$phe") || return 1
    saved=$(cat -- "$output.done") || return 1
    [[ "$saved" == "$signature" ]] || return 1
  else
    # Older runs recorded completion only after plotting succeeded. Adopt those
    # artifacts only when that record postdates BOTH the thin data and index.
    [[ -s "$legacy" && "$legacy" -nt "$output" && "$legacy" -nt "$output.tbi" ]] || return 1
    for f in "$hm3" "$pos"; do
      [[ -z "$f" ]] && continue
      [[ -s "$f" && ! "$f" -nt "$output" ]] || return 1
    done
    awk -F '\t' -v input="$input" -v output="$output" -v grch="$grch" -v cap="$cap" '
      $1=="thin_mode"&&$2=="TRUE"{t=1}
      $1=="thin_chr_max"&&$2==cap{c=1}
      $1=="grch"&&$2==grch{g=1}
      $1=="gwas_input"&&$2==input{i=1}
      $1=="display_input"&&$2==output{o=1}
      END{exit !(t&&c&&g&&i&&o)}' "$legacy" || return 1
  fi
  tabix -l "$output" >/dev/null 2>&1
}
