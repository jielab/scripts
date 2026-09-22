#!/usr/bin/env bash
# Prefer the activated environment; explicitly allow GRID_RSCRIPT overrides.
grid_select_r(){
  local packages=$1 candidate
  local -a runtime
  if [[ -n ${GRID_RSCRIPT:-} ]]; then
    GRID_R=("$GRID_RSCRIPT")
    # The system R must use its own startup configuration and package libraries.
    if [[ $GRID_RSCRIPT -ef /usr/bin/Rscript ]]; then
      GRID_R=(env -u R_ENVIRON_USER -u R_LIBS_USER "$GRID_RSCRIPT")
    fi
    "${GRID_R[@]}" -e "p<-strsplit('$packages',',',fixed=TRUE)[[1]];stopifnot(all(vapply(p,requireNamespace,logical(1),quietly=TRUE)))" || return
    return
  fi
  for candidate in "$(command -v Rscript || true)" /usr/bin/Rscript; do
    [[ -n $candidate && -x $candidate ]] || continue
    runtime=("$candidate")
    if [[ $candidate -ef /usr/bin/Rscript ]]; then
      runtime=(env -u R_ENVIRON_USER -u R_LIBS_USER "$candidate")
    fi
    if "${runtime[@]}" -e "p<-strsplit('$packages',',',fixed=TRUE)[[1]];stopifnot(all(vapply(p,requireNamespace,logical(1),quietly=TRUE)))" >/dev/null 2>&1; then GRID_R=("${runtime[@]}");return;fi
  done
  echo "ERROR: No usable R runtime; install packages: $packages, or set GRID_RSCRIPT" >&2
  return 1
}
