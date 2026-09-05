#!/usr/bin/env bash
# AA-aware VCF preparation shared by prep_gen and build.
prepared_record(){
  local c=$1 in=$IN_DIR/chr$1.vcf.gz index
  [[ -s $in ]] || { echo "ERROR: missing phased VCF: $in" >&2; return 1; }
  if [[ -s $in.tbi ]]; then index=$in.tbi; elif [[ -s $in.csi ]]; then index=$in.csi; else echo "ERROR: missing VCF index: $in" >&2; return 1; fi
  python3 "$HERE/arg.cache.py" --file "$in" --file "$index" --text-file "$HERE/arg.prep_vcf.sh" --text-file "$HERE/arg.validate_vcf.py" --value "missing_max=$MISS;mac_min=$MAC;build=$REFGEN_GRCH"
  if [[ $c == X ]]; then
    python3 "$HERE/arg.cache.py" --text-file "$HERE/arg_x.py"
    [[ ! -s ${REFGEN_PFILE_DIR}/chrX.psam ]] || python3 "$HERE/arg.cache.py" --text-file "${REFGEN_PFILE_DIR}/chrX.psam"
  fi
}
prep_one(){
  local c=$1 in=$IN_DIR/chr$1.vcf.gz out meta expected replace_word=FALSE
  out=$(prepared_vcf "$c"); meta=$out.input.meta.tsv; expected=$WORK/chr$c.prepared.meta.tsv
  prepared_record "$c" > "$expected"
  if [[ $REPLACE == FALSE && -s $out && -s $out.tbi && -s $meta ]] && cmp -s "$expected" "$meta"; then
    python3 "$HERE/arg.cache.py" --file "$out" --file "$out.tbi" > "$expected.outputs"
    if cmp -s "$expected.outputs" "$out.outputs.json"; then echo "SKIP verified ARG VCF: $out"; return; fi
  fi
  rm -f "$meta" "$out.outputs.json"
  [[ $REPLACE == FALSE && ! -e $out ]] || replace_word=TRUE
  local -a x_args=()
  if [[ $c == X && -s ${REFGEN_PFILE_DIR}/chrX.psam ]]; then x_args+=(--x-psam "${REFGEN_PFILE_DIR}/chrX.psam"); fi
  bash "$HERE/arg.prep_vcf.sh" --in "$in" --out "$out" --threads "$THREADS" --missing-max "$MISS" --mac-min "$MAC" --tmp-root "$(tmp_for "$c")" --replace "$replace_word" --chr "$c" --build "$REFGEN_GRCH" "${x_args[@]}"
  python3 "$HERE/arg.cache.py" --file "$out" --file "$out.tbi" > "$out.outputs.json.next"
  mv -f "$out.outputs.json.next" "$out.outputs.json"
  mv -f "$expected" "$meta"
}
