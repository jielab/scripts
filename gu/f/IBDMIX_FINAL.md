# IBDmix filters in `gu.sh final`

`normalize_results.py` calls `ibdmix_final_filter.py` before inserting IBDmix
segments into the Shiny database. Original chromosome results and `raw/` files
remain intact. Filtered files, exclusion BEDs and JSON audits are cached under
`final/normalize/ibdmix_filters/`. The database's `ibdmix_filter_runs` table and
`normalize/summary/ibdmix_filters.tsv` report which filters were applied.

## Denisovan refinement (enabled by default on GRCh37 autosomes)

1. Construct the union of original Altai calls across all populations, at
   LOD >= 4 and length >= 50 kb. Read `raw/`, not the Neanderthal final files
   already masked by African Denisova calls.
2. Separately for Denisova and Denisova25, construct regions carried by >=30%
   of tested ESN/GWD/LWK/MSL/YRI participants. Merge each person's intervals
   before counting. Use all tested members, including zero-call members, as
   the denominator. This implements a pointwise pooled-five-population carrier
   fraction, not a phased allele frequency. It is an explicit operational
   interpretation of the paper's regional-frequency wording.
3. Subtract both masks from that Denisovan reference's current final segments;
   discard residual pieces shorter than 50 kb. LOD remains the original
   parent call's score; fragment site counts are cleared, not fabricated.

Altai/Denisova is the paper reference pair. Denisova25 uses the same Altai
mask as an extension. Vindija and Chagyr calls remain available separately.
chrX is unchanged and marked as outside the Cell autosomal comparison.
Missing original calls or denominators are errors, not permission to publish
unfiltered Denisovan results as filtered.

## Highest-0.1% Altai derived-allele windows

The public paper describes this optional refinement but does not specify
sufficient window construction parameters in the inspected methods/workflow.
Do not estimate this mask from LOD scores, segment length, carrier frequency,
or by changing parameters to match the paper's reported totals.

Supply a GRCh37, 0-based half-open exclusion BED and its `.json` sidecar:

```bash
GU_IBDMIX_DAF_MASK=/absolute/path/altai_daf_top0.1pct.bed ./gu.sh final
```

Sidecar `/absolute/path/altai_daf_top0.1pct.bed.json`:

```json
{
  "genome_build": "GRCh37",
  "reference": "Altai",
  "percentile": 99.9,
  "source": "Author mask citation or reproducible method/source description"
}
```

Only autosomal coordinates are accepted. The mask is applied to Altai only;
applying an Altai DAF mask to other archaic individuals is not equivalent.
Each audit records the mask checksum and supplied provenance. An invalid
mask aborts final. If no mask is supplied, the DAF step is explicitly marked
`not_applied_missing_author_mask`, a warning is printed, and Shiny displays
that the strict DAF refinement has not been applied. Denisovan filtering still
runs. This is intentionally not described as an exact paper reproduction.

## Shiny comparison

The population map and Cell 2020 comparison table now use Altai-only exact
per-person autosomal interval unions, including tested zero-call people.
They no longer combine three Neanderthal references. The density overview
retains per-lineage unions for the extended reference panel. X remains outside
the Cell comparison. Input fingerprints, code and mask provenance invalidate
filtered caches. Database replacement remains atomic.

## Source

Chen et al., Cell (2020), DOI 10.1016/j.cell.2020.01.012, STAR Methods,
“Refining Neanderthal Callset by Using Denisovan Sequences as a Negative Control”.
https://rilab.ucdavis.edu/pdfs/IBDmix-cell.pdf
https://github.com/PrincetonUniversity/IBDmix/tree/main/j.cell.2020.01.012-workflow

The author's exclusion BED has not been located as of this implementation.
Exact window reconstruction and segment-by-segment reproduction remain open.
