# Shiny GWAS block browser

Start with `./gwas_compare.sh shiny`; its `--help` contains a complete height
example, including `trait`, `gwas_files` and `gwas_csv`. Defaults and CLI settings
remain in `gwas_compare.sh`. This directory contains:

- `app.R`: launcher, configuration, reference setup and cache worker entry point.
- `ui.R`: Overview / Blockview navigation, shared IGV header and panels.
- `server.R`: viewport state, interactions and asynchronous workers.
- `data.R`, `plot.R`: GWAS chromosome caches, BED lookup and Manhattan tracks.
- `ld.py`, `ld.R`: PRS-CSx HDF5 reader and population LD heatmaps.
- `ref.R`, `www/`: local FASTA range requests, JavaScript, CSS and bundled IGV.

Keep the launching terminal open; Ctrl+C stops the server. The default URL is
http://127.0.0.1:3841. Use the shell entry point to supply the inputs and defaults;
`ui.R` and `server.R` are sourced by `app.R`.

## Inputs and thinning

Inputs require standardized `CHR`, `POS`, `P`, with optional `SNP` and `LOG10P`.
Build metadata comes from `file.gz.grch` or the project's `qc/trait.grch`. All
inputs must have the same source build. Labels default to the filename, ignoring
`.thin.gz`; ancestry is its final dot-separated component. `--races CSV`
overrides that label explicitly.

When a full GWAS path is supplied, an adjacent `.thin.gz` takes precedence if it
exists. The viewer checks its generation metadata and refuses stale sidecars.
It still accepts full GWAS when no thin file exists. Generate sidecars with
`gwas_format.sh thin` or `thin,mplot`; see that script's complete batch example
for `4grid`, `main`, `met`, and `prot`.

`format` defaults to `--thin TRUE --thin-chr-max 10000`. Independently,
`--hm3 FALSE` keeps the standardized full `.gz`; `--hm3 TRUE` keeps the union of HM3
and P < 1e-3 in that file (threshold: `--p-hm3 1e-3`). Both modes support
`--thin TRUE` or `--thin FALSE`. Reference paths use `--hm3-file` and `--hm3-pos`.
The thin helper independently selects
HM3 SNPs (rsID or build-specific position) and all P < 1e-3 candidates, calls
the actual shared `0phe.f.R::thinP0` to prune P > 1e-3 background, and imposes the
per-chromosome total cap. If signals alone exceed that cap, strongest signals
win; `LOG10P` breaks underflow ties. The sidecar preserves original columns and
writes only `[gwas].thin.gz` and `[gwas].thin.gz.tbi`. Each requested thin run
regenerates both files to honor current inputs and settings, and removes legacy
thin `.grch`, `.n.tsv`, and `.meta.rds` sidecars after success. Build metadata
comes from the project’s `qc/trait.grch`; Shiny checks source modification time
for stale thin inputs. `--thin FALSE` disables generation; it does not delete
previous sidecars. Statistical analysis modules continue to use `[gwas].gz`, whose filtering is controlled by `--hm3`.
Mplot uses the sidecar for drawing, and the full GWAS for hit reports/cis flags.

At display time, `--max-points 12000` also limits each track using position-bin
peaks plus background samples. Track captions distinguish actual **shown** points
from all input variants in the current region. Zoom queries the current input
cache again; a thin input cannot restore variants removed during thinning.
Track y scales are independent, and finite `LOG10P` preserves P underflow.

## Overview and Blockview

Overview contains the GRCh dropdown and All / 1–22 / X chromosome bar. IGV,
region entry, zoom and pan share one viewport. Hover anywhere within a track's
BED-covered region to see `block: AFR 123`. Double-click opens Blockview:
IGV, all GWAS tracks, then EUR / AFR / EAS / SAS / AMR LD heatmaps.

The initial region is the selected block plus **one preceding and one following
block**, on that same chromosome. At chromosome ends it contains fewer blocks.
The selected block is shaded. Zoom and IGV changes update both GWAS and LD;
“所选 block ＋ 相邻 blocks” restores the initial region. Returning to Overview
restores its last viewport. Changing GRCh clears the selected block because
BED IDs are build-specific.

Only `[race].[grch].bed` supplies boundaries. Missing files/chromosomes produce
no boundaries and no inferred IDs; HIS is not automatically relabeled AMR.
BED intervals are zero-based half-open: START < GWAS POS <= END. IDs are
one-based row numbers after sorting the complete BED by chromosome/start/end,
across all chromosomes, independently per ancestry and build. All-chromosome
views omit boundary lines but still support block lookup. GRCh37 EAS/SAS use the
published ASN data. Reference provenance and excluded upstream intervals are
recorded under the block directory. `python3 f/ld_blocks.py` reproduces those
references without overwriting different existing BEDs.

## LD reference

`--ld-dir /mnt/i/refLD/csx` contains `ldblk_1kg_EUR`, `ldblk_1kg_AFR`,
`ldblk_1kg_EAS`, `ldblk_1kg_SAS`, `ldblk_1kg_AMR` and
`snpinfo_mult_1kg_hm3`. The SNP information file provides GRCh37 coordinates.
For GRCh38, reference SNP positions are actually mapped through hg19ToHg38;
unmapped/non-primary positions are excluded. LD values remain those of the
same 1KG reference genotypes.

The reader follows the [official PRS-CSx HDF5 schema](https://github.com/getian107/PRScsx/blob/master/parse_genet.py):
each `blk_*` contains `snplist` and `ldblk`. HDF5 block IDs do not correspond to
BED block IDs. Matrices are matched by SNP ID and source chromosome, normalized
by their diagonal and shown as r². Missing SNPs and cross-HDF5-block pairs are
white (unknown), never fabricated zeros. Allele sign flips do not change r².

All populations use the same physical coordinate axes and SNP subset. At most
`--ld-max-snps 300` SNPs are sampled evenly along the ordered reference positions
(range 2–1000). Captions give available and displayed counts. Missing population
files and autosome-only reference coverage are reported in their panels.

## Caches and dependencies

Outputs include `shiny.inputs.json`, `cache/grch37/track_NN/`, corresponding
GRCh38 caches, `cache/ld/` SNP indexes/region matrices, and worker logs. Build
switches convert GWAS viewing copies using the configured UCSC chain; source
GWAS are unchanged and all tracks switch together once ready. Keep separate
output directories for different GWAS collections.

R packages: `shiny`, `bslib`, `plotly`, `htmlwidgets`, `data.table`, `jsonlite`,
`processx`. `--ld-python` defaults to `~/anaconda3/bin/python` and requires
`h5py`, `numpy`, `pandas`. Shell tools: Python 3, gzip/pigz, flock, liftOver,
bgzip/tabix. IGV JS and its license are bundled in `www/vendor/`. Indexed
`GRCH37.fasta` / `GRCH38.fasta` and optional `glist.37.bed` / `glist.38.bed`
supply local sequence and gene tracks.

Checks from the project root: `Rscript --vanilla tests/bplot.R`,
`Rscript --vanilla tests/thin.R`, and `~/anaconda3/bin/python tests/ld.py`.
`tests/bplot_browser.cjs` checks a running six-height example; it requires Chrome
and `playwright-core`. Set `NODE_PATH` for an external package installation,
`BPLOT_URL` for a different server, and optionally `BPLOT_SHOT` for a screenshot.
