# PHOLE: proteome black-hole framework

PHOLE treats thousands of protein GWAS as one genome-ordered effect graph. It
does **not** read Manhattan PNG files. The numeric unit is a matrix whose rows
are source genes/loci on the coding genome (the houses/mothers) and whose
columns are measured proteins (the children).

## Upstream MAGMA blocker found on 2026-08-31

Do **not** start the formal PHOLE run from the MAGMA files currently being
generated. The PHOLE code runs mechanically, but the present MAGMA gene results
have two upstream input defects that can remove real cis signals.

1. The shared SNP-to-gene annotation was built from one protein GWAS and then
   reused for all proteins. The older flat cache records `source_gwas=A1BG`; the
   newer resource-tagged cache records `source_gwas=ALPI`. A protein-specific
   SNP that is absent from that one source GWAS is absent from the shared
   annotation even when it exists in another protein GWAS.
2. `gwas_post.sh` currently keeps only `P > 0` when it writes the MAGMA P-value
   input. Extremely strong P values stored as numerical zero are therefore
   discarded instead of reconstructed from `LOG10P` and floored to a small
   positive value.

The first defect is directly visible at ABCA2. The A1BG GWAS has no variants in
the ABCA2 coding interval, so the shared annotation contains no ABCA2/Entrez-20
entry. The ABCA2 GWAS itself has 116 variants in that interval, including 13 at
`P < 5e-8` (minimum `P=3.70e-15`), but ABCA2 is absent from
`ABCA2.genes.out`. Across the 2,958 PROT assays, 2,861 HGNC symbols map to the
NCBI gene-location file; 211 target coding genes are absent from the A1BG-based
cache and 206 are absent from the ALPI-based cache.

The second defect is also material. Within the existing ±100 kb cis files for
the alphabetically first ten proteins, AAMDC has 299 and ABO has 300 records
with `P=0`; their `LOG10P` values reach approximately 1,576 and 370,
respectively. These are the strongest signals, not invalid rows.

Before rerunning MAGMA:

1. Build the shared annotation from a canonical or union PROT SNP-coordinate
   universe, not from the first job acquiring the cache lock. Record and verify
   the universe/hash in `annotation.meta.tsv`.
2. For `P=0`, reconstruct from `LOG10P` where available and cap at a positive
   floor such as `1e-300`; never silently drop the row.
3. Invalidate both existing `window_0_0` annotation caches and all MAGMA outputs
   derived from them, then rerun MAGMA from clean outputs.
4. Verify that each assay's mapped target gene is present in its own
   `*.genes.out`, and compare its target-gene P value with the direct coding
   interval and ±100 kb cis SNP minima before PHOLE.
5. Keep `--window 0,0` as the source-gene-specific primary definition, then run
   a modest annotation-window sensitivity analysis rather than broadening the
   primary graph immediately.
6. The current summary-P run uses MAGMA's default `snp-wise=mean`. Add
   `multi=snp-wise` (mean plus top-SNP) as a sensitivity model for ABL1-like
   genes where only a small fraction of in-gene SNPs is strongly associated.

The primary unit is the complete gene (the mother), not a lead SNP. PHOLE reads
the full, completed `*.genes.out` files made by `gwas_post.sh magma`:

```text
<gwas-root>/common/<protein>/magma/<protein>.genes.out
```

For example, `common/ABO/magma/ABO.genes.out` supplies 17,540 whole-gene
mother-to-ABO-child scores. The local MAGMA files are not perfectly identical
gene grids: A1BG has 17,555 rows and ALPI contains genes absent from A1BG.
`prepare.py` therefore takes the union of the selected, completed MAGMA gene
IDs, orders it by chromosome/position, aligns every protein to that grid, and
keeps unavailable cells as `NaN`. The earlier
`<gwas-root>/<protein>/common/magma/` layout remains supported.

## Scientific model

For complete source gene `g` and protein child `j`, PHOLE writes the observed
standardized association score as

```text
W[g,j] = row[g] + child[j] + C[g,j] cis_alpha[j]
         + sum_h alpha[g,h] beta[h,j] + gamma[g,j] + epsilon[g,j]
```

- `C[g,j]` identifies the coding-gene cis window of protein `j`.
- `cis_alpha[j]` is the child's fitted cis excess.
- `alpha[g,h]` connects a genomic source to latent proteomic program `h`.
- `beta[h,j]` connects that program to protein child `j`.
- `alpha * beta` is the low-rank, mediation-inspired shared trans component.
- `gamma[g,j]` is a sparse direct/idiosyncratic residual edge.
- `epsilon[g,j]` is unexplained noise.

This is deliberately simpler and more identifiable than immediately fitting
the dynamical expression

```text
W = A (I - B)^(-1) + G + E.
```

The latter separates biological protein-to-protein propagation from genetic
injection, but summary statistics alone do not identify `A` and `B` without
strong anchors or interventions. PHOLE v1 therefore estimates the identifiable
product `alpha * beta`. It can later be upgraded when signed, allele-harmonized
effects, colocalization, perturbation, or longitudinal data are available.

This is inspired by DANDELION's use of a simple product structure, but it is not
the DANDELION/DACT test and is not a causal-mediation analysis. DANDELION combines
trans-regulatory and disease-association P values under a composite null. PHOLE
instead learns a multi-protein graph and tests out-of-sample prediction.

## Why the full MAGMA matrix is the primary input

Using a lead SNP as the graph node would split the mother into genetic pieces.
Clump/COJO variants remain useful evidence underneath a gene node, but they are
not the PHOLE node. Using only significant clump/COJO leads would also encode an
unreported cell as zero.
That is wrong: an absent lead can be a weak effect, low power, LD mismatch, or a
selection artifact. MAGMA supplies a nearly common annotated gene universe for
every protein; PHOLE takes the observed union and preserves trait-specific
missing cells rather than treating them as zero.

Important limitation: MAGMA `ZSTAT` is a gene-level association statistic. Its
sign is **not** the direction of protein-abundance change. With MAGMA, PHOLE is
therefore an association-strength graph. A future signed-effect adapter should
use common SNPs/loci, harmonize the effect allele across every protein, and
extract marginal `BETA/SE` from all full GWAS rather than merging protein-specific
conditional COJO `bJ` estimates.

## Cold-start prediction

The requested training/testing experiment is a true protein-level cold-start
test:

1. Fit row baselines, cis excess, and the `alpha * beta` graph in the training
   proteins. By default the training count is
   `floor(0.80 * number_of_retained_proteins)`.
2. Each training protein receives a learned `beta` embedding.
3. Predict the embedding of each held-out protein using only its coding-gene
   position and multi-scale codenome distance kernels.
4. Reconstruct its entire source-gene profile and evaluate **trans rows only**.

The kernel uses physical genomic distance within a chromosome plus a small
cross-chromosome background weight. This is the simplest test of the boulevard
hypothesis. A Hi-C/3D embedding can later replace this kernel without changing
the `alpha * beta + gamma` core.

Rank is chosen only inside the training set. Held-out proteins are not used to
choose it. The 80/20 split is calculated **after** MAGMA target-gene filtering.
`--n-train INT` remains available only as an explicit override of
`--train-frac 0.80`.

## Protein-level MAGMA gate

`prepare.py` now writes `matrix/protein.qc.tsv.gz` for every matched candidate,
including excluded proteins. The default `--protein-filter target-bonf`
requires the measured protein's mapped coding gene to be present in its own
MAGMA file and to satisfy

```text
P_target_gene <= --protein-alpha / N_finite_MAGMA_genes
```

with `--protein-alpha 0.05`. The target HGNC symbol is converted to the MAGMA
gene ID using the `gene_loc` recorded in `magma.meta.tsv`. This is stricter and
more relevant to PHOLE than merely requiring some significant gene anywhere in
the genome.

Alternative sensitivity gates are:

- `any-bonf`: at least one Bonferroni-significant MAGMA gene anywhere;
- `cis-bonf`: at least one significant MAGMA gene overlapping the protein
  coding interval expanded by `--cis-flank`;
- `coding-bonf`: at least one significant MAGMA gene whose interval overlaps
  the unexpanded coding interval;
- `none`: retain any protein with at least one finite gene P value.

The QC table distinguishes `target_gene_absent_from_magma_output` from
`target_gene_not_bonferroni_significant`. Do not interpret an absent target gene
as biological evidence of no cis effect; the current shared-cache bug produces
exactly that artifact.

In addition to a random split, run a whole-chromosome split as a harder spatial
stress test. With the current within-chromosome distance kernel, a held-out
chromosome has no local training anchors, so its query kernel reduces to the
global background. Position permutations can consequently be degenerate
(`NULL_SD = 0`, empirical P = 1); interpret chromosome mode as a conservative
generalization check, not as a second independent position-permutation test.

## Putative black holes

A latent factor becomes a candidate only when it has all of these features:

1. It improves held-out trans prediction (`CV_DELTA_R2 > 0`).
2. Its source-side `alpha` loadings localize to a genomic interval.
3. Its target-side `beta` loadings affect a coherent, non-trivial set of proteins.
4. The source peak is not already explained by a measured protein coding anchor.
5. Its predictive gain beats target-position permutations after FDR correction.

The output calls this a **putative latent proteome-regulatory factor**, never a
discovered unknown protein. A peak could instead reflect LD, an assay-binding
artifact, protein clearance, a protein complex, glycosylation, ABO/MHC/GCKR-like
hotspots, or other shared biology.

## Requirements

- Bash
- Python 3
- `numpy` and `pandas` for prepare/fit/validate
- `matplotlib` for plot/all
- `scipy` for synthetic selftest only

If `--python` is omitted, `phole.sh` tests `PHOLE_PYTHON`, `python3`, and common
user conda locations and selects the first environment with the packages needed
for the requested step. On this machine it selects
`/home/huangj/anaconda3/bin/python`; the system `/usr/bin/python3` is missing
`pandas`, `scipy`, and `matplotlib`.

The matrix is stored as float32 `score.npy`; 18,000 genes by 3,000 proteins is
about 216 MB on disk. Model fitting uses additional memory for randomized SVD.

## Input files

The stable PHOLE input is `*.genes.out` accompanied by `magma.done`;
variable-width `*.genes.raw` is retained for MAGMA diagnostics but is not
parsed. Requiring `magma.done` is the default and prevents a run from consuming
a file still being written. MAGMA files must contain:

```text
GENE CHR START STOP NSNPS NPARAM N ZSTAT P
```

Generic protein BED, without or with a header:

```text
CHR START END PROTEIN [GENE]
```

`PROTEIN` must match the MAGMA filename prefix. `GENE` is optional because cis
matching uses coordinates. An optional `--protein-map` can supply `PROTEIN GENE`.
The BED build must match the MAGMA gene-location build; filenames containing
`b37` or `b38` are checked automatically.

The actual local files are:

```text
/mnt/d/files/ppp_3k.38.bed   CHR START END Assay Panel
/mnt/d/files/ppp_3k.38.tsv   PROT metadata, including Assay and HGNC.symbol
```

The local BED's fifth column is an Olink panel label, not a gene symbol.
`prepare.py` detects this low-cardinality annotation and defaults `GENE` to the
Assay/PROTEIN name. Passing the TSV through `--protein-map` uses its
`Assay -> HGNC.symbol` mapping explicitly.

## Run

Do not start the formal output directory until corrected upstream MAGMA
generation is complete. File completion alone is necessary but no longer
sufficient. Every expected PROT directory must have both `<protein>.genes.out`
and `magma.done`, and the target-gene/annotation checks above must pass. The
current project has 2,940 PROT protein directories, so first check:

```bash
find /mnt/d/data.BIG/gwas/prot/common -mindepth 1 -maxdepth 1 -type d | wc -l
find /mnt/d/data.BIG/gwas/prot/common -mindepth 3 -maxdepth 3 \
  -type f -name '*.genes.out' | wc -l
find /mnt/d/data.BIG/gwas/prot/common -mindepth 3 -maxdepth 3 \
  -type f -name 'magma.done' | wc -l
```

For the planned complete random cold-start run, use:

```bash
/mnt/d/scripts/phole/phole.sh all \
  --gwas-root /mnt/d/data.BIG/gwas/prot \
  --category common \
  --protein-bed /mnt/d/files/ppp_3k.38.bed \
  --protein-map /mnt/d/files/ppp_3k.38.tsv \
  --out /mnt/d/data.BIG/gwas/prot/phole \
  --protein-filter target-bonf \
  --protein-alpha 0.05 \
  --train-frac 0.80 \
  --split random \
  --rank auto \
  --perm 100
```

Harder chromosome holdout:

```bash
/mnt/d/scripts/phole/phole.sh all \
  --out /mnt/d/data.BIG/gwas/prot/phole.chr \
  --gwas-root /mnt/d/data.BIG/gwas/prot \
  --category common \
  --protein-bed /mnt/d/files/ppp_3k.38.bed \
  --protein-map /mnt/d/files/ppp_3k.38.tsv \
  --protein-filter target-bonf \
  --train-frac 0.80 --split chromosome
```

Use a separate `--out` for a different split or parameter set. With
`--replace FALSE` (the default), existing step outputs are intentionally reused
even if a new command supplies different model parameters. Use `--replace TRUE`
only when you deliberately want to rebuild that output directory from scratch.

Synthetic end-to-end test:

```bash
/mnt/d/scripts/phole/phole.sh selftest --out /tmp/phole.selftest --replace TRUE
```

## Verified local tests

The following checks were run successfully. Test outputs are smoke/regression
artifacts, not scientific results to interpret.

- `bash -n`, `shellcheck`, and Python byte-compilation passed.
- On 2026-08-31, synthetic `selftest` completed prepare, fit, validation, and
  all three plots with the default target-gene gate. It retained 59/72 proteins
  and the default 80% rule split them into 47 training and 12 held-out proteins.
- A real mechanical test used the alphabetically first ten completed MAGMA
  files. The default target-gene gate retained six: A1BG, AAMDC, AARSD1, ABO,
  ACAA1, and ACADM. ABCA2, ABHD14B, and ABRAXAS2 lacked their mapped target gene
  in the MAGMA output; ABL1 target-gene `P=0.04694` was not Bonferroni
  significant. The isolated 17,555 x 6 run completed fit, validation, and all
  plots with a 4/2 split and three smoke-test permutations. This proves the new
  code path runs; it does **not** validate the flawed current MAGMA data.
- Historical pre-gate mechanical regressions used 72 proteins in random and
  chromosome splits, and a three-protein A1BG/ABO/ALPI test verified union-grid
  construction. They are retained only as software regression evidence because
  they used the now-disqualified upstream MAGMA files.
- The local `ppp_3k.38.tsv` Assay-to-HGNC mapping and no-`--replace` resume/skip
  behavior were both exercised successfully.

## Code flow

```text
phole.sh
  -> f/prepare.py   completed MAGMA discovery, union grid, BED/map matching, matrix/QC
  -> f/fit.py       split, inner rank selection, cis + alpha*beta, gamma, metrics
  -> f/validate.py  target-position permutations, empirical P/FDR
  -> f/plot.py      prediction, black-hole ranking, codenome-loading figures

f/common.py         atomic writes, genome helpers, kernels, SVD, metrics, BH-FDR
f/simulate.py       MAGMA-shaped synthetic data used by selftest
```

## Main outputs

```text
matrix/score.npy                 complete source-gene x protein score matrix
matrix/protein.qc.tsv.gz         all candidate proteins, MAGMA/target-gene gate and exclusion reason
matrix/source.tsv.gz             source gene coordinates
matrix/target.tsv.gz             protein child and coding-gene coordinates
model/cis.alpha.tsv.gz           fitted/predicted cis alpha
model/alpha.source.tsv.gz        source-to-factor alpha
model/beta.target.tsv.gz         factor-to-protein beta
model/gamma.edges.tsv.gz         sparse unexplained direct edges
model/black_holes.tsv            exploratory candidate factors
model/metrics.tsv                held-out trans metrics
validation/validation.tsv        position-permutation tests
validation/black_holes.validated.tsv
plot/prediction.png
plot/black_holes.png
plot/codenome_factors.png
```

## Rare LoF variants: separate analysis, not a MAGMA switch

Do not obtain a "rare-LoF MAGMA" result by filtering the present common GWAS
table. The standardized files contain SNP, alleles, EAF, N, BETA, SE, and P but
no functional consequence/LoF annotation. For A1BG only 27,925 of about 1.21
million variants have MAF below 1%, and only 2,015 are below 0.1%. More
importantly, summary-statistic MAGMA uses the 503-sample 1000 Genomes EUR panel
for LD; this is not an adequate carrier/LD representation for ultra-rare LoF
variants.

Use individual-level WES/WGS with ancestry/relatedness-aware gene-based tests:
REGENIE burden/SKAT-type masks, SAIGE-GENE+, or an equivalent validated
collapsing framework. Predefine masks such as high-confidence pLoF and damaging
missense at MAF thresholds (for example <1%, <0.1%, and singleton/ultra-rare),
apply minimum MAC/carrier QC, and report burden and variance-component tests.
The UK Biobank PROT rare-variant study used exome data from 49,736 participants
and gene-level collapsing, identifying 1,962 gene-protein associations; this is
the relevant design class, not summary-P MAGMA
([Dhindsa et al., Nature 2023](https://www.nature.com/articles/s41586-023-06547-x)).
Rare coding architecture is also much sparser than common-variant architecture,
with much of the burden contribution concentrated in ultra-rare pLoF variants
([Weiner et al., Nature 2023](https://pmc.ncbi.nlm.nih.gov/articles/PMC10614218/)).

For PHOLE, retain corrected common-variant MAGMA as the dense primary matrix.
Use rare-LoF burden associations as an orthogonal, sparse validation/anchor
layer for selected source-gene/protein edges rather than replacing the entire
dense matrix.

MAGMA summary-P analysis remains appropriate only after its annotation and P
inputs are fixed: it aggregates SNP associations while accounting for LD, and
the reference must overlap the SNP input and match ancestry
([MAGMA paper](https://pmc.ncbi.nlm.nih.gov/articles/PMC4401657/),
[MAGMA manual](https://ibg.colorado.edu/cdrom2021/Day10-posthuma/magma_session/manual_v1.09a.pdf)).

## Recommended analysis sequence

1. Fix the shared annotation and `P=0` handling, rebuild MAGMA, and pass the
   target-gene preflight audit.
2. Run PHOLE on the target-gene-filtered common-variant MAGMA matrix with an
   80/20 protein cold-start split.
3. Repeat after excluding MHC, ABO, GCKR, immunoglobulin regions, and known assay
   artifacts one block at a time.
4. Compare random and chromosome cold-start splits, plus `target-bonf` versus
   `any-bonf` sensitivity analyses.
5. Add rare-LoF burden results as orthogonal validation, not as a replacement
   for the dense common-variant matrix.
6. Rebuild a signed locus-by-protein matrix from harmonized full GWAS effects.
7. Add colocalization and fine-mapping before interpreting a source locus.
8. Replicate in an independent proteomic platform/cohort.
9. Only then nominate perturbation or biochemical experiments.
