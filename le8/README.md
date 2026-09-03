# LE8 5C analysis pipeline

This directory contains the LE8 supervised 5C omics pipeline.

## Run

```bash
cd /mnt/d/scripts/le8
./le8.sh -Y cvd_cad --biom prot,met
```

The LE8 wrapper passes bounded matching defaults to every sourced
`match_GRCH`/`match_SNP` call:

```bash
./le8.sh -Y cvd_cad --biom prot \
  --match-memory-mb 4096 --match-sort-memory-mb 512 \
  --match-tmp-dir /mnt/d/tmp
```

The same controls can be used directly. The 4 GiB `ulimit` is installed only
inside a function subshell, while the query/reference keys are external-sorted
on disk with a 512 MiB buffer:

```bash
source /mnt/d/scripts/0f/0phe.f.sh
match_GRCH --memory-mb 4096 --sort-memory-mb 512 --tmp-dir /mnt/d/tmp \
  --reference reference.tsv --output matched.tsv --audit matched.audit.tsv query.tsv
```

Create `/mnt/d/tmp` on a disk with enough free space. `--max-query-rows`
(default 5,000,000) catches a likely QUERY/REFERENCE reversal; use
`--allow-large-query` only after checking the argument order. For a strict
aggregate cap across AWK, sort, R and child processes, use a cgroup rather than
a shell-only per-process limit:

```bash
systemd-run --user --scope -p MemoryMax=16G -p MemorySwapMax=4G \
  bash /mnt/d/scripts/le8/le8.sh -Y cvd_cad --biom prot
```

Use `--replace TRUE` after upgrading this code or when every selected stage
must be recomputed. The large C1/C2/C5 caches now carry a code-version tag,
but not every auxiliary cache has a complete content hash of all upstream
individual and summary-statistic inputs.

## Incident/prevalent time contract

`0f/0phe.f.R::t2e()` defines incident events with `<Y>.Yt2e` and forward
follow-up in `<Y>.t2e`. Diagnoses before baseline are represented by
`<Y>.r2e`; `<Y>.b2e` is signed time around baseline (negative for prevalent
diagnosis-to-baseline time, positive for baseline-to-incident-diagnosis time).
`<Y>.bi2e` is attained age at diagnosis or censoring (birth to event/censor).
The valid birth-origin sensitivity model is
`Surv(age_at_baseline, <Y>.bi2e, <Y>.Yt2e)`: participants enter the risk set at
their baseline blood draw. It does not turn UK Biobank into a birth cohort and
does not project an adult protein/metabolite measurement back to birth. Merely
adjusting for age, sex and baseline diagnosis cannot identify that
counterfactual early-life omic level.

C1 uses incident cases for ProtWAS/MWAS and models baseline-prevalent cases in
separate logistic and case-only years-since-diagnosis analyses. The former
reverse-time Cox is retained as an explicitly exploratory figure and table: it compared
time-since-diagnosis in prevalent cases with age-since-birth in controls and
therefore does not define a coherent survival risk set and is excluded from
directionality/evidence grading. C1 instead performs
0.5-, 1-, 2-, 5- and 10-year incident landmark scans, plus formal
0/0.5/1/2/5/10/16-year diagnosis-window logistic models. Incident controls
must remain observed and disease-free through the upper edge of a window;
later cases can therefore be controls for an earlier window. C5
fits prediction models only in the incident-risk cohort, but applies the fitted
scores to prevalent participants for Yin–Yang score displays and temporal
diagnostics.

Protein groups are read from columns 4–5 of `ppp_3k.b38.bed`; metabolite
groups are read from the named `group` column of `met.lst`, or its penultimate
column when no header is available. `c1.input_feature_annotation_audit.csv`
records every assayed column and whether it matched that source. This must be
checked when the wide metabolite RDS reports 300 features but the intended
analysis panel is described as 249 metabolites.

The Yin–Yang figures are diagnosis-anchored, cross-sectional pseudo-trajectories
from one baseline omics measurement per participant. They are not longitudinal
within-person biomarker trajectories. In manuscripts, use
"incident/prevalent diagnosis-anchored triangulation" as the primary term and
define Yin–Yang only as a visual shorthand.

The third row of `c1.Fig1` and `c1.Fig2` is the attained-age delayed-entry
sensitivity analysis. `c1.Fig10.landmark_birthline_sensitivity.png` displays
0.5/1/2/5/10-year landmark persistence and baseline-time versus attained-age
effect concordance. `c1.Fig11.diagnosis_window_riskset.png` displays adjusted
window-specific betas and confidence intervals across -16 to +16 years.
`c1.Fig13.reverse_time_exploratory.png` preserves the legacy reverse-time
result with its incompatible-time-origin warning. If GO
annotation and g:Profiler are both unavailable, C1 Fig8 falls back to a clearly
labelled protein assay-group enrichment rather than producing an empty canvas.

## C2 figure contract

- `c2.Fig1.prots.top.png`: cis/local is always the left circle,
  trans/distal is always the top circle, and observational evidence is always
  the right circle. The panel-specific reference set is pale yellow.
- `c2.Fig2.pQTL_R2.png`: shows the distribution of per-instrument partial R²
  separately for cis/local and trans/distal instruments. The point is the
  median, the colored interval is the interquartile range, and the grey tail
  ends at the 90th percentile. It never sums R² and never applies a many-IV
  product approximation.
- `c2.Fig6.dandelion.png` contains the target and pair-level summaries. The
  dense path network is `c2.Fig9.dandelion_network.png`.
- `c2.Fig8.dandelion_mr_integration.png` displays at most 32 targets, ranked
  first by independent MR/observational support; its second panel summarizes
  support tiers across every DANDELION target.
- `c2.Fig12.genetic_decomposition.png` and
  `c2.Fig13.genetic_component_leadtime.png` separate the observed omic,
  its calibrated COJO-PGS component and the remaining non-genetic residual.
  `c2.Fig14.evidence_grades.png` assigns cis/local MR grades after explicitly
  auditing many-IV architecture, Cochran-Q heterogeneity, weighted-median
  concordance, Egger and Steiger directionality.


## Useful controls

```bash
./le8.sh --preflight
./le8.sh -Y cvd_cad --biom prot --steps c1_correlate,c2_cause
./le8.sh -Y cvd_cad --biom met --from c1_correlate --to c2_cause
./le8.sh -Y cvd_cad --biom prot,met --replace TRUE
```

The default locations and their CLI overrides are listed by
`./le8.sh --help`.

## GWAS directory contract

`le8.sh` consumes the trait-first layout written by `gwas_post.sh`:

```text
<project>/common/<trait>/gwas/<trait>.gz
<project>/common/<trait>/gwas/<trait>.cis.gz
<project>/common/<trait>/gwas/<trait>.jma.cojo
<project>/common/<trait>/gwas/<trait>.ldr.cojo
<main-project>/common/<trait>/magma/<trait>.genes.out
```

`LE8_GWAS_DIR`, `LE8_PQTL_IV_DIR`, and `LE8_MQTL_IV_DIR` are project roots,
not `common` or `clean` directories. `LE8_PROT_BED` defaults relative to the
pQTL project root. For a single-trait diagnostic run only, `--outcome-gwas` can select
an outcome file directly; the wrapper passes that exact resolved path to every
R job so preflight and runtime cannot diverge.

LE8 normalizes every GCTA `.jma.cojo` SNP column against the corresponding
full GWAS with `match_GRCH` from `0f/0phe.f.sh`. This permits a PLINK reference
to use `chr:pos:ref:alt` while the full outcome/QTL files use rsIDs; ambiguous
multi-allelic matches are excluded rather than joined by position alone. An
unchanged rsID can also bridge GRCh37/38 positions when its alleles agree; LE8
then uses the full GWAS coordinates for downstream cis/local classification.

Outcome GWAS genome builds are detected independently for each trait. The
default `--grch auto` accepts both GRCh37 and GRCh38 inputs in the same run
and passes the detected build to that trait's jobs. Use `--grch 37` or
`--grch 38` only when every selected trait must match a required build.

MR-link-2 is an expensive sensitivity analysis and uses the exact
case-sensitive control `RUN_MRlink2=Top|All|None` (CLI:
`--run-mrlink2 Top|All|None`). The default `Top` runs the prespecified
PCSK9/LPA/GDF15/NTPROBNP/MMP12 anchors when available, then the strongest
cis/local MR and C1 candidates, capped by `C2_TOP_MAX` (default 50). `All`
runs every eligible region and `None` skips the method. When selected, it
uses the same detected build to select
`<LE8_REFGEN_ROOT>/<37|38>/pfile` (default root:
`/mnt/d/data.BIG/refGen/1kg`). Both `.pvar` and compressed `.pvar.zst`
pfiles are accepted. The default population is EUR. If
`<1kg>/<build>/bfile/EUR/chr*` or `<pfile>/EUR.chr*` already exists it is used
directly; otherwise the runner loads `<pfile>/chr*` with PLINK2's `vzs`
modifier when needed and applies
`--keep` from the sibling `id/EUR.id.2col`. If that persistent ID file is
absent, the runner falls back to rows whose final `super_pop` column is `EUR`
in `samples.txt`. The filtered BED cache is written under each analysis
output, with multiallelic sites removed via `--max-alleles 2`, so separate
`EUR.chr*` pfiles do not need to be generated.

Regional inputs are mapped to the selected PLINK IDs by `match_GRCH`.
`MRLINK2_REF_PFILE_DIR`, `MRLINK2_REF_BED`, `MRLINK2_REF_POP`,
`MRLINK2_REF_ID_DIR`, or `MRLINK2_REF_SAMPLES` may override the defaults;
preflight rejects a build mismatch or an invalid population keep source.

## C2 MR and DANDELION input contract

MR takes the independent SNP set from each exposure's `.jma.cojo`, then reads
alleles, effects and coordinates from the matching full `.gz`.  `.cis.gz` is a
regional extract, not a replacement for a genome-wide independent instrument
set; `.ldr.cojo` is GCTA LD output and is retained as provenance rather than
treated as summary statistics.

For the protein-layer DANDELION analysis, disease exposures come from the
outcome `.jma.cojo`, each candidate protein's complete `.gz` supplies the
trans-pQTL P values at those SNPs, and a gene-level disease P-value file supplies
the second component.  A WES burden/gene-level file is preferred:

```bash
./le8.sh -Y cvd_cad --biom prot --steps c1_correlate,c2_cause \
  --dandelion-gene-p /mnt/d/data.BIG/gwas/main/common/cvd_cad/wes/cvd_cad.gene_p.tsv \
  --dandelion-gene-annotation /mnt/d/data.BIG/refGen/genes/GRCh37.genes.tsv
```

Without `--dandelion-gene-p`, C2 uses
`main/common/<trait>/magma/<trait>.genes.out` and labels the result explicitly as
an adapted GWAS gene-level sensitivity analysis, not WES burden evidence. It
is retained for backward compatibility but is excluded from primary C5
evidence grading. Disable this fallback with
`--dandelion-allow-magma FALSE`. C2 also flags a run when selected DPGs exceed
25% of tested targets (change with `--dandelion-max-target-fraction`); a broad
selection is not counted as primary evidence. Outputs now
include input/QTL coverage audits, all tested pairs, selected paths, DPG hubs,
MR/observational integration, a complete evidence landscape, a locus-to-DPG
heatmap, and the official `DANDELION::gen_fig` network PDF.

DANDELION uses the exact `RUN_Dandelion=Top|All|None` contract (CLI:
`--run-dandelion`). `Top` restricts candidate disease-proximal proteins to the
same audited C2 top set; `All` tests every eligible assayed protein; `None`
skips it. It is not run for metabolites: the method prioritizes gene/protein
nodes through trans-regulatory and gene-level disease evidence, and there is
no coherent metabolite analogue merely because mQTLs exist. Metabolite C2
Figures 6–9 are therefore retained as explicitly unavailable panels.

### Bidirectional MR and individual genetic decomposition

C2 also runs disease-liability-to-omic MR using independent outcome-GWAS lead
SNPs against each full pQTL/mQTL summary. This reverse direction estimates
genetic liability, not the effects of established disease, treatment or
survival. `c2.Fig11.bidirectional_mr.png`, `c2.reverse_MR_all.csv` and the
corresponding workbook sheets keep the two directions separate.

`c2.genetic_score_weights.tsv` is a PLINK-style manifest of cis/local weights.
If individual genetically predicted omic scores have been generated from UKB
genotypes, C2 auto-discovers `Rdata/prot.pgs.rds` or `Rdata/met.pgs.rds`.
Expected columns are `eid` plus `FEATURE.pgs`; `_pgs`, `.PGS`, `_PGS`,
`_GRS`, `GRS_FEATURE`, and an exact feature-name fallback are also recognized.
`C2_GENETIC_SCORE_FILE` may override the path. C2 calibrates each score among
participants known to be event-free at year ten (including later cases) and writes an
observed/genetic/residual comparison. The residual is explicitly a
non-genetic residual—it still contains lifestyle, environment, assay error,
treatment and latent disease and must not be called a pure preclinical-disease
component. External or cross-fitted QTL weights are preferred because weights
estimated in the same UKB proteomics sample can overfit and suffer winner's curse.

## C5 censoring contract

The glmnet panels fit penalized Cox models rather than any-event logistic
models. LightGBM remains a fixed-horizon classifier because upstream LightGBM
does not expose a native Cox objective; its training set excludes participants
censored before `C5_HORIZON`. The headline AUC is an IPCW cumulative/dynamic
AUC at that horizon, while Harrell's C-index uses the full follow-up.

`c5.Fig4.leadtime_prediction.png` reports score and individual-biomarker AUCs
at increasing minimum lead times, MASLD-style within-5/within-10/over-10/
over-12-year windows, and the eligible case/control counts. A red marker is
placed only at the farthest *consecutively* supported horizon: every tested
lead time through H must have lower 95% AUC CI > 0.50 and at least 100 cases
(change with `--c5-lead-min-cases`). Maximum follow-up alone never earns an
"up to N years" claim. These curves use a minimum-lead-time case/control AUC;
they are intentionally distinguished from the headline 10-year IPCW
cumulative/dynamic AUC.

`c5.Fig5.score_vs_topN_mechanism.png` compares a training-locked top-1
biomarker, an unweighted top-N mean, and a training marginal-Cox-weighted top-N score on
the same held-out participants. It shows AUC, validation Cox beta and CI,
case-control mean separation, pooled within-group SD, Cohen's d, bootstrap
stability, risk-set-matched diagnosis-anchored score separation, average biomarker correlation
and effective N. The grey band is the distribution of the N individual
biomarkers; the best validation biomarker is displayed only as an optimistic
reference and is never used for model selection. Raw SD is not itself the
mechanism: aggregation helps when mean separation grows relative to
within-group noise.

`c5.Fig12.attained_age_sensitivity.png` compares baseline-time and
attained-age delayed-entry Cox effects for each score. The large-metabolite
branch writes Figs 1–3 before optional matching/bootstrap diagnostics and caps
diagnostic resampling sizes, so an interrupted diagnostic cannot leave C5
without its primary figures.

## C3 credible-set diagnostics

`c3.Fig4.credible_sets.png` uses a log posterior axis so a one-SNP posterior
near 1 does not make all alternative variants visually disappear. It also
shows credible-set size against robust PP(H4), the full credible-set ECDF and
an explicit count of one-SNP posterior concentrations. Such concentration is
not automatically accepted as precise fine-mapping; verify the LD reference,
allele harmonization and single-causal-variant assumption using
`c3.credible_set_audit.csv` and `c3.credible_set_by_locus.csv`.

## Directionality and supervised modules

C1 contrasts incident, five-year landmark, baseline-prevalent logistic and
case-only years-since-diagnosis evidence. The resulting
`c1.Fig9.directionality_triage.png` and C2 integration are anchored by PCSK9,
LPA, GDF15, NTPROBNP and MMP12 by default; override them with
`--direction-anchors GDF15,PCSK9,...`.
C1 labels a signal as reactive-compatible rather than as a proven consequence:
prevalent associations can also reflect treatment, survival selection,
prevalence-incidence bias and preclinical disease. A valid forward MR plus
robust colocalization is therefore allowed to coexist with reactive evidence.
C4 retains its detailed figures and adds a consolidated supervised module atlas,
network globe, bootstrap module membership, and stable `YS_core` selection.
The atlas, network globe and selection/mediation composites are numbered
sequentially as `c4.Fig7`–`c4.Fig9`; `c4.Fig1`–`c4.Fig6` are the preceding
ordered figures.
Control module resampling with `--c4-module-boot` and
`--c4-module-stability`.

C5 reports both the original five-stage support count and a safeguard-aware
grade. Grade A requires observational association, forward MR and robust
colocalization. Grade B requires causal/locus evidence plus observational or
distal landmark support. Reactive evidence changes the assigned biomarker role
(for example, "causal + reactive mixed") but does not veto otherwise coherent
genetic evidence. C5 stability is shown as unavailable, rather than zero, when
nested cross-validation was not requested.

C5 has one upstream hard prerequisite: a readable C1 result for the same omic
layer. C2, C3 and C4 are optional. Missing, unreadable or explicitly
unavailable optional modules are recorded in `c5.upstream_availability.csv`;
their fixed figure columns/panels remain in place and are labelled unavailable
instead of terminating C5 or silently counting missing evidence as zero.

C5 additionally compares three prespecified score roles: a five-year distal
antecedent panel, a genetic-causal panel (the same feature must have cis/local
MR at FDR < 0.05 and robust colocalization), and their hybrid. Trans/distal-only
MR remains available as a separate sensitivity score. Cox linear predictors are retained as linear predictors and standardized
on the training set; they are not passed through `plogis()` or presented as
calibrated probabilities. The ROC panels use the same IPCW cumulative/dynamic
definition as the tabular headline AUC, and baseline prevalence is presented as
a cross-sectional proportion rather than a pseudo-survival curve.
