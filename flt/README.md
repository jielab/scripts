# Food-liking traits and brain health

This release is a replacement for `scripts/flt/`. It was prepared from
[jielab/scripts at 73be1f3](https://github.com/jielab/scripts/tree/73be1f3734b269028e991b9cb8dd871989113efa/flt)
and audited against [the current results and old draft at d9cccb4](https://github.com/jielab/analysis/tree/d9cccb43448e43228556a83cb7f01494a7351b9c).
Keep the repository's shared `0f/` and `diet/` directories alongside `flt/`.

## Run on HPC or locally

```bash
cd scripts/flt
./install.sh run
./flt.sh run --check-inputs --panel extended
./flt.sh run --panel extended --wait
./flt.paper.sh run --panel extended --wait
```

The extended panel uses all 28 configured outcome-side GWAS releases, including
education. Use the default `--panel primary` for a faster first pass with the
five draft diseases and education. Use the same panel on both entries. Both
entries show help when called with no arguments. `--dry-run` prints the job
graph without writing output. `--plan-only` writes the plan and command
files. `--dir-output DIR` overrides the entire analysis root.

The shared `0f/0os.sh` detects HPC and selects paths. The default maximum is
20 CPUs per job on HPC and 16 locally. Serial stages request one core.
The reverse MR, LDSC and PLACO scans use at most eight processes per job.
On HPC, independent stages submit together and dependent jobs wait on LSF job
IDs. On Linux/WSL without LSF, top-level jobs run sequentially with bounded
parallelism inside selected stages. LSF applies the queue's memory and walltime
defaults. If a heavy job needs over half a queue-33 node, `--whole-node`
reserves all 40 slots while limiting scientific threads to the selected budget.
Without `--wait`, HPC exits after submission and prints a receipt line per job.

To rerun a selected stage, use `./flt.sh run MODULE`. For example,
`./flt.sh run ldsc --panel extended` regenerates the LDSC products, while
`./flt.paper.sh run fig2 --panel extended` redraws Fig. 2 from saved estimates.
On a fresh analysis root, start with the full graph to build prerequisites.
The old layouts can stay alongside the new `main/` until you finish comparison.
The paper entry reads the new layout.

| Location | Contents |
|---|---|
| `main/results/01cox`, `01intake`, `01lm`, `01protein` | Aggregate cohort, dietary, behavioral, MRI and protein analyses |
| `main/results/02uvmr`, `02mvmr`, `02ldsc` | Genetic estimates, counts and diagnostics under the named panel |
| `main/results/03placo`, `03phewas`, `04smr` | Shared SNPs, candidate windows, functional lookups and SMR |
| `main/reports`, `qc` | Cohort flow, resource notes and per-stage status |
| `logs`, `jobs`, `tmp` | Logs, readable LSF plans and resumable intermediates, including individual-level cohort files under `tmp/cohort` |
| `paper/main/results` | PNG figures, the source workbook and separate plot data |
| `paper/main/reports`, `paper/logs`, `paper/jobs`, `paper/qc`, `paper/tmp` | Coverage, draft comparison and paper runtime |

## Study design

| Question | Analysis and reporting |
|---|---|
| Do liking scores predict subsequent disease? | The primary time origin is YGB field 20750, the completed food preferences questionnaire. Exclude a diagnosis on or before that day. End follow-up at the first of death, loss or the shared administrative cutoff. Fit per-point and grouped Cox models. Use age at the FPQ and the inherited BMI, socioeconomic and lifestyle adjustment. BH across the 182 planned traits within outcome and contrast. |
| Does the observation survive alternative choices? | Complete-case analysis, two-year landmark, and recruitment-time-origin history for comparison with the draft. The focal G732 analysis compares common samples with and without the mean liking of other raw items. It also controls for WebQ poultry intake when the matching fields are available. The comparator excludes chicken. The poultry proxy is the sum of fields 26069 and 26121, transformed by `log1p`. These comparisons address score specificity and recorded intake, with the limits of self-report and observational adjustment. |
| How much genetic sharing is present? | LDSC for all planned FLT–GWAS pairs. Report observed-scale h², its SE, genetic covariance, free intercept and genetic correlation with delta-method CI. The full panel-wide BH and within-GWAS BH are both exported. The main figure flags h² Z below 4 as low precision. Focal HDL is a separate sensitivity analysis if the reference is available. |
| Is there evidence for a direction? | UVMR independently selects and clumps each direction's exposure instruments at P < 5e-8, 10 Mb, r² < 0.001, and applies minimum per-SNP F = 10. Primary BH uses valid IVW tests within GWAS and direction. The planned-family correction remains available in a separate column. Wald ratio and the four sensitivity methods have no IVW FDR. Export Q, Egger intercept, Steiger eligibility and focal leave-one-out. |
| Does the focal signal persist with education and related exposures? | MVMR covers forward IVW discoveries for the five primary disease releases plus G732–MDD regardless of its IVW P. Include `covs/common/income`, `brain/common/edu` and `covs/common/physical_activity`. Report conditional F and Q, plus the assumption about SNP–exposure estimation error covariance. Education is also a continuous comparator in LDSC and both MR directions. |
| Which genomic regions could underlie sharing? | Focal PLACO retains the historical composite-null implementation and exports all SNPs, candidate SNPs at P < 1e-5 and FUMA upload data. Add predeclared pairs to `config/placo-pairs.tsv` after verifying both files' genome build. Coloc runs a single-signal Wakefield ABF in greedy, non-overlapping native windows around the strongest candidate SNPs. It uses marginal effects from both original files, including SNPs removed by PLACO's Z² cutoff. Export H0–H4 at three joint priors, SNP coverage, and the top conditional shared variant. |
| What does functional follow-up support? | PheWAS counts distinct FUMA loci per canonical phenotype. Show association records and source coverage separately. FUMA, MAGMA, deTS and Metascape provide mapping and enrichment evidence. SMR preserves the old nominal shared-probe screen and adds Bonferroni-adjusted results within tissue and GWAS. The pharmacological interpretation stays at the candidate-mechanism level. |

Prospective association, shared genetics and MR are distinct questions. The MR
workflow does not use Cox significance, LDSC significance or PLACO results to
select instruments. See the [shared study-design discussion](https://chatgpt.com/share/6b06a971-bf0c-83e8-a789-24a288b62c8e)
and [Han et al., Communications Biology, 2022](https://www.nature.com/articles/s42003-022-04070-9).
A reverse estimate with nominal P ≥ 0.05 is reported as an imprecise or null
test result, without treating it as evidence against reverse influence.
Binary-exposure reverse estimates are on the GWAS liability-proxy scale.
Food-liking MR estimates use the source GWAS beta units. They are distinct from
a change of one survey point and from an intervention on food consumption.

For the main MDD analysis, confirm that the local `mdd_2025` file is the PGC
European-ancestry release excluding YGB and 23andMe. Verify release, ancestry,
scale and cohort overlap before interpreting the extended panel, especially
alternate and non-`_eur` stroke releases. The supplied EUR LD reference is used
in these exploratory fits. Alternate releases are not automatically treated as
independent replication. For Steiger, optionally supply a `FLT_TRAIT_INFO`
CSV with `key,type,N,ncase,ncontrol,prevalence` using qualified IDs such as
`brain/mdd_2025`. Declared metadata and missing-information flags are exported. For MVMR, an optional
`FLT_MVMR_ERROR_COR_DIR` can supply one named error-correlation matrix per
food in `ID.csv`. Row and column names correspond to the food GWAS ID,
`income`, `edu` and `physical_activity`.

## Input locations

| Input | Default |
|---|---|
| FLT dictionary | `$DIR_FILES/diet/flt.xlsx`, including `item182` and `Sheet1` |
| YGB extracts and RDS | `$DIR_YGB/rap`, `$DIR_YGB/common`, `$DIR_YGB/Rdata` |
| Harmonized GWAS summary statistics | `$DIR_GWAS/{flt,brain,covs}/common/ID/gwas/ID.gz` with SNP, CHR, POS, EA, NEA, EAF, N, BETA, SE, P and LOG10P |
| MR LD | `$DIR_GEN/1kg/37/openGWAS/EUR.{bed,bim,fam}` via `FLT_MR_LD_PREFIX` |
| LDSC and HDL | `$DIR_LD/ldsc/hm3` and `$DIR_LD/hdl/YGB_imputed_SVD_eigen99_extraction` via `FLT_LDSC_REF` and `FLT_HDL_REF` |
| SMR | `FLT_SMR` executable and `FLT_GTEX`, `FLT_EQTLGEN` complete BESD prefixes |
| Imaging | `$DIR_YGB/rap/img.tab.gz` and `$DIR_FILES/diet/MRI_analysis_dictionary.csv`. `FLT_MRI_ATLAS` and `FLT_MRI_TOOLS` select anatomical-rendering resources. |
| Optional external returns | `FLT_FUMA_DIR`, `FLT_ATLAS_DIR` and `FLT_METASCAPE_RESULTS` default under `$DIR_FILES/flt/03placo` |

Use `FLT_YGB`, `FLT_GWAS`, `FLT_DICTIONARY`, `FLT_RSCRIPT` and
`FLT_PYTHON` to override input and interpreter paths. `./install.sh --check`
checks the starting environment. The installer uses the shared R and Python
environments and places extra R packages outside the analysis root. MRI
anatomical rendering uses the existing diet graphics runtime. It has a
coefficient-heatmap fallback when the renderer or atlas is unavailable.

### FUMA, PheWAS and Metascape

The focal PLACO stage generates
`main/results/03placo/G732-mdd_2025_PLACO_FUMA.txt.gz`. For a new FUMA job,
use the rsID-based upload and the intended EUR reference. Retain the original
return and parameters under `FLT_FUMA_DIR`. FUMA lead SNP IDs, chromosomes
and P values are checked against the newly generated PLACO candidate table before
archived annotations are used. The native PLACO plot uses the local GWAS
coordinates, while FUMA-based PheWAS uses rsID and the FUMA coordinates.

Put available Atlas CSV files in `FLT_ATLAS_DIR`. Read-only externally acquired
OpenGWAS and GWAS Catalog records are read from `FLT_PHEWAS_INPUT`.
New query responses and derived audit files are written under
`analysis/flt/cache/phewas` (`FLT_PHEWAS_CACHE`);
`--dir-out-phewas-cache` overrides that writable cache location. Set `FLT_PHEWAS_ONLINE=true` only when the worker can
reach the services and, for OpenGWAS, has an authorized JWT. Then run
`./flt.sh run phewas`. Missing database-source coverage is recorded.

`./flt.paper.sh run fig5` makes the PLACO Manhattan plot before checking the
external returns. Matching FUMA requests also prepare the Metascape candidate
lists in `paper/main/results/03fuma/metascape`. Supply the corresponding native
`Enrichment_GO/_FINAL_GO.csv` and rerun Fig. 5. Partial external coverage
is shown in `paper/main/reports/fig5-coverage.txt`.

## Manuscript outputs

| Output | Main panels |
|---|---|
| Fig. 1 | Cohort characteristics and the prospective Cox landscape |
| Fig. 2 | Forward IVW scan and comparison with MVMR |
| Fig. 3 | Cross-method evidence, actual intake, dietary profiles and cognition/mental health |
| Fig. 4 | MRI, limma liking–intake comparison and protein–disease follow-up |
| Fig. 5 | PLACO and PheWAS, mapping/SMR evidence, MAGMA/Metascape pathways and tissue results |
| Additional genetics figures | LDSC overview, focal chicken-liking Cox sensitivities and colocalization posterior probabilities |
| Source workbook | Aggregate data, methods/diagnostics, prior-draft comparisons and an explicit coverage sheet |

The current archive had Cox, intake, behavioral, protein, MRI, MR, focal PLACO,
PheWAS and SMR outputs. A combined paper directory, a completed LDSC-wide
scan, full focal PLACO summary statistics for plotting, and the native external
annotation returns were absent from that archive. The old draft's counts are
reported alongside newly computed values in `Draft-comparison.csv`. A new
time origin, harmonization, instrument availability or an updated reference can
change those numbers. The three named food pairs are kept for draft continuity
as exploratory follow-up, regardless of the new screen.

## Verification

The release was checked with Bash syntax and R parsing, the SMR regression
tests, archived aggregate results for publication rendering, synthetic cohort
and GWAS inputs for the new Cox sensitivity, PLACO and colocalization stages,
and mocked LSF submission. The LDSC adapter was exercised on a synthetic pair
using the upstream GenomicSEM regression, and Fig. 5 was rendered with a synthetic
FUMA return. Running all analyses on the actual YGB data,
GWAS and HPC reference files is required to obtain the revised numerical
results. No unpublished participant data is bundled.

## Generated results and external inputs

PLACO results are generated under `analysis/flt/main/results/03placo/G732-mdd_2025_PLACO_*`; all downstream consumers read these results. Computation checkpoints are stored under `analysis/flt/cache/placo/`, and PheWAS query and derived-audit caches under `analysis/flt/cache/phewas/`.

Only externally obtained inputs remain under `$DIR_FILES/flt/03placo`: FUMA_job787285, GWAS_Atlas, metascape, and read-only OpenGWAS/GWASCatalog query records in phewas_cache. These source records are not written by this workflow. Classification and summaries are regenerated. New API responses are written to the analysis cache. eQTLGen remains optional and unavailable eQTLGen analyses are skipped.
