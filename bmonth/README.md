# bmonth: UK Biobank birth-month and birth-location analysis

This pipeline mirrors `scripts/border` but treats month as a **circular** exposure rather than an ordinary rank.
Input variables already present in `all.rds`:

- `p52 -> birth_month`
- `p129 -> birth_north`
- `p130 -> birth_east`

## Main exposure definitions

| exposure | definition | interpretation |
|---|---|---|
| `bmonth_cos` | `-cos(2π month/12)` | winter (-1, Dec) to summer (+1, Jun) axis |
| `bmonth_sin` | `sin(2π month/12)` | orthogonal spring/autumn timing axis |
| `bmonth_winter_summer` | Dec-Feb=0, Jun-Aug=1; others NA | summer vs winter births |
| `bmonth_spring_autumn` | Mar-May=0, Sep-Nov=1; others NA | autumn vs spring births |
| `birth_north_100km` | birth north coordinate /100,000 | per 100-km northing |
| `birth_east_100km` | birth east coordinate /100,000 | per 100-km easting |
| `bmonth_rank` | Jan=1,...,Dec=12 | **secondary sensitivity only** |

`bmonth_relative` is deliberately not created: `(month-1)/11` is exactly a linear rescaling of month rank and gives no information beyond `bmonth_rank`. More importantly, Jan and Dec are adjacent biologically but far apart on a 1-12 linear scale, so rank/relative should not be the primary model.

## Joint PheWAS tests

In addition to scalar PheWAS figures, `3.phewas.R` runs four likelihood-ratio PheWAS tests:

1. `month11df`: unrestricted 12-month categorical test (11 df).
2. `cycle2df`: joint sine + cosine annual-cycle test (2 df).
3. `north_interaction2df`: whether the annual cycle differs by birth northing, adjusted for northing and the two seasonal main effects.
4. `east_interaction2df`: analogous interaction with birth easting.

The north interaction is particularly relevant to the hypothesis that the same calendar month can represent different prenatal temperature, photoperiod, infection, and food environments at different locations.

## Run

```bash
cd /mnt/d/scripts/bmonth
./bmonth.sh check
./bmonth.sh all
```

`./bmonth.sh` with no arguments also runs all stages and saves the full console
log under `/mnt/d/analysis/bmonth/logs/`. Unknown PheWAS dictionary codes remain
in the result tables and plots under the `other` category.

Every analysis invocation runs in its own systemd scope with a **24 GiB RAM
limit for the entire process tree** and **no swap**. R and its workers share
that limit; an OOM inside this scope terminates the analysis together. The
launcher verifies the kernel limits and refuses to run if the guard is
unavailable. Help does not require systemd.

PheWAS defaults to **automatic parallelism, up to eight workers**, reduced when
the available RAM or the 24 GiB process-tree budget is insufficient. Each worker
shares the current model's read-only design matrix using Linux fork and reads
one phenotype at a time from a packed disk cache. The full phenotype table is
released before forking. Workers communicate through local pipes, without TCP
sockets. BLAS/OpenMP threads are limited to one per worker.

If `raw/bmonth.cohort.rds` already exists, run only PheWAS with:

```bash
cd /mnt/d/scripts/bmonth
./bmonth.sh phewas
```

Run `./bmonth.sh prep` first only when that prepared cohort is missing.
Use `--phewas-cores 8` to request eight workers, or `--phewas-cores 1` for a
serial run. Explicit requests are also bounded by available memory and CPUs.

Sample alignment and covariate encoding are performed once. Regression uses
the same `stats::glm.fit` IRLS/QR solver as `glm`, with the same default
convergence settings. Scalar tests retain Wald P values; joint tests retain
likelihood-ratio P values and use the fitted model ranks for degrees of freedom.
The month and annual-cycle tests share a null fit when their sample rows and
covariate matrices are exactly identical. Disease-specific missing values,
seasonal subcohorts and minimum-case filters are preserved.

Every eight phenotypes per worker are checkpointed atomically under
`raw/phewas_cache_v2/`. Logs report completed phenotype counts, percentage and
an estimated remaining time **for the current group**. Repeating the command
resumes completed blocks, including after an interrupted group. Input-content
hashes, model code, covariates, contrast settings and R/PheWAS versions identify
the cache; older unversioned exposure caches are not used for computation.
`--replace` recomputes the requested groups. The packed phenotype cache uses
one byte per observation (about 0.7 GB for the current cohort).

Each completed scalar exposure is saved and plotted immediately; the two
calendar tests are published together, followed by each geographic interaction.
`phewas.rds/csv` and `phewas_joint.rds/csv` contain the groups completed so far
during a run. All four joint tests still run when `--exposures` selects only a
subset of scalar exposures. Complete-run output filenames remain unchanged.

For a reproducible test using generated data, run:

```bash
bash tests/run_smoke.sh /tmp/bmonth-smoke
```

This checks command dispatch, error propagation, exposure and survival-time
handling, all analysis stages, 100 bootstrap replicates, and the expected output
tables and figures. It also compares scalar estimates/standard errors/P values
with R PheWAS and joint tests with `glm`/`anova`, checks serial versus parallel
results, and exercises resume, corrupt checkpoints and worker failures. Test
outputs are kept in the specified directory.

Stages:

```bash
./bmonth.sh prep
./bmonth.sh phewas
./bmonth.sh assoc
./bmonth.sh med --med-boot 100
./bmonth.sh final
```

Default scalar exposures:

```text
bmonth_cos,bmonth_sin,bmonth_winter_summer,bmonth_spring_autumn,birth_north_100km,birth_east_100km
```

`bmonth_rank` can be added explicitly as a sensitivity analysis:

```bash
./bmonth.sh phewas --exposures bmonth_cos,bmonth_sin,bmonth_rank
```

## PheWAS outputs

- `raw/phewas.csv`: scalar PheWAS results.
- `raw/phewas_joint.csv`: 11-df month, 2-df annual cycle, and month-by-location interaction tests.
- `fig/Fig1_phewas_<exposure>.png`: scalar PheWAS Manhattan plots.
- `fig/Fig1_phewas_month11df.png`
- `fig/Fig1_phewas_cycle2df.png`
- `fig/Fig1_phewas_north_interaction2df.png`
- `fig/Fig1_phewas_east_interaction2df.png`

## Omics mediation

`prot.rds` and `met.rds` are handled in the same exploratory framework as `border`:

1. X -> mediator: BH FDR < 0.05.
2. mediator -> incident outcome, adjusting for X: P < 0.05.
3. coefficient-product + bootstrap for selected candidates.

For `bmonth_cos`, `bmonth_sin`, northing and easting, X is standardized within the complete-case mediation sample. The two binary seasonal contrasts remain 0/1.

A full 2-df seasonal cycle has no single scalar mediation coefficient, so mediation is performed on the sine and cosine components separately. This should be described as exploratory associational decomposition, not causal natural indirect-effect estimation.

## Important interpretation points

- Astrology/zodiac should not be used as the biological explanation. Month-of-birth associations are scientifically plausible through prenatal/perinatal seasonality: photoperiod, temperature, maternal diet, vitamin D, infections, air pollution and other exposures.
- Month is circular. The primary analysis is therefore `month11df` plus `cycle2df`; linear month rank is only sensitivity analysis.
- `birth_north`/`birth_east` are geographic coordinates, not climate variables. A stronger follow-up analysis should link birth location and month/year to historical temperature, daylight, pollution and possibly food availability.
- Month-of-birth effects can be confounded by birth cohort, geography, school-entry rules and secular changes. Consider adding year of birth / birth-cohort sensitivity analyses when the relevant variable is available in `all.rds`.

## Prior literature to compare against

- Boland et al. *JAMIA* 2015: **Birth month affects lifetime disease risk: a phenome-wide method**; introduced a high-throughput Season-Wide Association Study (SeaWAS).
- Didikoglu et al. *Am J Hum Biol* 2020: **Seasonality and season of birth effect in the UK Biobank cohort** (N≈502k).
- Lewis et al. *Chronobiology International* 2021: **Perinatal photoperiod associations with diabetes and chronotype prevalence in the UK Biobank**; explicitly combined time of year of birth with latitude.
- UK Biobank mental-health/neuroimaging work has used a cosine-transformed month-of-birth phenotype, supporting cyclic rather than ordinary linear month coding.
