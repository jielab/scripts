# MAHA — final single-CHNS publication pipeline

This package keeps exactly one CHNS analysis script:

- `f/MAHA_chns.R`

There are no parallel CHNS analysis fpreviousers or additional CHNS R scripts. CHNS outputs are written only to:

- `/mnt/d/analysis/maha/chns/` for cohort/intermediate outputs
- `/mnt/d/analysis/maha/` for publication-facing figures/workbooks

## Primary CHNS analysis

The internal primary identifier is `primary`, fixed in `f/MAHA_chns.R`.

The primary analysis is:

- 2011 baseline
- follow-up through 2015
- age >=40 years
- at least 2 diet-record days
- Cox regression
- covariates: age, sex, province, urban/rural residence, log income, education years, current smoking, BMI, and energy intake
- 2019 CHNS master/roster archive for mortality/follow-up linkage
- raw CHNS diet, PEXAM, and PACT inputs used by the intended 2011 analysis

The model identifier is metadata only. It is not printed on Figure 4.

A reproducibility guard stops the CHNS run if the fixed primary analysis no longer gives complete-case N=7,387 with 219 deaths or if DASH/MIND/MEDI do not have HR<1. The guard never changes estimates or switches models; it only prevents accidental publication of a result produced from the wrong inputs/configuration.

## Sensitivity analyses

Only same-cohort sensitivities are retained:

1. all 3 diet-record days
2. without BMI adjustment
3. excluding deaths in the first year
4. excluding baseline disease

No alternate multi-wave primary model is generated or compared in the publication figures.

## Figure 4

`f/MAHA_publication.R` creates the four-panel CHNS figure from the single primary analysis:

- a. Mortality associations
- b. Highest versus lowest quartile
- c. Direct coefficient comparison
- d. Equivalence and Bayes evidence

There is no model-description subtitle above Figure 4.

## Run

```bash
cd /mnt/d/scripts/maha
./maha.sh chns
```

This runs the one CHNS script and then rebuilds the CHNS publication outputs.

To rebuild publication figures from existing cohort outputs:

```bash
./maha.sh chns --pub-only
```
