# CHARLS AIPFI-CMM reproduction pipeline (v008)

This folder contains a reproducible workflow for the CHARLS paper:

**Associations of cumulative exposure and dynamic trajectories of the combined atherogenic and frailty index with incident cardiometabolic multimorbidity.**

**Transformer-based models for predictingcardiovascular risk in Chinese adults: developmentand validation**

Style is unchanged: one shell entry point, one R analysis file, and outputs split by paper figure/table sections marked with `# 🚩`.

## Files

```text
charls.sh      One entry point with usage()
charls.R       Main R analysis script
0README.md     This file
```

## Default paths

```text
Input  : /mnt/d/data/charls/clean
Output : /mnt/d/analysis/charls
```

Main required input:

```text
Harmonized/H_CHARLS_D_Data.dta
```

Raw Blood/Biomarker files used only to supplement biomarker variables:

```text
2011/Blood_20140429.dta       join by Harmonized ID_w1 <-> raw ID
2011/biomarkers.dta           join by Harmonized ID_w1 <-> raw ID
2013/Biomarker.dta            join by Harmonized ID    <-> raw ID
2015/Biomarker.dta            join by Harmonized ID    <-> raw ID
2015/Blood.dta                join by Harmonized ID    <-> raw ID
```

The file below is deliberately **not** used by default:

```text
2012/Biomarkers.dta
```

Your ID QC showed that `2012/Biomarkers.dta` has a 9-character `ID` and does not safely match the Harmonized wave-1/stable IDs. Therefore it is skipped rather than force-merged.

## Important v005/v006 fixes

v003/v004 caused the wave-1 mismatch because wave-1 2011 raw files were merged using the stable Harmonized `ID`. The ID QC shows the correct wave-1 join is:

```text
Harmonized ID_w1  <->  2011 raw ID
```

Specifically:

```text
2011/Blood_20140429.dta: raw ID exact overlap with H ID_w1 = 11,847 / 11,847 raw IDs
2011/biomarkers.dta    : raw ID exact overlap with H ID_w1 = 13,974 / 13,974 raw IDs
```

v005/v006 uses a fixed explicit ID plan. It does **not**:
- auto-test multiple ID transformations,
- zero-pad IDs,
- drop leading zeros,
- concatenate householdID + personID,
- choose the largest overlap automatically.

If an explicit join has low overlap, the code stops and writes QC.

## Run

```bash
cd /mnt/d/scripts/charls
chmod +x charls.sh

bash charls.sh clean
bash charls.sh prep
bash charls.sh
```

Or with explicit paths:

```bash
CHARLS_CLEAN=/mnt/d/data/charls/clean \
CHARLS_OUT=/mnt/d/analysis/charls \
bash charls.sh prep
```

## Modes

```text
all      prep -> Fig1 -> Tab1 -> Tab2 -> Tab3 -> Fig2 -> Fig3 -> Fig4 -> Fig5 -> Fig6 -> supp
clean    remove previous outputs under CHARLS_OUT
prep     build analytic datasets and QC files
fig1     flow diagram
tab1     baseline characteristics by AIPFI quartile
tab2     baseline AIPFI Cox models
tab3     cumulative AIPFI and trajectory Cox models
fig2     KM curves
fig3     K-means AIPFI trajectories
fig4     spline curves
fig5     subgroup forest plots
fig6     prediction / ROC / calibration / DCA / SHAP-like plots
supp     supplementary sensitivity analyses
qc       print expected QC and output paths
```

## QC files to check after `prep`

```text
qc/explicit_raw_id_plan.xlsx
qc/raw_id_match_qc.xlsx
qc/available_variables_after_raw_merge.xlsx
qc/var_map.xlsx
qc/fi_item_map.xlsx
qc/cohort_counts.xlsx
qc/manifest_file_choice.xlsx
```

Expected key lines in `raw_id_match_qc.xlsx`:

```text
2011/Blood_20140429.dta  : H ID_w1 <-> raw ID, overlap about 11847
2011/biomarkers.dta      : H ID_w1 <-> raw ID, overlap about 13974
2013/Biomarker.dta       : H ID    <-> raw ID, positive overlap
2015/Biomarker.dta       : H ID    <-> raw ID, positive overlap
2015/Blood.dta           : H ID    <-> raw ID, positive overlap
```

If TG/HDL-C still cannot be found, inspect:

```text
qc/debug_available_variables_wave1_TG.xlsx
qc/debug_available_variables_wave1_HDLC.xlsx
qc/var_map.partial_on_error.xlsx
```

Manual variable override is still allowed, for example:

```bash
CHARLS_VAR_W1_TG=raw2011_blood_20140429_newtg \
CHARLS_VAR_W1_HDLC=raw2011_blood_20140429_newhdl \
bash charls.sh prep
```

## Key outputs

```text
dat/charls_waves_2011_2020.rds
dat/baseline_2011_2020.rds
dat/longitudinal_2015_2020.rds

fig/Fig1.png ... fig/Fig6.png
tab/Fig1.out.xlsx
tab/Tab1.out.xlsx
tab/Tab2.out.xlsx
tab/Tab3.out.xlsx
tab/Supplementary.out.xlsx
```

## Method summary

AIP is calculated from TG and HDL-C using molar units:

```text
AIP = log10((TG / 88.57) / (HDL-C / 38.67))
```

FI is constructed from 32 deficits. AIPFI is `AIP * FI`. CMM is defined as coexistence of at least two of diabetes, heart disease, and stroke. Baseline analyses use the 2011/2012 baseline, with follow-up through 2020. Longitudinal analyses use repeated AIPFI at baseline and 2015, then follow-up after 2015.

## Important v006/v007 fixes

- Fixes the `centers_orig[, mean_level := ...]` error by converting the K-means center table from `data.frame` to `data.table` before using `:=`.
- `USE_MICE=1` now requires the `mice` package in the exact R environment used by `charls.sh`; if it is unavailable, the code stops with a clear message rather than silently falling back. Use `USE_MICE=0 bash charls.sh prep` only when you intentionally want median/mode fallback.



## Important v007 fix

- Fixes the Tab2 Cox extraction error: `object 'exp.coef.' not found`.
- The fix normalizes `summary(coxph)$coefficients` and `summary(coxph)$conf.int` column names using `make.names()` before extracting HR, 95% CI, and P values. This is more robust across `survival` package versions.
- Also suppresses a harmless data.table shallow-copy warning in Tab1 by explicitly copying the descriptive table before adding the P column.


## Important v008 fix

`cox_extract()` no longer depends on exact `summary(coxph)$conf.int` column names such as `exp.coef.`. Different `survival`/`data.table` versions may expose these names as `exp(coef)`, `exp.coef`, or `exp.coef.`. v008 extracts HR, lower 95% CI and upper 95% CI by the documented column order of `summary.coxph()$conf.int`, avoiding the v007 error:

```text
Items of 'old' not found in column names: [exp.coef.]
```
