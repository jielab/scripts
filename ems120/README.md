# EMS120 2026-08-12 — final figure renumbering

This package is based directly on the user's 2026-08-12 01:54 code set.

## What changed

### Main figures

1. **Fig1** — unchanged: call-volume-adjusted weekly EMS phenotype spectrum.
2. **Fig2** — NEW: **NLP / transformer phenotyping pipeline and validation** (3 x 2):
   - a. EMS text-to-phenotype workflow
   - b. expert-labelled OOF cohort composition
   - c. Keyword vs MacBERT overall OOF performance with bootstrap 95% CIs
   - d. per-phenotype F1 dumbbell / gain
   - e. largest OOF misclassification pathways
   - f. confidence-threshold utility (coverage, accuracy, Macro-F1)
3. **Fig3** — former Fig2: housing-price / spatial validation.
4. **Fig4** — former Fig3: 12-year phenotype gradients.
5. **Fig5** — redesigned to focus on **reopening**:
   - a. total EMS demand rebound
   - b. full-reopening phenotype-specific response
   - c. phenotype shifts across reopening stages
   - d. high-vs-low differential response / interaction
The former **Fig6** (on-scene time and circadian mechanism) is now **FigS2**.

The former PHSM panels from main Fig5 are moved into **FigS9**, together with detailed policy-window sensitivity analyses.

### Figure layout and supplementary numbering

- **FigS1** is unchanged.
- Former **Fig6** is **FigS2**.
- **Fig2** is a six-panel 3 x 2 layout with equal-width columns: overall/fold stability, precision-recall/threshold utility, and error pathways/calibration.
- The confusion matrix and largest error-pathway panels form the new **FigS4**.
- The previous **FigS4–FigS12** are renumbered **FigS5–FigS13**, respectively.
- New **FigS8** retains only the keyword-composition panel from the former FigS7; its left line-only panel was removed.
- PNG and `.out.xlsx` filenames follow the same mapping.
- TEST.Fig1 and TEST.Fig5 are no longer emitted.

## ML safety boundary

`f/ems120.py` and `f/ems120.f.R` are copied unchanged from the supplied latest version. The existing `ml_phone`, `ml_geo`, `ml_dx`, stratified 5-fold OOF training, cache logic, and model files are therefore not rewritten by this redesign.

The new Fig2 computes additional **post-processing summaries in R from existing OOF predictions**; it does not train a second model.

## Run from scratch

```bash
chmod +x ems120.sh
./ems120.sh
```

The shell still performs the CUDA check and will run MacBERT stratified 5-fold OOF CV when required. If the current OOF cache is valid, the existing pipeline may reuse it according to its cache logic.

## Step controls

The accepted processing steps remain `fig1` through `fig6`, plus `figs`; the `fig6` step now emits `FigS2.png` and `FigS2.out.xlsx`.
