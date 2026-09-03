# EMS120 analysis pipeline

The pipeline produces the following main figures:

1. **Fig1** — call-volume-adjusted weekly EMS phenotype spectrum.
2. **Fig2** — NLP/transformer phenotyping workflow and validation.
3. **Fig3** — housing-price and spatial validation.
4. **Fig4** — 12-year phenotype gradients.
5. **Fig5** — reopening analyses: total demand, phenotype-specific responses, stage shifts, and high-vs-low interactions.

Supplementary figures contain phone-score construction, on-scene and circadian analyses, confusion/error structure, geographic validation, policy-window analyses, transformer summaries, death-call validation, temporal spectra, robustness analyses, and dose-response results. PNG and `.out.xlsx` filenames use matching figure numbers.

## ML implementation

`f/ems120.py` contains disease-model training and inference. `f/ems120.f.R` contains statistical and plotting helpers. Fig2 derives post-processing summaries from out-of-fold predictions and does not train a second model.

## Run from scratch

```bash
chmod +x ems120.sh
./ems120.sh
```

The shell still performs the CUDA check and will run MacBERT stratified 5-fold OOF CV when required. If the current OOF cache is valid, the existing pipeline may reuse it according to its cache logic.

## Step controls

The accepted processing steps remain `fig1` through `fig6`, plus `figs`; the `fig6` step now emits `FigS2.png` and `FigS2.out.xlsx`.
