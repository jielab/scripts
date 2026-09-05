# GRID 2026.09.04 ARG-Needle release

- Replaced the large top-level script with a small dispatcher.
- Added independent `disco`, `arg`, and `grid` modules.
- Preserved the previous PCA/ancestry implementation as `grid.legacy.sh` and retained the bundled PRS-CSx and DiscoDivas reference assets.
- Added flat `[trait].[race].gz` discovery for height, LDL, and T2DM.
- Added robust GWAS normalization, PRS-CSx scoring, and zero-shot DiscoDivas comparison.
- Added a three-step large-sample ARG-Needle workflow for phased UKB Field 22438 BGENs.
- Added dated-local-genealogy feature extraction from tskit tree sequences.
- Added MAF/LD/global-PCA baseline versus allele-age/local-genealogy transportability modelling with blocked out-of-fold prediction.
- Added GRID effect construction and target scoring without phenotype-fitted ancestry-PRS weights.
- Added zero-shot and separately labelled phenotype-tuned evaluation.
- Added preflight, installation templates, synthetic tests, and overwrite/rollback installer.

This is research software. The evolutionary prior is a hypothesis to be tested, not a clinically validated PRS.
