# GRID ARG-Needle release validation

Release date: 2026-09-04

## Replacement compatibility

- The release was compared with the then-current public `jielab/scripts/grid`
  tree.
- All files from the public tree are represented in the replacement or
  preserved through an explicit legacy copy.
- The existing PRS-CSx source and both bundled DiscoDivas reference assets were
  retained byte-for-byte.
- `grid.legacy.unmodified.sh` preserves the previous public dispatcher.
- The installer was tested against a dummy destination containing an
  uncommitted local-only file; the complete previous directory was retained in
  the timestamped backup.

## Static and synthetic checks

- Bash syntax: all release shell scripts passed `bash -n`.
- Python syntax: all release Python files passed `py_compile`.
- R files were parsed when `Rscript` was available in the build environment.
- The small dispatcher passed its code-level preflight.
- Synthetic integration covered summary-statistic normalization, chromosome
  splitting, PRS-CSx HDF5 LD-score extraction, transportability-table
  construction, cross-validated prior fitting, GRID effect construction, score
  mixing, and evaluation-table creation.
- The ARG-Needle wrapper was checked against the official advanced three-step
  CLI argument names and converter function discovery.

## Limits of validation

No UK Biobank BGEN, phenotype RDS, full PRS-CSx reference panel, or real
ARG-Needle binary was present in the build container. Therefore, this validates
replacement integrity, syntax, dependency checks, and synthetic data flow; it
does not claim that a full UKB chromosome has already completed. The required
first real-data gate is the chromosome-22 pilot described in `RUNBOOK.md`.
