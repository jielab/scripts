# Bundled ArchaicSeeker3 source

This directory contains the minimal runtime source copied from the official
`Shuhua-Group/ArchaicSeeker3.0` GitHub `main` archive downloaded on 2026-08-27.
GU intentionally keeps the model checkpoints outside the scripts tree under
`/mnt/i/refGen/archaic/38/models`.

Two runtime-only compatibility edits are applied locally:

- NumPy 2 scalar conversion uses the first element of the one-element ancestry array.
- Pandas Copy-on-Write uses assignment instead of chained `fillna(..., inplace=True)`.

Upstream: <https://github.com/Shuhua-Group/ArchaicSeeker3.0>

Local source cleanup (2026-09-13): removed earlier, unconditionally overwritten
definitions of `find_introgression_segments` (3),
`inference_and_save_basemodel_overlap` (2), and `SmootherDataset` (1).
The final definitions, signatures, defaults and bodies are unchanged.
