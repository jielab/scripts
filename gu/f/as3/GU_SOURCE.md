# Bundled ArchaicSeeker3 source

This directory contains the minimal runtime source copied from the official
`Shuhua-Group/ArchaicSeeker3.0` GitHub `main` archive downloaded on 2026-08-27.
GU intentionally keeps the model checkpoints outside the scripts tree under
`/mnt/d/data.BIG/refGen/archaic/38/models`.

Two runtime-only compatibility edits are applied locally:

- NumPy 2 scalar conversion uses the first element of the one-element ancestry array.
- Pandas Copy-on-Write uses assignment instead of chained `fillna(..., inplace=True)`.

Upstream: <https://github.com/Shuhua-Group/ArchaicSeeker3.0>
