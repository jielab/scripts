Source: https://github.com/YunfengRuan/DiscoDivas
Commit: ee8e5d996e6fc6ffaf04255e1cf0fa8ecce68524
License: MIT (LICENSE alongside source).

The official distance-matrix inversion, distance interpolation, shrinkage and
PCA residualization are retained. Two local correctness fixes:
- Label target PC columns using the number actually selected (`npca`).
- Write IDs from the final merged `dat`, keeping each score attached to its ID.
  The upstream `IID` variable can refer to a different ordering/subset.

The wrapper validates finite PCs/scores, unique IDs, complete PRS sample overlap,
center order, positive distances and finite output before publishing.
