# GRID v1 method implemented in this release

For SNP j and discovery populations a and b, the observed cross-population effect heterogeneity is

    H_jab = log[1 + max((beta_ja-beta_jb)^2/(SE_ja^2+SE_jb^2)-1, 0)].

The baseline model contains mean/difference in MAF, mean/difference in ancestry-specific LD score,
global population distance in the fixed 1000 Genomes PCA coordinate system, and sampling uncertainty.
The evolutionary extension adds ARG mutation age, mutation-branch length, local root time, local
between-population branch divergence, lineage breadth, and lineage entropy. Coefficients are learned
with ridge regression and blocked out-of-fold prediction (leave-one-chromosome-out when at least three
chromosomes are available; genomic-block folds for a one-chromosome pilot).

The held-out predicted heterogeneity gives a conservation prior

    c_j = exp[-max(Hhat_j,0)/2].

Let beta_jk^CSX be the PRS-CSx posterior effect in population k and beta_j^shared the sqrt(N)-weighted
mean of aligned population effects. GRID uses

    beta_jk^GRID = c_j beta_j^shared + (1-c_j) beta_jk^CSX.

Thus a variant predicted to have a conserved effect is transported toward the shared effect; a variant
with weak transportability retains its ancestry-specific PRS-CSx effect. A missing evolutionary prior
sets c_j=0, so GRID falls back to PRS-CSx rather than imposing unsupported sharing.

The primary implementation is population-level and locus-specific. UKB high-confidence AFR/EAS/EUR/SAS
anchors are used to characterize local genealogies. Participant-specific, chromosome-averaged GNN
affinities are optional and are not required for the primary score. No UKB phenotype is used to construct
GRID weights. Phenotype-tuned stacking is evaluated only as a separately labelled comparator.
