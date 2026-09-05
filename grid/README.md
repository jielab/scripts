
# Methods notes: automatic ancestry

The preferred PRS-CSx replication uses a classifier trained with labelled reference individuals. 
The files supplied with DiscoDivas include projected-PCA weights and four population medians, but not the complete labelled individual-level 1KG PC training table. 
Therefore v9 implements an automatic, reproducible UKB-specific approximation rather than silently treating nearest-center assignment as formal ancestry.

The automatic classifier uses two independent signals:

1. released 1KG AFR/EAS/EUR/SAS centers in the DiscoDivas PC coordinate system;
2. UKB self-reported ethnicity, used only to identify high-purity training anchors when it agrees with the nearest center.

Balanced LDA is trained on projected PC1-PC20. Equal class priors prevent the very large EUR group from dominating. 
Posterior probability below 0.90 produces OTH, retaining ambiguous and admixed participants for DiscoDivas continuum analyses. 
High-confidence, genetic/self-report-concordant participants are used for fine-tuning cohorts in evaluation.

This is methodologically stronger than nearest-center grouping and requires no manually curated ancestry file. 
It is not claimed to be numerically identical to the random-forest classifier trained directly on labelled 1KG individuals in the original PRS-CSx paper. 
A user-supplied validated ancestry file still overrides the automatic procedure when exact external classification is available.

## ARG methods and data preparation

ARG construction is shared reference/target-data preparation and belongs in
`/mnt/d/scripts/0data/refGen.sh make-arg`. The two supported method names are:

- `needle`: the default, used by GRID with UKB phased genotypes. Durable output
  is stored beside the target genotype data, normally
  `/mnt/h/ukbGen/37/arg/{argn,trees}`.
- `tsinfer`: the original GU TRACE engineering backend. Both methods can be
  exported for TRACE with `--format trace`; method-specific trees and sample
  maps are stored under `arg/trace/needle` or `arg/trace/tsinfer`.

`grid.sh arg` and `gu.sh trace` validate the ARG files they consume. ARG data
is generated only by `refGen.sh make-arg`; GU no longer exposes a separate
`arg` module because `gu.sh trace` runs its read-only ARG validation internally.

The default Needle build is an ancestry-balanced pilot. Full UKB inference is
explicit because its time, memory, and disk requirements can be very large:

```bash
bash /mnt/d/scripts/0data/refGen.sh make-arg \
  --method needle \
  --dir-gen /mnt/h/ukbGen/37 \
  --dir-pfile /mnt/h/ukbGen/37/hap \
  --map-dir /mnt/d/data.BIG/refGen/maps/GRCh37 \
  --scratch /home/$USER/grid_arg_scratch \
  --chr 1-22 \
  --full TRUE
```

Use a native Linux/WSL ext4 path for `--scratch`; repeated ARG-Needle
intermediates on `/mnt/c`, `/mnt/d`, or `/mnt/h` are substantially slower.
Benchmark chromosome 22 before scheduling all autosomes. A full UKB run may
need an HPC node and substantially more native-Linux storage than the pilot.

ARG-Needle, rather than `tslmm`, performs GRID's ARG inference. `tslmm` is an
optional supervised ARG-BLUP benchmark and is not required for the primary
zero-shot GRID score. GRID uses local genealogy to learn a SNP-level
transportability prior; it does not require a participant-by-window local
ancestry matrix. GWAS used for performance evaluation must exclude UK Biobank
to avoid optimistic R2/AUC estimates.
