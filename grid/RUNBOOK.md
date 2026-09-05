# GRID v0.2 runbook

## Install as a safe replacement

From the directory containing this release:

```bash
bash grid/install/install_grid_release.sh /mnt/d/scripts/grid --yes
cd /mnt/d/scripts/grid
bash install/install_grid_env.sh
conda activate grid-argneedle
```

The installer first copies the complete existing destination to a timestamped
sibling backup and prints the exact rollback command.

## Reuse/build PCA and global ancestry

```bash
./grid.sh pca
./grid.sh ancestry
```

`ancestry` uses projected PCs, 1000 Genomes centers, and UKB anchors to produce
global AFR/EAS/EUR/SAS/OTH assignments and posterior probabilities. It does not
infer chromosome-local ancestry.

## Run PRS-CSx and zero-shot DiscoDivas

```bash
./grid.sh csx --trait height --dir-gwas /mnt/d/data.BIG/gwas/4grid
./grid.sh disco --trait height
```

Repeat with `--trait ldl` and `--trait t2dm`.

## ARG-Needle chromosome-22 pilot

Use fast native Linux scratch space:

```bash
mkdir -p "$HOME/grid_arg_scratch"

bash /mnt/d/scripts/0data/refGen.sh make-arg \
  --method needle \
  --dir-gen /mnt/h/ukbGen/37 \
  --dir-pfile /mnt/h/ukbGen/37/hap \
  --map-dir /mnt/d/data.BIG/refGen/maps/GRCh37 \
  --ancestry-file /mnt/d/data/ukb/phe/pca/ukb.ancestry.auto.tsv.gz \
  --scratch "$HOME/grid_arg_scratch" \
  --chr 22 \
  --max-individuals 20000 \
  --seed-haplotypes 4000

./grid.sh arg --chrs 22
```

This executes:

```text
balanced keep list
  -> phased BGEN to Oxford HAPS/SAMPLE
  -> map harmonization/QC
  -> ARG-Needle advanced step 1/2/3
  -> .argn to .trees conversion
  -> allele-age and local-genealogy features
  -> ancestry-anchor genealogical summaries
```

The pilot is a software and data-flow check. A paper-level transportability
model must use multiple chromosomes so that held-out-chromosome validation is
available.

## Build and evaluate GRID for the pilot

```bash
./grid.sh grid --grid-action all --trait height --chrs 22
./grid.sh eval --trait height
```

For an RDS with a nonstandard field name:

```bash
./grid.sh eval --trait ldl --phenotype-col YOUR_LDL_COLUMN
```

## Full autosomal analysis

After the chromosome-22 pilot has passed and resource use is acceptable:

```bash
bash /mnt/d/scripts/0data/refGen.sh make-arg \
  --method needle \
  --dir-gen /mnt/h/ukbGen/37 \
  --dir-pfile /mnt/h/ukbGen/37/hap \
  --map-dir /mnt/d/data.BIG/refGen/maps/GRCh37 \
  --scratch "$HOME/grid_arg_scratch" \
  --chr 1-22 \
  --full TRUE \
  --max-individuals 20000 \
  --seed-haplotypes 4000

for trait in height ldl t2dm; do
  ./grid.sh csx --trait "$trait" --dir-gwas /mnt/d/data.BIG/gwas/4grid
  ./grid.sh disco --trait "$trait"
  ./grid.sh grid --grid-action all --trait "$trait" --chrs 1-22
  ./grid.sh eval --trait "$trait"
done
```

Increase `--arg-max-individuals` only after documenting a prespecified scaling
analysis. Keeping one common ancestry-balanced ARG panel across all three
traits avoids constructing a different genealogical prior for each outcome.

## Interpretation of the methods

- `PRS-CSx component`: ancestry-specific posterior SNP effects.
- `Disco-zero`: phenotype-free whole-score interpolation based on global PCA
  distance; this is a deployment baseline, not the phenotype-tuned paper model.
- `GRID`: shared effect plus ancestry-specific residual effect, with residual
  borrowing controlled by a held-out prediction of effect conservation from
  MAF, LD, global distance, allele age and local genealogy.
- `CSx tuned stack`: phenotype-tuned out-of-fold comparator reported separately
  and never used to construct GRID weights.

A full scientific result requires effect/allele harmonization QC, GWAS–UKB
overlap exclusion, multiple-chromosome validation, sensitivity to ARG sample
size and maps, and independent-cohort replication. Passing the dispatcher's
internal checks and the software tests establishes executable consistency, not
biological validity.
