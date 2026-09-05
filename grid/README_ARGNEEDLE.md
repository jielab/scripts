# GRID: allele age and local genealogy for cross-ancestry PRS

This release implements the first testable version of the proposed article:

> **Allele age and local genealogy predict cross-ancestry effect transportability beyond LD and global ancestry**

The code is research software. It is designed to test whether ARG-derived mutation age and local genealogy
improve cross-ancestry effect-transportability prediction beyond MAF, ancestry-specific LD, and a fixed
1000 Genomes PCA distance. It is not a clinically validated risk score.

## 1. Dispatcher

```bash
./grid.sh csx|disco|pca|ancestry|arg|grid|eval|all [options]
```

The top-level `grid.sh` is only a dispatcher. Implementations are under `f/`. The previous PCA and ancestry
implementation is retained as `grid.legacy.sh`; the unmodified former dispatcher is also preserved as
`grid.legacy.unmodified.sh`.

## 2. Required input layout

### Multi-ancestry GWAS

Default directory:

```text
/mnt/d/data.BIG/gwas/4grid/
```

Required flat filenames, matched case-insensitively:

```text
height.AFR.gz  height.EAS.gz  height.EUR.gz  height.SAS.gz
ldl.AFR.gz     ldl.EAS.gz     ldl.EUR.gz     ldl.SAS.gz
t2dm.AFR.gz    t2dm.EAS.gz    t2dm.EUR.gz    t2dm.SAS.gz
```

Each file needs SNP/rsID or GRCh37 chromosome-position, effect and other alleles, beta or OR, SE or P,
and preferably EAF plus N. Case-control files can provide N_CASE and N_CONTROL. If N cannot be recovered,
pass, for example:

```bash
--n-gwas 'AFR=100000,EAS=200000,EUR=1000000,SAS=150000'
```

**Every GWAS used to evaluate UK Biobank must exclude UK Biobank.** Renaming a file does not establish
sample independence; verify it from the source documentation.

### Target genotypes for PRS scoring

Default:

```text
/mnt/d/data/ukb/gen/typ/chr1.{bed,bim,fam} ... chr22
```

or equivalent `chrN.{pgen,pvar(.zst),psam}` files.

### UKB phased genotyped data for ARG-Needle

The user's Windows location

```text
H:\ukbGen\37\hap
```

is seen from WSL as

```text
/mnt/h/ukbGen/37/hap
```

Accepted per chromosome:

```text
ukb22438_c1_b0_v2.bgen
ukb22438_c1_b0_v2_s*.sample
...
ukb22438_c22_b0_v2.bgen
ukb22438_c22_b0_v2_s*.sample

or the phased PLINK 2 representation:

```text
chr1.{pgen,pvar.zst,psam} ... chr22.{pgen,pvar.zst,psam}
```
```

These must represent phased Field 22438 data. The ARG builder exports a selected sample to Oxford HAPS while
preserving sample order, interpolates the map onto every retained variant, and runs the official
three-stage large-sample ARG-Needle workflow.

### Genetic maps

Provide one **GRCh37 sex-averaged cumulative genetic map** for each chromosome. Common map layouts are
accepted, including position/rate/cumulative-cM and four-column PLINK maps. Set either:

```bash
--arg-map-dir /path/to/GRCh37/maps
```

with recognizable chromosome filenames, or an explicit pattern:

```bash
--arg-map-pattern '/path/to/maps/genetic_map_GRCh37_chr{chr}.txt.gz'
```

Do not use GRCh38 maps with the GRCh37 UKB BGENs.

### Existing GRID reference assets

The release retains:

```text
csx/                                      modified PRS-CSx source
DiscoDivas/files/g1k_hm3_maf5_woamb_wolr.pca.weight
DiscoDivas/files/med.g1000.4pop.tsv
```

The existing PRS-CSx HDF5 panels remain required at the default `/mnt/d/data.BIG/refLD/csx`. GRID extracts
per-SNP ancestry-specific LD scores directly from those HDF5 matrices; no second LD reference download
is needed when all four panels are already present.

## 3. Software

Create the environment:

```bash
bash install/install_grid_env.sh
conda activate grid-argneedle
export PYTHONPATH=$HOME/software/arg-needle/arg-needle-scripts:${PYTHONPATH:-}
export ARG_NEEDLE_HOME=$HOME/software/arg-needle/arg-needle-scripts
```

The environment contains PLINK 2, Python numerical packages, tskit, and the R packages used by evaluation.
ARG-Needle includes compiled code, so its upstream build may still require following the README matching
the exact checked-out release. The installer reports failure unless both `arg_needle` and
`arg_needle_lib` import successfully.

`tslmm` is optional. It can be used later as a supervised ARG-BLUP benchmark, but it does not infer an ARG
and is not required by the primary GRID workflow.

## 4. Storage and compute safeguards

ARG-Needle input and intermediate files are large. `--arg-scratch` is mandatory and should point to a
native Linux/WSL ext4 filesystem, not `/mnt/c`, `/mnt/d`, or `/mnt/h`. Reading the source BGENs from
`/mnt/h` is unavoidable, but repeatedly writing HAPS/ARG intermediates to DrvFS is much slower.

The default is a balanced **20,000-person pilot**. Full UKB inference is never inferred from the sample
count and requires the explicit flag `--arg-full TRUE`. Measure time, peak RAM, and disk size on
chromosome 22 first. A full 22-chromosome run may require an HPC node and very large native-Linux storage;
the release does not claim that it will fit a particular desktop merely because the pilot succeeds.

## 5. Recommended first run

### Core checks and existing modules

```bash
cd /mnt/d/scripts/grid
./grid.sh pca
./grid.sh ancestry
./grid.sh csx --trait height
./grid.sh disco --trait height
```

### ARG-Needle chromosome-22 pilot

```bash
bash /mnt/d/scripts/0data/refGen.sh make-arg \
  --method needle \
  --dir-gen /mnt/h/ukbGen/37 \
  --dir-pfile /mnt/h/ukbGen/37/hap \
  --map-dir /mnt/d/data.BIG/refGen/maps/GRCh37 \
  --ancestry-file /mnt/d/data/ukb/phe/pca/ukb.ancestry.auto.tsv.gz \
  --scratch /home/$USER/grid_arg_scratch \
  --chr 22 \
  --max-individuals 20000 \
  --seed-haplotypes 4000

./grid.sh arg --chrs 22
```

By default the durable output is `/mnt/h/ukbGen/37/arg`, beside the target
genotype directories. If `--arg-dir` is changed, use the same location later:

```bash
./grid.sh grid --trait height --chrs 22 \
  --arg-trees-dir /custom/arg/trees
./grid.sh eval --trait height
```

The one-chromosome pilot uses genomic-block cross-validation only as a software/method smoke test. The
article's primary analysis must use leave-one-chromosome-out predictions across autosomes.

### Full method after the pilot

```bash
bash /mnt/d/scripts/0data/refGen.sh make-arg \
  --method needle --dir-gen /mnt/h/ukbGen/37 \
  --dir-pfile /mnt/h/ukbGen/37/hap --chr 1-22 --full TRUE \
  --map-dir /mnt/d/data.BIG/refGen/maps/GRCh37 \
  --scratch /native-linux/grid_arg_scratch

for t in height ldl t2dm; do
  ./grid.sh csx --trait "$t"
  ./grid.sh disco --trait "$t"
  ./grid.sh grid --trait "$t" --chrs 1-22
  ./grid.sh eval --trait "$t"
done
```

The dispatcher also supports `all`, but the staged commands above are preferable until chromosome 22 has
passed all QC checks.

## 6. Main outputs

```text
ARG_ROOT/panel/arg.panel.tsv                   selected UKB ARG panel
ARG_ROOT/argn/chrN.argn                        ARG-Needle output
ARG_ROOT/trees/chrN.trees                      tskit tree sequence
ARG_ROOT/trees/chrN.variants.tsv.gz            allele age/local genealogy features
OUTPUT/TRAIT/scores/csx.tsv.gz                 four PRS-CSx scores
OUTPUT/TRAIT/scores/disco_zero.tsv.gz          zero-shot distance-weighted comparator
OUTPUT/TRAIT/grid/transport.tsv.gz              pairwise transport observations
OUTPUT/TRAIT/grid/model/model_cv.tsv             baseline vs evolutionary blocked CV
OUTPUT/TRAIT/grid/model/variant_conservation.tsv.gz
OUTPUT/TRAIT/grid/weights/                       shared and population GRID weights
OUTPUT/TRAIT/scores/grid.tsv.gz                  GRID scores
OUTPUT/TRAIT/eval/performance.tsv                overall/by-ancestry evaluation
OUTPUT/TRAIT/eval/GRID_evaluation.xlsx           review workbook
```

Participant-level score and prediction files contain UKB identifiers and must remain in an approved secure
environment.

## 7. Replacing and rolling back the current directory

From the unpacked release:

```bash
bash install/install_grid_release.sh /mnt/d/scripts/grid --yes
```

The installer moves the complete current directory to a timestamped backup before copying the release and
writes an exact rollback command to `INSTALL_RECORD.txt`. This protects local uncommitted changes even
though the release is based on the current public GitHub `grid/` tree.
