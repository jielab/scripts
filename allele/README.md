# allele: allelic occupancy and functional geometry

## 正式运行：逐条复制命令

当前机器的 `allele` conda 环境和依赖已检查通过，**不需要再次运行 `install.sh`**，也不需要手动 `conda activate`。只有换电脑、删除环境或调整依赖时才需要重新安装。`./allele.sh --help` 内置同样的逐步命令。

正式输入为 `/mnt/f/annot/gnomad/gnomad.genomes.v4.1.sites.chr*.vcf.bgz`，参考序列为 `/mnt/f/gen/1hgp/38/GRCH38.fasta.gz`（已有 `.fai` 和 `.gzi`，无需解压）。本次开发测试使用 `/mnt/f/gen/1kg/38/vcf`，与正式 gnomAD 结果分开。

在 `/mnt/d/scripts/allele` 目录下，先检查下载是否完成，再按顺序继续。某一步报错先停止，不要继续后续步骤：

```bash
# 检查常染色体 chr1–22 的 VCF 和索引；检查成功后才继续
./allele.sh check-input --vcf-dir /mnt/f/annot/gnomad --chr autosomes

# 全部 chr1–22，不限制区间；输出目录第一次使用时必须没有 prepare 结果
./allele.sh prepare --vcf-dir /mnt/f/annot/gnomad --chr autosomes \
  --fasta /mnt/f/gen/1hgp/38/GRCH38.fasta.gz --out /mnt/d/analysis/allele_gnomad_v4.1

# 选择全部多等位位点及匹配对照；查看输出的 Atlas 查询位点数
./allele.sh select --out /mnt/d/analysis/allele_gnomad_v4.1

# 突变机会模型，不需要 API key
./allele.sh stage1 --out /mnt/d/analysis/allele_gnomad_v4.1

# 自动读取 /mnt/d/files/alphaG.api，无需手动输入或 export key
./allele.sh atlas --out /mnt/d/analysis/allele_gnomad_v4.1

# atlas 成功后继续
./allele.sh geometry --out /mnt/d/analysis/allele_gnomad_v4.1
./allele.sh stage2 --out /mnt/d/analysis/allele_gnomad_v4.1
./allele.sh final --out /mnt/d/analysis/allele_gnomad_v4.1
```

结果位于 `/mnt/d/analysis/allele_gnomad_v4.1/final/`。每一步显式写出相同的 `--out`，因此换终端后也不会误用其他输出目录。`atlas` 中断后再运行同一条 `atlas` 命令会续传；其他步骤成功后直接从下一步继续。`prepare` 拒绝覆盖已有候选表，更换输入或重复完整分析时用新的输出目录。

`check-input` 与 `prepare` 都会检查所选染色体是否齐全、是否还有 `.aria2` 标记、文件大小是否符合固定版本的官方记录、VCF 是否具有 BGZF 结束标记及配套索引。它们不会把已下载的少数染色体默认为全套数据。MD5 由下载脚本完成，这个快速检查不重新扫描数百 GB 文件。默认 gnomAD 分析为 chr1–22；`--chr 1,2,22` 可显式选子集，`--chr all` 才包括 X/Y。性染色体倍性需单独考虑。

完整 gnomAD 的 Atlas 查询可能很大，取决于 `select` 打印的位点数；处理使用全部多等位位点，不会静默抽掉病例。真实 API 已通过单个位点、3 个替换的连接与向量测试；全基因组吞吐和资源消耗尚未验证。具体测试范围见 [tests/README.md](tests/README.md)。

## 可选的本地 1KG 测试

不必重复已经完成的测试。要自己试跑小区间时，使用独立的新输出目录：

```bash
cd /mnt/d/scripts/allele
./allele.sh prepare --vcf /mnt/f/gen/1kg/38/vcf/chr22.vcf.gz \
  --region chr22:20000001-21000000 --out /mnt/d/analysis/allele_1kg_try
./allele.sh select --out /mnt/d/analysis/allele_1kg_try
./allele.sh stage1 --out /mnt/d/analysis/allele_1kg_try
./allele.sh test
```

`--region` 为 1-based、两端包含，须与索引齐全的 `--vcf` 一起使用。不带区间时读取整个文件。`--vcf` 与 `--vcf-dir` 互斥；普通队列目录读取其中的 VCF，gnomAD 目录只选择指定版本和染色体的标准文件名。流程只使用 VCF INFO 中的 AC/AN/AF，不解码个体 GT，也不需要 BAM/FASTQ。

## 配置与测试

修改分析路径时，先执行 `cp -n config.example.env allele.env`，再编辑 `allele.env`。输入目录变量为 `ALLELE_VCF_DIR`；未配置新变量时仍兼容显式设置的旧 `ALLELE_GNOMAD` 输入路径。安装脚本与分析入口读取同一配置，`allele.env` 中的赋值优先于同名环境变量；分析命令的 `--vcf`、`--vcf-dir`、`--fasta`、`--out` 参数可覆盖配置。若配置文件固定了 `ALLELE_OUT`，请将示例输出路径写进配置，或在每一步传入相同的 `--out DIR`。`./allele.sh --help` 会显示当前配置的路径。独立下载脚本 `gnomad_download.sh` 不读取此配置，用自己的 `--out`、`--proxy` 参数。

自定义 conda 可通过 `ALLELE_CONDA_EXE=/path/to/conda` 指定，配置文件可通过 `ALLELE_CONFIG=/path/to/file.env` 指定。环境缺失时会提示先安装；依赖导入失败时检查命令返回非零状态。数据路径和 API key 不影响依赖检查是否通过，但必须在对应分析阶段配置。

`./allele.sh test` 使用人工 VCF/FASTA 和仅在测试进程内替换的 Atlas 客户端，覆盖整个 `all` 命令、断点恢复、向量对应、几何数值和两阶段模型，不调用真实 API，不下载 gnomAD，不把模拟结果写进默认分析目录。**离线测试通过不代表真实 Atlas 授权或真实数据已验证。**

本次修复还包括混合 indel/SNV 记录中 AF/AC 的索引错位、AC/AF 为零的 ALT 被错误计入 occupancy、REF 不一致被静默丢弃、REF 放宽匹配失效、Stage 1 训练内分数泄漏、主 ALT 几何特征直接泄漏 observed 标签，以及未完成 Atlas 缓存被当作完整结果使用。

This pipeline tests a two-stage model of human multiallelic variation:

1. **Mutation generation / opportunity**: among the three possible single-base substitutions at a variable nucleotide, which alternatives are more likely to arise/be observed because of substitution class and local sequence context?
2. **Evolutionary retention**: after accounting for mutation opportunity, does the AlphaGenome molecular-effect geometry help explain whether one, two, or all three alternative alleles are retained in human populations?

The central state variable is **allelic occupancy**: for a reference nucleotide, `O=0..3` counts how many of the three possible alternative bases are observed. In a variant-only VCF, this workflow analyzes `O=1..3`: biallelic, triallelic, and tetra-allelic loci.

## Four research questions

**Q1. Why do some nucleotide positions become multiallelic?** Compare all triallelic/tetra-allelic loci with context-matched biallelic controls. Stage 1 captures mutation opportunity; Stage 2 asks whether molecular tolerance / functional geometry adds predictive information.

**Q2. Are additional alleles functionally redundant or divergent?** AlphaGenome signed molecular-effect scores are treated as a high-dimensional vector. For the secondary observed ALT at a triallelic locus, the pipeline calculates distance to the reference (zero-effect vector), distance/cosine to the primary observed ALT, an affinity score, and perpendicular novelty.

**Q3. Does functional geometry predict which alternative is retained?** A conditional logistic model compares observed versus unrealized alternatives *within the same nucleotide position*. This controls the genomic position by design and adjusts for the Stage-1 mutation-opportunity score.

**Q4. What distinguishes complete mutational occupancy?** Tetra-allelic sites (`REF + 3 ALT`) are the extreme class. The pipeline reports their functional magnitude/diversity and contrasts them with tri- and biallelic sites.

## Why AlphaGenome Atlas fits this design

AlphaGenome Atlas contains precomputed effects for over 9 billion possible human SNVs. The public SDK exposes the Atlas through an API. This pipeline deliberately does **not** download the ~1 PB Atlas. Instead, it queries all three possible alternatives only at selected loci. The default vector uses fixed-dimensional, signed molecular-effect scorers: `ATAC, DNASE, CHIP_TF, CHIP_HISTONE, CAGE, PROCAP`. You may override this list, but gene-centric scorers have variable row structure and are not enabled by default.

Official resources:
- AlphaGenome SDK: https://github.com/google-deepmind/alphagenome
- AlphaGenome documentation: https://www.alphagenomedocs.com/
- AlphaGenome Atlas overview: https://deepmind.google/science/alphagenome/
- gnomAD v4.1: https://gnomad.broadinstitute.org/

AlphaGenome outputs are research predictions, not clinical pathogenicity calls. Respect the AlphaGenome terms of use; Atlas/API access is intended for research and is subject to its current licensing and rate limits.

## Directory layout

```text
allele/
├── allele.sh
├── install.sh
├── gnomad_download.sh
├── environment.yml
├── config.example.env
└── f/
    ├── common.sh
    ├── vcf_inputs.py
    ├── gnomad_v4.1_genomes.json
    ├── parquet_io.py
    ├── extract_sites.py
    ├── select_atlas_sites.py
    ├── atlas_query.py
    ├── geometry.py
    ├── stage1_mutation.py
    ├── stage2_retention.py
    ├── final.py
    └── validate.py
```

## Input details

The analysis requires GRCh38 because AlphaGenome Atlas is queried on GRCh38 coordinates. Do not mix genome builds.

### Choosing a cohort

- **Local 1KG WGS:** the existing GRCh38 chr22 VCF contains 3,202 samples and INFO counts/frequencies. It is suitable for developing and testing the workflow. The [30x collection](https://www.internationalgenome.org/data-portal/data-collection/1000genomes_30x) comprises 2,504 unrelated samples plus 698 relatives; use an unrelated subset or account for relatedness in formal population analyses. After subsetting samples, recompute AC/AN/AF before this pipeline; `prepare` trusts the existing INFO fields.
- **gnomAD genomes:** 76,215 individuals in the [v4 release](https://gnomad.broadinstitute.org/news/2023-11-gnomad-v4-0/). The sites-only VCF provides aggregate counts/frequencies with extensive annotations, not sample GT columns. More samples improve detection of rare additional alleles, but gnomAD is not a methodological requirement. The complete chr22 file is about 8.73 GB and is unnecessary for the local pilot.
- **UKB:** WGS can support the same analysis if available. The local directories inspected here contain `typ`, `hap` and `imp` chip/imputation data. These should not be treated as an unbiased inventory of observed rare third/fourth alleles: array content and the imputation reference panel constrain which variants appear. Cohort size alone does not resolve that limitation.

Occupancy means **observed within the chosen cohort**, not present or absent throughout humanity. Smaller cohorts will miss more rare alternatives; compare sources with attention to sample size, ancestry composition, relatedness, coverage and filtering. Do not mix VCFs from different cohorts or duplicate releases in one input directory.

### gnomAD 独立下载与多机分工

[`gnomad_download.sh`](gnomad_download.sh) 是单文件脚本，内置 chr1–22、X、Y 的 Google Cloud/AWS 地址和官方 MD5（2026-10-15 核对）。下载的是 gnomAD v4.1 genomes **sites summary** 及 `.tbi` 索引，没有个体 genotype；全量 48 个文件约 563 GB。原 aria2 清单仍可单独使用，但新脚本不依赖它。

其他电脑只需复制 `gnomad_download.sh`，安装 `aria2c`，在 Bash/WSL 中运行；不需要 conda、Hail、Google 登录、项目配置或其他代码文件。Ubuntu/WSL 可用 `sudo apt install aria2` 安装下载工具。`install.sh` 现在仅负责环境安装和检查。

按染色体分工，例如分别在两台电脑运行：

```bash
# 电脑 A
bash gnomad_download.sh --proxy http://127.0.0.1:7897 --chr 1,2,3,4,5 --out /mnt/f/annot/gnomad

# 电脑 B；用 --out 改成本机保存目录
bash gnomad_download.sh --proxy http://127.0.0.1:7897 --chr 6,7,8,9,10 \
  --out /data/gnomad/4.1/genomes
```

其余电脑按同样方式分配未下载的染色体。每个染色体会同时加入 VCF 和索引。之前复制的下载脚本默认目录仍为 `/mnt/f/annot/gnomad/4.1/genomes`；现在请显式用 `--out /mnt/f/annot/gnomad`，在该目录运行也可用 `--out .`。默认代理为 `http://127.0.0.1:7897`。可以用 `--proxy none` 直连，`--jobs N` 调整同时下载的文件数（默认 2）。命令行参数优先于 `GNOMAD_OUT`/`GNOMAD_PROXY` 环境变量；脚本不读取 `allele.env`。

单染色体续传或全量下载：

```bash
./gnomad_download.sh --chr 22 --out /mnt/f/annot/gnomad
./gnomad_download.sh --chr all --proxy none --out /mnt/f/annot/gnomad
```

每文件最多 16 个连接、每服务器最多 8 个，完成后自动校验 MD5。Ctrl+C 退出后重跑同一命令可续传，包括先前 curl 下载的 chr22 部分文件。若旧命令仍在下载同一文件，先等它退出；保留中断下载旁边的 `.aria2` 控制文件。文件显示完整大小不代表下载已完成。网络速度仍受线路和代理节点影响。

各电脑下载并校验完成后，将不同染色体的 `.vcf.bgz` 和 `.tbi` 汇总到分析机的 `/mnt/f/annot/gnomad`，再显式指定该输入，例如：

```bash
./allele.sh prepare \
  --vcf /mnt/f/annot/gnomad/gnomad.genomes.v4.1.sites.chr22.vcf.bgz \
  --region chr22:20000001-21000000 \
  --out /mnt/d/analysis/allele_gnomad_chr22_pilot
```

下载 X/Y 不等于自动将它们纳入统计分析；纳入前需处理性染色体倍性。更换队列或区间后使用新的分析输出目录，保持 Atlas 缓存与对应位点一致。

### Reference FASTA

The configured local FASTA is already indexed. For a different, uncompressed GRCh38 FASTA, create its index:

```bash
samtools faidx /path/GRCh38.fa
```

Set `ALLELE_FASTA=/path/GRCh38.fa` in `allele.env`. BGZF-compressed FASTA is also supported when both `.fai` and `.gzi` indexes are present; ordinary gzip FASTA does not support random access.

### AlphaGenome

The API key is stored as one plain-text line in `/mnt/d/files/alphaG.api`, outside the scripts directory. `allele.sh` reads it automatically. Use `ALLELE_API_KEY_FILE` to change this path; an explicitly exported `ALPHA_GENOME_API_KEY` takes precedence. The file is read as data, never executed. `./allele.sh inspect` prints only whether the key is set. Keep the key itself out of source/config files.

## Analysis stages

### `prepare`

Streams coordinate-sorted SNV VCFs, merges same-position same-REF alternatives across records/files, verifies REF against GRCh38 FASTA, and expands every variable locus to all three possible substitutions. Candidate and site tables are written in bounded batches through temporary files, committed only after a successful scan. Unsorted input and existing output tables are explicit errors. Important columns include `observed`, `occupancy`, `af`, trinucleotide `context`, transition/transversion status, and CpG-transition flag.

AF/AC retain the original ALT indices for mixed indel/SNV records. Split records at the same position/REF are merged, preserving additional alleles. Known zero-count alternatives and loci with an absent reference allele (summed ALT AF=1) are excluded. ALT AF sums above one are rejected; reference AF and entropy remain missing when ALT frequencies are incomplete. REF mismatches stop the run; QC counts are in `work/candidates.parquet.qc.json`. Source file paths, sizes, header sample counts and the requested region are recorded in `work/candidates.parquet.inputs.json`. A zero header sample count denotes a sites-only VCF, not a zero-person source cohort. An existing FASTA `.fai` index is required.

### `select`

Keeps every locus with `occupancy >= 2` and samples context-matched biallelic controls. Scans the compact site table in batches; random-priority sampling retains only the needed controls per stratum instead of loading all candidates. Default control ratio is 2:1. This controls Atlas cost while preserving the main comparison.

Sparse context strata are filled only from controls with the same REF. The achieved ratio can be smaller if controls are insufficient. No multiallelic loci is an explicit error. A ratio of zero is allowed but makes the site-level binary model unavailable.

### `atlas`

Queries every possible ALT at the selected loci. Scores are stored as a Zarr matrix and aligned using variant IDs returned by Atlas. The analysis fails rather than silently using a scorer whose output is not fixed-dimensional per variant.

Only subsets of the six default signed scorers are supported. Response matching uses chromosome, position, REF and ALT, including `1`/`chr1` normalization; duplicate, missing or nonfinite responses are errors. Each chunk commits a checkpoint. Re-running resumes from the last committed chunk; a complete cache needs no network calls. Changed sites, scorers or SDK version require a new output directory, preserving the old cache. `geometry` rejects incomplete caches. The manifest records the SDK version, request hash and row counts; it never stores the API key.

### `geometry`

Feature-wise RMS scaling is followed by equal weighting of AlphaGenome modalities, preventing a modality with many tracks from automatically dominating Euclidean distance. For a candidate vector `v` and the primary observed ALT vector `b`:

```text
magnitude       = ||v||
distance        = ||v-b||
cosine          = v·b / (||v|| ||b||)
affinity        = (||v|| - ||v-b||) / (||v|| + ||v-b||)
novelty         = perpendicular distance from the A->B functional axis / ||v||
```

`affinity > 0` means the candidate is functionally closer to the primary observed ALT than to the reference zero-effect vector. The `primary ALT` is operationally defined as the observed ALT with the highest present-day AF; it is **not** claimed to be the historically first derived allele. Ancient DNA can later replace this with temporal ordering.

For exactly zero vectors, cosine is defined as 0; novelty is 0 for a zero candidate and 1 for a nonzero candidate with a zero anchor. These are numerical conventions, not evidence for a biological direction.

### Stage 1: mutation opportunity

A chromosome-blocked model predicts `observed vs unrealized` substitution using only mutation-context variables. It is intentionally called **mutation opportunity**, not mutation rate: because a variant-only VCF does not contain monomorphic positions, the estimand is conditional on the locus already being variable.

Stage 2 receives out-of-fold mutation scores for O<3. With multiple chromosomes, folds are chromosome-blocked; a single-chromosome pilot falls back to stable hashed locus folds and records this in `stage1_metrics.json`. All three alternatives at a locus stay together. O=3 scores come from the full model, which is trained on O<3. Balanced class weights make the scores relative, not calibrated probabilities. Repeated feature/outcome patterns are aggregated within each fold; frequency weights reproduce expanded-data logistic training without sampling away observations. Numeric scaling uses the same frequency weights. Predictions are streamed to compact tables; these retain candidate IDs, occupancy, scores and validation provenance, while the original attributes remain in `work/candidates.parquet`.

### Stage 2: evolutionary retention

Two main models are generated:

1. **Within-locus conditional logistic regression** has two outputs. `stage2_conditional_all_candidates.tsv` compares all three ALT candidates at O=1/2 loci, using mutation opportunity and magnitude. `stage2_conditional_candidate.tsv` focuses on O=2 loci and compares the secondary observed ALT with the unrealized ALT, using mutation opportunity plus magnitude/affinity/cosine/novelty. The observed primary ALT is excluded from this second model: its affinity/cosine=1 and novelty=0 encode its observed status by construction. This is an exploratory association conditional on an AF-selected observed anchor, not a causal test of evolutionary retention.
2. **Site-level multiallelic model**: `occupancy >=2` vs `occupancy=1` using mutation-opportunity summaries and functional tolerance/diversity across all three theoretical alternatives.

Scaling is fitted within each site-model validation fold. Single-chromosome or one-class validation is explicitly marked unavailable (null AUC and NaN held-out predictions), never replaced by a fabricated baseline prediction. Nonidentifiable/nonconverged conditional models are reported as unavailable in metrics and have empty coefficient tables. The site AUC remains exploratory: global geometry scaling and upstream cross-fitting are not a nested end-to-end validation, and context-selected case/control results are not population-calibrated.

The code also writes a focused table for the secondary observed ALT at triallelic sites and a table of tetra-allelic extreme loci.

## Primary outputs

```text
work/candidates.parquet
work/sites.parquet
work/atlas_sites.parquet
atlas/vectors.zarr
atlas/rows.parquet
results/geometry_candidates.parquet
results/geometry_sites.parquet
results/stage1_candidates.parquet
results/stage2_conditional_candidate.tsv
results/stage2_site_multiallelic.tsv
results/stage2_triallelic_secondary.tsv
results/stage2_tetra_sites.tsv
final/Fig1_occupancy.png
final/Fig2_observed_unobserved.png
final/Fig3_triallelic_affinity.png
final/Fig4_multiallelic_model.png
```

## Interpretation guardrails

- Multiallelism is not automatically balancing selection. Recurrent mutation and sample size can generate many multiallelic loci.
- `unobserved` means absent from the analyzed population resource, not biologically impossible or globally extinct.
- AlphaGenome is a molecular-effect prediction model. A large effect vector is not synonymous with deleterious fitness or clinical pathogenicity.
- Present-day AF does not establish allele age or order of emergence.
- The first paper can focus on present-day global occupancy and functional geometry. Ancient DNA should be added as an independent temporal layer only after the static model is stable.
- Tetra-allelic calls deserve stringent QC because an apparent fourth base is especially sensitive to sequencing/calling artefacts.

## Runtime scope

Extraction, Stage 1, occupancy summaries and Atlas-vector normalization use bounded batches. Selection, Atlas request metadata, geometry output tables and Stage 2 still scale with the number of selected loci. A live Atlas probe succeeded for one public 1KG locus (three alternatives, six scorers, 3,763 features). Full-genome gnomAD and sustained real Atlas throughput have not been benchmarked; inspect the selected-locus count before querying the API. Regional and full-chromosome 1KG validation does not establish whole-genome resource requirements. Do not combine overlapping VCF cohorts: duplicate allele records retain the largest AF rather than pool counts. Variant-only inputs cannot estimate O=0 or establish that absent alleles are biologically impossible.

## Recommended next extensions

1. Add mutation-rate maps or de novo mutation tables so Stage 1 becomes a calibrated biochemical mutation-opportunity model rather than a context proxy.
2. Repeat occupancy in 1000 Genomes/HGDP/SGDP and test whether the same secondary ALT appears across ancestries.
3. Add ancient DNA to estimate emergence and persistence times of the second/third ALT.
4. Add direct UKB phenotype vectors as orthogonal validation of AlphaGenome functional geometry, without redefining the primary model around PheWAS significance counts.
