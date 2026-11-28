# METAG · 双端宏基因组分箱与证据复核

## 从原始 FASTQ 开始

将包内的 `metag/` 安装到 `scripts/metag/`，与项目的 `scripts/0f/` 并列。切换前请确认仍指向旧安装路径的 split、refine 与 paper 作业均已结束或取消。将上一次分析的 `/work/sph-huangj/analysis/metag/` 和 `/work/sph-huangj/data/metag/species/` 移出原路径。从这两个空输出位置开始，程序自行创建所需目录。`data/metag/raw/` 提供原始读对，`data/annot/metag/lib/` 提供已经建好的 full 索引、taxonomy 和 `group/optimized/references.tsv`。竞争比对还需要该表指向的原始参考 FASTA。

在 `/work/sph-huangj/scripts/metag/` 执行：

```bash
./metag.sh run --dat-list england,Wxiamen --batch 1-2 --batch-size 10
```

每个队列的第 1、2 批各取最多 10 个完整双端样本，总数最多 40 个。显式分批的 split worker 对整批查询一次 Kraken，再逐样本竞争比对和写 FASTQ。每批完成后，依赖它的 refine 作业自动运行。两个队列完成后汇总：

```bash
./metag.paper.sh run --dat-list england,Wxiamen
```

split 输出格式为 `metag-ranked-v2`，refine 输出格式为 `metag-refine-v2`。首次运行从原始 FASTQ 计算。中断后再次执行同一命令，默认 `--replace FALSE` 复用本版已经发布且输入、参数和记录匹配的样本。尚未发布的样本从查询开始重算。已存在的其他格式结果会触发明确的错误。

## 运行监测

ranked panel 的 MMI 建库和逐样本比对中使用 minimap2 的 `-v 2`，省去遍历索引哈希表以打印统计的步骤。使用短读长 preset 和同一面板内的竞争比对。缓存命中仍需要将 MMI 加载进每个新 mapper 进程。缓存不足时，mapper 会从临时 FASTA 当场建索引。

用 `bpeek JOBID` 查看批量查询计数和随后的逐样本 `RANKED` 行。阶段切换立即输出，长阶段每五分钟输出一次。读对进度每十万对更新，完成时记录最终数。新发布样本的 `summary.json` 保存 `phase_seconds` 和 `mapping.cache`。

| 日志阶段 | 含义 |
|---|---|
| `panel_lookup`、`panel_fasta`、`build_mmi` | 查找缓存或等同一面板的锁、整理参考 FASTA、构建可复用的 MMI |
| `mapper_setup` | 记录选中面板的碱基数、缓存命中与回退原因 |
| `mapper_to_first_sam` | 启动 mapper 到观察到第一行 SAM，包含加载 MMI 或从 FASTA 建索引，以及可能的输出缓冲 |
| `sam_to_first_pair`、`mapping_pairs` | 等待可处理的首对结果、随后处理比对记录。长阶段附 mapper 的累计 CPU 秒、主缺页数、RSS 与存储层累计读取字节（Linux 可用时） |
| `write_bins`、`flush_bins` | 顺序读取原 FASTQ、核对摘要、写 assignments 和压缩分箱，然后排空压缩缓冲并写计数表 |

相邻心跳的 `cpu_s` 差值可与实际经过秒数比较。操作系统可能从页缓存供数据，因此 `read_bytes` 不增长也可能在正常工作。单凭一行日志或一次进程采样无法区分缓慢与死锁。

## 可选的海豹复核

海豹复核需要计算节点上能执行 NCBI BLAST+ 的 `blastn`、`makeblastdb`、`blastdb_aliastool`、`blastdbcmd`。可使用现有 module 或单独的 Conda 环境。项目的 Python 环境继续由 `f/env.sh` 配置。把 BLAST+ 的 bin 加入提交与计算节点共同继承的 PATH，也可通过已有的 METAG_ENV_FILE 加载。以下安装步骤仅在集群缺少 BLAST+ 时使用。

```bash
conda create -n metag-blast -c conda-forge -c bioconda blast
export PATH="$(conda info --base)/envs/metag-blast/bin:$PATH"
for x in blastn makeblastdb blastdb_aliastool blastdbcmd; do command -v "$x"; done
```

在 `scripts/metag/` 运行。可以直接提交整队列。

```bash
./metag.sh refine --dat-list Wxiamen --seal-audit
# 等上面所有复核作业完成
./metag.paper.sh run --dat-list Wxiamen
```

若只想先看读对数最高的样本，可先运行 `./metag.sh refine --sample Wxiamen/SRR39271697 --seal-audit`。等待完成后再提交整队列，默认 `--replace FALSE`。若目录中已发布超过 20 个 Wxiamen 样本，`refine --dat-list` 和 `paper --dat-list` 会纳入当时所有已完成样本，限制范围可重复指定 `--sample Wxiamen/SRR...`。

新样本可运行 `./metag.sh run --dat-list Wxiamen --batch 1-2 --batch-size 10 --replace FALSE --seal-audit`。已有 split 的样本用 refine 更直接。普通 `run`、`split` 和 `refine` 仍按原方式工作。添加 `--seal-audit` 才启用这项全量复核。一般食源 bin 的抽样 `--blast-db` 和这个选项可分别使用。首次建库与 BLAST 请给磁盘留足空间，所需空间随可用的整套参考基因组而变化。

## 自动建的参考库

第一个 worker 从现有 `lib/group/optimized/references.tsv` 或 `lib/group/references.tsv` 选择实际存在的组装，创建单独的 `analysis/metag/refine/references/weddell/`。不同样本共用这一份。进程锁和临时目录保护并发首次构建。源分类库仍由 `1.lib.sh` 管理。

自动面板纳入威德尔海豹、最多四个可用的其他海豹科物种，优先共同祖先最近者，并从事先列出的候选中纳入人、常见家养动物、细菌等背景。参考范围记录在 `seal_panel.json`，原 FASTA 路径、taxid 和解压后 SHA-256 记录在 `references.json`。对于每条命中，`contigs.tsv.gz` 把内部 subject 名映射回原 FASTA accession 和描述。若本地目录没有威德尔海豹组装，补充有来源记录的组装并使用下方自定义面板。若缺少其他海豹科、人、其他非海豹科哺乳动物或细菌背景中的任一类，程序仍报告观察到的命中，并将解释标为 `panel_limited`。常规小面板不能代表整个 nt 数据库。

若怀疑现有海豹组装本身含有污染，可另选独立组装，同时扩充近缘种、污水优势微生物或其他实际疑似来源。建立新版本目录，并保留两次报告供比较。manifest 的字段为 `taxid`、`fasta`、`role`，每个物种一个组装，FASTA 可 gzip，支持相对于 manifest 的路径。每个 taxid 须在项目 taxonomy 中。自定义建库命令请在计算节点运行。为作业提供 frozen panel 的目录，而非某个 FASTA。

```bash
python3 refine/seal_panel.py --manifest seal_refs.tsv --out /work/sph-huangj/data/annot/metag/seal_panel_v2
./metag.sh refine --dat-list Wxiamen --seal-db /work/sph-huangj/data/annot/metag/seal_panel_v2
# 等复核完成
./metag.paper.sh run --dat-list Wxiamen --out /work/sph-huangj/analysis/metag/paper_seal_v2
```

默认每样本保留最近一次发布的海豹复核；需要比较不同版本时，在切换面板之前先用独立的 `--out` 发布旧版本报告，并另存需长期对照的 `refine/seal/` 逐对明细。组装污染报告可从 [NCBI FCS](https://www.ncbi.nlm.nih.gov/datasets/docs/v2/data-processing/policies-annotation/quality/contamination/fcs-contamination/) 核对。

## 判读

- 对原物种 9713 与原属 9712 各自核对全部 R1/R2、读对总数及内容摘要。每个结局均以对应 bin 的全部读对为分母。没有该 bin 时为 `no_bin`。审查这些阳性分箱能评估标签的支持程度，估计漏检还要在其他分箱和未分类读对中引入独立已知标签。
- BLAST 使用 dust，逐端保留 identity ≥95%、单段 query coverage ≥80% 的命中。bit score 在合格最佳分数 98% 以内的命中共同参与判定。达到候选截断、低复杂度、只有一端支持、跨物种并列或同一物种多位置并列分别计数。两个端都指向同一替代物种时记为 `alternative_taxon`，逐端分类和得分见明细。
- 强核证据 `target_nuclear` 要求面板类别齐备、两端均仅支持目标种单一位置、同一已标注核染色体、方向相向、外侧跨度不超过 1,000 bp。该跨度适用于本项目短片段双端文库。线粒体、重复、未定位 scaffold、未能确认是染色体的 contig 单列。保守的 FASTA 标题规则可能低估染色体证据。不同分箱始终分开统计。
- 铺在多少个 10 kb 窗口、去除相同起止后还剩多少片段、相距至少 10 kb 的位置数分别输出。`distributed_panel_support` 要求至少三个相互分离的位置（同一 accession 内的片段中点相距至少 10 kb，不同 accession 单列）。这是面板内用于复查的描述性门槛。共享参考、基因组污染、重复和实验混入仍可能产生一致命中。即使支持真实 DNA，也不能据此推断当地活体或食用来源。
- 论文级验证建议用不同组装重新审查支持度最高的区段，检查是否重叠重复元件或 NCBI FCS 可疑记录，并与同批提取空白、建库空白及组成已知的阳性对照比较；使用独立标签校准物种级 precision、recall、F1，按样本或文库分组。可对有独立真值的样本运行 `./metag.paper.sh run --sample COHORT/SAMPLE --raw-check --truth truth.tsv`。FASTQ 原样保留与物种判定准确率由不同指标衡量。

## 输出

| 输出 | 查看内容 |
|---|---|
| `analysis/metag/refine/seal/Wxiamen/SRR.../seal_pairs.tsv.gz` | 每个候选读对的原始 bin、bin 内从 1 开始的顺序、互斥结局、两端候选 taxids、最佳及最佳非目标得分、可定位目标的 accession 与 1-based 坐标 |
| `seal_hits.tsv.gz` | 返回的各端 BLAST 命中、bit score、identity、覆盖率和坐标。subject 对应面板 contigs.tsv.gz |
| `seal_bins.tsv`、`seal_loci.tsv`、`receipt.json` | 逐样本复核分母、证据强度、窗口、面板范围、输入与软件版本 |
| `analysis/metag/paper/main/ranked/seal_status.tsv` | 每个样本已复核、尚未复核或 split 已变化 |
| `seal_bins.tsv`、`seal_loci.tsv`、`seal_recurrent_loci.tsv` | 全队列分箱结果、各样本原始 accession/窗口、同一面板版本内跨样本重复出现的区段 |
| `analysis/metag/paper/main/index.html` | 筛选队列和样本，查看种级候选与独有核染色体证据的比较图，导出 SVG 并查阅证据表 |

通用食源、群落 AMR 和常规 BLAST 的参数见 [refine/README.md](refine/README.md)。缺少使用 BLAST+ 的环境时，原有 split 和常规无序列搜索的汇总仍然可运行。
