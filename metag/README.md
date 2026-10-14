# 环境宏基因组整理与物种分流 v2

本目录是在 `metag_split_20261014` v1 代码上审核、修改的完整版本。
目标：ENA 清单 → `raw/<sample_key>/<run>` → fastp → 同一 Kraken2 库竞争分类 → 互斥 FASTQ 分流。
病毒按**病毒自身的基因组类型**分到 `virus_DNA/`、`virus_RNA/`、`virus_unknown/`。
按 `le8` 的入口方式整理：所有命令从 `metag.sh` 调用，内部模块统一放在 `f/`。

```text
metag/
├── metag.sh
├── install.sh
├── README.md
├── environment.yml
└── f/
    ├── conda_env.sh
    ├── conda-linux-64.lock
    ├── metag.py
    ├── datasets.py
    ├── taxonomy.py
    ├── splitter.py
    ├── virus_metadata.py
    ├── k2_compat.py
    ├── build_db.sh
    ├── run_full.sh
    ├── align_candidates.sh
    └── datasets.example.tsv
```

`tests/`、验证产物、重复实施报告和旧说明已清理。原始实施依据为用户当时提供的
`codex.txt` 与 `metag_split_20261014` 代码包。
数据、数据库、冻结映射和结果保存在各自的工作目录，不放入代码目录。

## 换电脑：安装与准备

把整个 `metag/` 代码目录复制到新电脑，在 **Linux x86_64 或 WSL2 x86_64** 中执行。
先把已下载的原始 FASTQ 放到 `/mnt/f/metag/各数据集名/`；路径不同就修改第二条命令的 `--data-root`。
以下只安装环境和准备参考、清单，不运行全量分析。参考下载和建库可很大，成功后再按下节试跑。

```bash
(
set -e

# 1. 安装并检查软件环境；不需要手动 source 或 activate。
./install.sh

# 2. 整理六个已下载数据集，下载/建立共享参考及病毒映射，并检查覆盖。
./install.sh --prepare --data-root /mnt/f/metag --datasets emp,wbe,长江,厦门污水,全球污水,纳米布沙漠
)
```

`./install.sh --check` 只检查已安装的环境；`./install.sh --dry-run --prepare` 只显示计划。
失败后可重跑相同命令，成功的参考步骤与有效下载缓存会复用。
若本机代理失效，可给安装命令加 `--direct`；只影响该脚本及子进程，无需手动 unset。
本机参考任务已经在运行，目前不要再次执行 `--prepare`。

脚本优先复用已有 Anaconda/Miniforge/Miniconda；没有 Conda 时下载固定的 Miniforge **26.7.2-0**，
验证 SHA256 后安装到 `~/miniforge3`，不修改 shell 启动文件。
`f/conda-linux-64.lock` 固定本次 **181 个软件包**的版本、build 和下载 MD5，包含 Python 3.11.17、
Kraken2 2.17.1、fastp 1.4.0；`environment.yml` 保留可读依赖说明。安装默认使用锁文件，避免下次重新求解出不同版本。
已有环境与锁文件不一致时停止，不自动升级正在使用的环境。
非标准 Conda 路径可通过 `METAG_CONDA_ROOT` 指定，之后调用 `metag.sh` 时也需保留该变量。

`--reference-root` 可单独指定共享参考父目录；默认是 `DATA_ROOT/reference`。
建库默认 16 线程、4 下载线程、24 GiB hash 上限，可通过 `--threads`、`--download-threads`、`--max-db-gb` 修改。
后续分析若使用非默认数据/参考路径，也要相应设置 `--data-root`、`--db` 和 `--virus-map`。
ICTV 固定为 `VMR_MSL41.v1.20260729.xlsx` 并校验 SHA256；accession 缓存保留原始响应与校验信息。

**软件环境复现与参考数据复现是两件事。** 新建库会获取当时的 NCBI 参考与 taxonomy，不能保证与旧库逐字节一致。
要复现旧分析，应一起迁移原始 FASTQ、人工整理的 metadata，以及整个已完成的 `reference/`，
包括索引、同库 taxonomy、ICTV 文件、accession 缓存、病毒映射和来源收据；迁移后重新审计。
数据路径变化时不要原样搬用指向旧路径的 `raw/` 绝对软链接，应从原始 FASTQ 重建；
给准备命令加 `--refresh-metadata` 会备份并重建清单，随后须从备份合并人工注释，再执行 pilot。
`install.sh` 不下载新的原始样本；ocean/sub 下载完成后再纳入，neon 仍需核对文件映射。

## Usage examples：参考准备完成后运行

`metag` 环境已安装，emp、wbe、长江、厦门污水、全球污水、纳米布沙漠的清单与 raw 软链接已准备完成。
`metag.sh` 自动激活环境，保留当前代理设置，无需手动 source、activate 或 unset。
**当前完整共享参考仍在下载/建库，正式库审计和真实数据分类试跑尚未完成。**
下面列出参考准备完成后要执行的命令，与 `./metag.sh --help` 一致。
**准备完成后，完整复制第 1–3 步的整个代码块。** 括号内任一步失败就停止，不会关闭当前终端。
第一步会检查参考；未完成时停止，不会重复启动正在运行的建库任务。

```bash
(
set -e

# 1. 选择六个已下载的数据集，并检查共享参考是否准备完成。
DATASETS=emp,wbe,长江,厦门污水,全球污水,纳米布沙漠
./metag.sh datasets "$DATASETS"
REF=/mnt/f/metag/reference
DB="$REF/kraken2_broad_v1"
MAP="$REF/virus_map_v1/virus_genome_types.sqlite"
if [[ ! -f "$DB/metag_build_state/all.complete" || ! -s "$DB/hash.k2d" || ! -s "$MAP" ]]; then
  echo '共享参考尚未准备完成，请等待；不要重复启动 build-db。' >&2
  exit 1
fi
./metag.sh db-audit --db "$DB" --virus-map "$MAP" --threads 16

# 2. 六个数据集逐个试跑；每个最多 2 个 run，每输入流最多 100000 records。
./metag.sh pilot "$DATASETS" --db "$DB" --virus-map "$MAP" --threads 16 --execute

# 3. 试跑全部成功后，全量处理这六个数据集。
./metag.sh full "$DATASETS" --db "$DB" --virus-map "$MAP" --threads 16 --execute
)
```

**ocean 和 sub 均在下载，前三步只处理上述六个数据集。两个目录均下载完成后，才运行以下命令，并先完成共享参考准备。**
若人工修改过 samples.tsv，刷新后需从备份合并人工注释，再试跑。

```bash
(
set -e

# 4. 两个目录均下载完成后，重建清单并整理 raw；旧 metadata 会先备份。
./metag.sh prepare ocean,sub --refresh --execute
./metag.sh pilot ocean,sub --execute
./metag.sh full ocean,sub --execute
)
```

`manifest` 建立 FASTQ 清单，记录 run、样本、单双端、ENA 大小/MD5 和下载状态；
`prepare` 将建清单和 raw 软链接整理串联起来，已有清单默认保留。
`pilot` 完整校验少量 run 后试跑，结果位于各目录的 `pilot_virus_v2/summary/`；
`full` 完整校验本地清单后进行 QC、分类、分流、汇总，结果位于 `results_v2/summary/`。
清单检查通过不等于已经通过 MD5：实际 MD5 校验由 pilot/full 执行。

批量参数支持英文逗号和中文目录名；各目录独立保存结果，按列表顺序运行。
某个数据集失败会继续其余已选项，最后报告失败项并返回非零状态。
`prepare/pilot/full` 不带 `--execute` 时只打印计划。

本次检查有 9 个 FASTQ 数据目录，`0metadata`、`url` 是辅助目录。
sub 仍有 `.aria2`，neon 有 142 个 submitted 文件名（例如 `4637809.3.050.upload.fastq.gz`），
需先确认 run 归属及处理、合并、配对关系，再加入处理列表。ocean 仍在下载。
`datasets` 显示实时状态；FASTQ 路径数包含 raw 链接，不代表独立样本数。

共享参考默认使用上述 DB/MAP 路径；`--data-root` 仅修改数据集父目录。
脚本自动激活环境；上述变量在括号代码块内设置，整块复制即可。
`prepare --refresh` 会先备份再重建清单，不能用来保留人工注释；元数据补充用 `metadata-upgrade`。
全量仅处理本地清单，不自动下载项目里的其他样本。更新数据库或映射须使用新版本目录。

## 中断与重试

建库时按 Ctrl+Z 会暂停进程，暂停的进程仍持有建库锁。若提示锁被占用，在启动建库的原终端
用 `jobs -l` 查看，`fg` 恢复；不要删除锁文件或重复启动。需要终止时先恢复到前台，再按 Ctrl+C。
建库任务已停止且没有索引文件时，可重跑原 build-db 命令，已完成且产物非空的参考步骤会跳过。
不要在现有建库进程仍运行时重复启动；部分索引须先核查，不能直接覆盖。

`virus fetch` 支持失败后留下的空目录，并校验/复用已有完整 XLSX 和来源收据。
不会覆盖非空但不完整或校验失败的目录；这类目录保留供检查。
已有病毒映射由 `db-audit` 核验，损坏或不兼容时停止，不自动覆盖冻结版本。

## 环境与入口说明

入口通过 `f/conda_env.sh` 自动查找 Anaconda、Miniforge 或 Miniconda，并激活其 `metag` 环境，
Python 和 fastp、Kraken2 等工具均从该环境加载；激活仅影响脚本及其子进程。
可用 `METAG_PYTHON` 或 `PYTHON_BIN` 显式指定解释器。帮助页无需激活环境即可显示。
入口可从其他目录用绝对路径调用；输入、输出相对路径以调用目录为准。
不带参数时显示帮助；入口关闭 Python 字节码缓存，代码目录不会因此生成 `__pycache__/`。

依赖 `openpyxl` 解析 ICTV XLSX。分类、分流和汇总不联网。
Kraken2 2.17.1 的源码使用有序输出队列；多线程模式限定该已核对版本。
其他版本须先使用 `--threads 1`，完成顺序兼容验证后再扩展版本支持。
保留每条 read 的 ID、长度、C/U、taxid 和配对检查；不以删除检查解决软件不兼容。
fastp 多线程可能改变输出次序；分类和分流均读取同一份 QC 后文件。
fastp 1.4.0 的 failed_out 不包含所有双端均失败记录。QC 计数交叉核对报告中的主输出读数、
过滤计数与实际孤立 mate；failed FASTQ 仅作诊断，另报其条数及未写出的失败条数。原始 raw 完整保留。

## 清单、raw 与校验

`manifest` 只扫描根目录和 `raw/`，检查 ENA 文件数组、大小、缺 mate 和 `.aria2`。
ENA 返回 HTTP 200 也可能包含服务端错误。程序会检查整张表的必要列、run 编号和行完整性，
拒绝错误或截断响应并重试，校验通过后才写入 ENA 缓存；`--ena-tsv` 同样检查表格。
网络错误反复出现时可改用完整的本地 ENA 报告。清单中的 `valid=no` 表示文件缺失、
大小不符或仍有 `.aria2` 等待处理项，详见 `metadata/inventory.audit.json`。
其 `samples.tsv` 保留 sample → run 两层，`files.tsv` 记录原路径、raw 路径和 ENA MD5。
同一样本的多个 runs 不自动合并。缺 `_2` 不能把 `_1` 当单端；额外无后缀文件保留为 `ARCHIVE_SINGLE`。
单文件但元数据为 PAIRED 时需审核 merged/interleaved，不能直接当原始单端。

`organize` 默认软链接，不修改原 FASTQ；软链接依赖原文件，不能随后删除原件。
首次也可选择 `--mode hardlink` 或同文件系统的 `--mode move`，搬移后保留操作收据。
MD5 必须完整读取压缩文件；`--gzip` 另做全解压与四行 FASTQ 检查。存在/大小匹配不能代替 MD5。
不要对人工整理后的清单使用 `manifest --refresh` 来做元数据升级。

## 非破坏性元数据升级

升级先创建 `samples.before_v2.<时间>.<唯一值>.tsv`，追加缺失列，初始值 `unknown`。
逐字段变更和冲突分别记录于 `metadata_upgrade.*.changes.tsv`、`*.conflicts.tsv`。
保留人工值、额外自定义列、include/valid、sample_key 和路径。更新表不允许修改这些控制字段。
非空旧值与新值不同会进入待审表，不覆盖。更新表必须有 `metadata_evidence` 和
`metadata_review_status=reviewed/provisional`。既有值要更正时，先人工审阅冲突再明确编辑原清单。

新增列：

```text
input_nucleic_acid   DNA / RNA / total_nucleic_acid / unknown
library_method      shotgun / capture / amplicon / unknown
reverse_transcription, dna_retained, virus_enriched, is_control   yes / no / unknown
processing_protocol, metadata_evidence
metadata_review_status   reviewed / provisional / unknown
```

保留 assay_class、ENA 原字段、human_read_masking、pooling、source_matrix 等。
新建 ocean WGS 清单仅根据原 README 初填 DNA/shotgun，标记 provisional 并记录依据。
CASPER 不统一推定纯 RNA、去 DNA 或总核酸；RNA-Seq 标签不能代替 run 级实验协议。
这些字段只用于解释和分层，**不参与病毒路由**。

## 冻结病毒基因组映射

先确定实际 Kraken2 库，保留它建库时的 `taxonomy/nodes.dmp`、`names.dmp`，以及已有的
`merged.dmp`、`delnodes.dmp`。不能拿新版 taxonomy 覆盖旧库。

支持离线本地 XLSX、包含 `Species / Genome composition / GenBank accession` 的 TSV；
实际 ICTV 文件中的 `Genome` 列也支持，按工作表与表头定位，不写死列号。
只解析正式 accession，多个片段归并到同库物种；不同片段指向不同物种时记录一对多异常。

交叉表支持 `accession + taxid`、NCBI `accession.version + taxid`，或以 seqid/taxid 开头的 seqid2taxid.map（可含附加来源列）。
构建器优先检查实际存在的库内 `seqid2taxid.map`、`taxonomy/seqid2taxid.map`、`accession_taxid.tsv`；
多个候选时要求显式指定。按 VMR 需要的 accession 流式过滤本地大表。
显式版本号不自动替换为其他版本。
无交叉表时仅允许**唯一且严格匹配**同库科学名/官方同义名，方法记入审计，不做模糊匹配。
可用 `virus cache-accessions` 单独批量获取 VMR 所需 accession 的 NCBI 交叉引用。

每批最多 100 个，限速、重试、保存带哈希的原始响应，支持中断续跑。
接口出错、访问限制或返回非预期 schema 时停止；不把错误页面当映射。
`mol-type` 只描述序列记录，绝不作为病毒基因组类型。

每个新的映射目录包含：

```text
virus_genome_types.sqlite / virus_genome_types.tsv
virus_metadata.source.json / virus_metadata.build.json
virus_metadata.evidence.tsv / virus_metadata.conflicts.tsv
virus_metadata.unmatched.tsv / virus_metadata.unknown_labels.tsv
```

记录源版本、URL/本地来源、原文件 SHA256、交叉表/规则/override 哈希、代码规则版本和 taxonomy 身份。
读入时校验内容与来源；同一运行中的映射不可变，内容变化即停止。构建拒绝覆盖已有目录。
同库 taxonomy 缓存位于 `.metag_taxonomy/<内容身份>.sqlite`，加锁构建、原子发布，不热更新旧缓存。

可选 `--rules rules.tsv` 和 `--overrides overrides.tsv`：

```text
taxid  evidence_scope  genome_composition_raw  source_ref  source_version  reason  review_date
```

文件为制表符分隔；scope 为 taxon/species/subtree，review_date 为 YYYY-MM-DD。
subtree 必须有权威的全分支适用依据，不能用数据库碰巧包含的几个后代推断。
override 是明确审核结论，保留被覆盖证据；祖先 override 不自动抹掉更具体的冲突。
这些表必须填写适用于实际数据库的已审核证据，不将合成示例用于正式分析。

## 类型、分类层级与输出

- dsDNA/ssDNA/dsDNA-RT → DNA；dsRNA/ssRNA/ssRNA-RT → RNA。
- Hepadnaviridae 是 DNA；Retroviridae 是 RNA。复制中间体/前病毒的分子类型不改变该规则。
- Riboviria/Viruses 等混合高层节点默认 unknown；缺失、新标签、未解决冲突均保留 unknown。
- 菌株允许继承物种；属/科仅可继承有来源的全分支规则，目录仍 `_unresolved/实际节点`。
- 病毒宿主为人/动物/植物时，病毒仍在病毒分支。viroid/satellite 等实体仅沿实际 taxonomy，entity_type 无依据则 unknown。

```text
ROOT/raw/<sample_key>/<run>/<原始文件名>
OUT/human|animal|plant|fungi|bacteria|archaea/<sample_key>/<run>/<species>/
OUT/virus_DNA|virus_RNA|virus_unknown/<sample_key>/<run>/<species 或 _unresolved/taxid_name>/
OUT/other_eukaryote|ambiguous|technical|other|unclassified/...
OUT/reports/<run>/
OUT/summary/
```

默认 OUT=`ROOT/results_v2`。不再生成聚合 `virus/` 或重复数据软链接。
非病毒语义沿用原版：human 优先于非人动物；plant 是 Viridiplantae；无法排除人的 Metazoa 等高层结果为 ambiguous。
物种目录始终依据同库 NCBI species，不把属/科结果提升为物种。
保留 PE.R1/PE.R2、QC_U1、QC_U2、ARCHIVE_SINGLE 流；PE 以 fragment 联合分配，两端同路由。

## 数据库覆盖与执行

`db-audit` 使用实际 `kraken2-inspect` 的 direct minimizers，每个节点只计一次；不能按 VMR/taxonomy 中存在物种就声称库覆盖它。
`metag.coverage.tsv` 保留七大基础类；`metag.virus_coverage.tsv` 给出 DNA/RNA/unknown、节点、物种、比例与分母。
必须满足基础 virus direct minimizers = 三类之和。未知节点不是必须存在的覆盖条件。
默认要求 DNA 和 RNA 都有实际支持；`--allow-incomplete-virus` 为显式探索选项，run 时也须提供并记录不足。
元数据未注释与无相应参考支持分别记录。七大类不足另由 `--allow-incomplete-db` 控制，不能替代病毒类型检查。
`db-audit --base-only` 仅用于新建库的基础检查，**不能据此启动 v2 run**。

已有索引优先复用。`./metag.sh build-db --db "$DB"` 默认仅打印计划；已有完整索引时仅准备审计，不重建 hash。
新建库需要显式 `--execute`，并使用新目录；含 archaeal/bacterial/viral/human/fungal/plant/protozoal/vertebrate/invertebrate/UniVec 参考。
参考下载通过 `f/k2_compat.py` 调用固定版本的 `k2 download-*`，修正 NCBI 路径尾部斜杠、404 空文件和并行错误传播；断线重试有次数限制，续传文件核对 NCBI MD5。
每个参考类别必须产生非空 FASTA 和映射文件，才标记完成。
索引构建使用 `kraken2-build --build`，规避 k2 2.17.1 的容量上限及失败返回码问题。
构建后检查三个索引文件再标记完成。
脚本默认 hash 上限为 48 GiB；本次后台构建显式设为 24 GiB，为本机其他任务预留内存。
该上限不限制参考下载、全部磁盘或实际内存，且下采样会降低敏感性。
保留 taxonomy、建库来源与审计；不要运行会删 taxonomy 的 clean。

批量 full 接受逗号列表；单目录形式用 `--root` 指定目录。执行时先审计和完整 MD5 校验。
不要盲扫 E/F 盘，不下载新推荐的 PRJEB87273，不自动循环 14 组。

## 汇总、断点与解释

`assignments.tsv.gz` 保留原列并追加 canonical_taxid、base_group、病毒类型、状态和 evidence_taxid。
详细证据在 taxid 级表，读数量不会导致全量 ID 缓存或反复扫描 FASTQ。
GzipBins 保留有限缓冲/句柄以及合法的拼接 gzip members。

- `species_by_run.tsv` 兼容原视图，历史语义包括 `_unresolved` 行。
- `groups_by_run.tsv`：互斥组计数；`qc_by_run.tsv`：QC 与守恒。
- `virus_taxa_by_run.tsv`：所有病毒 assigned taxid，包括属/科/unknown。
- `virus_species_by_run.tsv`：仅 species_taxid>0，按 run/species/group/stream 聚合。
- `virus_genome_by_run.tsv`：每个 run 的 DNA/RNA/unknown，含零计数和未知比例。
- `virus_annotation_audit.tsv`：使用的节点、来源、继承、冲突、缺失。

这些是不同视图，不能再次相加。分节段或不同 assigned 节点不能重复计为新物种。
所有视图包含 schema 和映射身份。汇总只使用有 done.json 的完整结果与每个 run 的 input.json 元数据快照，
检查签名、报告校验和、计数，并拒绝 v1/v2 或不同映射混合。

```text
raw_reads_processed = clean_reads + qc_failed_reads
clean_reads = 所有互斥输出组 reads 之和
virus_all_reads/fragments = virus_DNA + virus_RNA + virus_unknown
PE: reads=2*fragments；单端/孤立 mate: reads=fragments
```

空分母写 NA。未知 read 比例与节点比例分开。每百万 clean reads 不是颗粒、细胞、活性或绝对浓度。
普通 DNA、RNA、总核酸、病毒富集、捕获、扩增子、对照不可无分层比较。
DNA 文库未检出 RNA 病毒不能直接称为阴性；mask 人源后的零计数也不能证明原样本无人源。
候选分类不等于 confirmed species、感染或完整活病毒；必要时用 `./metag.sh align` 做独立竞争比对与覆盖核查。

旧 virus/、旧 schema、输入/参数/脚本/映射变化均要求新的 OUT；同路径换映射内容也会失效。
缺少 done.json 的部分结果不汇总、不覆盖。完成前检查全部输出，done 最后写入。
跨类别提交不是单一文件系统事务；中断产生的文件须按 input/FAILED 收据人工核查，不盲目清理。
原 raw 不删除、不重新编码、不先去人源。

## 各数据集的约束

14 组项目编号、示例路径和实验字段见 `f/datasets.example.tsv`；批量入口从该表读取显式指定数据集的项目编号，不会自动遍历执行全部条目。
特别保留：ocean 的 DNA 原核粒径子集（含 454 SE）；WBE 禁止再自动剪接头；
Oxford 每个 run 使用 UNRESOLVED_<run>；London 捕获/扩增子/对照/合池默认待审；
CASPER 按 run 协议补充且保留已报告 human mask。其他未知实验属性保持 unknown。
F: 目录名已做本次存在性检查；E: 在交付检查时不可访问，示例路径不能当成已核实路径。

<a id="ocean-source"></a>

### Ocean 数据依据

原始附件 `SOURCE_README.md` 第 76–78 行说明：指定的 `ocean.url` 对应 Tara Oceans
原核生物粒径组分的 DNA shotgun 数据，包含双端 Illumina 和少量单端 454，
不含宏转录组，也不代表 Tara 的全部项目。这段证据已保留在本节；完整原件位于
`/mnt/d/Downloads/metag_split_20261014/SOURCE_README.md`。
程序对 PRJEB1787 的 WGS 行仅据此预填 DNA/shotgun，并保持 `provisional`，仍需按实际协议审核。

## 已完成验证与执行状态

2026-10-14 已安装本机 conda `metag` 环境（Python 3.11、Kraken2 2.17.1、fastp 1.4.0），
依赖来自 conda-forge/bioconda；environment.yml 用 nodefaults 排除全局 defaults。
本次包下载中断后重试发现本机代理 SSL 错误，安装改用直连成功；未修改系统代理或 base 环境。

已实际生成上述六个数据集的清单与 raw 软链接，清单检查无错误。
用临时的小型参考库和两个合成数据集，实际验证了建库、病毒映射、审计、prepare、pilot、full、
QC、DNA/RNA 分流、汇总、重复运行跳过及 MD5 错误拦截；双端/单端及孤立 mate 的读数守恒通过。
ICTV 官方 XLSX 下载与解析也已实测。测试文件放在临时目录，未增加 tests/。

六个数据集各选 2 个 run，共 12 个 run 的原始 FASTQ 已通过完整 MD5 校验；选择与文件状态记录在
`/mnt/f/metag/reference/pilot_runs_20261014.json`。这不等于完成分类试跑。
正式共享参考库仍在构建，真实数据的分类 pilot、全量分类及其生物学结果仍未验证。
完整参考的下载/建库在独立任务中继续；上面的分析代码块不会重复启动它。
24 GiB 限制仅作用于 hash 表，下采样会降低敏感性。

## 主要依据

- ICTV VMR：https://ictv.global/vmr；本次实际获取 `VMR_MSL41.v1.20260729.xlsx`。
- 基因组类型：https://ictv.global/report/genome
- Hepadnaviridae：https://ictv.global/report/chapter/hepadnaviridae/hepadnaviridae
- Retroviridae：https://ictv.global/report/chapter/retroviridae/retroviridae
- Kraken2 手册：https://github.com/DerrickWood/kraken2/wiki/Manual
- 2.17.1 顺序实现：https://github.com/DerrickWood/kraken2/blob/v2.17.1/src/classify.cc
- fastp 1.4.0 配对与失败输出实现：https://github.com/OpenGene/fastp/blob/v1.4.0/src/peprocessor.cpp
- ENA taxonomy API：https://ena-docs.readthedocs.io/en/latest/retrieval/programmatic-access/taxon-api.html
