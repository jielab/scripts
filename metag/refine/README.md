# METAG refine

本模块读取本版发布的 ranked 分箱，汇总人源和食源序列证据、可选的群落 AMR 标记及原始分类表中的公共卫生病毒标签。

海豹疑似标签的全量读对复核、所需软件和输出见[运行说明](../README.md)。

## 运行顺序

先完成本版 split。使用 `metag.sh run` 时，refine 已列入同一批次的依赖作业；手动针对已完成样本补充证据可使用下列 refine 命令。

```bash
cd /work/sph-huangj/scripts/metag
./metag.sh refine -h
./metag.sh refine --dat-list england,Wxiamen
# 等本轮 refine 作业全部结束
./metag.paper.sh run --dat-list england,Wxiamen --raw-check
```

上面的 refine 命令可直接利用现有计数和证据表。BLAST 和 AMR 各需要外部参考，配置后会随同一个 refine 命令运行。可先指定一个样本检查参考和运行成本。每次请求显式给出参考都会重新执行对应搜索。在相同代码版本和同一份 split 记录下，后续仅汇总时会保留已完成的序列复核。候选表、病毒标签或 AMR 未执行均可从状态表区分。

## 启用碱基复核与 AMR

BLAST 核酸复核需要本地 NCBI BLAST+（blastn、blastdbcmd）。建立自选参考库还需 makeblastdb 和 blastdb_aliastool。AMR 需要 DIAMOND 2.x。可在独立环境安装这两个程序。原有 METAG Python 解释器及环境可继续使用，只需让登录节点和计算节点的 PATH 能找到新增程序。

```bash
conda create -n metag-refine-tools -c conda-forge -c bioconda blast diamond
conda activate metag-refine-tools
blastn -version
diamond version
```

**核酸参考。** 已有本地广谱、带 subject taxid 的核酸库可直接传给 --blast-db。检查人、植物和哺乳动物的核 DNA 时，尤其要查库中有哪些核染色体组装及其近缘竞争参考。NCBI core_nt 省略部分大型真核染色体，仅用它复核全基因组短读段会丢失候选。最好使用与 split 独立的组装加近缘和背景。复用同一参考可验证另一种比对方法是否得到一致结果，两种算法的共同偏差仍可能保留。

可在有磁盘空间的计算节点，用随包脚本从本地组装清单建立有版本的 BLAST 库。manifest 是含 taxid、fasta、role 三列的 TSV，每行一个已知来源的核酸组装，支持 FASTA 或 FASTA.gz。同一物种可有多个组装。程序保留物种清单、原路径和解压后内容摘要。输出目录可位于原 lib 旁边。

```bash
python refine/build_blast_db.py --manifest source_assemblies.tsv \
  --out /work/sph-huangj/data/annot/metag/source_blast_v1
```

也可以直接从现有 lib 的代表基因组目录选择物种。下面是演示人源及两个稻种的较小面板，包含黑猩猩和细菌背景。若某个 taxid 在本地缺少组装，程序会列出缺项。猪、牛、海豹及其近缘物种的复核需要将相应组装纳入面板。出于统计解释的需要，未被覆盖的物种会明确列为无法检验。

```bash
python refine/build_blast_db.py --library-taxids 9606,9598,4530,4538,1423 \
  --out /work/sph-huangj/data/annot/metag/human_rice_blast_v1
```

**AMR 参考。** 使用已验证的 ShortBRED 抗生素耐药蛋白 marker FASTA。软件只使用长度至少 25 aa 的 TM 标记。Family、ID_Orig、Merged.ID 三列的可选注释可把部分参考 accession 翻译为基因或功能家族。下面给出一个可复现的公开参考版本。根据建库文档中的读长选相应 marker 文件。150bp 文件仅是例子。该历史面板覆盖有限，不能当作所有当前耐药机制的清单。

```bash
git clone https://github.com/biobakery/shortbred.git /work/sph-huangj/data/annot/metag/shortbred_markers
git -C /work/sph-huangj/data/annot/metag/shortbred_markers checkout ef6d8b5b8489ccb38cb08434ae965fe92548b3ec
marker_dir=/work/sph-huangj/data/annot/metag/shortbred_markers/sample_markers
export METAG_AMR_MARKERS="$marker_dir/ShortBRED_ABR_150bp_markers.faa"
export METAG_AMR_ANNOTATION="$marker_dir/metadata_and_consensensus_sequences/ShortBRED_ABR_Metadata.tab"
export METAG_BLAST_DB=/work/sph-huangj/data/annot/metag/source_blast_v1/competition

./metag.sh refine --sample Wxiamen/SRR39271683
# 单样本检查成功后
./metag.sh refine --dat-list england,Wxiamen
# 等 refine 完成
./metag.paper.sh run --dat-list england,Wxiamen --raw-check
```

若建的是 human_rice_blast_v1，请将 METAG_BLAST_DB 指向该目录内的 competition。参考参数也可直接通过 --blast-db、--amr-markers、--amr-annotation 传入。未提供某项参考时该项记为 not_run。默认纯表格汇总每队列一个作业，启用搜索后每个作业按顺序处理同队列最多四个样本，可用 --batch-size 调整。项目共享调度器支持 --plan-only 和 --submit local。32 GiB 为 refine 作业默认预留；与 split 的大内存预加载是不同的任务。

## 样本元数据与输出

元数据是制表符分隔的文本。sample 为 DATASET/SAMPLE。sample_type 可用 stool、wastewater、environment、blank、mock、other、unknown。library_type 可用 DNA、RNA、DNA_RNA、unknown。DNA_RNA 指已有记录证明实验覆盖两者。传 --metadata 时应包括该命令选中的所有样本。单独的 unknown 不含“阴性”的含义。未提供来源记录的新样本默认为 unknown。此前同一 split 已登记的样本元数据在无参数重汇总时沿用。

analysis/metag/refine/main/DATASET/SAMPLE/ 保存逐样本结果，jobs/ 和 logs/ 保存提交文件和日志。临时查询只在本地使用。完成后再运行 metag.paper.sh，报告会在 analysis/metag/paper/main/ranked/ 中自动汇入选中样本的结果，并比对每个样本的 split receipt。split 已更新的样本标为 split_changed，需重新 refine。不同样本类型、建库方案的来源统计分别汇总。

| 文件 | 可回答的问题 |
|---|---|
| refine_status.tsv | 哪些样本及可选分析已完成，使用什么参考和工具，AMR 测了多少家族、多少对无法区分家族 |
| refine_sources.tsv、refine_sources_by_dataset.tsv | 每个样本和各队列分到哪些人源、植物和哺乳动物分类单元，哪些与维护的小型食源词表匹配；鸟类仅含词表内的分类单元 |
| refine_blast.tsv、refine_blast_by_bin.tsv | 各 bin 随机抽了多少对，其中多少支持已报告层级、更高层级、其他谱系，或缺乏可判别证据 |
| refine_amr.tsv | 所给 marker 集合中，每个可测家族获多少读对支持，每百万原始输入读对的支持数。没有相应注释时保留参考的家族 ID |
| refine_virus.tsv | 既有种级分类表中是否出现 HBV 相关物种、HIV-1/2 相关物种或流感 A/B/C/D 标签。当前与历史物种名按精确名称匹配 |
| target_panel_qc.tsv、target_panel_by_dataset.tsv | 原 split 中人及已检出植物的面板机会、窗口门槛和检出分母 |

## 方法和证据边界

**来源。** 食源词表按原始 rank 和学名标记来源可能性，保留未列入词表的全部植物和哺乳动物。牛源 DNA 的可能来源包括肉、乳、加工原料和环境。鸡属于鸟类。稻属和具体稻种是互斥的原始分箱，汇总到属时应分别加各 bin 一次。粪便的食源 DNA 研究可用已知膳食及阴性对照校准。混合环境材料的群体信号无法与某位贡献者对应。人源结果仅报告 Homo sapiens 层级，当前流程没有估计个体身份、社会种族或遗传祖源。

**BLAST。** 各个非空 human、plant、mammal 种级 bin，词表内食源种及食源属级 bin 默认固定种子、无放回均匀抽取最多 100 对。读取被抽查的 bin 时核对整个 bin 的内容摘要和读对数。每端在本地搜索。保留单段覆盖至少 80%、identity 至少 95%、bitscore 至少为合格最佳分数 98% 的匹配。先求单端共同祖先，再对两端求共同祖先。缺失 taxid、任一端缺乏支持或达到每端 100 个候选的上限时单列。两个端的共同谱系不等同于对同一基因组位置的成对比对。以上启发式阈值用于复核，在独立 mock 上校准后才能当作准确率。抽取阳性 bin 能评价条件一致性，估计漏检还需对 shared、unclassified 及其他 bin 有代表性地抽样。

用随包工具生成的数据库有明确物种清单。原始种不在库中、缺少近缘对手、缺少该分支外背景均记为 panel_not_evaluable，实际抽查分母为零。植物种级要求同属对手，哺乳动物和鸟类种级要求同科，属级复核要求属内至少两个参考种。对其他广谱库，reference_scope 会写 external_inventory_unverified。即使库中有一个同科物种，也可能遗漏真正的最近缘种。相同 DNA 在不同建库策略中的检出率不能直接换算成摄入量。

**AMR。** 从原始两端核对 split 的读对总数和完整 FASTQ 摘要，再以 DIAMOND blastx 搜索全部 marker。按单端保留 identity ≥95%、比对至少 25 aa、长度覆盖读段可翻译长度和标记长度较小者的至少 95% 的匹配。取合格命中中距最佳 bitscore 2% 以内的全部家族。两端证据合并，一对最多给一个家族计数，跨家族歧义单列在状态表。一个有清晰支持的端也可计数。both_mates_pairs 给出两端均支持该家族的读对数。uniquely_supported_markers 给出至少被一端唯一指向的不同 marker 数，避免同一家族中多个并列 marker 虚增覆盖证据。输出以所有原始输入读对为分母，混合样本中宿主或其他组分的份额会影响该指标。该计数器使用公开 marker 思路并采用自己的双端判据，不是 ShortBRED 软件的 RPKM。它不能判定某株细菌携带该基因，也不覆盖依赖特定位点突变或表达调控的全部机制。菌株归属需要独立的物理连锁证据，表型需要药敏证据。

**病毒。** 公共卫生表是既有种级标签的回顾性核对，阳性显示 classifier_signal_requires_confirmation。零表示在这批 split 的种级清单中无匹配名称。普通 DNA 建库对游离 HIV RNA 和流感 RNA 缺乏检测能力，HIV 前病毒 DNA 是不同的检测对象。HBV 相关物种级名称也不足以确定人来源。公共卫生结论另需与样本材料相符、经过验证的检测和对照。

**验证设计。** 将提取空白、建库空白、已知单一物种和近缘混合材料与研究样本一起分析。先在这些材料上评估 species 和 genus 两个层级，再比较研究样本。对有逐对独立真值的样本，可在 paper 中提供 --truth 计算 precision、recall、F1。盲测参考基因组尽量与建库组装分开。食物摄入问题还需要膳食记录或控制喂养结果；技术一致性无法代替这种验证。

## 文献与软件

- Diener 等。Metagenomic estimation of dietary intake from human stool，Nature Metabolism，2025。[研究](https://doi.org/10.1038/s42255-025-01220-1) · [MEDI](https://github.com/Gibbons-Lab/medi)
- Tomofuji 等。Reconstruction of the personal information from human genome reads in gut metagenome sequencing data，Nature Microbiology，2023。[研究](https://doi.org/10.1038/s41564-023-01381-3)
- Kaminski 等。High-Specificity Targeted Functional Profiling in Microbial Communities with ShortBRED，2015。[研究](https://doi.org/10.1371/journal.pcbi.1004557) · [标记文件及注释](https://github.com/biobakery/shortbred/tree/ef6d8b5b8489ccb38cb08434ae965fe92548b3ec/sample_markers)
- Buchfink 等。Sensitive protein alignments at tree-of-life scale using DIAMOND，Nature Methods，2021。[研究](https://doi.org/10.1038/s41592-021-01101-x)
- 关于短读长 AMR 与表型的解释。[Boolchandani 等，Nature Reviews Genetics](https://doi.org/10.1038/s41576-019-0108-4)
- 关于核酸库覆盖和污水数据解释。[NCBI BLAST 核酸库说明](https://blast.ncbi.nlm.nih.gov/doc/blast-help/blastdatabases.html) · [CDC 污水监测说明](https://www.cdc.gov/nwss/about-data.html)
