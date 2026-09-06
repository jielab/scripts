# GU：现代人古人类基因片段分析流程

GU 用统一入口运行四类分析，并把可比较的结果标准化到 SQLite/TSV，供 R Shiny 查询和可视化。当前 PhyML 分层流程版本为 `2026-09-04.4`。

已安装 `2026-09-05.2` 复核扩展：原生 Shiny 的 PhyML 页面支持输入 lead SNP 与最佳单倍型 tag SNP 对照，并新增 Evidence review 页面。先运行下列只读后处理生成展示数据，再重启 Shiny：

```bash
python3 /mnt/d/scripts/gu/f/prepare_review.py \
  --normalize /mnt/d/analysis/gu/normalize \
  --database /mnt/d/analysis/gu/gu.sqlite \
  --output /mnt/d/analysis/gu/review
./gu.sh shiny
```

运行实现直接位于 `f/`：ARG 为 `arg.sh` / `arg.py`，共用逻辑为 `comm.sh` / `comm.py`，Shiny 为 `shiny.R`；PhyML 的大型阶段保留为平铺的 `phyml_*.py`。AS3 第三方源码位于 `f/as3/`。复核入口为 `prepare_review.py`，调用 `dual_lead.py` 并使用 `review_template.html`。`Rscript f/shiny.R --review --data /mnt/d/analysis/gu/review` 可仅显示独立复核页面。测试源码已移除，不需要旧安装补丁或归档。`/mnt/d/analysis/gu/review/` 是结果目录，与源码位置无关。

这次结构整合保留科学计算参数；ARG 的源码指纹现引用合并后的文件，旧准备缓存可能需要重新验证或生成。ARG/PhyML 的对外命令仍为根目录 `arg.sh` / `gu.sh`。

| 模块 | 适合的分析范围 | archaic reference VCF | 关键限制 |
|---|---|---:|---|
| `phyml` | 候选 loci；也可按固定窗口扫描整条染色体 | 需要 | region-first、anchor-optional；先保存完整区域相似性，再做 derived/outgroup/复现性与树过滤；不是正式的 genome-wide tract caller |
| `ibdmix` | 全染色体/全基因组；loci 可作为有侧翼的局部运行 | 需要 | 常染色体是原生路径；男性 non-PAR X 需 pseudo-diploid 接口适配，结果标记 experimental |
| `trace` | 已有整染色体 ARG 后做全染色体联合推断；loci 在结果端裁剪 | 不需要 | UKB 单倍型必须进入同一套 ARG；不能像 PCA 一样投影到只含 1KG 的 ARG |
| `as3` | 官方支持的 GRCh38 常染色体 chr1–22 | 需要 | loci 模式仍按整条染色体预处理/推断，最后裁到 loci；GRCh37 和 chrX 都在任何大文件预处理前停止 |

Zeberg 与 Pääbo 2020 Nature 文章没有发布一个以作者命名的 caller。其 Methods 是构建 1000 Genomes 单倍型，并以 PhyML 3.3、HKY85 模型和 bootstrap 重建系统树。因此本流程使用 `phyml` 作为唯一模块名。由于完整 1KG 数据在一个 locus 内也可能产生数千种不同单倍型，GU 默认不运行包含全部 unique H types 的超大探索树；`--plot-phy TRUE` 只运行常见单倍型区域树和过滤后的 candidate-focused trees，除非显式设置 `PHYML_RUN_ALL_UNIQUE_TREE=TRUE`。

## 可复制的运行顺序

TRACE 需要先有整染色体 ARG。ARG 是可复用的数据准备产物，统一由
`gu/arg.sh build` 构建。用户只需在 ARG 生成后运行 `gu.sh trace`；
TRACE 会在启动分析前自动完成只读 ARG 检查，不再提供单独的 `gu.sh arg`
用户模块：

```bash
./gu.sh phyml --loci /mnt/d/files/gu.37.bed --plot-phy TRUE --jobs 1

./gu.sh ibdmix --chr X --grch 37 --target 1kg --target-dir /mnt/d/data.BIG/refGen/1kg/37/pfile/chr --jobs 1

bash /mnt/d/scripts/gu/arg.sh build --method tsinfer --dir-gen /mnt/d/data.BIG/refGen/1kg/37 --chr 22,X

./gu.sh trace --chr 22,X --grch 37 --target 1kg --target-dir /mnt/d/data.BIG/refGen/1kg/37/pfile/chr --jobs 1

bash /mnt/d/scripts/gu/arg.sh build --method needle --format trace --dir-gen /mnt/h/ukbGen/37 --dir-pfile /mnt/h/ukbGen/37/hap --chr 22

./gu.sh trace --chr 22 --grch 37 --target ukb --target-dir /mnt/h/ukbGen/37/hap/chr --arg-dir /mnt/h/ukbGen/37/arg/trace/needle --jobs 1

./gu.sh as3 --chr 3,22 --grch 38 --target 1kg --target-dir /mnt/d/data.BIG/refGen/1kg/38/pfile/chr

./gu.sh final

./gu.sh shiny
```

AS3 仅支持 GRCh38 chr1–22；指定 chrX 会在预处理前报错。

## 命令文件与并行调度

`phyml`、`ibdmix`、`trace` 和 `as3` 使用同一个外层调度规则。每次请求先生成永久、可独立重跑的 `.cmd` 文件和 `<request>.cmd.list`；命令只引用正式 target prefix、永久结果目录与必要的永久 BED，不引用 `run.*`、临时 VCF 或临时 BED。TRACE 在进入该调度器前通过内部 `f/arg.sh trace_check` 调用 `arg.sh check`，只读验证 ARG，不写入或转换 ARG 数据。

重新提交 IBDmix 请求时，若某个染色体目录同时存在 `.complete` 和对应阈值的非空 `final/all_archaic_refs.*.segments.tsv.gz`，且未指定 `--replace-ibdmix TRUE`，外层调度器会直接报告 `[GU CMD] SKIP`。该判断发生在 target 预检和临时 PFILE→VCF 转换之前；`.cmd.list` 仍保留全部请求单元，便于审计和独立重跑。

- 默认 `--run-cmd TRUE --foreground FALSE`：请求调度器通过 `setsid`/`nohup` 在后台运行，入口立即报告 launcher 的日志、PID 和状态文件。`--jobs N` 最多并行 N 个 `.cmd`（默认 N=4）。AS3 是例外：当前每个独立 chromosome 命令使用同一组 GPU assignment，本地调度会把 effective jobs 固定为 1，避免多个染色体同时争用 GPU 0；`--jobs` 仍影响其他方法。每个任务的 stdout/stderr 写入同目录 `<unit>.log`，失败时另写 `<unit>.err`。
- 后台提交成功后会按 `gwas_post.sh` 风格打印可直接复制的 `Progress: tail -f ...`、永久 `Command list`、`Job: pgrep -af 'gu.sh METHOD'` 和 `Kill PID ...: kill -- -PID`，并显示 PID/status 文件。负 PID 表示按进程组停止整个请求及其 worker；只执行 `kill PID` 不保证已经启动的后代进程一并退出。
- 指定 `--foreground TRUE` 可让整个请求保持连接终端并等待完成，适合交互调试或由外部作业调度器管理。永久的单元 `.cmd` 是 leaf worker，直接执行时始终以前台方式运行，不会再次把自己送到后台。
- 显式指定 `--run-cmd FALSE` 时只生成 `.cmd`，不启动分析。
- 以后在 HPC 上可将 `.cmd.list` 中的命令逐个交给 `bsub`；当前版本不猜测集群参数，也不内置 `bsub`。

不同模块的后台请求彼此独立，可以同时提交，例如同时运行 `phyml --foreground FALSE` 与 `ibdmix --foreground FALSE`。launcher 日志/PID/状态文件按 method、scope、时间和启动进程编号唯一命名；每个请求的 target 转换也使用独立的 `run.*` 临时目录。`--jobs` 是单个请求内部的并行度，因此同时提交多个模块时应按总 CPU、内存和 GPU 容量分配各自的值。默认 `--auto-final FALSE` 可避免多个模块完成时同时刷新同一个 SQLite；需要自动刷新时建议只给其中一个请求设置 `--auto-final TRUE`。

调度单元遵循各方法的原生分析边界：PhyML/IBDmix 的 `--loci` 每行一个命令；TRACE/AS3 把同一染色体的 loci 聚合成一个命令，以保留整染色体模型上下文；所有 `--chr` 请求都是每条染色体一个命令。因此 `--chr 22,X` 会生成两个完全独立的任务；全基因组 IBDmix/TRACE 生成 23 个命令，AS3 生成 22 个。`--jobs` 只控制任务间并行；PhyML 的 `--phyml-jobs` 仍控制单条染色体任务内的窗口并行，两者不要同时设得过高。外层只改变任务组织方式，不改变四个模块内部的科学计算。

## 输出、缓存与分享目录

`/mnt/d/analysis/gu` 是可整体删除并重跑的分析目录；顶层文件夹与 GU 模块对应。最大的 SQLite 数据库直接放在根目录；`normalize/` 只放跨方法 summary、各模块标准化表和 Shiny 真正读取的轻量附件。数据库内的附件路径相对分析根目录保存。

```text
/mnt/d/analysis/gu/
  gu.sqlite             # Shiny 主数据库（大文件）
  as3/                 # AS3 本次分析结果
    [target]/          # 例如 1kg/ 或 ukb/
      chr22/
      log/
  ibdmix/              # 同样使用 [target]/<analysis-unit> 与 [target]/log
  phyml/               # 同样使用 [target]/<analysis-unit> 与 [target]/log
  trace/               # 同样使用 [target]/<analysis-unit> 与 [target]/log
  normalize/           # 轻量分享目录；启动 Shiny 时与 gu.sqlite 一起提供
    summary/            # 跨方法标准化 TSV
    as3/               # 模块 segments.tsv.gz
    ibdmix/            # 模块 segments.tsv.gz
    trace/              # 模块 segments.tsv.gz
    phyml/              # loci、树和 Shiny 序列文件
```

`./gu.sh final` 会刷新根目录的 `gu.sqlite` 和轻量 `normalize/` 分享内容；大型输入、参考面板和日志不会被复制进去。运行 `./gu.sh shiny` 默认读取根目录 `gu.sqlite`。若对方需要实际启动 Shiny，应把 `gu.sqlite` 与 `normalize/` 一起提供并保持二者相邻；只分享表格 summary 时，发送 `normalize/` 即可。请求 BED 的规范化中间文件在系统临时目录创建并在退出时删除；稳定的 request BED/map 会随各模块结果保存。PFILE 是首选输入；只在下游程序必须读取 VCF 时，GU 才在 `/mnt/d/tmp/gu/<target>/<build>/run.*/chrN/` 生成单次运行文件，例如 `1kg/37/run.*/chrX/chrX.vcf.gz`。混合请求按 `chr22/`、`chrX/` 分开，不生成 `chr22_X` 目录；`chrX` 默认就是 male-only non-PAR。已有原生 VCF 时只建立软链接，不复制数据；`run.*/vcf/chr*.vcf.gz` 也只是供下游读取的统一软链接接口。ARG 的大型 VCF-Zarr/LMDB 工作文件位于 WSL ext4 的 `/tmp/gu/arg/<target>/<build>/run.*/chrN/`，不放到 `/mnt/d`；任务退出时删除。`analysis/gu` 只放正式结果，source reference tree 不保存 target conversion 或 AS3 preprocessed 文件。`normalize/reference-callsets.sqlite` 是跨不同 chromosome/locus 的多次 `final` 与 Shiny 查询共同使用的官方 published-callset 标准化数据库。

## BED、坐标与 build

`--loci` 文件为 BED：`chrom  start0  end  locus_id`，0-based、half-open。`--loci` 与 `--chr` 互斥。第四列对 PhyML 是可选的 exact-marker 标签，只用于完整区域扫描之后的二级 zoom；缺失、写错或在 VCF 中找不到都不会删除区域候选。SNV、MNV 和 `rs8176719` 一类 indel 都按确切 VCF record 单独记录，不用邻近 SNP 代替。GRCh37/38 分别使用 `/mnt/d/files/gu.37.bed` 与 `/mnt/d/files/gu.38.bed`；PhyML 未给 `--loci/--chr` 时按 `--grch` 自动选择其中一个，不会改写源 BED。当前默认文件分别包含 8 个 GRCh37 loci 和 5 个 GRCh38 loci，ABO 与 APOE 是 negative controls。

`--loci-flank` 默认 `100kb`，支持 `bp/kb/mb`。GU 保存原始 `*.core.bed`、实际计算的 `*.bed` 和同时列出 core/analysis 边界的 `*.map.tsv`。侧翼会在染色体端点截断，并进入请求签名，因此不同 flank 不会误用同一结果。

`--target` 是稳定、适合做文件夹名的数据集名称（例如 `1kg`、`ukb`），`--target-dir` 是按染色体拼接编号的 VCF/pfile 前缀（例如 `/mnt/h/ukbGen/37/hap/chr`）。两者必须一起提供；都省略时默认使用当前 build 的 1KG pfile。输出统一写入 `<module>/<target>/`，因此同一模块的 1KG 与 UKB 结果不会互相覆盖。

标准 chrX 路径固定使用无 PAR 的 `<PREFIX>X.male`：1KG reference 和 target 中每名男性都只有一个明确 phased haplotype。代码仍保留 PAR/`<PREFIX>XY` 兼容路径，但它不是默认 chrX 分析的一部分。chrX 会先按实际 analysis interval（已包含 `--loci-flank`）分类；`--chr X` 以及完全位于两个 PAR 之外的 X loci 使用 `<PREFIX>X.male`。跨越 PAR/non-PAR 边界会拒绝运行。PLINK2 pfile 已经足够，`.pvar.zst` 会自动以 `--pfile PREFIX vzs` 打开；`chrX.male` 会校验 non-PAR、male-only、haploid。

也支持预先生成的 `<PREFIX>X.male.vcf.gz`（加 `.tbi/.csi`），但不是必需步骤。原生 VCF 路径会额外校验所有 records 都在 X、两个 PAR 区域没有记录、所有非缺失 GT 都是单倍体，并要求 `--sample-panel`/邻近 `samples.txt` 中的男性 ID 与 VCF samples 完全一致。`arg --target-vcf-dir` 的内部 `chrX.vcf.gz` 也执行同一校验。无论输入是 pfile 还是 VCF，需要 VCF 的下游只看到单次运行临时接口 `chrX.vcf.gz`。不要用含 PAR、女性或二倍体 chrX 的普通 VCF 代替。

输出布局服从各方法的模型边界。PhyML/IBDmix 的染色体或 BED 行是独立分析单元；PhyML 的单-locus 结果直接写在运行目录的 `loci/` 下，整染色体分窗结果写在 `loci/<window_id>/` 下。AS3 把同一染色体的所有 loci 聚合为一次整染色体推断；TRACE 的单个 `.cmd` 也限定为一条染色体，避免调度上把多染色体绑成一个作业。`./gu.sh final` 再递归合并为 Shiny 使用的 SQLite/TSV。每个结果保留自己的 BED、map、日志、provenance 和 build；预检/运行日志位于对应模块的 `<target>/log/` 下。IBDmix 对 pfile 输入记录原始 `.pgen/.pvar/.psam` 转换签名，而不把可删除临时 VCF 的 mtime 当作科学输入；完成标记、provenance 和全部压缩产物通过校验时，重复命令只报告已完成并跳过计算。

```text
3   45859651   45909024   rs35044562
12  113350796  113425679  rs10735079
```

每个输出保存自己的 `genome_build`。标准化数据库按 build、来源类别、染色体和边界分别建目录，GRCh37 与 GRCh38 不会再被聚成同一个 `segment_code`。

若 PLINK2 pvar 使用 `chr:pos:ref:alt` 一类 coordinate ID、没有 gold-standard rsID，`check_GRCH` 的 rsID 检查失败后会自动进入严格的坐标 sentinel fallback。本机真实 pvar 集成测试得到 GRCh38 `0/4`、GRCh37 `5/1` 个 build-specific sentinel，分别正确推断为 GRCh38/GRCh37。

## PhyML locus workflow

分析顺序为区域发现 → 精确 marker 注释 → diagnostic/outgroup 筛选 → 候选树 → GWAS 风险核对。各 `PHYML_MIN_*` 等环境变量优先于 profile。

默认候选筛选要求 YRI 频率 ≤0.05、单倍型拷贝数 `n>10`、至少一个匹配 diagnostic 位点，并对**实际 diagnostic marker span**要求 ILS 概率 `<0.1`；不额外限制 diagnostic 匹配比例或相邻位点最大距离。ILS 参数为 0.53 cM/Mb、分化 550000 年、古样本 50000 年、29 年/代，用 Gamma(shape=2) survival 计算。`PHYML_MAX_ILS_PROBABILITY=1` 可关闭此门槛。树使用 HKY85 和 100 次 bootstrap。

区域发现使用完整区域 SNP 集、archaic-derived 定义与 male-only X。GWAS risk allele 不参与候选发现。树确认要求 bootstrap ≥70，更高置信度还需满足更严格的 bootstrap/纯度条件；高相似度本身不等于渐渗。

`final/evidence_parameters.json` 记录 profile、实际门槛与 ILS 模型；`candidate_ils_probability` 是模型筛选量，不是渐渗概率。`f/phyml_markers.tsv` 保存 build/坐标/REF/ALT/风险等位基因。anchor 查找在 rsID 缺失时只使用该表的**精确**四元组，不猜最近 SNP。可用 `PHYML_MARKER_MAP` 指定自己的表。`final/gwas_risk_links.tsv` 在发现和树分析后统计候选携带风险等位基因的比例；`named_anchor_matches_lineage` 则表示是否等于古参考的 anchor allele，两者均不代表 GWAS 共定位或因果证明。

`gu.sh final` 自动运行 `f/phyml_report.py`，从每个谱系的 `evidence_haplotypes.tsv` 和 `evidence_trees.tsv` 构建多行报告。Shiny **PhyML evidence** 展示全部位点的逐谱系汇总、逐单倍型明细、各谱系候选树以及 IBDmix/TRACE 验证。不会用一个 locus 的主谱系代替另一个谱系的结论。无候选和无 diagnostic 位点也保留行。

每个候选均有指定 lead SNP、最佳 tag SNP、等位基因和 SNP–haplotype r²。单倍型选择不依赖 GWAS；最佳 tag 在保存的 core 双等位 SNP 中寻找，要求该 SNP 的目标等位基因与最匹配古参考一致，按 r²、F1、call rate 和诊断位点/区间位置排序。汇总行的 tag 标记整个谱系候选集合，明细行的 tag 标记单个 H 类型；两者的 r² 不可混用。单个位点不能良好标记目标时保留低 r² 或缺失值，不把 tag SNP 当作 GWAS 因果 SNP。另保留全区相似度最高的比较单倍型及指定 marker 古参考等位基因的常见比较单倍型，其角色不构成阳性判定。

验证同时报告同一个体的任意片段重叠，以及单条片段覆盖候选实际区间 ≥80% 的人数。IBDmix 要求同谱系；TRACE 是 Ghost/Unknown，只提供位置与染色体拷贝编号的支持。IBDmix 无相位信息，chrX 仅作探索。未运行、样本/区间未被分析、有部分重叠、已运行但无重叠分别标记。分母是参与比较的候选携带者，不能解释成全基因组可检出率。结果写入数据库旁的 `review/phyml_locus_report.tsv`、`phyml_haplotype_report.tsv`、`phyml_copy_validation.tsv`、`phyml_lineage_trees.tsv`，可在 Shiny 下载。完整运行顺序见 [八个位点运行说明](RUN_EIGHT_LOCI.md)。

```bash
./gu.sh phyml check --loci /mnt/d/files/gu.37.bed --grch 37
./gu.sh phyml --loci /mnt/d/files/gu.37.bed --plot-phy TRUE --jobs 1
./gu.sh phyml --chr 3,9,12,19,X --grch 37 --run-cmd FALSE
```

### 分层顺序

对 `--loci`，BED 的完整 core 永远先于 anchor 分析。`--phyml-region-mode ld` 和 `--phyml-ld-r2` 为旧命令兼容及请求 provenance 保留，但不会再用 anchor-SNP LD block 缩小 Stage 1。

1. **区域发现（anchor-independent）**：保存每个 core 内全部现代 H types、全部可用古人类参考投影、完整 H × reference pairwise matrix、局部高相似片段和按 overlap 聚类的 regional families。Stage 1 不使用 BED 第 4 列、祖先等位基因、YRI/AFR 频率、复现性或 bootstrap。
2. **常见单倍型区域树（anchor-independent）**：`--plot-phy TRUE` 时按旧 workflow 的 `n > 10` 规则（默认 `PHYML_COMMON_TREE_MIN_COPIES=11`）建立 all-archaic、Neanderthal 和 Denisovan common-H trees。
3. **精确 anchor zoom（可选）**：读取 BED 第 4 列指定的 exact VCF record，把 allele 分配到 sample-haplotype copies 和 SNP-defined H types。未提供、未找到或 phase 无法判定时写出明确的 `not_evaluable` 状态，不把它解释成阴性。
4. **derived/outgroup filter**：要求古人类谱系 consensus allele 相对 `INFO/AA` 为 derived，并在 YRI（不可用时 AFR）中低频；同时要求 diagnostic markers 形成连续 cluster。
5. **复现性过滤与 focused tree**：兼容 profile 要求至少 11 个总 copies、11 个 non-outgroup copies，并通过上述 ILS 筛选；strict profile 保留此前 2 copies、3 diagnostic matches 的门槛。低频 H types 仍永久保留在 Stage 1，并可作为 `candidate_context`，但不能单独产生高置信度 locus call。
6. **独立方法验证**：PhyML 是 haplotype candidate layer；正式生物学结论应再与同一个体的 IBDmix、TRACE 或 AS3 tract overlap 结合。

`introgression_call` 报告区域 tract，`anchor_linked_introgression_call` 单独报告 named marker 与该 tract 的连接。anchor 与区域 candidate 不一致时，两者不会被强行合并。

### 关键输出

所有分层结果位于每个运行目录的 `final/`：

| 层 | 主要文件 | 含义 |
|---|---|---|
| Stage 1 区域清单 | `region_haplotypes.unfiltered.tsv.gz`、`region_archaic_references.unfiltered.tsv.gz`、`region_matches.unfiltered.tsv.gz` | 全部现代 H types、实际 chromosome-copy 分母、全部古参考以及所有 H × reference 比较 |
| Stage 1 候选与片段 | `region_haplotype_candidates.unfiltered.tsv.gz`、`region_match_segments.unfiltered.tsv.gz`、`region_candidate_families.unfiltered.tsv.gz` | 宽松相似性候选、观察到的 matching-marker span 和重叠 family；均不是 introgression call |
| 常见 H trees | `region_common_tree_inputs.tsv`、`region_common_trees.tsv`、`region_tree_clades.unfiltered.tsv.gz` | `n >= 11` 的区域树输入、状态和全部满足基本结构的 archaic-modern clades |
| exact anchor | `anchors.tsv`、`anchor_copies.tsv`、`anchor_haplotypes.tsv`、`anchor_allele_groups.tsv` | exact marker 的 REF/ALT、古参考 allele、现代 copies、H groups 和 common-H 归属 |
| filtered evidence | `evidence_sites.tsv`、`evidence_haplotypes.tsv` | 完整 fate tables，包括通过和未通过原因 |
| pass-only views | `diagnostic_sites.filtered.tsv`、`diagnostic_candidates.filtered.tsv`、`anchor_zoom_candidates.tsv` | 通过 derived/outgroup/连续性/复现性规则的便捷视图 |
| focused trees / 结论 | `evidence_trees.tsv`、`evidence_loci.tsv` | bootstrap、candidate purity、candidate sensitivity、区域 call 与 anchor-linked call |
| provenance | `archaic_references.selected.tsv`、`region_scan_parameters.json`、`evidence_parameters.json` | 实际参考面板、阈值和运行参数 |

所有 `*.unfiltered.tsv.gz` 清单直接以 gzip 格式写出，无需手动压缩；证据汇总、dual-lead 和 normalize 打包均支持读取压缩文件，并兼容旧的 `*.unfiltered.tsv`。normalize 打包旧 TSV 时也会自动压缩。新文件写入完成后才替换旧结果并移除同名未压缩清单。

`region_scan_loci.unfiltered.tsv.gz` 会为失败、skip 或缺少 sequence artifact 的 locus 保留 `stage1_status=not_evaluable:<raw status>`；不能把“未完成/不可评价”当成“未检测到 introgression”。`direct_match_pass` 只表示原始 A/C/G/T identity，保留用于兼容，不是 carrier 判定或 introgression call。诸如 `798/847` 的值表示 847 个双方可调用 alignment columns 中有 798 个碱基相同，不表示 carriers 数量。

### 默认阈值与树行为

```bash
# Stage 1 local similarity
export PHYML_REGION_SCAN_MIN_SITES=12
export PHYML_REGION_SCAN_MIN_PROP=0.90
export PHYML_REGION_SCAN_MAX_GAP_BP=50000
export PHYML_REGION_FULL_MIN_SITES=10
export PHYML_REGION_FULL_MIN_PROP=0.80
export PHYML_REGION_FAMILY_MIN_OVERLAP=0.50

# Common-H trees
export PHYML_COMMON_TREE_MIN_COPIES=11
export PHYML_COMMON_TREE_MIN_SITES=10
export PHYML_REGION_TREE_BOOTSTRAP_MIN=70

# Optional strict-profile overrides (omit these to use legacy_compatible defaults)
export PHYML_THRESHOLD_PROFILE=strict
export PHYML_OUTGROUP=YRI
export PHYML_OUTGROUP_FALLBACK=AFR
export PHYML_MAX_OUTGROUP_ALLELE_FREQ=0.05
export PHYML_MAX_OUTGROUP_HAPLOTYPE_FREQ=0.05
export PHYML_MIN_OUTGROUP_COPIES=20
export PHYML_MIN_DIAGNOSTIC_SITES=3
export PHYML_MIN_DIAGNOSTIC_MATCH_PROP=0.80
export PHYML_MIN_CANDIDATE_COPIES=2
export PHYML_MAX_DIAGNOSTIC_GAP_BP=50000

# Candidate-focused tree
export PHYML_EVIDENCE_BOOTSTRAP_MIN=70
export PHYML_EVIDENCE_BOOTSTRAP_STRONG=90
export PHYML_EVIDENCE_MIN_PURITY=0.80
export PHYML_EVIDENCE_MIN_SENSITIVITY=0.50
export PHYML_EVIDENCE_CONTROLS=24
export PHYML_EVIDENCE_TREE_FLANK_BP=5000
```

不要为让特定 locus “通过”而反复调阈值；修改阈值时应预先声明，并用固定 sensitivity grid 报告。PhyML bootstrap 是重采样 alignment columns 后某拓扑边的出现比例，不是“该 locus 有多少概率来自 introgression”。

`--plot-phy FALSE` 仍完成区域清单、anchor zoom、evidence fate tables 和显式的 `not_requested` tree 状态；`--plot-phy TRUE` 才运行分层树与 QC PNG。全 unique-H tree 默认关闭，可显式设置 `PHYML_RUN_ALL_UNIQUE_TREE=TRUE`。单棵树失败时已完成比较和其他树会保留；默认 `PHYML_TREE_FAILURE_POLICY=warn`，需要把任一树失败视为作业失败时改为 `error`。

默认自动寻找本地完整的 `Altai Chagyr Vindija Denisova Denisova25` 面板，并把实际选择写入 `archaic_references.selected.tsv`。可显式固定面板：

```bash
export PHYML_REFS='Altai Chagyr Vindija Denisova Denisova25'
```

比较缓存使用 `final/comparison.input.meta.tsv`；原始输入签名一致时可复用旧 raw comparison，但新版 region scan、anchor zoom、outgroup filters 和 layered summaries 每次重建。需要强制重做 raw comparison 时使用 `--replace-phyml TRUE`。

`--chr X` 按 `--phyml-window-bp`（默认 500 kb）在 pure non-PAR `chrX.male` 上分窗，每名男性只贡献一个 X haplotype；因此 `total_haplotype_copies` 是实际男性 X chromosome 数，不是常染色体式的二倍样本数。PAR-scoped loci 使用男女样本的 diploid phased haplotypes。全染色体模式默认只做原始分窗比较；如确需开启分层 evidence，可设置 `PHYML_EVIDENCE_WHOLE_CHROMOSOME=TRUE`，但它仍是探索扫描，不能与 IBDmix/TRACE/AS3 tract calls 等同。

共享函数 `match_HAP`、`match_SNP` 和 `check_GRCH` 位于 `0f/0phe.f.sh`，不是 `0phe.f.R`。路径默认相对本仓库解析，也可用 `PHE_F` 覆盖。

## IBDmix

IBDmix 原生适合染色体/基因组扫描。GU 的 loci 模式默认计算扩展后的 `BED ± 100 kb`，再把 segment 裁回 analysis 区间，以降低局部窗口边缘效应；未经改动的 core 保存在 request map 中：

```bash
./gu.sh ibdmix --loci /mnt/d/files/gu.37.bed --loci-flank 100kb --grch 37
./gu.sh ibdmix --chr X --grch 37
```

代码会按所请求的染色体在 VCF header 中选择 `1/chr1` 或 `X/chrX/23`，不再错误使用 header 中第一条 `##contig`。whole/pure X 建立 `X_MALE` 单倍体分析单元；PAR loci 建立 `X_PAR` 全样本二倍体分析单元。IBDmix 需要二倍体 GT 接口时，只对男性 non-PAR 单倍体临时复制为同型二倍体，sample 仍只计算一次。pure chrX 需要 `samples.txt` 中的 `sample`、`sex`；PAR/常染色体分析使用全样本列表。

IBDmix 上游模型和常用评估以二倍体常染色体 genotype 为对象，并不原生处理男性 hemizygous X。GU 的复制只解决输入接口，不等于完成 chrX 的 FDR/LOD 校准；因此 normalize/Shiny 将该结果标为 `experimental male-X pseudo-diploid`，保留 LOD 供人工比较，但不把它纳入 chrX 的统一 evidence score。

## ARG 与 TRACE

TRACE 的输入不是 VCF，而是覆盖整条染色体的 mutation-bearing tskit `.trees/.tsz` ARG。`tsinfer` 与 `needle` 都可生成该格式；`--format trace` 把方法原生输出转换到互不冲突的 `arg/trace/<method>/`。Needle 转换会把样本节点表改成 TRACE 的 `tree_node_id/sample/haplotype` 契约，并明确将归一化节点时间标为 generations。

1KG VCF 的 `INFO/AA` 用作 tsinfer 的 ancestral allele；准备步骤会规范 `A|||` 等值，只保留 AA 等于 REF/ALT、低缺失、双等位 SNP 和 phased GT：

```bash
bash /mnt/d/scripts/gu/arg.sh build --method tsinfer --dir-gen /mnt/d/data.BIG/refGen/1kg/37 --chr 22,X
./gu.sh trace --chr 22,X --grch 37 --target 1kg --target-dir /mnt/d/data.BIG/refGen/1kg/37/pfile/chr

bash /mnt/d/scripts/gu/arg.sh build --method needle --format trace --dir-gen /mnt/h/ukbGen/37 --dir-pfile /mnt/h/ukbGen/37/hap --chr 22
./gu.sh trace --chr 22 --grch 37 --target ukb --target-dir /mnt/h/ukbGen/37/hap/chr --arg-dir /mnt/h/ukbGen/37/arg/trace/needle
```

执行顺序是 `arg.sh build → gu.sh trace`；ARG 检查是 TRACE 启动时的内在步骤。`--method tsinfer` 先把
`vcf/chr*.vcf.gz` 规范为持久的同级 `vcf.4arg/chr*.vcf.gz`，只保留有效
`INFO/AA` 与 `FORMAT/GT`，再生成 ARG。Needle 的原生产物仍位于
`arg/argn` 和 `arg/trees`；只有显式 `--format trace` 才会在
`arg/trace/needle` 建立 TRACE 文件，因此不会与 tsinfer 或 GRID 文件混淆。
ARG 只接受整染色体 `--chr`，Needle 与 tsinfer 均支持常染色体及 non-PAR X。

也可以使用 `--method threads`（PalamaraLab/threads，独立环境，不修改 Needle 环境）：

```bash
bash f/install_threads.sh
./arg.sh prep_gen --dir-gen /mnt/d/data.BIG/refGen/1kg/37 --method threads --format trace --threads 8 --jobs 1
./arg.sh build --dir-gen /mnt/d/data.BIG/refGen/1kg/37 --method threads --format trace --threads 8 --jobs 1
./arg.sh check --dir-gen /mnt/d/data.BIG/refGen/1kg/37 --method threads --format trace
```

Threads 默认覆盖 chr1–22、non-PAR X；加 `--chr 22` 可先跑一条染色体。
准备数据保存在 `pfile.4arg.threads/chrN.pgen/pvar/psam`，不使用 Needle 的 HAPS。
默认输出根目录是 `arg.threads`，原生 threading instructions 位于 `threads/chrN.threads`，
转换后树位于 `trees/chrN.trees`，TRACE 输入位于 `trace/threads`；支持 `--arg-dir` 自定义根目录。
`--threads 8` 实际启动最多 8 个 Threads 推断 worker（私有 Ray runtime），
整个任务仍受 `--memory-limit-gb` 限制；build 的 `--jobs` 仍限制为 1，避免同时运行多个整染色体。
准备、推断和转换都有输入/输出指纹检查点。转换失败后重跑会复用已经验证的推断结果。
`--action infer` 只推断，`--action convert` 从匹配的推断结果继续转换；普通 build 包含完整检查。

模型默认 `--threads-mode wgs --mutation-rate 1.25e-8`，常数 haploid Ne=20000
（等价 diploid Ne=10000）；这是明确的模型假设，不是从 1KG 数据拟合的历史。
`--demography FILE` 可指定从第 0 代开始的两列 time/haploid-Ne 文件；同一次请求的染色体共用该文件，
X 若需独立人口历史，应单独用 `--chr X --demography FILE` 构建。
`--threads-fit-to-data TRUE` 可打开上游拓扑一致性后处理，默认关闭；转换始终映射输入突变。
X 中每名男性保留一个真实单倍型、女性两个；上游 PGEN 存储要求偶数个单倍型，
所以沿用 `--x-odd-male drop-last|error`，默认必要时排除一个男性并写出 selection.qc.json。
真实样本身份保存在 haplotypes.tsv / sample_map.tsv，成对存储名不用于生物学分析。
不允许缺失基因型或未相位杂合子进入 Threads。GU 的 `--features`/`--affinity` 仍为 Needle 专用。
Threads → TRACE 的工程接口已接入；这不等同于验证古老混合检出的科学准确性。
上游文档：<https://palamaralab.github.io/software/threads/manual/>。

```bash
./gu.sh trace --chr 22 --grch 37 --target 1kg \
  --target-dir /mnt/d/data.BIG/refGen/1kg/37/pfile/chr \
  --arg-dir /mnt/d/data.BIG/refGen/1kg/37/arg.threads/trace/threads
```

TRACE 启动时会先验证：每条请求染色体存在 ARG；不是误标成 posterior 的坐标 chunk；site 数和覆盖跨度像整染色体；每个 chromosome-scoped run 内 node/sample/haplotype map 一致；缓存 provenance 未漂移。tsinfer ARG 要求 target 与 ARG 样本严格相等；Needle 的固定 panel 允许是 target 的子集，但 ARG 中每个样本必须存在于 target，TRACE 也只分析 ARG panel 中的节点。`f/arg.py infer-tsinfer` 从 VCF-Zarr 自动读取 ploidy，因此 pure male X 的每名男性是一个 ARG sample node。由于 X 是单倍体、常染色体是二倍体，`--chr 22,X` 之类的混合请求会自动拆成独立的 `chr22` 与 `chrX` TRACE 输出，完成后由 normalize 合并，不会把不兼容的 node map 强行联合。GU 随后依次执行：按 node 分组的 `trace-extract`、`trace-infer`、逐 haplotype 的 `trace-summarize`、标准化合并。

TRACE 默认保留完整染色体上下文，`--loci` 只在最终 segment 表中做 post-hoc clipping。若要做局部提取敏感性分析，可显式使用 `--trace-loci-mode extract`，此时才把 BED 传给 `trace-extract --include-regions`。未配置 HapMap genetic map 时，上游 TRACE 会按 `1e-8/bp/generation` 的均匀重组率运行。
对于 UKB，应把已相位的 UKB batch 与固定的 1KG anchor haplotypes 一起建 ARG，再对 UKB node 做 TRACE。只用 1KG 建树后“投影”UKB 没有树序列/祖先路径上的等价操作。

GRCh37 archaic callable masks 统一位于 `/mnt/d/data.BIG/refGen/archaic/37/mask/<sample>/chrN_mask.bed.gz`；`mask` 是永久 reference 数据目录，不是运行 cache。

## ArchaicSeeker3

标准运行使用独立的 `as3_mamba` 环境、GRCh38 和 chr1–22。GU 直接使用 AS3 官方预处理好的 Ref1028：每条染色体一个联合 reference VCF，配套共享的 `Ref_Panel.map.txt`。该面板含 146 个 AFR、3 个 Neanderthal 和 1 个 Denisovan 样本；GU 不再从分离的 archaic VCF 自建 reference panel，也不兼容旧布局。

本机布局由 `gu.env` 固定如下：

```text
/mnt/d/data.BIG/refGen/archaic/38/vcf/
  Ref_Panel.chr{1..22}.vcf.gz
  Ref_Panel.chr{1..22}.vcf.gz.tbi
  Ref_Panel.map.txt

/mnt/d/data.BIG/refGen/archaic/38/mask/
  chr{1..22}_mask.bed
  chr{1..22}_mask.wchr.bed

/mnt/d/data.BIG/refGen/archaic/38/models/
  Basemodel_Mamba_4096/{args.pckl,best_model.pth}
  Smoother_512_Kernel_8192/{args.pckl,best_model.pth}

/mnt/d/data.BIG/refGen/archaic/38/database_introgression_call/
  hg38_Neanderthal_introgression.bed
  chm13_Neanderthal_introgression.bed

/mnt/d/analysis/gu/normalize/
  reference-callsets.sqlite # normalize/Shiny 共用的官方 callset 标准化数据库

/mnt/d/tmp/gu/<target>/<build>/run.*/
  chr22/                  # chr22 自己的 VCF、source.tsv、ARG 临时 VCF
  chrX/                   # chrX 自己的 VCF、source.tsv、ARG 临时 VCF
  vcf/                    # 仅含指向各 chrN/ 文件的下游接口软链接
  as3-input/              # 仅限本次任务的 manifest
  as3-runtime/            # 仅限本次任务的 AS3 working directory

/tmp/gu/arg/<target>/<build>/run.*/
  chr22/                  # chr22 VCF-Zarr/LMDB；退出时删除
  chrX/                   # chrX VCF-Zarr/LMDB；退出时删除

/mnt/d/scripts/gu/f/as3/
  ArchaicSeeker3.1-mamba # GU 内置的最新版官方推断源码
  src/
```

| 文件 | 用途 |
|---|---|
| `vcf/Ref_Panel.chrN.vcf.gz`、`Ref_Panel.map.txt` | 最新 AS3 推断的直接 reference 输入 |
| `mask/chrN_mask.bed`、`chrN_mask.wchr.bed` | 官方 3N1D callable-mask 两种 contig 格式；GU 做完整性校验 |
| `models/` | 与 Ref1028 同级保存的固定模型参数和 checkpoints；启动时做 SHA-256 校验 |
| `f/as3/` | GU 内置的 AS3 源码；不读取 `/mnt/d/software/gu/as3` |
| `database_introgression_call/*.bed` | 已发表 callsets，仅用于外部比较 |

正式命令是：

```bash
./gu.sh as3 check --chr 22 --grch 38 --target 1kg --target-dir /mnt/d/data.BIG/refGen/1kg/38/pfile/chr
./gu.sh as3       --chr 22 --grch 38 --target 1kg --target-dir /mnt/d/data.BIG/refGen/1kg/38/pfile/chr
```

`as3 check` 会验证内置源码与模型 checksum、Ref1028 VCF/索引/contig、map 的 150 个样本及分类，以及两种 mask 文件；缺失或损坏时直接报出文件路径。官方 Ref1028 已完成 mask 相关预处理，因此推断直接传入 reference VCF 和 map，mask 只作官方数据清单校验，不会再次过滤 Ref1028。`--loci` 按染色体运行一次 AS3，再把结果裁到 loci。

最新版标准结果读取 `introgression.bed`，并同时保留 Denisovan、Neanderthal、Mosaic 分类文件。该文件虽以 `.bed` 命名，其 `start/end` 是 1-based inclusive VCF positions；GU normalize 会转换成统一的 0-based half-open 坐标。标准后处理按 10 万行分块过滤 raw BED，并单遍流式解压 SNP 明细，只为通过阈值的区段累计合并统计；不再把整条染色体的 SNP 表及逐区段 NumPy 数组同时装入内存。日志每扫描 25 万条 SNP 明细报告一次进度。AS3 只支持 GRCh38 chr1–22，`--chr X`（包括混在 `3,22,X` 中）会在 target 转换和 GPU 运行前明确拒绝。

每条染色体在全部 GPU chunk 落盘后原子写入推理检查点。若此后的 canonical post-processing 失败或被中断，GU 会保留 `introgression.raw.bed`、`introgression.raw.snps.gz` 和预测文本；用相同输入、模型、stride、ancestry 参数及 chunk size 重跑同一命令时，会显示 `RESUME ... stage=postprocess`，只重做流式 CPU 汇总而不重复 GPU 推断。检查点 provenance 不包含 `--merge`，因此可在保留同一份 raw inference 的前提下修改 merge distance。参数或输入签名不一致的旧临时目录不会被复用。

### AS3 官方 callset、haplotype 与 burden

`database_introgression_call/*.bed` 带有 population 与 sample 字段。标准化后只进入 `reference_callsets`，并标记为 `external_reference`；不会进入模型输入、`segments` 或 `carriers`。

若该可选目录丢失，可重新下载、校验并解压；无需重新下载 archaic VCF：

```bash
cd /mnt/d/data.BIG/refGen/archaic/38
wget -c -O database_introgression_call.tar.gz 'https://zenodo.org/records/14552025/files/database_introgression_call.tar.gz?download=1'
echo '86f0fd27fb123139fdaefa62b1682541  database_introgression_call.tar.gz' | md5sum -c -
mkdir -p database_introgression_call
tar -xzf database_introgression_call.tar.gz -C database_introgression_call
```

该目录缺失不会阻止 AS3，只会使 published overlay/overlap 为空。

官方 callset 约有两百万条 carrier-level 记录。只需运行统一的 final 整合命令：

```bash
./gu.sh final
```

`final` 会先快速检查官方参考缓存；缓存缺失或无效时自动构建，已有有效缓存时直接复用，然后继续处理各分析模块并刷新 Shiny 数据。参考缓存默认写到 `/mnt/d/analysis/gu/normalize/reference-callsets.sqlite`，保存 raw carrier 表和按 population 合并的 union 表。


- `reference_callsets`：`dataset_id, population, genome_build, chr, start, end, source_class, reference_role, raw_file`；官方 BED 的 carrier-level 重叠行会先按 population/build/chromosome 做 interval union，避免匿名个体行重复 overlay 或成倍放大 overlap；
- `reference_callset_overlaps`：仅比较 GU-AS3 Neanderthal calls；官方 BED 有 individual ID 时严格匹配同一个 1KG sample，再按 population 分层，保存 overlap bp、result overlap、reference overlap 和 reciprocal overlap。官方 carrier calls 只在标准化期间存在于 TEMP 表，不会进入 GU 的 `segments/carriers`；
- Shiny Segments 页的 published-callset 虚线 overlay 与 population overlap 汇总。

AS3 输出的 `Haplotype` 是当前 chunk 内的索引，不是稳定的 1/2。标准化器始终从 `SampleID_HapID` 拆出 `sample_id` 与 `haplotype`，原始 chunk 索引另存为 `method_haplotype_index`。

Individuals 页同时给出三种 burden：`raw_call` 是原始 calls 直接求和；`nonredundant_union` 按 sample/method/source_class/chromosome 做 interval union，避免 Altai/Chagyrskaya/Vindija 对同一 Neanderthal tract 重复计长；`consensus_catalog` 基于 `segment_catalog/carriers`，并另存 dosage-aware bp。

## UKB 的推荐顺序

1. 验证 UKB hardcall/相位状态与样本 sex 表。
2. 按可承受规模生成 batch；每批保留同一套 1KG anchor。
3. 为每批生成带 build-matched AA 的 phased VCF，把 UKB 与 anchor 一起建 ARG；ARG build 自己的派生 VCF 只在单次 `run.*` 中存在。
4. IBDmix、TRACE、AS3 分批输出；保留 `batch_id` 与 `genome_build`。
5. 运行 `./gu.sh final`，再运行 `./gu.sh shiny`。

标准化器先建立临时 SQLite、执行 `PRAGMA quick_check`，成功后原子替换正式数据库。各分析模块默认不再在每个 chromosome/batch 完成后重扫全部历史结果；应在一批任务结束后显式运行一次 `./gu.sh final`。小规模交互运行可使用 `--auto-final TRUE` 恢复自动重建。除原有 segment/burden 表外，final 还生成以下 PhyML-centered 表：

- `phyml_haplotype_carriers`：每一个被 PhyML 实际测试的 sample/haplotype copy，同时保留 pairwise `direct_match_pass` 与独立的 `tree_candidate_pass`；
- `locus_sample_support`：每个 locus × sample 的 tree-defined PhyML、IBDmix、TRACE、AS3 布尔支持；IBDmix/AS3 必须与 candidate tree lineage 一致，TRACE 作为 positional、lineage-agnostic 支持；
- `segment_catalog` 与 `coord_catalog` 按 target dataset 独立聚类并生成代码，1KG 与 UKB 不共享 cohort-level segment code；`sample_populations` 也使用 `(dataset_id, sample_id)` 复合键，避免跨 target 的同名 sample 串联；
- `locus_method_support`：按 `ALL` 和 1KG population 汇总 prevalence、Wilson 95% CI，以及每个 tract caller 对 PhyML-positive 个体的交集支持率；
- `method_runs`：从完成标记、provenance 和最终文件记录方法 × chromosome 的完成状态；即使一个合法运行产生零条 segment，也会显示为“已运行、阴性”，不会误报为“未运行”；
- `locus_evidence`：每个 locus 一行，保留各方法原生统计，但不再把 sequence identity、bootstrap、LOD、posterior 混成一个分数；`phyml_qc` 与 `evidence_state` 明确区分 core 太窄、缺祖先 outgroup、无支持 clade、未运行和独立方法 overlap；
- `locus_trajectory`：在 analysis interval 内按固定 bins 预计算两套曲线：全体 PhyML-tested cohort prevalence，以及以 PhyML-positive 个体为分母的 conditional tract support。分母不是全数据库样本数，也不是 raw segment 行数。

Shiny Overview 以五行 validation matrix 为主：预期角色、PhyML QC/tree carriers、IBDmix、TRACE、AS3 和结论；`not run`、`unsupported`、`exploratory` 分开显示。单-locus 解释与独立方法表保持在首屏，trajectory 和 individual rows 收入折叠面板。PhyML 页按“古参考 → tree candidate haplotypes → outside-clade controls”显示序列和树。

`f/prep_ukb.sh` 专门检查、导出 UKB 原始 BGEN/芯片数据，通过 `gu.sh ukb` 调用，提供 `inspect-hap`、`make-panel`、`batches`、`hap-vcf`、`hap-arg-vcf` 等动作。已经有 phased PGEN 的数据直接使用 `arg.sh prep_gen`，不需要经过这套 BGEN 导出流程。

## 代码入口

- `arg.sh`：ARG 的 `prep_gen`、`build`；`build` 自动验证输出，验证失败则返回非零状态。`check` 仅用于事后只读复查。运行 `./arg.sh -h` 查看 1KG/UKB 的准备和构建命令，`--help-all` 查看完整选项。
- `f/arg.sh` / `f/arg.py`：整合 ARG 参数解析、日志、准备、Needle/tsinfer、逐染色体调度及 TRACE 导出；`f/arg.sh trace_check` 是内部只读检查入口。
- `f/comm.sh` / `f/comm.py`：多模块共享的区域解析、目标数据检查、染色体长度、loci 扩展和宽 TSV 兼容逻辑。
- `f/shiny.R`：完整 Shiny 界面、服务端和双 SNP 展示，支持可选 `--review --data DIR` 独立复核页。
- `f/as3/`：第三方 AS3 源码及 GU 的环境检查、模型检查和输入准备脚本。
- `f/phyml.sh`：PhyML 流程调度；所有 Python 阶段统一通过 `f/phyml.py` 的 `compare`、`region`、`anchor`、`evidence`、`tree`、`plot` 子命令调用，阶段实现为 `f/phyml_*.py`，树绘图为 `f/phyml.tree_plot.R`。
- `f/prep_ukb.sh`：UKB 原始数据准备；`gu.sh ukb` 的命令和参数保持兼容。

`gu.env` 是本机运行配置，`gu.sh` 和 `install.sh` 都会读取，不能作为无用文件删除。目前它固定 AS3 Python 3.11；删除后安装与健康检查的默认版本会退回 3.9，还会丢失本机环境路径、模型路径及内存设置。

## 安装与检查

```bash
./install.sh
./install.sh --check
./install.sh --check-references
./install.sh --repair-as3
```

软件安装与大体量参考数据检查已解耦。通用 `gu` 环境承载 bcftools、tsinfer、TRACE、R/Shiny 等；AS3 保留独立的 `as3_mamba` 环境，并使用 `gu.env` 固定的 Python 版本，避免 PyTorch/CUDA/mamba-ssm 与通用环境互相污染。AS3 推断源码固定在 `f/as3/`，模型固定在 GRCh38 `models/`，安装和运行都不依赖外部 AS3 checkout。非 Blackwell GPU 默认验证上游 Torch 2.4/CUDA 12.1 profile；本机 RTX 5090（compute capability 12.0）使用兼容的 Torch 2.8/CUDA 12.8 profile。`./install.sh --check` 会明确打印实际 Python、Torch、Torch CUDA、GPU 和期望版本，而不是只检查 import 是否成功。AS3 运行日志中的 `physical_gpu_id=0` 或旧版 `gpu=0` 表示选择物理 GPU 0，并非关闭 GPU；该卡在设置 `CUDA_VISIBLE_DEVICES=0` 的子进程中映射为 PyTorch 的 `cuda:0`。GU 默认传递 `--target-chunk-size 64`，避免大队列一次性载入全部样本而耗尽 WSL 内存；可用 `--as3-target-chunk-size INT` 显式覆盖。每个 AS3 Python 子进程还会运行在独立 cgroup 中，默认 `MemoryHigh=20G`、`MemoryMax=24G`、`MemorySwapMax=8G`。即使出现异常内存增长，OOM 也限制在该 AS3 scope 内，不会再成为 WSL 全局 OOM；这些值可分别通过 `gu.env` 中的 `AS3_MEMORY_HIGH`、`AS3_MEMORY_MAX` 和 `AS3_MEMORY_SWAP_MAX` 调整。

若 AS3 日志出现 `numpy.compat`、Dask array requirements、`torchvision` 或 h5py/HDF5 不匹配，说明独立环境已发生依赖漂移。网络可用时在仓库目录运行 `./install.sh --repair-as3`；该动作不会更新通用 `gu` 环境，会保留匹配的 Torch/CUDA profile，并按 `f/as3/requirements-gu.txt` 修复数据栈。`gu.env` 中的 `GU_CONDA_BASE` 与 `AS3_PYTHON` 应指向同一 Conda 安装，避免机器上同时存在 Anaconda/Miniforge 时修复与运行使用不同环境。随后先运行 `./install.sh --check`，再重跑原 AS3 命令；target 转换和 manifest 会按本次任务重新生成并在退出时删除。

## 主要上游

- Zeberg & Pääbo 2020 Nature：<https://www.nature.com/articles/s41586-020-2818-3>
- PhyML：<https://github.com/stephaneguindon/phyml>
- IBDmix：<https://github.com/PrincetonUniversity/IBDmix>
- TRACE：<https://github.com/YulinZhang9806/trace>
- ArchaicSeeker3：<https://github.com/Shuhua-Group/ArchaicSeeker3.0>


## ARG 内存保护

`./arg.sh` 默认通过 user systemd / cgroup v2 限制整个任务树为 32 GiB RAM + 4 GiB swap。`build` 每次只运行一个染色体；即使传入 `--jobs 4` 也会明确提示并降为 1。准备阶段仍可使用多个 jobs，但受同一总内存上限约束。缺少 cgroup v2 或 user systemd 时拒绝启动，不退回单进程 ulimit。

使用 `--memory-limit-gb N`（正整数）和 `--memory-swap-gb N`（非负整数）调整；环境变量 `ARG_MEMORY_LIMIT_GB` / `ARG_MEMORY_SWAP_GB` 设置默认值。保护适用于公共入口 `./arg.sh`，内部 `f/arg.sh` 不应直接用于启动任务。单染色体仍可能超过上限，届时任务会失败；32 GiB 是上限，不是完成推断所需内存的保证。该保护不改变样本、SNP、推断模式和已有检查点，也不约束修改前已启动的任务。

Needle 的调用器 `../grid/f/run_argneedle_advanced.py` 自动加载同目录的
`argneedle_memory.py`，将 hashed threading 的 posterior 缓存由
`(scaffold-1) × SNP` 的 float64 矩阵改为每个 SNP 一个选中值。例如 chr1
的 4000 scaffold / 2,054,102 SNP，这一项从约 61.2 GiB 降到 15.7 MiB。
平滑区间内选中 cousin 不变，因此仍对同样的 float64 值求均值；不改变
样本、SNP、哈希候选、ASMC padding 或随机种子。另将原生 ASMC `batchSize`
从 64 改为 16，保留全部候选对和 SIMD。独立构建的 ASMC 1.4.0 后端
在每个窗口输出结果后，用 Linux `MADV_DONTNEED` 释放前向/后向 scratch
缓冲区内部的完整页，避免窗口推进时整条染色体的工作区逐渐常驻内存。
此优化不消除输入哈希和 ARG
本身的内存需求。上游函数不兼容时明确报错，不静默退回大矩阵。
修改仅作用于本次 Python 进程，不改写 site-packages；调用器的检查点指纹
包含优化模块、独立后端二进制和 decoder 源码，准备好的基因型输入仍可复用。
后端还将 SIMD 大数组索引改为 64 位，避免 chr1 尾部的索引超过 int32 范围后崩溃。
后端位于 `~/.cache/gu/asmc-memory/v1.4.0-p2/<Python ABI>/`。迁移机器或
更换 Python 时，在 grid 环境运行
`python3 ../grid/f/build_asmc_memory.py --source /path/to/ASMC-1.4.0` 重建。
源代码使用上游 `PalamaraLab/ASMC` 的 `v1.4.0` tag；构建需要 CMake、
C++17 编译器、Boost 和 zlib，其余依赖由上游 CMake 按固定版本获取。
缺少修补后端时会明确报错，不自动回退到原版的大内存实现。

使用 grid 环境的 Python 运行 `../grid/f/tests/test_argneedle_memory.py
--haps /path/to/chr1.haps.gz --map /path/to/chr1.map`，可对真实输入的小子集
比较原实现与优化后的 sequence/array 三阶段树表（包含 normalization），
并验证断点复用。`../grid/f/tests/probe_argneedle_memory.py` 用完整输入及
指定 scaffold 大小测试最初几个样本；它不写检查点，不代表整条染色体已完成。
也可加 `--probe-windows 12`，反复解码覆盖染色体首、中、尾的 12 个代表窗口，
检查窗口切换后的缓冲页释放及重复计算的一致性；此模式不完成整个样本的 threading。

### final 模块

`./gu.sh final` 是最后的分析整合模块：整合各方法的证据，生成供 Shiny 展示及论文图表使用的最终结果。`./gu.sh shiny` 是展示入口，不是分析模块。现有 `normalize/` 输出目录、参考缓存及 `gu.sqlite` 路径保持兼容。
