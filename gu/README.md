# GU 代码目录

公共入口保持 `gu.sh`、`arg.sh`、`install.sh`；内部代码按方法合并，准备和公用工具以 `0.` 开头。

| 模块 | shell | Python / R |
|---|---|---|
| 公共配置、区域和表格工具 | `f/0.common.sh` | `f/0.common.py`、`f/0.common.R` |
| UKB、祖先等位基因准备 | `f/0.prep.sh ukb` / `ancestral` | 共用 `f/0.common.py` |
| ARG 公共步骤及 Needle/tsinfer | `f/arg.f.sh`，通过根目录 `arg.sh` 运行 | `f/arg.py` |
| SINGER | 通过根目录 `arg.sh --method singer` | `f/arg.singer.py` |
| Threads / 私有 Ray worker | 通过根目录 `arg.sh --method threads` | `f/arg.threads.py`（含 `worker` 子命令） |
| ArchaicSeeker3 | `f/as3.sh`（含 `prepare` 子命令） | `f/as3.py`（`health`、`model-check`） |
| COJO | 主入口调用 | `f/cojo.py`、`f/cojo.R` |
| PhyML | `f/phyml.sh` | `f/phyml.py`、`f/phyml.R` |
| IBDMIX | `f/ibdmix.sh` | `f/ibdmix.py`（含 `provenance`） |
| TRACE | `f/trace.sh` | `f/trace.py`（`output`、`combine`） |
| 结果归一化、外部参考缓存、密度 | 主入口调用 | `f/normalize.py`（含 `reference-cache`、`density`） |
| 结果复核及双 lead 分析 | 主入口 / Shiny 数据准备 | `f/review.py` |

进程组内存限制是跨项目的功能，GU 和 GRID 共用 `../0f/memory_cap.sh`。GU 专用的区域、染色体及结果约定仍放在本项目 `f/0.common.*`；已有 `../0f` 手写工具保持原有接口。

需要 GRID 的 ARG 辅助功能时，调用 `../grid/f/0.arg.py` 的相应子命令。第三方 AS3 程序、模型来源、许可证和依赖锁定文件保留在 `f/as3/`。SINGER / Threads 的安装统一通过 `install.sh --singer` / `--threads`。

正式结果按方法、数据集和位点／染色体组织。每张 PhyML PNG 配同名 XLSX，记录该图实际显示的单倍型计数、分支长度、节点支持度和判定结果；不再生成同图 PDF。树的统计推断和显示筛选方法保持不变。

所有主要结果永久保存在 `/mnt/d/analysis/gu/`。格式按用途决定：供查看的结果表使用 XLSX（包括含 ID 的表）；可复用原生数据、拟合结果使用紧凑归档；大规模关系数据库另有 RDS 恢复副本。不会仅因为有 `sample_id`、`eid` 等列就强制保存 RDS。

| 内容 | 永久文件 |
|---|---|
| PhyML | 每个位点的 `phyml.haplotypes.xlsx`、树图 PNG 和同名 XLSX、`phyml.raw.tar.gz` |
| IBDmix | 每个范围的 `ibdmix.tracts.xlsx`、`ibdmix.raw.tar.gz`、简要运行参数和实际调用样本名单 |
| TRACE / AS3 | `trace.segments.xlsx` / `as3.tracts.xlsx`，以及对应的 `*.raw.tar.gz` |
| Final / Shiny | `final/gu.sqlite`、汇总和复核 XLSX、`final/normalize/`、`final/review/` |
| 数据库恢复 | `final/gu.results.rds` 保存大型关系数据，供数据库恢复复用 |
| UKB / COJO | `ukb/` 的准备结果和 `phyml/<数据集>/inputs/` 的 lead 转换结果 |

计算中的原生结果先写入永久目录。完成后，将可复用序列、树、片段和检查点整合进每次分析的 `*.raw.tar.gz`，逐文件校验后移除分散副本。继续分析时按需恢复到永久目录；Shiny 使用 `/tmp/gu-native-view/` 中可随时重新解压的读取副本。删除 `/tmp` 不会删除唯一的分析结果，重新运行 `./gu.sh shiny` 会恢复所需读取副本。

临时日志、锁、生成的 mask、预处理 VCF、命令列表和计算工作区放在 `/tmp`。PhyML 位点的完成日志和命令保留在永久目录；每个位点的参数和来源记录保存在原生归档内。XLSX 使用多个工作表整合结果，并保留精确表格导出供代码读取；超出 Excel 行数上限时拆分工作表，过长文本按顺序拆到附加工作表，不截断内容。发布使用暂存文件校验后替换。

`GU_ANALYSIS_ROOT` 与 `GU_PUBLISHED_ROOT` 使用同一永久目录；结果路径不能指向 `/tmp`、`/var/tmp`、`/dev/shm` 或 `/run`。跨项目的 XLSX/RDS 读写共用 `../0f/results.R` 和 `results.py`。格式整理不重新拟合模型，也不改变科学结果或伪造分析完成状态。

代码格式遵循公用 R 代码：运算符两侧空格、tab 缩进、使用 `# 🚩` 划分模块。第三方源码保留上游格式。

本次整理核对了全部已保存的科学结果和 258 棵完整树，未重新推断树或调用 IBDmix／TRACE。TRACE 目录当前保留的是可继续分析的提取数据；整理不会把未完成的阶段标成已完成。

## PhyML 运行与续跑

```bash
PHYML_TREE_CPUS=4 PHYML_TREE_TIMEOUT=7200 ./gu.sh phyml \
  --loci /mnt/f/gwas/main/common/bald0/gwas/bald0.jma.cojo \
  --loci-format cojo --grch 38 --target 1kg \
  --target-dir /mnt/f/gen/1kg/37/pfile/chr \
  --jobs 6 --memory-cap 32G --replace-phyml FALSE --foreground TRUE
```

默认整批共用 32 GiB 内存上限，每棵树最多计算 7200 秒（2 小时）。参考版本检查使用 `f/0.build-check.sh` 的带锁缓存：同一源文件、哨兵表、目标 build 和校验代码只扫描一次；文件路径、大小、纳秒修改时间或校验代码变化时重新检查，失败不缓存。

每个位点先完成输出校验、写入 `.phyml.locus.complete.json`，再发布和压缩。完成记录、位点 `.cmd` 和 `.log` 保留在永久目录及原生归档中；下次运行会核对输入和输出后跳过已完成位点。

超时继续处理其余位点，记录到 `phyml/<数据集>/phyml.timeouts.tsv`，包含位点、序列数、位点数、转入清单前已保存的 bootstrap 树数、耗时、输入和命令路径。输入 alignment、日志和超时记录保留；不完整的树不进入有效树结果。普通续跑跳过已记录的超时位点。已中断的旧树若输出时间证明它已超过本次预算，也直接列入待处理清单，`reason` 标明旧运行中断。

以后补跑时在同一入口命令前设置 `PHYML_RETRY_TIMEOUTS=TRUE`，并按需要覆盖 `PHYML_TREE_TIMEOUT`、`PHYML_TREE_CPUS`；成功后会从清单移除。仅增加超时或 CPU 数不会自动重跑清单，也不会替换已完成的树。不完整 bootstrap 暂不支持接着续算。

运行控制回归测试：`python3 tests/test_phyml_runtime.py -v`。

## 2026-10-09 PhyML 全量科学复核

对照 [Zeberg & Pääbo, Nature 2020](https://www.nature.com/articles/s41586-020-2818-3)，逐一读取并校验本批 709 个位点的原生归档。原始 COJO 有 711 条记录，2 条因 SNP ID 与 Chr/bp 冲突被排除；成功转换为 GRCh37 的 709 条均有状态记录。保存结果中 352 条高 LD 标记不足、77 条风险/非风险序列无法区分、9 条 lead 缺失或不唯一、9 条建树位点不足、1 条重复单倍型不足、18 条建树失败（17 超时、1 数值失败），只有 243 条完成树。未完成和不可评估不能算作阴性。

原有 `tree_pass` 检验的是全部重复风险单倍型与该谱系全部参考共同形成排他分支，BS ≥70；这是一个严格的特定假设，并不是完整的渗入发现算法。论文在 COVID 位点观察到了这样的风险分支，不意味着所有性状中的渗入都必须增加表型，或所有带同一 lead 等位基因的序列都必须共同成支。新增 `nonrisk_tree_pass`、`nonrisk_tree_bootstrap`、`supported_allele_role`，在同一完整树上按同样规则复核另一等位基因；原有风险方向调用和携带者验证口径保留。非风险方向不自动解释为临床保护效应。Overview 可查看任一方向支持、两个方向各自支持或全部输入 lead，并显示输入/完成分母。

本批已存树中，Neanderthal 原有风险方向支持 1 个 locus，非风险方向另有 4 个；Denisovan 风险方向 1 个。4 个新增非风险信号的输入 ID（GRCh38）和 BS 是：`20:47731721:C:T` (72)、`5:40059671:G:A` (81)、`7:47009169:A:G` (100)、`7:522772:G:A` (93)。这些是树支持，不能直接改名为已确认渗入。审计中的单参考/部分现代单倍型分支另列为探索性亲缘证据，不进入严格调用，也不把扫描分支的 BS 当作已校准的全基因组 P 值。

已修复 GRCh38→37 后 REF/ALT 互换造成的精确 lead 匹配失败：比较同一坐标上的无序双等位集合，随后按目标 VCF 顺序解码 GT，保持效应等位基因方向；第三等位基因和多个匹配记录仍不接受。`8:116228954:A:G` 在 GRCh37 的真实记录是 `8:117241193:G:A`，独立复核恢复了 lead、6 个高 LD 标记和 28 个建树位点，但其序列仍无法区分两种 lead 等位基因，因此没有被提升为阳性。另修复 Denisovan 报告沿用 Neanderthal 阶段状态的问题，以及系统 PhyML shell 启动器在串行模式下自行启动 MPI 的问题。

过滤影响也已量化：有逐位点 QC 的核心共 139,784 个 SNP（不同 lead 的重叠区间重复计数），三 Neanderthal 共调用为 80,369，五参考共调用为 79,910。额外两个 Denisovan 参考减少 459 个位点；77 个序列无法区分的 loci 中，仅 2 个的 lead 能在三个 Neanderthal 中共调用。352 个高 LD 不足的 loci 中，218 个在 r² >0.8 时有至少两个标记；这里只报告敏感性数量，没有降低主分析的 r² >0.98 阈值或重新拟合这些位点。243 棵完整树中，20 棵不足 10 个 SNP，76 棵五个古参考序列完全相同，不能把难以区分谱系当作无渗入证据。

论文使用的 0.53 cM/Mb 是其 chr3 区域的局部重组率，不能用于所有 COJO loci 的正式 ILS 显著性判断。当前固定重组率的结果仍仅作 LD 跨度敏感性指标；严格渗入判断还需要共享衍生变异、实际共享片段边界、局部重组/ILS 和独立片段证据。1–4% 是个体基因组碱基的祖源比例，不是经 GWAS/COJO 筛选后 loci 为渗入来源的概率，也不提供最少阳性数量；即使假定 709 次独立、同概率抽样，7–28 也只是期望值。

复核输出在 `/mnt/d/analysis/gu/final/review/phyml_audit_20261009/`，包含工作簿、阶段统计、全部 loci、双方向树支持、探索性分支、逐参考缺失原因、输入排除记录和校验清单。独立的论文阳性对照序列重建得到 GRCh37 `3:45859651–45909024`、450 个 SNP 和 253 种现代单倍型，与论文报告一致。原始 GWAS 归档不改写；普通续跑会重新检查旧规则下报告 lead 缺失的 9 个位点一次，其余有效完成缓存继续复用。

可重复执行审计；`--refresh-report` 可选，仅刷新现有概览的诊断列，先备份旧表，不重新推断树或改写原始风险调用：

```bash
python3 f/phyml.py audit \
  --dataset-dir /mnt/d/analysis/gu/phyml/1kg \
  --lead-table /mnt/d/analysis/gu/phyml/1kg/inputs/bald0.jma.cojo.dbac4a5b8be8/gwas_leads.GRCh37.tsv \
  --output /mnt/d/analysis/gu/final/review/phyml_audit_20261009 \
  --refresh-report /mnt/d/analysis/gu/final/review/phyml_locus_report.tsv
python3 tests/test_phyml_science.py -v
python3 tests/test_phyml_runtime.py -v
/home/huangj/anaconda3/envs/gu/bin/Rscript --vanilla tests/test_shiny_evidence.R
```

## 2026-10-06 chrX 复核更新

先运行 `./gu.sh final`，从正式原生结果重建跨方法证据、完整 PhyML 位点状态、X 汇总和工作簿。此步骤会从保存的 Newick 刷新树图及同名 XLSX，不重新推断树；完整图和剪枝图分别校验无根分裂的 bootstrap。展示代码的变更会刷新导出缓存。旧格式工作簿也可单独迁移：

```bash
python3 f/0.common.py results upgrade-workbooks --published /mnt/d/analysis/gu
```

`Altai.2013` 和 `Denisova.2013` 仅保留在 IBDmix 原生结果及复现汇总。规范化数据库保留 `source`、`reference_role`，跨方法片段目录、携带者统计、PhyML validation 和 review 使用较新参考白名单。五参考 PhyML 的原有整体拓扑检验及阈值保持不变。

X 默认使用独立的 `nonpar-v2` 评分配置：modern 和 archaic 输入统一采用对应 build 的非 PAR 区间，检查 genotype 表的所有位置，各连通区间分别调用。常染色体的 caller 参数、mask 算法和 `pipeline_version` 不变。新的 X 结果写入 `chrX.nonpar-v2.mac1`（位点任务也加同样后缀），保留旧 `chrX` 基线。完整新结果出现后，final 优先采用新结果，不与旧基线取并集。

`IBDMIX_X_PROFILE=legacy` 可显式复用旧 X 基线。`IBDMIX_X_MINOR_ALLELE_COUNT=2` 或 `3` 是独立敏感性运行，写入 `.mac2` / `.mac3`，要求先有相同 calling cohort 的完整 MAC1 基线；这些敏感性结果不进入默认 final。需要自定义基线位置时设 `IBDMIX_X_BASELINE_DIR`。这些设置仅用于 X，不改变常染色体参数，也不代表 X LOD 已校准。

PhyML 概览保留所有输入 lead，并分开显示输入、LD、alignment、建树、拓扑及序列 QC。新增运行保存输入人数、lead 可调用人数／拷贝、重复及单次拷贝、modern 筛选原因和逐参考位点可调用性。历史结果中能从已存序列恢复的指标在 final 阶段重算；没有保存的输入分母或缺失原因显示未知，不能用后续保留人数代替。

`final/gu.ibdmix.summary.xlsx` 保存各参考／人群的检测、零片段和未检测人数、Mb、物理分母、染色体范围、过滤状态及同一男性队列的 X／22 条常染色体比较。全景和区域放大的 X 热图均采用 bin 与 non-PAR 的交集。仅含汇总的审核包也能运行 `./gu.sh shiny`；完整环境可用 `GU_SUMMARY_ONLY=1 ./gu.sh shiny` 强制汇总模式。没有数据库或序列的功能会明确标为不可用。

LD span 与实际匹配 tract 边界分别保存；固定 0.53 cM/Mb 的 ILS 数值只标为 LD 跨度敏感性统计。可选 `PHYML_X_GENETIC_MAP=/path/map.tsv` 接受 GRCh37、1-based 的 `chr`、`pos`、`cM` 列，同时必须设 `PHYML_X_MAP_SEX_CONVENTION=sex_averaged`、`female_meiosis` 或 `historical_X`。程序只在图谱覆盖内插值，保存来源 SHA-256 和性别惯例；不将男性历史 X 重组率置零，不把该长度当已证实的古人类共享连续片段。

运行 metadata 的变长记录会完整矩形化导出；mask 组件状态、允许范围 bp、mask 校验值与 genotype QC 保存到原生结果 `qc/` 并随归档保留，清理 `/tmp` 后可恢复。历史未保存的 mask 执行状态仍为 unknown。

本地基本验证已覆盖：Shell/Python/R 语法、`gu.sh --help`、真实 `generate_gt` 的 PAR 输入不变性、真实 IBDmix 分区调用、2013 片段隔离、36 个保存的 X 位点归一化、15 棵保存树的完整／剪枝重定根、旧 XLSX 迁移与校验、完整及汇总模式 Shiny 启动。EUR 的 X Altai 和 Altai.2013 片段并集合计分别保持 115,785,748 bp 和 78,936,509 bp，分母均为 240。测试数据与脚本在 `/tmp/gu-update-20261006/`；没有替用户执行全量 final、重新推断树或重跑全染色体 caller。实际 X 灵敏度、局部图谱及模型校准仍需正式运行和研究验证。

## 2026-10-07 UKB PGEN 与样本筛选

PhyML 使用已定相的 `hap`，IBDmix 使用较密集的 `imp`。统一参数是 **`--keep FILE --keep-males`**：名单按 PLINK2 规则传入 `--keep`；`--keep-males` 原样传给 PLINK2，性别来自输入 `.psam`，女性和未知性别均排除。没有 `--keep-males` 时保留所选常染色体样本的全部性别；男性非 PAR X 仍独立选取男性。不会将 `/` 伪造为 `|`。

```bash
./gu.sh phyml \
  --loci /mnt/f/gwas/main/common/bald0/gwas/bald0.jma.cojo \
  --loci-format cojo --grch 38 \
  --target ukb --target-dir /mnt/f/gen/ukb/37/hap/chr \
  --keep /path/batch01.txt --keep-males

./gu.sh ibdmix --grch 37 \
  --target ukb --target-dir /mnt/f/gen/ukb/37/imp/chr \
  --keep /path/batch01.txt --keep-males
```

旧命令中的 `--target 1kg` 若搭配明确的 `/ukb/` 路径，也会自动识别为 UKB，不会覆盖 1KG 结果。建议显式写 `--target ukb`。名单支持 IID 单列或 PLINK FID/IID 格式；带表头时使用 `#IID` 或 `#FID IID`。保留字符串 ID 的前导零。常染色体缺少请求样本、或 FID/性别不一致时会停止。直接 UKB 入口允许 chrX 的样本数少于常染色体：仅在 X 上排除 PSAM 中不存在的男性，保留原常染色体名单；日志、准备 QC 和 `inputs/ukb.excluded.samples.tsv` 明确记录排除情况，`ukb.scored.samples.tsv` 记录实际分析者。缺失者不当作零片段。

不同 keep 内容或性别筛选生成独立的 `ukb-cohort-<hash>` 数据集和默认结果根目录 `/mnt/d/analysis/gu-ukb-cohort-<hash>/`。同一名单及样本 metadata 的 `hap`/`imp` 使用同一 cohort 身份，分别写入 PhyML/IBDmix 子目录，以便 final 按同一人交叉核对；每个方法仍分别核对 source/build/QC，不能用另一套 PGEN 覆盖已完成结果。可通过 `GU_ANALYSIS_ROOT` 显式选择独立结果根目录。1KG 使用筛选参数时也写入独立 subset，未使用筛选参数的原路径和 caller 不变。

每个结果根目录保存供汇总使用的 `gu-target.env`。结束分析后，在子 shell 中加载日志所示的对应文件：

```bash
(
  source /mnt/d/analysis/gu-ukb-cohort-<hash>/gu-target.env
  ./gu.sh final
  ./gu.sh shiny
)
```

`<hash>` 替换为实际目录名。这个持久配置仅用于 final/Shiny；重新推断仍使用上述带 PGEN 路径和名单的命令。已归档的实际样本名单和准备 QC 支持 `/tmp` 清理后的恢复。

直接 PGEN 入口使用现有 hardcalls，不额外重设 dosage 阈值、INFO 或缺失率阈值；QC 明确记录 `existing_PGEN_hardcalls;INFO_not_reassessed`。它检查样本、唯一坐标、双等位 SNP、FASTA REF、实际 GT、倍性及 PhyML 所需 phase，但不声称重新验证了原 imputation 质量。IBDmix UKB 入口只将现代数据实际保留的 SNP 送入古参考评分，现代或古参考缺失记录不补为参考纯合；这项适配位于准备层，原 IBDmix 调用文件及常染色体参数、mask 代码没有修改。

UKB PVAR 中重复的 rsID 使用临时 PVAR 行编号提取，输出 ID `gu_row_N` 对应源 PVAR 的第 N 条变异记录；保留不同坐标的同名位点，仍排除同一坐标存在多个双等位 SNP 的情况。原 PGEN/PVAR/PSAM 不改写。若 PVAR 的 ALT 才是 FASTA 参考碱基，用 PLINK `--ref-from-fa force` 同步调整 REF/ALT 与 GT 编码，再执行严格 FASTA 校验。男性非 PAR X 的无效杂合 hardcall 通过 `--set-invalid-haploid-missing` 置缺失；QC 记录原有缺失数和本步骤新增缺失数，不将其补为任何等位基因。此步骤不应用于常染色体。

UKB PhyML 的 LD 由本次所选 UKB 样本计算，仍使用原来的 phased r² > 0.98、重复单倍型和五参考整体树检验。祖先等位基因从 GRCh37 1KG VCF 的 `INFO/AA` 按 POS/REF/ALT 精确匹配，仅转移注释，不借用其现代基因型；未匹配保持 N。可用 `PHYML_ANCESTRAL_VCF_DIR` 指定另一套同 build 的 AA 注释 VCF。没有 ancestry panel 时标注 UKB/ALL，不标为 EUR。`hap` 位点较稀疏，原始 lead 缺失、LD 位点不足等状态如实保留，不以邻近 SNP 替代。

UKB IBDmix 默认使用 Altai/Vindija，可通过 `--archaic-gen` 选择现有五个较新参考；2013 参考保留给原复现入口。UKB 不具备原 AFR5 对照，因此此入口关闭 AFR Denisova background subtraction，并把结果标为未校准的参考匹配候选，不计入严格跨方法证据。男性 X 从 `chrX.pgen` 按 PSAM 选出男性，检查 non-PAR，独立调用不相连区间。未评估范围和女性不会补为零。

也可用 **`--keep-psam FILE,CHUNK_SIZE,CHUNK_INDEX`** 自动分批，例如：

```bash
./gu.sh ibdmix --grch 37 \
  --target ukb --target-dir /mnt/f/gen/ukb/37/imp/chr \
  --keep-psam /mnt/f/gen/ukb/37/imp/chr1.psam,2000,3 --keep-males
```

chunk 编号从 1 开始，按 PSAM 原始样本顺序计数，表头、空行和 `##` 注释不计入人数。上例先取第 4001–6000 人，再保留其中男性，因此实际分析人数可能少于 2000。最后不足 2000 人的一块照常处理；编号越界、非正整数、没有剩余男性均报错。`--keep` 与 `--keep-psam` 互斥，二者均可与 `--keep-males` 同用。FID/IID 保留前导零；按 ID 匹配目标染色体，不要求不同染色体的 PSAM 行顺序相同。

UKB 的 `--memory-cap` 是整个任务及其所有 worker 共用的内存上限，默认 `32G`；PLINK 的 `--memory` 仅控制单个进程的工作区。调度器为其他进程保留总限额的 25%（至少 1 GiB），并额外预留每个 worker 2 GiB，按剩余预算限制 `--jobs`。默认 32G 下最多同时运行 2 条染色体，每个 PLINK 工作区 8192 MiB；16G 下最多 1 条，8G 下使用 1 条及 4096 MiB 工作区。`--jobs 1` 仍保持串行，PhyML 默认 1 个任务。资源预算及实际并发数写入详细日志；生成的 `.cmd` 保存相同的 PLINK 内存设置。可用 `GU_PREP_MEMORY_MB` 指定工作区 MiB，超出可用预算会在启动分析前报错。这是并发预算，不能保证任意样本数下的完整推断都不会超限。loop 中依次等待各批完成，避免多次独立调用的内存消耗叠加。

IBDmix 原有的每个人群至少 10 人要求仍然生效，按最终保留人数检查；尾块筛选后不足 10 人时应调整分块大小。

PhyML 将同一参数加到上面的 `hap` 命令即可；若要在 PhyML/IBDmix 间使用同一批人，应使用同一份 PSAM 分块名单。生成的 worker 保留分块参数及源 PSAM 校验值；调度后源文件改变会要求重新生成命令。准备 QC 归档记录源文件、校验值、chunk 编号和样本范围。

可以在 loop 中替换最后的 chunk 编号；`--foreground TRUE --jobs 1` 让每一批完成后再继续下一批。例如先试跑第 1–3 批：

```bash
for chunk in 1 2 3; do
  ./gu.sh ibdmix --grch 37 \
    --target ukb --target-dir /mnt/f/gen/ukb/37/imp/chr \
    --keep-psam "/mnt/f/gen/ukb/37/imp/chr1.psam,2000,$chunk" \
    --keep-males --foreground TRUE --jobs 1 || break
done
```

不同批次的实际样本名单分别生成独立结果身份，避免相互覆盖。IBDmix 的 AF/MAC 在每批中计算，分批结果不等价于全队列一次调用；程序保留 cohort 身份，不自动把不同批次拼成同一次统计分析。

下载方案中的严格准备入口也已保留：`./gu.sh ukb inspect-pgen`、`make-panel-pgen`、`pilot`、`pgen-vcf`。`pgen-vcf` 默认要求真实 INFO/MFI 或明确的预先 QC 位点名单，对有 dosage 的 imp 先重新生成 hardcalls，再做缺失率过滤和完整 GT 审核。详见 `./gu.sh ukb pgen-vcf --help`；它与上述“直接使用现有 PGEN”入口的 QC 语义不同。支持 `--ukb-info-file`、`--ukb-info-min`、`--ukb-hardcall-threshold` 等下载方案中的参数。

本次小规模验证覆盖真实 PLINK2 的名单与男性筛选、PSAM 分块边界及互斥参数、常染色体和非 PAR X 的 PGEN 导出、PhyML 的 hap 输入／LD／序列准备／归档、原版 IBDmix 输入检查、UKB 结果归一化和已有 1KG 回归测试。未运行 UKB 全量推断。本次 `ibdmix.sh`、`ibdmix.py` 与 `0.common.sh/.py/.R` 均与更新前逐字节一致。
