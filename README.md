# 代码与分析结果约定

本文件统一记录 `/mnt/d/scripts/`（Windows：`D:\scripts\`）下所有项目的整理和开发约定。后续修改沿用这些规则，不另建 `AGENTS.md`。

## 通用下载校验

`download_check.py` 下载缺失文件，用现有的权威 MD5 检查已有文件，并修补混入二进制数据的 `ConnectionClosedException` 错误文本。仅使用 Python 标准库。例如从本目录运行：

```bash
./download_check.py --dir /mnt/f/metag --url 全球污水.url --folder 全球污水
./download_check.py --dir /mnt/f/metag --url 厦门污水.url --folder 厦门污水
./download_check.py --dir /mnt/f/metag --url 纳米布沙漠.url --folder 纳米布沙漠
./download_check.py --dir /mnt/f/metag --url 长江.url --folder 长江
```

相对路径均相对于 `--dir`；默认读取 URL 清单同名的 `.md5`，也可用 `--md5` 指定。MD5 必须来自可信来源，不可用损坏的本地文件生成期望值。`ena-md5.metadata.tsv` 存在时还会核对其中的校验值与文件长度。普通 URL 清单可含 aria2 的 `out` 和 `checksum` 选项；其他逐文件选项会明确报错。

默认流程：缺失文件下载到同目录的 `.download_check.part`；已有文件及新下载文件均先检查 MD5，不符时尝试局部修补，完整 MD5 一致才算成功。新下载文件通过校验后才发布为正式文件，并删除临时文件；中断时保留临时文件供下次续传。完整但校验失败的临时文件也会检查和尝试修补，不凭文件大小当成成功。网络请求默认直连，可加 `--use-proxy`。已有文件无法局部修复时明确报告失败，不自动重下整个文件。`--check-only` 只校验，不下载；`--plan` 只验证清单，不创建文件；`--only 文件名` 可选取少量文件。请避免同时下载和修复同一目录；程序会阻止并行检查及旧修复任务，跳过带 `.aria2` 控制文件或检测到正在写入的文件。

每份 URL 清单维护两份记录：`.download_check.jsonl` 合并历史结果和修补前后的字节备份；`.download_check.fail.tsv` 保存未解决的失败，列为 `file`、`url`、`reason`。每处理完一个文件就原子更新失败清单；再次通过的文件自动移除，用 `--only` 处理部分文件时保留其他失败项。全部解决后失败清单只剩表头，不生成单独的 repair 清单、锁文件或逐样本日志。重跑仍会重新计算 MD5，已正确的文件不会重复下载或修补。退出码：`0` 本次选取的文件全部通过，`1` 本次仍有文件需要处理，`2` 参数、清单或运行冲突错误，`130` 用户中断。

## 代码风格

- 参考 `0f/` 中手写代码的风格：运算符两侧留空格，使用 tab 缩进。不同功能模块之间留两行空行，再用 `# 🚩 模块说明` 分隔。
- 变量、函数和文件使用直观、能说明用途的名称；避免过短的缩写，例如使用 `add_rsid.py`、`hf_download.py`。
- 只保留一个最终版本。不要保留旧代码、兼容别名、转发文件、重复实现或永久备份。确有整理需要时，暂存文件放 `/tmp`。
- 删除函数前检查整个脚本目录的调用关系，包括 shell 调用、动态加载、配置和命令入口；确认没有用途后再删除。改名或合并时同步修改调用方、路径和说明。
- `github_copy.sh le8,grid` 直接接收项目列表，不再使用 `scan`、`sync` 子命令；也支持空格分隔多个项目。先扫描并生成 `github_files.lst`，显示警告和汇总后等待确认，输入 `y`/`yes` 才复制；回车、`n` 或输入结束则取消复制并保留清单。扫描仅在两条同时满足时判定 UKB ID：列名包含 `id`（忽略大小写，例如 `id`、`eid`、`IID`、`participant_id`），且该列所有数据行均为七位整数。空表不命中，缺失、空白或非七位值会使该列不符合条件；数值导出的 `.0`、科学计数法按等值整数识别。无表头数据、零散七位数字、坐标、GO 编号、元数据及文件名不单独触发 ID 排除。用户已确认 `gu` 项目没有 UKB 个体数据，1KG 个体数据可以公开，因此仅对 `gu` 的正式图表跳过内容及 ID 检查；不能仅凭目录名包含 `1kg` 就豁免。无法检查的内容仍排除。
- UKB 范围内的 XLSX 会检查所有 worksheet（包括隐藏表）和附加表格。识别 `le8-workbook-v1`、`analysis-tables-v1` 的内嵌导出，核对清单、大小（若提供）及 SHA-256 后，按原始 CSV/TSV/JSON/JSONL/文本格式检查内容（支持 gzip）；不能因 `.bin` 后缀直接排除，也不能跳过检查。未知二进制附件仍排除。Excel 的 `~$` 临时锁文件不复制。
- 扫描按顺序先排除大于 50 MB（50,000,000 字节，等于上限可纳入）的文件，再排除 RDS，最后检查 UKB ID 列。某 XLSX worksheet 命中时，屏幕日志打印 `WARNING`，列出文件路径和 worksheet 名称，不打印个体值；该工作簿不进入复制清单。扫描全部结束后在同一次运行中等待用户确认，再直接按清单复制。
- RDS 可能含有个体 ID，现有扫描器不能可靠检查其内容，因此所有项目（包括 `gu`）的 `.rds`、`.rds.gz` 均在扫描阶段排除。
- `gu` 的扫描清单只保留正式 XLSX/XLSM、图片、SVG 和命名 PDF，不收集展开的 BED/TSV、原生输入等文件；跳过 `raw`、`mask`/`masks`、`cache`、`extract`、`infer`、`inputs`、`tmp`/`temp`、`normalize`、`storage-migration` 目录。单个非 UKB 项目达到 10,000 个候选文件时扫描直接报错，保留原清单，不截断清单。
- 扫描继续排除 `.log`、`.done`（含 gzip 压缩形式）以及没有归属名称的 `plots.pdf`、`Rplots.pdf`。分析源目录内用于复用的完成标记不因扫描或复制而删除。
- 确认后的复制阶段严格按 `github_files.lst` 复制，不再次筛选清单，不复查大小、扩展名或 UKB 内容；只保留路径安全检查。覆盖清单内同名目标文件，未列出的目标文件不变。源文件已经不存在时跳过并汇总，不改写清单。复制暂存和回滚备份均在 `/tmp`，复制失败时恢复已改动目标文件。

## 目录、归属和合并

- 按模块或方法组织代码，尽量每个模块只保留一个 `.R`、一个 `.py`、一个 `.sh` 文件。通用逻辑集中复用；确有独立运行或依赖要求时，可保留必要的少量文件。
- LE8 的公共调度、读写和临时续跑核验集中在 `f/0.common.py`，队列输入及运行衔接在 `f/0.host.py`，独立资源调度保留 `f/0.resources.py`。ABM、分子偏离和汇总分别由 `abm.py`、`deviation.py`、`final.py` 实现。`config/abm.json` 与 `config/deviation_{prot,met}.json` 是默认运行参数，整理代码时保持其内容不变。已完成 C1 的核验记录只放 `/tmp/le8-completed/`，文件状态改变或记录不存在时重新执行原有结果检查。
- 项目 `f/` 中的前期处理和公用代码使用 `0.` 前缀，例如 `0.common.R`。其他模块使用对应前缀，例如 `c1.`、`c2.`；不用 `0_`、`c1_` 作为代码文件的模块前缀。
- `0f/` 只放确实常用、跨项目复用的工具；目录名已有 `0`，其中的文件不再加 `0.`。项目专用功能留在对应项目，例如 GWAS 专用功能放 `gwas/`。
- GRID 的依赖统一由 `grid/environment.yml` 定义，使用 `grid/install_grid.sh` 安装；不另保留 twin 的重复依赖清单。

## LE8 可选输入

`le8/config/` 只保留默认运行读取的 `abm.json`、`deviation_prot.json` 和 `deviation_met.json`。多路线分析由各模块代码实现，不另存未被程序读取的路线说明 JSON。

- C4 自定义时钟通过环境变量 `C4_CLOCK_MANIFEST` 指定 TSV，列为 `clock`、`source`、`variable`。同一 `clock` 只能使用一种来源：`all`、`img` 或 `omics`，且至少有三个实际可用变量。例如肾脏时钟可列出 `all` 中的 `creatinine`、`cystatin_c`、`eGFR`，前提是这些列确实存在。不设置该变量时，仍按原流程处理全局组学、可用的表型系统及脑影像时钟。
- 显式 PGS 映射由 `--pgs-mapping` 指定，列为 `pgs_column`、`feature`，例如 `L_VLDL_TG.pct.SCORE_SUM` 对应 `L_VLDL_TG.pct`。
- 客观终点数据由 `--objective-file` 指定。`known_cad_at_index` 表示采样时已知 CAD（0=无，1=有，缺失=未知）；`index_date` 为分子采样日期；`objective_date` 为独立判读的冠脉检查日期；`objective_cad` 为按预先定义标准判定的 0/1 客观病变终点；`objective_modality` 使用 `CCTA` 或 `INVASIVE_CORONARY_ANGIOGRAPHY`，`CAC_CT` 仅用于独立钙化任务。两次日期不同时需要 `intervening_cad` 和 `intervening_treatment`，0 分别表示期间无 CAD 事件、无相关新治疗或治疗改变。

## 分析结果

- 正式输出只保留重要的分析结果、必要的质量控制结果，以及确实可复用的数据或模型对象。不要输出零散的运行参数、旧代码版本记录、输入路径清单等管理文件。
- 不生成或保留 `analysis_options.rds`、`revision_manifest.rds`、`analysis_manifest.rds`、`run_options.json`、`fit_options.json`，也不要换成同用途的其他格式文件。程序需要的最少缓存一致性检查信息放在相应结果对象内部。
- 每张正式 `.png` 配一个同目录、同名的 `.xlsx`，例如 `c3.Fig2.regional_top_loci.png` 与 `c3.Fig2.regional_top_loci.xlsx`。工作簿只放该图对应的分析结果；不同 panel 的结果可放不同 worksheet。不要把整个模块的所有表塞进一个巨型工作簿。
- 同一主题分页时，各图的工作簿应对应各自展示的数据范围。没有对应图片、但确有独立用途的分析结果，按方法或主题另存少量工作簿；不要为了减少文件数而制造难以浏览的大量 worksheet。
- 工作簿按分析内容组织，不按“每 N 张表”机械拆分，也不使用含义不明的 `.2.xlsx`、`.3.xlsx` 分卷。先把结果归到对应图片或分析主题，再合并同结构的表，用结果类型列区分；确需拆分时使用具体主题名称，例如 `c1.pgs.comparison.xlsx`、`c1.pgs.decomposition.xlsx`、`c1.pgs.temporal.xlsx`。时间窗口、landmark 及各自对比放在一起，富集结果不混入关联分析工作簿。图片重新编号后，结果表按主题继续归属该图。
- Worksheet 使用简短、直观的结果名称。删除空表、无结果的占位表、重复表、过期的跨模块注释，以及运行参数、输入路径、版本和文件清单等管理性 worksheet。例如 C3 不保留无用的 `CIGMA_annotation`。必要的分析质量控制结果可以保留。
- 不在 XLSX 内重复嵌入旧 XLSX，也不同时保留同内容的 CSV、XLSX、汇总 RDS。区域／SNP 等完整大数据若已保存在可复用的分析对象中，图表工作簿只导出该图涉及的数据，不复制整个底层数据集。
- XLSX 必须能被 Microsoft Excel 正常打开，不能只验证 R/Python 能读取。ZIP 包不写入独立目录项；内部精确数据导出使用明确的内容类型，不能继承 `.bin` 的打印设置类型。修改打包逻辑后，用 Excel 正常打开模式验证代表性工作簿，并核对单元格内容与内部数据未改变。
- 图片可保留多个 `.png`；同一主题尽量用多个 panel 组织。已有 PNG 时不再生成同内容 PDF。只用于拼图的中间 panel 放 `/tmp`，正式目录保留合成图。
- 图片文件名应说明归属和内容，不生成或保留通用名称的 `plots.pdf`、`Rplots.pdf`；R 绘图应显式指定输出文件，避免默认图形设备产生这些文件。
- 图号只出现一次，采用“模块＋图号＋主题”，例如 `c3.Fig5.pgs_coloc_triangulation.png`；不生成 `c3.Fig1.c3.Fig7...` 这样的重复前缀。同一主题下内容完全相同的图片只保留一份，改名时同步更新工作簿、网页和 Shiny 索引。
- 带 UKB 等个体 ID 的数据保存为 `.rds`，不能混入汇总工作簿或公开查看界面。文件名按数据内容确定，如 `test_individuals.rds`、`individual_explanations.rds`、`review_split.rds`，不要统一命名为 `individual.rds`；嵌套个体解释也不用 JSONL 留在正式目录。
- 可复用的模型或分阶段分析对象也可以保留 RDS。文件名简短、固定、说明用途，例如 `c2.dandelion.rds`、`c2.instruments.rds`、`c2.mr.rds`，不把日期、模式和版本串堆进文件名。
- GRID 的个体 PRS-CSx posterior 使用 `1csx.posterior.rds`，不以 `individual_posterior.tsv.gz` 留在正式结果目录。迁移个体文件时更新所有生产方、读取方和报告链接，并核对行数、ID、数值和类型。
- GRID 的评分分别使用 `1csx.scores.rds`、`2disco.scores.rds`，已有 COJO 评分继续读取 `Yeval.cojo.rds`；统一评估入口为 `grid.sh final`，图表使用 `final.<主题>.png/.xlsx`，个体结果只保存为 RDS。已有 MCMC 后验 HDF5 和 SNP 模型权重保留原生格式供复用，评分交换文件及缓存进入 `/tmp`。
- GU 的个体／单倍型及原生拟合结果按方法、数据集和位点／染色体保存为 `phyml.haplotypes.rds`、`ibdmix.tracts.rds`、`trace.segments.rds` 或 `as3.tracts.rds`；综合个体结果保存为 `final/gu.results.rds`。SQLite 仅作为 `/tmp` 中的查看缓存；`gu.validation.rds` 保留需供 Shiny 读取的个体验证数据。正式结果不重复复制到 `final/normalize`。
- 旧图片若未保存某项作图统计量，应明确说明缺失，不能把其他汇总统计量当作原图数据。是否重新拟合须单独获得授权；允许重算时应同步更新 PNG 和 XLSX。
- 为减少文件数量，可以迁移现有结果，并同步修改读写路径。迁移后验证数据完整、数值不变、下游可读取，再删除被替代的文件。
- LE8 的综合结果直接放在 `final/`，不另建 `overview/`，不生成结果目录的 README。C4 固定预算验证并入 `c4_connect/`，按对应图片分配工作簿；验证对象使用 `c4.validation.rds`，避免覆盖连接分析的 `c4.res.rds`。

## 临时文件与复用

- 所有临时数据、交换用 CSV/TSV、测试及验证脚本、临时备份、工具缓存均放 `/tmp`。`tests/`、`__pycache__/`、`.ruff_cache/` 等不留在 `/mnt/d/scripts/` 或正式结果目录。
- `/home/huangj` 不用作分析输出、临时工作目录或工具缓存目录；正式结果写入各项目指定位置，测试和缓存写入 `/tmp`。主目录中的现有软件和运行环境可继续读取。
- 独立回归测试、性能比较脚本及其日志不留在正式代码目录；临时需要时放 `/tmp`。主流程实际调用的数据校验、统计检验、训练和预测代码按功能命名并保留，不按文件名是否含 `test` 机械删除。
- 清理正式结果时，先确认哪些结果可以复用、哪些可以从保留的数据重新生成。可重新生成的旧 PNG、CSV 可以删除；不要删除唯一的分析结果或仍需复用的数据。
- 不因整理文件而重新拟合模型、改变分析方法、修改估计值或伪造缓存信息。
- LE8 通过 `le8.sh` 在 `/tmp` 展开交换数据、执行分析和重绘，成功后发布结果。当前结果工作簿同时保存可复用的汇总导出，因此应保留；个体数据仍单独保存为按内容命名的 RDS。
- 结果保留“疾病／数据层／模块”的必要归属，减少方法名、运行版本、导入副本等重复目录层级。同一份 ABM 结果直接读取，不另建 `imported_*` 副本。中间扫描、GPU 格式转换和自动输入表使用 `/tmp/le8-cache/`。
- 分享 Shiny 使用专门的汇总查看包，按已登记的表格和图片收集依赖，保留相对路径；不直接发送整个分析目录，不打包个体 RDS、模型权重、原始数据或临时缓存。
- LE8 工作簿可以正常查看。Excel 等软件重新保存可能移除程序用于恢复精确数据的内部内容；需要编辑展示版时另存到其他目录。具体入口和选项见 `le8/le8.sh --help`。
