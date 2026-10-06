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

日志、锁、生成的 mask、预处理 VCF、命令列表和计算工作区放在 `/tmp`。每个位点的参数和来源记录保存在原生归档内。XLSX 使用多个工作表整合结果，并保留精确表格导出供代码读取；超出 Excel 行数上限时拆分工作表，过长文本按顺序拆到附加工作表，不截断内容。发布使用暂存文件校验后替换。

`GU_ANALYSIS_ROOT` 与 `GU_PUBLISHED_ROOT` 使用同一永久目录；结果路径不能指向 `/tmp`、`/var/tmp`、`/dev/shm` 或 `/run`。跨项目的 XLSX/RDS 读写共用 `../0f/results.R` 和 `results.py`。格式整理不重新拟合模型，也不改变科学结果或伪造分析完成状态。

代码格式遵循公用 R 代码：运算符两侧空格、tab 缩进、使用 `# 🚩` 划分模块。第三方源码保留上游格式。

本次整理核对了全部已保存的科学结果和 258 棵完整树，未重新推断树或调用 IBDmix／TRACE。TRACE 目录当前保留的是可继续分析的提取数据；整理不会把未完成的阶段标成已完成。
