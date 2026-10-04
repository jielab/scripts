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

- `phyml.haplotypes.rds`：单倍型、个体携带关系、保存的树及可复用算法结果。
- `ibdmix.tracts.rds`、`trace.segments.rds`、`as3.tracts.rds`：各方法的个体片段和可复用结果。
- `final/gu.results.rds`：综合关系数据表；使用 `readRDS()` 可按表名读取。它含个体数据，不属于公开分享包。
- `final/gu.validation.rds`：个体携带者交叉验证与 Shiny 报告数据。
- `final/gu.*.xlsx`：独立的汇总证据与质量控制结果。轨迹结果按染色体划分，避免超过 Excel 行数限制。

`./gu.sh` 自动在 `/tmp/gu-cache/` 展开算法交换文件与 SQLite 查看缓存，完成后发布 RDS 和图表。Shiny 仍通过 `./gu.sh shiny` 启动，直接运行 `shiny/app.R` 也会从 RDS 准备缓存。`final/normalize` 的重复 PhyML 副本不再保留；临时查看目录通过链接读取同一份方法结果。日志、测试、参考调用缓存与运行参数交换文件均在 `/tmp`。

跨项目的 RDS／XLSX 写入由 `../0f/results.R` 和 `results.py` 共用；GU 的原生结果展开、树图数据和关系数据库结构留在本项目。整理现有结果不重新拟合模型，也不伪造缓存签名。

代码格式遵循公用 R 代码：运算符两侧空格、tab 缩进、使用 `# 🚩` 划分模块。第三方源码保留上游格式。

本次整理核对了全部已保存的科学结果和 258 棵完整树，未重新推断树或调用 IBDmix／TRACE。TRACE 目录当前保留的是可继续分析的提取数据；整理不会把未完成的阶段标成已完成。
