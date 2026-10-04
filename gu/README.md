# GU 代码目录

公共入口保持 `gu.sh`、`arg.sh`、`install.sh`；内部代码按方法合并，准备和公用工具以 `0.` 开头。

| 模块 | shell | Python / R |
|---|---|---|
| 公共配置、区域和表格工具 | `f/0.common.sh` | `f/0.common.py` |
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

结果目录和科学定义不变。代码签名仍用于检查缓存；整理前的缓存可能需要重新计算，不能仅凭文件存在就跳过。旧日志是历史记录，不修改其命令内容。`f/phyml.py cache adopt` 保留原有受验证的缓存接纳机制。

代码格式遵循公用 R 代码：运算符两侧空格、tab 缩进、使用 `# 🚩` 划分模块。第三方源码保留上游格式。
