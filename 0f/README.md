# 公用工具

这里只存跨项目工具或可直接使用的独立工具。目录名已有 `0`，文件名不再加数字前缀，使用能说明用途的名称；调用方同步使用这些路径。

| 文件 | 用途 |
|---|---|
| `phenotype.R` / `phenotype.sh` | 表型、基因型、GWAS 公用数据处理与配置 |
| `association.R` | 关联模型 |
| `prediction.R` | 预测模型 |
| `results.R` / `results.py` | LE8、GRID、GU 的结果工作簿、可复用原始对象与临时交换数据 |
| `abm.py` | GRID 使用的精确分块近邻检索、Transformer 基础层、家系 ESS、掩码重建损失及梯度审计 |
| `plotting.R` | 分析图表与主题 |
| `manhattan_plot.R` | GWAS Manhattan 图 |
| `configure_ml.R` | R/reticulate 与机器学习环境 |
| `memory_cap.sh` | GRID、GU 共用的进程组内存限制 |
| `console.sh` / `console_run.py` / `process_lifetime.py` | 统一命令行日志、前台任务与子进程生命周期 |
| `join_file.py` | 按键合并文件 |
| `add_rsid.py` | rsID 补充 |
| `hf_download.py` | Hugging Face 下载 |
| `img2pdf.sh` | 图像转 PDF |

GWAS 格式化及 BGZF/tabix 索引实现位于 `../gwas/f/format.f.sh`，LE8 通过 `source .../format.f.sh --index` 只加载索引函数。清理依据覆盖 `/mnt/d/scripts` 下源代码、命令模板及配置中的引用；直接入口和动态调用使用的名称保留。已删除 68 个无用函数，剩余函数/类的引用复查通过。运行中仍需各分析项目自己的依赖环境。

结果格式按用途决定：查看和交付的表格使用 XLSX，含 ID 的表也可以写入工作表；RDS 用于后续代码需要复用的大型原始数据或拟合对象。XLSX 保留精确表格源数据，按 Excel 容量拆分工作表，长文本拆分到附加工作表；转换校验成功后才能删除旧 RDS。日志和临时交换文件放在 `/tmp`。

GRID 的检索/注意力数值核心维护在 `abm.py`，由 `f/attention.py` 加载。
LE8 当前使用 `../le8/f/abm.py` 的 twin confidence 实现和 `../le8/f/deviation.py` 的遗传校正分子偏离实现。两项目各自维护输入、结局定义和训练/验证策略。

默认前台运行。普通日志、`SCRIPT_VERBOSE=1`、Shiny 和原始输出模式使用同一进程清理机制。
Ctrl+C、关闭终端或启动器退出会停止整棵任务进程树；启动器被 SIGKILL 时，Linux 父进程死亡信号触发子监督进程清理。
内部并行 worker 不等于后台任务：入口持续等待它们，并清理被中间进程遗留的 worker。
只有用户明确指定 `--foreground FALSE`，或自行使用 `nohup`/后台作业命令时，才允许任务脱离交互终端。
不要用 `SCRIPT_VERBOSE`、重定向日志或创建 systemd service 来隐式切换后台运行；内存限制使用同步 `systemd-run --scope`。
