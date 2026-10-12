# 公用工具

这里只存跨项目工具或可直接使用的独立工具。目录名已有 `0`，文件名不再加数字前缀，使用能说明用途的名称；调用方同步使用这些路径。

| 文件 | 用途 |
|---|---|
| `phenotype.R` / `phenotype.sh` | 表型、基因型、GWAS 公用数据处理与配置 |
| `association.R` | 关联模型 |
| `prediction.R` | 预测模型 |
| `results.R` / `results.py` | LE8、GRID、GU 的结果工作簿、可复用原始对象与临时交换数据 |
| `abm.py` | LE8、GRID 共用的精确分块近邻检索、Transformer 基础层、家系 ESS、掩码重建损失及梯度审计 |
| `plotting.R` | 分析图表与主题 |
| `manhattan_plot.R` | GWAS Manhattan 图 |
| `configure_ml.R` | R/reticulate 与机器学习环境 |
| `memory_cap.sh` | GRID、GU 共用的进程组内存限制 |
| `console.sh` / `console_run.py` | 统一命令行日志 |
| `join_file.py` | 按键合并文件 |
| `add_rsid.py` | rsID 补充 |
| `hf_download.py` | Hugging Face 下载 |
| `img2pdf.sh` | 图像转 PDF |

GWAS 格式化及 BGZF/tabix 索引实现位于 `../gwas/f/format.f.sh`，LE8 通过 `source .../format.f.sh --index` 只加载索引函数。清理依据覆盖 `/mnt/d/scripts` 下源代码、命令模板及配置中的引用；直接入口和动态调用使用的名称保留。已删除 68 个无用函数，剩余函数/类的引用复查通过。运行中仍需各分析项目自己的依赖环境。

结果格式按用途决定：查看和交付的表格使用 XLSX，含 ID 的表也可以写入工作表；RDS 用于后续代码需要复用的大型原始数据或拟合对象。XLSX 保留精确表格源数据，按 Excel 容量拆分工作表，长文本拆分到附加工作表；转换校验成功后才能删除旧 RDS。日志和临时交换文件放在 `/tmp`。

ABM 的数值核心只维护在 `abm.py`，LE8 的 `f/c1.abm.py` 和 GRID 的
`f/attention.py` 加载它。各项目保留自己的输入、结局定义和训练/验证策略。
LE8 使用 1e-10 距离并列容差，GRID 使用精确距离排序；两者保留 FP64 检索、
自身/家系/训练折排除，不作近似检索或自动 CPU 回退。结果签名包含共用文件的哈希。

2026-10-12 优化：身份编码常驻设备，局部 top-k 后只排序小候选集；LE8 默认检索块
256×8192（原 32×2048），固定 global/tail 邻域在一次拟合内复用，每轮一次训练候选
检索，温度网格共用编码和检索，参考键/标签在 GPU 复用。梯度审计每轮剔除已证实有
梯度的参数，避免每步重复检查。训练批量、折数、重复次数和早停设定保持原值。
缓存不进入序列化模型，不在不同队列或拟合之间共享。LE8 的 ABM 按结局/组学分别发布。

共用数值核心已通过 CPU/CUDA 数值一致性、两项目短程训练、温度缓存等价、
冻结状态、模型重载和 GRID 测试标签隔离验证。独立测试脚本已移出正式代码目录；
临时验证代码和运行产物统一放在 `/tmp`。
