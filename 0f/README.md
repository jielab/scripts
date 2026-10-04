# 公用工具

这里只存跨项目工具或可直接使用的独立工具。目录名已有 `0`，文件名不再加数字前缀，使用能说明用途的名称；调用方同步使用这些路径。

| 文件 | 用途 |
|---|---|
| `phenotype.R` / `phenotype.sh` | 表型、基因型、GWAS 公用数据处理与配置 |
| `association.R` | 关联模型 |
| `prediction.R` | 预测模型 |
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
