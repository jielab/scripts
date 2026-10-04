# LE8

统一入口为 `le8.sh`。`f/` 中公用文件以 `0.` 开头，分析文件以 `c1.` 至 `c5.` 标明归属。全项目整理约定统一见 [根目录 README](../README.md)。

| 模块 | 文件 |
|---|---|
| 公用、调度、索引、验证和导出缓存 | `0.common.R`、`0.common.py`、`0.engine.sh` |
| C1 相关与 ABM | `c1.correlate.R`、`c1.abm.py` |
| C2 因果分析、MR-link-2、cis 索引 | `c2.cause.R`、`c2.cause.sh`、`c2.mr_link2.py` |
| C3 共定位 | `c3.coloc.R`、`c3.coloc_GPU.py`、`c3.coloc_GPU.sh` |
| C4 连接分析 | `c4.connect.R` |
| C5 细胞分析 | `c5.cellulation.py` |
| 最终汇总、表格、科学问题报告 | `final.R`、`final.py` |

C1 使用 ABM（agent-based modeling）命名。模型运行、注释及图形代码集中在 `c1.abm.py`，注释和重绘分别使用 `annotations`、`figures` 子命令。新模型统一使用 `c1_abm` 序列化身份；不注册已经合并掉的旧 Python 模块名。

`final.py` 合并原先分散的最终报告、表格和科学问题处理。Shiny 入口为 `shiny/app.R`。正式执行统一使用 `le8.sh`；辅助文件是运行目录内的内部步骤。

## 结果文件

- 每张 PNG 配同名 XLSX，只放对应图的分析结果；同类 panel 用 panel 列区分，尽量合并为少量 worksheet。例如 `c3.Fig2.regional_top_loci.png` 对应 `c3.Fig2.regional_top_loci.xlsx`。重要的未绘图结果按主题另外保存。
- 未绘图结果按分析内容组织，不按表数量分卷：C1 使用 `c1.association.xlsx`，PGS 配对比较、分解、时间分析分别为 `c1.pgs.comparison.xlsx`、`c1.pgs.decomposition.xlsx`、`c1.pgs.temporal.xlsx`；C4 验证按表现、面板、模型和遗传证据区分。Final 和 Shiny 同样按主题组织，不再生成 `.2.xlsx`、`.3.xlsx`。
- 工作簿内部同时保存原始汇总导出及校验信息，供重绘、缓存复用和 Shiny 读取；不再另存通用汇总 RDS。完整区域／SNP 结果留在可复用分析 RDS，图表只导出实际展示的范围，不重复嵌入旧工作簿。程序可以在 `/tmp` 完整恢复这些导出。
- 工作簿是可复用结果，需保留。Excel 等编辑器重新保存时可能移除内部来源记录；需要编辑展示版时，另存到其他目录。
- 带个体 ID 的表按具体数据名称保存为 `.rds`，例如 `test_individuals.rds`、`individual_explanations.rds`、`c4.focus.roles.rds`，不进入汇总工作簿或 Shiny。嵌套的个体解释和 attention 数组也保存在 RDS 内，不另留 JSONL、NPZ 文件。
- 不单独输出分析参数、代码版本或输入路径清单，不生成 `analysis_options`、`revision_manifest`、`analysis_manifest` 文件。缓存复用所需的最少一致性信息保存在相应分析结果对象内部。
- 可复用的 R 模型和分阶段结果保留为 RDS。名称简短且固定，例如 `c2.dandelion.rds`、`c2.instruments.rds`、`c2.mr.rds`、`c2.reverse_mr.rds`；版本信息放在对象内部。
- `le8.sh` 在 `/tmp/le8-run-*` 中展开交换表、执行分析和重绘。成功后发布结果工作簿、PNG 和个体/模型 RDS；运行失败时保留 `/tmp` 中的诊断现场。

## 目录与阅读顺序

| 目录 | 内容 |
|---|---|
| `<疾病>/<prot或met>/c1_correlate/` | C1 工作簿、PNG 和可复用分析对象；`abm_reference/`、`abm_tabicl/` 各保存一个 ABM 方法的结果 |
| `<疾病>/<数据层>/c2_cause/` | 按图配对的 XLSX、PNG 和命名明确的 RDS |
| `<疾病>/<数据层>/c3_coloc/` | 共定位汇总、PNG 和可复用结果；GPU 格式转换和逐区域计算缓存在 `/tmp` |
| `<疾病>/<数据层>/c4_connect/` | LE8 连接、代理、交互、非线性及固定预算验证结果；每张 PNG 配同名 XLSX，验证拟合保存在 `c4.validation.rds`；Yin、YinYang 区分代理发现人群 |
| `<疾病>/<数据层>/c5_cellulation/` | 细胞注释汇总和 PNG；自动构建的输入表在 `/tmp` |
| `final/` | Fig1–8、补充图及各自的 XLSX、`report.html`、`index.html` 和图注；研究问题综合结果也在本层，不另建 `overview/`、`tables/`、`panels/` 或结果 README |
| `shiny/` | 按主题拆分的索引工作簿，含查看器索引和综合视图；表格详情、图像仍通过相对路径引用各模块结果 |

先看 `final/report.html` 或 Shiny 总览，再按疾病、数据层和模块查阅结果。工作簿使用简短的结果页名，不保留空表、重复表和无用的管理页。已有 PNG 时不生成同名 PDF；中间 panel 图、PGS 扫描缓存、MR-link-2 工作文件、GPU 转换文件和日志均留在 `/tmp`。ABM 的输入副本、逐折神经网络检查点也放在 `/tmp/le8-cache/`；正式目录保留自包含的最终 `model_bundle.joblib` 和必要的个体 RDS。ABM 的冻结状态、校准信息等复用所需信息收进方法工作簿内部，不再散落为 JSON 文件。不因整理重新拟合模型。

## 分享 Shiny

只发送 `shiny/` 和 `final/` 不足以打开全部详情与图片。生成完整的汇总查看包：

```sh
./le8.sh share --share-out /mnt/d/analysis/le8-share.zip
```

发送这个 ZIP 即可。包内只有两份查看代码（`shiny/app.R`、`f/0.common.R`）、说明和 `results/` 中所需的汇总工作簿与 PNG；不包含个体 RDS、训练权重或原始 UKB 数据。创建时会核对来源校验值、检查个体 ID 列，并在独立目录验证 Shiny 读取。

接收方安装 R，以及 `shiny`、`DT`、`data.table`、`ggplot2`、`digest`、`openxlsx`、`jsonlite` 包。解压后在包的根目录运行：

```sh
Rscript shiny/app.R
```

浏览器打开 `http://127.0.0.1:3839`。查看包不需要 Python、WSL 或原作者的绝对路径；请保留包内相对目录结构。Windows 的 RStudio 也可切换工作目录到包内 `shiny/`，然后运行 `source('app.R')`。

自己的原结果仍可这样打开：

```sh
./le8.sh shiny --no-reindex --analysis-root /mnt/d/analysis/le8
```

运行 `./le8.sh`，或使用 `./le8.sh --prepare-only` 生成结果后退出、不启动 Shiny。测试、临时验证、备份和 Python 缓存均放在 `/tmp`，不保留在代码目录。
