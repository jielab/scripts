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

## 2026-10-05 方法更新

- C1 保留常规关联及 PGS 分析，reference ABM 更新为 selective Yang pipeline：开发集内选择、校准和冻结门控，验证集评估覆盖率与风险增益。默认 C1–C5 流程包含 reference ABM，无需添加 `--run-abm --abm-backend reference`；TabICLv2 仍可通过 `--abm-backend tabicl` 或 `both` 选择。显式选择 `c1_abm` 也会执行训练；已有匹配结果按缓存规则复用。
- C2 将交叉拟合的遗传预测部分 G 与残差 R 同时放入调整模型，直接检验二者差异；残差不解释为纯环境作用。MR 使用边际效应、独立工具和坐标/等位基因身份；冲突重复变异排除并记录质量控制。MR-link-2 按完整预定检验家族校正。DANDELION 原生结果、全局多重校正、逐位点剔除结果分别展示；真实 WES 与 GWAS 基因适配证据分开。冻结的 Yang 状态投影及近远期风险属于补充分析。
- C3 固定同一组、有序且方向一致的变异供 CPU/GPU 比较，稳定计算 H3，并保留先验敏感性与输入条件。未通过敏感性检查的 H4 不标记为稳健共定位。SuSiE 的可信集与各 MR 工具逐一对应，只有部分工具信号获得支持时标记 partial_signal_support，不升级整个 MR 汇总结果。SuSiE 可选分支要求有符号 LD、参考样本量、版本和祖源元数据；缺失条件时给出不可用状态。
- C4 分开基本调整的总连接与条件特异性，允许一项分子连接多个 LE8 域。YS 按无疾病结局参与的连接强度和冗余排序；固定预算下比较 YS、NS、YSplus、原始分子及嵌套交叉拟合的概念表示。保存临床与分子线性预测值贡献、惩罚路径和非零系数，避免把零分子贡献解释为模型优势。路径乘积及其有符号比例只作描述；bootstrap 不足时不输出推断。模块稳定性是固定候选分子集合下的重聚类稳定性。
- C5 以全部已测分子为背景，保留未知注释，使用等预算且排他基因的细胞对比。原生 CIGMA、外部导入结果和不可用状态严格区分；只有符合输入条件时才运行对应分析，不从缺失证据推导衰老或细胞谱系结论。
- Final/Shiny 增加上述方法、有效性状态、验证比较及下载表。缺失分析明确显示缺失，汇总报告不把不可用结果补成阳性证据。

有家系或相关个体时，可用 `LE8_GROUP_FILE` 指定含 `eid,group` 的文件，或用 `LE8_GROUP_COLUMN` 指定个体表中的分组列；交叉验证、开发/验证划分和相关 bootstrap 按组进行。未提供时以个体为组，程序不会猜测家系。疾病 PRS 是可选临床基线输入，需显式提供 `C4_DISEASE_PRS_COLUMN` 或 `C4_DISEASE_PRS_FILE`，同时设置 `C4_DISEASE_PRS_TRAIT` 与 `--Y` 一致；不会拿单分子 PGS 替代疾病 PRS。

分析缓存检查代码、输入文件及相关设置；输入或方法变化时重新计算。同一结果根目录下，不同 `Y/biom` 的命令可以同时计算；相同 `Y/biom` 的命令排队，避免重写同一组结果。读取结果快照和发布结果时使用短期共享目录锁，每个分析任务只发布自己的疾病/数据层。`final,shiny` 等汇总任务等待正在运行的分析结束，再合并已有疾病/数据层；已有 Shiny 服务时不重复启动。测试请使用独立的 `--analysis-root`，不要将缩小特征数或 bootstrap 次数的测试结果发布到正式目录。

## 分开分析，最后汇总

可以在不同终端运行以下命令，使用同一个默认结果目录：

```sh
./le8.sh --Y cvd_cad --biom prot
./le8.sh --Y cvd_cad --biom met
./le8.sh --Y ra --biom prot
```

全部分析完成后，统一汇总多个疾病和数据层：

```sh
./le8.sh final,shiny --Y cvd_cad,ra --biom prot,met
```

默认流程为 `c1_correlate,c1_abm,c2_cause,c3_coloc,c4_connect,c4_panel_validation,c5_cellulation`。`final,shiny` 只汇总已有结果，不自动训练 ABM；显式模块列表只运行所选模块及其必要前置分析。需要跳过默认 ABM 时使用 `--skip-abm`，需要在其他模块组合中加入 ABM 时仍可使用 `--run-abm`。

C1 PGS 按特征原子保存检查点，终端显示已完成数、复用数和估计剩余时间。相同输入、模型和参数下重新执行原命令即可续跑；输入或方法变化不会混用检查点，`--replace TRUE` 会重新计算。每个特征仍计算完整的六类模型，完成整个特征家族后统一计算 FDR。`--cores N` 现在也控制 PGS 特征并行；默认仍为 1，以控制完整 PGS 矩阵与模型数据的内存开销，多终端并发时需合计各任务的内存。已完成的 PWAS/MWAS 扫描同样保留在 `/tmp/le8-cache/` 供重启复用，检查原始输入、协变量、数值方法及相关设置；修改调度或报告代码不会单独触发这些扫描重算。清理 `/tmp` 或系统重启后，未发布的检查点可能丢失。

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

运行 `./le8.sh` 默认只执行 C1–C5（包含 reference ABM 和 C4 面板验证），完成后退出。汇总与查看单独使用 `./le8.sh final,shiny`；只生成汇总而不启动服务，可加 `--prepare-only`。测试、临时验证、备份和 Python 缓存均放在 `/tmp`，不保留在代码目录。

测试覆盖包括安装后的 C1/C5 检查、七类遗传/反向疾病模拟、家系隔离与嵌套选择防泄漏、模型重载、原生 MR-link-2、DANDELION、CPU/GPU 后验一致性、SuSiE 多信号和输入不匹配。实际 UKB 验证使用独立目录中的缩小样本/特征子集；不等于全队列结果已经验证。SuSiE 的真实分析仍需提供与各 GWAS 的 `ancestry` 一致的 LD 元数据，未知效应尺度或缺失强信号不会升级为完整 MR–coloc 支持。Final 默认展示 10-assay 预算，若未运行该预算则展示最近的已配置预算，选择不依据预测结果。

可选扩展的边界：本版未启用组织 eQTL/sQTL 扩展或 SuSiE 信号 LBF 的 GPU 批处理；GPU 主线比较的是相同区域边际 BF。风险尺度因果中介和年龄断点搜索也未启用，现有输出分别为描述性路径乘积和基线年龄平滑曲线。原生 CIGMA 仍需匹配的单细胞表达和基因型输入。
