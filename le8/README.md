# LE8 5C — disease PRS + proteome + metabolome

版本：2026-09-16。完整 LE8 模块代码；本包不含分析结果或 patch。

## 1. 安装与最简运行

将包内 `le8/` 的内容合并到现有 `scripts/le8/`，覆盖同名代码文件，保留自己的 `revision.conf.sh` 和本地路径配置。使用 `le8.sh` 运行，不需要应用补丁。

本模块继续使用原项目 `scripts/0f/` 公共函数、UKB RDS 输入以及 R/Python/conda 环境。默认 `/mnt/d` 路径沿用原项目；不同机器应核实路径。

在生成 `all.rds` 前为 CAD 增加同名字段即可：

```r
dat0$cvd_cad.pgs <- dat0$cad.6m.pgs
# 随现有流程保存 all.rds；不必删除旧的 cad.6m.pgs 或 cad.pgs。
```

其他疾病使用完全相同的规则，例如 `ra.pgs`、`amr.pgs`。必须是该疾病的遗传风险分数，不是某个蛋白或代谢物的 PGS。该命名不会自动保证发现样本独立。

```bash
cd /mnt/d/scripts/le8
bash le8.sh c5_consolidate --Y cvd_cad --biom prot,met --preflight
bash le8.sh c5_consolidate --Y cvd_cad --biom prot,met --replace TRUE
```

上述命令保留并运行原单组学 C5，同时运行联合比较和全框架汇总。若这次只重跑联合 C5，可选：

```bash
C5_RUN_REFERENCE=FALSE bash le8.sh c5_consolidate --Y cvd_cad --biom prot,met --replace TRUE
```

完整框架可省略模块参数：`bash le8.sh --Y cvd_cad --biom prot,met`。C2/C3 需要相应 GWAS、QTL、LD 等输入；没有输入的分析不会凭空生成。`--replace TRUE` 会重算相应模块，旧结果需要留存时请先保留副本。

日常运行继续使用短命令 `./le8.sh --Y ra --biom prot,met`（或换成 `cvd_cad`）；默认使用 1 个 worker、32 GiB RAM 上限，不必重复填写。联合 C5 已有完整正式结果且汇总设置一致时，会在读取原始组学前直接复用，即使折缓存已经清理也不会因此重新训练。修改输入或模型设置要重算时，显式使用 `--replace TRUE` 或新的输出目录。其余模块按各自缓存规则复用或生成输出。

`C5_JOINT_SUMMARY_ONLY=TRUE` 仅用于所有折已保存后的故障恢复，不是日常运行必需选项。中途停止、尚未完成所有折时仍用普通短命令；只会复用与当前签名匹配的折，代码/输入/设定变化导致签名不同的旧折会重算，不会强行混合新旧版本。

如果 C5 的所有 landmark/fold 已保存，但在最后汇总时因内存不足退出（137），使用：

```bash
C5_JOINT_SUMMARY_ONLY=TRUE ./le8.sh c5_consolidate --Y cvd_cad --biom prot,met --cores 1 --memory-limit-gb 32
```

这条命令重试单组学 reference，再从 `_c5_cache/L*.fold*.rds` 恢复联合汇总。恢复会核对保存的签名、每个 landmark 的完整折数、参与者和原始分折；联合阶段不读取新的原始表型/组学，也不重训模型。缺失、损坏或混合签名的缓存会明确报错。汇总按 landmark 和 ablation 分批处理，预测导出仍保留全部行；无需提高内存上限。不要加 `--replace TRUE`，汇总成功前不要手动删除折缓存。仅处理 joint 时可再设置 `C5_RUN_REFERENCE=FALSE`。

若要同时接着运行原完整流程的后续步骤，将模块参数换成 `c5_consolidate,s1_interact,s2_nonlin,final`。恢复模式沿用缓存中的模型和队列；修改输入、协变量或训练设定后要重新拟合时，不应设置 `C5_JOINT_SUMMARY_ONLY`。汇总使用当前 `C5_BOOT`、`SEED` 和主比较设定，记录在 `c5.summary_provenance.csv`。

C5 不再生成或保存 `_c5_private/L*.rds` 中的逐模型副本。原始 Cox 对象的公式环境会连带序列化整张训练/测试组学表，使一个小模型膨胀到约 1 GB；现在只返回预测、系数、预处理参数和基线风险等必要数据。新折缓存同时保存预处理参数及基线风险，汇总后导出相应 CSV。

每折保存成功，以及从旧折缓存恢复汇总时，自动清理该折可核验的旧模型副本，记录到 `c5.temporary_cleanup.csv`。`_c5_cache` 中的折缓存仅供运行和失败恢复使用；汇总、图表、注释、完整预测导出及 `c5.res.rds` 全部成功后，自动删除本次用过且未被其他进程更新的折缓存，记录到 `c5.fold_cache_cleanup.csv`。汇总失败时保留折缓存；分折角色、完整预测导出和正式结果始终保留。相同设置下再次执行汇总恢复会复用已完成结果；折缓存清理后，修改汇总设置需要重新运行分析。

不清理未完成折、未知文件、符号链接或比折缓存更新的模型副本。汇总临时分片压缩保存，每批用完即删除，正常结束或报错时清理本次临时目录；被强制杀死的进程遗留目录在下次汇总时通过进程身份核验后清理。已运行的 R 进程不会自动加载修改后的代码。

Reference 的年龄时间尺度敏感性分析在完整样本/事件不足或 Cox 拟合失败时会输出带列定义的空表、不可估计图和原因日志，不再因 `se` 列缺失而退出。基线用药字段缺失导致的 LE8 协变量缺失仍需检查 `c0.baseline_provenance.txt` 和 `c0.baseline_rebuild_audit.csv`；不会自动把未知用药当作未用药或更换协变量。

C3 coloc 默认只处理 C2 中 `FDR_all < 0.05` 的 MR 分析，并只用这些分析实际保留的工具变量定位 loci；不再用 C1 排名或指定蛋白补足候选。`C3_MR_FDR` 可调整阈值，仍保留 `C3_MAX_FEATURES=200` 和 `C3_MAX_LOCI_PER_FEATURE=3` 的默认上限，截断会写入日志。仅 distal/trans 显著时，不会额外加入 local/cis 区域。单独运行 C3 时会自动检查/运行 C2 依赖。

C3 优先按 C2 保存的 SNP 坐标查询 QTL 的 tabix 索引，必要时只做一次缺失 SNP 的 ID 扫描；提取结果保存在分析根目录的 `.c3_qtl_cache`，可通过 `C3_QTL_CACHE_DIR` 改位置。每个特征和 locus 都记录进度、缓存命中与耗时。`c3.selected_mr.csv` 和 `c3.mr_locus_selection.csv` 分别记录入选 MR 和工具变量/locus 选择。

重跑原命令即可启用新规则，无需 `--replace TRUE`：旧的 C3 整体缓存不再直接复用，条件一致的旧单 locus 缓存可迁移到不依赖候选排序的缓存。GPU-coloc 只接收新清单，结果保存在 `c3_coloc/gpu_coloc/runs/<输入签名>/`，避免复用旧的全候选 GPU 结果；统一汇总仍写入 `c3.res.rds` 和 C3 输出表。


## 2. Disease PRS 自动识别

C5 将 `[Y].pgs` 纳入 `all.rds` 的字段读取，按 eid 对齐，且只做精确匹配。不会把 `cad.pgs`、`cad.6m.pgs` 或任意 biomarker PGS 猜测为 `cvd_cad.pgs`。

- 有有效数值的 `[Y].pgs`：自动加入 PRS 比较。
- 字段缺失、非数值、没有变异或共同事件人群不足：写明原因，继续不含 PRS 的比较。
- PRS 缺失不插补为零。PRS 可用时，含/不含 PRS 的全部主比较限制在同一个具有 PRS 的共同 prot/met 人群；缺失人数写入 `c5.prs_cohort_flow.csv`。
- 存在 `C5_PRS_MANIFEST` 时，明确配置优先于自动字段。原 manifest 输入仍要求独立发现样本声明。
- 自动字段默认 `sample_overlap=unknown`，允许探索性分析并明确标注。外层交叉验证不能消除 PRS 发现 GWAS 与当前队列重叠造成的信息泄漏。

确认来源后可选填元数据；`none` 必须有真实依据，不应为了运行而填写：

```bash
export C5_DISEASE_PRS_SOURCE='实际 GWAS / PRS 来源及版本'
export C5_DISEASE_PRS_BUILD=b38
export C5_DISEASE_PRS_ANCESTRY='实际开发人群'
export C5_DISEASE_PRS_OVERLAP=unknown  # none / unknown / yes
```

`yes` 表示已知发现样本重叠：这类预先计算的 PRS 不进入验证比较，但不含 PRS 的模型继续运行。

Disease PRS 在输出中为 **G**。原来的 `prot.pgs.rds`、`met.pgs.rds` 是 **biomarker PGS**，其测量重建、独立筛选、PGS 捕获部分/剩余部分和预测分析全部保留。

## 3. 两类互补的比较

### A. 原有等总检测数量比较，继续保留

`Clinical_Protein_NS`、`Clinical_Metabolite_NS`、`Clinical_ProtMet_NS` 在同一个总 assay budget k 下比较。例如 k=10 时为 10 proteins、10 metabolites、5 proteins+5 metabolites。原 LE8 YS/YSplus、Yin/YinYang 比较保持这一规则。

### B. 新增固定每层 panel 的完整组合

C=clinical，P=protein，M=metabolite，G=disease PRS。PRS 可用时拟合以下 15 个非空组合：

| 范围 | 模型 |
|---|---|
| 单层 | C、P、M、G |
| 两层 | C+P、C+M、C+G、P+M、P+G、M+G |
| 三层 | C+P+M、C+P+G、C+M+G、P+M+G |
| 四层 | C+P+M+G |

输出模型名为 `F_C`、`F_P`、…、`F_CPMG`。PRS 缺失时保留 7 个不含 G 的模型，其他模型记为 unavailable。

**这一组 k 指每个组学层的 assay 数。** k=10 时，P 为 10 个蛋白，M 为 10 个代谢物，P+M 为 20 个 assay；G 始终是一个已计算的 disease PRS。真实数量随表输出。

每个外层训练折内分别选定 P 和 M panel，随后含该层的所有组合均使用同一 panel。加入 M 时保留 P 中全部蛋白，加入 G 时不改变 P/M panel。所有组合使用同一批人、相同外层分折、同一临床变量及时间窗。

输出全部 28 个“增加一层”的可配对比较。重点包括：

1. C+P+M+G vs C+P+M：已有实测双组学后，PRS 是否仍有增益？
2. C+P+M+G vs C+P+G：已有 proteome 和 PRS 后，metabolome 是否仍有增益？
3. C+P+M+G vs C+M+G：已有 metabolome 和 PRS 后，proteome 是否仍有增益？

固定 panel 不表示从全队列选定一次：**每个外层训练折重新筛选，验证折结局不参与筛选或拟合**。同一折内的重复模型复用拟合对象，减少重复计算。主临床模型为原定义的临床变量模型，不冒称为 SCORE2、PREVENT 或 FRS。

## 4. 评价指标与图

- 同一预测时间窗的 IPCW AUC、Brier、校准截距/斜率、校准图与决策曲线。
- 限制在预测时间窗内的 Harrell C；跨折汇总只比较各折内部的可比较个体对，避免不同模型 LP 尺度影响比较。
- 配对 ΔC、ΔAUC、ΔBrier、IDI、连续 NRI（同时报告 event/nonevent 部分）。
- 指定风险阈值后才报告分类 NRI，例如 `C5_RISK_CUTS=0.05,0.10`。这些阈值必须与具体 Y、预测年限和用途匹配；CAD 结局不能直接套用复合 CVD 指南阈值。
- 同一组 bootstrap 抽样同时计算各模型、所有增益和层贡献，相关参与者按组重抽样。区间条件于已训练模型，**不包含整套重新筛选/拟合的不确定性**。
- 对 C 基础上的 P/M/G 全部六种加入顺序计算平均增益分配（层级 Shapley decomposition），三层贡献之和等于 full model minus Clinical。它是预测增益的分解，不是因果贡献或遗传/环境效应比例。
- 蛋白、代谢物、PRS 单层评分的折内 Spearman 相关；训练集 75% 分位界定的三层高/低风险组合，报告 N、事件数、IPCW 风险和 KM 净风险区间。
- 保留 0/2/5 年 landmark、去 GDF15/利钠肽、去直接脂质重建指标的分析和所有之前的多面板图。

新增两张六 panel 图，各输出 PNG/PDF：

| 图 | Panel 内容 |
|---|---|
| `c5.Fig23.PRS_prot_met_prediction` | C-index、AUC、每层额外 ΔC、校准、决策曲线、平均层贡献 |
| `c5.Fig24.PRS_prot_met_complementarity` | ROC、层评分相关、三层风险分组、landmark、NRI/IDI、实际 assay 数–表现曲线 |

图默认显示预先设定的主预算和主 landmark；没有可估计结果时明确留空，不自动换成最漂亮的模型。点估计可为负增益，不预设联合模型一定胜出。

所有新汇总表均位于 `analysis/le8/[Y]/`：

| 文件 | 内容 |
|---|---|
| `c5.prs_provenance.csv` | Disease PRS 精确字段、来源、样本重叠及使用状态 |
| `c5.prs_cohort_flow.csv` | 共同人群、PRS 可用人数、最终比较人数 |
| `c5.factorial_design.csv` | 15 个计划模型和 PRS 缺失状态 |
| `c5.factorial_metrics.csv` | 每种组合的 AUC/C/Brier、区间和实际 assay 数 |
| `c5.factorial_contrasts.csv` | 所有一层增益及配对 ΔC/ΔAUC/ΔBrier/NRI/IDI |
| `c5.factorial_shapley.csv` | 相对于临床模型的平均层贡献 |
| `c5.factorial_status.csv` | 完整/缺失/失败预测，未运行不写成无增益 |
| `c5.factorial_strata.csv` | 三层评分高低组合及风险 |
| `c5.factorial_score_correlations.csv` | 单层评分在各验证折的相关 |
| `c5.factorial_ROC.csv` | IPCW ROC 绘图坐标 |
| `c5.metrics.csv`, `c5.calibration.csv`, `c5.decision_curves.csv` | 包括新增组合及原有模型的性能与绘图数据 |
| `c5.coefficients.csv`, `c5.fit_diagnostics.csv` | 显式系数、拟合状态及 measured / biomarker PGS / disease PRS 数量 |

默认 CI 计算覆盖 `all_assays` 下的主预算及各 landmark；其他预算/ablation 仍保留点估计，CI 为 NA 而不是伪造为零。`C5_FACTORIAL_BOOT_ALL=TRUE` 可拓展 CI 到全部设置，但计算时间明显增加。

## 5. 常用设置与保留的分析

```bash
export C5_ASSAY_BUDGETS=5,10,50
export C5_PRIMARY_BUDGET=10
export C5_LANDMARKS=0,2,5
export C5_PRIMARY_LANDMARK=5
export C5_END_YEARS=10
export C5_OUTER_FOLDS=5
export C5_BOOT=200
# 如仅关注常规基线起算的十年预测，可事先设置 C5_PRIMARY_LANDMARK=0。
# endpoint 固定在基线第10年：L=5 是存活/无事件者的第5至10年风险，不是另一个十年。
```

默认 ridge Cox 稳定处理共线性，同时输出显式系数和预处理参数；可用 `C5_RISK_SOLVER=cox`，病态矩阵会标记失败。已有 LightGBM、单组学 glmnet、轨迹、lead-time 等 reference 分析均保留。没有在这次改动中另行套用全队列残差后再交叉验证的流程。

LE8 支持 1–8 个有数据支持的领域，不强迫凑齐八个。`C4_FOCUS_COMPONENTS`、`C5_DOMAINS` 控制候选领域；`C5_PROXY_TARGETS=bmi,nonhdl` 保留 BMI/non-HDL 重建、加入临床模型和替换原指标的分开比较。没有验证 proxy 改善前，不称其为更好的临床指标。

有真实亲缘分组时使用 `C5_GROUP_COLUMN`；缺失时按 eid 分折并报告限制。已知记录覆盖期可提供 `LE8_ENDPOINT_MANIFEST`；模板中的空日期不可直接使用。AMR 仍须核对具体耐药结局、感染/药敏试验记录及记录覆盖，不能把药敏表成员直接当作耐药病例。

C1–C4、TwoSampleMR、MR-link-2、Dandelion、coloc/GPU-coloc、CIGMA、CellAge 注释和原图均保留。Dandelion 原生检验与项目聚合检验分开标注；CIGMA 仍要求真实 donor×cell-type 表达、SEM²、细胞比例和 kinship，不能用 PGS 或单位矩阵替代。CellAge 富集使用全 assay 基因背景，表达标签不证明血浆蛋白的细胞释放来源。

## 6. 文献依据与适用边界

1. Du et al. *Multi-omics integration predicts the incidence of 17 diseases in the UK Biobank*. Nature Communications (2026). [全文](https://www.nature.com/articles/s41467-026-73017-z)。已核对全文：采用共同人群，比较单/双组学，Cox martingale residuals + mixOmics `block.spls`，五次重复五折验证及 C-index。该研究发现双组学相对 proteome 的额外增益通常较小。这里借鉴其比较问题，使用现有训练折内筛选和生存预测框架；不是对 mixOmics 算法的逐项复现。
2. Ding et al. *Integrating clinical, proteomics, and polygenic scores to improve cardiovascular risk prediction: a prospective cohort study*. Journal of Advanced Research (2026). [论文](https://doi.org/10.1016/j.jare.2026.06.034)，[PubMed](https://pubmed.ncbi.nlm.nih.gov/42386074/)。已核对摘要及公开方法片段：LASSO Cox/稳定性筛选的蛋白评分，与临床评分和 disease PRS 比较，并评价 NRI、IDI、DCA；完整补充方案尚未取得。本模块借鉴增量评价设计，不能称为该文的复现。

这些新增分析回答“各数据层是否提供额外预测信息”，不据此认定蛋白是因果因素，也不把 GDF15/BNP 的预测能力本身判定为数据泄漏。Disease PRS 也不等于出生时已知的成人血浆蛋白水平。所有内部结果仍需要来源独立性核对及外部验证。

死亡目前作为删失，预测概率和 KM 曲线为净风险，不是考虑竞争死亡的累计发生风险。IPCW 使用边际独立删失假设。不同人群、结局和时间窗的 AUC/C-index 不能直接作论文间高低比较。

## 7. 代码检验

从 `scripts/` 目录运行（即本文件所在目录的上一层）：

```bash
python3 -m unittest discover -s le8/tests -p 'test_*.py' -v
Rscript le8/tests/test_prediction_audit.R
Rscript le8/tests/test_joint_prediction.R
Rscript le8/tests/test_final_joint.R
Rscript le8/tests/test_factorial_prediction.R
Rscript le8/tests/test_c5_attained_age_empty.R
Rscript le8/tests/test_c5_streaming_summary.R
Rscript le8/tests/test_c5_storage_cleanup.R
```

测试涵盖 PRS 精确字段/缺失/方差/来源优先级、模型组合、配对样本一致性、NRI/IDI、增益分解、预测相同则增益与区间为零，以及原先的基线时间和 proxy 验证。测试使用合成数据，不产生 UKB 结果。本次未在 UKB 个体数据上重跑，不能宣称 prediction 已改善。

本次实际通过：9 个 Python 测试、55 个 R 文件的语法解析，以及上述四组 R 合成数值测试（WebAssembly R 4.6.0，含 survival/glmnet）。这不等于在你的原生 R 环境和 UKB 数据上完成全流程验证。

新增两张六 panel 图也已使用明确标注的合成数据渲染为 PDF 并检查布局；这些测试图不包含在交付包中。真实 PNG/PDF 由你重跑 C5 后生成。

`data/cellage/` 和 `f/assets/` 是代码所需的外部参考注释、GO 名称及脑区轮廓，不是 AMR/CAD 分析结果。压缩包只保留当前 README；外部资源来源和许可收录在下方。


### f/assets/GO_SOURCE.md

Generic GO term-name dictionary, copied from the existing local ontology cache on 2026-09-12. This is the full term dictionary, not the article protein list or statistical results. Used for labels only; enrichment estimates are unchanged. Gene Ontology: https://geneontology.org/docs/download-ontology/



### f/assets/ggseg/SOURCE.md

Atlas polygons: ggseg v1.6.5, https://github.com/ggsegverse/ggseg/tree/v1.6.5/data . Desikan-Killiany and FreeSurfer aseg, MIT package licence. Region values are generated from the current LE8 analysis.



### f/assets/ggseg/LICENSE.md

# MIT License

Copyright (c) 2018 Athanasia Mowinckel & Didac Vidal Pineiro

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.



### data/cellage/LICENSE

MIT License

Copyright (c) 2026 Daisy Ding

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.


## 本机更新说明（2026-09-16）

- 本机 GWAS 路径保留为 `/mnt/e/gwas/{main,prot,met}`；`revision.conf.sh` 为 Python 辅助脚本指定已有的 `le8` 环境，补齐 `matplotlib` 依赖。
- 正常终端仅显示模块/组学层及主要步骤的 `[LE8] START / DONE`、失败和诊断信息；PGS 不再逐特征报告进度。详细输出仍保存到 `analysis/le8/logs/` 和疾病目录内的 `logs/`。`SCRIPT_VERBOSE=1` 可查看详细输出。
- `/mnt/d/scripts/0f/console_run.py` 的过滤调整仅适用于 LE8，其他脚本的显示策略保持原样。
- 分析目录清理保留可读取的 RDS、MR-link-2 标准化压缩输入和可续算的单任务结果。旧图表、日志、可重建的审计表及新版不再读取的疾病目录内重复 BED 文件删除；外部 GWAS、1KG 源数据和共享参考缓存不受影响。
- 重新运行 `./le8.sh --Y cvd_cad --biom prot,met` 会按现有缓存规则复用数值结果并生成输出；未设置 `--replace TRUE`。
