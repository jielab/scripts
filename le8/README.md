# LE8：Twin confidence、遗传校正分子偏离与汇总证据浏览器

当前目录只保留一套实现，统一入口为 `le8.sh`。`c1_correlate` 运行 PWAS/MWAS、PGS、时间分析、富集及常规临床/biomarker预测基准；`c1_abm` 运行 twin confidence、molecular deviation 及与C1基准的比较；`final` 汇总当前分析；`shiny` 浏览汇总证据。代码位于 `f/`，配置位于 `config/`，页面位于 `shiny/`。不提供历史模型入口、迁移分支或旧结果对比。

常规预测现在随 `c1_correlate` 自动运行，输出到 `<analysis-root>/<Y>/<biom>/c1_correlate/prediction/`，包括 clinical、biomarker-only、clinical+biomarker 的性能、增量及固定面板评估。特征筛选仅使用训练集；C2–C4结果不进入这个基础基准。关联分析与预测分别提交，已有C1关联结果会复用，只补缺失的预测阶段。`final --fit-reference` 已移除；原有预测配置的 `FINAL_*` 环境变量名暂时保留，`final --fit-joint` / `--fit-genetic` 仍用于另外的跨组学联合分析。

`c1_abm` 自动读取C1冻结的测试集预测，在两者测试集交集内、同一人群子集上重新计算 IPCW AUC 与 Uno C，输出到 `c1_abm/baseline_comparison/C1_ABM_comparison.xlsx`。比较使用同一个训练期删失模型，报告实际人数及覆盖率；时间窗或结局定义不一致时输出不可比较状态。C1的风险排序分数不当作概率，因此不据此比较Brier/LogLoss。差值为探索性点估计，不表示统计显著优于基准。两次运行传入相同的 `--outer-roster` 可保留相同的外部测试名单；未传入时使用各自原名单的交集，不改变已有划分。已有ABM拟合会复用，重跑 `c1_abm` 可刷新比较而无需 `--replace TRUE`。缺少客观终点的deviation仍不作诊断性能比较。

每个阶段完成（或复用已完成结果）后都会输出累计性能表和图，目录为 `<analysis-root>/<Y>/<biom>/<模块>/performance/`；`final` 的累计报告位于 `<analysis-root>/final/<Y>/<biom>/performance/`。主要文件是 `performance_summary.xlsx`、`performance_comparison*.png` 和同一测试人群的 `C1_ABM_comparison_*.png`。重跑 `c1_correlate` 也会读取已完成的ABM并刷新这些图，无需等到final，也不会重拟合ABM；如果ABM尚未运行则标记未评估。C2、C3、C5当前没有独立预测模型，其阶段报告保留先前已完成模型的比较，并将本阶段标为未评估；C4已有的验证模型加入累计报告。各图保留各自的终点、时间窗、样本量和验证口径，不把Harrell C和Uno C混为一个指标。累计报告只读取已有模型，不根据阶段证据另建模型。

从空结果目录重新运行：

```bash
cd /mnt/d/scripts/le8
./le8.sh --Y cvd_cad --biom prot,met --preflight
./le8.sh --Y cvd_cad --biom prot,met
./le8.sh final --Y cvd_cad --biom prot,met
./le8.sh shiny --Y cvd_cad --biom prot,met --no-reindex
```

默认完整流程包含 C1–C5；需要分析 RA 时设 `--Y cvd_cad,ra`。当前完成结果仍支持断点续跑，`--replace TRUE` 明确重拟合。家系与客观冠脉终点缺失时明确标记。

ABM twin 的 elastic-net 遇到迭代上限时，会保留当前系数续算，默认每次预算为 2000、4000、8000 次（`max_iter=2000`、`max_iter_retries=2`）；最多累计 14000 次，不放宽原收敛容差，也不改变 solver、惩罚、样本权重或数据划分。日志显示重试和收敛状态，训练审计记录实际迭代数与尝试次数。仍未收敛则明确失败，不发布未收敛模型；可通过 `--abm-config` 调整预算，`max_iter_retries=0` 禁用续算。失败输入保留在结果目录的 `_abm_staging_*` 中；重跑会重新准备未完成模块，已完成模块继续跳过，无需 `--replace TRUE`。终端失败摘要会显示 ABM 原始异常，而不是仅显示外层命令失败。

ABM 与 C1 统一默认读取 `<ukb-phe>/Rdata/prot.rds` 和 `met.rds`，表型来自同目录的 `all.rds`。已处理矩阵沿用原分子名和数值尺度，不重复做原始字段映射或 `log1p`。默认协变量为 `age,sex,tdi,PC1,PC2`，分类变量为 `sex`；`plate`、`well`、`center`（含 `met.plate`、`prot.well` 等带前缀形式）从协变量、额外调整和分子输入中排除，共享配置带入时会记录排除项。显式原始代谢物输入需同时指定 `--omics-file` 和 `--met-input raw`，才执行原始字段映射与相应默认变换。

已有 ABM 若使用原始输入或旧调整变量，不能直接当作新输入的完成结果跳过。执行 `./le8.sh c1_abm --Y cvd_cad --biom prot,met --replace TRUE` 会按新配置拟合并在成功发布时保留旧结果快照；此命令只运行 ABM，不重做 C1 关联分析。日志会显示实际输入路径、变换和调整变量。

## 2026-10-14 最小质控修复

此次仅修复分子缺失质控、偏离效应回归的冗余协变量，以及展示/下载口径。默认10分子预算、twin 距离、profile_concordance 筛选、Full/Clean/Soft/Random 和三类效应分解保持原设置。

- 分子列由 `always_include_indices` 明确区分。每个 OOF teacher 在自己的训练数据上冻结可用分子列、样本QC、面板及预处理；最终模型和 donor bank 同时检查冻结的分子集合与实际面板。临床哑变量不进入分子缺失率的分母。技术排除和 profile 不一致排除分别记入审计。
- query 的 `molecular_missing_fraction` 指实际冻结面板；`frozen_molecular_missing_fraction` 指选面板前冻结的可用分子集合。任一超过预设阈值均拒绝。全分子缺失者 confidence=0、selected=false、policy risk=NA，保留原始行与人数分母；技术合格但无 twin 时仍可使用 full-model fallback。
- 解释性状态均值差回归保留截距、用满秩 nuisance basis 去掉冗余协变量，再检查状态项能否识别；保留 cluster-robust SE/FDR，记录保留/删除列、设计秩和自由度。参考校准及七个疾病预测模型不重编码。
- Shiny 的图、PNG 和“当前筛选 TSV”共用 DT 全部分页的过滤行；零行保持为空，切换表/分子后旧行号失效。实现依据 [DT Shiny §2.2 的 rows_all 语义](https://rstudio.github.io/DT/shiny.html)。已知CAD、冠脉钙化、客观冠脉病变按 task 分别标注；无 task 的旧分子表展示为“记录状态/结局状态”，不改写历史源表。

**旧 ABM 拟合不能原地升级。** 新投影遇到旧模型时会要求显式重拟合。独立 `abm fit` 应复用原 development/test 文件和保存的配置（尤其临床列索引），输出到新的 run 目录。宿主流程可执行 `./le8.sh c1_abm --Y cvd_cad --biom prot,met --replace TRUE`，同时沿用原运行参数及原 `--outer-roster`；成功替换前会将旧模块保存到 `_history`。该命令会运行 c1_abm 的 twin/deviation 两部分，不重跑 C1关联及 C2–C5。仅重做解释性回归或展示时，可用旧 deviation 模型与预测执行 `deviation evaluate` 到新目录，再更新 evidence/展示，无需重拟合七个疾病模型。

宿主显式配置 `task=known_cad_state` 或 `coronary_calcium` 时，分别存入 `c1_abm/<task>/molecular_deviation/`；默认客观病变仍使用原 `c1_abm/molecular_deviation/`。这些独立任务可在证据浏览器中查看，不混入默认客观病变汇总。

没有家系信息时，宿主流程及独立 `abm prepare` / `deviation prepare` 均按个人 ID 分组并标记 `individual_unverified`，不自动添加 `assume_unrelated=true`。输出明确记录 `family_status=unavailable_individual_ID_groups_only`、`independent_family_support_verified=false`；无客观检查仍可校准，诊断评价为 unavailable。原测试名单保持不变，已查看/未核实的测试结果继续标记 exploratory。

回归检查（在本目录运行）：

```bash
python -m pytest -q tests/test_qc_repairs.py tests/test_prepare_optional_inputs.py -p no:cacheprovider
Rscript --vanilla tests/test_shiny_filters.R
```

R 检查需要 Shiny、DT、DBI、RSQLite、ggplot2、plotly、testthat；使用 SQLite 合成表和 `shiny::testServer` 实际验证过滤、分页、下载、切换范围、零行与 task 绘图。合成验收覆盖20名全缺失 query、30名全缺失 donor、dummy数量不变性、同人单独/批量投影、家系排除、原 outer roster、共线状态不可估计等。此次修复前后合成完整样本的 ABM 预测/选择/confidence 一致，参考校准和七模型预测也一致；有缺失样本的模型会因排除不合格训练者、重新冻结预处理/半径等发生必要数值变化。这里不声称完成真实 UKB 重跑或整个仓库的历史测试。

同日使用 `le8_audit_evidence_20261014` 中未经改写的 `reproduce_codex.py` 再验收：20名全缺失 query 的QC通过数、支持数、入选数与最大confidence均为0，30名全缺失训练者无一保留；三个分子的完整one-hot状态效应与参考水平编码一致。另补齐两个独立prepare入口的可选家系处理。26项当前Python测试及R `testServer` 检查通过；原附件的22项旧测试仅提供日志，未冒充本轮已重跑。原证据文件未覆盖，当前源码哈希、复现前后数字和检查范围保存在 [核查记录](tests/audit_recheck_20261014.json)。

## 1. 本版落实的研究目标

### ABM 的主线

1. 在 build 数据内选出模型使用的分子面板（默认10个分子；明确列出的临床协变量不占这10个名额）。
2. 每个训练 donor 用重复、家系分隔 OOF 模型评价 **profile 与观察结局是否符合该模型**。同样是模型预测低风险的 profile，观察到 case 的 donor 会获得较低的模型一致性。
3. 主方案按模型一致性过滤 donor；保留全训练、软权重、病例×风险构成匹配的随机过滤对照。`class_risk_stratified` 是预声明敏感性方案，不替代主假设。
4. 在冻结面板的整个标准化特征向量中找 twin。系数只调整距离权重，不把整个 profile 压成一个总风险。固定预算面板中若有pilot系数为零或最终系数被压缩为零的项，会在training_audit记录数量；不把“名义10项”冒充10个非零作用。每一维保留非零距离权重。
5. query 结局未知时计算 0–100 `pre_prediction_confidence`。证据包含距离、独立家系支持度、邻居 donor 的 OOF 一致性，以及匹配半径。没有充分支持时得分为0且专门预测不适用；0不代表健康。
6. confidence 主方案不取决于 attention 是否胜过 baseline。预期误差、AUC排序误差、普通风险置信度、低风险选择、随机选择是补充或对照。
7. 分别检验子集选择收益、同一子集中的训练过滤收益、合并收益。总体和拒绝者仍保留在评价分母中。

`confidence = 100 × (proximity × independent_support × donor_concordance × profile_radius)^(1/4)`，且必须满足技术QC及最小 twin 支持。这是透明的参考支持分数，**不是“预测正确的百分比”**；独立 audit/test 的 `confidence_validation` 表检验其与实际误差的关系。

“一次测量的10个蛋白无明显异常”在代码中是相对于训练分布和模型的 profile，不冒充同一个人的纵向“没有变化”。低模型一致性不证明原始病例错误。默认仅删除5%；若删除后病例/对照不足，流程报错，不偷偷将病例加回。

### 分子遗传偏离的主线

按每个分子的自己的 PGS，拟合 `M_j ~ clinical/ancestry/technical covariates + G_j`，与 `M_j ~ covariates` 对比。参考资格默认“采血时无已知CAD”；不使用未来发病排除参考人群，也不称其为影像正常健康人。

以家系交叉拟合得到参考残差尺度，计算 `Z_j=(M_j−expected_j)/residual_SD_j`。训练疾病模型使用 OOF Z；独立 tune/calibrate/test 使用冻结参考模型。

主终点是无已知CAD史者中、与采血相近时间的客观冠脉病变。已知CAD状态与冠脉钙化是两个**独立任务**，不能替代主终点。没有客观终点时只完成参考校准，状态是 `not_validated_for_diagnosis`，不生成虚构ROC/AUC。

七个比较模型：Clinical、Clinical_G、Clinical_M、Clinical_M_G、Clinical_Z0（不含PGS的残差）、Clinical_Z、Clinical_Z_G。残差化不创造信息，最后一组检验与原始M+G的重参数化关系。弱遗传捕获保留标记，不能把Z解释成纯环境或疾病成分。

### 页面组织

`evidence` 只建证据索引，不重新拟合。可读取原有 `.xlsx` 的 `le8/manifest.json` 与精确内嵌导出，也可读普通Excel sheet及新增聚合表。不会要求先手工导出CSV。

`viewer` 按 gu 的模块化方式拆成 data/ui/plots/server。全局下拉：疾病、组学、模块、Table、protein/metabolite。Overview → molecule explorer → twin confidence → genetic deviation → audit。选择一张表、一个分子才动态出图；不循环生成3000张静态图片。表格分页/过滤，分子选项 server-side selectize，点击行可定位分子。小单元格抑制的分布信息保留为缺失，不改成0。

共享页面不导入个人ID、个体预测、donor、model或原始矩阵。每个testing sample的confidence保存在本地 `LOCAL_ONLY` 文件中；默认没有把个体信息暴露到公共网页的开关。

## 2. 依赖与入口

依赖由 `environment.yml` 统一管理：

```bash
./install.sh --env le8
# 或只补充现有 Python/R 环境所需依赖
./install.sh --abm
```

所有常规分析和独立模型工具均通过 `le8.sh` 调用，无需另外的安装器或启动脚本。下面的 `abm` / `deviation` 命令用于显式输入文件的独立研究；日常全队列运行使用上面的 `c1_abm` 调度。

## 3. ABM 使用

数据准备：`cohort.rds` / TSV / Parquet 必须含 `eid`、原宿主已推导好的随访time/event、基线CAD状态；家系编号可选。不要为使用本包重新定义原CAD结局。可通过 `--time-col`、`--event-col`、`--prevalent-col` 指向原列名。

```bash
cd /mnt/d/scripts/le8
./le8.sh abm prepare \
  --cohort-file /path/to/derived_cad_cohort.rds \
  --omics-file /path/to/prot.rds \
  --config config/abm.json \
  --time-col time --event-col event --prevalent-col known_cad_at_index \
  --clinical age,sex,PC1,PC2 \
  --output /mnt/d/analysis/le8/cvd_cad/prot/abm_twin_prepared
```

`/path/to/...` 是需要替换的实际数据文件，程序不假装能访问或自动推导尚未确认的临床字段。临床分类变量必须沿用宿主固定编码；稀有类别的编码不能用test结局决定。有家系资料时可加 `--group-col family_group`；没有时直接运行，按个人 ID 分组并明确亲缘未核实。`--assume-unrelated` 仅兼容显式用户声明，运行不需要此选项，声明也不等于家系支持已核实。已有角色名单应通过 `--role-col` 保留，或在宿主入口继续使用原 `--outer-roster`。

```bash
P=/mnt/d/analysis/le8/cvd_cad/prot/abm_twin_prepared
O=/mnt/d/analysis/le8/cvd_cad/prot/abm_twin
./le8.sh abm fit --development "$P/_LOCAL_ONLY/development.npz" \
  --config "$P/config.json" --output "$O/fit"
./le8.sh abm project --model "$O/fit/model_LOCAL_ONLY.joblib" \
  --features "$P/_LOCAL_ONLY/test_X.npz" --output "$O/test_predictions_LOCAL_ONLY.joblib"
./le8.sh abm evaluate --model "$O/fit/model_LOCAL_ONLY.joblib" \
  --predictions "$O/test_predictions_LOCAL_ONLY.joblib" \
  --truth "$P/_LOCAL_ONLY/test_truth.npz" --output "$O/test_report"
```

投影输出旁边有 `test_predictions_LOCAL_ONLY.confidence.tsv.gz`，含每人的confidence、支持分量、是否入选和理由。完整概率和诊断数组在该本地joblib中。加载joblib只限可信本地模型文件。

默认10年净风险目标；死亡按宿主给出的删失处理，非竞争风险累积发生率。输出分别保留 IPCW AUC、Uno C、Brier、LogLoss、人数覆盖率、事件覆盖率。不能把Uno C改称AUC，也不把文献0.75当作本实验固定对照。

`model_concordance` 的筛选与普通数据QC不是同一件事。原始标签残差筛选容易削减病例，务必看 `training_audit.json` 和随机对照。默认只删5%是运行起点，不是已证明最优阈值；所有阈值在看新test前确定。

## 4. PGS读取与精度

输入识别与当前 `gwas/format.sh pgs` 一致：

- 合并表：`/mnt/f/gwas/met/pgs/met.pgs.gz`，`eid + <assay>.SCORE_SUM`；prot同理。
- 单分子表：`#IID + <assay>.ALLELE_CT + <assay>.SCORE_SUM`。
- 只按eid对齐，不能依赖原行顺序。
- 不把 `.SCORE_AVG` 自动替换成 `.SCORE_SUM`。
- 不模糊匹配不同assay，不把百分比映射成绝对量，不合并 NTPROBNP/NPPB。
- 无工具变异的header-only输出与全常数分数标记不可用，不补0。

当前上游合并脚本使用 `%.3f` 并将很小的分数写为0；新模块报告 `rounded_3dp_fraction`，但它只是精度风险提示，不证明每个分数都已严重损坏。可以不重跑GWAS，直接从保留的逐分子结果重合并：

```bash
./le8.sh pgs-merge \
  --manifest /mnt/f/gwas/met/pgs/met.pgs.files.tsv \
  --output /mnt/f/gwas/met/pgs/met.pgs.fullprecision.gz
```

输出create-only、17位有效数字；无法恢复逐分子原文件自己已丢失的精度。大量蛋白合并需要相应内存。PGS输入读取分块，仅保留研究所需eid。

PGS总分和trait级ALLELE_CT不能验证每个人的SNP组合、等位基因方向、祖源距离或每个SNP缺失。ABM可评价传入的特征profile；原始基因型适用性需要另有variant级数据，当前会明确记录未评估。

## 5. 遗传校正分子偏离使用

配置默认met；prot使用 `deviation_prot.json`。`gwas_overlap` 必须按真实来源记录，交叉拟合参考模型不消除原始GWAS样本重叠。

```bash
./le8.sh deviation prepare \
  --cohort-file /path/to/cad_cohort.rds \
  --omics-file /path/to/met.rds \
  --pgs-root /mnt/f/gwas/met/pgs/met.pgs.fullprecision.gz \
  --clinical age,sex,PC1,PC2 \
  --config config/deviation_met.json \
  --output /mnt/d/analysis/le8/cvd_cad/met/deviation_prepared
```

对名称不一致的指标传 `--pgs-mapping`，格式见示例。PGS目录有多个候选时必须给确切文件；绝不猜用哪个。

仅做客观病变诊断评价时，所需终点字段见 `config/objective_endpoint_schema.tsv`；这些字段不是参考校准与偏离计算的前置条件，家系列同样可选。`objective_cad` 是研究者事先定义并独立判读的0/1病变终点，程序不擅自把ICD未来诊断当影像病变。默认采血与冠脉检查相差不超过90天；日期不同需确认间隔内没有相关疾病事件或新治疗。90天是预声明默认，需按研究方案确定，不是通用诊断标准。

```bash
P=/mnt/d/analysis/le8/cvd_cad/met/deviation_prepared
O=/mnt/d/analysis/le8/cvd_cad/met/molecular_deviation
./le8.sh deviation fit --data "$P/_LOCAL_ONLY/development.npz" \
  --config "$P/config.json" --output "$O/fit"
./le8.sh deviation project --model "$O/fit/_LOCAL_ONLY/model.joblib" \
  --data "$P/_LOCAL_ONLY/test_X.npz" --output "$O/prediction_LOCAL_ONLY.joblib"
./le8.sh deviation evaluate --model "$O/fit/_LOCAL_ONLY/model.joblib" \
  --predictions "$O/prediction_LOCAL_ONLY.joblib" --truth "$P/_LOCAL_ONLY/test_truth.tsv.gz" \
  --anchors L_VLDL_TG.pct,Total_TG,ApoB --output "$O/results"
```

主要输出：遗传校准OOF R²，以及有可用终点时的同人群七模型比较、四组配对增量对比、ROC、概率校准、固定特异度阈值下表现、探索性decision curve、所有分子的状态分布及校正均值差。`feature_state_association` 的β是调整后的Z均值差，不是疾病OR；标准误按提供的家系或未核实亲缘的个人 ID 分组计算，FDR按全指标族。

共享图不使用“最显著feature”自动选主角。默认只画预声明锚点；全指标在Shiny下拉查询。检测样本是被选择去做冠脉检查者时，PPV和net benefit只对应该样本，不自动代表一般人群筛查。

没有影像终点仍可完成fit_reference并输出校准，诊断模型为空、状态明确；不能把缺失影像者当阴性。若研究者已有已知CAD状态分析，可将task设为 `known_cad_state`，但必须是另一结果目录，不能叫隐匿冠脉病变诊断。

## 6. Final / Shiny

```bash
./le8.sh evidence --analysis-root /mnt/d/analysis/le8 --Y cvd_cad --biom prot,met
./le8.sh viewer --db /mnt/d/analysis/le8/evidence/evidence.sqlite --port 3839
```

浏览本机 `http://127.0.0.1:3839`。默认只绑定localhost。索引完成后重启或刷新浏览器新会话以读取新数据库；运行中的会话保持自己的只读快照。

`evidence/evidence.sqlite`、`catalogue.xlsx`、`index.html` 分别存放浏览器索引；`final/` 保存当前汇总报告。SQLite每张源表保持独立来源、run、单位和hash，不自动跨运行合并、不选择最小P/最大H4、不把一条placeholder算一项发现。超过展开大小限制的表明确写入audit，界面不能声称索引完整。默认128MB/源表，可配置 `--max-table-mb`。

动态图包括单分子效应、遗传—实测关联、OOF PGS捕获、偏离分布/均值、ROC、calibration、decision curve、twin confidence误差曲线、coverage与配对对比。其他表提供自选X/Y轴。相关时间窗/HR/OR先筛选范围再比较，页面不会自动把它们合成一个因果效应。

## 7. 合成演示与结果分享

```bash
./le8.sh abm demo --output /tmp/le8_abm_demo_unique
./le8.sh deviation demo --output /tmp/le8_deviation_demo_unique
```

合成演示用于检查程序运行，不代表真实队列效果。

`./le8.sh share --share-out /path/le8-share.zip` 导出汇总 SQLite、索引工作簿和完整 Shiny 页面。个体预测、矩阵和模型保存在本地 `LOCAL_ONLY` 路径，不进入浏览器或分享包。
