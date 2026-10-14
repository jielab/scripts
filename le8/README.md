# LE8：Twin confidence、遗传校正分子偏离与汇总证据浏览器

当前目录只保留一套实现，统一入口为 `le8.sh`。`c1_correlate` 运行 PWAS/MWAS、PGS、时间分析及富集；`c1_abm` 运行 twin confidence 与 molecular deviation；`final` 汇总当前分析；`shiny` 浏览汇总证据。代码位于 `f/`，配置位于 `config/`，页面位于 `shiny/`。不提供历史模型入口、迁移分支或旧结果对比。

从空结果目录重新运行：

```bash
cd /mnt/d/scripts/le8
./le8.sh --Y cvd_cad --biom prot,met --preflight
./le8.sh --Y cvd_cad --biom prot,met
./le8.sh final --Y cvd_cad --biom prot,met
./le8.sh shiny --Y cvd_cad --biom prot,met --no-reindex
```

默认完整流程包含 C1–C5；需要分析 RA 时设 `--Y cvd_cad,ra`。当前完成结果仍支持断点续跑，`--replace TRUE` 明确重拟合。家系与客观冠脉终点缺失时明确标记。

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

数据准备：`cohort.rds` / TSV / Parquet 必须含 `eid`、原宿主已推导好的随访time/event、基线CAD状态、家系编号。不要为使用本包重新定义原CAD结局。可通过 `--time-col`、`--event-col`、`--prevalent-col` 指向原列名。

```bash
cd /mnt/d/scripts/le8
./le8.sh abm prepare \
  --cohort-file /path/to/derived_cad_cohort.rds \
  --omics-file /path/to/prot.rds \
  --config config/abm.json \
  --time-col time --event-col event --prevalent-col known_cad_at_index \
  --group-col family_group --clinical age,sex,PC1,PC2 \
  --output /mnt/d/analysis/le8/cvd_cad/prot/abm_twin_prepared
```

`/path/to/...` 是需要替换的实际数据文件，程序不假装能访问或自动推导尚未确认的临床字段。临床分类变量必须沿用宿主固定编码；稀有类别的编码不能用test结局决定。无家系字段时，仅在已经排除亲属后显式 `--assume-unrelated`。

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
  --cohort-file /path/to/cad_objective_cohort.rds \
  --omics-file /path/to/met.rds \
  --pgs-root /mnt/f/gwas/met/pgs/met.pgs.fullprecision.gz \
  --clinical age,sex,PC1,PC2 \
  --group-col family_group \
  --config config/deviation_met.json \
  --output /mnt/d/analysis/le8/cvd_cad/met/deviation_prepared
```

对名称不一致的指标传 `--pgs-mapping`，格式见示例。PGS目录有多个候选时必须给确切文件；绝不猜用哪个。

必要主终点字段见 `config/objective_endpoint_schema.tsv`。`objective_cad` 是研究者事先定义并独立判读的0/1病变终点，程序不擅自把ICD未来诊断当影像病变。默认采血与冠脉检查相差不超过90天；日期不同需确认间隔内没有相关疾病事件或新治疗。90天是预声明默认，需按研究方案确定，不是通用诊断标准。

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

主要输出：遗传校准OOF R²、同人群七模型比较、四组配对增量对比、ROC、概率校准、固定特异度阈值下表现、探索性decision curve、所有分子的偏离分布及校正均值差。`feature_state_association` 的β是调整后的Z均值差，不是疾病OR；标准误按家系聚类，FDR按全指标族。

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
