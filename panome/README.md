# panome

**person → molecular state → disease**。默认用 UKB-PPP 蛋白数据，预测新发 `cvd_cad`；同一套算法也接受代谢组、转录组等“人 × 数值特征”矩阵。每次运行分析一个 omics layer，不把独立运行误称为联合多组学学习。

入口为 [panome.sh](panome.sh)，研究设计见 [DESIGN.md](DESIGN.md)。不依赖 LE8 的全局 R 环境，不修改 LE8 及 UKB 输入文件。

## 环境与运行

本机已经配置项目 `.venv`，继承现有 `ai` 环境的 PyTorch，并在项目内安装其他依赖。跨机器推荐独立环境：

```bash
cd /mnt/d/scripts/panome
python3 -m venv .venv
.venv/bin/python -m pip install -r requirements.txt
# R 需要 survival 和 data.table；不在分析时自动安装包。
./panome.sh -Y cvd_cad --biom prot --preflight
./panome.sh -Y cvd_cad --biom prot --dry-run
./panome.sh -Y cvd_cad --biom prot
```

本机默认输入：

- 表型：`/mnt/d/data/ukb/phe/Rdata/all.rds`。
- PPP：`/mnt/d/data/ukb/phe/rap/raw/prot.tab.gz`。保留原始缺失值，在 split 内填补。
- 默认输出：`/mnt/d/analysis/panome/cvd_cad/prot/main/`。
- `PANOME_PYTHON=/path/to/python` 可覆盖 Python；`--r-bin` 可覆盖 Rscript。

`prot.rds` 经 `/mnt/d/scripts/ukb/f/phe.R` 调用 `clean_biom(..., imp=TRUE)`，已经全人群填补。因此它仅作为兼容输入，不能据此声称从原始数据开始无信息泄漏：

```bash
./panome.sh --omics-file /mnt/d/data/ukb/phe/Rdata/prot.rds --run-name legacy_rds
```

快速真实数据试运行与软件模拟验证：

```bash
./panome.sh --run-name pilot6000 --max-samples 6000 --epochs 30 --hidden 256 \
  --stability-repeats 3 --bootstrap 30 --attribution-samples 200
./panome.sh --demo --run-name synthetic --epochs 8 --hidden 64 --latent 8 \
  --graph-dims 5 --stability-repeats 2 --bootstrap 5 --attribution-samples 30
```

试运行结果不替代全量分析；`--max-samples` 是按固定随机种子、在疾病筛选前取子样本。

## 六个阶段

| 阶段 | 内容 | 核心输出 |
|---|---|---|
| s1_prepare | 对齐 person ID，读取基线表型/原始分子数据 | phenotype.csv、omics.f32、features.txt |
| s2_preprocess | 基线排除、随访、60/20/20 split、训练内 QC/填补/校正 | cohort_audit.json、preprocessor.joblib、x.npy |
| s3_representation | Denoising autoencoder + PCA 参考 | autoencoder.pt、AE/PCA 坐标、重构损失 |
| s4_graph | 训练 AE 空间 kNN 图、Leiden、diffusion、留出者投影 | sample_graph.npz、graph_node_ids.csv、状态权重、稳定性 |
| s5_predict | 临床 / PWAS / PCA / AE / state / diffusion / elastic-net Cox | 独立测试指标、个体预测、bootstrap、状态 HR/PH 检查 |
| s6_report | 分子状态解释、个体非线性贡献与报告 | REPORT.md、图、person_molecular_states.csv、attributions |

支持 `--from s4_graph`、`--to s3_representation`、`--steps s5_predict,s6_report`。恢复时必须保留相同的分析参数，例如 pilot 恢复时仍需传入相同的样本数/epoch 等。`--replace` 必须从 s1 开始，会清理该 run 的阶段目录。使用新 `--run-name` 保存另一套分析。

Manifest 记录输入路径/字节数/修改时间、源码 SHA256、参数及依赖版本。缓存 v2 以 Python AST 分析代码指纹判断兼容性，忽略注释、帮助文字和帮助排版；参数默认值、分析逻辑、输入或依赖变化仍会阻止旧缓存复用。原版本的已验证运行可自动升级缓存元数据，原 manifest 保存在 `provenance_history/`，原阶段结果和完成标记保持不变；升级严格限定于 `provenance/cache_v1_to_v2.json` 记录的已审查版本，不会接受任意旧代码。输入大文件用 size+mtime 指纹，不是完整内容哈希。阶段成功才写 `DONE.json`。并发写同一个 run 会被 `.lock` 阻止。若进程异常退出，先检查 lock 中 PID 是否仍运行，再清理遗留锁。

## 人群、协变量和疾病

默认时间原点 `date_attend`：需要确认是所用 omics 的**基线采样对应日期**。PPP 原始表按上游代码作为 baseline proteomics 使用；若输入重复测量，必须先选定 visit 并提供对应采样日期。

- 结局日期：`fod_icd10_cvd_cad`，也可 `--diagnosis-col` 指定。
- 诊断 ≤ 基线：排除（同日也不是新发）。
- 截尾：`min(date_death, date_lost, 2023-04-01)`，默认截止沿用 `0phe.f.R`；可通过 `--end-date` 调整。
- 截尾当天诊断计为事件；截尾后诊断不是事件。
- 只有基线日期有效且随访 > 0 才纳入。
- 如果存在“诊断存在但首次日期缺失”的计数标记，可通过 `--disease-evidence-col` 提供，代码将排除这些人。
- 默认 **CAD-free**，不是论文“所有预定义疾病均未发生”。通过 `--healthy-date-cols` 提供额外首次诊断日期列，按基线排除；不会自动把所有 ICD 字段都解释为疾病。

```bash
# 下面的额外列需确实存在；代码对缺列报错。
./panome.sh --run-name broader_free \
  --healthy-date-cols fod_icd10_cvd_stroke_i,fod_icd10_cvd_htn
```

临床协变量默认 `age,sex,tdi,PC1,PC2,center`，对应 LE8 的 basic + center。这不是充分调整的临床 CAD 风险模型；正式论文需加入可靠的基线吸烟、血压、糖脂代谢及药物等，再以 `--covariates` 指定。

分子矩阵默认按 `age,sex,prot.plate` 做训练内 ridge 残差化（alpha=1），其中 plate 明确按分类处理。残差化使主分析回答“年龄/性别之外的分子差异”。推荐技术校正-only 敏感性分析：

```bash
./panome.sh --run-name technical_only --residualize prot.plate
```

`--categorical` 是分类列名列表，默认 `sex,center,prot.plate`，**不能把 plate ID 当连续量**。新 batch 不会凭空得到可靠校正：未知类别按全零 dummy 投影，并产生 sklearn 警告；外部新平台需重新验证。

可用 `--group-col family_id` 对亲缘连通分量分组划分。需事先把该列加入表型输入；每个亲缘组的成员不能跨 split。默认随机 person split 没有自动处理亲缘关系，也不自动限制祖源。原始 PPP 选择性入组会限制外推。

## 其他 omics

```bash
./panome.sh --biom met --run-name met_cad
# generic: phenotype 和 omics 可均为 CSV，也可 phenotype RDS + omics RDS/分隔文本。
./panome.sh --biom transcriptome --phe-file /path/phenotype.csv \
  --omics-file /path/normalized_expression.csv --run-name rna_cad \
  --residualize age,sex,rna_batch --categorical sex,center,rna_batch
```

Omics 必须有唯一、非空 `eid`，其余列都是数值特征；metadata 放在表型表。CSV phenotype 的日期使用 ISO `YYYY-MM-DD`，允许空值。代码拒绝重复 ID 和缺少所需列。

蛋白 NPX 不再次 log；其他组学必须先使用适合其测量尺度的处理，例如 RNA 的 library-size normalization/log-expression、甲基化的合适变换。不能把原始计数直接当作跨组学可比较的 NPX。跨层共同表征、MOFA、multi-view AE 尚未实现。

## 结果阅读

1. 先看 `s2_preprocess/cohort_audit.json` 与 split 的样本量/事件数。
2. 检查 `s3_representation/ae_history.csv`，AE 不保证优于 PCA。
3. 检查 `s4_graph/resolution_search.csv`、`graph_perturbation_stability.csv` 和连通分量；状态数不固定为 5。
4. `s6_report/person_molecular_states.csv` 中同时保留连续 AE/diffusion 坐标、离散状态和邻域权重。图节点顺序由 `s4_graph/graph_node_ids.csv` 定义。
5. `s5_predict/test_metrics.csv` 比较相同测试集。`metric_limitations.json` 明确记录随访支持不足的指标。
6. `s5_predict/test_state_cox.csv` / `test_state_ph_diagnostics.csv` 为探索性独立测试状态关联；检查少事件和 PH 假设。
7. `s6_report/person_ae_cox_attributions.csv` 给测试者的 AE-Cox 个体 integrated gradients；`attribution_completeness.csv` 给近似误差。参考点是经过预处理后的 0 向量，临床协变量固定。它描述模型 log-hazard 差异，不是每个蛋白的因果效应。

所有风险均为 death-censored Cox 的 **net risk**；不是死亡竞争风险下的实际累计发生率。所有 `.joblib` / `.pt` 只应加载可信本地产物。人级输出仍应按 UKB 访问要求保存，不上传到公共站点。

## 对新个体使用冻结模型

完成六阶段后，`f/project.py` 可读取不含疾病结局的外部 CSV，以已有 QC/校正/AE/graph/Cox 模型投影，完全不重新训练：

```bash
.venv/bin/python f/project.py \
  --run-dir /mnt/d/analysis/panome/cvd_cad/prot/main \
  --phenotype /path/new_baseline_metadata.csv \
  --omics /path/new_omics.csv --output /path/new_person_predictions.csv
```

外部蛋白列名需与训练时大写名称一致；必须覆盖训练保留的特征。需要临床/技术 metadata，不需要 Y/随访列。工具拒绝覆盖既有输出。输出状态和各模型 log-hazard/net risk；新平台、新 batch 的可靠性仍需外部验证。

## 本次已完成的全量运行

真实队列、模型比较及需要保留的解释边界见 [VALIDATION.md](VALIDATION.md)。全量输入 53,013 人，最终分析 42,190 人 × 2,919 蛋白。AE 已实际运行，但本次预测没有优于 PCA/elastic-net，候选社区的扰动稳定性也偏低。

额外中心分层诊断可单独执行：

```bash
.venv/bin/python validation/audit_results.py /mnt/d/analysis/panome/cvd_cad/prot/main
```

## final：最终汇总、论文图与配套 Excel

正常运行 `./panome.sh -Y cvd_cad --biom prot` 完成六阶段后，会自动生成 **8 对 400 dpi PNG 与同名 XLSX**。默认 `main` 的图表直接位于 `/mnt/d/analysis/panome/cvd_cad/prot/`，同时保存在 `main/publication/`；其他 run-name 的图表仅保存在对应运行目录的 `publication/`，避免覆盖主分析。

与 `gu.sh final` 一样，`final` 是独立汇总模块；不带参数时默认 PPP / cvd_cad / main。已有结果仅执行汇总，无需重新拟合：

```bash
./panome.sh final -Y cvd_cad --biom prot
```

| 文件前缀 | 内容 |
| --- | --- |
| Fig1.cohort_representation | 队列流程、数据划分、AE 训练曲线 |
| Fig2.molecular_landscape | 训练/测试分子坐标、状态蛋白谱 |
| Fig3.state_survival | 测试集净风险曲线、风险人数、状态 Cox HR |
| Fig4.prediction_benchmarks | 模型 C 指数、相对临床基线的配对差异及区间 |
| Fig5.calibration | 5/10 年校准、Brier 分数 |
| Fig6.individual_molecular_profiles | 相似预测分数个体的状态权重和分子贡献 |
| FigS1.state_diagnostics | 社区扰动稳定性、resolution、novelty |
| FigS2.landmark_sensitivity | 固定预测模型的 landmark 敏感性 |

每份 Excel 含 `Readme`、实际绘图值及相关原始统计表。`Fig_index.csv` 是文件目录，`Fig_legends.md` 提供英文图注。图中状态是探索性候选社区；净风险为 1−KM（死亡删失），不是竞争风险累计发生率。个体例子按模型预测分数选择，不按观察结局挑选；配套图表不含原始参与者 ID。

导出器 `f/final.py` 只读取已保存结果，具有独立代码/输入/输出校验；调整排版不会触发模型重训。核验已生成文件：

```bash
.venv/bin/python validation/check_publication.py /mnt/d/analysis/panome/cvd_cad/prot/main
```
