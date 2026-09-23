# Panome 5.1：学习个体匹配，再向真实参考人借用结局

新增独立入口 `panome_TF.sh`：使用已下载的 TabICLv2 数值表格 Transformer，对 prot/met → 疾病 Y 做监督微调并预测测试样本。模型、隔离策略、运行命令和限制见 [PANOME_TF.md](PANOME_TF.md)。输出位于单独的 `/mnt/d/analysis/panome_TF/`。

本版围绕你的 **row-based、individual-based、imputation-like** 思路重写。主模型输出的风险由真实 reference 的观察结局及其权重组成；Transformer 学习蛋白组合和个体之间的匹配。原始“挑选最合适的 100 人、找最近的人、COPY Y”仍作为独立模型完整保留。

本次是 v5 的 attention 深化版，内部版本 **5.1.0**。主模型包含两层含义：个人内部的分子模块 Q/K/V self-attention，以及个人对真实 reference 的多头 Q/K attention（Value 为真实 Y）。`ATTENTION_DESIGN.md` 说明哪些部分来自原 Transformer，哪些为了保留 COPY 解释而改变。方法细节见 `METHODS.md`。

主流程保留 top-fit、随机面板、距离匹配及独立训练的模型对照，用于检验 reference 选择和 attention 是否有效；运行时间本身不构成方法优劣的证据。

## 安装与正式运行

当前目录使用 v5.1 运行代码。`panome.sh` 优先使用 `PANOME_PYTHON`，否则复用 `$HOME/venvs/panome-v3/bin/python`，最后回退到 `python3`。本项目不创建 `.venv`，现有外部环境已满足依赖。若使用其他已配置环境，可先设置 `export PANOME_PYTHON=/path/to/python`。`panome.sh` 默认使用 CUDA、16 核、`quality-teacher all` 和运行名 `v5_attention_allteachers`，命令行显式参数可以覆盖这些默认值。不需要下载预训练模型。

```bash
cd /mnt/d/scripts/panome
bash panome.sh --check-device
./panome.sh --Y cvd_cad --biom prot,met --preflight

# 完整推荐分析：依次运行蛋白与代谢物，包括神经 OOF 榜样评分
./panome.sh --Y cvd_cad --biom prot,met

# 多个结局 × 多个分子层：依次运行四个独立分析
./panome.sh --Y cvd_cad,ra --biom prot,met
```

`--biom prot,met` 按顺序执行两个独立分析，各自使用对应的原始矩阵、映射、变换和批次列。蛋白默认残差化 `prot.plate`，代谢物默认残差化 `met.plate` 并使用原始代谢物映射及 `log1p`。结果分别写入 `/mnt/d/analysis/panome/cvd_cad/prot/v5_attention_allteachers/` 和 `/mnt/d/analysis/panome/cvd_cad/met/v5_attention_allteachers/`，首个分析失败时停止。多个 biom 时用 `--analysis-root`、`--run-name` 修改输出位置；需要 `--run-dir`、`--omics-file` 或 `--output` 时分别运行单个 biom。

原始代谢物映射配合 `log1p` 时，将负丰度作为无效测量转为缺失，并在 `input/metabolite_mapping.csv` 的 `negative_to_missing` 列逐指标记录。该固定规则在样本/特征 QC 之前应用，不使用结局；缺失值随后由训练集拟合的中位数填补。原始代谢物的无 Y 投影使用相同规则，并输出 `*_metabolite_mapping.csv`。`--transform none` 保留有符号值；已变换的蛋白或代谢物矩阵应使用此选项。直接提供的命名矩阵仍保留 `log1p` 非负检查。

这里的 Transformer 是从分子数据随机初始化并训练的自定义网络，BERT 的文本词表和权重不适配它，无需 Hugging Face 模型或 token。遮蔽重建预训练阶段不计算分类 logloss，因此日志显示 `N/A (reconstruction only)`；联合训练阶段才显示实际 logloss。

`--Y cvd_cad,ra --biom prot,met` 按 `cvd_cad/prot`、`cvd_cad/met`、`ra/prot`、`ra/met` 顺序运行，分别读取 `fod_icd10_cvd_cad` 或 `fod_icd10_ra`，结果目录为 `<analysis-root>/<Y>/<biom>/<run-name>/`。可先加 `--dry-run` 核对配置。多个结局时不能指定共用的 `--run-dir`、`--output` 或 `--diagnosis-col`；单个 biom 的 `--omics-file` 可跨结局复用。列表不接受空项或重复项，任何一次分析失败即停止后续运行。直接调用 `f/panome.py` 时仍须指定单个结局和单个 biom。

`--device cuda` 检测不到 GPU 会报错；`--device auto` 自动选择。若 `--check-device` 显示 CPU 版 PyTorch，请按 [PyTorch 官方安装器](https://pytorch.org/get-started/locally/) 选择与你机器相符的 CUDA wheel。不要仅凭安装了 `transformers` 就判断 GPU 已可用。

入口默认启用 `--quality-teacher all`。显式设置 `--quality-teacher ensemble` 时，榜样评分用重复 OOF elastic net + 小树模型，主流程仍会训练 Transformer、MLP、独立的无 retrieval Transformer、无预训练 Transformer，固定均匀 self-attention 对照、距离匹配对照，并保留自监督 encoder。`all` 把严格 OOF 的神经网络预测也加入评分，默认额外训练 3 × 5 个网络，计算量明显增加；每个 OOF 网络有自己的早停样本，其外层 OOF 对象从未参与该网络的预处理、token 构建或训练。

默认路径和终点：

| 项目 | 默认值 |
|---|---|
| 表型 | `/mnt/d/data/ukb/phe/Rdata/all.rds` |
| 蛋白 | `/mnt/d/data/ukb/phe/rap/raw/prot.tab.gz` |
| 首诊 / 基线 | `fod_icd10_cvd_cad` / `date_attend` |
| 死亡 / 失访 | `date_death` / `date_lost` |
| 行政截止 | `2023-04-01`，须核对是否适合你的数据版本 |
| 风险时间窗 | 10 年；死亡按删失处理的 net risk |
| 临床基准 | `age,sex,tdi,PC1,PC2,center` |
| 蛋白残差化 | **仅 `prot.plate`**；本版保留蛋白中的年龄、性别相关生物变化 |
| 主参考库 | 100 个真实人；300/1000 人及全部已知结局 build 样本作对照 |

默认残差化不同于 v4 的 `age,sex,prot.plate`。因此不能把跨版本差异全归功于 Transformer。本版同一次运行的各模型使用相同预处理；重复运行器可同时检验两种残差化策略。蛋白风险也不能据此解释为独立于年龄、性别的因果效应。

如果已经准备了亲缘连通分量，增加 `--group-col family_component`；外层拆分、OOF、神经训练 donor fold 排除、查询排除和 bootstrap 都使用亲缘组。未提供时不能声称排除了亲缘泄漏。可使用 CSV、TSV、gzip、RDS 输入；Parquet 另需安装 `pyarrow`。RDS 保留最新版代码中二进制传递和极小数值分类编码的修复。

## 本机检查与修复

2026-09-22 已完成 CPU/CUDA 合成端到端检查（1,600 人、800 人测试集），包括六种训练方案、42 个预测输出、CUDA 混合精度、LightGBM 和神经 OOF。15 项数学/attention 测试、CPU 断点恢复、无 Y 投影、查询顺序/批次一致性及测试结局隔离均通过。GPU 训练的冻结预测在 CPU 投影的最大差约 1.1e-6；固定拆分改变 800 人的测试结局后，重新训练的 30 个输出、OOF 评分、参考名单与神经参数保持一致。

修复 CUDA autocast 下概率 BCE 报错：编码器继续使用混合精度，reference 匹配与概率损失使用 float32。保留 RDS 二进制中转修复，并将 R 读取桥接与 shell 入口纳入断点恢复的代码指纹。

真实输入抽样 1,600 人已完成读取；排除既往病例后 1,521 人，QC 后 1,273 人、2,915 assays，按新版默认仅残差化 `prot.plate`，处理后的分子和临床矩阵全部为有限值。这是抽样预处理核查，不是全量训练。

历史审计材料、发布包验证记录和演示动画已从运行目录移除。训练时生成的 calibration-audit、attention 检查与输入/输出审计属于模型流程，继续保留。完整真实队列训练、全量资源消耗及科学有效性仍需正式运行验证。

## 多次训练、残差化敏感性与真正的负对照

```bash
"${PANOME_PYTHON:-$HOME/venvs/panome-v3/bin/python}" run_experiments.py --output-root /mnt/d/analysis/panome_suite \
  --seeds 2026,2027,2028 --residualization-sensitivity --include-null -- \
  --Y cvd_cad --biom prot --device cuda --cores 16 --quality-teacher all
```

默认固定第一个运行的外层拆分，再改变训练 seed，避免把不同测试集的结果误当成模型稳定性。所有运行先冻结，然后统一评估测试集。`--vary-splits` 可进一步研究拆分敏感性。输出按残差化/负对照方案分别汇总 AUC、Brier、覆盖率及其跨 seed 波动，同时比较 reference Jaccard 和同一批人的风险相关性。

`--include-null` 额外运行 development 结局打乱：build、tune、calibration-fit、calibration-audit 内分别打乱成对的随访时间和事件，测试集结局保持原样。它是完整拟合的负对照，**不是**置换检验的 p 值。普通主流程中的 `permuted_values_same_panel` 只扰动已学模型的 donor 结局，回答另一个问题：预测是否实际依赖 reference 的 Y。二者不能混称。

## 无结局信息的新个体预测

```bash
bash panome.sh project \
  --run-dir /mnt/d/analysis/panome/cvd_cad/prot/v5_attention_allteachers \
  --phe-file /path/query_baseline.csv --omics-file /path/query_proteins.csv \
  --output /path/query_predictions.csv --device cuda
```

这里只需要 ID、训练所需的基线协变量/批次，以及蛋白数据，不需要 diagnosis、event、time。训练中使用的 assay 列可以部分缺失；缺失比例计入支持判断。每位 query 只匹配 build reference；不会把同一批的其他待预测者作为 donor。输出预测、是否有支持、拒绝理由，以及旁边的 `*_references.csv.gz`。没有支持时保留探索性结果，但 `released_net_risk` 为空。

## 运行控制与计算量

- 默认 48 个分子 token、宽度 64、4 个头、2 层，batch 128。先遮蔽重建预训练 20 epochs，再联合训练最多 60 epochs，早停 patience 10。每个 epoch 记录损失、内层验证值并保存 checkpoint。
- `--tokens` 调整模块 token 数，不丢弃保留的蛋白。token 内对各 assay 使用不同的可学习参数，再跨 token 做 self-attention；这是分层压缩，不是声称实现了 3,000 个蛋白的全配对 attention。
- `--quality-neural-epochs 20 --quality-pretrain-epochs 5` 控制每个 OOF 神经网络的预算。
- `--experiments transformer` 可先检查主网络；默认会运行六个训练方案：`transformer,mlp,direct,no_pretrain,uniform,metric`。`ssl` 对照自动来自主 Transformer 的预训练快照。
- 显存不足先降低 `--batch-size`，再考虑 `--width`、`--tokens`；主库检索分批计算，不构造 50,000 × 50,000 的矩阵。完整数据仍需保存若干蛋白矩阵和 reference bank，请按实际机器观察内存日志。
- `--train-only` 冻结模型但不计算测试集指标；随后 `bash panome.sh evaluate --run-dir ...`。
- 中断后使用完全相同参数增加 `--resume`。会核对代码、输入指纹、依赖和配置，重跑准备/经典模型阶段，并从神经 epoch checkpoint 继续。`--full-input-hash` 可对大输入做完整 SHA256；默认大输入记录大小和修改时间。若强制终止留下 `.lock`，先核实其中 PID 不再运行再移除。
- 可以直接重复原命令：已完成的分析显示 `SKIP completed` 并跳过，未完成或失败的分析显示 `RESTART incomplete_or_failed`，清空原运行输出后从头训练。完整完成要求 `DONE.json` 和模型、报告、指标等核心文件存在；仅有 checkpoint 或 `MODEL_FROZEN.json` 不算完成。`--train-only` 使用训练结束后写入的 `TRAIN_DONE.json`。此规则对单个和批量分析均适用。
- 默认按输出目录判断完成状态，不因代码或参数变化自动重训已完成结果；需要更新它们时使用新的 `--run-name` 或显式 `--replace`。`--resume` 仍按严格指纹检查续训。覆盖操作全程持有 `.lock`，拒绝覆盖正在运行的目录，也不清理没有有效 Panome manifest 的非空目录。

v5.0 的 checkpoint/bundle 与本次 Q/K/V 参数结构不兼容；请用新目录重新训练，不要将旧模型复制到本版 run 目录。既往结果保留用于审计。默认超参数只是预设起点，不因为训练更久就自动更好。若要增加容量，可在开发集计划中比较 `--tokens 48/96 --width 64/128 --layers 2/4`，固定最终 test，不依据 test 反复挑模型。

## 主要输出

| 文件 | 回答的问题 |
|---|---|
| `reference_candidates.csv` | 每个人的 OOF 预测、相对常数风险的拟合增益、稳定性；不是人的“好坏” |
| `panel_*.csv` / `panel_specs.json` | 实际选了哪些人，病例覆盖、邻居数、温度、收缩量是多少 |
| `neural/*/learning_curve.csv` | 是否真正训练、早停在哪、重建和预测损失如何变化 |
| `embedding_diagnostics.csv` | 表征是否趋于常数、有效维度是否塌缩 |
| `test_reference_matches.csv.gz` | 每个 target 的真实参考人、Y、权重、风险贡献、删除/翻转 Y 的影响 |
| `individual_risk_decomposition.csv` | reference 贡献 + 先验贡献能否重构原始风险和校准后的风险 |
| `individual_explanations.jsonl.gz` | 抽样个体的共同异常蛋白、主要不匹配蛋白和解释边界 |
| `masked_reconstruction.csv` / `masked_feature_metrics.csv` | 遮蔽位置的 RMSE、R² 是否优于均值、PCA和随机参考人 |
| `module_perturbations.csv.gz` | 移除某模块的信息后，风险与匹配人是否改变 |
| `test_mosaic_matches.csv.gz` | 不同分子模块是否匹配到不同参考人 |
| `attention/self_attention.npz` | 真实逐人 × 层 × head × query token × key token 的权重矩阵 |
| `attention/reference_head_contributions.csv.gz` | 每个 head 对实际 donor Y 的贡献，能否重构总风险 |
| `attention/attention_interventions.csv.gz` | 改变 attention / 移除 head / 遮蔽模块后，风险及参考人如何改变 |
| `attention/attention_vs_sensitivity.csv` | attention rollout 与实际模块遮蔽敏感性是否一致 |
| `same_risk_pairs.csv` | 相似预测风险、不同分子组合的个体对 |
| `reference_readiness.json` | 独立 audit 中是否优于常数风险；参考库是否完整 |
| `support_error_audit.csv` / `subgroup_metrics.csv` | 支持规则是否仅排除了高风险人或特定人群 |
| `approach_comparison.csv` / `paired_contrasts.csv` | 所有方法及相同个体上的成对比较 |
| `raw_test_metrics.csv` / `test_calibration.csv` | 区分原始 COPY/borrowing 和校准的影响 |
| `REPORT.md` / `Fig_*.png` | 本次运行的汇总和科学图 |

默认另对 128 人导出完整 attention 和机制干预，使用 `--attention-samples` 调整或设为 0 关闭。逐 donor 表及无 Y 投影输出都包含各 head 的权重与贡献列。默认对 1,000 人做重复遮蔽重建，对 256 人做全部模块移除实验；reference 贡献表覆盖所有测试者。可用 `--explanation-samples`、`--perturbation-samples` 修改。任何模块名 `data_token_*` 都只是数据驱动分组。提供 `--module-file`（TSV：`feature`,`module`）才使用外部命名；重叠通路按模块名字排序取第一个，未覆盖蛋白另分组，映射完整写出。

本版为固定时间窗疾病分析，不实现基因型 phasing、连续表型、个体因果效应、竞争风险 CIF 或治疗数字孪生。代码能够计算并检验个体证据，不能预先保证 Transformer 或 100 人方案会胜出。
