# Panome TF：用预训练数值 Transformer 学习 prot/met 与疾病 Y

`panome_TF.sh` 是独立的监督预测入口。原 `panome.sh` 已有自行训练的分子 Transformer，主模型的风险来自真实参考人的结局加权；本方案使用下载的 TabICLv2 权重、实际梯度微调和学习到的分类头，输出测试样本的疾病风险与二元预测。原程序和正在运行的任务继续保留，新结果写入 `/mnt/d/analysis/panome_TF/<Y>/<biom>/tabiclv2_finetuned/`。

## 从 EMS120 借鉴的部分

- `ml_phone`：从特征和人工标签训练模型、比较候选模型、保存后批量预测。其当前实现为 elastic net／随机森林／boosting，并不是 Transformer。
- `ml_dx`：加载预训练 MacBERT，使用训练标签做梯度微调，独立验证，保存权重后批量推断。
- Panome 的输入是连续分子丰度，不是中文诊断文本。因此借鉴训练与验证流程，使用数值表格模型，不把数千项测量拼接成可能截断的文本。

## 模型和下载

采用 [TabICLv2 官方实现](https://github.com/soda-inria/tabicl) 与 [官方分类权重](https://huggingface.co/jingang/TabICL)。模型先在合成表格任务上预训练，再适配本地研究数据；不是预先训练好的 UKB 疾病诊断模型。[论文](https://arxiv.org/abs/2602.11139)

- 本地文件：`/mnt/e/AI_model/TabICL/tabicl-classifier-v2-20260212.ckpt`，即 `E:\AI_model\TabICL\...`。
- 大小 110,368,038 bytes，27,552,258 个参数。
- 仓库 revision：`4dcd344ece2c00be9e831fdd35bed57b5ad83e19`。
- SHA256：`bdc7dbd5e4ff21f8f0456fcf90c6b7cdf72dbea960f2d05b19bec19f9b3d4ed0`。
- 本次已通过 Git + Git LFS 下载真实权重并验证 SHA256；目录内其余 checkpoint 可能仍是 LFS 指针，本程序只使用上述分类文件。
- `./panome_TF.sh --download-model` 可下载或核验后退出。先尝试 Git LFS，失败后使用与 `0hf_download.py` 相同的 `snapshot_download` 方法。公开权重无需 token；如需要认证，从 `HF_TOKEN` 环境变量读取，不写进代码、命令参数或输出文件。

继续使用 `PANOME_PYTHON` 或现有外部 Python 环境，不创建项目 `.venv`。本次在现有解释器安装了 `tabicl==2.2.0`，没有升级原分析依赖。其他机器可在既有环境安装 `requirements-tf.txt`。微调适配层调用与官方 fine-tuning 相同的底层 forward，并固定包版本；升级 TabICL 前应重新验证。

## 训练设计

1. **结局**：沿用已有的首诊／基线／死亡／失访日期，默认预测基线之后 10 年内首次发病的 net risk。排除基线已患病者。未达到时间窗且被删失的人不当作阴性；其监督权重为 0。删失分布仅由 build 估计。需要普通横断面 0/1 标签时，应另行定义任务，而不是混用这套随访标签。
2. **外层隔离**：默认 build 30%、tune 10%、calibration 10%、test 50%。可传入已有 `--split-file`，文件包含 `eid,split`。提供 `--group-col` 时按亲缘连通分量分组。初始分层可用结局，但 test 的 Y 不参与任何拟合、早停、阈值或特征选择。
3. **预处理**：特征 QC、截尾、中位数填补、批次残差化和标准化均只拟合 build。原始代谢物负丰度在 `log1p` 前记为缺失并记录。默认蛋白保留有符号值。
4. **输入容量**：默认按 build 内的 IPCW 加权单变量关联选 256 个分子指标，保持原特征顺序；保存全部筛选分数和入选名单。这样控制蛋白的 attention 计算量。单变量筛选可能漏掉纯交互信号，`--max-features 0` 可保留全部 QC 合格指标，但内存和运行时间显著增加。容量比较只能在开发集进行。主网络默认只输入分子数据；`--with-clinical` 可加入临床协变量。
5. **实际梯度微调**：从官方权重初始化，AdamW、学习率 `1e-5`、最多 10 epochs、patience 3，CUDA 混合精度和梯度裁剪。默认冻结列／行编码器，更新 ICL Transformer 与分类头；仍约有 2,628 万个可训练参数。`--unfreeze-encoder` 可微调全网络。每个已知结局 build 样本在每个 epoch 作为 query 参与一次加权交叉熵。
6. **上下文隔离**：TabICL 仍是 in-context Transformer：输入一组带 Y 的 build context 和待预测样本，经过学习的分类头输出概率。训练 query 的 context 必须来自其他内层 fold，亲缘组也不能跨 fold。默认最多 256 个不重复 context，按 IPCW 和类别加权比例抽取，而非人为配成 50:50。query 损失显式使用 IPCW；context 抽样只是近似加权，不能称为精确的 IPCW attention。最终预测使用固定的 build context；test Y 从不进入它。
7. **选择与校准**：默认在最多 2,048 个已知结局 tune 样本早停，保存这些 ID；`--validation-samples 0` 使用所有已知 tune 样本。未微调的 epoch 0 同样参与选择，若微调无改善就保留基础模型，不强行宣称微调有效。独立 calibration 集做单调 logistic 风险校准；二元阈值用 tune 集的加权 Youden J 确定并映射至校准后的概率。
8. **对照和测试**：临床 logistic、同一入选分子特征的 elastic net、常数风险对照。各自的超参数只在 tune 选择，校准只用 calibration。冻结模型、context、阈值和测试预测以后，才计算 test 的 IPCW AUC、AUPRC、Brier、logloss、灵敏度、特异度和分箱校准。预测不是临床诊断；本版没有输出置信区间，也不声称排除了未提供的亲缘关系。

TabICL 的普通 `fit()` 本身不做梯度训练。本代码在 `tf_model.py` 显式执行反向传播，`training_summary.json` 记录梯度步数、选中 epoch、与预训练权重相比改变的 tensor 数。`selected_weights.ckpt` 是微调后可重载的权重，`model_bundle.joblib` 同时封装预处理、固定 context、分类权重、校准和阈值，可离线预测。

## 运行

```bash
cd /mnt/d/scripts/panome
./panome_TF.sh --Y cvd_cad,ra --biom prot,met --preflight
./panome_TF.sh --Y cvd_cad,ra --biom prot,met

# 一个较短的合成端到端检查，不代表真实数据效果
./panome_TF.sh --demo --max-samples 1200 --demo-features 12 \
  --epochs 1 --max-features 12 --context-size 32 --batch-size 32 \
  --predict-batch-size 64 --validation-samples 80 --folds 2 --min-events 5 \
  --run-name synthetic_check

# 调整计算量；这属于开发阶段的预设比较
./panome_TF.sh --Y cvd_cad --biom prot --max-features 128 \
  --context-size 128 --batch-size 32 --run-name smaller

# 使用既定拆分便于与原模型公平比较，单次指定一个 outcome/layer
./panome_TF.sh --Y cvd_cad --biom prot \
  --split-file /mnt/d/analysis/panome/cvd_cad/prot/v5_attention_allteachers/split_before_qc.csv
```

完成的分析自动跳过；未完成／失败的分析清理后重新训练。`--replace` 强制重训已完成结果。运行锁防止并发覆盖，不清理无有效 manifest 的非空目录。该入口不实现 epoch checkpoint 续训，失败重跑从下载的基础权重开始。`--train-only` 冻结训练结果后停止，可随后执行 `evaluate --run-dir ...`。

四个组合依次运行，首次失败时停止，修复后重复原命令即可。原 `panome.sh` 正在占用 GPU 时，建议待其结束再启动全量 TF 批次；本次只进行了短小验证，没有启动第二个全量批次。

## 新样本预测

```bash
./panome_TF.sh project --biom met --met-input raw \
  --run-dir /mnt/d/analysis/panome_TF/cvd_cad/met/tabiclv2_finetuned \
  --phe-file /path/query_baseline.csv --omics-file /path/query_raw_met.csv \
  --output /path/query_predictions.csv
```

无需 Y／诊断／随访列，仅需训练所需的 ID、协变量／批次和分子测量。命名代谢物矩阵使用 `--met-input named`。context 中已出现的 ID／亲缘组不能冒充独立新样本。缺失指标计入输入 QC；QC 不通过时仍保留探索性分数，`released_risk` 为空。QC 通过不代表临床有效性验证。

## 主要文件

| 文件 | 内容 |
|---|---|
| `manifest.json` | 输入、代码、依赖、预训练文件 hash 和运行配置 |
| `sample_qc.csv` / `cohort_audit.json` | 样本剔除与队列信息 |
| `feature_selection.csv` | build 内筛选分数和入选指标 |
| `inner_folds.csv` / `context.csv` | 内层标签隔离、最终上下文名单 |
| `learning_curve.csv` / `training_summary.json` | 微调损失、早停与权重更新证据 |
| `selected_weights.ckpt` / `model_bundle.joblib` | 选中权重、可复用完整预测模型 |
| `test_predictions.csv` | 冻结的逐人原始／校准风险和预测 Y，不含真实 Y |
| `test_individuals.csv` | 评估后追加真实 Y 和 IPCW；删失未知的 Y 留空 |
| `test_metrics.csv` / `test_calibration.csv` / `REPORT.md` | 测试指标和报告 |
| `DONE.json` / `FAILED.json` | 完成／失败状态 |

论文、模型的通用 benchmark 优势不等于在本队列有效；本地全量独立测试、跨 seed 稳定性和外部验证仍然需要实际运行检验。

## 本次验证（2026-09-23）

22 项回归测试通过；四个真实输入组合的路径、CUDA 和权重 checksum 预检通过。1,200 人合成数据完成蛋白与原始代谢物的 GPU 微调、冻结和测试评估；微调进行了 8 次反向更新，选中的模型有 248 个 tensor 与原始权重不同。代谢物检查包含仅出现在 test 的负丰度与亲缘组隔离。

固定拆分后改写 595 个测试样本的 Y，重训所得权重、筛选特征、context、阈值及冻结预测完全一致。重载模型对 599 个不含 Y 的 query 成功预测，4 个缺失率超标样本的 `released_risk` 为空。由于 CUDA FP16 和预测批次不同，595 个可比样本的校准风险最大差为 0.000382；需要更高数值精度时可在训练和预测中使用 `--no-amp`。验证记录位于 `/mnt/d/analysis/panome_checks/tf_20260923/validation.json`。这些是功能验证，尚未重跑真实队列的完整 TF 训练。
