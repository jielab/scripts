# 在 grid 中运行 PRSformer：height、LDL、T2DM

`3.prsformer.sh` 把 UKB 个体基因型和表型接入 **23andMe 官方 PRSformer 网络**，训练一个含 height、ldl、t2dm 三个输出的模型，并在预先隔离的 test 样本中评分。PRSformer 软件由 **23andMe, Inc.** 开发；这个 grid 数据准备、训练和结果发布适配器独立编写。

PRSformer 的输入是个体水平的 genotype 和 phenotype。这里不读取 CSx 的四套 posterior SNP effects，也没有可直接下载应用的预训练 height/LDL/T2DM 权重。保留了官方网络的 genotype encoder、attention 和 phenotype head；为这三个任务新增训练流程、缺失表型 mask、协变量基线、连续/二分类损失及验证集选择。**这是在你的数据上训练官方架构的适配实现，不是论文原始训练实验的数值复现。**

## 1. 文件和运行环境

把以下文件放到现有 scripts 仓库的 `grid` 目录，保持相对位置：

| 文件 | 职责 |
|---|---|
| `3.prsformer.sh` | 路径、阶段、模型参数、依赖预检及任务编排 |
| `f/3.prsformer_data.py` | 对齐样本/SNP、划分数据、训练集 QC、生成矩阵及发布报告/RDS |
| `f/3.prsformer.py` | 调用官方架构，训练、验证选择、checkpoint、test 预测和评估 |
| `README.prsformer.md` | 本说明 |

脚本用 `--python` 指定解释器，不读取会切换现有 grid Python 环境的 `f/0.common.sh`，也不自动安装软件。

### 单独建立 Python 3.11 环境

以下命令适用于有兼容 NVIDIA 驱动的 Linux/WSL。PyTorch/NATTEN 的版本组合必须匹配；这里锁定官方 PRSformer 使用的旧 NATTEN API，对应安装组合见 [NATTEN 官方安装页的 0.17.5 部分](https://natten.org/install/#0175)。

```bash
python3.11 -m venv "$HOME/.venvs/grid-prsformer"
"$HOME/.venvs/grid-prsformer/bin/python" -m pip install --upgrade pip
"$HOME/.venvs/grid-prsformer/bin/python" -m pip install \
  numpy pandas scipy pgenlib pyreadr matplotlib openpyxl scikit-learn
"$HOME/.venvs/grid-prsformer/bin/python" -m pip install \
  torch==2.6.0 --index-url https://download.pytorch.org/whl/cu124
"$HOME/.venvs/grid-prsformer/bin/python" -m pip install \
  natten==0.17.5+torch260cu124 -f https://whl.natten.org
"$HOME/.venvs/grid-prsformer/bin/python" -m pip check
```

也可以用自己的环境，并在每次调用时传 `--python /你的环境/bin/python`，或设置 `GRID_PRSFORMER_PYTHON`。正式的 neighborhood attention 路径要求 CUDA；找不到 GPU、torch/CUDA/NATTEN 不匹配或缺少 native attention API 时会退出。准备数据和发布报告两种阶段不要求 GPU。

### 获取官方源码

官方仓库的 [LICENSE.txt](https://github.com/23andMe/PRSformer/blob/main/LICENSE.txt) 包含内部非商业研究和再分发条款；本适配包提供到原作者的链接，运行时读取你从原作者取得的源码。

```bash
mkdir -p /mnt/f/software
git clone https://github.com/23andMe/PRSformer.git /mnt/f/software/PRSformer
git -C /mnt/f/software/PRSformer checkout 7dea4e1bb27975885c937f2be82f9243bc82bda7
```

已有 checkout 时只需核对版本，勿再次 clone 到非空目录。路径可用 `--upstream-dir` 或 `PRSFORMER_UPSTREAM_DIR` 更改。源码中无用的 `pytorch_lightning` import 在加载时经过检查后跳过，不改写磁盘上的官方文件；模型参数和训练 checkpoint 记录官方源码摘要以校验重载一致性。

## 2. 默认输入及表型定义

| 参数 | 默认值/内容 |
|---|---|
| `--dir-gen` | `/mnt/f/gen/ukb/37/hap`，`chr1` 至 `chr22` 的 PGEN/PVAR/PSAM；支持压缩 `.pvar.zst` |
| `--pheno-file` | `/mnt/d/data/ukb/phe/Rdata/all.rds`，含 `eid`、表型和数值协变量 |
| `--ancestry-file` | `/mnt/d/data/ukb/pca_proj/ukb.ancestry.auto.tsv.gz` |
| `--group-col` | `genetic_ancestry`；用于分层划分和各祖源结果展示 |
| `--snp-list` | `/mnt/f/refLD/csx/snpinfo_mult_1kg_hm3`；限制到既定 HM3 候选 SNP 集合 |
| `--covariates` | `age,sex,PC1,PC2`，与当前 grid 比较的默认协变量一致 |
| `--remove` | `/mnt/d/files/ukb.exclude.id`；文件缺失会报错，显式 `none` 可关闭 |
| `--keep` | `none`；可传仅保留样本的 ID 文件 |

`eid` 必须唯一。表型、祖源表和各染色体基因型按 ID 对齐，不假定输入行序相同；只保留所有选定染色体中都存在、协变量齐全且至少一个任务有表型的人。单个任务缺失的标签不计入该任务 loss。祖源字段已存在于 phenotype 文件时，直接使用该字段；否则与 ancestry 文件合并。

height 和 ldl 默认读取同名列，可用 `--height-col`、`--ldl-col` 指定其他列。数据准备不自动进行 LDL 用药校正、逆正态变换或临床异常值清理；输入列应当已经是你打算比较的表型定义。连续任务在训练样本内拟合协变量基线，把基线残差按训练集残差标准差缩放后训练，输出时恢复原表型单位。

T2DM 使用**基线患病状态**，不是随访发病的生存模型。`--t2dm-col auto` 沿用现有字段的时间含义：

- `t2dm.Yr2e == 1`：基线病例，标签为 1。
- 其余样本中，`t2dm.Yr2e == 0` 或 `t2dm.Yt2e` 为 0/1：基线对照，标签为 0；随访新发病例在基线仍属于对照。
- 无法据此确定者：缺失标签。

如果另有经过核查的基线 0/1 表型列，优先用 `--t2dm-col 你的列名` 明确指定。程序先用 train 拟合协变量 logistic 基线，再把神经网络输出作为 log-odds 增量，以 binary cross entropy 训练。它不使用发病时间、删失时间，也不产生 Cox HR。

如果希望与论文的 20 PCs 调整更接近，可传完整列表；所有对照方法也应使用相同列表和样本：

```bash
--covariates age,sex,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10,PC11,PC12,PC13,PC14,PC15,PC16,PC17,PC18,PC19,PC20
```

协变量必须是数值列，分类变量应事先明确编码。`--covariates none` 可关闭协变量基线。

## 3. 数据划分和泄漏控制

默认在祖源与 T2DM 标签的分层内，按固定 seed `20260904` 划分约 60% train、20% validation、20% test。小层的实际人数会受取整影响。每个任务的三个 split 必须有足够有效标签；T2DM 各 split 必须同时有病例和对照。

SNP 的 MAF 和 call rate **仅用 train 个体计算**，默认阈值分别为 0.01 和 0.98。这是整个 train cohort 的频率阈值。保留的 SNP 按数值染色体、位置排序；候选清单若含坐标，会核对 SNP ID 与 CHR:BP 的一致性。保留 ALT 剂量 0–2，以 `-1` 表示缺失，避免把等位基因计数误当成已经标准化的 GWAS effect。

协变量拟合、残差尺度估计和网络权重学习均只使用 train；validation 决定 early stopping 和最优 epoch。test 不参与这些拟合及模型选择。多任务 loss 对各任务有效观察的平均 loss 等权，而不是让缺失更少的表型自然获得更多权重。

默认按个体随机分割**不保证亲属分离**。分析前已排除相关个体时可直接使用；保留亲属时，传入 `--split-group-file family_groups.tsv`，文件列为 `eid,group`，每个 kinship connected component 使用相同 group，单独个体使用自己的 group。程序按 group 划分，并检查同一 group 没有跨 split。此模式需要 scikit-learn，且祖源/病例比例未必精确保持。

若已有固定划分，使用 `--split-file splits.tsv`，其两列是 `eid` 和 `split`；split 只能是 `train`、`validation`、`test`，并须覆盖所有保留个体。也可同时提供 family group 文件来验证亲属隔离。**公平比较 CSx、DiscoDivas 和 PRSformer 时，应固定这同一批最终 test IDs。**

## 4. 推荐运行方式

### 与新 GRID 使用共同测试集

先按本文配置独立的 PRSformer 环境及官方源码，再从 `grid.sh all --run-prsformer` 启动共同划分比较，完整选项见 [README.md](README.md)。这条路径使用 GRID 冻结的外层 50/50 名单，在 training half 内划分 PRSformer 的 train/validation，整体约为 40%/10%/50%。下方直接调用 `3.prsformer.sh` 的默认 60%/20%/20% 仍适用于独立实验。

导出的 scores RDS 额外保存 `covariate_names`、`endpoint_definition` 和 `preparation_signature`，供共同测试集比较时核对来源。

### 先做输入和运行环境检查

```bash
cd /mnt/d/scripts/grid
bash 3.prsformer.sh --traits height,ldl,t2dm --check
```

`--check` 的 all/prepare 模式先检查标签、样本 ID、split 和候选 SNP 元数据，并报告预计矩阵大小；尚未解码基因型时不能声称已完成 train-only MAF/call-rate QC。all 模式还检查生产训练依赖；若 cache 已准备好，会进行使用真实官方模型的小型 forward/backward 检查。这个模型检查验证接口与 kernel 可调用性，不保证全规模模型能装入 GPU。

在没有 GPU 的节点上只检查输入时使用 `--mode prepare --check`。仅查看将执行哪些命令，用 `--dry-run`，它不读取数据，也不验证依赖兼容性。

### 正式训练三个任务

先选一块容量充足的本地 SSD 放 cache，下面的路径按你的磁盘修改：

```bash
bash 3.prsformer.sh \
  --traits height,ldl,t2dm \
  --cache-dir /mnt/f/prsformer-cache/three_traits \
  --device cuda:0
```

默认 `all` 一次完成 prepare → train → report。`train` 在最优 validation epoch 选定后，已经生成 test 预测，不需要再调用一次 predict。也可以分阶段执行，`--stage` 是 `--mode` 的别名：

```bash
bash 3.prsformer.sh --mode prepare \
  --cache-dir /mnt/f/prsformer-cache/three_traits
bash 3.prsformer.sh --mode train \
  --cache-dir /mnt/f/prsformer-cache/three_traits --device cuda:0
bash 3.prsformer.sh --mode report \
  --cache-dir /mnt/f/prsformer-cache/three_traits
```

如果指定了不同的 trait 列表、协变量、输入、split 或模型参数，分阶段运行时保留相应选项。已有缓存会校验输入/配置和已保存数据清单；不匹配时要求新 cache 或显式 `--replace`。

每次实际运行都会对规范化后的 cache 路径取得 `flock` 锁，并保持到整个 wrapper 退出；同一个 cache 上的 prepare、train、predict、report 不能同时运行，避免 `prepare --replace` 覆盖正在训练的数据。不同路径指向同一 cache 时也使用同一把锁。`--check` 和 `--dry-run` 不取得写入锁。

### 很小的 CPU 流程检查

只有显式使用 `--attention global --device cpu --amp off` 才能在 CPU 上检查流程；此分支最多允许 4,096 SNP。以下运行三个任务的少量 SNP/样本示例，结果与正式分析分开保存：

```bash
bash 3.prsformer.sh \
  --chrs 22 --max-samples 2000 --max-variants 256 \
  --device cpu --attention global --amp off \
  --embed-dim 8 --heads 2 --layers 1 --ff-dim 16 \
  --batch-size 4 --accumulation-steps 4 --epochs 2 --patience 2 \
  --cache-dir /tmp/grid-prsformer-smoke \
  --output-root /mnt/d/analysis/grid/prsformer-smoke \
  --score-dir /mnt/d/data/ukb/pgs/prsformer-smoke
```

小样本若某个 T2DM split 没有病例，会在预检阶段报错。这个检查分支使用官方 global attention，验证数据和训练/发布流程；它不构成正式 neighborhood attention 模型的性能测试。

### 中断恢复和重算结果

训练中断时，用**原来的全部选项**加 `--mode train --resume`。输入、split、seed、源码摘要及训练参数必须一致；训练 epoch 总数也属于检查项。正常训练完成后会删除可恢复的临时 `training.pt`，保留 validation-selected `model.pt`。

```bash
bash 3.prsformer.sh --mode train --resume \
  --cache-dir /mnt/f/prsformer-cache/three_traits

# 对已有模型重算同一 cache 中的 test 预测
bash 3.prsformer.sh --mode predict \
  --cache-dir /mnt/f/prsformer-cache/three_traits

# 重新发布正式结果，不重新训练
bash 3.prsformer.sh --mode report --replace \
  --cache-dir /mnt/f/prsformer-cache/three_traits
```

`--checkpoint` 可指定 predict 使用的模型；预测使用 checkpoint 的架构和训练时预处理，并核对 SNP 顺序/REF/ALT 及官方源码。`--replace` 允许重建 preparation cache 或覆盖发布结果，**不删除已有模型**。重训另一配置时选择新的 `--run-dir`；`--resume` 与 `--replace` 不能同时使用。

## 5. 资源规模

最终 float16 genotype matrix 约占 `N × V × 2` 字节。准备阶段从 variant-major 转为 sample-major，要同时保留两份布局，峰值约为：

`2 × N × V × 2` 字节，另加至少 512 MiB 预留，以及现有缓存/模型所占空间。

例如 400,000 人 × 1,000,000 SNP，最终矩阵约 0.8 TB，准备期间约需 1.6 TB 空间。缓存默认位于 `/tmp/grid-prsformer`；不要在磁盘容量不足时直接启动全量解码。脚本会先检查可用空间，按 `--chunk-variants` 控制准备阶段的工作 buffer，以 memmap 存储矩阵。

GPU 需求也随 SNP 数线性增长，attention 的计算量与窗口相关；官方表型 head 同样包含与 SNP 数相关的参数。默认 embedding 64、2 层、窗口 385、micro-batch 1、梯度累积 64，并开启 FP16 和 gradient checkpointing。模型预检报告参数量、参数/梯度/Adam states 的最低显存和一份 activation 的大小；这些是分项估算，全训练峰值还有 attention workspace 等开销。小 batch 能减少激活显存，不能消除模型和 optimizer 的固定开销。

先根据实际 GPU 与候选 SNP 数确认是否可行，再决定正式规模。`--max-samples`/`--max-variants` 只用于显式的小规模开发运行；改变正式 SNP panel 时，应使用预先定义的候选集合，避免用完整 cohort 的表型筛选 SNP 后再声称 test 独立。

## 6. 正式输出与指标含义

默认正式结果如下：

| 文件 | 内容 |
|---|---|
| `/mnt/d/analysis/grid/prsformer/3.prsformer.pt` | 最优模型、训练时预处理、模型配置和输入/官方源码身份信息 |
| `/mnt/d/analysis/grid/prsformer/3.prsformer.split.rds` | 每人的 `eid,target,split`，提供分组文件时还含 group |
| `/mnt/d/analysis/grid/prsformer/3.prsformer.performance.xlsx`、同名 `.png` | 各任务、ALL 及各祖源的 test 评估指标和样本数 |
| `/mnt/d/analysis/grid/prsformer/3.prsformer.training.xlsx`、同名 `.png` | 每 epoch 的 train/validation loss 及最佳 epoch 信息 |
| `/mnt/d/data/ukb/pgs/<trait>/3.prsformer.scores.rds` | 该任务 test 个体的结果；包含 eid、target、split、outcome、baseline、prediction、prsformer 和 prsformer_scale |

`prsformer` 是模型的遗传预测分量：height/LDL 使用原表型单位，`prediction = baseline + prsformer`；T2DM 使用 log-odds 增量，`prediction` 是把基线 logit 与此增量相加后得到的概率。

临时的 `genotypes.npy`、`data.tsv.gz`、`variants.tsv.gz`、`prepare.json` 留在 cache；`history.tsv`、`metrics.tsv`、`test_predictions.tsv.gz` 和恢复 checkpoint 留在 run 目录。正式个体输出采用 RDS，正式汇总采用 XLSX/PNG，不把临时 JSON/个体 TSV 作为分析结果发布。默认 `/tmp` 可能被清理；要继续训练或原样重算，应保留 cache/run，正式模型和 split 文件另有发布副本。

### 连续表型的几种 R²

设 `b` 是仅在 train 拟合的协变量基线，`g` 是神经网络输出的遗传分量，所有下式在 test 样本上计算：

| 指标 | 数学含义 |
|---|---|
| `prediction_R2` | `cor(Y - b, g)^2`，残差与遗传分量的相关平方，忽略预测尺度校准误差 |
| `SSE_partial_R2` | `1 - sum((Y - b - g)^2) / sum((Y - b)^2)`，固定模型相对协变量基线减少的残差平方和比例 |
| `full_R2` | `1 - sum((Y - b - g)^2) / sum((Y - mean(Y))^2)` |
| `baseline_R2` | 用上式的 `b` 替代 `b + g` |
| `full_RMSE`、`prediction_bias` | 绝对预测误差及平均偏差，保留原表型单位 |

固定模型的 SSE 指标可为负值，这反映其在 test 上比基线差；不截断到零。这里的 `SSE_partial_R2` 不在 test 上重新拟合 PRS 系数，因此不能未经说明等同于论文在 test 上再拟合 `Y ~ covariates + PRS` 得到的 partial R²。比较时应同时匹配具体指标公式、协变量、表型变换和测试人群。

T2DM 报告 AUC、baseline AUC、ΔAUC、AUPRC、Brier、log loss、病例数及患病比例，不用连续表型的相关平方代替疾病预测评价。小祖源样本出现单一类别时，相应 AUC/AUPRC 是缺失值。图表给出点估计；它们本身不提供增益显著性的检验。

### 与现有 Yeval 的衔接

`3.prsformer.scores.rds` **只含固定 test 样本**，不是全 cohort 的 out-of-fold PRS。现有 Yeval 如果重新把样本随机分折、拟合组合系数或校准系数，其任务和本报告的固定模型评估就发生了变化。不要把这份 test-only 输出当作已有全样本 OOF 输入。

公平的比较流程是：预先固定 train/validation/test，CSx 组合权重和 DiscoDivas anchors 只在相应训练/验证数据中确定，在同一 test ID 交集上计算三个方法的指标。若需要“所有人都有 OOF PRS”，必须在外层各折重新完成该折的 genotype QC、预处理、模型训练、validation 选择与 test 评分，随后合并各外层 test 预测。

## 7. 来源和适用范围

- [23andMe/PRSformer 官方仓库](https://github.com/23andMe/PRSformer)，本适配验证的源码版本为 `7dea4e1bb27975885c937f2be82f9243bc82bda7`。
- [官方 PRSformer 模型源码](https://github.com/23andMe/PRSformer/blob/7dea4e1bb27975885c937f2be82f9243bc82bda7/src/model.py) 与 [LICENSE.txt](https://github.com/23andMe/PRSformer/blob/7dea4e1bb27975885c937f2be82f9243bc82bda7/LICENSE.txt)。
- [NATTEN 官方安装与旧版本 wheel 说明](https://natten.org/install/)。

适配器可以在具备上述数据和依赖的环境中启动训练；交付代码的流程检查不等于已经在你的 UKB genotype 上训练完成。真实 GPU 的 full-scale neighborhood attention 吞吐、显存与预测结果，须以你的运行日志和独立 test 输出为准。
