# GRID — Genetic Risk based on Individual Distance

GRID 是一个新的研究实现：先把 PRS-CSx 的 SNP 后验权重按有来源的人类变异注释拆成遗传特征，再在训练数据里为每个待预测者寻找相近的参考个体，借用他们的交叉拟合残差，并输出匹配证据与预先冻结的筛选结果。

公共入口是 **`grid.sh`**。它替代旧 `3grid.sh` 的方法设计，不沿用旧 transportability 再收缩、频率推年龄或“古老变异较安全”的规则。研究目标是检验年龄信息和逐人匹配能否改善 prediction 与 explainability；能否超过 PRS-CSx、DiscoDivas 或 PRSformer，必须由同一批独立测试对象的结果决定。

## 1. 先运行起来

以下命令使用本项目已有的 `/mnt/d`、`/mnt/f` 默认目录。先准备好 PRS-CSx 四组分数、对应 SNP 权重及 1KG 空间中的目标样本投影；`grid.sh all` 不会自动重跑这些上游分析。

### 安装共享 CPU 环境

```bash
cd /mnt/d/scripts/grid
bash install_grid.sh
conda activate grid
bash grid.sh --python-help
```

安装器建立或更新 `grid` 环境，包含 Python 3.11、R、PLINK 2 和 GRID 的 CPU 依赖，不再安装或编译 ARG-Needle。`rdata==1.1.0` 的 RDS 读写要求 Python 3.11 或更新版本；旧 Python 3.10 环境需要更新，不能仅补装该包。`grid.sh` 启动 CPU 流程时会自动执行 `conda activate grid`，同时让 PLINK 2 进入 PATH；可用 `--python` 或 `GRID_PYTHON` 指定另一套已经具备依赖的 Python。若要保留原有环境，可用 `ENV_NAME=grid-evolution bash install_grid.sh` 安装到独立环境，再用 `conda activate grid-evolution` 激活，并传 `--python "$CONDA_PREFIX/bin/python"`。

PRSformer 的 GPU 环境独立安装，步骤见 [README.prsformer.md](README.prsformer.md)。当前适配器对应 PyTorch 2.6.0、CUDA 12.4 与 NATTEN 0.17.5；运行时不自动安装依赖。

### 下载真正的变异年龄注释

```bash
python f/grid.evolution.py download \
  --geva-dir /mnt/f/ref/GEVA --chromosomes 1-22
```

下载的是 [GEVA Atlas 的染色体级 summary](https://human.genome.dating/download/index)，GRCh37、年龄单位为 generations，全部压缩表约 5.5 GB。下载器核对官方 MD5、计算 SHA256、完整读取 gzip 后才发布文件；已验证的文件直接复用。数据留在参考目录，不放入 Git。已有 GEVA 目录可直接使用；提供符合第 4 节规范的 canonical 年龄表也可以替代 GEVA。

### 检查并运行三个性状

```bash
bash grid.sh --traits height,ldl,t2dm --check
bash grid.sh all --traits height,ldl,t2dm
```

`--check` 检查可用输入、样本对齐、划分与配置，不训练模型，不产生正式预测报告。逐 SNP 的坐标、等位基因和注释匹配在 `evolution` 阶段检查。只打印将执行的命令可用 `--dry-run`。

### 在相同测试半区运行 PRSformer

安装好独立 GPU 环境及官方 PRSformer checkout 后：

```bash
bash grid.sh all --traits height,ldl,t2dm --run-prsformer \
  --prsformer-python "$HOME/.venvs/grid-prsformer/bin/python" \
  --prsformer-device cuda:0
```

执行顺序是 `prepare → evolution → PRSformer → fit → report`。需要覆盖其训练参数时使用 `--prsformer-epochs 30` 等前缀参数；详见 `bash grid.sh --help`。GRID 会强制传入同一份 keep、family 和 split 表，避免 PRSformer 另分测试集。

## 2. 输入、性状和共同队列

### 默认数据位置

| 参数 | 默认路径或取值 | 用途 |
|---|---|---|
| `--pheno-file` | `/mnt/d/data/ukb/phe/Rdata/all.rds` | 表型与协变量 |
| `--score-dir` | `/mnt/d/data/ukb/pgs` | 各 trait 的 `1csx.scores.rds`，以及可选 `2disco.scores.rds` |
| `--pca-file` | `/mnt/d/data/ukb/pca_proj/ukb.discodivas.pca.tsv.gz` | 同一参考空间的个体匹配 PCs |
| `--ancestry-file` | `/mnt/d/data/ukb/pca_proj/ukb.ancestry.auto.tsv.gz` | 祖源分层和基线校准标签 |
| `--dir-gen` | `/mnt/f/gen/ukb/37/hap` | 染色体基因型 |
| `--gwas-dir` | `/mnt/f/gwas/4grid/common` | 已完成 CSx 推断的来源目录 |
| `--snpinfo` | `/mnt/f/refLD/csx/snpinfo_mult_1kg_hm3` | CSx SNP 信息 |
| `--geva-dir` | `/mnt/f/ref/GEVA` | 真实年龄 summary |
| `--covariates` | `age,sex,PC1,PC2` | 所有主比较一致的协变量 |
| `--distance-pcs` | `10` | 参考投影的匹配 PCs，表内命名为 `match_PC1…` |
| `--out-root` | `/mnt/d/analysis/grid/GRID` | 正式结果 |
| `--cache-dir` | `/tmp/grid-cache/grid` | 私有临时数据、阶段缓存和日志 |

`height`、`ldl` 使用指定连续表型列的原有单位；代码不自动进行 LDL 药物校正。`t2dm` 是**基线 0/1 结局**，不把事件时间或 incident indicator 当成二分类风险。默认 `--t2dm-col auto` 使用 `t2dm.Yr2e`、`t2dm.Yt2e`：基线病例为 1，明确基线非病例或具有有效随访状态者为 0，其余缺失。也可显式指定已经定义好的基线列。

准备阶段求三个性状输入、基因型、分数和 PCs 的共同样本交集，保留至少一个性状有结局的人，并冻结一个共同 outer roster。某性状缺失结局者不进入该性状的拟合和测试指标，缺失情况不用于选择高性能样本。各性状实际 n 与角色人数会单独报告。

### 家系与固定划分

默认以固定 seed `20260904` 做约 50/50 划分，不使用 Y 的数值。若有亲缘/家系连通分量，给出 `eid,family_id` 或 `eid,group` 表；同一分量必须留在同一半区和内部角色。

例如，先把实际家系表保存为 `/mnt/d/files/ukb.family.tsv.gz`，随后每个阶段使用同一参数：

```bash
bash grid.sh all --traits height,ldl,t2dm \
  --split-group-file /mnt/d/files/ukb.family.tsv.gz
```

没有家系表时，程序以单个 ID 为分组并记录这一事实；它不会因此宣称样本彼此无亲缘。已有固定划分可用 `--split-file`，字段 `eid,split`，split 为 `train` 或 `test`；仍需满足约 50/50 和家系隔离。家系大小与结局缺失会使比例略有偏差，代码会检查过大的偏离。

### 外部准备好的个体特征

`--data-file` 可以接入已有特征表，支持 RDS 或 TSV/gzip TSV，以及路径中的 `{trait}`。至少包含 `eid,split,y`、四个 `csx.EUR/AFR/EAS/SAS` 和指定协变量；匹配 PCs 使用 `match_PC1…`，演化模块使用 `evo.*`。可提供 `family_id,ancestry,disco`。

这是特征接口，不会自动证明外部年龄来源、SNP 对齐或无结局泄漏。该路径在结果中标记 `external_prepared_features`、`evolution_verified=False`；完整复现优先使用本项目的原生基因型评分路径。

## 3. 四个 Python 模块分别做什么

| 文件 | 责任 | 接触哪些结局信息 |
|---|---|---|
| `f/grid.data.py` | 统一 ID、性状、家系、共同划分；核对 CSx 权重；调用 PLINK 生成遗传模块 | 读取表型并确定是否有观测；不按 Y 值学习划分或注释 |
| `f/grid.evolution.py` | 合并外部年龄/选择/频率注释，输出逐染色体评分权重和 QC | 不读取任何个体 Y 或 train/test 标签 |
| `f/grid.abm.py` | 开发集内训练基线、监督距离、donor bank、gate 和审计；冻结预测 | 只学习开发半区结局；test Y 立即从输入视图移除 |
| `f/grid.py` | 调度阶段、核验缓存身份、保存模型、连接测试结局并评价 | 冻结预测后，在 report 阶段连接 test Y |

`grid.sh` 负责命令入口、共享参数、阶段顺序、并发锁及 PRSformer 对接。代码的分工对应不同数据和内存需求，不是四套重复方法。

## 4. 演化信息：输入什么，解释到哪里

### 人类群体历史的准确含义

这里研究的是人类群体内变异的年龄、分布和选择证据。现代 AFR、EUR、EAS、SAS 是同时代、内部仍有异质性的祖源概括，**现代 AFR 样本不是其他三组的直接祖先**；群体历史包含迁移、分化和混合。[人类起源原始研究](https://doi.org/10.1038/s41586-023-06055-y)支持这种网络式历史，而非四级演化阶梯。

GEVA 从单倍型共享、突变与重组信息推断年龄，不能只从某个 SNP 的变异频率推导年龄。频率和 LD 可以携带历史的结果，不能据此把频率差直接称为自然选择或某个起源时间。[GEVA 原文](https://doi.org/10.1371/journal.pbio.3000586)与 [Atlas FAQ](https://human.genome.dating/info/faq)说明其估计和假设。

年龄也不等于 pathogenicity。古老变异仍可影响现代常见病；例如祖先型 APOE ε4 与 LDL 升高的关系可见[原始人群研究](https://doi.org/10.1371/journal.pgen.1012285)。本实现不令 old SNP 自动变安全、不按年龄翻转 beta，也不额外强制衰减 beta。

GPN-Star、PrimateAI-3D、Evo 2 是启发来源，当前代码不调用这些跨物种模型。它们分别使用系统发育/比对、灵长类容忍变异与结构、跨生命域序列学习；相应任务上的成功不能直接当成 GRID 的 PRS 增益。原文见第 10 节。

### 年龄来自哪个等位基因

原生 GEVA 读取 `atlas.chr*.csv.gz`，优先 `Combined > TGP > SGDP`，使用 `AgeMode_Jnt`、相应区间与质量字段。坐标固定 GRCh37，按位置和完整等位基因对匹配，拒绝仅凭 rsID、互补链猜测或隐式 liftover。

GEVA Atlas 将 ALT 视作待定年的 derived allele。本实现只有在外部祖先碱基 `AlleleAnc == AlleleRef` 时才接受这一方向；未知或冲突的年龄留作原始审计值，评分分组进入 unknown。**derived allele 的起源年龄不一定是评分 A1 或 disease-risk allele 的年龄**，输出 `effect_allele_is_dated` 标明区别。

所有年代是依据遗传数据**推断**的 generations，不是直接观察到的历史年代、进入某个现代群体的年代或四组人群的分化日期。GEVA 的复合后验区间也不能当作完全校准的历史置信区间。

### 固定年龄模块与负对照

| 年龄模块 | 默认边界与要求 |
|---|---|
| `young` | `<1000` generations |
| `middle` | `1000–4000` generations，含 1000、不含 4000 |
| `old` | `>=4000` generations |
| `uncertain` | 有相应年龄信息但低质量、缺少界限、界限跨箱或点估计不在界限内 |
| `unknown` | 没有可用年龄或 derived/ancestral 方向不能确认 |

前三组要求质量达到 0.8，并且记录的界限落在同一个年龄箱且包含点估计。这些边界是预先固定的特征，不是人群分化时间。每个祖源的五组 beta 相加必须精确恢复该祖源原始 beta；同一批基因型得到的模块分数也必须重构对应总分。

`perm_age_*` 只置换 qualified 的 young/middle/old 标签，**每个 unknown/uncertain SNP 的位置、标签和权重保持不变**。置换在 CHR × reference MAF bin 内进行：`.01,.05,.10,.25` 为分箱边界；没有参考 AF 的位点单列一层。参考 MAF 由可用源群体 ALT AF 的等权均值折叠得到，与目标测试样本无关。

固定 seed 的负对照保留 SNP 身份、beta、模块数和每层年龄标签计数；输出可置换、固定及实际改变的 SNP 数。它是一次 prespecified negative control，不产生 permutation P value，也没有排除所有 LD、重组率和注释质量相关信息。

真正年龄与随机年龄预测臂具有相同的 frequency、selection、proxy 输入，只替换年龄标签。仅比较两条各有不同附加注释的模型，不能把增益归因于年龄。

### Canonical 注释和权重

原生主流程读取已完成的 CSx 结果；`--weights-file` 可替换为 TSV/gzip TSV，字段为：

```text
CHR BP SNP A1 A2 beta_EUR beta_AFR beta_EAS beta_SAS
```

A1 是评分等位基因。主流程要求四个 beta 列，允许某 SNP 缺少部分祖源 beta；不能一行全缺失。独立 builder 也支持只指定可用群体。仅接受 GRCh37、1-based、常染色体、双等位 A/C/G/T SNP，并检查重复 ID/位点。

`--annotation-file` 接受下表字段；至少有 `CHR BP REF ALT`，同一位点/等位基因对只能一行。可用 `{trait}` 为不同性状指定文件。

| 类别 | 字段 | 约定 |
|---|---|---|
| 坐标 | `CHR BP REF ALT`，可选 `BUILD` | GRCh37；精确等位基因匹配 |
| 年龄 | `ANC AGE_GEN AGE_LO AGE_HI AGE_QUAL` | generations、质量 0–1；ANC 必须为 REF |
| 年龄来源 | `AGE_SOURCE AGE_METHOD AGE_UNCERTAINTY` | 有限年龄必须有来源、方法及不确定性说明；未知界限要明确说明 |
| 选择 | `SEL_LOG10P_EUR/AFR/EAS/SAS SEL_SOURCE SEL_METHOD` | 已是 `log10(P)<=0`；仅提供有数据的群体 |
| 频率 | `AF_EUR AF_AFR AF_EAS AF_SAS FREQ_SOURCE` | ALT AF，范围 0–1；允许与年龄来自不同来源 |
| 可选代理 | `PROXY_<name> PROXY_SOURCE PROXY_METHOD` | 非负数值，需显式 `--include-proxy` |

Selection 采用来源方法报告的证据，默认按 `P<=0.001` 分 detected/not_detected/unknown；这一特征阈值不是经多重检验校正的全基因组发现，也不表示疾病作用方向。Relate 的原生 RData 需先按实际字段和等位基因转换；其 `pvalue` 已为 log10(P)，不要再次取负对数。

Frequency differentiation 使用至少两个参考群体的 `max(AF)-min(AF)`，按 `.05,.20` 分 low/middle/high，另有 unknown。它是频率分化特征，**不是 FST、变异年龄或选择检验**。代码不自动以目标测试集 AF 补参考注释。

缺少真实年龄时默认停止。只有显式 `--allow-proxy-only` 才允许 `proxy_only`；下游关闭真实年龄及随机年龄臂，主模型退回 frequency/proxy 或 no-evolution 分支。保留的 age_unknown 交换列只用于总分重构，不能据列名宣称已加入真实演化信息。

## 5. 从“4 列拟合”到逐人借用的实际算法

GRID 仍保留一个经校准的 CSx 基线 `b_i`，再加逐人匹配得到的修正；它不声称完全去掉 PRS 的全局校准。当前基线在样本足够的目标祖源中分别拟合四个 PRS 加协变量，数据不足或新祖源使用 pooled fallback；连续结局用 ridge，T2DM 用 logistic。

在 build 数据中按家系进行嵌套交叉拟合，donor `j` 的残差为：

$$
r_j=Y_j-\widehat b_{-\mathrm{fold}(j)}(X_j).
$$

该 donor 自己所在家系的 Y 不参与它的基线拟合或超参数选择。匹配距离的 ridge 监督权重使用 build 的这些 OOF 残差学习，再对 ancestry、CSx、frequency、evolution 各 block 标准化和压缩，最终最多 12 维；因此整个距离学习不是 outcome-blind，但只使用开发集。

KD-tree 默认精确检索；可选 HNSW 会在训练数据抽样核对 recall，低于阈值时停止。对 query 排除相同 ID 和家系，在最近 k 个候选中保留距离小于固定 caliper 的 donor。基础半径 `rho` 是训练参考对象到第 5 个可用邻居距离的中位数，使用最多 2048 个确定性抽样的 build 对象估计。

令有效近邻集合为 `N_i`，权重为：

$$
w_{ij}=\frac{\exp[-d_{ij}^2/(2\rho^2)]}{\sum_{l\in N_i}\exp[-d_{il}^2/(2\rho^2)]}.
$$

同一家系的权重先相加，以家系质量计算 `ESS_i=1/sum_f W_if^2`。实际借用比例与预测为：

$$
s_i=\exp[-d_{i,\min}^2/(2\rho^2)],\qquad
q_i=\alpha s_i\frac{ESS_i}{ESS_i+20},\qquad
\widehat Y_i=b_i+q_i\sum_{j\in N_i}w_{ij}r_j.
$$

少于 3 个匹配 donor、家系 ESS 小于 2 或没有可靠半径时，修正为 0；T2DM 最后截到 `[0,1]`，并单独记录 clipping 修正。默认在 tune_model 中比较 `k=16,32,64`、`alpha=0,.25,.5,.75,1`、半径倍数 `1,2,4`，始终允许 alpha=0 回退。

解释输出包含实际匹配 ID、最多 k 位候选中过 caliper 的人数、独立家系数、ESS、分块距离、每位 donor 的 OOF 基线/残差及加权贡献。这回答“哪些参考个体、以多大权重、使这次预测改变了多少”；不是统计学或生物学意义上的双胞胎，也不是因果机制。候选预测和最终 policy 都可以逐项重构。

## 6. 50/50 设计、预先筛选和公平对照

### 内部样本预算

同一个 outer test 半区在所有模型冻结后才用于评价。训练半区默认按以下角色分配，比例均以完整共同队列为分母：

| 角色 | 约占全队列 | 用途 |
|---|---:|---|
| build | 30% | 全局模型、监督距离及 donor bank；默认 5 折 donor OOF |
| tune_model | 7.5% | 选 ridge/HGB 参数和 k、alpha、caliper |
| tune_gate | 7.5% | 学习哪些人有更大预期增益、哪些人误差更低 |
| calibration_fit | 2.5% | 从已训练 gate 的分值确定固定 coverage 阈值，不用其 Y 再调模型 |
| calibration_audit | 2.5% | 对固定候选策略做独立内部审计 |
| test | 50% | 最终评价，不参与上述学习 |

只有 build 是 donor bank，所以匹配库约为总队列的 30%，不是整个训练半区。家系作为整体划分；实际人数及每种方法用过的独立结局数记录在 `grid.training.xlsx`。

PRSformer 在同一 outer development 半区内划为 **40% fit / 10% validation / 50% test**。`CSx_full_training`、`DiscoDivas_full_training` 使用开发期选定的超参数，在整个 **50% training** 重拟合；还提供全训练 Ridge/HGB 控制，避免把数据预算不足的基线当作主要对手。它们不回流修改已经完成独立审计的 matching policy。

### 两种筛选问题分别回答

Gate 的 gain 目标是 `(Y-CSx)^2-(Y-GRID)^2`；低误差目标是 `(Y-GRID)^2`。默认 `selected` 按预期 gain 选人，`selected_absolute_error` 按预期平方误差选人。前者表示可能相对内部 CSx 改善，后者表示相对容易预测。

阈值在 calibration_fit 中固定，默认目标 coverage=0.5；test 中只应用阈值，因此实际 coverage 可以不是 0.5。额外保留 `.2,.4,.5,.6,.8,1` 覆盖率，以及 random、support_only、absolute_error 控制；T2DM 另有 clinical_lowrisk 控制。

`candidate` 还要求预期 gain>0 且匹配支持足够。独立 calibration_audit 对候选人计算配对 MSE/Brier 差值，按家系 bootstrap；只有信息量足够且 95% 区间上界<0 才标记 supported_gain。随后 test 的 `released` 才可采用匹配候选，其余 `GRID_policy` 保留内部 CSx 基线。内部审计没有通过是一个真实结果，不能转而在 test 挑阈值。

这种筛选不保证某个人预测准确，也不保证 selected 子集 R² 高于全样本：较低 MSE 与较高 R² 是不同条件，子集 Y 方差缩小就可能降低 R²。

### 主比较与消融

| 模型/实验臂 | 检查的问题 |
|---|---|
| `CSx`、可选 `DiscoDivas` | build 校准的内部模型 |
| `CSx_full_training`、`DiscoDivas_full_training` | 充分使用同一开发半区的强基线 |
| `GRID_no_evolution` | 仅 ancestry + CSx 的逐人匹配 |
| `GRID_frequency_only` | 在共同匹配特征上加入频率分化 |
| `GRID_evolution` | 加入真实年龄及可用 selection/proxy，默认主 candidate |
| `GRID_permuted_evolution` | 相同附加信息，仅替换 qualified 年龄标签 |
| `Ridge_evolution`、`HGB_evolution` 及 `_full_training` | 同类特征是否仅靠全局回归/非线性拟合就足够 |
| `GRID_policy` | 预先冻结的实际输出策略，主报告对象 |
| `PRSformer` | 在同一测试集合上的独立官方架构训练适配器 |

主配对损失参考是 **CSx_full_training**。另输出 primary 对 no-evolution、frequency-only、permuted-age、全训练 Ridge/HGB 的直接差值，以及 primary/policy 对 PRSformer 的直接差值；不能用两条各自相对 CSx 的区间替代两模型直接比较。

所有方法在相同 all、gain-selected、rejected、low-error-selected/rejected 对象上比较，并分祖源报告。连续性状报告 MSE/RMSE/MAE、`total_R2`、`SSE_partial_R2` 以及 `Prediction_R2=cor(Y-covariate_baseline,prediction-covariate_baseline)^2`；T2DM 报告 Brier、log loss、AUC、average precision 和 case coverage。

这里的 `Prediction_R2` 是相对共同冻结的 build 协变量基线的描述性相关。完整 training 重拟合模型及 PRSformer 的增量中也包含协变量系数重拟合差异，不能把它解释为纯遗传增量，或直接等同于旧 Yeval、论文中的同名 R²。主要性能判断采用同一批人的配对 MSE/Brier。

`CSx_full_training` 的 full training 指目标人群组合/校准模型使用完整 development 半区；这里没有重新运行 CSx MCMC 或选择上游 phi。如果输入后验来自固定 phi，比较仍以那套后验为条件。要主张超过充分调优的 PRS-CSx，还应在 development 内完成上游参数选择，并核实 discovery GWAS 与目标 test 的样本重叠。

主结论首先看完整 test 的配对差值。区间是条件于已冻结模型的家系 bootstrap，不覆盖重新训练全部流程的不确定性。三个性状、多个祖源和 coverage 的探索结果不能自动当作多重检验后的显著发现。

PRSformer 比较还核对完整 outer roster、相同 test ID、家系隔离、结局数值/单位和协变量定义。缺少相同划分的结果时明确记录未运行；旧实验分数不自动混入。相同测试集也不代表外部 GWAS、个体训练量或预训练信息量完全相同，结果表保留实际方法预算。

## 7. 分阶段运行与模型复用

独立运行各阶段时，trait、输入、协变量、家系、seed 和模型参数保持一致：

```bash
bash grid.sh --stage prepare --traits height,ldl,t2dm
bash grid.sh --stage evolution --traits height,ldl,t2dm
bash grid.sh --stage fit --traits height,ldl,t2dm
bash grid.sh --stage report --traits height,ldl,t2dm
```

也可在 evolution 后运行 `bash grid.sh prsformer --traits height,ldl,t2dm --prsformer-device cuda:0`，再继续 fit/report。`prepare` 产生临时共同 keep、family 与 PRSformer split 表；它们必须与这次配置一起复用。改变输入或模型源码后，身份核验会拒绝陈旧缓存。确需重建时选新缓存/输出目录，或显式使用 GRID 的 `--replace`；PRSformer 自身恢复训练的约定见其 README。

### 新个体预测

`predict` 接收已经按训练方案计算好的**个体特征表**，不从原始基因型自动开始，也不需要 Y。表需含 `eid`、模型使用的协变量、`csx.*`、`match_PC*` 和 `evo.*`；有已知家系和祖源时一并提供。新人的 ID 不得出现在开发集。

必须沿用训练保存的 SNP 权重、A1、GRCh37 坐标、相同 PCA 投影空间、模块定义及**训练中心化频率**，不能用新人或测试样本重新估计中心。原生评分用 PLINK `center no-mean-imputation` 并只读取 SUM：缺失基因型贡献 0；各年龄分区继续重构同批评分总分。

原生流程发布的 `<trait>/scoring/` 保留可复用权重、模块和频率。外部 prepared-feature 路径没有这些原生评分材料，必须由特征提供者保留。模型中保存训练期数值转换、metric、donor bank 和 gate；不能只拿一个总 PRS 代替全部输入特征。

例如，将遵守这些条件的新 height 特征表保存为 `/mnt/d/data/ukb/pgs/height/grid.new_features.rds` 后：

```bash
bash grid.sh --stage predict --trait height \
  --model-file /mnt/d/analysis/grid/GRID/height/grid.model.rds \
  --data-file /mnt/d/data/ukb/pgs/height/grid.new_features.rds
```

结果写入该 trait 的 `projection/`，包括 `grid.predictions.rds` 和 `grid.matches.rds`。没有提供 family_id 时以新 ID 自身分组；提供了与开发集相同的真实家系时会排除该家系 donor，并标记 overlap，不能把亲属预测当作无亲缘外部验证。

## 8. 正式结果和临时材料

默认 `/mnt/d/analysis/grid/GRID/` 下有共同 `grid.split.rds`，每个 trait 一目录。

| 文件 | 内容 |
|---|---|
| `grid.performance.png/.xlsx` | 完整 test、同一筛选子集、祖源分层、主配对损失和 coverage |
| `grid.contrasts.png/.xlsx` | 对 PRSformer 和各消融/全局控制的直接配对比较 |
| `grid.coverage.png/.xlsx` | 冻结阈值下的实际 coverage 与预测误差；不同 selector 控制 |
| `grid.support.png/.xlsx` | 实际匹配人数与有效家系支持 |
| `grid.training.xlsx` | 内部角色、调参、donor OOF、检索质量、metric、标签预算和独立审计 |
| `grid.model.rds` | 可复用冻结模型和最少必要 provenance；Python 对象以无损 payload 保留 |
| `grid.test_individuals.rds` | 测试对象、结局、所有模型预测及筛选诊断 |
| `grid.scores.rds` | 主要分数和 selected/released 标记 |
| `grid.individual_explanations.rds` | 逐人匹配支持、预测修正、策略来源及实际匹配输入特征 |
| `grid.matches.rds` | packed donor 证据；0-based index，-1 表示无匹配 |
| `grid.training_roles.rds`、`grid.donor_residuals.rds` | 开发期角色，以及 donor 的交叉拟合残差和匹配输入特征 |
| `scoring/` | 原生流程可复用的 SNP/模块权重与训练频率 |

每张正式 PNG 对应同目录、同名 XLSX，工作簿只保留该主题的汇总结果。带个体 ID 的结果与匹配证据保存为 RDS，不混入汇总工作簿或公开仓库。交换 TSV、PLINK 工作目录、临时验证和日志留在 `/tmp/grid-cache/`；正式模型与评分材料应保留用于复现和新人预测。

PRSformer 的共同划分比较默认单独写入 `GRID/benchmark/prsformer/`，个体分数位于 `GRID/benchmark/scores/<trait>/3.prsformer.scores.rds`。用 `--prsformer-root` 可整体更换这一位置。

## 9. 保留的上游入口和旧 Yeval

```bash
bash 0.pca.sh --check
bash 1.csx.sh --traits height,ldl,t2dm --check
bash 2.disco.sh --traits height,ldl,t2dm --check
```

检查通过且确需计算时去掉 `--check`；同样可以使用 `grid.sh pca/csx/disco` 转发。`0.pca.sh` 把目标样本投影到已有 1KG PCA 空间；`1.csx.sh` 运行已有 PRS-CSx 和评分流程；`2.disco.sh` 生成基于该空间的 Disco 分数。新 GRID 复用这些产物，并在自己的开发半区重新校准。

`1.csx.sh` 的 MCMC、posterior、auto/meta 和内存限制仍由其入口负责；新 GRID 不要求另生成逐人 posterior covariance，也不将 beta shrinkage 当作变异年龄。GWAS 效应等位基因频率 EAF 必须与 A1 对应，不能把 MAF 无条件充作 EAF。已有原生 CSx SNP 权重和 MCMC 数据保留供复用。

旧评估仍可独立调用，例如：

```bash
bash Yeval.sh --trait height --type ct --covar-name age,sex,PC1,PC2
```

**旧 Yeval 的 OOF 预测、旧 10/90 或其他分组结果不能混入新 GRID 的 50/50 测试评价。** 旧 T2DM 生存结局也不能与这里的 baseline 0/1 Brier/AUC 比较。只有同一测试 ID、结局、协变量、评分信息和预算清楚的结果，才构成这里的方法比较。

## 10. 验证状态、原始文献与作者资源

本地整合的修复、运行环境和 10 项可复跑测试见 [VALIDATION.grid.md](VALIDATION.grid.md)；下方为下载包原有验证记录，不能视为本机重新执行过 PRSformer 训练。本机整合保留 `0.pca.sh`、`1.csx.sh`、`2.disco.sh` 的点号命名及已有 CSx 缓存锁修复。旧 `f/3grid.f.sh`、`f/3grid.py` 已移出运行目录，`f/0.arg.py` 继续保留供 GU 使用。

本实现提供研究流程和可审计输出；本代码编写会话没有进行真实 UKB 的完整训练，也没有产生“已经超过论文”的实证结论。已有 height 报告约 0.3 的 R² 或 Disco 无提升，不能直接当成新 GRID 的基线值；首先要核对 R² 定义、测试祖源组成、GWAS 重叠、样本量、校准与划分。

演化 builder 已用合成 SNP/注释完成 23 项检查，包含真实 GEVA 格式解析、等位基因方向、分数重构、固定未知位置的负对照、来源和 gzip/MD5 异常。独立 ABM 小例验证了 test Y 改动不影响拟合预测、无 Y 新预测、非零 donor 贡献与 policy/距离分解重构。这些验证不代替真实队列外部验证。

最终端到端联调使用 600 人、300 个家系、2 条染色体、24 个 SNP 的模拟数据和真实 PLINK 2/PGEN，同时运行三个性状。共同名单为 300/300；官方 PRSformer 小模型实际训练 2 epochs，在 240/60/300 名单上完成前向、反向及预测，经过 QC 保留 23 SNP，三个性状均通过共同 test、结局和协变量核对。这是 CPU global-attention 的接口验证，不是论文规模或 CUDA/NATTEN 性能复现。

从正式个体 RDS 独立重算 1,200 行指标和 1,575 行配对比较，最大绝对差为 `8.88e-16`。12 张 PNG 均有可读同名 XLSX，15 本汇总工作簿未混入个体 ID；三个 model RDS 重载后均可对无 Y 新人重现预测。27 份永久评分材料与模型内 SHA256 一致，用保存的权重和 training AF 给 300 个新 ID 重打 80 个模块，数值与原分数一致。输入来源变化、同大小权重篡改、错用 trait 模型、PRSformer 测试家系重划均被拒绝。

该小型联调中三个性状的 audit 均未放行个体修正，`released=0`，策略正确回退 CSx；独立匹配核测试另外覆盖了非零借用与精确贡献分解。未实际运行 HNSW 可选路径、完整 UKB 训练或全基因组 GPU PRSformer；默认 KD-tree 路径已验证。

| 原始工作 | DOI / 作者资源 | 与本实现的关系 |
|---|---|---|
| PRS-CS，2019 | [10.1038/s41467-019-09718-5](https://doi.org/10.1038/s41467-019-09718-5)；[作者代码](https://github.com/getian107/PRScs) | 连续收缩的多基因预测 |
| PRS-CSx，2022 | [10.1038/s41588-022-01054-7](https://doi.org/10.1038/s41588-022-01054-7)；[作者代码](https://github.com/getian107/PRScsx) | 多祖源后验权重及主要比较基线 |
| DiscoDivas，AJHG 2026 | [10.1016/j.ajhg.2026.05.006](https://doi.org/10.1016/j.ajhg.2026.05.006)；[作者代码](https://github.com/YunfengRuan/DiscoDivas) | 祖源连续空间中的 PRS 插值启发；不是局部祖源推断 |
| PRSformer，NeurIPS 2025 | [会议原文](https://proceedings.nips.cc/paper_files/paper/2025/hash/b9f2c7f6cc690434047ec546c83270dc-Abstract-Conference.html)；[预印本 DOI 10.1101/2025.10.26.684578](https://doi.org/10.1101/2025.10.26.684578)；[作者代码](https://github.com/23andMe/PRSformer) | 原始个体基因型的多任务监督模型；本地三个性状是训练适配，不复现原论文私有大队列 |
| GEVA，2020 | [10.1371/journal.pbio.3000586](https://doi.org/10.1371/journal.pbio.3000586)；[Atlas](https://human.genome.dating/)；[代码](https://github.com/pkalbers/geva) | 当前可直接接入的 GRCh37 变异年龄 |
| Relate，2019 | [10.1038/s41588-019-0484-x](https://doi.org/10.1038/s41588-019-0484-x)；[公开人类年龄/选择数据](https://doi.org/10.5281/zenodo.3234689) | 可转换为 canonical 注释；分支上下界不等于 95% CI |
| CLUES，2019 | [10.1371/journal.pgen.1008384](https://doi.org/10.1371/journal.pgen.1008384)；[代码](https://github.com/standard-aaron/clues) | 基于 genealogy 的选择/轨迹推断；不是代码内置现成的全基因组年龄表 |
| 1000 Genomes Phase 3 | [官方数据资源](https://www.internationalgenome.org/data-portal/data-collection/phase3/) | 祖源 PCs 和可选参考 AF；AF/FST 不直接等于年龄或选择 |
| GPN-Star，Nature 2026 | [10.1038/s41586-026-11005-5](https://doi.org/10.1038/s41586-026-11005-5)；[代码](https://github.com/songlab-cal/gpn) | 系统发育与比对建模的启发，不是当前调用模块 |
| PrimateAI-3D，Science 2023 | [10.1126/science.abn8197](https://doi.org/10.1126/science.abn8197)；[原文](https://pmc.ncbi.nlm.nih.gov/articles/PMC10713091/) | 灵长类容忍变异是概率证据，不是对人类无害的保证 |
| Evo 2，Nature 2026 | [10.1038/s41586-026-10176-5](https://doi.org/10.1038/s41586-026-10176-5) | 跨生命域序列规律，不是 AFR→EUR 的定年器 |
| SBayesRC，Nature Genetics 2024 | [10.1038/s41588-024-01704-y](https://doi.org/10.1038/s41588-024-01704-y) | 注释参与 PRS 已有先例；GRID 的研究点在个体匹配、证据分解和筛选验证 |

方法的新意应由消融与独立测试支持；公开已有注释参与预测，并不自动构成新的生物学机制发现。
