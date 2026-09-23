# GRID — PRS-CSx posterior evaluation

基于 `jielab/scripts@584b4c023c62acd11f67b8ddf8ade321108ea165`。保留最新 GRID/ARG 实验模块，扩展 PRS-CSx 后验评分并重写 Yeval 指标、四面板。没有引入 LDpred2。

## 改动

- height/LDL 主指标改为 **Prediction R²**：`cor(y − baseline_OOF, full_OOF − baseline_OOF)^2`。协变量调整、标准化、四分数组合都在训练折完成。旧 SSE-based partial R² 另存 `SSE_partial_R2`，不是给旧数字换名。
- T2DM 生存模型用 Harrell C；四面板 b 和经验距离图突出 **ΔC**。二分类模式以 AUC 为主。
- a：原始祖源标签的 PCA；b：实测群体预测；c：四个训练中心、代表性个体及距离连线；d：真实逐人后验点。
- d 不使用分箱、复制或抖动的点。经验距离分箱放到单独辅助图。
- SNP 后验抽样在每条染色体内对齐四个祖源；逐人保存均值与完整协方差。组合方差为 `w'Σ_i w`，包括交叉协方差。染色体间按 PRS-CSx 的独立模型相加，不按抽样序号制造关联。

**方法边界：** 堆叠 PRS-CSx 的个体可靠性是条件于拟合组合系数和外部遗传方差的模型估计，是 Nature 方法的探索性扩展，尚不代表已经校准的个体准确度，也不等于表型 Prediction R²。T2DM 没有“单个人的 C-index”；d 默认展示 PRS log-hazard 的后验 SD。

## 安装和重跑

将完整 `grid` 替换到 `/mnt/d/scripts/grid`，避免混用内部模块。根目录没有重复的旧 csx/disco/pca 入口；`f/csx.sh`、`f/disco.sh`、`f/pca.sh` 是必要的内部模块。

```bash
cd /mnt/d/scripts/grid
chmod +x *.sh
conda env update -n grid -f environment.yml
```

- **0pca.sh 本次没有修改，无需因此重跑。** 训练中心必须在相同坐标系；如果改用训练样本自己的 PCA 基底，则目标必须重新投影。
- **需要重跑 1csx.sh。** 旧平均权重不能还原后验抽样；新代码签名防止复用旧推断。
- 可以复用已保存的 Disco/COJO 分数作为比较输入；Yeval 必须重跑。
- 本包没有 UKB 数据，也没有声称已重算的真实结果。

## 1. 生成个体后验

```bash
./1csx.sh --traits height,ldl,t2dm --check
./1csx.sh --traits height,ldl,t2dm --posterior TRUE --jobs 4 --threads 4
```

`--posterior TRUE` 是默认值。默认 4000 次迭代、burn-in 2000、thin 5，共 400 个保留抽样。后验评分按染色体顺序执行，四个人群分别评分；每个 PLINK 默认 8192 MB，可用 `--posterior-memory` 调整。比平均权重评分耗时明显增加。

GWAS 必须提供正确的效应等位基因频率 **EAF**。已移除 MAF 自动映射 EAF 的错误别名。缺失或等位基因不明确会报错，不用 1KG 频率静默填补。可用外部 discovery EAF：

```bash
./1csx.sh --trait height --posterior TRUE \
  --posterior-frequency-dir /path/to/height_discovery_eaf
```

目录包含 `AFR.tsv.gz EAS.tsv.gz EUR.tsv.gz SAS.tsv.gz`，字段 `SNP A1 A2 EAF`，覆盖全部后验 SNP。EAF 为 A1 频率；代码支持无歧义翻转/互补链，拒绝歧义和缺失。

输出位于 `/mnt/d/data/ukb/pgs/<trait>/`：

- `csx.pgs.gz`：原有模型的平均分数。
- `csx.posterior.tsv.gz`：四个 discovery-centred 后验均值和十个协方差元素。
- 配套 `.json` 和 `.metadata.tsv`：来源、染色体、校验值；必须一起保留。

SNP 抽样在 `.../csx/<trait>/<signature>/raw/chr*/joint_posterior.h5`；如以后需要重新评分而不重跑 MCMC，应保留。完成后的个体 posterior 表可直接供 Yeval 复用。

PLINK 使用 `center no-mean-imputation`，只读 **SUM**。中心化后缺失位点贡献 0，等价于用相同 discovery EAF 均值填补；不能读取随缺失改变分母的 AVG。这也避免旧 PLINK 将未中心化均值加到缺失位点的问题。Yeval 使用 posterior 文件中的四个均值进行组合，使均值、方差和尺度对应同一个预测器。

## 2. 训练中心：必须有真实来源

默认 `--distance-source discovery`，要求 `<score-dir>/<trait>/training_centers.tsv`，或显式 `--training-centers FILE`。不能从最终 PRS 或 1KG 标签推断真实 GWAS 训练中心。

**有发现样本的同坐标投影：** 输入至少含 `POP PC1 ... PC10`，POP 为 AFR/EAS/EUR/SAS。

```bash
python f/training_geometry.py --mode individuals \
  --input /path/to/discovery_projected_pcs.tsv --pcs 10 \
  --pca-space ukb_1kg_projection_v1 \
  --source 'Actual discovery cohorts projected with the same SNP weights and scaling' \
  --output /mnt/d/data/ukb/pgs/height/training_centers.tsv
```

计算均值，不使用中位数。如果只提供发现人群的子样本，中心只能称为该子样本估计；输出 `N_GWAS` 应按实际 GWAS 人群权重核对，不能把代表性子样本人数当成完整 GWAS N。

**只有 discovery EAF：** 对现有 0pca 的线性剂量和投影，可以计算 `mean_PC = Σ 2 × EAF × PC_weight`，但必须覆盖目标投影的全部实际位点，并保持等位基因及尺度完全一致。

```bash
python f/training_geometry.py --mode eaf \
  --input /path/to/discovery_frequency_manifest.tsv \
  --pca-weights /mnt/d/files/DiscoDivas/g1k_hm3_maf5_woamb_wolr.pca.weight \
  --projection-snps /path/to/exact_target_projection_snps.txt \
  --frequency-allele A1 --pcs 10 --pca-space ukb_1kg_projection_v1 \
  --source 'Discovery EAF on the complete target-projection SNP set' \
  --output /mnt/d/data/ukb/pgs/height/training_centers.tsv
```

manifest 字段 `POP FILE N_GWAS`；FILE 指向 `SNP A1 EAF` 表。A1 必须已与 PCA 权重第 6 列完全一致，不猜测翻转。权重沿用 0pca 格式：SNP 第 2 列，效应等位基因第 6 列，PC 权重从第 7 列开始。SNP 列表取本次目标投影 `chr*.sscore.vars` 的实际位点，不是整个权重表。

EAF 中心对应无缺失基因型的期望投影；缺失严重或投影有额外变换时，应使用实际发现样本的投影均值。原始 GWAS EAF 比标准化后仅保留 HM3 的文件可能覆盖更多 PCA 位点。

也可直接填写 `templates/training_centers.tsv`。`pca_space` 是实际坐标系的身份声明，字符串相同本身不能证明坐标可比。不同 trait 的训练样本组成不同，应分别准备。

多训练群体距离为 `sqrt(Σ (N_g/ΣN) × d_ig²)`，并保存到每个源人群的距离。它避免简单混合中心因不同人群位置互相抵消而误导，但仍是描述性扩展，不保证准确度单调下降，也不是 Nature 的单训练群体距离。

## 3. 连续表型的遗传方差

复制并填写 `templates/genetic_variance.tsv`，不可直接运行含 NA 的模板：

| 字段 | 含义 |
|---|---|
| trait / target | 如 height / EUR；每个 trait-target 一行 |
| scale | ct 必须为 residual_phenotype |
| h2 | 与表型、协变量及 SNP 范围匹配的外部 SNP 遗传力，0<h2<1 |
| genetic_variance | 或直接给同一表型单位的遗传方差；与 h2 二选一 |
| source | 估计来源、样本和尺度说明 |

提供 h2 时，用 `h2 × var(training phenotype residuals)` 换算，每折只用训练样本。个体 model-based R² 为 `1 − posterior_variance/genetic_variance`。负值保留，提示尺度、先验或校准问题。它不等于 panel b 的实测 Prediction R²。

T2DM 默认画后验 SD；若探索 model-based R²，需绝对 `genetic_variance`，尺度为 `log_hazard`（t2e）或 `log_odds`（dt）。**不能将病例对照 liability h2 直接代入。**

## 4. Yeval

准备每个 trait 的训练中心，以及包含 height/LDL 遗传方差的表后：

```bash
./Yeval.sh --trait height --type ct --covar-name age,sex,PC1,PC2 \
  --pca-space ukb_1kg_projection_v1 \
  --genetic-variance-file /path/to/genetic_variance.tsv

./Yeval.sh --trait ldl --type ct --covar-name age,sex,PC1,PC2,drug.lipid \
  --pca-space ukb_1kg_projection_v1 \
  --genetic-variance-file /path/to/genetic_variance.tsv

./Yeval.sh --trait t2dm --type t2e \
  --covar-name age,sex,PC1,PC2,drug.dm,drug.htn \
  --pca-space ukb_1kg_projection_v1
```

LDL 不自动除以 0.7。T2DM 使用 `t2dm.Yt2e` / `t2dm.t2e`，可用 `--event-col` / `--time-col` 覆盖。已有 COJO 可显式用 `--pt-file`，避免自动评分。phi 沿用上游固定值/auto，没有额外 phi 网格调参。

暂时没有训练中心或遗传方差时，可明确选择探索模式：

```bash
./Yeval.sh --trait height --type ct \
  --distance-source reference --individual-metric sd
```

此时 c/d 标明四个 1KG 参考中心的代理距离，d 为真实后验 SD，不声称是准确度。没有 posterior 文件则不会画伪散点；仅做原有评分基准可显式 `--posterior-mode off`，关闭个体后验分析。

## 输出顺序

默认 `/mnt/d/analysis/grid/Yeval/<trait>/`：

1. `comparison.png`：含 COJO 的整体方法比较。
2. `combined_scores.png`：auto/fixed meta、四分数回归及 Disco 组合策略；与第一图共享的 bar 数值相同，两图有颜色图例。
3. `distance_performance.png`：PRS-CSx 四面板。
4. `paired_improvement.png`：DiscoDivas-tuned vs PRS-CSx，放在四面板之后。
5. `distance_bins.png`：单独的经验距离分箱图；t2e 显示 ΔC。

`performance.tsv` 保存主指标、增益、旧 SSE 指标及 RMSE 等；`individual_posterior.tsv.gz` 保存全部个体后验 SD、model-based R²（有尺度时）、对角/交叉协方差贡献、组合权重及源人群距离。`distance_centers.tsv` 和 `fold_coefficients.tsv` 保留可复核定义。

默认每个祖源最多显示 5000 个真实个体点，全部估计保存在表中。抽样只用于显示，不根据结局或准确度选点。经验分箱使用全部合格样本，并根据事件数减少稀疏箱。

`--write-predictions TRUE` 另存 OOF 预测。`--bootstrap 0` 仅用于快速检查。`--out-root` 可以避免不同结局类型覆盖同一个 trait 报告。

## 其他输入与验证

默认 GWAS `/mnt/e/gwas/4grid/common`；LD `/mnt/e/refLD/csx`；基因型沿用现有 hap/typ 路径；表型 `/mnt/d/data/ukb/phe/Rdata/all.rds`；PRS `/mnt/d/data/ukb/pgs/<trait>/`。T2DM AFA 仍映射 AFR。详细参数见各主入口 `--help`。

`grid.sh`、`f/arg*` 等实验模块保留当前 GitHub 实现；本次没有重新对其性能作结论。

验证包括：实际 1csx.sh 小数据端到端运行及缓存续跑；实际小规模 PRS-CSx MCMC；真实 PLINK 评分中的缺失基因型、等位基因翻转和字符串 ID；两个染色体独立后验合并；连续、二分类、生存结局的模拟端到端运行；独立重算 Prediction R²、AUC、C、训练折组合、协方差与遗传方差尺度。尚未在完整 UKB 数据上运行，真实数据上的 MCMC 收敛及个体可靠性校准仍需评估。

参考：PRS-CS https://doi.org/10.1038/s41467-019-09718-5；PRS-CSx https://doi.org/10.1038/s41588-022-01054-7；Ding et al. https://doi.org/10.1038/s41586-023-06079-4。
