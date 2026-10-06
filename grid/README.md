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

将完整 `grid` 替换到 `/mnt/d/scripts/grid`，避免混用内部模块。公共入口为 `0pca.sh`、`1csx.sh`、`2disco.sh`、`3grid.sh` 和 `Yeval.sh`。

```bash
cd /mnt/d/scripts/grid
chmod +x *.sh
conda env update -n grid -f environment.yml
```

- PCA 算法未因目录整理改变。 训练中心必须在相同坐标系；如果改用训练样本自己的 PCA 基底，则目标必须重新投影。
- 如果只有旧平均权重，生成逐人后验仍需保存抽样的推断结果。已有后验表可继续供 Yeval 使用。
- 可以复用已保存的 Disco/COJO 分数作为比较输入；Yeval 必须重跑。
- 本包没有 UKB 数据，也没有声称已重算的真实结果。


## 代码目录（2026-10-03）

| 模块 | 入口 | 辅助文件 |
|---|---|---|
| 公共环境、参数、数据准备、评分文件 | `f/0.common.sh preflight` | `f/0.common.py` |
| PCA、距离、祖源 | `0pca.sh` | `f/0.pca.R`、`f/0.pca.py` |
| ARG 数据准备工具 | 由 `../gu/arg.sh` 调用 | `f/0.arg.py` |
| PRS-CSx | `1csx.sh` | `f/1csx.py` |
| DiscoDivas | `2disco.sh` | `f/2disco.R` |
| GRID | `3grid.sh grid` | `f/3grid.f.sh`、`f/3grid.py` |
| 评估 | `Yeval.sh` | `f/Yeval.R`、`f/Yeval.py` |

Python 辅助文件通过子命令选择步骤，`--help` 列出子命令。CSX 的 populations、auto/meta 和评分 shell 逻辑全部位于根目录 `1csx.sh`。共享内存限制实现位于 `../0f/memory_cap.sh`，保持原有 cgroup 校验。

`f/csx/` 保存 PRS-CSx 第三方实现；DiscoDivas 许可证和来源记录仍在 `f/disco/`。这些第三方文件保留其代码结构。

结果路径及数据表结构不变。代码签名仍参与缓存校验，重构前的缓存可能失效；不会仅通过改写签名把未经验证的旧结果认定为新结果。临时目录小数据验证覆盖三个 CSX 模式、后验评分、缓存续跑和 GRID 拟合，并对照整理前代码逐表比较。

## 1. 生成个体后验

`1csx.sh` 默认对整次调用及全部并行子进程设置 **16 GiB RAM + 2 GiB swap** 的 cgroup v2 上限，超限无法回收时仅终止该任务组。需要 systemd user manager；限制无法建立或核实时拒绝无上限运行。多个独立调用各自计额；这不是整台 WSL 的总上限，也不限制 GPU 显存。内部 `f/0.common.sh` 和 Python 模块不应作为有内存保护的入口直接调用。

```bash
./1csx.sh --traits ldl,t2dm --jobs 4 --threads 2 --memory-cap-gb 16 --swap-cap-gb 2 --score-memory 2048
```

可通过 `GRID_MEMORY_CAP_GB` / `GRID_SWAP_CAP_GB` 设置默认上限（整数 GiB，swap 可为 0）。普通 population/auto/meta PLINK 评分统一用 `--score-memory`，默认每个进程 **2048 MiB workspace**；`--posterior-memory` 仍默认 8192 MiB，按染色体顺序运行。最多用 RAM 上限的 75% 分配评分 workspace，超过时自动降低**评分**并发数，保留 `--jobs` 控制 MCMC 并发。PLINK workspace 不等于进程总内存，最终总上限由 cgroup 执行；可使用更小的 `--jobs` 降低 MCMC 峰值。

启动时打印 scope 名称和 PID，可用 `systemctl --user status <scope>` / `journalctl --user -u <scope>` 查询。2026-09-24 的错误日志已带 `--memory 4096`，当时整机仅余约 465 MiB；单靠 PLINK 参数不能保护 WSL。

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

- `1csx.scores.rds`：原有模型的平均分数。
- `1csx.posterior.rds`：四个 discovery-centred 后验均值和十个协方差元素。
- posterior 的必要输入校验、染色体和中心化信息保存在 RDS 的 `model_info` 属性内，不再生成 JSON 或 TSV 边车文件。

SNP 抽样在 `.../csx/<trait>/phi-1e-2/raw/chr*/joint_posterior.h5`（自动 phi 为 `auto/`）；如以后需要重新评分而不重跑 MCMC，应保留。完成后的个体 posterior 表可直接供 Yeval 复用。

CSx 的 MCMC 推断结果仍按性状、phi 和染色体保存，以便重新评分时复用。原始 HDF5 后验抽样、SNP 权重及工具恢复运行必需的状态属于模型数据，保留其原生格式。普通评分、combined 评分、posterior 评分及 PCA、DISCO 工作目录统一位于 `/tmp/grid-cache/`；日志、命令和交换表也在 `/tmp`。代码根据正式路径计算稳定的缓存位置。仅修改注释或排版不触发重算；确需重新计算时使用 `--replace TRUE`。

PLINK 使用 `center no-mean-imputation`，只读 **SUM**。中心化后缺失位点贡献 0，等价于用相同 discovery EAF 均值填补；不能读取随缺失改变分母的 AVG。这也避免旧 PLINK 将未中心化均值加到缺失位点的问题。Yeval 使用 posterior 文件中的四个均值进行组合，使均值、方差和尺度对应同一个预测器。

## 2. 训练中心：必须有真实来源

默认 `--distance-source reference`，使用探索性的 1KG 参考中心距离。显式选择 `--distance-source discovery` 时，要求 `<score-dir>/<trait>/training_centers.tsv`，或通过 `--training-centers FILE` 指定，并提供匹配的 `--pca-space`。不能从最终 PRS 或 1KG 标签推断真实 GWAS 训练中心。

**有发现样本的同坐标投影：** 输入至少含 `POP PC1 ... PC10`，POP 为 AFR/EAS/EUR/SAS。

```bash
python f/0.pca.py training-centers --mode individuals \
  --input /path/to/discovery_projected_pcs.tsv --pcs 10 \
  --pca-space ukb_1kg_projection_v1 \
  --source 'Actual discovery cohorts projected with the same SNP weights and scaling' \
  --output /mnt/d/data/ukb/pgs/height/training_centers.tsv
```

计算均值，不使用中位数。如果只提供发现人群的子样本，中心只能称为该子样本估计；输出 `N_GWAS` 应按实际 GWAS 人群权重核对，不能把代表性子样本人数当成完整 GWAS N。

**只有 discovery EAF：** 对现有 0pca 的线性剂量和投影，可以计算 `mean_PC = Σ 2 × EAF × PC_weight`，但必须覆盖目标投影的全部实际位点，并保持等位基因及尺度完全一致。

```bash
python f/0.pca.py training-centers --mode eaf \
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

所有表型默认 `--individual-metric sd`，展示后验 SD，无需遗传方差。连续表型如需展示个体 model-based R²，显式选择 `--individual-metric reliability` 并提供下面的遗传方差表。

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

`./Yeval.sh -h` 的三条示例可直接复制运行，使用上述默认输入路径和已完成的后验表：

```bash
./Yeval.sh --trait height --type ct --covar-name age,sex,PC1,PC2
./Yeval.sh --trait ldl --type ct --covar-name age,sex,PC1,PC2,drug.lipid
./Yeval.sh --trait t2dm --type t2e --covar-name age,sex,PC1,PC2,drug.dm,drug.htn
```

默认使用探索性的 1KG 参考中心距离，个体面板展示后验 SD，因此不需要训练中心、`--pca-space` 或遗传方差文件；height/LDL 的群体预测指标仍为 Prediction R²。可在任一命令末尾追加 `--check`，只检查输入、不拟合或生成 COJO 分数。T2DM 也支持 `--type dt`，未指定 `--phenotype-col` 时从 `t2dm.Yr2e` / `t2dm.Yt2e` 推导基线 T2DM；更改类型会覆盖同一输出目录中的先前报告。

需要使用真实发现样本距离及连续性状的个体 model-based R² 时，使用下面的完整配置。

准备每个 trait 的训练中心，以及包含 height/LDL 遗传方差的表后：

```bash
./Yeval.sh --trait height --type ct --covar-name age,sex,PC1,PC2 \
  --distance-source discovery --individual-metric reliability \
  --pca-space ukb_1kg_projection_v1 \
  --genetic-variance-file /path/to/genetic_variance.tsv

./Yeval.sh --trait ldl --type ct --covar-name age,sex,PC1,PC2,drug.lipid \
  --distance-source discovery --individual-metric reliability \
  --pca-space ukb_1kg_projection_v1 \
  --genetic-variance-file /path/to/genetic_variance.tsv

./Yeval.sh --trait t2dm --type t2e \
  --covar-name age,sex,PC1,PC2,drug.dm,drug.htn \
  --distance-source discovery --pca-space ukb_1kg_projection_v1
```

LDL 不自动除以 0.7。T2DM 使用 `t2dm.Yt2e` / `t2dm.t2e`，可用 `--event-col` / `--time-col` 覆盖。已有 COJO 可显式用 `--pt-file`，避免自动评分。phi 沿用上游固定值/auto，没有额外 phi 网格调参。

默认探索模式下，c/d 标明四个 1KG 参考中心的代理距离，d 为真实后验 SD，不声称是准确度。没有 posterior 文件则不会画伪散点；仅做原有评分基准可显式 `--posterior-mode off`，关闭个体后验分析。

## 输出顺序

默认 `/mnt/d/analysis/grid/Yeval/<trait>/`：

1. `Yeval.comparison.png`：含 COJO 的整体方法比较。
2. `Yeval.combined_scores.png`：auto/fixed meta、四分数回归及 Disco 组合策略；与第一图共享的 bar 数值相同，两图有颜色图例。
3. `Yeval.distance_performance.png`：PRS-CSx 四面板。
4. `Yeval.paired_improvement.png`：DiscoDivas-tuned vs PRS-CSx，放在四面板之后。
5. `Yeval.distance_bins.png`：单独的经验距离分箱图；t2e 显示 ΔC。

`Yeval.comparison.xlsx` 保存主指标、增益及 RMSE 等；`Yeval.individuals.xlsx` 保存全部个体后验 SD、model-based R²（有尺度时）、对角/交叉协方差贡献、组合权重及源人群距离。参考中心保存在 `Yeval.distance_performance.xlsx`，逐折系数与入组计数保存在 `Yeval.models.xlsx`。

默认每个祖源最多显示 5000 个真实个体点，全部估计保存在表中。抽样只用于显示，不根据结局或准确度选点。经验分箱使用全部合格样本，并根据事件数减少稀疏箱。

`--write-predictions TRUE` 另存 OOF 预测。`--bootstrap 0` 仅用于快速检查。`--out-root` 可以避免不同结局类型覆盖同一个 trait 报告。

## 其他输入与验证

默认 GWAS `/mnt/f/gwas/4grid/common`；LD `/mnt/f/refLD/csx`；基因型沿用现有 hap/typ 路径；表型 `/mnt/d/data/ukb/phe/Rdata/all.rds`；PRS `/mnt/d/data/ukb/pgs/<trait>/`。T2DM AFA 仍映射 AFR。详细参数见各主入口 `--help`。

`3grid.sh` 和 `f/3grid.*` 的结果读写已接入 RDS／XLSX；分析方法保持原样。本次用合成结果验证发布入口，没有重新评估其预测性能。

验证包括：实际 1csx.sh 小数据端到端运行及缓存续跑；实际小规模 PRS-CSx MCMC；真实 PLINK 评分中的缺失基因型、等位基因翻转和字符串 ID；两个染色体独立后验合并；连续、二分类、生存结局的模拟端到端运行；独立重算 Prediction R²、AUC、C、训练折组合、协方差与遗传方差尺度。本次只检查现有 height、ldl、t2dm 数据的输入及迁移完整性，没有重跑这些模型；上述既有方法测试不替代真实数据上的 MCMC 收敛与个体可靠性校准评估。

参考：PRS-CS https://doi.org/10.1038/s41467-019-09718-5；PRS-CSx https://doi.org/10.1038/s41588-022-01054-7；Ding et al. https://doi.org/10.1038/s41586-023-06079-4。

## 结果文件与临时目录

- `analysis/grid/Yeval/<trait>/`：每张 `Yeval.*.png` 对应同名 XLSX；`Yeval.models.xlsx` 保存逐折系数和入组计数，`report.html` 包含方法说明。
- 个体评估结果使用 `Yeval.individuals.xlsx`；可选的所有方法预测使用 `Yeval.predictions.xlsx`。有 ID 的结果表照常进入工作簿。
- `data/ukb/pgs/<trait>/`：`1csx.scores.rds`、`1csx.posterior.rds`、`2disco.scores.rds` 和 `Yeval.cojo.rds` 保存后续评估会读取的大型评分数据；`2disco.coefficients.xlsx` 保存供查看的个体系数表；`Yeval.cojo.xlsx` 保存变异匹配质量结果。
- 新的 GRID 方法输出位于 `analysis/grid/3grid/<trait>/`：`3grid.model.xlsx` 为验证和系数表，`3grid.model.rds` 为模型，`3grid.conservation.xlsx` 和 `3grid.variant_predictions.xlsx` 为完整变异结果，`3grid.scores.rds` 为个体评分。
- PLINK 评分、Disco 输入、PCA 交换文件、拆分的 GWAS 输入、日志和运行锁位于 `/tmp/grid-cache/` 或 `/tmp/grid/Yeval/`；正式目录不保留这些临时副本。已有 PRS-CSx MCMC 后验和 SNP 权重仍保留供复用，不重新运行推断。
- Python/R 结果读写共用 `../0f/results.py`、`../0f/results.R`。结果工作簿保存精确表格导出；是否含 ID 不决定文件格式。
- 旧评估未保存 `paired_improvement` 的 bootstrap 区间。整理时保留原 PNG，工作簿明确标记该限制；新运行会保存完整的配对比较结果，不从图片估算区间。
