# GRID — audit v3

基于 `jielab/scripts@8c49af6` 的二次审计修订。详细问题、方法边界与验证记录见 **AUDIT_V3.md**。

## 安装

将本包的整个 `grid` 目录放到 `/mnt/d/scripts/grid`，先备份原脚本目录，避免把不同版本的 `f/` 混在一起。数据、GWAS、PCA 和分析结果目录不包含在本包中。

```bash
cd /mnt/d/scripts/grid
chmod +x *.sh
# 已有 grid 环境时更新依赖：
conda env update -n grid -f environment.yml
# 首次完整安装（也安装 ARG-Needle 依赖）：
# bash install_grid.sh
```

入口默认激活 `$HOME/miniforge3/envs/grid`。其他位置设置 `GRID_CONDA_ENV`。
R 优先使用当前环境；必要时设置 `GRID_RSCRIPT=/path/to/Rscript`，并确保该 R 有报告所需的包。

## 三个 trait 的 CSx / Disco 比较

```bash
# 队列级 PCA，只需准备一次；检查通过不等于已计算。
./0pca.sh --check
./0pca.sh

./1csx.sh --traits height,ldl,t2dm --check
./1csx.sh --traits height,ldl,t2dm --jobs 4 --threads 1

# 保存无表型调优的官方 Disco 插值，作为诊断比较。
./2disco.sh --traits height,ldl,t2dm

# 默认额外在每个训练折内调优 Disco 输入模型。
./Yeval.sh --trait height --type ct
./Yeval.sh --trait ldl --type ct
./Yeval.sh --trait t2dm --type t2e
```

上述 LDL 命令沿用输入表型。只有在你的分析方案要求、且表型上游尚未重复处理时，才显式加 `--covar-name age,sex,PC1,PC2,drug.lipid`。本程序不自动除以 0.7。

二分类 T2DM：

```bash
# 使用已有、明确的 0/1 列：
./Yeval.sh --trait t2dm --type dt --phenotype-col t2dm
# 不指定 phenotype-col 时，按 Yr2e/Yt2e 明确构造基线 ICD10 T2DM。
```

**dt 的 R² 现在是 observed-scale partial R²，另有 AUC/Brier；不是 liability R²。**
`ct` 报告 OOF partial R²、总 R² 和 RMSE；`t2e` 报告考虑删失且在折内比较的 Harrell C。
`--prevalence` 仅保留为患病率背景描述，不参与一个未经验证的 partial-R² liability 转换。

## 输入和输出

| 内容 | 默认位置 |
|---|---|
| GWAS | `/mnt/e/gwas/4grid/common/<trait>.<POP>/gwas/<trait>.<POP>.gz`；也支持目录下 `<trait>.<POP>.gz` |
| 非欧洲 T2DM 别名 | `t2dm.AFA` 映射 AFR LD；输出保留输入名称 |
| CSx target | 按原环境查找 `/mnt/d/data/ukb/gen/typ` 或 `/mnt/e/ukbGen/37/hap`；建议显式 `--dir-gen` |
| CSx SNP 列表 | `<dir-gen>/ukb_array.bim`；可显式 `--csx-bim-prefix` |
| LD/SNPINFO | `/mnt/e/refLD/csx`；`--csx-ref-dir` / `--csx-snpinfo` |
| PCA imputed target | `/mnt/e/ukbGen/37/imp/chr1..22`；`--dir-imp` |
| PCA 权重/中心 | `/mnt/d/files/DiscoDivas/` |
| 投影与祖源结果 | `/mnt/d/data/ukb/pca_proj/` |
| 全队列表型 | Yeval 默认 `/mnt/d/data/ukb/phe/Rdata/all.rds`；`--pheno-file` |
| CSx 永久 SNP 权重 | 原 GWAS 名加 `.csx.gz`，另有签名和元数据 |
| CSx 永久 GWAS 预处理 | 原 GWAS 名加 `.csx.sumstats.gz`，配套 `.csx.sumstats.json` |
| 个体分数 | `/mnt/d/data/ukb/pgs/<trait>/csx.pgs.gz`；6 个 CSx 模型列 |
| untuned Disco | 同目录 `disco.pgs.gz` 和 `disco.coef.tsv.gz` |
| Yeval | `/mnt/d/analysis/grid/Yeval/<trait>/report.html` 等，成功后有 `SUCCESS` |

切换 `ct/dt/t2e` 默认会覆盖同一 trait 的报告；需要同时保留时设置不同 `--out-root`。
`--write-predictions TRUE` 可增加一个含 IID、折分及各模型预测的压缩表。默认保持报告目录简洁。
`--bootstrap 0` 仅用于快速流程检查，没有置信区间。

### Yeval 图表和指标

报告依次展示以下四张图；前两图均有方法颜色图例，重复方法使用相同估计值。

1. **Prediction by target ancestry**：包括 COJO 的整体方法比较。
2. **Combined-score comparison**：聚焦 CSx 分数如何组合，比较 auto-meta、fixed-meta、四分数回归和 Disco 插值；去掉 COJO，增加 fixed-meta。
3. **Prediction along genetic distance**：只展示 PRS-CSx，采用 2×2 四面板。a 为原始祖源标签着色的 PCA；b 为各祖源表现；c 为相同样本按到 1KG EUR 中心的距离着色；d 为各祖源内部距离分箱的预测表现。
4. **DiscoDivas-tuned versus PRS-CSx**：补充差值图。

这里的“原始祖源”指 `--group-col` 指定的已有标签（默认 `genetic_ancestry`），不是根据遗传距离重新分组，也不自动等同于自报族群。OTH/UNASSIGNED 保留原标签。
四面板结构参考 [Ding et al., Nature 2023, Fig. 1](https://www.nature.com/articles/s41586-023-06079-4/figures/1)。原论文图是示意图；本报告用实际群体/分箱估计，不生成个体 R²，也不强加随距离下降的曲线。
PRS-CSx 仍沿用每个目标祖源训练折拟合的四分数组合，因此不同祖源并非共享同一套组合系数。

**height/LDL 的 partial R²**：`1 − SSE(full)/SSE(covariates)`，其中两个 SSE 都来自留出预测。
也等于 `(full_R2 − baseline_R2)/(1 − baseline_R2)`。它衡量加入 PRS 减少了多少协变量模型剩余误差；`delta_R2` 则以表型总方差为分母。
报告开头用本次 PRS-CSx 结果代入公式，`performance.tsv` 保存两个 SSE、两个总 R²、增量 R² 及 RMSE。partial R² 不是按参数数目校正的 adjusted R²。
PRS-CS 论文描述的是协变量调整后的 observed-versus-predicted R²；仅凭该描述不能断言与本程序的 OOF SSE 定义数值完全一致，尤其预测校准不同时。

**T2DM 生存分析**：主图显示 `covariates + PRS` 的 Harrell C-index；前两图虚线和四面板 b 的黑色刻线显示 covariates-only C。
概览表同时显示 `baseline_C`、`delta_C` 及配对区间，以及 N、事件数和观察到的事件比例。事件比例不是固定时间风险；C-index 也不是 R² 或 PRS 单独的贡献。
二分类 `--type dt` 继续使用 observed-scale partial R²，并单列 AUC、baseline_AUC、delta_AUC 和 Brier。

距离使用投影的前 10 个 PC 和 1KG EUR 参考中心，不是发现 GWAS 的中心。图中平方根刻度只改变显示间距，刻度值保持原始距离单位。
`--distance-bins` 默认最多 10 个祖源内等人数箱，`--min-n` 默认每箱 100 人；`--min-bin-events` 默认要求 dt/t2e 每箱至少 20 个事件/病例和 20 个非事件/对照，必要时减少箱数。无法形成两个合格箱的祖源只显示总体结果。
所有分箱使用已有留出预测，不在箱内重新拟合。`distance_performance.tsv` 保存实际箱数、边界、中位距离、N、事件数和指标，便于复查图中点。

替换代码后需在 UKB 数据环境重跑 Yeval；旧的汇总表和图片不足以恢复个体距离分箱或协变量基线 C。重跑可显式使用已有 `--pt-file`，无需重新运行 CSx/Disco 或 COJO 评分。

`--chrs 22` 等子集有独立的权重文件名与评分目录 `<trait>/chr22/`，不会覆盖全基因组结果。Yeval 对子集运行需显式指定对应的 `--pgs-file` 等输入。

## 各方法准确含义

| 标签 | 内容 |
|---|---|
| COJO | 同祖源 `.jma.cojo` 的 SNP/refA/bJ；不是标准 P+T 网格搜索 |
| PRS-CSx-auto-meta | 自动学习 φ，再合并各祖源后验效应 |
| PRS-CSx-fixed-meta | 使用固定 φ，再合并各祖源后验效应 |
| PRS-CSx | 目标祖源训练折内拟合四个 CSx 分数的组合 |
| DiscoDivas-tuned | 四个祖源开发模型各自组合所有 CSx 分数，再基于训练中心插值；所有表型调优仅用外层训练折 |
| DiscoDivas-untuned | `2disco.sh` 保存的原始四分数官方插值，作为补充诊断 |

`--disco-tune FALSE` 只评价保存的 untuned Disco，不计算 tuned 分支。
默认每个祖源每折至少 100 个开发样本；不足时明确报错，可用上述选项做诊断。
Yeval 的 `--disco-a` 与 `2disco.sh --a-list` 顺序统一为 **AFR,EAS,EUR,SAS**。默认全为 1。

没有 φ 网格选择，也没有独立实现 PRS-CS-mult；标签不会混用。
官方权重几何不是局部祖源推断，插值系数也不是 ancestry proportion。

## 续跑与缓存

- 固定 φ 的 populations 和 fixed-meta 共享同一条联合 MCMC 链，auto 使用学习 φ 的另一条链。
- 标准化 GWAS 永久保存在原始 GWAS 旁，例如 `ldl.AFR.csx.sumstats.gz`，包含原始效应 BETA、SE、N 等；`ldl.AFR.csx.gz` 则是 MCMC 后的 SNP 权重。
- populations、auto/meta 和 GRID 共用上述预处理缓存。原始 GWAS、SNPINFO 和预处理代码签名匹配时直接复用；改变 φ、MCMC 参数、目标基因型或删除分析临时目录不要求重新预处理。`--replace TRUE` 可强制重新预处理。
- 发布前校验压缩文件完整性和行数，并加锁、原子发布。分析目录按需生成独立工作副本及染色体拆分文件；每次推断使用的 φ/N 记录在各自运行目录，不会写入共享缓存。
- 对仍保存在旧分析目录的匹配预处理表，运行时可自动导入永久缓存。迁移工具 `f/sumstats_cache.py --input GWAS --snpinfo SNPINFO --trait TRAIT --pop POP --work /mnt/d/analysis/grid/csx --migrate-only --remove-legacy` 会在校验成功后移除旧压缩表，保留运行元数据；使用前应停止对应性状的任务。
- 重新发布 populations 不会删除已有 auto/meta 列。若样本集合改变，使用独立 `--score-dir`，避免混合不同队列。
- 缓存包含配置、输入文件路径/大小/修改时间及相关代码；不对大型基因型文件逐字节计算 hash。
- 旧版本无签名或签名不一致时重算。PCA 旧缓存也可能触发一次完整重建，随后按投影/距离/祖源三个阶段复用。
- `--replace TRUE` 强制重算所选模块。改变 target、keep/remove、φ 或参考 LD 后无需靠手工删除缓存来更新。
- CSx 的 SNPINFO/BIM/LD 必须来自正确的 genome build 和祖源，参考 SNP 匹配/评分数应从日志及 `.sscore.vars` 核查。不要把软件成功退出等同于 SNP 覆盖充分。
- GWAS 预处理与预检查共用 98% 坐标匹配下限；预处理在完整输入上再次检查。少量冲突位点会被剔除，并在预处理日志和元数据的 `coordinate_mismatch` / `coordinate_mismatch_examples` 中记录；低于下限则停止，不发布部分结果，也不将冲突坐标改写为参考坐标。
- 顺序运行 CSx 和 Disco 时可用 `./1csx.sh ... && ./2disco.sh ...`；这样 CSx 失败后不会继续用评分目录中已有的旧分数运行 Disco。

## GRID 实验流程

ARG 构建仍由已有 `refGen.sh` / `gu/arg.sh` 数据准备流程负责。GRID 消费已有 `.argn`、`.trees`、sample map、anchors 和 `.variants.tsv.gz`。

```bash
./1csx.sh --trait height --chrs 1-22 --jobs 4 --threads 1
./grid.sh grid --trait height --chrs 1-22 \
  --arg-dir /mnt/e/ukbGen/37/arg --jobs 4 --threads 1
./Yeval.sh --trait height --type ct \
  --grid-file /mnt/d/analysis/grid/height/scores/grid.tsv.gz
```

`grid.sh grid` 先检查 ARG，然后衔接新版永久 CSx 输入，依次计算 LD、运输模型、收缩权重和分数。
GRID 中间输入目录 `/mnt/d/analysis/grid/<trait>/csx/` 由桥接步骤按当前输入生成，不需要用户手工复制。

`--grid-transport-model evolutionary_full` 是预先指定的默认模型；`baseline` 是消融模型。
在查看 OOF 结果后切换到更好的模型，会引入选择偏差；需要额外独立验证。
单染色体仅适合作为区块验证的 pilot；数据必须覆盖足够多区块。

## 验证说明

当前脚本目录不包含开发测试和测试产物；此前已进行合成数据及小型真实程序测试。
本次没有重跑真实 UKB 三个 trait，也没有证明 GRID 或 tuned Disco 的真实预测性能已经提高。

统一使用 `0pca.sh / 1csx.sh / 2disco.sh`；详细参数用各主入口 `--help`。根目录不再提供同功能的无编号入口。`f/` 下的同名文件是内部执行模块，必须保留。
