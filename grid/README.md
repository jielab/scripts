# GRID

Genealogy-informed effect transport for multi-ancestry PRS.

## 文件结构

```text
pca.sh          队列级准备：投影到 1KG PC 空间、参考距离、祖源推断
csx.sh          独立运行 PRS-CSx，生成权重和个体分数
disco.sh        独立运行现有 DiscoDivas zero-shot 比较方法
grid.sh         GRID 计算与 eval；保留 pca 快捷入口
README.md       使用说明
install_grid.sh 环境安装脚本
environment.yml Conda 环境定义
f/              实现代码、PRS-CSx 源码、安装脚本和测试
```

参考文件独立放在 `/mnt/d/files/DiscoDivas/`（Windows：`D:\files\DiscoDivas`）：

- `g1k_hm3_maf5_woamb_wolr.pca.weight`
- `med.g1000.4pop.tsv`

## 安装与使用

首次使用运行 `install_grid.sh` 安装 `grid` 环境。四个入口会自动激活
`$HOME/miniforge3/envs/grid`，无需手动执行 `conda activate`；其他位置通过
`GRID_CONDA_ENV` 指定。

```bash
cd /mnt/d/scripts/grid
bash install_grid.sh
./grid.sh -h
```

## PCA：同一队列只准备一次

```bash
./grid.sh pca # project PCA to 1KG, calculate genetic distance, infer ancestry
# 等价的独立入口：
./pca.sh
```

PCA 是使用已有 1KG 权重做投影（PLINK2 `--score`），不重新拟合 PCA。
UKB 自带的 PCs 不能直接替代此处与 1KG 中心一致的坐标。
默认 PC1–PC20 用于协变量和祖源分类，PC1–PC10 用于到四个人群中心的欧氏距离。

三个内部阶段分别续跑：

1. 已有完整投影表时跳过 PLINK 投影；否则复用已完成的染色体分数。
2. 投影不变且距离/QC 输出完整时跳过；缺失时直接从投影表计算。
3. 距离更新或祖源/QC 输出缺失时推断祖源；完整结果直接复用。

不存在独立 `ancestry` 模块。投影、距离和祖源结果供不同性状共用。
输入或参数改变后，用 `--replace TRUE` 重建；程序不自动识别全部参数变更。

默认文件：

```text
/mnt/d/data/ukb/phe/pca/ukb.discodivas.pca.tsv.gz
/mnt/d/data/ukb/phe/pca/ukb_reference_distances.tsv.gz
/mnt/d/data/ukb/phe/pca/ukb.ancestry.auto.tsv.gz
/mnt/d/analysis/grid/pca/chr*.sscore
```

可用 `--pca-file`、`--ancestry-file`、`--med-file`、`--pca-weight`、
`--dir-imp`、`--phe-file` 覆盖路径。PCA 始终投影 22 条常染色体。

祖源推断用自报族群与最近 1KG 中心一致的 UKB 样本作为训练锚点，训练均衡 LDA，
输出 AFR/EAS/EUR/SAS 概率；最高概率低于 0.90 时分为 OTH。
锚点不足或模型不可用时退回距离方法。这是当前 UKB 近似实现，
并非原始 PRS-CSx 论文中以带标签 1KG 个体训练的分类器；也不推断染色体局部祖源。

## 独立运行 CSX 和 Disco

```bash
./csx.sh --trait height --chrs 22
./disco.sh --trait height
```

PRS-CSx 源码位于 `f/csx/`，输出在 `分析目录/性状/csx` 和 `scores/csx.tsv.gz`。
Disco 沿用现有的 `f/disco_zero.R` 比较实现，使用投影 PCs 与参考中心计算连续人群权重，
输出 `scores/disco_zero.tsv.gz`。这次拆分不改变其方法实现。
两个入口不再属于 `grid.sh` 的模块。

## GRID 模块

用户只需要两个命令：`grid` 完成整套 GRID 计算，`eval` 单独进行性能评估。
检查 ARG、提取 LD、组装特征、拟合模型、生成权重和评分是 `grid` 内部连续执行的
步骤，不是需要用户分别调用的模块。

以染色体 22 为例：

```bash
./grid.sh grid --trait height --chrs 22
./grid.sh eval --trait height
```

GRID 输出在 `性状/grid/` 和 `性状/scores/grid.tsv.gz`，评估输出在 `性状/eval/`。
旧的跨方法 `all` 入口已移除。运行多个性状可写显式循环：

```bash
for trait in height ldl t2dm; do
  ./csx.sh --trait "$trait" --chrs 1-22
  ./disco.sh --trait "$trait"
  ./grid.sh grid --trait "$trait" --chrs 1-22
  ./grid.sh eval --trait "$trait"
done
```

## 数据及 ARG

默认 GWAS：`/mnt/e/gwas/4grid/{trait}.{AFR,EAS,EUR,SAS}.gz`。
默认分析目录：`/mnt/d/analysis/grid`；可用 `--output-root` 覆盖。
其他选项见 `./grid.sh -h`。`--replace FALSE` 复用缓存；更换输入/配置时明确重建。

ARG 建树属于共享数据准备，在 `refGen.sh` 中完成；GRID 只消费已有结果：

```bash
bash /mnt/d/scripts/gu/arg.sh build --method needle \
  --dir-gen /mnt/e/ukbGen/37 --dir-pfile /mnt/e/ukbGen/37/hap \
  --map-dir /mnt/e/refGen/maps/GRCh37 --chr 22
```

默认 ARG 位于 `/mnt/e/ukbGen/37/arg/{argn,trees}`，可用 `--arg-dir` 覆盖。
全队列建树的资源需求应先通过小规模实验评估；模型验证应使用多个染色体。
用于性能评估的 GWAS 应排除 UKB，以避免样本重叠导致的乐观估计。

## UKB 多方法预测评估（Yeval）

```bash
./Yeval.sh --trait height --type ct
./Yeval.sh --trait ldl --type ct --covar-name age,sex,PC1,PC2,drug.lipid
./Yeval.sh --trait t2dm --type dt --covar-name age,sex,PC1,PC2 --prevalence cohort
./Yeval.sh --trait t2dm --type t2e --covar-name age,sex,PC1,PC2
```

输出统一为 `/mnt/d/analysis/grid/Yeval/<trait>/`，不再创建 `csx/ct` 等子目录。
默认使用完整队列 `phe/Rdata/all.rds` 与 `ukb.ancestry.auto.tsv.gz`，避免 White-only
`ukb.phe` 无法评估其他祖源的问题。四种主方法为同祖源 COJO 的 PT、`csx.auto`
（按约定显示为 PRS-CS-multi）、训练折内组合四个祖源评分的 PRS-CSX、已保存的
DiscoDivas。`csx.meta` 和各单独祖源评分也在同一批样本及相同折上评估。

打开 `report.html` 查看四方法主图、auto/meta 比较、DiscoDivas 相对 PRS-CSX 的配对
差异、PCA/遗传距离分布和沿距离变化的性能。目录只保留网页、5 张配套 PNG、网页
链接的 `plots.pdf`（全部图）、`performance.tsv`、`cohort.tsv`、`methods.md`，以及
记录命令与运行过程的 `eval.log`、防止并发运行的 `run.lock`。
不再导出单图 PDF、逐人预测/折分/距离、折内系数及其他中间 TSV；输入配置、评分
定义、跳过的项目和患病率假设直接放在网页中。成功运行后会清理旧版本留下的这些
已知冗余文件，其他文件不受影响。`--check` 只把校验结果写入日志，不改已有报告。

PT 默认使用完整 imputed 基因型 `/mnt/e/ukbGen/37/imp/chr*`，按 COJO 位点提取，
使用 `SNP/refA/bJ` 和 `SCORE_SUM`；缓存为 `/mnt/d/data/ukb/pgs/<trait>/pt.pgs.gz`。
按 rsID、`Chr/bp` 和 `refA` 匹配；重复 rsID 只有唯一兼容记录时才参与评分。
`refA` 支持 SNP 和由 A/C/G/T 组成的完整插入/缺失序列（如 `CAA`、`TC`）；序列
必须与基因型等位基因完全匹配，不截取首个碱基，也不把符号编码猜测成序列。
临时 PVAR 使用逐记录唯一 ID，保留 PGEN 的记录顺序。无法区分的重复记录明确排除。
评分缓存旁的 `pt.pgs.gz.variants.tsv` 记录每个祖源/染色体的请求、评分、缺失、
坐标/等位基因不匹配、歧义排除及重复 ID 解析数量；`pt.pgs.gz.matches.tsv` 提供
每个请求位点的候选与匹配明细。这些缓存文件不再复制到评估目录。

`dt` 输出 Liability R²（Lee 转换）及 AUC/Brier；默认的 t2dm 二元表型是**基线已确诊
ICD10 T2D**。`--prevalence cohort` 用评分和协变量筛选前的各祖源 UKB 队列比例作为
工作假设，具体 K 和分母显示在报告的患病率表格中，不声称是普通人群患病率。可以指定
`--phenotype-col` 和 `--prevalence EUR=...,AFR=...,EAS=...,SAS=...`。
`t2e` 保留删失信息，报告 Harrell C，不把 incident event 指标转换为 Liability R²，
因此不需要也不使用 `--prevalence`。

缺少 `csx.auto/csx.meta` 时默认明确报错；仅在显式提供 `--allow-missing-scores` 时
生成标明缺项的部分报告。所有误差线是固定 OOF 预测的配对个体 bootstrap 区间，
不包括 GWAS 和模型估计的不确定性。
