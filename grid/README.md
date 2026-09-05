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

默认 GWAS：`/mnt/d/data.BIG/gwas/4grid/{trait}.{AFR,EAS,EUR,SAS}.gz`。
默认分析目录：`/mnt/d/analysis/grid`；可用 `--output-root` 覆盖。
其他选项见 `./grid.sh -h`。`--replace FALSE` 复用缓存；更换输入/配置时明确重建。

ARG 建树属于共享数据准备，在 `refGen.sh` 中完成；GRID 只消费已有结果：

```bash
bash /mnt/d/scripts/gu/arg.sh build --method needle \
  --dir-gen /mnt/h/ukbGen/37 --dir-pfile /mnt/h/ukbGen/37/hap \
  --map-dir /mnt/d/data.BIG/refGen/maps/GRCh37 --chr 22
```

默认 ARG 位于 `/mnt/h/ukbGen/37/arg/{argn,trees}`，可用 `--arg-dir` 覆盖。
全队列建树的资源需求应先通过小规模实验评估；模型验证应使用多个染色体。
用于性能评估的 GWAS 应排除 UKB，以避免样本重叠导致的乐观估计。
