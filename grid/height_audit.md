# GRID height 结果与代码审计

审计日期：2026-10-07。

本文件保留下载包的审计记录及原始来源链接。其中 `0pca.sh`、`1csx.sh`、`2disco.sh`、`f/1csx.py`、`f/2disco.R` 在当前目录分别命名为 `0.pca.sh`、`1.csx.sh`、`2.disco.sh`、`f/1.csx.py`、`f/2.disco.R`；引用行号对应原审计快照。

## 核心结论

当前上传结果中的 **EUR PRS-CSx Prediction R² 为 0.30234497861424，DiscoDivas-tuned 为 0.297219924901021**。后者的点估计下降约 0.005125，约为前者的 1.70%。AFR、EAS、SAS 也都没有观察到 Disco-tuned 相对四分数组合的点估计提升。这个结论来自已有工作簿和 HTML；本次没有重新运行真实 UKB 数据，也没有修改远端仓库。

0.302345 是特定协变量调整和 OOF 评估下的预测相关系数平方。它不能直接表述为“PRS 单独解释了总身高方差的 30.2%”。这个数量级本身也不能证明程序出错：PRS-CSx 原论文补充表 11 中，EUR 的对应 height 数值为 0.3187356909；原论文与本地分析在发现 GWAS、SNP 覆盖、协变量和评估设计上均需进一步对齐。

正式发表的 DiscoDivas 论文没有报告 height 实验，因此目前不存在一个可以直接拿来对照的“论文中 DiscoDivas 的 height R²”。已有 GRID 代码保留了作者的核心插值算法，但本地 untuned 分数和 Yeval 的折内 tuned 实现对应不同的锚点构造，不能合称为完全复现了论文实验。

本报告依据：代码快照 `jielab/scripts/grid/`，结果快照 `jielab/analysis/grid/`，PRS-CSx 原论文和官方补充工作簿，以及 DiscoDivas 正式论文和作者仓库。末尾列出原始链接。仓库只上传了部分体积较小的结果文件；未上传原始后验目录、完整个体预测或大型基因型文件，不等于相应分析没有运行。

## 1. 现有结果实际上是什么

### 1.1 四个主要目标群体

以下点估计读取自 `analysis/grid/Yeval/height/Yeval.combined_scores.xlsx` 的 `performance` 表。差值和相对变化从未四舍五入的单元格数值计算。

| 目标群体 | N | PRS-CSx 四分数组合 | DiscoDivas-tuned | 差值（Disco − CSx） | 相对变化 |
|---|---:|---:|---:|---:|---:|
| EUR | 462,511 | 0.302344979 | 0.297219925 | -0.005125054 | -1.70% |
| AFR | 8,764 | 0.086167832 | 0.084018591 | -0.002149240 | -2.49% |
| EAS | 2,656 | 0.239759269 | 0.217341124 | -0.022418145 | -9.35% |
| SAS | 11,249 | 0.208207521 | 0.197378254 | -0.010829267 | -5.20% |

这四组结果支持“本次分析中 Disco-tuned 没有提升点估计”，尚不足以进一步下结论说“DiscoDivas 显著更差”或“这个方法对身高无效”。要判断差异，需要同一批个体、同一次重抽样得到的配对差值区间。

HTML 还显示，EUR 中 `csx.EUR` 单套祖源权重的 Prediction R² 约为 **0.300600**，四套 PRS 组合后为 **0.302345**。组合增加约 0.001745，约为 0.58% 的相对增加。这说明在当前 EUR 目标样本和当前四套发现 GWAS 下，额外三套分数的增量有限；不能把“四祖源联合推断”和“四套 PRS 的表型组合”理解成必定显著提高 EUR 的预测。

HTML 中 OTH 的 N 为 127，PRS-CSx 和 Disco-tuned 分别约为 0.099099 与 0.151870。该小组的正向差值值得记录，但当前代码把低祖源分类置信度者归入 OTH，**没有把 OTH 定义为一个经过验证的 admixed 队列**。因此，这个结果不能作为已经验证了 DiscoDivas 混合祖源优势的证据。

### 1.2 工作簿中的全部 16 行组合模型结果

下表保留工作簿数值精度。区间是各方法本身的区间，不是方法间配对差值的区间。

| 目标群体 | 方法 | N | Prediction R² | 下限 95% | 上限 95% |
|---|---|---:|---:|---:|---:|
| EUR | PRS-CSx-auto-meta | 462511 | 0.299997975504701 | 0.297574922298897 | 0.302296989757441 |
| EUR | PRS-CSx | 462511 | 0.30234497861424 | 0.300067708273389 | 0.304693275837102 |
| EUR | PRS-CSx-fixed-meta | 462511 | 0.297729841459392 | 0.295301166679483 | 0.300020928904638 |
| EUR | DiscoDivas-tuned | 462511 | 0.297219924901021 | 0.294878979056281 | 0.29950335149029 |
| AFR | PRS-CSx-auto-meta | 8764 | 0.0892755194415303 | 0.0771992131628097 | 0.100359647254962 |
| AFR | PRS-CSx | 8764 | 0.0861678315647692 | 0.0746667244883744 | 0.0965468459680438 |
| AFR | PRS-CSx-fixed-meta | 8764 | 0.0890701173385803 | 0.0774163582380618 | 0.100355612405998 |
| AFR | DiscoDivas-tuned | 8764 | 0.0840185913211264 | 0.0726978672073572 | 0.0945474757124773 |
| EAS | PRS-CSx-auto-meta | 2656 | 0.207542615642992 | 0.181752953637718 | 0.233102534461906 |
| EAS | PRS-CSx | 2656 | 0.239759269280184 | 0.212634395917718 | 0.267968430716106 |
| EAS | PRS-CSx-fixed-meta | 2656 | 0.204151562984421 | 0.178834109928687 | 0.230026161275314 |
| EAS | DiscoDivas-tuned | 2656 | 0.217341124276985 | 0.189011511818207 | 0.245103289342724 |
| SAS | PRS-CSx-auto-meta | 11249 | 0.207524828497636 | 0.194213364627067 | 0.223776498894279 |
| SAS | PRS-CSx | 11249 | 0.208207521147393 | 0.195546895635352 | 0.224169361928997 |
| SAS | PRS-CSx-fixed-meta | 11249 | 0.206189471538178 | 0.192913288831431 | 0.222188533097478 |
| SAS | DiscoDivas-tuned | 11249 | 0.197378253886737 | 0.183869334498218 | 0.212573487796201 |

`PRS-CSx-auto-meta` 在 AFR 的点估计 0.0892755194415303 高于这里固定 phi、再拟合四分数组合的 0.0861678315647692。这也是现有数据不支持“使用目标表型组合一定优于 auto/meta”的一个具体例子。

来源：结果工作簿 [R1]，HTML [R2]。

## 2. Prediction R²、partial R² 和 full R² 要分开解释

令 \(y_i\) 为实际身高，\(\hat y_{0i}\) 为仅含协变量模型的 OOF 预测，\(\hat y_{1i}\) 为协变量加 PRS 模型的 OOF 预测。当前代码的主要指标为：

\[
R^2_{\mathrm{prediction}}
=\operatorname{cor}\!\left(y-\hat y_0,\;\hat y_1-\hat y_0\right)^2.
\]

也就是说，它衡量“去掉协变量预测后的表型”与“PRS 模型带来的预测增量”之间的相关程度。现有 height 报告的协变量为 age、sex、PC1、PC2。协变量模型、PRS 标准化和四分数回归均在训练折拟合，再预测留出折。

代码同时单独保存以下指标：

\[
R^2_0=1-\frac{\mathrm{SSE}_0}{\mathrm{SST}},\qquad
R^2_1=1-\frac{\mathrm{SSE}_1}{\mathrm{SST}},\qquad
\Delta R^2=\frac{\mathrm{SSE}_0-\mathrm{SSE}_1}{\mathrm{SST}},
\]

\[
R^2_{\mathrm{SSE,partial}}=1-\frac{\mathrm{SSE}_1}{\mathrm{SSE}_0},\qquad
\mathrm{SST}=\sum_i(y_i-\bar y)^2.
\]

EUR 的 PRS-CSx 四分数组合给出如下实际数值：

| 指标 | 数值 | 解释 |
|---|---:|---|
| Prediction R² | 0.30234497861424 | 协变量调整后的 OOF 预测相关系数平方 |
| Prediction r | 0.549859053407542 | 上述相关系数，保留正负号 |
| Baseline R² | 0.525077706502023 | age、sex、PC1、PC2 的 OOF SSE 指标 |
| Full R² | 0.668668076854101 | 协变量加四套 PRS 的 OOF SSE 指标 |
| ΔR² | 0.143590370352078 | 相对于总身高方差的预测增量 |
| SSE partial R² | 0.302344977942565 | 相对于协变量剩余 SSE 的改善比例 |
| RMSE | 5.32450534626537 | 按输入身高的原单位计算 |

在这个大 EUR 样本中，Prediction R² 与 SSE partial R² 极其接近。这种接近不能证明指标只是换了名字；二者定义不同。相关系数平方对预测的比例缩放不敏感，而且会消除符号，SSE 则会反映校准偏差，因此仍应保留 signed r、SSE 指标及 RMSE。小 OTH 组在 HTML 中的 Prediction R² 约为 0.099099，而 SSE partial R² 约为 0.055579，已经显示出二者可以明显不同。

证据位置：`scripts/grid/f/Yeval.R` 第 530–545 行计算相关系数平方，第 774–780 行计算其余 SSE 指标；HTML 的方法说明也给出了主指标公式。[S1, R2]

## 3. 数据拆分和可能的泄漏边界

### 3.1 已检查到的折内处理是合理的

`Yeval.R` 第 668–683 行按祖源分组随机分配五折。第 719–744 行用 `fold != k` 的样本拟合模型，并在 `fold == k` 的样本上预测；第 734–740 行用训练折的均值和标准差标准化 PRS。

Disco-tuned 的锚点也在第 617–658 行按外层折重建，使用 `eligible & fold != k` 的样本。因此，所检查代码中没有看到用外层测试折身高直接拟合四分数组合或 Disco 锚点的问题。

### 3.2 五折不等于家系隔离

上述旧评估代码的折分配是个体层面的随机分配，没有读取亲缘组或家系连通分量。如果输入样本仍包含相关个体，亲属就可能跨训练折和测试折。是否在上游已经去除了相关者，现有上传材料不能确认。

推荐在明确保留/排除样本后，固定同一份划分，并确保相关者的连通分量不跨集合。后续比较 PRSformer、PRS-CSx 与 Disco 时，测试集、表型可用性、协变量和亲缘处理都应一致。新的单次留出测试结果不能直接与旧五折全样本汇总结果计算方法增益。

### 3.3 发现 GWAS 与 UKB 是否重叠，目前无法证实

height 配置只给出了 `height.AFR.gz`、`height.EAS.gz`、`height.EUR.gz`、`height.SAS.gz` 的本地路径；缺少足以追踪具体论文、版本、组成队列及 UKB 纳入情况的完整 manifest。`n_gwas` 留空也不能直接说“没有样本量”：入口代码支持从标准化 GWAS 的保留 SNP 推导样本量，并允许覆盖。

下游 OOF 只能隔离本地模型拟合使用的表型。如果发现 GWAS 本身已经包含待评估 UKB 个体，或与其有未处理的近亲关系，下游五折不能消除这一来源的偏差。当前证据支持“重叠未知，需要核实”，不支持“已确认发生重叠”或“已确认独立”。

证据位置：`config.json` 的四个 GWAS 输入路径及样本量选项；`1csx.sh` 第 285–303 行的样本量处理。[S2, R3]

## 4. PRS-CSx 的推断和 SNP 覆盖

### 4.1 四套权重是联合后验，不是四次互不相干的 PRS-CS

PRS-CSx 利用不同祖源的 summary statistics 和相应 LD 参考联合推断，跨祖源共享收缩信息，同时保留各祖源的后验效应。之后再在目标个体上分别计算四套 PRS，并用目标训练表型拟合其组合。

本地封装中有三个容易混淆的名称：

| 本地输出 / 模式 | 实际含义 |
|---|---|
| `PRS-CSx` / populations | 固定 phi 的四套祖源后验分数，再在目标训练折拟合组合系数 |
| `PRS-CSx-auto-meta` / auto | 从 GWAS 学习 phi，并用后验方差的逆数对祖源后验效应进行 meta 合并 |
| `PRS-CSx-fixed-meta` / meta | 使用给定 phi，并进行同样的后验效应 meta 合并 |

这里的命名是 GRID wrapper 的命名。官方 PRS-CSx 的 `auto` 和 `meta` 是两个维度：auto 负责学习全局收缩参数；meta 负责合并后验效应。auto 本身并不自动把四套权重变成一套。上述单套 meta 分数在 Yeval 中仍会接受折内回归校准，用于保证评估口径一致。[S2, P1, P3]

### 4.2 当前 populations 分支固定 phi = 1e-2

上传配置明确为 `phi = "1e-2"`，`1csx.sh` 的参数说明也明确指出这是固定值，没有用目标表型调优。不能把它称为经过 phi 网格选择的最优 PRS-CSx。若要公平比较，应在训练/验证范围内比较合理的 phi 候选或 auto，锁定后再评估测试集；不能根据最终测试结果选 phi。

### 4.3 目标 BIM 在后验推断之前就限制了 SNP 集合

配置把目标列表设为 `/mnt/f/gen/ukb/37/hap/ukb_array.bim`，并使用 `snpinfo_mult_1kg_hm3`。`f/csx/parse_genet.py` 第 49–62 行读取目标 BIM，第 90–102 行把可用 SNP 和等位基因限制到目标 BIM、参考 SNP、GWAS 的共同集合。

这意味着当前使用的权重集合在推断时已经受到 `hap/ukb_array.bim` 的限制。即使 GWAS 原文件含有数百万 SNP，也不能据此推断当前 PRS 真的用了数百万 SNP。改变最后 PLINK 评分的基因型目录，不能恢复前面已经被过滤掉的变异；如果希望扩大到 imputed HapMap3 可用集合，需要提供对应的完整目标变异列表、确认 build 和等位基因，再重新生成后验权重。

目录名 `hap` 本身不是性能差异的解释。加性 PRS 主要使用剂量，是否相位通常不是这里的关键；需要量化的是实际变异覆盖、剂量质量、等位基因匹配和 LD 参考。当前 PCA 入口另外使用 `imp` 目录，不能把 PRS 权重的 `hap` 限制推广为“整个项目都只用了 hap”。

已直接数过上传的 EUR 五条染色体权重文件：chr13 为 4,812，chr14 为 4,477，chr15 为 4,541，chr21 为 2,103，chr22 为 2,567，合计 18,500。这个数字只覆盖所检查的五条染色体，不能外推成全基因组 SNP 总数。正式性能排查应输出每条染色体在 GWAS、参考、目标基因型、等位基因协调及最终评分各步骤的数量。

证据位置：`config.json` 第 24、29、502 行；`1csx.sh` 第 35–39 行；`f/csx/parse_genet.py` 第 49–62、90–102 行；PCA 使用目录见 `f/0.common.sh` 第 120 行及 `0pca.sh` 第 119–122 行。[S2, S3, S7, R3]

## 5. DiscoDivas 的两条本地路径

### 5.1 保存的 untuned 分数

`2disco.sh` 默认使用 1KG 参考中心、10 个 PC 和四组相同的质量参数 A=1。`f/2disco.R` 保留了作者方法的中心间距离矩阵求逆、个体到中心的距离插值、收缩及 PC 残差化。`f/disco/UPSTREAM.md` 指向作者仓库的 commit `ee8e5d996e6fc6ffaf04255e1cf0fa8ecce68524`，并记录本地修复了 PC 数目标签和最终输出 ID 对齐。

作者 README 要求输入 PRS 已在各自的 validation cohort 中调优，并使用实际 validation cohort 的 PC 中位数作为中心；示例 1KG 中心需要替换。直接把四套未完成这种目标调优的 CSx 分数和 1KG 中心送入插值，属于一个 untuned 诊断分支。可以用它理解行为，但不应据此宣称完整复现了论文的验证设计。

该分支对目标 PRS 做 PC 残差化和尺度处理时会使用目标样本的 X 信息，没有使用身高 Y；这应描述为使用目标分布的预处理，不能直接称为测试表型泄漏。[S4, S6, P4]

### 5.2 Yeval 的 DiscoDivas-tuned

Yeval 的 tuned 路径在每个外层训练集合中：

1. 分别在 AFR、EAS、EUR、SAS 训练样本中，以协变量和四套 CSx 分数拟合目标表型，生成四个锚点模型。
2. 用各组训练样本的 PC 中位数作为中心。
3. 从四组训练样本中均衡抽样，拟合 PRS 对 PC 的残差化和尺度参数。
4. 根据作者核心距离插值公式为个体组合分数，再在外层留出折评估。

因此当前 `DiscoDivas-tuned` 已经处理了前面 untuned 分支缺少目标训练锚点的问题，不能仅用“Disco 输入没调优”解释它的全部下降。它仍然是对作者核心方法的一种本地 OOF 实现；锚点队列、中心、质量参数、祖源分布等与论文实验是否一致，需要逐项核实。

另外，PCA 输出的 `fine_tune_eligible` 包含更严格的置信度和自报祖源一致性条件，但当前 `Yeval.R` 的锚点构建按 `target` 分组，没有使用该字段。应明确这是有意的队列定义，还是需要进一步缩小锚点训练队列；不要根据测试性能反向挑选锚点。

证据位置：`f/2disco.R` 第 31–43、169–179、210–231 行；`Yeval.R` 第 617–658 行；`f/0.pca.R` 第 444–464 行。[S1, S4, S5, S6]

## 6. “遗传距离”和 OTH 的实际定义

当前评估默认使用 10 个 PC，把个体到四个 1KG 参考中心的欧氏距离合成等权 RMS：

\[
d_i=\sqrt{\frac{1}{4}\sum_{k\in\{\mathrm{AFR,EAS,EUR,SAS}\}}
  \left\|\mathrm{PC}_i-\mathrm{center}^{1KG}_k\right\|_2^2}.
\]

这不是单独“到 EUR 的距离”，也不是“到实际 GWAS 发现样本的距离”，更不是推断一个人每个染色体片段祖源的 local ancestry。它是同一 PC 空间中的四参考中心距离代理。代码已经提供单独的 discovery 模式，要求记录发现中心来源、N_GWAS 以及 PCA 空间；缺少这些信息时不应把参考中心重新命名为发现中心。

当前祖源分类从四类中取最大概率，低于默认 0.90 阈值就标为 OTH。OTH 可以包含分布边缘、置信度不足或难以归入四类者，不能自动解释为 admixed。验证连续祖源方法的核心问题，需要一个预先定义且样本量足够的混合祖源目标集，并报告其祖源连续分布；仅分析四个主要群体和 127 个 OTH，不能完整回答这一问题。

证据位置：`Yeval.R` 第 106、205–240 行；`0.pca.R` 第 444–447 行；HTML 距离说明。[S1, S5, R2]

## 7. 与论文数字怎样比较

### 7.1 PRS-CSx 原论文的 height 数值

PRS-CSx 正式论文的官方补充工作簿给出以下结果。表 11 使用 UKB White British + BBJ 发现 GWAS；表 14 加入 PAGE。以下全部取各表的 PRS-CSx 列。

| 目标群体 | 当前 GRID，四分数组合 | 论文 Supplementary Table 11：UKB + BBJ | 论文 Supplementary Table 14：UKB + BBJ + PAGE |
|---|---:|---:|---:|
| EUR | 0.30234497861424 | 0.3187356909 | 0.3226759453 |
| AFR | 0.0861678315647692 | 0.0831843257 | 0.1071273685 |
| EAS | 0.239759269280184 | 0.2223786000 | 0.2258614402 |
| SAS | 0.208207521147393 | 0.2067800935 | 0.2121811292 |
| AMR | 本地没有同定义的可比组 | 0.2715423416 | 0.2858984561 |

这个表说明，本地结果并非各群体都普遍低于论文；AFR、EAS、SAS 与表 11 的相对关系也不同。不能只看到 EUR 的一个差距，就推断所有后验收缩或评分都失败了。

同时，上表也不是一个方法优劣检验。原论文使用与 UKB discovery 无亲缘关系的 unrelated UKB 目标个体，调整 age、sex 和前 20 个 genotype PC；在各目标人群中多次随机分成 50% validation 与 50% test，并汇总 100 次划分。当前报告调整前 2 个 PC，并采用五折 OOF。发现 GWAS 的版本与队列、变异覆盖、样本关系和指标估计方式都尚未统一，不能据此断言本地实现比论文方法低多少。[P1, P2]

### 7.2 PRS-CS 2019 的基准数据也不同

PRS-CS 原论文 Table 1 的 height 分析使用 Yengo 等人的 EUR GWAS（N=693,529），目标是 Partners HealthCare Biobank 的 3,957 人，模型使用共同可用的 750,888 个 HapMap3 SNP。Figure 2 报告 height 的预测结果；每次用约三分之一目标样本调参、三分之二测试，重复 100 次。它可以提供历史量级参照，但不能把该目标队列、SNP 集合或拆分设计视为当前 UKB 分析的条件。[P6]

### 7.3 DiscoDivas 正式论文没有 height

正式论文 DOI 为 `10.1016/j.ajhg.2026.05.006`。Figure 3 报告的连续性状为 BMI、DBP、SBP、TC、HDL、LDL 和 lg(TG)，Figure 4 涉及 CAD、DM2；没有 height。

因此，论文支持的是其实际评估性状和人群中的发现，不能挪用成 height 上必然提升的承诺。当前 height 的下降应按当前发现 GWAS、锚点、目标队列和插值行为来排查。[P5]

## 8. 区间、工作簿与可复核性

### 配对不确定性

当前代码说明各方法共用 bootstrap 抽样索引，默认 200 次。但已上传 HTML 同时说明旧 exporter 没有保留方法间配对差值的数值区间，相关图被保留。此次没有取得可重算的真实逐人预测，也没有从图像反推精确区间。

因此本报告只给点估计差异；不会依据两个方法各自 95% 区间是否重叠来判断差异显著性。后续需直接导出相同测试个体上的 paired bootstrap 差值，涉及相关个体时以家系/亲缘组为重抽样单位，并保留抽样设计说明。[S1, R2]

### XLSX 元数据缺陷不等于数据缺失

`Yeval.combined_scores.xlsx` 的工作表 XML 把 dimension 错写为 `A1`，但 `sheetData` 中实际有 A1:R17，即 18 列、17 行，包含表头和全部 16 行结果；autoFilter 的范围也是 A1:R17。直接用 openpyxl 的 read-only 模式迭代会误以为只有一个单元格，调用 `reset_dimensions()` 后可正确读取全部数据。

```python
import openpyxl

wb = openpyxl.load_workbook(path, read_only=True, data_only=True)
ws = wb["performance"]
ws.reset_dimensions()
rows = list(ws.values)
```

此外，该文件含有指向缺失 drawing/VML 文件的关系，普通加载路径可能报错。这里属于导出兼容性问题；本次已从原始 sheet XML 和上述 read-only 读法交叉确认数值完整。不能因此把这个工作簿判为“没有结果”。当前 exporter 是否仍会产生这些问题，需要用当前版本另行验证，不能仅凭旧文件归因。

HTML 大表中 EUR 的 N 显示为 462510，而概览和工作簿为 462511，源于 `Yeval.R` 第 1097 行把所有数值列用 `signif(..., 5)` 显示。两处不一致是显示精度问题，不能据此推断评估时使用了两个队列。[R1, R2, S1]

## 9. 推荐处理顺序

1. **先补 GWAS 来源和重叠证据。** 为每个性状、每个祖源记录原始论文/下载版本、发现队列、有效 N、UKB 纳入情况及可用 SNP 数。若存在 overlap，应更换独立发现 GWAS 或使用真正独立的目标样本。
2. **量化并对齐变异覆盖。** 明确当前 `ukb_array.bim` 限制下实际保留多少 SNP；如果转向 imputed HapMap3，协调坐标和等位基因后重新做后验推断与评分。保存每一步 SNP 数和翻转/排除原因。
3. **统一评估设计。** 固定目标人群、样本、表型、协变量及亲缘分组，训练和验证完成所有选择后锁定测试集。增加到 20 个 PC 可作为预先规定的敏感性分析；不要根据测试集改善程度选协变量。
4. **在验证范围内调参。** 为 PRS-CSx 比较固定 phi 候选与 auto；Disco 的质量参数、锚点定义和需要的预处理同样只在训练/验证范围内决定。四套组合、meta 和 Disco 都保留为有明确定义的比较方法。
5. **针对连续祖源目标设计验证。** 使用预先定义的 admixed 队列，检查锚点中心是否来自实际训练/验证群体且位于同一 PCA 空间，再报告连续分布上的表现。OTH 只是当前分类器的一个输出类别。
6. **报告可解释的同口径结果。** 同时导出 Prediction R²、signed r、ΔR²、full R²、RMSE 及配对差值区间。新 PRSformer 模型须与在相同留出样本上评估的 CSx/Disco 比较，不能与旧五折点估计直接拼接成改进百分比。
7. **最后修复结果展示与兼容性。** N 保留整数；修复工作簿 dimension 和 dangling drawing 关系；完整保存配对区间及可追踪的参数/输入标识。

前六项决定科学可解释性和公平比较，第七项影响可读性和读取兼容性。此次代码与合成数据验证不构成真实 UKB 性能实验；PRSformer 的实际身高结果必须等用户环境中的真实训练和独立测试完成后才能给出。

## 原始来源

以下链接对应本次所读材料；`main` 链接后续可能变化，行号指本次下载的代码快照。

- [R1] [精确组合模型工作簿](https://github.com/jielab/analysis/blob/main/grid/Yeval/height/Yeval.combined_scores.xlsx)
- [R2] [现有 height 报告](https://github.com/jielab/analysis/blob/main/grid/Yeval/height/report.html)
- [R3] [height phi=1e-2 配置](https://github.com/jielab/analysis/blob/main/grid/csx/height/phi-1e-2/config.json)
- [S1] [评估和 Disco 折内调优](https://github.com/jielab/scripts/blob/main/grid/f/Yeval.R)
- [S2] [PRS-CSx 入口及模式](https://github.com/jielab/scripts/blob/main/grid/1csx.sh)
- [S3] [PRS-CSx SNP/等位基因交集](https://github.com/jielab/scripts/blob/main/grid/f/csx/parse_genet.py)
- [S4] [本地 Disco 核心实现](https://github.com/jielab/scripts/blob/main/grid/f/2disco.R)
- [S5] [祖源分类与置信度](https://github.com/jielab/scripts/blob/main/grid/f/0.pca.R)
- [S6] [Disco 来源记录](https://github.com/jielab/scripts/blob/main/grid/f/disco/UPSTREAM.md)
- [S7] [PCA 入口](https://github.com/jielab/scripts/blob/main/grid/0pca.sh)
- [P1] [PRS-CSx 正式论文](https://www.nature.com/articles/s41588-022-01054-7)
- [P2] [PRS-CSx 官方补充工作簿（Supplementary Tables 11、14）](https://media.springernature.com/original/springer-static/esm/art%3A10.1038%2Fs41588-022-01054-7/MediaObjects/41588_2022_1054_MOESM4_ESM.xlsx)
- [P3] [PRS-CSx 官方软件与参数说明](https://github.com/getian107/PRScsx)
- [P4] [DiscoDivas 作者仓库](https://github.com/YunfengRuan/DiscoDivas)
- [P5] [DiscoDivas 正式论文](https://doi.org/10.1016/j.ajhg.2026.05.006)
- [P6] [PRS-CS 正式论文](https://pmc.ncbi.nlm.nih.gov/articles/PMC6467998/)
