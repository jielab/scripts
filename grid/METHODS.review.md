# GRID：科学前提、方法边界与验证设计

**当前 GRID 实现的是 PRS-CSx 之后的个体预测层：用注释拆分后验评分，通过监督距离寻找参考个体，再借用其交叉拟合残差。它没有重新拟合加入进化先验的贝叶斯 SNP 效应模型。** 本会话未完成真实 UKB 训练，不能宣称性能已超过 PRS-CSx、DiscoDivas 或 PRSformer。

## 1. 进化信息支持什么推断

控制频率后，较年轻变异平均具有较大表型效应，已有 [Gazal 2017](https://www.nature.com/articles/ng.3954) 和 [CoalNN 2023](https://academic.oup.com/mbe/article/40/10/msad211/7279051) 的证据；这属于效应分布规律，不能据此令每个古老变异的权重衰减。自然选择作用于繁殖适合度，不等同于现代疾病风险；晚发疾病、环境变化及多效性可使风险持续存在。祖先型 APOE ε4 仍与 LDL 升高相关，见 [2026 年原始人群研究](https://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.1012285)。身高也不宜用“致病性”描述。

现代 AFR、EUR、EAS、SAS 是同时代的祖源概括，不能排列为进化阶梯；历史包含分化和基因交流，见 [Ragsdale 2023](https://www.nature.com/articles/s41586-023-06055-y)。频率分化、变异年龄、选择证据应分别标记；频率或 LD 不能直接充当年代。由 [GEVA](https://journals.plos.org/plosbiology/article?id=10.1371/journal.pbio.3000586) 推断的 derived allele 年龄，也不一定属于 GWAS 的 effect/risk allele。

[GPN-Star，Nature 2026](https://www.nature.com/articles/s41586-026-11005-5) 使用跨物种比对和系统发育树预测功能约束；[PrimateAI-3D，Science 2023](https://pmc.ncbi.nlm.nih.gov/articles/PMC10713091/) 利用灵长类常见容忍变异与结构信息，但原文亦记录在人类中致病的例外及序列背景补偿。二者启发特征设计，不能证明人类变异越古老越安全，也没有直接证明 GRID 的个体预测收益。当前代码不调用这些模型。

## 2. 当前算法及解释对象

对祖源来源 \(a\)、注释模块 \(c\)，由原始 CSx 后验权重计算：

$$
S_{iac}=\sum_s\widetilde G_{is}\widehat\beta^{\mathrm{CSx}}_{sa}\mathbf 1(s\in c).
$$

\(\widetilde G\) 使用保存的训练期评分中心化。年龄模块含 young、middle、old、uncertain、unknown；分区之和重构使用同样中心化及缺失值规则重算的 CSx 总分，不要求等于原先未中心化的评分列。年龄改变特征分组，不翻转 beta，也不施加单调衰减。

基线 \(b_i\) 由四组 CSx 及协变量校准，连续结局用 ridge，二分类用 logistic；目标祖源样本不足时使用 pooled fallback。build 内按家系嵌套交叉拟合得到 donor 残差：

$$
r_j=Y_j-\widehat b_{-\mathrm{fold}(j)}(X_j).
$$

监督距离用这些残差学习 ridge 特征权重，经分块标准化及 PCA 压缩到最多 12 维。它使用开发集结局，不能称为无监督；待预测者的结局不参与距离学习。排除相同 ID、家系后，在最近 \(k\) 个候选中保留满足 caliper 的集合 \(N_i\)：

$$
w_{ij}=\frac{\exp[-d_{ij}^2/(2\rho^2)]}{\sum_{l\in N_i}\exp[-d_{il}^2/(2\rho^2)]},\qquad
n_{\mathrm{eff},i}=\frac{1}{\sum_f(\sum_{j\in f}w_{ij})^2},
$$

$$
q_i=\alpha\exp[-d_{i,\min}^2/(2\rho^2)]\frac{n_{\mathrm{eff},i}}{n_{\mathrm{eff},i}+20},\qquad
\widehat Y_i=b_i+q_i\sum_{j\in N_i}w_{ij}r_j.
$$

\(\rho\) 来自训练参考距离；\(k,\alpha,\mathrm{caliper}\) 在开发期选择。少于 3 位 donor 或有效家系数不足 2 时回退基线；\(\alpha=0\) 始终可选。二分类预测另截到 \([0,1]\)。

这解释了“哪些参考人、各自权重及残差、使预测改变多少”。人数须区分：全库排除同 ID/家系后、距离满足 caliper 的合格人数，以及其中实际用于借用的最多 \(k\) 位 donor；后者才对应权重与贡献列表。有效家系数另体现权重集中程度，不能与前两者互换。参考人不是生物学双胞胎，贡献分解也不等于疾病机制或因果效应。

## 3. 两个需要实测的假设

**假设一：** 有来源的年龄/选择注释，在频率及 LD 信息之外，能够改善预测相关特征或残差结构。功能注释用于 PRS 已有 [LDpred-funct 2021](https://www.nature.com/articles/s41467-021-25171-9) 与 [SBayesRC 2024](https://www.nature.com/articles/s41588-024-01704-y)；CoalNN 亦已有多群体年龄注释。创新需落实到增量信息及个体适配，不能只称“首次加入进化”。

**假设二：** 参考个体的局部经验能改善预测，并在不查看新人的结局时识别较低误差或较大增益者。祖源接近或匹配人数多本身不足以证明这一点：[Wang 2026](https://www.nature.com/articles/s41467-026-68565-3) 中，遗传距离的灵活拟合仅解释身高个体平方预测误差方差的 **0.51%**。

以下 2×2 对照已在代码中配置，尚未完成真实 UKB 实验。全局同类算法的有/无演化版本使用相同超参数网格；基础特征均含 CSx、匹配 PCs 及可用频率信息，全局预测另含共同协变量。

| 特征条件 | 全局预测 | 个体局部校准 |
|---|---|---|
| 基础特征，不用真实年龄 | `Ridge_no_evolution` / `HGB_no_evolution` | `GRID_frequency_only`；无频率输入时为 `GRID_no_evolution` |
| 基础特征，加真实年龄及可用演化注释 | `Ridge_evolution` / `HGB_evolution` | `GRID_evolution` |

四个全局模型另有 `_full_training` 重拟合版本。列间估计器及实际标签用途不同，应分别报告 build 与完整开发期的预算，不能将交互差值直接归因于匹配机制；selection/proxy 同时变化时，行间差值也不是纯年龄作用。`GRID_no_evolution` 对 `GRID_frequency_only` 检查频率增量，真实年龄对 `GRID_permuted_evolution` 则保持其余注释相同。

年龄置换只改变合格年龄标签，保留 unknown/uncertain 位置。当前**没有 LD 控制**：仅按染色体分层，在外部参考 AF 可用时增加 MAF 分层；缺少 AF 的分层须标为染色体内探索对照，不能称 MAF 匹配。一次置换不是 permutation P 值。正式年龄归因仍需 LD 匹配/分层敏感性分析及 SBayesRC 等注释方法基线。

## 4. 数据隔离、筛选与评价

默认完整队列约为 build 30%、模型调参 7.5%、筛选器训练 7.5%、阈值确定 2.5%、独立审计 2.5%、test 50%。**donor bank 仅从目标 build 建立，并非通过 discovery GWAS 获得的个体库，也不是全部开发半区。** 其中人员是否与 discovery GWAS 重叠仍需根据实际来源核对。 家系表有效提供时按家系隔离；无家系表不能声称无亲缘。

筛选器分别学习相对增益 \((Y-b)^2-(Y-\widehat Y)^2\) 与绝对误差 \((Y-\widehat Y)^2\)，二者回答不同问题。阈值独立固定后，再审计配对增益及低误差筛选是否获得经验支持；test 只应用冻结规则。内部审计不是个体准确性的保证，也不能在审计或 test 失败后重新挑 coverage。

必须核对每个实际 GWAS 文件是否包含目标 UKB test 及其亲属，记录版本、排除 UKB 情况和上游调参；目标内部重新划分不能消除 discovery 重叠。还应报告各方法实际个体结局预算，比较使用完整开发半区校准的 CSx 强基线。

连续性状的冻结预测 \(R^2\) 定义为：

$$
R^2_{\mathrm{pred}}=1-\frac{\sum_{i\in T}(Y_i-\widehat Y_i)^2}{\sum_{i\in T}(Y_i-\overline Y_T)^2}.
$$

test 结局仅用于评价。可以拟合测试集的校准截距和斜率作为诊断，但不能用所得系数重写预测或重新报告校准后的测试性能。残差平方相关 \(\mathrm{cor}(Y-\widehat Y_{\mathrm{cov}},\widehat Y-\widehat Y_{\mathrm{cov}})^2\) 是另一指标，不能替代上述误差定义；协变量系数也变化时，增量不能全称纯遗传贡献。

主结论先比较完整 test 的配对 MSE/Brier，并辅以 RMSE、MAE、校准及二分类 AUC、log loss、average precision。再在同一组预先选定者中比较全部模型，报告实际 coverage、病例覆盖和祖源分层。低 MSE 子集未必有较高 \(R^2\)，因为其结局方差可能更小；不能用 GRID 的选中子集与别的方法的全样本制造增益。
