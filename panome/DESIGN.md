# Panome 设计依据

## 原问题与研究假设

原讨论：[探索Panproteome📍](https://chatgpt.com/c/6aa01c8e-2ce0-83ec-981e-fdca59accfbd)。将分析单位从蛋白列扩展为个体及其分子环境，寻找“相似 CAD 风险可以对应不同分子组合”的证据。

参考：[Zhang et al., Proteomic health archetypes identified in disease-free adults enable risk assessment for diverse chronic diseases, Genome Medicine 2026](https://doi.org/10.1186/s13073-026-01696-w)。本实现借鉴其疾病前的样本相似图与状态迁移思路；引入 AE 与严格独立测试评估。没有复现其完整基线疾病排除清单、蛋白比值迁移分类器、GWAS/MR 或外部队列，因此不使用论文的状态名称，也不宣称重复得到相同原型。

需要纠正原沟通中的过强推论：PWAS 是关联分析，并不假设疾病只有一种原因；AE 及状态聚类也不自动识别因果结构。单次 baseline 分子数据不支持真实的个体疾病轨迹或 molecular ancestry。

## 与原 Step 1–6 的对应

1. 读取 **person × omics**，保存唯一 person ID，按实际原始列数处理，而非硬编码 2,900。
2. 先按 person（或亲缘组）划分，再在 train 学习缺失过滤、0.5%/99.5% winsorization、median 填补、年龄/性别/技术项 ridge 残差化及标准化。
3. 先学习 AE，再在 AE 空间构建个体 kNN 图，避免在 3,000 维直接构造完整 N×N 距离矩阵；保留稀疏邻接。
4. 同时得到 AE 连续坐标、Leiden 状态与图的 diffusion 坐标。test 只通过 train 邻居和冻结模型投影。
5. 固定这些表示后，训练并验证 Cox 生存预测。Y 不进入 AE、图权重、状态数选择或状态命名。
6. 同一 split、相同临床协变量下比较 PWAS score、PCA、elastic-net Cox 和 panome 模型，再解释个体分子组合。

原步骤“先图后 embedding”在此拆成 **AE 表征 → 图 → diffusion embedding**。AE 降噪用于距离计算，diffusion 刻画人群局部结构；两者的贡献可分别检验。

## 训练/验证/测试边界

- 60% train、20% validation、20% test；默认按新发事件分层随机划分。事件仅用于划分和下游预测，不用于无监督学习的损失或结构选择。
- 所有分子预处理参数只从 train 拟合。样本 missingness 用 train 确定的特征集合计算，并保持原 split。
- AE 用 observed-entry MSE：原始缺失条目不作为重构目标；随机 mask 输入 10% 形成 denoising。AdamW，validation 重构损失 early stopping，保存最佳模型。
- 固定默认结构 `P → 512 → 64 → 20`，decoder 对称。不给某 latent 轴预设“炎症”含义。latent 坐标具有旋转/非可识别性，需后续 feature/pathway 解释。
- PCA 同样只用 train 拟合，提供线性表征对照。
- Leiden 在 train 的 AE 坐标上运行。候选 resolution 按最小社区规模及最大状态数约束后，以 train silhouette 选取；不看 Y。不存在合格多状态解时明确报错，不强行输出五类。
- validation 同时用于 AE early stopping 和各生存模型的正则强度选择，test 始终留出。默认不做 train+validation refit，从而保持图、状态与模型定义一致。
- test bootstrap 是已训练模型下的条件不确定性；不是完整模型重训 bootstrap。

## 图、连续状态与新个体

训练 AE 坐标按 train 的坐标 SD 标准化后计算欧氏 kNN。权重为 `exp(-distance / kth_neighbor_distance)`；有向边以最大值对称化。训练点去掉 self-loop。N>5,000 使用 NNDescent；小数据使用精确 sklearn kNN。

状态编号按训练社区人数降序确定，S1 仅是最大社区的数值参照，没有“健康”或“低风险”语义。留出者的 `state_weight_k` 是其 train 邻居权重落入状态 k 的比例，最大者为预测状态；它不等于炎症等生物机制的真实比例。训练点的硬状态使用 Leiden 标签，权重仍使用 leave-self-out 邻居计算。

以训练邻接归一化矩阵求非平凡正特征值/特征向量，建立 time=1 diffusion 坐标；留出者用相对于 train 的随机游走权重做 Nyström 延拓。多连通分量会造成谱退化，因此输出分量数并要求审查。

训练点 nearest-other distance 的 99% 分位数为分布外提示阈值。它只是分子空间的距离标记，不是保证检测到所有 batch/platform shift 的检测器。

稳定性首先检查 10% graph-edge dropout 后重聚类的 ARI。这不覆盖 AE 重训与样本重抽样；发表前应重复 random seeds / family splits，并比较共同样本的 ARI 和 out-of-sample predictions。

## 疾病预测与传统比较

九组模型：clinical；clinical+PWAS score；clinical+PCA；clinical+AE；clinical+hard state；clinical+soft state；clinical+diffusion；clinical+AE+soft state；clinical+all-features elastic-net Cox。

PWAS 基准是在 train 经相同分子校正后的矩阵上，按单变量 Breslow Cox 估计系数和 Wald P 值、BH FDR，选训练 P 值最小的固定 top 30 构造加权和。后续 Cox 再加 clinical covariates。它是**marginal prediction benchmark**，不是全面临床校正的病因 PWAS。所有模型在相同测试集比较，不用 test 选 top proteins。

默认年龄校正和临床年龄项都是线性的，不能保证消除非线性年龄效应；正式临床比较宜增加年龄样条及更完整风险因子。

低维 Cox 在 validation 上从 alpha=0.1/1/10/100 选择；全特征 elastic-net 用 l1_ratio=0.5 的 alpha path。生存终点保留截尾时间，不把随访不足的无事件者都当成固定期限阴性。

评估：Harrell C、5/10 年 IPCW time-dependent AUC、Uno C、Brier score、预测风险分位组的 KM net-risk calibration；时间范围不支持的指标写入 limitations。IPCW 测试随访截断在训练的 99% follow-up 分位数以内，避免超出训练删失分布支持。死亡按截尾，因此所有概率都是 net risk。

独立 test 的状态 Cox 是探索性关联，输出 HR/95%CI 和 Schoenfeld PH 检查。固定 2 年 landmark 排除在此之前事件或失访者，使用已有 baseline 预测分数评估剩余随访，避免把排除的早期病例重标成对照。

## 如何体现 mosaic / 个体异质性

仅以 `Cox(AE(X))` 预测，终层系数虽然固定，但 encoder 非线性，因此原分子尺度上 `∂log h_i/∂X_ij` 可以因人而异。加入状态权重可反映局部人群环境。输出三层证据：

1. person × AE/diffusion 坐标：连续状态。
2. person × state weights：邻域混合状态。
3. person × molecular contributions：AE-Cox integrated gradients，相对于预处理后的零向量，固定临床项。输出每人的 top 20 绝对贡献及完整积分和误差，不把模型归因当作因果效应。

分子状态表型与具体蛋白/通路的关联需要独立注释、外部队列与时间信息支持。当前不硬编码 LPA/APOB、IL6/TNF 或肾脏衰老状态名。

Mixture-of-experts 是后续可检验扩展：如果每个状态有足够 CAD 事件，可研究状态依赖的风险映射。但在第一版就加入多个疾病监督专家会提高模型自由度并混淆“无监督发现”的含义，因此本版先提供冻结状态 + AE 的可审查基础与非线性个体归因。

## 正式研究需要补齐的验证

预先定义 disease-free endpoint panel；核实分子采样日期与结局数据覆盖；处理 PPP 选择性抽样及亲缘；临床强基准；plate-only 对照；多次完整训练/AE seed 稳定性；外部队列及平台校准；竞争死亡的累计发生率；必要时非比例风险；独立 pathway 富集与疾病分型验证。横断面 latent manifold 不自动等于病程。
