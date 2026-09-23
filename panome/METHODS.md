# 方法、估计目标与验证边界

## 研究对象与拆分

每一行是一个基线个体。首先排除基线已患病、无有效随访和不满足 QC 的样本；QC 阈值、缺失填补、winsorization、批次残差化、标准化均由 build 学习。外层约 50% test 完全保留，另 50% 分成约 30% build、10% tune、5% calibration-fit、5% calibration-audit。后两个 5% 来自原 calibration 按 ID/亲缘组的结局盲拆分。没有把“50% reference”全部用于最终 donor；其中一部分必须用于开发验证，避免自我证明。

所有 build 内部神经训练还划出约 20% 作早停；这部分不进入梯度训练。主模型冻结后，已知 10 年结局的 build 人都可加入最终候选 bank，包括神经内部早停者。内部早停损失只用于选择 epoch，不作为泛化成绩。预处理/token 划分使用 build X；严格 OOF 评分则对每个外层 OOF fold 重新学习预处理及 token。

split 文件可以预先指定 build/tune/calibration/test。未指定时使用事件分层，因此 test Y 仅用于初始分层，manifest 明确记录这一点。验证测试标签隔离时固定同一个拆分。

## 为什么与基因型 imputation 类似、又不能照搬

参考 haplotype matching 的价值在于借用真实群体结构。蛋白丰度不是二倍体的相位序列，没有与 phasing 对应的必然分解，也没有基因组邻近关系保证可以“复制”一段连续蛋白模块。同一组可测 X 对应的未来疾病 Y 仍有随机性、未测量因素和随访删失。

因此同时实现三种不同输出：

1. literal COPY1：最近 reference 的已观察 0/1 结局，忠实检验原提案。
2. soft borrowing：多个真实人的 Y 按 learned matching 权重组合；这是主 reference 风险的来源。
3. molecular reconstruction：遮蔽若干实际测得的蛋白，用 references 或 decoder 重建，检验“类似 imputation”的信息借用是否成立。

COPY1 的 0/1 不是一个人的确定未来。所有概率输出在 0–1 之间；“150% 风险”若指相对风险 1.5 倍，必须另定义参照群体，不能当成 150% 的发病概率。

## 编码器与损失

每个保留蛋白有各自的数值 embedding 参数和缺失 embedding 参数，先在互不重叠的模块内聚合成 token。数据模块按 build X 的 PCA loading 聚类；外部模块采用明确的 assay 映射，未映射 assay 自动分组。每个蛋白恰好进入一个 token，不默默只取 top 100 proteins。token 内聚合是有损表示，应通过重建与对照检验。

每层显式计算 Q/K/V 的 scaled dot-product multi-head self-attention，再经过残差、pre-LayerNorm 和逐 token FFN；最后 CLS 输出个体表示。不同层仍逐层计算，同一层的 token 可以并行更新。模块的顺序没有基因组位置含义，因此不加入序列位置编码，而使用蛋白身份和模块身份 embedding；MLP 对照使用相同 token 输入、训练划分和损失，但不同参数量会写入 training summary。没有把所有人的 Y 或待预测人群作为 Transformer 输入。没有额外下载外部蛋白预训练权重。

预训练仅遮蔽 observed assays，在被遮蔽的位置计算重建 MSE。联合训练损失为：

`L = IPCW_BCE(borrowed risk, Y) + direct_weight * IPCW_BCE(auxiliary head, Y) + reconstruction_weight * masked_MSE`。

默认权重为 1、0.3、0.2。独立 `direct` 模型去掉 retrieval loss；`no_pretrain` 保留相同监督训练预算而跳过预训练，因此总更新次数较少，不能说是严格计算量匹配的对照。`ssl` 使用预训练快照的欧氏距离而不使用未训练的疾病 head/gate。

训练 donor 的 **row embeddings** 由当前网络每个 epoch 刷新，缓存后 stop-gradient；当前 batch query 保持梯度，Q/K 投影矩阵在每个 batch 均接收梯度。这是交替更新的近似，不是每个 batch 同时反向传播所有 donor 的编码器。query 所属整个 inner fold 的 donor 被屏蔽；亲缘组同 fold，因此自己的 Y 和亲属的 Y 都不作为该 query 的直接 value。不同 fold 共享训练网络参数是通常的监督学习，不应将训练预测叫 OOF。

## 个体匹配与可重构风险

主模型采用真正的 **多头 Q/K/V reference cross-attention**。每个人的 `z` 是归一化 CLS 投影，先乘 `sqrt(width)` 恢复尺度；每个 head 用不同的可学习投影产生 `Q_i^h`、`K_j^h`：

`S_ijh = (Q_i^h · K_j^h) / [sqrt(d_head) * learned_temperature_h * tune_temperature]`。

query 的 `g_ih = softmax(gate(z_i))` 为 head 权重。首先屏蔽同一 ID / 同一亲缘组，再对所有可用 references 按 donor 维度做每 head softmax。将各 head 的概率按 g 混合，选择共同 top-k reference 集合。这个离散检索步骤本身不可微；已选集合中的 soft attention、Q/K、gate 和 query encoder 可微。

对于选中集合 J，每个 head 单独重新归一化：

`A_ijh = softmax_j(S_ijh + log(IPCW_j * class_correction_j)), j in J`；

`head_value_ih = sum_j(A_ijh * observed_Y_j)`；

`w_ij = sum_h(g_ih * A_ijh)`。

这里 **Value 就是实际观察到的 reference Y**。没有把 outcome 输入 Q/K，也不把预测的 Y 伪装为真实 donor 的 Y。不同 head 可以对同一候选集里的不同人赋予较大权重。每 head 的分布先归一化，再混合，不是把几个距离简单相加后称为 multi-head attention。

这是为可解释风险借用设计的 Transformer 变体。个人内部是标准 Q/K/V self-attention + FFN block；跨个人的 Value 不做任意线性变换，head 输出也使用非负凸组合，**有意区别于原论文的任意输出投影**，以保留真实结局的贡献恒等式和 0–1 概率。它不是原机器翻译 Transformer 的完整复现。

仅 horizon 结局可识别的 build 样本可作 label donor。IPCW 是校正删失的正测度；病例富集 panel 的 class correction 只能修正类别比例，不能修正类别内部的非随机筛选。主模型 support 的 `nearest_distance` 是最高匹配分数 donor 与 query 的 embedding 平方距离，而不是未归一化 Q/K logit 的负值；它不声称该 donor 同时是欧氏空间最近者。

令 `ESS_i = 1 / sum_j w_ij²`，build 全局 IPCW 风险为 `pi`，先验强度为 `alpha`：

`p_raw(i) = [ESS_i * sum_j(w_ij * Y_j) + alpha*pi] / (ESS_i + alpha)`。

即 `p_raw(i) = sum_j(c_ij * Y_j) + c_prior_i*pi`，所有系数非负且总和为 1。随后使用 calibration-fit 的单调 logistic 校准器。`individual_risk_decomposition.csv` 检查贡献总和与最终概率完全一致；不把非线性校准后的风险强行拆成相加的蛋白因果效应。

COPY1 固定 k=1、alpha=0、不做校准。soft borrowing 的 k、温度和先验强度仅在 tune 选择。同一 panel 的 encoder、距离、随机 donor、Y 扰动对照分别保留。`uniform_same_panel` 独立训练固定均匀 self-attention，仍保留 V/输出投影、FFN、相同深度和跨人 attention；`metric_same_panel` 独立训练原多头平方距离匹配；`equal_heads_same_panel` 冻结主 encoder 后令跨人 head gate 均匀，其 k/温度/先验和校准仍分别在开发集确定。检索排除相同 ID、相同 family；可用参考人少于 k 时在该 query 内降低有效人数并给予零 padding 权重，不受其他 query 批次组成影响。最终 ESS 不足会被 support 规则拒绝。

## Reference 的拟合质量与覆盖

默认采用 3 次重复、5 折 OOF 的 elastic net 与小树模型，计算真实标签 log loss 相对 fold 内常数风险的增益，并要求增益和重复稳定性达到阈值。`--quality-teacher neural` 使用 fold-specific Transformer outcome borrowing；`all` 平均三种 teacher 的 OOF 概率。每个 OOF 神经网络的预处理、分组、早停均只在该 fold 的训练部分发生，OOF 查询者的 Y 从未进入其 donor bank。

主 panel 在可靠候选中按结局分层、质量加权的分子覆盖选择真实人。默认保留接近已知结局 bank 的病例比例，同时确保有病例和对照。原始全局 top-fit 不强行修正，允许暴露其退化。病例富集、只按多样性、随机及全部 bank 是敏感性方案。

`utility` 另外利用 tune 对象，评估删去某个 donor 后 loss 如何变化，并对低曝光 donor 的估计向 0 收缩，再用于覆盖选择。它是使用 tune Y 的开发选择，不能声称是独立 OOF 质量证据。所有最终评价依赖未参与这些选择的 audit/test。

主 panel 候选不足时保留实际数目；少于两人时显式标记 diversity fallback；readiness 不通过。不会把更大 panel 的表现悄悄当成“100 人成功”。

## Mosaic、支持规则及解释验证

模块 mosaic 是额外的简单模型：每个模块在同一 100 人中按该模块的标准化蛋白距离直接借用 Y，再用 tune 拟合非负、和为 1 的模块混合权重。它没有伪装成与主 Transformer 等价的模型，模块距离仍依赖填补后的测量空间。

支持规则使用 matching distance、测得蛋白的 reconstruction RMSE、ESS、缺失比例、未知类别；另外用 tune 上的误差拟合一个小型 expected-error 诊断模型，阈值由 tune 冻结。独立 calibration-audit 检查整体 Brier/log loss 是否优于冻结的常数风险、AUC 是否高于 0.5、覆盖是否足够。不会用 test 调规则。

支持不是置信区间。audit pass 也不是经置信区间验证的个人可靠性保证。test 中分别汇报 accepted/rejected 人群的事件率、同一人群上的所有模型表现、误差分层和性别/年龄/中心亚组，防止靠排除高风险者制造“准确”的表象。

解释检验包括：

- 原始贡献恒等式、donor Y 翻转的解析影响。
- 固定检索邻居集合，删除一个 donor 并重新归一化，重新计算 ESS、先验和校准风险；这不同于重新搜索第 k+1 近邻。
- 相同 query、相同 masks 下比较 reference/随机 reference/decoder/PCA/零 build 均值的 RMSE 和逐 assay R²。未测量的 donor assay 不作为已测值，分母按 assay 的可用 donor 权重计算；无 donor 值回退到零中心均值。
- 所有模块逐一遮蔽的敏感性、参考人 Jaccard 及风险变化；遮蔽模块可能属于分布外扰动，不能解释为治疗效果。
- 同风险但分子组合不同的个体对、重复 seed 的 panel Jaccard 和风险一致性。

新增 `attention/` 输出逐人、逐层、逐 head 的实际 softmax 矩阵及 reference head 贡献。采样仅由 ID hash 决定。矩阵导出使用 eval 模式，代码核对显式矩阵路径与部署 SDPA 路径的输出一致性。完整张量有明确轴名，不能把 head 轴误作人群亚型。

解释检验还包括：冻结参数时将 self-attention 变成均匀/仅自身，分别只修改 query 或同时重新编码 bank；固定邻居集合删除某个 cross-attention head；遮蔽 rollout 得分最高、最低和随机模块。attention rollout 只是平均 head + 残差的传播启发式，忽略 FFN、LayerNorm 和 Value/输出投影，**不等同于实际风险归因**。所有模块遮蔽敏感性与 rollout 的逐人相关也单独输出。随机模块与 top/bottom 匹配模块数量，不匹配 assay 数量，模块大小差异是混杂因素。推理时破坏结构的实验可能导致分布外输入，其敏感性不能替代重新训练对照的泛化成绩。

模块 attention 或权重不是机制发现。prototype ID 可读不等于解释可信。尤其当校准器把 reference 系数压到零时，必须承认模型不再使用这一证据来源。

## 终点、选择与局限

当前只实现固定 horizon 的死亡删失 net risk。病例是 horizon 前发病；对照是随访超过 horizon 而未在该窗发病；此前被删失且未发病者不伪装为健康标签。reverse Kaplan–Meier 只由 build 拟合；IPCW 假设独立删失。没有实现条件删失模型、竞争死亡 CIF、多时间点 hazard 或外部中心迁移自动证明。

最终图表排名是探索性的。选择“最佳 approach”要结合区分度、校准、coverage、重建增益、reference fidelity 和重训稳定性；不能看一次 test 排名后继续调参再把同一 test 称为外部验证。成对 bootstrap 区间只条件于当前拟合，且多重比较未校正。重复训练 SD 不是独立队列的标准误。

与 v4 相比，本版改变了默认残差化、校准/audit 拆分和模型结构；应依靠同一 v5 运行内部对照及显式残差化敏感性来判断增益来源。没有实现基因型意义上的 phasing，也没有证明给每个人估出了独立的 causal beta；共享训练模型仍来自群体，个体化体现在参考集合、权重、风险分解和适用范围。

## 相关原始来源

以下用于方法定位，不把本代码宣称为这些论文的复现：

- Vaswani et al., **Attention Is All You Need**：[原文](https://arxiv.org/html/1706.03762v7)。self-attention、Q/K/V、多头、残差和 FFN 的依据；原始 decoder 在生成时仍自回归。
- Gorishniy et al., **TabR: Tabular Deep Learning Meets Nearest Neighbors**, ICLR 2024：[会议原文](https://proceedings.iclr.cc/paper_files/paper/2024/hash/4ef594af0d9a519db8fb292452c461fa-Abstract-Conference.html)。从训练对象的特征和标签检索信息，说明 retrieval 本身已有系统研究。
- Gorishniy et al., **Revisiting Deep Learning Models for Tabular Data**：[作者论文](https://arxiv.org/abs/2106.11959)。FT-Transformer 和 MLP/ResNet 基准提示不能预设 attention 必胜。
- Somepalli et al., **SAINT: Improved Neural Networks for Tabular Data via Row Attention and Contrastive Pre-Training**：[作者论文](https://arxiv.org/abs/2106.01342)。与 row attention 相关；本实现刻意使测试 query 之间互不影响。
- Qiu et al., **ProLM: a plasma proteomics pretrained model for the general population**：[UK Biobank 论文记录](https://www.ukbiobank.ac.uk/publications/prolm-a-plasma-proteomics-pretrained-model-for-the-general-population/)，DOI `10.1038/s41467-026-75507-6`。蛋白自监督预训练已有相关工作；本包不使用其权重、不借其结果宣称本方法有效。
- [PyTorch 官方安装说明](https://pytorch.org/get-started/locally/) 与 [Transformer 实现资料](https://docs.pytorch.org/tutorials/intermediate/transformer_building_blocks.html)。
