# 从个人组合到 reference borrowing：Panome v5.1 的 attention 设计

本次改进保留原问题：**这个人的分子组合最接近哪些真实人？这些人观察到了什么结局？借用这些人的证据是否有支持？** 神经网络学习匹配规则，而真实 reference 的 Y 仍是主风险的 Value。

## 个体化的含义

线性加性模型的系数在人群中共享；交互模型、非线性模型和局部模型也可以使某个变量的关联依赖个人背景。Panome 的模型参数同样在人群中共享，个体化体现在内容相关的 attention、reference 集合、权重和支持范围，不能宣称摆脱了群体信息或为每个人识别了独立因果 beta。

## 为什么这次不只是改名

| 环节 | 上一版 v5.0 | 本次 v5.1 |
|---|---|---|
| 个人内部 | PyTorch TransformerEncoder，模块 self-attention | 显式 Q/K/V block、残差、pre-LN、FFN；部署 SDPA 与矩阵导出路径互相核对 |
| 人与 reference | embedding 分块距离，经 gate 合成一个分数 | 学习 Q/K 投影，每个 head 独立 softmax，再读取真实 Y；用凸组合融合 head |
| Value | 单一匹配权重加权 Y | 每 head 都有可核对的 Y 加权值，汇合后仍能精确分解 |
| 对照 | MLP、无预训练、无 retrieval 等 | 增加独立训练的均匀 self-attention、距离匹配，以及冻结主模型的均匀 head gate |
| 解释证据 | reference 贡献、模块遮蔽 | 加上真实 attention 张量、逐 head 贡献、移除 head、破坏 attention、attention 与遮蔽敏感性的对应 |

这不是首次提出 attention 或样本检索；有关方法已有 FT-Transformer、SAINT、TabR 等研究。值得探索的研究主张是：将可审核的真实人群参考、疾病结局借用、个体匹配支持与分子重建组合成系统，然后以独立数据检验。

## 第一层含义：理解个人内部的分子上下文

对每位受试者，所有 QC 保留的 assays 都进入模型。每个 assay 有自己的数值和缺失 embedding。先聚合成约 48 个模块 token，再加入 CLS token，经过默认两个 Transformer block。每个 block 的计算是：

```
U = LayerNorm(T)
Q_h, K_h, V_h = U W_Qh, U W_Kh, U W_Vh
A_h = softmax(Q_h K_h^T / sqrt(d_head))
T' = T + Dropout(Concat_h(A_h V_h) W_O)
T_next = T' + Dropout(FFN(LayerNorm(T')))
```

同样高的某个蛋白，在不同的其他蛋白背景下，可以产生不同的上下文表示。这是一种可学习的条件关系，不自动意味着生物相互作用已被证明。

模块是计算上的分层压缩，有损而非完整保存 3,000 个蛋白的所有两两关系。数据驱动模块不能自动命名为“炎症通路”。外部 pathway 映射可由 `--module-file` 提供；请审查映射及重叠处理。assay 身份、缺失模式和模块 embedding 提供身份信息，没有把 CSV 列顺序当作蛋白序列或基因组位置。

## 第二层含义：向真实人的结局提问

个人 CLS 表示形成 Query，build reference 的 CLS 表示形成 Key；每个 head 有不同 Q/K 投影。每 head 可关注不同的参考相似性，但不能预先把某个 head 宣称为脂质或炎症。

```
Query  = target 分子组合的投影
Key    = reference 分子组合的投影
Value  = reference 的已观察疾病 Y
```

每个 head 对选中 reference 单独做含 IPCW 正测度的 softmax，得到一组可审核的权重。head gate 是 target X 的函数。最后凸组合各 head，再加显式先验收缩和独立校准。校准后的概率不能按简单蛋白加法拆开，但原始风险有精确的 donor 分解。

跨个人这一层**不是**原论文完整的任意 Value/输出投影：这里有意保留 `Value=真实 Y`，输出必须为非负凸组合，以满足用户要求的 COPY / borrowing 可解释性。直接神经疾病 head 作为辅助损失和独立对照，不能悄悄混入主 reference 风险。`panome_clinical` 等混合模型单独命名。

训练会排除 query 所属 inner fold 的 donors，并保持亲缘组同 fold；测试者只访问冻结 build bank，不互相注意，也不访问自己的 Y。缓存 donor 表示每 epoch 刷新，Q/K 投影每 batch 更新。这是可扩展的近似优化，不声称对完整 bank 编码器做每 batch 全量反向传播。

## 把 explainability 变成可检验的问题

按以下顺序阅读新输出，避免只看漂亮的 heatmap：

1. `individual_risk_decomposition.csv`：donor 贡献加先验能否还原概率？
2. `test_reference_matches.csv.gz`：这些 donor 是谁？真实 Y、IPCW、每个 head 权重分别是什么？
3. `attention/reference_head_contributions.csv.gz`：每个 head 是否使用不同 reference？移除该 head 后风险是否改变？
4. `attention/self_attention.npz`：真实每层、每 head 的连接是否有个体差异，还是退化成近乎均匀？
5. `attention/intervention_summary.csv`：把 attention 置均匀/仅自身、遮蔽高低 attention 模块，会发生什么？
6. `attention/attention_vs_sensitivity.csv`：rollout 与真正模型敏感性是否对应？若不对应，保留这个负结果，不把 attention 排名当蛋白重要性。
7. `approach_comparison.csv` / `paired_contrasts.csv`：独立训练的完整 attention 是否优于均匀 attention、距离匹配、MLP；校准、覆盖、重建是否同时可接受？

rollout 只是传播启发式，未包含 FFN、LayerNorm 和 Value 投影的完整效应。删 head 的测试固定邻居集合；结构干预保持原校准器；模块遮蔽可能造成分布外输入。它们是忠实性和依赖性诊断，不是干预治疗效果或因果归因。

## 计算量不是研究深度的替代品

默认复杂度约为 assay embedding 的 `O(batch × assays × width)`、模块 self-attention 的 `O(batch × tokens² × width)`，以及 query-reference matching 的 `O(batch × bank × width)`。分批读取不需要完整人群 × 人群矩阵。层数和 token 数可以增加，但是否需要更大模型应由开发集学习曲线、稳定性和公平对照决定。

建议先跑默认完整流程及 `--quality-teacher all`，再用同一外层 split 做 seed、token 数和深度敏感性。评估一经看过，就不应在同一 test 上继续挑参数并宣称它仍是独立确认。特别是 100 个 reference 是否足够，需要与更多 reference 和 full bank 一起评估。

## 这次实际看到的界限

本次工程验证及本机修复见 `README.md`；尚未生成真实队列的新版模型结果。短预算样例中某些 head/attention 干预会改变风险，说明计算确实参与预测；但 attention rollout 高分模块并不保证更敏感，且独立 audit 可以拒绝参考库。代码保留拒绝和负结果，不自动将 Transformer 或 100 人方案宣布为最佳。

原始来源：[Vaswani et al., Attention Is All You Need](https://arxiv.org/html/1706.03762v7)，尤其 3.1、3.2、3.3 节。源码 `f/attention.py`、`f/neural.py`、`f/borrowing.py` 与这里的公式对应；其他方法来源见 `METHODS.md`。
