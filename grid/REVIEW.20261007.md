# GRID 代码审阅与新版说明 · 2026-10-07

审阅基线：[jielab/scripts@b849964](https://github.com/jielab/scripts/tree/b84996411653bc4e3b190d9cb997578396da4f29)。范围包括 `grid/` 全部入口与核心实现、PRSformer 适配器、旧 Yeval，以及作为设计参考的 `le8/f/c1.abm.py`。本次修改集中在 `grid/`，未修改 le8 的分析流程。

**结论：保留“有来源的人类变异注释 + 逐人匹配 + 冻结筛选”的研究路线，先修复数据身份、评价和解释上的缺口，再用真实固定测试集决定是否有增益。新版没有真实 UKB 性能结果，也没有重写 PRS-CSx 的贝叶斯推断。**

原代码已经具备不少正确的基础：共同 outer roster、开发集内部角色、家系 OOF donor 残差、真实年龄和 unknown/uncertain 分区、未知年龄位置固定的负对照、保存模型以及测试结局隔离。这次是在现有实现上补强，不能把这些既有设计全称为本次新增。

运行入口和完整参数见 [README.md](README.md)；科学依据、算法方程与文献见 [METHODS.review.md](METHODS.review.md)。

## 1. 已更正的问题

| 优先级 | 原问题及影响 | 本次修改 | 主要文件 |
|---|---|---|---|
| 高 | 已监督训练的 GRID 分数被旧 Yeval 重新随机分折后，不能成为这些新折的 OOF 输入；原训练可能见过新测试折的结局 | Bash 和 R 两个入口拒绝 `--grid-file`，引导使用 GRID 的冻结 test report | `Yeval.sh`、`f/Yeval.R` |
| 高 | ALT_AF 直接当作效应等位基因频率，会在 EA=REF 时方向相反；无标识的 AF 也可能含义不明 | 明确 EAF 优先；ALT_AF 根据 EA/REF/ALT 翻转；generic AF 需要频率等位基因标识，否则保留缺失；MAF 不充作 EAF | `f/0.common.py` |
| 高 | 一个既有 CSx 分数表和另一轮后验权重可能被同时输入 GRID，表面上性状和列名均一致 | 同锁生成分数来源 sidecar；GRID 核对分数哈希、各祖源权重哈希、联合推断签名、染色体集合与评分规则 | `1.csx.sh`、`f/0.common.py`、`f/grid.data.py` |
| 高 | PRSformer 新预测时，仅检查当前 split 标签不足以阻止把已训练对象重新标成 test | checkpoint 保存开发期 ID/家系注册表，预测时核对身份；旧 checkpoint 只接受其原有输入身份 | `f/3.prsformer.py`、`f/3.prsformer_data.py` |
| 中 | `age_usable` 可能把质量通过但区间跨箱、实际进入 uncertain 的年龄算作可用 | 分开 quality-pass 与 usable；后者与实际年龄箱一致；增加非零权重可用年龄的最低数量要求 | `f/grid.evolution.py` |
| 中 | 缺少参考 AF 时，年龄负对照事实上只能按染色体置换，不能声称已匹配 MAF；也未控制 LD | 增加 auto/chr-maf/chr-only 模式；严格 chr-maf 缺失即拒绝；报告实际方案、AF 覆盖及有效置换对比，明确 LD 未控制 | `f/grid.evolution.py`、`f/grid.data.py` |
| 中 | `Prediction_R2` 名称掩盖了实际计算的是残差平方相关，可能在方向或尺度错误时仍很高 | 区分 `predictive_R2=1-SSE/SST` 与 `residual_correlation_R2`；同步 GRID、PRSformer、Yeval 的输出与图注 | `f/grid.py`、`f/3.prsformer.py`、`f/Yeval.R` |
| 中 | 缺少直接的预测校准诊断；二分类完全分离时，数值优化器可能返回貌似有限的斜率 | 增加校准截距、斜率和状态；常量预测、缺少类别、完全/准完全分离不报告可识别估计；测试系数不改写预测 | `f/grid.py` |
| 中 | “有多少匹配者”最多只能到 k，不能回答全参考库有多少合格者 | 精确计算 `eligible_match_count`，同时保留实际 `matched_count`、截断标志、donor/family ESS；排除同 ID/家系 | `f/grid.abm.py` |
| 中 | 已有 donor 修正分解，但无法直接重构基线全部贡献 | 增加实际 pooled/祖源模型的逐特征基线贡献；连续尺度或 logit 尺度可重构，再接 donor 修正和 clipping | `f/grid.abm.py` |
| 中 | “相对 CSx 增益更大”和“本身误差更低”不是同一筛选目标；低误差选择缺少独立审计 | 保留原 gain gate，增加独立 low-error audit；报告 pooled/祖源结果、样本量、病例覆盖、配对家系区间 | `f/grid.abm.py`、`f/grid.py` |
| 中 | 仅有全局演化模型时，难区分注释增量与估计器复杂度 | 增加 Ridge/HGB 的 no-evolution 对照及完整开发集重拟合；同估计器、同超参数网格直接配对 | `f/grid.abm.py`、`f/grid.py` |
| 中 | PRSformer 与 GRID 协变量同名，不证明使用了相同实际数值 | 导出原始 `covariate.<name>`，按 eid 核对数值；同时核对共同 roster、test、家系、结局和单位 | `f/3.prsformer.py`、`f/grid.py` |
| 中 | 新人的预计算特征即使列名相同，也可能使用不同权重、PCA 或中心化；pickle 版本兼容也未充分约束 | 保存特征契约；predict 必须提供契约和精确输入哈希声明；加载前检查算法源码及运行环境 | `f/grid.py`、`grid.sh` |
| 低 | PCA 最后汇总阶段没有收到 score-dir，无法读回 `.sscore.vars` 做实际使用 SNP 数核查 | 接通该目录与最终 QC；保留与官方配套 loadings/centers 一致的投影公式 | `0.pca.sh`、`f/0.pca.R` |
| 低 | CSx preflight 拒绝只有 BETA/OR+P 的 summary statistics，但准备阶段其实支持 | 前后一致地接受 SE 或 P；有 SE 时优先使用其精度 | `f/1.csx.py` |
| 低 | 显式指定 Python 仍会触发默认 conda 激活，使独立可用环境被不必要地阻断 | 显式解释器直接执行，默认路径继续使用 grid 环境 | `grid.sh` |

例如 `REF=A, ALT=C, EA=A, ALT_AF=0.2` 时，正确 EAF 为 0.8；这类方向问题不能靠后续回归解释掉。另一方面，年龄模块分数重构的是按同样中心化与缺失规则计算的 `evo.total_<POP>`，不要求与原始未中心化的 `csx.<POP>` 相等。

## 2. 必须准确标注的比较对象

### PRS-CSx

`CSx_full_training` 的 full training 指目标队列的组合/校准回归使用完整开发半区，**没有重新选择 discovery GWAS 的 phi**。当前 `1.csx.sh` 的四祖源默认后验配置是固定 phi；auto/meta 是另行提供的配置。要主张超过充分调优的 PRS-CSx，需要在开发数据中选择 phi，并纳入 auto/meta 等适当基线。

来源 sidecar 证明“这套分数与这套后验对应”，不证明原始 GWAS 没有 UKB、没有祖源/坐标错误或满足研究设计。实际 discovery 样本重叠仍须逐个来源文件核对。

### DiscoDivas

`2.disco.sh` 默认将四个原始 ancestry CSx 分数、1KG PC 中位数及 quality=1 送入官方距离矩阵插值，并对全部目标 X 做不读 Y 的 PC 残差化。它采用官方公式，但没有在四个实际验证人群中拟合表型调优的 anchor PRS。

新版 GRID 将其下游回归标为 `DiscoDivas_calibrated_reference` 及 `_full_training`。原始 `disco` 输出列继续兼容旧代码。该标签假设分数来自现有 reference 流程；外部 prepared features 的实际来源仍由提供者负责。

旧 `Yeval::build_disco_folds` 已在每个 outer 训练折内拟合四个表型 anchor、训练期中心与 PC harmonization，确实可称为 phenotype-tuned 的自定义 CV 实现。但它固定 phi、quality、PC 数，每折每祖源至少需要 100 个 eligible training 样本，病例/事件也须充足；其 folds 按个人而非家系分组，不是 GRID 的固定 roster。因此这两个结果不能直接做配对差值，也不能据此宣称复现了整篇论文。

### PRSformer

仓库代码调用官方架构，但 height/LDL 的连续残差任务与 baseline T2D 的 logistic-offset 任务属于本地三性状适配，不等于原论文的大规模多疾病训练。新版统一其外部 test 和来源核查，并记录标签预算；没有在本次会话训练真实 UKB 的 PRSformer。

## 3. 新版 GRID 的实际研究对象

当前方法可表达为：

\[
\widehat Y_i=b_i+q_i\sum_{j\in N_i}w_{ij}
\{Y_j-\widehat b_{-\mathrm{fold}(j)}(X_j)\}.
\]

这里 `b_i` 是四个 CSx 加协变量的校准基线，距离和参考库在 build 集建立，donor 残差来自排除其家系的嵌套交叉拟合，`q_i` 根据距离和有效家系支持进行收缩。二分类还要作概率边界处理。最终 `GRID_policy` 由开发期冻结的 gain gate 和独立审计决定是否采用候选修正。

这仍然包含全局基线；增加的是逐人借用和选择。它也没有从 GWAS summary statistics 中恢复原始 GWAS 个体。参考库来自目标 build，其中人员是否与 discovery 重叠是另一项待核实事实。

基线项、距离块、donor 权重和实际借用量支持**预测过程解释**，不直接给出 SNP 或匹配人的因果效应。`eligible_match_count=1000`、`matched_count=32`、`family_ess=8` 分别表示全库支持、实际使用、有效独立支持，不能全叫“1000 个 twin”。

当前 exact 支持计数使用批量 KD-tree radius count，避免巨大 query×donor 矩阵。HNSW 可近似检索，但为了报告精确全库人数，仍需要额外的精确计数索引；这部分可能成为大样本计算瓶颈，本次没有 UKB 规模的速度或内存基准。

## 4. 筛选和消融的解释规则

默认完整队列约为 build 30%、模型调参 7.5%、gate 7.5%、阈值确定 2.5%、独立审计 2.5%、test 50%，家系完整留在同一角色。PRSformer 则共享同一 outer development/test，在开发半区内作 40/10 的训练/验证划分。完整开发期基线在 matching audit 后重拟合，不回流改变 gate。

- gain-selected 回答“相对于内部 build CSx 是否可能改善”；low-error-selected 回答“候选 GRID 本身的误差是否可能较低”。独立审计和最终 test 都应分别看这两个问题。
- 当前 release 仍由 pooled gain audit 决定；祖源审计行属于诊断，pooled 通过不等于每祖源都得到支持。
- 被选人群的 MSE 低，不等于其 R² 高；所有方法必须在同一组预先选定者上比较，再报告完整 test 的总体损失。
- family bootstrap 区间条件于冻结模型，不包含重新训练全部流程的不确定性；多个性状、祖源、coverage 的探索不能自动当作多重检验后的确认性结论。
- no-evolution 与 evolution 全局模型使用相同估计器/调参网格，但后者若同时增加 selection/proxy，差值是注释组合增量。真实年龄与置换年龄臂保留其余注释，仍未排除 LD 相关结构。

## 5. 进化前提需要怎样调整

现代 AFR、EUR、EAS、SAS 是同时代的祖源概括，不能排列成 AFR→EUR→EAS→SAS 的演化阶梯。应使用变异年龄、选择、频率分化等有明确来源的变量，分别表述它们。

“古老变异较安全”不能成为逐 SNP 单调衰减规则。选择约束、现代疾病风险、GWAS 效应、tag SNP 与因果变异是不同对象；derived allele 年龄也不一定是 EA 或 risk allele 的年龄。GPN-Star 和 PrimateAI-3D 的跨物种结果提供启发，未证明人类内的年龄注释一定改善 PRS。原始证据和限定见 [METHODS.review.md](METHODS.review.md)。

较稳妥的主张是检验：**可追溯的年龄/选择信息，在频率、LD 和已有统计效应之外，是否提供增量；这种增量是否帮助逐人的局部校准或预测误差筛选。** 当前代码把它写成可失败、可回退的假设，而不是预先写进 beta 方向。

## 6. 迁移和运行

在仓库 `grid` 目录、已安装共享环境和上游输入的前提下：

```bash
# 重新评分并补齐旧 CSx 分数的权重来源；不为此重跑 MCMC。
bash 1.csx.sh --traits height,ldl,t2dm --stage score

# 若要比较 reference Disco，使用更新后的 CSx 分数重新生成。
bash 2.disco.sh --traits height,ldl,t2dm

# 采用新目录，避免混用旧模型和缓存。
bash grid.sh all --traits height,ldl,t2dm \
  --split-group-file /mnt/d/files/ukb.family.tsv.gz \
  --cache-dir /tmp/grid-cache/review-20261007 \
  --out-root /mnt/d/analysis/grid/GRID-review-20261007
```

家系文件必须来自实际亲缘连通分量；没有时不应生成一份“人人无亲缘”的假数据来满足参数。可省略该参数运行，但结果只能说明未提供家系信息，不能宣称亲缘已经隔离。

真实年龄可通过 `f/grid.evolution.py download` 获取 GEVA，或使用 canonical 注释。参考 AF 齐备时推荐显式 `--age-permutation-mode chr-maf`；缺少 AF 时 auto 仍能运行，并记录该对照只覆盖到哪一层。`--min-age-variants` 需要根据预期注释覆盖和研究设计设定；默认 1 仅防止完全没有有效年龄。

需要同 roster PRSformer 时，在上面的同一配置加 `--run-prsformer` 和对应独立 GPU Python/device 参数。旧 PRSformer 输出需使用新版适配器重新发布原始协变量列，不能仅改列名混入比较。

新人预测要求沿用模型保存的评分材料、PCA、中心化和 feature contract，并提供 `--feature-manifest`。这是一份可核验身份的生成者声明，不是对外部计算正确性的独立证明。完整格式和示例见 README 第 7 节。

## 7. 本次验证与实际边界

测试只使用合成数据。统一运行包含 baseline、evolution、matching、evaluation/provenance 以及原有集成测试；其中真实 PLINK 2 在 600 人、300 个家系、2 条染色体、24 个 SNP 的 PGEN 上执行评分。LDL 特意缺失一部分训练结局以核对共同 roster 不被重新划分；合成 Disco scalar 用于测试下游校准连接，不冒充官方 R 插值输出。

**本次最终统一运行：57 项通过，0 失败、0 错误、0 跳过。**

| 测试组 | 通过数 |
|---|---:|
| baseline 输入/身份/指标契约 | 17 |
| 演化注释与负对照 | 10 |
| 全库计数与实际借用 | 9 |
| 评价、校准、来源与特征契约 | 11 |
| 三性状集成流程 | 10 |

另有 20 个 Python 文件及 8 个 Bash 文件通过语法检查，`git diff --check` 通过。生成的代表性 performance/contrasts PNG 已检查布局；这些图完全来自测试 fixture，不作为 UKB 性能展示。

重点验证包括：等位基因翻转、分数与权重哈希、年龄区间和未知分区、MAF 缺失的负对照、全库匹配计数的独立暴力计算对照、家系排除、基线/donor 贡献重构、二分类校准分离、无 Y 重载预测、错误特征声明、PRSformer 训练身份及实际协变量、三性状报告，以及改变 test Y 不会改变模型预测和筛选。

运行方法：

```bash
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 \
  python -m unittest discover -s f -p 'test_*.py' -v
```

PLINK 2 不在 PATH 时设置 `GRID_TEST_PLINK2=/实际路径/plink2`。测试环境为 Python 3.12、rdata 1.1.0 和 CPU 数值库；正式安装器所用 Python 3.11 与此不同，不能把本次测试理解成所有平台/依赖版本均验证。

**未执行：** 原生 R PCA/Yeval/Disco 程序、PRS-CSx MCMC、GPU/NATTEN PRSformer 训练、实际 HNSW 近似检索、完整 UKB 或外部队列分析。R 修改做了接口/静态契约检查；实际端到端集成走 Python、RDS 和 PLINK 路径。原 `VALIDATION.grid.md` 保留为历史记录，本次不继承其中 GPU 训练的“已运行”状态。

## 8. 下一轮真正决定方法价值的实验

1. 固定 discovery 来源并核查 UKB/test 重叠；给出每个祖源的 GWAS N、build/test N、病例数、SNP/年龄覆盖和实际标签预算。
2. 在开发集完成 CSx phi 选择，纳入 auto/meta，并把 phenotype-tuned Disco anchors 放到同一个冻结 roster；原 Yeval 分折结果不能替代这一步。
3. 先做完整 test 的配对 MSE/Brier，再做预先确定的 gain/low-error 分层；同时报告校准、覆盖与各祖源的结果。
4. 做 MAF+LD 匹配的年龄置换/敏感性分析、重复种子与外部人群验证，并加入适当的注释 PRS 基线（例如 SBayesRC）。
5. 若进一步研发“演化先验版 PRS-CSx”，需要把年龄/选择注释真正写入 SNP 效应先验，并重新估计 posterior；当前后验分区和 ABM 层不能宣称已经完成这项工作。

这五步中哪些成立、哪些不成立，应由预先定义的实验结果报告。新代码的贡献是让这些判断对应可追溯的输入、冻结的预测和可核对的个体证据。
