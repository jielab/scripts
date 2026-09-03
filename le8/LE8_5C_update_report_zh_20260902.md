# LE8 / 5C 框架复核与 2026-09-02 代码更新说明

## 1. 复核范围

- 代码基线：`jielab/scripts` commit `329c1d17c0392feed4be36243c52b04ba4c75b3e`。
- 结果基线：`jielab/analysis` commit `f7bf3e1c887185f0e8bda55fb85f270b6a559c44`。
- 重点阅读：`le8/f/c1_correlate.R` 至 `c5_consolidate.R`、`c2_genetic_decompose.R`、`0f/0phe.f.R::t2e()`、`0data/gwas_post.sh` 的 COJO PGS 构建，以及 CAD 的现有 protein/metabolite 输出。
- 本版没有直接写入 GitHub；交付的是完整可替换的 `le8/` 代码目录。

## 2. 先给结论

1. “Yin–Yang”在概念上确实对应 baseline-prevalent 与 incident disease，但它不是标准流行病学术语，也不能把两类样本当成同一个风险集直接合并估计。更准确的论文表述是 **incident–prevalent–lead-time triangulation**；“Yin–Yang”可以保留为图形隐喻。
2. 你的核心质疑是成立的：一次成人血浆测量关联未来医院诊断，可能反映潜伏期疾病、诊断延迟、治疗、共病或器官储备，而非病因。现有 CAD 结果中 GDF15、NPPB/NT-proBNP、MMP12 的“实测值强、COJO PGS 弱、残差强”与 reactive/preclinical-consequence 解释相容。
3. 但“因此不应该用于 incidence prediction”过强。非因果 consequence biomarker 可以用于早期发现或短期风险分层；它不适合被解释为病因或一级预防靶点。必须把 **etiologic score、distal primary-prevention score、early-detection score** 分开报告。
4. PGS 是受精时已经确定的 inherited propensity，不是“出生时测到的 protein/metabolite level”。它不能把一个成人蛋白浓度反事实地外推到出生。
5. 现有分析仓库尚未完成真正的 5C：protein/metabolite 有 C1，protein 有较完整 C2，metabolite 只有部分 C2 原始表；没有可审计的 C3、C4、C5 结果。因此目前不能宣称完成 consolidation、网络重塑或优于临床模型的最终结论。

## 3. 现有 CAD 结果发现了什么

### 3.1 C1 蛋白质

- 蛋白 assay 样本 38,576 人，incident CAD 3,419 例，baseline-prevalent 1,756 例，2,910 个蛋白。
- adj2 complete-case incident 模型为 23,665 人、2,057 个事件。
- incident CAD：146 个 Bonferroni 显著，434 个 FDR<0.05。
- 代表性实测蛋白：NT-proBNP HR 1.326，`P=1.26e-41`；GDF15 HR 1.275，`P=1.07e-26`；MMP12 HR 1.275，`P=1.39e-26`；NPPB HR 1.243，`P=1.08e-20`。
- 方向性分类不是简单的 “causal/noncausal”：1,984 unresolved，466 established-disease associated，232 distal + disease-state mixed，128 distal antecedent，56 near-diagnosis/reactive-compatible，44 incident only。
- NT-proBNP、GDF15、MMP12、NPPB、PCSK9、LPA 在当前规则下均属于 mixed：它们在 5 年 landmark 后仍有关联，同时也与 prevalent disease 有关联。这能排除“纯粹只在诊断前几个月升高”的极端解释，但不能证明病因。

### 3.2 C1 代谢物

- 459,254 人，incident CAD 38,929 例，baseline-prevalent 18,439 例，300 个代谢特征；adj2 incident complete-case 为 283,437 人、23,767 事件。
- incident 分析 198 个 Bonferroni 显著、243 个 FDR<0.05。
- prevalent 结果中 LDL、non-HDL 等出现很强的反向关联，最合理的首要解释是既往诊断后的 statin/生活方式治疗与 survivor selection，而不是“低 LDL 导致既往 CAD”。这正是 prevalent 证据必须独立解释的例子。

### 3.3 C2 MR

- protein MR 共 957 行：460 cis、497 trans；FDR_analysis 显著 12 个 cis、6 个 trans。
- 当前 evidence grade 的 A 级 robust 候选为 TIMP4、LAYN、NOS3、IL1RN、HYAL1；B 级 heterogeneous 包括 GALNT2、LPA、FGF5、FURIN、PCSK9、CFB。
- 这说明强 observational predictor 与 causal/locus candidate 不是一回事；最终 panel 不应只按 ProtWAS P 值排序。

### 3.4 现有 individual PGS decomposition

当前 C2 只校准了 100 个匹配 score，而且是探索性的、非 cross-fitted 分解。几个锚点很有信息量：

| Feature | 实测值 incident | PGS incident | residual incident | 解释边界 |
|---|---:|---:|---:|---|
| GDF15 | beta 0.228, `P=9.07e-27` | beta -0.004, `P=0.855` | beta 0.266, `P=9.01e-28` | 强烈符合 acquired/reactive 成分，但不能单独证明 disease→GDF15 |
| MMP12 | beta 0.238, `P=8.47e-27` | beta -0.018, `P=0.422` | beta 0.264, `P=1.56e-29` | 同上；还需与 MR 方向不一致一起解释 |
| PCSK9 | beta 0.121, `P=1.86e-7` | beta 0.027, `P=0.231` | beta 0.119, `P=3.74e-7` | 弱 PGS 不能否定 PCSK9 因果性；score 强度、样本重叠与蛋白测量均会影响结果 |
| LPA | beta 0.096, `P=1.13e-5` | beta 0.042, `P=0.0569` | beta 0.087, `P=7.12e-5` | incident PGS 边缘；prevalent PGS `P=1.07e-5`，遗传成分较清楚 |

“residual”包含环境、药物、测量误差、未建模遗传成分和潜伏期疾病，不能命名为 “preclinical disease component”。现有代码已经正确写了这条警告，应保留。

### 3.5 DANDELION 与 MR-link-2

- 当前 CAD DANDELION 是 MAGMA gene-level sensitivity，不是基于合格 disease rare/WES gene evidence 的 primary DANDELION。
- 24/40 targets 被选择（60%）过宽，且所有结果的 `consolidation_eligible=FALSE`。这些结果可以画 sensitivity landscape，但不能计入主要 C2 或 C5 causal evidence。
- `c2.Fig5` 有 MR-link-2 结果，而 `c2.mrlink2_job_audit.csv` 为空，构成 provenance contract 不一致。新版会在原 job manifest 缺失但 completed results 存在时恢复最小审计表。

## 4. 没有发现什么、目前缺什么

### 4.1 结果层面的缺口

- 当前 GitHub analysis 中没有 C3 coloc/fine-mapping 结果，不能确认 MR hit 与蛋白 QTL/疾病 GWAS 是否由同一 causal signal 驱动。
- 没有 C4 LE8 connection/module/mediation 结果，更没有 disease-state network remodeling 结果。
- 没有 C5 held-out prediction/consolidation 结果，因而没有可复核的 C-index、time-dependent AUC、Brier、calibration、lead-time performance 或 compact panel。
- metabolite C2 不完整，protein 与 metabolite 不能做平行 5C 比较。
- 没有外部 cohort/platform replication；UKB 内部结果不能回答 assay transportability 与 ancestry generalizability。

### 4.2 设计层面的缺口

- `prot.pgs.rds`/`met.pgs.rds` 的权重来自 UKB-PPP GWAS，再回到 UKB 个体做 association/decomposition，存在 sample overlap 与 winner's curse。应优先使用外部 pQTL/mQTL 权重，或在 UKB 内 split/cross-fit 构建 score。
- 每个 PGS 必须报告 score→实测 omic 的增量 R²、F statistic/有效强度、cis-only 与 cis+trans 版本；trans score 不能自动解释为该蛋白的因果作用。
- full genetic cohort 与 omic assay subset 的样本构成不同。新版同时输出 full-genetic PGS scan 和 same-omic-sample PGS scan；Figure 1/2 用 same-omic-sample 版本，避免把 40 万人的 P 值与 3.9 万人的蛋白 P 值直接比较。
- 成人 baseline omic 与疾病诊断之间没有重复测量，不能把个体内临床前斜率与个体间差异分开。真正验证 consequence biomarker 最有力的是 serial samples、影像/亚临床表型以及诊断后的轨迹。
- hospital-record incident date 是“首次被记录的诊断”，不是生物学 onset。landmark 与 risk-set window 能降低但不能消除这个误差。
- prevalent cases 是入组时仍存活并能参加 UKB 的 survivors，受 Neyman prevalence–incidence bias、治疗与疾病严重度选择影响。不能把 prevalent OR 与 incident HR 当作可交换效应。
- 若目标是临床部署，还需要独立验证、校准、decision-curve/net benefit、成本/可检测性、跨平台稳定性以及预先规定的少量 panel。

## 5. Yin–Yang 提法与已有先例

### 5.1 是否等于 incidence–prevalence

概念上是，但建议在 Methods 中用标准术语：

- `baseline-prevalent disease`：诊断日在 baseline 之前；
- `incident disease`：baseline 时处于风险集、随访后首次诊断；
- `near-diagnosis incident` 与 `distal incident`：按真正的 risk-set/landmark 定义；
- `attained-age Cox with delayed entry`：成人采血后才进入风险集，不能称为从出生开始观察蛋白。

“Yin–Yang”只用于总览图和叙事，不用于变量名、模型名或因果结论。最稳妥的框架名称是 **Yin–Yang visualization of incident–prevalent–lead-time triangulation**。

### 5.2 有没有类似文章

有，而且与你的思路非常接近，但尚未等同于本版的个体级融合：

- Lind 等 2024 同时比较 2,919 个 **measured** protein 与 MR 所代表的 **genetically predicted** protein，发现大多数 observational protein–CVD 关联不支持因果；这是“实测 vs 遗传代理”的直接先例，但不是 40 万个体的 PGS–实测–prevalent–lead-time 联合分析。<https://www.nature.com/articles/s44161-024-00545-6>
- Henry 等在 incident HF 中明确讨论未识别/无症状 prevalent disease 和 subclinical disease 引起 reverse causation；NT-proBNP、GDF15 等 observational association 没有强 MR 支持。<https://pmc.ncbi.nlm.nih.gov/articles/PMC9010023/>
- 主动脉瓣狭窄的 prospective proteomics 研究发现，排除采血后 5 年内发病者后若干 protein 关联衰减，并指出 concurrent atherosclerosis 可能驱动结果，是 lead-time/lag sensitivity 的清楚先例。<https://pubmed.ncbi.nlm.nih.gov/29487139/>
- Pradeep/Schuermans 的四种 cardiac disease 研究与 Yu/Chen 的系统 CVD 研究属于现有的 observational association→enrichment→MR/coloc→prediction 路线。<https://www.nature.com/articles/s44161-024-00567-0>；<https://pmc.ncbi.nlm.nih.gov/articles/PMC12987571/>

本项目比较有辨识度的贡献应表述为：**把 baseline-prevalent、incident lead-time、个体级 omic PGS、实测 adult level、MR/coloc 及 held-out compact prediction 放进同一套可审计证据分层**，而不是声称首次注意到 reverse causation。

## 6. 两篇 2026 论文怎样处理

### 6.1 Nature：CIGMA

DOI `10.1038/s41586-026-10577-6` 对应的是 *Cell-type-specific eQTLs underlie the genetic architecture of complex traits*。CIGMA 是 population-scale scRNA-seq 的 variance-component 模型，用 donor、cell type、pseudobulk/cell-to-cell variance、cell proportion/count 与 kinship 分解 shared/specific eQTL；它不是新的 cis/trans MR。

- 仅有 UKB-PPP plasma protein、GWAS summary statistics 与 PGS 不能运行 CIGMA。
- 不应把它作为新的 C2 causation test，也不能因为 cell-type-specific eQTL 就增加一个 causal evidence count。
- 合理用法是：在外部 scRNA 数据跑完 CIGMA 后，把 cell-type specificity 作为 cis instrument/coloc locus 的生物学注释，帮助解释 trans 或 tissue-of-action。
- 新版 C2 增加 `c2.method_scope_audit.csv`，明确记录 CIGMA 当前不具备输入资格。

论文：<https://pubmed.ncbi.nlm.nih.gov/42649298/>；方法代码：<https://github.com/Minhui-Chen/CIGMA>。

### 6.2 Science：autism PPI rewiring

该研究用 AP-MS 对 100 个 ASD risk proteins 及 54 个患者 missense variants 进行物理互作测量，并用结构预测、Xenopus 与 forebrain organoid 验证 mutation-induced differential PPI。<https://www.science.org/doi/10.1126/science.ady4523>

UKB-PPP 的 plasma abundance correlation 不是 PPI，更不是 mutation-induced rewiring。新版借用的是分析结构，而不是证据名称：

1. distal disease-free participants 作为 reference state；
2. near-diagnosis incident 与 baseline-prevalent 作为 disease states；
3. basic-covariate residualized omic correlations；
4. Fisher-z differential-correlation test；
5. FDR + effect-size 双阈值的 remodeled edges；
6. convergent hubs；
7. 可用 `C4_PPI_FILE` 注释已知 PPI backbone，但图注仍明确写“co-abundance remodeling, not physical PPI rewiring”。

## 7. 新版代码改动

### C1 correlation

- 新增 `f/c1_pgs.R`。
- `c1.Fig1`/`c1.Fig2`：左列改成 same-omic-sample 的 inherited protein/metabolite PGS；右列保留 measured adult omic adj2。
- 原 basic observed models 没有删除，继续保留在 CSV/RDS/Excel。
- 同时输出 full-genetic PGS 与 same-omic-sample PGS，避免样本量混淆。
- 新增 `c1.Fig14.pgs_actual_concordance.png` 和 Excel worksheet `pgs_actual_concordance`。
- 保留 `as.character(eid)`，并在 PGS/phenotype 两侧都显式转换。
- observed scan cache 版本保持不变，只新增 PGS cache，避免无意义地重跑 ProtWAS/MWAS。

### C2 causation

- 新增 method-scope audit，明确 CIGMA 不是当前 C2 causal method。
- 保留 cis/trans MR、reverse MR、MR-link-2 与 DANDELION。
- 修复 completed MR-link-2 results 与空 job audit 的审计不一致。
- estimator stage cache 版本与 final-output 版本分离，避免只改审计输出就重算 MR/DANDELION。

### C3 colocalization

- 新增 `f/c3_pgs.R` 和 `c3.Fig7.pgs_coloc_triangulation.png`。
- 生成 PGS、实测值、residual 与 robust PP(H4) 的 feature-level triangulation table。
- PGS 只作正交的 feature-level support，绝不进入 `coloc.abf()`，避免把全基因组 score 冒充 locus summary statistic。

### C4 connection

- 新增 `f/c4_state_network.R`。
- 保留原 Figure 1–9，新增 Figure 10–11 与同一 `c4.out.xlsx` 中的多 worksheet。
- 网络以 distal controls 为 reference，比较 near-diagnosis incident 与 prevalent state。
- 可选 `C4_PPI_FILE` 只做 known-PPI annotation。

### C5 consolidation

- 保留所有原有 prediction paradigms。
- 新增 `Mechanism-weighted compact`：候选只来自 outer-training screen、training-only 5-year distal screen 和 MR+coloc genetic set；相对 glmnet penalty factor 为 genetic+distal 0.40、genetic 0.55、distal 0.70、training-ranked only 1.50。
- 这比按 GDF15/NT-proBNP 名称手工降权更可复现，也避免使用 validation outcome 选择 feature。
- evidence consolidation 另给数据驱动的 `evidence_weight_prior`；reactive-compatible 且没有 locus evidence 的 feature 权重最低，但该列是可审计 prior，不冒充 causal effect。

## 8. 建议的运行顺序

先确认自动发现的文件存在：

```bash
ls -lh <UKB_PHE>/Rdata/prot.pgs.rds <UKB_PHE>/Rdata/met.pgs.rds
```

预检与全流程：

```bash
cd /mnt/d/scripts/le8
./le8.sh -Y cvd_cad --biom prot,met --preflight

export C1_PGS_MAX=0
export C4_NETWORK_TOP=80
export C4_NETWORK_DELTA=0.20
# 可选：两列 feature1/feature2 的物理 PPI reference
# export C4_PPI_FILE=/path/known_ppi.tsv

./le8.sh -Y cvd_cad --biom prot,met \
  --steps c1_correlate,c2_cause,c3_coloc,c4_connect,c5_consolidate \
  --run-mrlink2 Top --run-dandelion Top --nested-cv TRUE
```

若没有合格的 CAD WES/rare-burden gene P values，DANDELION 只能按 sensitivity 解释：

```bash
./le8.sh -Y cvd_cad --biom prot --steps c2_cause \
  --run-dandelion Top --dandelion-allow-magma TRUE
```

有合格输入时才做 primary DANDELION：

```bash
./le8.sh -Y cvd_cad --biom prot --steps c2_cause \
  --dandelion-gene-p /path/cvd_cad.WES_burden.tsv \
  --dandelion-gene-annotation /path/genes.GRCh37.tsv
```

## 9. 运行后必须核对的 headline 表

1. `c1.pgs_status.csv`：PGS 文件、匹配 feature 数、full 与 same-omic N。
2. `c1.pgs_actual_concordance.csv`：不能只看 P 值，要看方向与 same-sample effect。
3. `c2.individual_genetic_decomposition.csv`：R²、cross-fitting 状态、observed/PGS/residual。
4. `c2.method_scope_audit.csv` 与 `c2.mrlink2_job_audit.csv`。
5. `c3.pgs_observed_coloc_triangulation.csv`：PGS support 不等于 PP(H4)。
6. `c4.state_network_status.csv`、counts、edges、hubs；不得把无 PPI reference 的 edge 称作 physical interaction。
7. `c5.mechanism_weight_prior.csv`、held-out performance、calibration、lead-time curves 与 score coverage audit。

## 10. 验证状态

- 已完成 Git diff whitespace 检查和所有 R 文件的引号/括号词法平衡检查。
- 当前交付容器没有 `R`/`Rscript`，也没有 UKB phenotype/QTL/GWAS 私有数据，因此没有伪造运行结果；仍需在 `/mnt/d` 的 `le8` conda/R 环境执行 `--preflight` 与真实数据 rerun。
- 新版保留原 Figure 1–13/1–14、原 raw outputs 和 Excel 结构，只增加 panels/worksheets；没有删除原分析逻辑。
