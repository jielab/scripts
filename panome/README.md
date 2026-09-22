# Panome 4 — 个体 reference matching 与分子组合解释

本版按照“建 reference → 找匹配 → COPY / 借用 reference 的结局信息”重写。以真实个体作为榜样，保留整体匹配、模块 mosaic、原始 COPY Y 和多参考风险估计，不以 AUC 决定唯一赢家。默认主面板预先固定为 100 人；300、1000 人只作为并列敏感性分析。

这是研究实现。真实 UKB 新版结果须在你的本地数据上重新运行；本包不包含这些结果。旧版结果的核查见 `AUDIT.md`。本版聚焦有首诊日期和随访的疾病，按指定时间窗预测 net risk；不提供旧版 height、Transformer、Leiden 或 PGS 接口。

本机替换前的修复、真实输入抽样核查及运行记录见 [LOCAL_REVIEW.md](LOCAL_REVIEW.md)。

## 运行

将整个 `panome/` 目录作为新版使用，避免混入旧 `f/`。默认输出到新的 `v4_reference` 子目录，不覆盖旧结果。

```bash
cd /mnt/d/scripts/panome
python3 -m venv .venv
.venv/bin/python -m pip install -r requirements.txt
chmod +x panome.sh
./panome.sh --Y cvd_cad --biom prot --dry-run
./panome.sh --Y cvd_cad --biom prot --preflight
./panome.sh --Y cvd_cad --biom prot --cores 16
```

不读取真实数据的快速检查：

```bash
./panome.sh --demo --tree hist --max-samples 1600 --run-name synthetic_check --analysis-root /tmp/panome_check
```

默认输入保持现有 UKB 布局：

- 表型 `/mnt/d/data/ukb/phe/Rdata/all.rds`。
- 蛋白 `/mnt/d/data/ukb/phe/rap/raw/prot.tab.gz`。
- CAD 首诊 `fod_icd10_cvd_cad`，基线 `date_attend`，死亡 `date_death`，失访 `date_lost`。
- 行政截止日期沿用已有分析的 `2023-04-01`，请按实际数据覆盖期设置 `--end-date`。
- 临床协变量 `age,sex,tdi,PC1,PC2,center`；蛋白残差化 `age,sex,prot.plate`。
- 死亡作为删失，输出是死亡删失下的 net risk，**不是竞争死亡下的实际累计发病概率**。

RDS 优先用已有 `Rscript` 提取字段，再由 `pyreadr` 读取二进制中间文件；`pyreadr` 已列入依赖，没有 R 时直接使用它读取。这样可保留批次等分类元数据中的极小数值。CSV、TSV、gzip 文本也可直接输入。不会 source LE8 或修改原始数据。`--tree lightgbm` 是默认树基准；未安装会报错，不会默默替换。可明确选 `--tree hist` 使用 sklearn 的 histogram gradient boosting。

亲缘连通分量列已准备好时，请加 `--group-col family_component`。这样外层划分、内层 OOF、匹配排除和 bootstrap 都以亲缘组处理。没有亲缘信息时是个人随机划分，不能声称已排除亲缘泄漏。

```bash
# 先完整训练并冻结，稍后才打开测试结果
./panome.sh --Y cvd_cad --biom prot --run-name cad_reference --train-only
./panome.sh evaluate --run-dir /mnt/d/analysis/panome/cvd_cad/prot/cad_reference

# 保留蛋白的年龄/性别相关部分，另做明确命名的敏感性分析
./panome.sh --Y cvd_cad --biom prot --residualize prot.plate --run-name cad_plate_only_residual

# 5 年风险是另一个明确时间窗，须重新拟合
./panome.sh --Y cvd_cad --biom prot --horizon 5 --run-name cad_5y

# 多次划分用于重训稳定性，不能挑选其中 AUC 最高的一次作为最终结果
for seed in 2026 2027 2028; do
  ./panome.sh --Y cvd_cad --biom prot --seed "$seed" --run-name "cad_seed_${seed}"
done
```

默认外层 50% test，另外 50% 分为 build/tune/calibration，占总样本约 30%/10%/10%。QC 后比例可能略变。所有模型都采用同样的划分。可通过 `--split-file` 提供含 `eid,split` 的 CSV，split 必须是 `build,tune,calibration,test` 之一。用固定 split 比较多个方案，避免样本变化造成混淆。

## 并列实现的方案

| 输出名 | 实际行为 |
|---|---|
| `copy1_fullproteome` | 最低 OOF log loss 的 100 个真实人；在全部保留蛋白的标准化空间寻找 1 人；直接 COPY 其已知 10 年 Y=0/1。原始方案的直接对照，不做概率校准。 |
| `copy1_topfit` | 同样的榜样，改用训练获得的分子距离匹配，再 COPY 1 人 Y。 |
| `copyk_fullproteome`, `copyk_topfit` | 相同榜样，K 人距离加权 COPY，加入 IPCW；保存原始值和独立校准值。 |
| `random_panel` | 按已知结局比例抽样的随机 100 人，检验是否真的需要挑榜样。 |
| `diversity_panel` | 优先覆盖分子空间的真实 100 人，不按 fit 挑选。 |
| `reliable_100` / `panome` | 重复 OOF 质量、结局类别配额和分子多样性共同选择真实个体；在其局部邻域估计风险，再进行加权匹配。 |
| `reliable_300`, `reliable_1000` | 面板大小敏感性分析；不会自动取代预设的 100 人主面板。 |
| `all_reference` | 使用全部已知时间窗结局的 build 人员，普通 IPCW 加权 kNN；检查压缩到 100 人损失了多少信息。 |
| `unsupervised_matching` | 与主方案相同参考身份，用不利用 Y 的 PCA 距离重新估计邻域和匹配。 |
| `mosaic_equal` | 不同分子模块分别匹配不同真实参考人；各模块风险等权组合。 |
| `mosaic_weighted` | 相同模块匹配，根据 tune 的 Brier 改善给予模块权重，保留 10% 等权收缩。 |
| `panome_clinical`, `mosaic_clinical` | 参考风险与临床预测在独立 calibration 中组合。 |
| `panome_elasticnet_hybrid` | 参考风险与 elastic net 组合；用于评价预测增量，不把全部预测归因于参考人。 |
| `clinical`, `protein_elasticnet`, `elasticnet`, `clinical_pca`, `lightgbm` | 临床、蛋白 EN、临床+蛋白 EN、临床+PCA、临床+蛋白树模型基准。 |
| `clinical_technical` | 临床+残差化使用的元数据+蛋白缺失比例；检查明显技术信息是否也能预测。不是完整的逐 assay 缺失模式负对照。 |

`--panel-size` 改变主面板大小；`--panel-sizes` 决定额外大小。主面板不够人数时会记录实际人数及 `not_ready`，不会偷偷补成“合格 100 人”。所有可以计算的探索性结果仍输出。

## 如何定义“榜样”

在 build 内重复 3 次、每次 5 折。每折重新拟合分子预处理和模型。默认质量教师为蛋白 elastic net 与浅层非线性树的平均预测，避免仅挑选“符合线性模型的人”。`--quality-teacher elasticnet` 可单独运行线性教师作为敏感性分析。

每个人记录 OOF 概率、OOF log loss、相对该折常数风险的 log-likelihood gain、重复之间的变化，以及正 gain 比例。主候选要求平均 gain>0 且至少 2/3 重复为正。未知时间窗结局者不能被当作 Y=0，也不能成为 COPY Y 的来源，但其分子数据可以参与 build 表征。

原始 `topfit` 对照按最低 log loss 全局取人；这很可能主要选出低风险非病例。改进方案在结局类别内挑选可靠且互相不同的人，配额接近 build 已知标签的比例。100 人本身仍不是疾病总体的无偏随机样本，因此独立校准以及未筛选邻域的风险估计很重要。

一个人只有一次观察结局。这里的“可靠”指给定模型下重复 OOF 的一致性，不是证明这个人的蛋白和疾病之间有确定关系，也不意味着其他人是坏数据。

## 个体解释与 copy 稳定性

| 文件 | 要回答的问题 |
|---|---|
| `reference_candidates.csv` | 每个候选人的 fit 依据、稳定性、两种教师的 OOF 预测是什么？ |
| `reference_panel.csv` | 最终参考人是谁？观察标签、局部风险、邻域有效人数和事件支持是多少？ |
| `test_reference_matches.csv` | 这个人整体匹配了谁？各参考人的权重是多少？ |
| `test_mosaic_matches.csv.gz` | 这个人的不同分子模块分别匹配了谁？各模块支持怎样的风险？ |
| `test_feature_explanations.csv.gz` | 哪些蛋白共同偏高/偏低，哪些蛋白与参考组合最不一致？ |
| `test_person_cards.jsonl.gz` | 每人一行完整 JSON 解释卡，可用于后续个体浏览界面。 |
| `molecular_profiles.csv`, `molecular_blocks.csv` | 每个人的模块坐标，以及模块包含哪些蛋白、如何计算。 |
| `same_risk_pairs.csv` | 预测风险相近、分子组合不同的个体，不要求不同离散 subtype。 |
| `masked_copy_stability.csv` | 遮住已测蛋白后，能否重构它们？参考身份和原始风险变化多大？ |
| `matching_diagnostics.csv.gz` | 各匹配方案的距离、有效参考人数、重构误差与局部事件支持。 |
| `reference_utilization.csv` | 有哪些参考人真正被使用？是否大部分目标都落在极少数人上？ |

默认分子模块从 build 蛋白的 PCA loading 相似性获得，并在每个模块中学习 PC1。它们是数据驱动的坐标，不会自动被命名为“炎症”“脂质”或被解释为通路活动。若已有审定的通路映射，可用 `--module-file pathways.tsv`，两列必须为 `feature,module`；特征名须精确匹配。允许一个蛋白属于多个外部模块，但这种相关性使模块风险不能作为可相加的病因占比。

mosaic 借鉴“不同部分匹配不同参考”的思想，没有蛋白 phasing、染色体顺序、重组率或 HMM。单个模块仅用 PC1 是可检验的首版假设，不能声称保留了该模块全部信息。

## 无合适匹配与面板评估

`reference_readiness.json` 分开记录结构/支持不足与预测性能提示。AUC 没有改善不会自动否定其解释价值。`provisional_pass` 只表示预设的内部支持检查通过，仍不是外部有效性的认证。

个体判断依赖冻结的 tune-X 阈值：参考距离、标准化蛋白重构误差、参考风险分歧、有效参考人数、未知分类水平及缺失比例。它们可在不知道该人的 Y 时计算。`supported_match` 和 `rejection_reason` 对所有 test 人员输出；面板和个体都通过时，`released_net_risk` 才有值。所有人的探索性原始/校准风险仍保留。这个字段仅用于研究流程，不是临床使用许可或个体正确概率保证。

`approach_comparison.csv` 同时给出预测指标和参考使用/匹配描述。`coverage_curve.csv` 在完全相同的被保留人员中比较所有模型，防止只报告容易预测的人而夸大优势。遮蔽重构默认抽取 1000 个 test 人、重复 3 次、每次遮住 10% 的已测蛋白；对 Panome、全蛋白 COPY 1、随机面板和多样性面板使用相同的人员与遮蔽位置，不使用其疾病标签。可用 `--mask-models` 选择其他匹配方案；可用 `--explanation-samples`、`--mask-repeats`、`--mask-fraction` 修改。

## 给新个体预测

只需要分子测量和基线元数据，不需要 Y。所有预测、校准、参考面板和阈值固定，目标批次不会成为彼此的 reference。

```bash
./panome.sh project \
  --run-dir /mnt/d/analysis/panome/cvd_cad/prot/v4_reference \
  --phe-file /path/new_baseline.csv --omics-file /path/new_proteins.csv \
  --output /path/new_people.csv
```

输出各方案预测、整体匹配支持，并另存 `new_people_mosaic.csv.gz`。保留 assay 必须齐全，个别测量可缺失。若训练使用亲缘组，新个体也要提供完整连通分量编号。跨平台测量单位或系统偏移未被此接口自动解决。joblib 仅加载你信任的本地训练产物。

## 执行、验证与边界

日志只记录主要步骤 START/DONE 和 OOF 重复结束。`--cores` 同时限制 BLAS/OpenMP 线程。保存冻结模型后，可独立 `evaluate` 重画图和计算指标；不提供折训练中断的自动续跑。已有结果目录默认拒绝覆盖；需要重跑用新 run-name，或对该 run 明确 `--replace`。

主要输出包括 `Fig_prediction_coverage.png/pdf`、`Fig_individual_profiles.png/pdf`、`REPORT.md`、`test_metrics.csv`、`paired_contrasts.csv`。Bootstrap 是对冻结模型的条件不确定性；跨 seed 重训与外部队列验证另行进行。当前小型合成数据验证见 `VALIDATION.md`。

本版疾病基准使用同时间窗的 IPCW logistic/树模型，与旧版 Cox 的 Harrell C 不同，不应把新旧不同划分的数值直接相减。IPCW 使用 build 的边际删失分布，依赖独立删失假设；病例抽样权重、协变量条件删失、竞争死亡及经验证的完整临床风险模型仍需按正式研究设计扩展。本版未执行动态数字孪生、干预模拟、MR 或因果机制发现。
