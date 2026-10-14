# GRID 新版：genetic-profile twin、事先 confidence、Final 与 Shiny

本包用于交给 Codex 核查和安装到现有 `jielab/scripts/grid`。**这是完整的新 ABM 数值实现、调用入口、结果索引和 Shiny 源码，不是只包含修改意见的文档。** 为避免覆盖已经反复检查的上游，不重复打包 PCA、PRS-CSx、DiscoDivas、PRSformer 的原实现；安装器保留这些文件，只新增独立模块并替换公共分发入口和 Shiny。

核对的上游基线：`jielab/scripts` commit `c2872afde581bec243f118717496a164a833a69c`。该提交的 GRID 文件与前一提交 `086be9413a93…` 相同。本包没有推送到 GitHub，也没有更改已有分析结果。

本目录已于 **2026-10-14** 完成本机整合。旧入口完整保存在 `grid.upstream.sh`，被替换文件备份在 `_twin_install_backups/20261014T132630896588/`。上游受保护文件的 SHA256 验证均通过，记录在 `verification/upstream-integrity.json`。

整合修复了旧入口日志跳转、子进程退出码、Python/R 环境选择、冻结人员及家系名单校验、PRSformer 缺失来源标记校验、损坏 manifest 的导入回退，以及 Shiny 的同名函数遮蔽、多子集符号和大表分面绘图问题。索引会在遍历前跳过缓存、工作目录和旧 twin run，并避免对挂载盘重复解析路径。数据库在临时文件内构建，完成后统一同步到磁盘并原子替换；失败保留原索引。

## 1. 这次解决什么

主路径为：

**冻结的原始 train/test 名单 → build 内重复 OOF 评价参考者 → Full/Clean/Soft/Random 训练对照 → 多维 genetic profile 匹配 → 每个 query 的事先 reference-support confidence → 单独读取 test Y 作评价。**

训练参考者的筛选默认采用 `profile_concordance`。二分类用 OOF log loss；连续性状用按 teacher 训练方差标准化的 OOF 平方误差。误差较大者可被排除，但结果不会把其命名为已证实的错误或“噪声”。默认只排除 1%（`keep_fraction=0.99`）；这是预声明、保守的执行起点，不是观察 test 后找到的最优值。可以在看结果前修改配置做预声明的 95% 等保留比例；若病例被删到不足，程序报错，绝不偷偷补回病例或重新分 test。

- **Full**：所有技术 QC 合格的 build 人员。
- **Clean**：通过技术 QC 且 OOF profile 一致性较好的参考者。
- **Soft**：同一个 OOF 合格集合按一致性连续降权。
- **Random**：随机保留与 Clean 相同规模、相同病例×OOF 风险分层构成的对象。连续结局匹配 OOF 预测分层。

技术排除与 profile 不一致排除分开记录。全遗传特征缺失的人不能通过临床变量稀释缺失率，不能成为 donor，也不会获得正 confidence。所有 testing 人员保留在原始人数和 coverage 的分母中。

## 2. Twin 指的是什么

`twin` 是特征相似的参考者，不是必须找到真实双胞胎。默认沿用现有 GRID prepared cache 的 `csx / ancestry / frequency / evolution` 多维遗传摘要。临床协变量不进入遗传缺失率或匹配距离的分母。分子版的“10 个蛋白”预算**没有**机械复制到 GRID；GRID 默认保留当前传入的遗传特征。

距离是在标准化后的多维遗传向量中计算，每个保留维度都有正距离权重；不是只按 `beta'X` 或一个 PRS 总分匹配。高于 32 维的输入可用 build-only PCA，并保留一个重构残差坐标。压缩后的相似性不能证明原始基因型完全相似。

**重要范围：若只输入多个 PRS、PC 或遗传模块分数，本程序不能还原或检查组成它们的所有 SNP 组合。** 输出如实标记 `profile_kind=genetic_summary`。只有明确提供 dosage 向量并设置 `variant_dosage` 时，才执行 0–2 范围的剂量质控；它仍不是全基因组的穷举单倍型检验。推荐复用现有低维遗传特征，不要直接传入未经设计的百万列 dense genotype 矩阵。

参考半径只用开发阶段 tune_model 的 X 冻结；Full/Clean 使用相同半径，不因删人而扩大匹配范围。半径内没有足够参考者时，不强行补齐 k 个邻居。

## 3. 每个 testing sample 的事先输出

本地 `_LOCAL_ONLY/confidence.tsv.gz` 包含：

- `pre_prediction_confidence`：0–100 参考支持分数。
- `genetic_missing_fraction`、`covariate_missing_fraction`。
- `eligible_match_count`、`matched_count`、`effective_reference_groups`。
- `nearest_radius_ratio`、`radius_ratio`、`donor_oof_loss`。
- `selected`、`reason` 和各模型预测；`prediction::ABM_policy` 为选择＋回退结果。

confidence 是距离、有效参考支持和 donor OOF 一致性的组合，**不是“预测正确的概率”**。其与实际误差的关系在 audit/test 中另外评价。

无足够 twin 但技术合格：回退到 `CSx_reference`（没有 CSx 输入时回退 `ABM_full`）。技术 QC 不合格：策略预测为 NA，confidence 为 0；不能把依靠大量填补的 fallback 包装成可靠预测。

没有家系资料时按个人 ID 分组，并标记 `individual_unverified`。存在明确家系分组则排除同家系参考和跨组泄漏，但“有一列 family_id”不自动证明相关个体已全部识别。**不需要冠脉检查表，不新增这类运行前提。**

## 4. 怎样区分真正的收益

令 T 为原测试中有结局且技术合格的人群，S 为冻结规则事先选出的人群，B 为相同特征／相同架构的 `ABM_full`，A 为过滤训练的 `ABM_clean`：

```
selection_effect           = Performance_B(S) - Performance_B(T)
training_effect_same_subset= Performance_A(S) - Performance_B(S)
combined_effect            = Performance_A(S) - Performance_B(T)
```

CSx、Disco、PRSformer 等有可核实逐人预测的方法也在同一个 S 中比较。不能将 `A(S)-B(T)` 全部归因于模型改善。

连续性状主误差是 MSE，并报告 RMSE、MAE、`prediction_R2 = 1-SSE/SST` 和以冻结训练方差为分母的 MSE。**prediction_R2 允许为负，不等于相关系数平方；不同子集的表型方差不同，R² 升高不自动意味着误差下降。**

二分类报告 Brier、LogLoss、AUC、AUPRC；同时显示人数和病例 coverage。不能只保留易预测对照而不展示漏掉多少病例。

区间通过同一次原测试群组重抽样计算，保持 S 嵌套于 T。bootstrap 默认 500；少于 20 个有效抽样不提供区间，少于 200 个有效抽样不触发“支持”状态。区间只条件于已拟合模型和冻结规则，未计入完整训练不确定性。主 coverage 默认 50%，其他 coverage 为预声明探索，不能事后挑最好的一档。

额外对照包括 baseline-confidence、低风险、随机、同风险×祖源构成随机子集和 gate 预测误差。主 selector 仍是 `support_quality`，不是“attention 相对 CSx 的 gain”。

## 5. 原方法与新参考模型的名称不能混淆

| 名称 | 本包实际含义 |
|---|---|
| `CSx_reference` | 使用现有分来源 CSx 分数，在本实验 build/tune/calibration 角色中拟合的 pooled 参考模型 |
| `Disco_reference` | 对传入 Disco 分数作本实验内的训练限定校准；不是重新运行 DiscoDivas |
| 历史 `PRS-CSx` / `DiscoDivas-untuned` | 原 final 的已保存结果，数值与名称原样展示，标记 legacy |
| `PRSformer` | 严格核实原 cohort、test 名单、结局单位及协变量后的原保存预测 |
| 外部 benchmark manifest 中的方法 | 已冻结的真实结局尺度预测；文件和人员身份核实，训练来源仍需提供者核查 |

**`CSx_reference` 不冒充原 final 中按祖源拟合的 PRS-CSx；`Disco_reference` 不冒充 untuned 结果。** 如果缺少原方法逐人预测，只能显示其已有整体汇总，不能从原来的 AUC/R² 反推出新 S 上的结果。

原 PRSformer 的默认自动读取位置是：

```
<analysis-root>/prsformer/split.rds
<analysis-root>/prsformer/scores/<trait>/predictions.rds
```

不满足原测试、结局、协变量的严格一致性时，该比较标记 unavailable，不缩小 ABM 的测试人群来迁就它。

其他已经保存的官方 benchmark 预测可以用 `--comparators benchmarks.json` 接入。`examples/benchmark_manifest.example.json` 提供格式；`f/build_benchmark_manifest.py` 可自动计算文件及 test 名单 hash。只允许**结局尺度预测**（height 的实际单位，T2DM 的 0–1 概率），不能将原始 PRS 值直接当预测值计算 MSE/Brier。

PGS 原始 GWAS 与目标 UKB 的样本重叠不能被下游 OOF 自动修复。本包将该状态默认标为 unknown；这是解释边界，不是索要新数据的运行门槛。

## 6. 安装：保留原上游

在解压后的 `grid_new` 目录：

```bash
conda activate grid
python install.py --grid-dir /mnt/d/scripts/grid --check
python install.py --grid-dir /mnt/d/scripts/grid
```

安装器：

1. 将原 `grid.sh` 完整保留为 `grid.upstream.sh`。
2. 保留原 `0.pca.sh`、`1.csx.sh`、`2.disco.sh`、`f/main.py`、`f/data.py`、`f/abm.py`、PRSformer 等文件，检查其 SHA256 没有变化。
3. 加装新 `f/twin_*.py`、浏览器索引及测试，替换公共入口和 Shiny；被替换文件备份在 `_twin_install_backups/`。
4. 不接触 `/mnt/d/analysis/grid` 的任何旧结果。

新代码**不依赖 LE8 目录中的 Python 模块**。新 Python 文件可以独立交给 Codex 检查。原模型仍可使用 `./grid.sh legacy-abm`，原 final 可使用 `./grid.sh legacy-final`；`./grid.sh prsformer` 继续委托原入口。

本机 `grid.sh` 默认使用现有 `$HOME/miniforge3/envs/grid` 的 Python，可通过 `GRID_CONDA_ENV`、`GRID_PYTHON` 或 `ABM_PYTHON` 指定其他环境。Shiny 自动选择依赖完整的 R，支持 `GRID_RSCRIPT` 或 `--r-bin`。本次为 GRID Python 补装了 pytest，为系统 R 的用户库补装了 RSQLite；依赖声明也已同步更新。分析命令不会自动安装软件。

## 7. 直接复用你的 prepared cache

```bash
cd /mnt/d/scripts/grid
./grid.sh abm --check --traits height,ldl,t2dm

./grid.sh abm --traits height,ldl,t2dm
```

该命令默认使用 **4 个线程**，自动读取 `grid.sh` 同目录下的 `config/twin.json`，复用 `/tmp/grid-cache/grid/prepared.joblib`，将新结果写入 `/mnt/d/analysis/grid/abm`，并从 `/mnt/d/analysis/grid` 读取已保存的比较方法。无需再传这些参数；显式提供参数仍可覆盖默认值。配置路径按安装目录定位，从其他工作目录调用也有效。bootstrap 默认次数取自配置文件（当前为 500）。

不更改原 outer test。若原缓存位于别处，改 `--prepared`；缓存完全不存在时，应先明确运行保留的原准备步骤，不能凭空生成一份新测试名单。也可用 `--data-file`＋`--schema` 读取明确的数值特征表，见 examples。

每次拟合新建不可覆盖的 run：

```
analysis/grid/abm/twin/<trait>/<UTC-run-id>/
  protocol.json
  fit/aggregate/       # role、OOF、训练筛选、遗传特征参数
  audit/               # 独立开发审计
  test/                # 全体/入选/拒绝、三类效应、各benchmark、confidence
  _LOCAL_ONLY/         # 个人预测、真实结局、模型及临床校准
  COMPLETE.json
```

全部成功后才更新 `CURRENT.json`。失败 run 不进入新 Shiny 的默认当前结果；旧 run 保留。三种性状独立执行，一个失败不会删除前面已经完成的 run。

`run` 会保存原 dev 和 test_X/test_truth 文件；以后的单独运行可用：

```bash
./grid.sh fit --data DEV.npz --config CONFIG.json --output NEW_FIT
./grid.sh project --model NEW_FIT/_LOCAL_ONLY/model.joblib \
  --data TEST_X.npz --output LOCAL_PREDICTION.joblib
./grid.sh evaluate --model NEW_FIT/_LOCAL_ONLY/model.joblib \
  --predictions LOCAL_PREDICTION.joblib --truth TEST_TRUTH.tsv.gz \
  --endpoint height --output NEW_REPORT
```

单独 project 严格拒绝含 y/time/event/outcome 等额外字段的 NPZ。既有 test 已被多次查看，结果始终标为 exploratory；重新分 split 不能让已参与设计的数据变成未触碰验证。

## 8. Final 和 Shiny

```bash
./grid.sh final --analysis-root /mnt/d/analysis/grid
./grid.sh shiny --analysis-root /mnt/d/analysis/grid --port 3839
```

- `final` 只读取现有结果：历史 benchmark `.xlsx` 和本包新 aggregate `.tsv.gz`，不拟合模型。
- SQLite：`analysis/grid/evidence/evidence.sqlite`。
- 历史 dashboard 快照不重复导入；新 twin 只默认显示各性状 `CURRENT.json` 指向的完成 run。
- 不加载任何个人 RDS、pickle 或 joblib；个人 confidence 不会被自动公开。
- 表格、图、TSV/PNG 下载使用同一个筛选后的全部分页行；零行不回退全表，切换表格后旧行号失效。
- 页面以 trait/module/table/ancestry/subset/method/metric 下拉选择为主，按选择画一张图；没有数千张静态图循环。
- 原 `heldout_residual_correlation_R2` 保留其原名称，不包装成新预测 R²。
- Shiny 默认只绑定 localhost。启动后是前台进程，Ctrl+C 停止。

本机已执行 R `testServer` 验收和真实历史结果的 `--review`，并使用本地 Chrome 实际检查 trait/源表切换、多方法选择、列过滤、零行、分页和 TSV/PNG 下载。浏览器记录见 `verification/browser-acceptance.json`。这是本地 Shiny 应用，未部署到公网。

## 9. 测试和核查

```bash
python -m pytest -q tests -p no:cacheprovider
Rscript --vanilla tests/test_shiny.R
./grid.sh shiny --analysis-root /mnt/d/analysis/grid --review
```

`--review` 检查 DB 可读及有界绘图构建，不替代浏览器交互测试。

合成演示（不代表 UKB 效果）：

```bash
./grid.sh demo --output /tmp/grid_twin_height_smoke
./grid.sh demo --binary --output /tmp/grid_twin_t2dm_smoke
./grid.sh final --analysis-root /tmp/grid_twin_height_smoke
```

自动测试包括：全遗传缺失拒绝／donor 排除、临床列数量不稀释 QC、家系排除、冻结 test 保留、query 结局隔离、单人／批次不变、模型序列化、audit 标签不改变拟合、随机训练对照、三类效应、错误 R² 防护、比较人群检查、索引隐私／checksum／历史隔离和安装保护。

本机 Python 回归测试为 **58 passed、1 skipped**；跳过项是只在原分发目录执行的安装器验收，该项另行执行且通过。连续及二分类 demo、独立 fit/project/evaluate、R 测试均已通过。真实缓存的三个性状均保留 **243,372** 名原测试人员；现存 PRSformer 的三个性状保存预测均通过完整名单、结局与协变量核查。官方 CSx/Disco 逐人 final 预测文件未保存，因此其新子集比较仍需显式提供冻结预测，历史汇总照常显示。

精确执行记录见 `verification/integration.json` 与同目录日志。**本次未重新训练真实 UKB，未声称新版已胜过 CSx/Disco/PRSformer。** 对代码正确性的验收，不要求制造一个阳性预测结果。
