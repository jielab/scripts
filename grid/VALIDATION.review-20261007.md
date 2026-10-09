# GRID 审阅版本机核查与更新 · 2026-10-07

## 来源与更新范围

已读取用户指定的 [ChatGPT 对话](https://chatgpt.com/c/6ac62614-6400-83ec-bfd2-5f7c6b6bd6e0)，并核查下载目录 `/mnt/d/Downloads/GRID/GRID-reviewed-20261007`。

- `SHA256SUMS.txt` 所列文件全部通过校验。
- 当前目录不是 Git checkout。将交付包复制到临时目录并逆向应用 `GRID-20261007.patch`，重建审阅基线；本机更新前的 38 个文件与该基线逐字节一致。补丁 dry-run 也全部通过。
- 完整应用交付包的 24 个新增或修改文件，应用后 44 个文件与包内 `scripts/grid` 逐字节一致，再实施下述本机补充修复。
- 包内共享 `0f/memory_cap.sh` 与本机一致；`results.py`、`results.R` 仅末尾空行不同。保留本机共享 helper，未覆盖其他项目文件。
- 保留原有 `VALIDATION.grid.md`、`VALIDATION.prsformer.md` 作为历史记录。本文件记录此次本机实际验证；交付包原有的 57 项记录见 `REVIEW.20261007.md`。

完整更新前备份：`/mnt/d/Downloads/GRID/grid-before-review-20261007-200816/grid`。同级 `before-sha256.json` 保存更新前文件哈希，`validation/` 保存本轮测试、预检、运行环境和入口检查日志；`local-review.patch` 保存本机相对交付包的补充改动，`after-sha256.json` 保存最终目录哈希。

## 核查发现与补充修复

1. **旧 CSx 权重迁移会被预处理代码哈希阻断。** 旧 `preparation_signature` 包含整份 `0.common.py` 的内容哈希。交付包虽然要求仅重新评分补齐来源，`native_weights` 却要求旧权重元数据匹配新版代码哈希，两者不兼容。本机的 12 组既有权重也存在这一情况。

   现在保留旧预处理元数据，通过原始联合推断配置核对 GWAS/SNP 参考的路径、大小和修改时间，并核对性状、祖源、原始输入路径；随后仍强制验证重新评分的分数哈希、权重哈希、联合推断签名、染色体及评分规则。通过的旧权重标记为 `historical_preparation_inference_inputs_verified`。输入改变、无法识别的旧签名或缺失/错配的评分来源均被拒绝。此兼容检查不表示旧预处理逻辑已重新认证，也不免除发现上游错误后重新推断的要求。

2. **无协变量的 PRSformer 比较误把 `none` 当成列名。** GRID 数据准备与 PRSformer 均支持 `none`；新版比较层原本会要求 `covariate.none`。现在使用统一协变量解析，并仅在实际存在协变量时核对数值。

3. **缺少年龄注释时预检仍先加载巨大表型文件。** 首次本机预检在读取表型时已消耗数分钟，而缺少默认 GEVA 目录可以立即判定。已停止该次只读预检，并将相同检查提前到表型加载之前。修复后的真实入口立即返回明确缺失注释错误，无需读取整份表型。新增回归测试确保这一分支不加载数据模块。

## 本机实际验证

交付包先在本机环境完成 **57 项测试，全部通过**。补充修复后，在最终项目目录统一运行 **65 项测试，0 失败、0 错误、0 跳过**，耗时 60.875 秒。

| 测试组 | 通过数 |
|---|---:|
| baseline 输入、等位基因频率、身份与指标 | 17 |
| 演化注释与负对照 | 10 |
| 三性状集成及无表型预测入口 | 11 |
| 精确匹配人数、实际借用及家系排除 | 9 |
| 旧权重迁移与输入变更拒绝 | 4 |
| 原生 R PCA QC、Yeval 入口保护 | 2 |
| 评价、校准、来源、特征契约与提前预检 | 12 |
| **合计** | **65** |

集成测试使用 600 人、300 个家系、两条染色体和 24 个 SNP 的合成数据，通过真实 PLINK 2 执行 PGEN 评分，再运行 prepare → evolution → fit → report。新增测试进一步经过保存的模型、特征契约和精确输入哈希调用 predict，在没有 Y 的输入上重现三个性状的保存预测。

原生 R 验证执行 `0.pca.R projection` 的投影后 QC，使用独立 `--score-dir`，核对工作簿中全部 22 条染色体的实际 SNP 行数；另直接执行 `Yeval.R`，确认监督训练分数在读取输入及创建结果目录之前被拒绝。R openxlsx 写入的工作表 dimension 提示可能是 `A1`，测试通过重置流式读取维度核对实际完整单元格。

另外通过：22 个 Python 文件语法、8 个 Bash 文件语法、3 个 R 文件解析、ShellCheck error 级检查、`pip check`，以及 10 项帮助入口、默认/显式 Python 环境与含空格路径的编排 dry-run 检查。

使用既有 `/home/huangj/miniforge3/envs/grid` 环境，未更新或重装依赖：Python 3.11.6、NumPy 1.26.4、pandas 2.3.3、SciPy 1.17.1、scikit-learn 1.9.0、rdata 1.1.0、joblib 1.6.0、PLINK 2 a6.9、R 4.2.3。完整版本见备份目录 `validation/runtime.json`。

复现统一测试时使用共享环境入口，以免 R 继承与当前解释器不兼容的个人包库：

```bash
cd /mnt/d/scripts/grid
source f/0.common.sh
grid_activate_environment
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 \
  python -m unittest discover -s f -p 'test_*.py' -v
```

## 真实数据就绪状态与后续运行

只读核对确认 height、ldl、t2dm 共 12 组祖源权重的原始推断输入身份均可验证。但三个性状目前均没有 `1csx.scores.provenance.json`，默认 `/mnt/f/ref/GEVA` 也不存在。最终真实入口 `bash grid.sh --traits height,ldl,t2dm --check` 返回 exit 2，提示缺失 GEVA/canonical 年龄注释；它不是完整数据预检通过。

正式分析之前，按 README 使用相同原始推断参数执行 `1.csx.sh --stage score`，为旧分数补齐来源；使用 Disco 时同步重新生成相应分数。另需提供有来源的年龄注释，存在真实家系连通分量表时传入 `--split-group-file`，并为新版 GRID 使用独立缓存与结果目录。

本次没有重新评分或修改现有真实 CSx/Disco 结果，没有运行真实 UKB 拟合、PRS-CSx MCMC、完整 R PCA/Yeval/Disco 分析、GPU PRSformer 训练或实际 HNSW 检索。代码及所覆盖契约已经验证，不构成真实队列预测优势的证据。
