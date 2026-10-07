# PRSformer 适配器验证记录

后续整合说明（2026-10-07）：新 GRID 包为 `f/3.prsformer_data.py` 增加了评分结果的协变量、结局定义和准备签名字段；入口及训练模块保持原样。以下“三份代码逐字节一致”描述首次导入时的状态。本轮未重跑 PRSformer 神经网络训练，GRID 侧验证见 [VALIDATION.grid.md](VALIDATION.grid.md)。

## 当前目录落地验证（2026-10-07）

来源：`D:\Downloads\grid_prsformer_height_audit\grid`。`3.prsformer.sh`、`f/3.prsformer.py`、`f/3.prsformer_data.py` 按原包导入，三份代码逐字节一致。运行说明及 height 审计报告一并保留。既有入口及调用引用统一为 `0.pca.sh`、`1.csx.sh`、`2.disco.sh`，相关辅助文件为 `f/1.csx.py`、`f/2.disco.R`。

本次在本机实际执行：

- Shell/Python 语法检查、PRSformer ShellCheck，以及四个入口的 `--help`。
- PRSformer 五个模式的 dry-run，含空格路径、`--option=value`、冲突参数拒绝；确认 dry-run 不创建缓存。
- `f/1.csx.py` 与 `f/3grid.py` 的 `--help`，确认改名后的动态模块导入可用。
- 使用现有 `gu-threads` Python 环境，在临时目录生成 62 人、两个染色体、24 SNP 的 BED/BIM/FAM 和 TSV 表型。两个染色体使用不同样本顺序，变异倒序，包含缺失剂量、缺失表型、负 ID、撤回 ID、仅在非训练样本中变异的 SNP。
- 准备结果为 60 人、23 SNP；完整剂量矩阵与独立生成的真值逐元素一致。固定 train/validation/test 为 36/12/12，家系组没有跨集合；只检查输入时不创建缓存，再次 prepare 复用已有缓存。
- 修改 validation/test 的协变量和连续表型不改变训练集拟合的基线参数；保存的预处理参数重用后基线一致；含并列预测值的 AUC、average precision 与手算值一致。

临时合成数据已清理。本次未安装或修改现有 Python 环境。默认 `~/.venvs/grid-prsformer/bin/python` 与 `/mnt/f/software/PRSformer` 尚不存在，默认 `--check` 按预期报告缺失解释器。本次没有重跑真实 UKB 准备/训练，也没有重做神经网络、RDS 发布或 CUDA/NATTEN 验证；下方这些项目属于下载包原有验证记录。

## 下载包原有验证记录

日期：2026-10-07。以下记录是代码验证，不是 UKB 实证结果。

## 执行环境与范围

使用 PyTorch 2.6.0+cpu、pgenlib 0.95.1、pyreadr 0.5.7，实际加载官方 PRSformer 的 `g2p_transformer_ExplicitNaNDose2`。官方源码对应 `7dea4e1bb27975885c937f2be82f9243bc82bda7`。CPU 验证显式选择官方已有的 global-attention 小数据分支；没有把它当作 neighborhood attention 的性能测试。

未访问用户本地 `/mnt/d`、`/mnt/f` 的 UKB 基因型和表型，未执行真实 UKB 训练，未执行 NATTEN/CUDA 的全基因组训练。正式结果须以用户环境中的独立测试集输出为准。

## 端到端验证

构造 182 人、两个染色体、48 个候选 SNP 的合成数据，分别写成真实的 PGEN/PVAR/PSAM 和 BED/BIM/FAM 文件，并把表型写成 RDS。

故意在两个染色体中采用不同样本顺序，把变异顺序倒置，加入缺失剂量、训练集单态 SNP、仅在 held-out 样本中变异的 SNP、撤回 ID、负数 ID 和部分缺失的 height/LDL 标签。T2DM 使用 Yr2e/Yt2e 字段测试基线标签构建；固定 split 并提供 family group。

验证结果：

- PGEN 和 BED 两条读取路径都保留 180 人、46 SNP。
- 两个读取路径的最终矩阵逐元素相等，并与独立保留的合成基因型真值相等。
- 验证 ALT 剂量方向、染色体/BP 排序、跨染色体 IID 对齐、缺失值 `-1`、撤回/负 ID 排除，以及训练集单态 SNP 筛除。
- 训练、验证、测试分别为 108、36、36 人；所有 family group 均只属于一个 split。
- 三性状联合训练实际完成三个 epoch，使用两层官方 Transformer、gradient checkpointing、batch 7 和 accumulation 5，覆盖最后一个不足完整累积组的 batch。
- `3.prsformer.sh` 的 prepare → train → publish 全流程实际成功。
- 三个正式 RDS 均只含同一批 36 个 test ID，保留 outcome、baseline、prediction、prsformer 等列。
- 模型及 split RDS 成功发布；performance、training 的 PNG/XLSX 配对输出成功，工作簿可重新读取，PNG 已人工查看。
- 从正式发布的 `.pt` 经 wrapper predict 重新生成预测并再次 publish 成功。

## 训练模块独立验证

另外使用 144 人、48 SNP 的合成输入，检查训练控制和数值实现：

- 同一 checkpoint 重载后的 test 预测与首次输出一致。
- 在 epoch 1 后模拟中断，再恢复到原定 epoch 总数，最终预测与不中断运行一致。
- 修改 held-out 协变量和标签不改变训练集拟合的基线与标准化参数。
- AUC、average precision（含预测 ties）与 sklearn 的独立实现一致。
- 从测试输出独立重算 height/LDL 的 prediction R²，与 metrics 一致。
- `.pt` 内含完整 ordered SNP/REF/ALT、预处理、数据准备配置、训练历史和源码摘要。
- 临时 `evaluation.json` 的输入/输出摘要与模型文件信息均核对一致。
- `--covariates none --check` 成功；不生成模型输出目录。

## 防错检查

已实际验证以下输入会明确报错并停止：

- 把输出表型或其原始来源列作为协变量。
- SNP list 与目标基因型的 CHR:BP 不一致。
- 同一 family group 跨 train/validation/test。
- 改动已完成 run 的 metrics 后尝试发布。
- 两个 wrapper 同时操作同一个规范化 cache 路径，包括通过含 `..` 的路径别名访问。

另外验证：TSV 中空字段不会把后续列向前错位；prepare `--check` 不创建缓存目录；Bash 语法、带空格路径、`--option=value`、可选 TRUE/FALSE 参数和各阶段 dry-run 正常。

## 验证边界

上述检查覆盖数据对齐、剂量方向、训练/测试隔离、训练控制、保存重载及结果发布。它们不能代替真实 GPU 的 kernel 兼容性和内存验证，也不能证明 PRSformer 在 height、LDL、T2DM 上优于现有 CSx/Disco。完整比较必须使用固定、相同且与 discovery 不重叠的测试样本；当前交付未产生这一实证比较结果。
