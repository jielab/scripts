# GRID 默认无需外部年龄数据库：验证记录（2026-10-09）

## 修改范围

删除了原数据库的默认目录、命令行参数、输入预检、染色体读取器、源优先级合并、下载入口、URL 和下载校验代码。默认 `./grid.sh --traits height,ldl,t2dm` 用已有 CSx 四组分数、协变量和祖源投影 PCs 训练 CUDA selective-attention ABM，主模型名为 `GRID_no_evolution`。

默认 `evolution` 阶段只核对 CSx 评分与权重的来源，记录 `reference_matching`，不重算注释模块分数，不创建年龄或置换特征，不报告年龄模型或重复的年龄增益比较。gate 使用当前实际存在的 Ridge/HGB 对照。

显式 `--annotation-file` 仍可读取其他符合既有规范的注释表，相关来源、等位基因方向、质量和区间约束保留。这是一个单独的可选实验，不参与默认运行。模型版本为 `1.3-optional-annotations`；旧模型和旧缓存仍受源码及特征契约校验，不能混用。

没有修改 PRSformer 脚本或它的运行进程。修改前文件备份位于 `/tmp/grid-remove-geva-before-20261009`。

## 55 分钟测试窗口

北京时间 09:23:45 开始，10:19:03 结束，实际持续 **55 分 18 秒**，达到用户要求的 55 分钟。用户的一小时延迟任务开始于 09:22:28，预计 10:22:28 启动；测试结束时它仍在等待。所有测试使用独立缓存和输出，默认正式缓存没有被测试占用。

最终状态：**96 项自动检查通过；8,000 人真实子集三性状端到端通过；全量 height 训练、预测和独立重载校验通过。** 全量长测试未发现异常、非有限训练损失或显存不足。测试窗口结束时仅对确认属于临时目录的测试 Python 进程发送 SIGTERM；包装脚本退出码 143 是主动限时结束的结果，不是程序报错。已确认该进程退出、GPU 资源释放、默认缓存锁可获取；用户的延迟任务和原有 PRSformer 进程均未被终止。生产代码的 SHA256 与测试开始时完全一致。

时间、进程、源码校验和结束阶段记录见 `/tmp/grid-no-annotations-validation/final-window-result.json`。本次并未完成全量 LDL、T2DM 以及全量三性状汇总报告，不能把限时测试表述为整个正式分析完成。

| 检查 | 结果 |
|---|---|
| 默认 GRID 环境 CUDA attention 单元检查 | 6/6 通过 |
| 无注释三性状 GPU 端到端 | 13/13 通过 |
| 显式 canonical 注释三性状 GPU 端到端 | 13/13 通过 |
| 无注释三性状 CPU reference 对照 | 13/13 通过 |
| 匹配、来源、注释和原生 R 回归 | 34/34 通过 |
| 基线完整性 | 17/17 通过 |

端到端检查覆盖真实 PLINK/PGEN 合成夹具、RDS、prepare/evolution/fit/report/predict、新进程加载模型和缓存、测试结局隔离、donor 贡献重构。无注释模式还验证：不存在 `evo.*` 特征、年龄模型和模块评分产物；原 CSx 分数保持不变；权重来源被篡改时仍拒绝。

执行环境：`/home/huangj/miniforge3/envs/grid/bin/python3`，PyTorch `2.12.0+cu130`，RTX 5090 Laptop GPU，默认 `cuda:0`。

## 真实数据测试

全量三性状默认流程已启动，参数仅改变缓存和输出目录。日志：`/tmp/grid-no-annotations-validation/full-real-run.log`。真实队列准备通过：486,745 人，训练 243,373 人、测试 243,372 人。与前一版保存的名单逐项比较，参与者及训练／测试归属完全一致。

全量 height 已完成 fit、242,649 人的冻结测试预测及 `fit.joblib` 保存。另用新 Python 进程加载该模型，抽取 1,000 人重新预测，CSx、GRID 主模型及最终 policy 均与保存值一致（绝对与相对容差 `1e-10`）；全部测试预测有限，全部 donor 贡献可重构。实际设备 `cuda:0`，选定 epoch 3，共 3,983 次梯度更新、45 个状态张量发生变化，峰值已分配显存 182.6 MiB。记录见 `/tmp/grid-no-annotations-validation/full-height-verification.json`。

全量 height 的独立增益审计状态为 `insufficient_audit_information`，policy 按既有规则回退到 CSx，没有报告已证实的预测增益。

全量 LDL 在 `cuda:0` 完成 8 个训练轮次、3,824 次梯度更新，匹配参数选择、可靠性筛选及全训练对照拟合均完成。窗口结束时正在执行 203,630 人的冻结测试预测，尚未保存完整 fit 结果；全量 T2DM 拟合和全量三性状 report 尚未开始。上述未完成部分需要正式任务继续验证。

另从共同队列中按固定种子、完全不读取结局，分别抽取训练和测试各 4,000 人；保持其原有外层训练／测试归属。这个 8,000 人子集使用默认模型、训练及报告参数，三个性状均完成 prepare/fit/report，进程退出码 0。日志：`/tmp/grid-no-annotations-validation/real-subset-run.log`。输出用于验证代码，不能作为正式 UKB 预测性能结论。

随后三个模型在新的 Python 进程中运行 predict，并与保存的预测逐项比较（绝对和相对容差均为 `1e-10`），donor 贡献重构也通过；T2DM 概率均在 `[0,1]`。没有年龄特征或年龄模型。

| 性状 | 测试人数 | 选定 epoch | 实际梯度步数 | 变化的状态张量 | 模型 GPU 已分配显存峰值 |
|---|---:|---:|---:|---:|---:|
| height | 3,984 | 4 | 80 | 45 | 106.76 MiB |
| LDL | 3,369 | 2 | 48 | 45 | 109.89 MiB |
| T2DM | 4,000 | 5 | 90 | 45 | 107.72 MiB |

各性状测试人数差异来自结局缺失；共同外层名单仍为 8,000 人。详细结果：`/tmp/grid-no-annotations-validation/real-subset-verification.json`。

子集的 12 个工作簿均可读取、9 张图像均有效；汇总工作簿未出现个人 ID 表头。记录见 `/tmp/grid-no-annotations-validation/real-subset-artifacts.json`。

全部测试日志位于 `/tmp/grid-no-annotations-validation`。源码中不再存在原数据库的读取或下载入口，当前 README 已删除其下载步骤；历史验证记录保留当时发生的事实，并注明已被新流程取代。
