# GRID GPU ABM 验证记录（2026-10-09）

> 历史记录：本记录早于当天移除默认外部年龄数据库依赖的修改。当前流程及验证见 [无外部注释验证记录](VALIDATION.grid-no-geva-20261009.md)；下文的缺失数据阻挡仅描述当时版本。

本次将 GRID 默认 ABM 改为借鉴 LE8 c1 的模块 Transformer 和可学习 donor Q/K attention。遗传特征提供 token，CSx 的家系 OOF 残差提供 donor value；独立模型选择、gate、审计及外部测试集继续保留。PRSformer 仍由独立脚本运行，GRID 使用 GPU 不需要 `--run-prsformer`。

## 实际运行环境

- 默认解释器：`/home/huangj/miniforge3/envs/grid/bin/python3`，Python 3.11。
- PyTorch：`2.12.0+cu130`，CUDA runtime 13.0。
- GPU：NVIDIA GeForce RTX 5090 Laptop GPU。
- `bash grid.sh --check-device`：真实模块 Transformer、donor Q/K 和 decoder 的前向/反向计算通过，关键参数存在非零梯度，设备为 `cuda:0`。
- `python -m pip check`：No broken requirements found。
- 安装包从官方 PyPI 下载，下载器校验各 wheel 的 SHA256 后离线安装。未修改 LE8 或 PRSformer 环境。

## 自动测试

| 检查 | 结果 | 说明 |
|---|---|---|
| 默认 GRID 环境 `f/test_grid_attention.py` | 6/6 通过 | 实际 CUDA；连续/二分类训练、梯度、冻结和重载、结局隔离、贡献重构；GPU 检索对照 NumPy，涵盖 ties、家系及 OOF 折排除；空候选的有限梯度；禁止隐式 CPU 回退 |
| 默认 GRID 环境 `GRID_TEST_ABM_BACKEND=selective_attention ... f/test_grid_integration.py` | 13/13 通过 | 三性状合成数据，真实 PLINK/PGEN、RDS、prepare/evolution/fit/report/predict；新进程加载模型及复用 fit/report 缓存 |
| 显式 CPU reference 集成流程 | 13/13 通过 | 保留旧对照算法兼容性 |
| `test_*review.py` | 34/34 通过 | 原有匹配、数据及注释约束，以及含 RAW 属性的原生 RDS 表读取 |
| `f/test_baseline_integrity.py` | 17/17 通过 | 基线完整性与来源约束 |
| Bash/Python 语法 | 通过 | 修改的脚本及新增模块 |
| shell 错误路径 | 通过 | 无 CUDA 时在读大表前拒绝；缺注释时提前报错；保存运行日志并保留非零退出状态 |
| 大参考库检索 | 通过 | 100,000 donors × 1,024 queries，k=64，1.93 秒，峰值已分配显存 83.2 MiB；抽样独立 NumPy 对照一致；在本机 LE8 的相同 CUDA PyTorch 版本上测量 |
| 默认网络规模与 dropout | 通过 | GRID 环境；width=64、heads=4、layers=2、dropout=0.1、batch=256，真实 CUDA 执行 3 次 AdamW 更新；重复运行参数逐项完全一致 |

另在已有 LE8 CUDA 环境中独立执行了同样的 6 项 attention 与 13 项端到端测试，均通过。测试中发现并修复了动态模块反序列化失败，以及非确定性 CUDA 计算导致的微小重复运行漂移；没有放宽原有端到端结局隔离断言。

默认 GRID 环境的详细日志位于 `/tmp/grid-default-gpu-preflight.log`、`/tmp/grid-default-attention-tests.log`、`/tmp/grid-default-gpu-integration.log` 和 `/tmp/grid-pip-check.log`。其他回归日志位于 `/tmp/grid-review-tests.log`、`/tmp/grid-baseline-integrity.log`、`/tmp/grid-native-integration.log`。

## 真实输入与执行边界

大型表型 RDS 改为通过 R 读取后只向 Python 传入需要的列。实际 `/mnt/d/data/ukb/phe/Rdata/all.rds` 测得 502,371 行、9 列，19.89 秒，返回的 Python 表约 61.3 MiB；这不是 R 子进程的总内存峰值。

真实 `prepare` 首次检查还发现已存在的 CSx RDS 含有 `rdata` 不能解析的 RAW 属性。现在所有输入 RDS 表优先走原生 R 读取，并保留字符串列类型（包括带前导零的个人及家系 ID）；已增加原生 R 生成这类文件的回归测试。

修复后，真实三性状 `prepare` 已完成：共同队列 486,745 人，训练集 243,373 人、测试集 243,372 人。三个性状来源记录补齐后，使用最终代码重新执行该检查，仍成功，退出码 0。使用独立缓存 `/tmp/grid-abm-real-prepare-20261009` 和输出目录 `/tmp/grid-abm-real-output-20261009`，不占用延迟任务的默认缓存。日志为 `/tmp/grid-real-prepare.log`。这是输入和队列检查，尚未进行真实 UKB 的年龄特征构建和 ABM 训练。

为补齐真实 CSx 分数与后验权重的来源记录，使用原权重执行 `1.csx.sh --traits height,ldl,t2dm --stage score --models populations --posterior FALSE --jobs 4 --threads 2`，三个性状均完成，进程退出码 0。本次没有重新执行 MCMC。原分数备份位于 `/tmp/grid-abm-before-20261008/data`，日志位于 `/tmp/grid-csx-provenance-rescore.log`。

三个性状的 `canonical_weights()` 均通过来源验证，覆盖 22 条常染色体；height/LDL/T2DM 分别有 142,024 / 147,751 / 147,703 个变异。逐性状按 `eid` 对齐全部 487,162 人后，新旧 `csx.EUR/AFR/EAS/SAS/auto/meta` 六列均完全一致（零绝对和相对容差，保留缺失值）。摘要见 `/tmp/grid-all-real-input-check.json`。重新评分只补齐可验证的来源链，没有改变这些数值。

**完整 UKB 年龄分析仍受 GEVA 注释缺失阻挡。** 默认 `/mnt/f/ref/GEVA` 中缺少 `atlas.chr1.csv.gz` 至 `atlas.chr22.csv.gz`；本机访问[官方 bulk 下载端点](https://human.genome.dating/download/index)失败。作者发布于 [Figshare 的早期 SGDP 表](https://figshare.com/articles/dataset/Atlas_SGDP/7098680/1)只含年龄点估计等字段，缺少当前方法需要的区间和质量分数（已核查实际 TSV 表头），因此未用它替换当前 Atlas summary，也没有开启 `--allow-proxy-only` 或改变分析定义。

自动测试通过证明上述代码路径可执行，不代表已跑完真实 UKB 的 evolution/fit/report，也不证明模型改善预测性能。实际 `prepare/evolution` 的表格处理及 PLINK、部分基线和 gate 使用 CPU；attention 训练、donor attention 与默认分块检索使用 CUDA。

最终 `bash grid.sh --traits height,ldl,t2dm --check` 的退出码为 2：提前指出缺少 GEVA chr1–22，不再报告评分来源记录缺失，也未载入大型表型表。详细错误保存于 `/tmp/grid-final-prerequisites.log`。同时独立 `bash grid.sh --check-device` 仍通过。再次探测官方 bulk URL 仍返回 HTTP/2 PROTOCOL_ERROR，不能把 GPU 检查通过当作全量输入检查通过。

每次正式执行的日志自动写入 `CACHE_DIR/logs/grid.*.log`。代码原始备份位于 `/tmp/grid-abm-before-20261008`。用户提交的 `sleep 2h` 进程没有被终止或替换；其后续命令将读取启动时磁盘上的新版代码，缺少 GEVA 时会明确失败并留下日志。
