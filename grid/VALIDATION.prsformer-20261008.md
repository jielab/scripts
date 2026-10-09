# PRSformer 本机安装与验证（2026-10-08 晚至 2026-10-09 凌晨）

修复起点：默认 `~/.venvs/grid-prsformer/bin/python` 与官方源码目录不存在。机器为 RTX 5090 Laptop GPU，旧的 CUDA 12.4 / NATTEN 0.17.5 组合不适配该 GPU。

## 安装与代码

- 独立环境：`/home/huangj/.venvs/grid-prsformer`，Python 3.12.13，torch 2.12.0+cu130，NATTEN 0.21.7+torch2120cu130，CUDA 13.0。
- 官方源码：`/mnt/f/software/PRSformer`，commit `7dea4e1bb27975885c937f2be82f9243bc82bda7`；源码文件未修改。
- NATTEN wheel 与官方 GitHub release 的 SHA256 一致：`702d435bcdf84145c2407af7d2267b127f77249d83e9cb6d381ec0b3f8a2cb6c`。
- 为避免重复慢速下载，从已有独立 AI 环境的已安装文件重建了本地 torch/CUDA 依赖 wheel，逐文件核对 RECORD 哈希后安装；没有让新环境继承 AI 环境的 site-packages。新的 setuptools 与数据依赖来自 PyPI。
- `pip check` 通过；57 个软件包的 25,132 个安装文件（5.42 GB）逐一通过 RECORD 哈希复核。所有直接依赖版本记录在 `requirements.prsformer.txt`，完整版本冻结文件见日志目录。
- 新增可重复安装入口 `install.prsformer.sh`；加载器显式使用新版 NATTEN `cutlass-fna`，保持官方 QKV/投影参数键与局部注意力语义。
- 默认 `cuda` 规范化为 `cuda:0`；新增 `--check-runtime`；all 模式在读取真实数据前运行真实 GPU kernel 检查。
- 提前丢弃无关表型列，减少内存复制；基因型准备、转置与训练过程增加进度输出。
- 修复整仓回归中 `f/grid.data.py` 的 Rscript 多行 `-e` 参数截断问题：改用临时 R 文件，保留既有按列读取逻辑。修改前版本另存为备份目录内 `grid.data.before-rscript-fix.py`。

## 已完成验证

1. **默认 GPU 配置**：64 维、4 heads、2 层、FF=128、窗口 385、FP16、gradient checkpointing，真实 forward/backward 通过。显式 CPU/global 路径也通过。
2. **4 项原生 GPU 回归测试全部通过**：含 12 组窗口/dilation/精度组合的独立局部注意力输出和梯度比对；原生 PGEN 跨染色体样本重排；三性状训练/预测/报告；覆盖保护；进程中断后恢复及历史一致性。
3. **默认训练控制**：合成数据使用默认 30 轮上限、patience=5，第 20 轮正常早停，选择第 15 轮，完成报告发布。
4. **全队列只读预检**：默认原始路径全部可用；486,745 人、147,792 个候选 SNP，22 条染色体与 HM3 坐标零不一致。矩阵两种布局加保留空间的上界为 268.49 GiB，当前缓存磁盘充足。
5. **真实全 SNP 端到端试跑**：2,048 人，全部 22 条染色体，train-only MAF/call-rate QC 保留 145,136 SNP；两轮训练与验证完成，413 人的独立测试集预测、RDS 评分和 XLSX/PNG 均生成成功。其余模型参数与默认命令一致。
6. **独立剂量抽查**：PLINK 原生 `--export A --export-allele`，22 条染色体、30 人、66 个 SNP，共 1,980 个 ALT 剂量与缓存一致（float16 允许 0.001 绝对舍入误差）。
7. **全长 GPU 压力测试**：147,792 SNP、12 次 forward/backward，47,360,515 个参数；峰值 allocated 1.18 GiB、reserved 1.64 GiB。此为 batch=1 的局部压力测试，不是全队列运行峰值保证。
8. 缺少 Python、缺少官方源码、CPU/neighborhood、不兼容的 resume/replace 均能明确退出；Bash 语法、ShellCheck、Python AST 检查通过。

## 默认参数持续训练与最终回归

- 真实 2,048 人 × 145,136 SNP 使用默认 30 轮上限、patience=5、学习率、模型、FP16 和 batch/累积设置，在第 11 轮正常早停，选用第 6 轮模型；所有训练和验证损失均为有限值，未出现 OOM。
- 单轮训练加验证约 47–50 秒。按人数比例外推，正式全队列每轮约 3 小时量级；这是小样本吞吐量外推，真实磁盘、温度与其他负载会影响耗时。
- 重载上述保存的模型，重新预测全部 413 人；1,239 条“个体 × 性状”结果与训练结束时的输出一致，再次完成 RDS/XLSX/PNG 发布。
- 独立复核真实评分文件的 test IDs、源协变量、endpoint、preparation signature、预测分解，并用 sklearn 重新计算 R²/AUC，与报告一致。
- 整仓最终回归在独立 PRSformer 环境中运行：78 项中 74 项通过，4 项显式 GPU 测试在该次运行中跳过；这 4 项已单独启用并全部通过。日志为 `final-stable-suite.log` 与 `gpu-regressions-v3.log`。最终测试前后源码哈希一致，详见 `final-stable-suite-source-hashes.json`；此前一次受并发源码修改影响的缓存保护退出不属于最终结果。
- 从已激活的 GRID 环境调用默认 PRSformer 入口仍能使用独立 Python 并通过 GPU 检查。
- 最后重跑原始默认路径的 `--check` 成功，再次确认 486,745 人、147,792 候选 SNP、268.49 GiB 空间上界。
- 全人数补充实测通过：486,745 人 × 512 候选 SNP，train-only QC 保留 504 SNP，完成完整矩阵准备及全部协变量基线拟合，模型 forward/backward 预检通过。

这些验证没有启动或完成 486,745 人 × 全 SNP 的正式神经网络训练；正式运行由用户已提交的延时命令启动。

验证完成时间：2026-10-09 00:23（Asia/Shanghai）；自用户指定的约 55 分钟测试窗口开始持续验证。用户的延时任务仍保留，预计约 00:27:45 启动。

## 文件位置

- 修改前备份与安装/测试日志：`/mnt/d/Downloads/GRID/prsformer-runtime-before-20261008-231848/`。
- 通过的合成测试：`/tmp/grid-prsformer-validation-20261008/synthetic-v3/`。
- 真实数据试跑：`/tmp/grid-prsformer-validation-20261008/ukb-full-snps/`，使用独立 cache、run、reports、scores。
- 正式默认 `/tmp/grid-prsformer`、`/mnt/d/analysis/grid` 与 `/mnt/d/data/ukb/pgs` 未被这些测试写入。

## 解释范围

这些检查验证安装、设备兼容性、数值操作、数据方向、训练恢复和输出合同，不代表短期试跑模型具有可用于研究结论的预测性能。真实小样本试跑的验证损失仍较大，不能把运行成功解释为模型已达到合格预测性能；正式全队列模型需要独立评估。pyreadr 导入原始 RDS 的无关日期列时会产生 NA 转换 RuntimeWarning；所需数值表型、协变量及最终输入检查通过。
