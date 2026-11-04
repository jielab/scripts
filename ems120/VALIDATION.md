# EMS120 迁移验证记录

基线：`jielab/scripts@bf7368b5fd72c9b85b55ea8e942062531e165d59`。
新位置：`scripts/ems120`。测试针对随包代码，不涉及患者数据或真实训练。

## 已实际执行

| 检查 | 结果 | 覆盖范围 |
|---|---|---|
| Python 标准库 unittest | 28 项全部通过 | 15 项 CLI/调度/worker，13 项模型策略/缓存保护 |
| Bash 语法 | 2 个入口全部通过 | `ems120.sh`、`f/run_worker.sh` |
| Python 语法 | 包内 Python 文件全部通过 | 静态编译，不导入 Torch 或科学模块 |
| R 4.3.3 语法 | 7 个文件全部通过 | 主分析、两个 R runtime helper、R 测试及三个共享 R helper |
| 基础 R 测试 | 3 组全部通过 | 路径优先级、RDS 完整发布/失败保留旧缓存、缺少包时明确失败 |
| Python 科学函数 AST 对比 | 14 个关键函数与基线相同 | 评分、特征、文本、关键词、统计指标和 geo 训练数据构造 |

28 项 unittest 使用模拟 `bsub`、Conda 激活脚本、记录参数和环境的 R 替身，以及假的 XGBoost/Torch 接口。它们实际运行了 shell → planner → worker → 子进程链，但没有将测试替身当作真实 LSF/CUDA 运行。

具体检查包括默认 HPC、显式 local、路径带空格、资源校验、禁止自动回退、GPU 检查必须调度、提交回执、失败退出、禁止登录节点直接执行 HPC `.cmd`、分析/模型目录锁、缓存绘图重跑、临时文件位置、设备可见性和研究参数冻结。模型测试检查 CPU 回退检测、旧/非 CUDA XGBoost 拒绝、缓存 CPU 赢家保留、线程限制、模型重复加载消除、失败不会变成全 NaN，以及文件/目录发布失败时旧版本恢复。

基础 R 测试通过已有 R 4.3.3 执行；没有安装或加载完整科研分析环境。

## 没有完成或不能据此保证的事项

- 未连接启明账号、未实际向启明提交作业，无法验证账号队列权限、收费、实时可用节点或站点配置差异。
- 此环境没有真实 CUDA GPU / Torch / XGBoost 分析环境，实际 GPU 检查需在目标 worker 执行。
- 未提供 EMS Excel、预训练模型或旧模型，因此未运行完整研究流程，也未测得实际峰值内存和总耗时。
- GPU 重训不保证与 CPU 结果逐位相同；模型候选的 CUDA histogram 后端及浮点差异可能影响数值很接近的比较结果。
- 旧缓存按原缓存规则复用，不包含全量输入内容指纹；更换数据或标签时使用新分析目录。

## 在目标仓库复核

```bash
cd /path/to/scripts
bash -n ems120/ems120.sh
bash -n ems120/f/run_worker.sh
python3 -m unittest discover -s ems120/tests -p 'test_*.py' -v
Rscript --vanilla ems120/tests/test_runtime.R
bash ems120/ems120.sh --check-gpu --wait
```

最后一条是真实 GPU 验收任务，默认检查 R/reticulate CUDA 运算和 XGBoost CUDA 拟合。通过后再提交完整分析，并核对代表性评分、类别分布、CV 表和图形。

打包验证还检查迁移补丁能应用于基线、补丁应用结果与完整代码目录一致、ZIP 无旧 `0LOCAL/ems120` 代码目录、归档 CRC 和逐文件 SHA-256 校验。实际日志随 ZIP 放在 `verification/`。
