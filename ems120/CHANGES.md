# 本次代码修改

## 目录与调度

完整项目从 `0LOCAL/ems120` 移到 `ems120`，共享代码仍在 `0f`。新版 shell 不再调用强制 local 的 `local_paths.sh`，也不再从 `/mnt/d/scripts/0f` 取 helper。默认模式明确为 HPC，显式 `--submit local` 使用同一个 worker。

`f/dispatch.py` 只做参数、路径和资源校验，再生成可检查的 `plan.json` / `.cmd`；计算依赖在 `f/run_worker.sh` 中激活。调度 helper `f/workflow_jobs.py` 是共享 helper 的项目内快照，带前次改进的主机内存和 GPU 按 host 请求修复，并支持 EMS 独立临时目录。因此本包无须覆盖已经配置好的 `0f/0os.sh` 或其他 project。

## 文件职责

| 文件 | 修改作用 |
|---|---|
| `ems120.sh` | 新位置引导、读取共享目录表、明确 HPC/local、轻量 planner |
| `f/dispatch.py` | 路径/资源/阶段验证，冻结运行参数，建立单 GPU 作业 |
| `f/run_worker.sh` | LSF 分配检查、双目录锁、Conda/R/Python 环境绑定、实际运行 |
| `f/workflow_jobs.py` | 提交回执、日志、状态、退出码和信号传递，未知提交结果不盲目重试 |
| `f/runtime.R` | R 运行路径、依赖检查、reticulate 设备约束、RDS 完整发布 |
| `f/check_runtime.R` | 实际 R/reticulate Torch CUDA 运算；显式 GPU 检查时附加 XGBoost 探测 |
| `f/ems120.R` | 去除本机目录、自动安装与旧 Python 查找路径；保留原分析步骤 |
| `f/ems120.py` | 模型目录、GPU 策略、线程/批量控制、电话模型缓存复用、输出保护 |
| `tests/` | 调度、运行环境、GPU 策略、错误传播及缓存保护回归测试 |

## 保留的分析方法

研究年份、疾病映射、phone 规则评分、训练样本来源、五折 CV 定义、随机种子、Transformer 训练轮数及其批量设置保留。仍按交叉验证 RMSE 选择 phone 赢家。XGBoost CUDA 是原候选的执行后端调整，并没有强制将所有 CPU 模型替换成另一种算法。

完整分析仍按原依赖顺序进行，避免把有共享状态的 phone/geo/dx 阶段同时提交后覆盖同一训练文件或 RDS。不同数据分析可以使用不同输出目录分别提交。

## 关键可靠性修复

- Torch 与 XGBoost 分别检查实际 CUDA 能力；显式 CPU 调试只允许 local。
- 保留分配的 GPU 可见性，要求一张逻辑可见卡；不在 HPC 指定物理 GPU 编号。
- 原 `ml_phone` 训练失败会返回全 NaN；新版将失败传给 R 和调度状态，避免把错误结果缓存成成功。
- 原 phone 预测批次重复加载同一模型；新版在进程内复用，且约束已加载模型的 CPU 线程数。
- RDS、JSON、phone 模型单文件完整发布；phone 成对文件有未完成标记；Transformer 模型先完整保存后发布。
- 分析目录和外置模型目录均加锁；不同计划不能同时写相同结果。
- HPC 分析/模型默认放持久 `/work`，临时文件使用共享 scratch 配置。

## 官方接口依据

- XGBoost `device=cuda` 与 `tree_method=hist`：<https://xgboost.readthedocs.io/en/release_2.1.0/gpu/>
- `device` 参数自 2.0.0 引入：<https://github.com/dmlc/xgboost/blob/master/doc/parameter.rst>
- LSF GPU 请求参数：<https://www.ibm.com/docs/en/spectrum-lsf/10.1.0?topic=o-gpu>

队列硬件规格来自用户附件；实时可用性仍以启明账号查询为准。环境、数据和地图相关依赖见 `requirements-notes.md`，实际测试边界见 `VALIDATION.md`。
