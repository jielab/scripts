# flux

基于 jielab/scripts 的 0bdec70（2026-11-26）。完整替换仓库中的 flux/，沿用相邻的 0f/ 和已有的 ai Python 环境。在 flux/ 下运行 ./flux.sh -h 和 ./flux.paper.sh -h 可查看示例。

## Python 环境

默认复用已有的 `ai` Conda 环境（Python 3.11–3.13），不创建独立环境，也不替换 Python。
`./install.sh run --plan-only` 显示 KEEP/NEEDS；`./install.sh run` 只补充缺失或版本范围不满足的依赖。
兼容的 Torch/CUDA 保持原样；安装时约束其他已安装包，依赖冲突则停止，不自动升级共享环境。
`./install.sh run --check` 只检查版本、导入和 CPU/GPU 前向反向运算。
版本范围是兼容性初筛，检查通过不代表完整数据分析已经验证。

```bash
./install.sh run
./flux.sh run
# 非标准位置：安装和运行使用相同覆盖值
export FLUX_ENV=/path/to/conda/envs/ai
```

`--cuda` 在已有 Torch 上用于显式校验，默认 auto 接受能运行的 CUDA 版本；
缺少或版本不兼容的 Torch 才选择安装 wheel。GPU 驱动或二进制加载错误会报告，不自动重装。
共享环境中与 FLUX 无关的 `pip check` 冲突仅报告，不自动修复。

## 从零运行

先用 bkill 停止本批仍在运行或排队的 flux 作业。用 bjobs 确认它们全部退出后，删除输出目录 `/work/sph-huangj/analysis/flux`，解压本包替换 `scripts/flux/`。在 HPC 上运行

```bash
cd /work/sph-huangj/scripts/flux
./flux.sh run
```

程序会新建默认输出目录，从原始 YGB 输入构建队列，训练主模型，再完成评估与两项消融。输入由相邻仓库目录 `0f/` 中的共用配置解析。

每次全新主分析共 8 个 LSF 作业。作业依赖如下，符合依赖后仍需等待队列分配资源。

| 作业 | 前置作业 | 资源和用途 |
|---|---|---|
| flux.prepare | 无 | CPU，准备队列 |
| flux.train | prepare | GPU，拟合主模型 |
| flux.readout | train | GPU，提取状态、训练预测头 |
| flux.dynamics | train | GPU，评估下一事件及时间 |
| flux.ablate | train | GPU，依次重训并评估 no-ODE 和 no-update |
| flux.evaluate-baseline | readout | CPU，先评分主模型，再拟合及评估基线模型 |
| flux.blackout | readout | GPU，评估无新观测时的预测 |
| flux.rollout | readout | GPU，模拟后续事件 |

单独运行 evaluate 和 baseline 的入口仍然可用。对只需主模型的运行使用 ./flux.sh run --no-ablation。预览调度用 ./flux.sh run --plan-only。

## 汇总

等所有主分析作业成功完成后，另行提交：

```bash
./flux.paper.sh run
```

打开输出目录中的 paper/main/reports/index.html 查看汇总。

## 性能与验证

- 训练及冻结状态推断按模态打包实际有数值测量的访视。诊断时间、缺失掩码、同日事件顺序、训练窗口和目标沿用原定义。
- 每个微批一次性取回历史区间的 RK4 步数。每个区间把相同的年龄和背景投影用于多个 RK4 子阶段。隐藏状态仍逐步积分，保持 float32。全局批量为 3,072，默认最多 2,000 个优化步，每 50 步用完整验证集评估。
- 训练日志每 10 步给出平均秒数（s/update），验证行另列 validation_seconds。比较两版每 50 步总耗时，用相邻两条验证完成日志的时间差。
- 在 PyTorch 2.6.0、2 个 CPU 线程下，以相同合成 SQLite 数据、完整模型宽度、32 人、96 个时间点和 5,389 个数值特征做前向和反向对照。交替测量、剔除各自首轮后取 5 次中位数。读取与组批由 0.134 秒降至 0.023 秒，读取、组批、前向和反向合计由 0.794 秒降至 0.556 秒。测试批的数值载荷从 82,775,040 字节降至 149,936 字节。实际比例取决于每人的观测密度。
- 核对了同一随机种子下的损失与参数梯度、两种消融和跨块历史推断。320 人合成队列已运行主分析和独立论文汇总。批量矩阵运算可能产生浮点舍入差异。真实 YGB 的指标与 V100 用时须以 HPC 日志为准。
