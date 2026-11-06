# QiMing 提交资源约定

按管理员 2026-11-05 的建议维护。先执行 `queueinfo` 核对队列；本表描述单节点请求，不能把多个节点的内存相加供单进程使用。

## 33 队列

每节点 40 核、192 GiB。令 `f = max(所需核数 / 40, 所需内存 GiB / 192)`，按较大的占比向上取整：

| 占比 | 申请核数 `-n` |
| --- | ---: |
| 0 < f ≤ 25% | 10 |
| 25% < f ≤ 50% | 20 |
| 50% < f ≤ 75% | 30 |
| 75% < f ≤ 100% | 40 |

24 核、160 GiB 的任务应申请 40 核，内存仍申请 160 GiB。4 核、64 GiB 的任务应申请 20 核。代码保留此前为操作系统预留空间的预算：33 队列默认允许单个任务最多 160 GiB；超过预算需选择更大队列。

独立 CPU worker 的参数示例（将 `worker.sh` 换成实际 worker）：

```bash
bsub -q 33 -n 40 -W 24:00 -M 160GB \
  -R 'span[hosts=1] rusage[mem=160GB/host]' bash worker.sh
```

`/host` 表示整台节点上的任务总内存。不要把 160 GiB 当作每核内存再乘以 40。需要 per-task 配额的现有提交器，会用取整后的核数计算每核份额。

## 队列选择和自动处理

`qiming_resources.py` 统一处理新 CPU 计划：普通已知 CPU 队列按各自节点规格取整；管理员明确推荐的 `v4-3t`、`v3-6t` 保留实际所需核数。其他已知 CPU 队列暂按分档规则保守处理，除非管理员进一步确认例外。GPU 任务及未知队列不套用虚构的 CPU 节点规格。

```bash
python3 /work/sph-huangj/scripts/0f/qiming_resources.py \
  --queue 33 --cores 24 --memory-gib 160 --hours 24 --format bsub
# 输出 -q 33 -n 40 ...

python3 /work/sph-huangj/scripts/0f/qiming_resources.py \
  --queue v4-3t --cores 24 --memory-gib 160 --hours 24 --format bsub
# 输出 -q v4-3t -n 24 ...
```

自动选择 `--queue auto` 仍遵守各队列保守内存预算和最长时限。资源分档只改变 LSF 预约核数，不自动增加分析线程、并行样本数、内存额度或 GPU 数。

## 各项目的使用方式

下面是加到已有分析命令中的资源参数，数据、模块和输出路径仍使用该项目原来的参数。优先用原生提交入口，它们管理依赖和提交记录；不要在其外层再套一次 `bsub`。

| 项目 | 160 GiB CPU 任务的资源参数 / 说明 |
| --- | --- |
| allele3 | `--submit hpc --queue 33 --job-cores 40 --job-memory-gb 160 --job-hours 24 --write-cmd` |
| grid | `--submit hpc --queue 33 --job-cores 40 --job-memory-gb 160 --plan-only`；GPU 阶段使用单独配置 |
| gwas | `--submit hpc --queue 33 --cores 40 --memory-gb 160 --hours 24 --plan-only` |
| gu / IBDmix | `--GPU FALSE --submit hpc --queue 33 --slots 40 --memory-cap 160G --run-cmd FALSE` |
| gu / hpc plan | `--queue 33 --threads 24 --memory-gb 160`；申请 40 slots 后再用原来的 submit 命令 |
| ld | `LD_MEMORY_GIB=160 bash ld.sh --grch 38 --chr 1 --submit hpc -q 33 -n 40 --write-cmd` |
| le8 | `--submit hpc --queue 33 --cores 40 --memory-limit-gb 160 --hours 24 --plan-only`；GPU 阶段单独配置 |
| metag | `submit=hpc bash 1.library.sh --config hpc.yml --submit FALSE`；按类别生成不同资源的命令 |
| ems120 | GPU 工作流：`--submit hpc --queue 4a100-40 --cores 8 --memory-gb 96 --gpus 1 --plan-only` |
| dose / maha / nextstrain / pfi / phole / ukb | 主入口仍初始化本地路径；已补资源说明，不把它们宣称为已验证的 HPC 工作流 |

查看生成的 `.cmd`、计划和路径之后，通过原入口移除 `--write-cmd` / `--plan-only`，或按各项目的 receipt-aware 提交方式提交。GU 使用 `--run-cmd TRUE`；metag 使用 `--submit TRUE`。对于已保存的旧计划，应创建新计划才能采用新资源规则；读取旧计划不会悄悄改写它的资源。

`ukb/gen.sh --background` 是本地入口中保留的特殊旧提交分支，三个分支已明确使用 `33 / 40 核 / 160 GiB / 24 小时`；仍需显式提供 HPC 数据、输出及软件路径。其科学运行未因此被验证。

直接执行的安装脚本与分析计划器不同。其文件头给出独立安装 worker 的 `bsub` 示例；资源按安装所需量设置。GPU 安装检查应使用相应 GPU 队列。

## metag 当前配置对应的预约

| 阶段 | 队列 | 申请核数 | 总内存 GiB |
| --- | --- | ---: | ---: |
| setup | 33 | 10 | 32 |
| prepare | 38 | 16 | 64 |
| build viruses | 33 | 10 | 48 |
| build microbes | 38 | 64 | 448 |
| build plants / invertebrates | 1t81c | 120 | 896 |
| build vertebrates | v3-6t | 64 | 1792 |

这是 `hpc.yml` 现有资源配置经分档后的结果，不是对实际峰值内存的测量。2026-11-05 的 prepare 失败是 `Disk quota exceeded`；增加核数不能解除磁盘配额限制。先解决组配额，再生成新计划复用已有有效准备结果。
