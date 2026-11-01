# GU FINAL 2026 验证

本次合并回归：**66 项测试通过，18 项子测试通过**；HPC 另有 **29 项测试通过**。没有将可选测试跳过计入通过数量。测试数据为合成输入，不含私人基因型。

| 范围 | 数量 | 实际覆盖 |
| --- | ---: | --- |
| IBDmix | 23 | 官方原生 prepare 对照、实际 CUDA 源码的主机算法对照、缺失/单倍体/重复位点、默认五参考、显式参考、精确 DAF provenance、严格续跑、真实 Tab 与批次记录 |
| PhyML / UKB | 13 | 染色体精确匹配、重复 archaic 位点、lead 等位基因翻转、男性 X、LD 样本、重复变异、stderr 管道、实际 pruning 源码对照、当前 CLI |
| ARG / TRACE | 13 | 真实固定上游的合成 ARG 端到端，chr1+chr2 联合 HMM，男性非 PAR X 独立模型，bedtools mask，输入/地图/样本身份，严格续跑，inference 文件完整性 |
| 汇总 / 主 CLI | 8 | 并集分母、成功零调用与未检测区别、完整同男性比较、染色体范围、CPU/CUDA 统一目录、默认五参考、显式额外参考、单一模块名、TRACE 联合请求 |
| 新结果数据库 | 9 | 实际创建 schema 2026，三方法 SQL，完成但零片段，动态参考集合，实际 roster，当前 PhyML 字段，X 非 PAR 分母 |
| HPC | 29 | mock LSF、JOBID 依赖、输出隔离、冻结样本与配置、批次标签、HPC 路径、ARG ancestry、失败/零调用区分、严格汇总、归档去重、活动作业期间拒绝重试 |

## 可复现运行

从包顶层（`gu/` 同级）运行。先安装 GU 科学依赖，并为测试安装 `pytest`、`pytest-subtests`、`msprime`。真实 TRACE 集成还需要 `trace-extract`、`trace-infer`、`trace-summarize` 与 `bedtools` 在 PATH。缺少依赖的本地测试会明确 skip，不应报告为已完成相应验证。

```bash
# 指向正确版本的上游源码与已编译准备程序。
export IBDMIX_UPSTREAM=/absolute/path/to/IBDmix
export IBDMIX_PREPARE_BIN=/absolute/path/to/ibdmix-prepare
export PHYML_UPSTREAM=/absolute/path/to/phyml

# planning fixtures 写到单独的持久临时目录，不混进代码目录。
mkdir -p /absolute/path/to/gu-test-output
export GU_TEST_OUTPUT_ROOT=/absolute/path/to/gu-test-output
python3 -m pytest -q gu/tests

# HPC 包也可单独运行：
python3 -m unittest discover -s gu/tests -p test_hpc.py -v
```

单独运行非 TRACE 子集也可使用 `python3 -m unittest discover -s gu/tests -p test_ibdmix_core.py -v`；新结果测试与 ARG/TRACE 测试使用 pytest。

官方对照版本：

- IBDmix：`6a9d4bc67dcaf602e65cac0fa6191e1e7a1291ef`。
- PhyML：`42bb4b52f479938a3427dfe4357b861f0e306686`。
- TRACE：`70cab6c6687900826eb0b2ded69c39554b1b4146`。

原生 prepare 使用 HTSlib 1.22.1 构建，含严格 C++ 编译检查。算法对照包含 16,384 个站点 LOD/AF、1,024 组随机分段 trial ×3 样本及 2,048 个 PhyML pruning 配置。

主 CLI planning 测试只模拟昂贵的 mask 预检查，执行真实命令生成逻辑。CUDA 算法对照编译实际 `.cu` 源码中的数值代码为主机 C++，**未执行 GPU device kernel**，因此不涵盖驱动、显存、设备浮点、CUDA 并发或性能。HPC 测试模拟调度器，不代表真实集群提交。

## 文件与界面验证

所有交付 Bash 脚本执行 `bash -n`；Python 执行 AST 语法检查；R 文件使用 R 4.3.3 的 `parse()` 检查。Shiny 的动态参考排序、选择和分母计算另有小型实际 R 语义检查，但没有在完整依赖环境中交互启动整个 Shiny 页面。

两个包各自检查公开入口、实际命令生成、共享文件依赖和 SHA256 清单；HPC 包检查打包后的调度测试。科学核心文件在两包中逐字一致，差别是部署配置及 HPC 调度附加文件。顶层 `RELEASE.json` 记录版本、路径、默认参考及验证边界。

## 尚未执行

完整 1KG / UKB 重跑、真实 GPU 设备数值和性能测试、真实南科大 LSF 作业、全新 Conda/R 环境完整安装、完整 Shiny 交互测试。正式运行前，可先选 chr22、固定同一批样本及新结果根，确认本机或 HPC 实际程序、输入和资源需求。
