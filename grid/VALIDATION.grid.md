# GRID individual evolution 本地核查与整合

日期：2026-10-07。来源为 `D:\Downloads\grid_individual_evolution\grid_individual_evolution`，结合原 ChatGPT 对话最新要求核查。当前目录不是 Git checkout，因此依据 `GRID.patch` 的改动清单逐项合并；核对了补丁内 9 个新增文件与包内源码完全一致。

## 实际合并范围

- `grid.sh` 改为新方法入口，新增 `f/grid.py`、`f/grid.data.py`、`f/grid.evolution.py`、`f/grid.abm.py`。
- 更新 CPU 环境、安装器、共享帮助入口；PRSformer 数据导出增加用于共同测试集比较的来源字段。
- 保留 `0.pca.sh`、`1.csx.sh`、`2.disco.sh` 及当前辅助文件命名，修正新代码中的旧调用名称；保留本地 CSx 运行锁改动及现有 PCA、Disco、Yeval 实现。
- 旧 `f/3grid.f.sh`、`f/3grid.py` 退出运行目录；保留 GU 所需的 `f/0.arg.py`。完整旧目录在应用前另行备份。
- 原有 `1csx.*`、`2disco.*` 结果文件名继续使用；没有移动或覆盖已有 UKB 分析结果。

应用前完整代码、Conda 环境导出、显式包清单及 pip 清单保存在 `/mnt/d/Downloads/grid_individual_evolution/grid_before_integration_20261007_174931`。该目录也保存安装和验证日志。

## 核查发现及修复

1. **Python 依赖冲突。** 包内 `environment.yml` 固定 Python 3.10，但 `rdata==1.1.0` 要求 Python ≥3.11。实际安装复现了失败；较旧的 rdata 0.11.2 没有本流程所需的 `write_rds`，不能直接降级替代。CPU 环境改为 Python 3.11，入口在旧解释器下给出明确错误。安装器默认更新命名环境 `grid`。
2. **缺失表型误触发 50/50 错误。** 原实现对单个 trait 删除缺失 Y 后再次强制 50/50，可能拒绝已经正确冻结的共同名单。现在只在共同 roster 阶段验证 50/50；各 trait 保留原归属，按实际可分析人数拟合和评价。
3. **缓存未跟踪名单文件。** 原始 `keep`、`remove`、`split_file`、`split_group_file` 未进入输入身份清单。现在纳入来源检查，改变撤回或家系等名单后拒绝继续使用旧缓存；允许合法的空撤回文件。
4. **原生输入预检漏验模型选项。** 原生 `--check` 只验证默认 seed 配置，非法的 `k_grid` 等参数可能延迟到评分后才报错。现在在无需已构建注释特征的情况下核验实际模型选项。

## 本次实际执行的验证

测试文件：`f/test_grid_integration.py`，按用户要求集中放入 `f/`。它生成可独立计算预期结果的模拟样本，用来检查程序逻辑，不读取真实 UKB 数据，也不用于评价真实队列的预测性能。真实输入预检使用 `bash grid.sh --traits height,ldl,t2dm --check`；正式评分、拟合和报告阶段还会执行相应的一致性检查。

先用隔离临时 Python 3.12.3 环境（rdata 1.1.0、NumPy 2.5.3、pandas 3.0.6、scikit-learn 1.9.1）完成 10 项测试，再依用户要求更新已有命名环境 `grid`，在项目目录和最终环境中重新执行，10 项全部通过（19.105 秒）。

最终环境为 `/home/huangj/miniforge3/envs/grid`：Python 3.11.6、rdata 1.1.0、NumPy 1.26.4、pandas 2.3.3、SciPy 1.17.1、scikit-learn 1.9.0、Matplotlib 3.9.1、openpyxl 3.1.5、PLINK 2 a6.9、R 4.2.3。已从用户现有 Anaconda 初始化脚本实际执行 `conda activate grid`，确认名称能正确定位该环境。

更新依据项目 `environment.yml`；完整索引下载较慢，实际安装改用 Mamba 的本地有效索引缓存，再单独执行安装器中的 `pip install rdata==1.1.0` 步骤。未修改全局 Conda 配置。`pip check` 无依赖冲突，GRID/CSx Python 帮助入口可用，共享流程所需的 15 个 R 包均能加载。升级前后环境导出、包清单、运行版本及两轮测试日志均保存在上述备份目录。Matplotlib 的第三方 pyparsing 弃用提示不影响本次图表生成。

生成 602 人、两个染色体、24 SNP 的 BED 文件，再用真实 PLINK 转换为 PGEN。排除一个负 ID 和一个撤回 ID 后，共同名单为 600 人、300 个家系，训练和测试各 300 人。加入跨染色体样本乱序、缺失剂量、效应等位基因反向，以及 80 个仅在 training half 缺失的 LDL 标签。LDL 仍使用同一外层名单，其可分析数据为 220 train / 300 test。

三个性状均完成 prepare → evolution → fit → report，10 项回归测试全部通过：

- 训练频率中心化的四套总分与独立 NumPy 真值一致，模块分数能重构总分。
- 共同名单不变，家系没有跨 train/test，缺失 LDL 不重分组。
- unknown/uncertain 年龄标签及其 SNP 位置在置换对照中固定，各层标签计数不变。
- 撤回、保留、家系、外层拆分文件变更均触发缓存拒绝；空撤回表可记录。
- 原生预检能拒绝非法模型参数。
- 修改全部测试集 height Y 后，重新拟合的预测和筛选列均保持一致；donor OOF 审计中家系不重叠。
- 三个正式模型 RDS 重载后，对不含 Y 的测试输入重现预测。
- 匹配核强制非零借用时，参考人贡献加基线及截断修正能重构预测，并排除同家系 donor。
- 生成的 PNG 均有可读同名 XLSX，汇总工作簿不含个体 ID；配对损失差与独立公式一致。
- 同大小评分权重篡改被拒绝，PCA/CSx/Disco 转发使用点号命名。

另通过 Python/Bash 语法、ShellCheck 和 GRID＋PRSformer 编排 dry-run，含空格路径及强制共同 keep/family/split 参数。合成数据和结果仅在临时目录产生，并在测试结束后清理。

复跑：

```bash
conda activate grid
cd /mnt/d/scripts/grid
python f/test_grid_integration.py
```

解释器需要安装 GRID 的 CPU 依赖；PLINK 默认从 PATH 或 `~/miniforge3/envs/grid/bin/plink2` 选择，也可设置 `GRID_TEST_PLINK2`。

## 验证边界

本轮没有下载完整 GEVA、运行真实 UKB 分析、重新训练 PRSformer 或验证 HNSW/GPU 全基因组路径。原包 README 中 PRSformer CPU 训练及更多指标复算属于交付方原有记录。当前 `grid` CPU 环境已就绪，三个性状的既有 CSx 分数均存在；正式运行仍需有来源的年龄注释（默认 GEVA 目录尚不存在）。联合 PRSformer 比较还需要其独立环境及官方源码，当前默认位置尚未配置。

这些检查证明所覆盖的程序行为，不能证明 GRID 的预测性能优于 CSx/PRSformer；发现 GWAS 与测试样本的独立性、真实年龄注释覆盖和外部验证仍需在实际分析中核实。
