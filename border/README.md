# border: UK Biobank birth-order analysis

代码已从 `D:\Downloads\border\border` 核查并修订。正式输出默认写入 `/mnt/d/analysis/border`；输入只读，不修改 `all.rds` 或其他 UKB 数据。

## 正式运行

```bash
cd /mnt/d/scripts/border
./border.sh all
```

随访截止日直接读取 `/mnt/d/scripts/0f/phenotype.R` 第一行的 `date_follow_end <- as.Date("2023-04-01")`，无需额外参数；可用 `--follow-end YYYY-MM-DD` 覆盖。脚本仅解析这行日期，不执行共享文件的其他代码；日期和来源保存到队列元数据及日志。它用于缺少随访变量的结局（当前为 RA）。CAD 使用 `all.rds` 已有的 `cvd_cad.t2e`、`cvd_cad.Yt2e`；指定日期不会改写已有随访。上游截止日更新后，重新 `prep` 即可读取新值。

先单独检查输入和依赖，可运行：

```bash
./border.sh check
```

分阶段运行：

```bash
./border.sh prep
./border.sh phewas --phewas-cores 4
./border.sh assoc
./border.sh med --med-boot 100
./border.sh final
```

屏幕日志显示各阶段的开始/结束和耗时、每个 PheWAS 暴露的拟合/缓存复用、中介筛选数量及已完成候选数。中介进度最多每完成 5 个候选输出一次，不输出每次 bootstrap；详细模型日志仍保存在输出目录的 `logs/`。日志更新在新启动的任务中生效。

`all` 的顺序为 `prep → phewas → assoc → med → final`。PheWAS 是第一个分析阶段，完成后立即保存四张独立图，无需等待中介分析。

`prep` 保存已构建的结局，后续阶段无需重复指定随访日期。更新 `all.rds`、结局或截止日后须重新 `prep`，随后重跑下游阶段。`--in DIR` 指定包含 `Rdata/` 的输入目录，`--out DIR` 指定独立输出目录。`--outcomes cvd_cad,ra`、`--exposures border12,border1n,border_rank,border_relative` 可缩小分析范围。

## 暴露及模型

使用 `cnt_older`（p5057）、`cnt_brother`（p1873）、`cnt_sister`（p1883）以及已由上游生成的 `adopted`（p1767）。家庭儿童数 `sibship_size = cnt_brother + cnt_sister + 1`，包含本人。

| 暴露 | 定义 | 单位 |
|---|---|---|
| `border12` | 家庭恰有两个孩子；长子女=0、次子女=1 | 次子女 vs 长子女 |
| `border1n` | 家庭至少两个孩子；最大=0、最小=1，中间孩子缺失 | 最小 vs 最大 |
| `border_rank` | 非独生子女中 `cnt_older + 1` | 排行每增加一位 |
| `border_relative` | 非独生子女中 `cnt_older / (sibship_size - 1)` | 最小(1) vs 最大(0) |

排除 `adopted==1`，缺少整个 `adopted` 字段会报错。个人收养状态缺失者按原设计保留并在 QC 中计数。负数、非整数、无穷值 sibling counts 记为缺失；`cnt_older` 大于兄弟姐妹总数者排除。

主模型保持参考对话的 `age + sex + edu.sco + tdi + factor(center)`，不加入 birth year。`edu.sco` 必须是上游已派生的分数，不能用原始多选 `edu` 编码直接替代。扩展暴露不额外调整家庭大小，因此保留潜在 family-size confounding；这里没有实现家系内比较或因果识别。

- **CAD/RA**：Cox 回归，使用基线后的随访时间（年）；排除既往病例、无效时间及非 0/1 事件。需重建结局时优先使用 `fod_icd10_<outcome>`，并在失访、死亡和行政截止日中取最早日期。`outcome_qc.csv` 记录具体来源。
- **PheWAS**：对 `PheWAS.rds` 中的记录疾病状态做 logistic 回归，不是 incident Cox PheWAS。强制转换 0/1 为逻辑型，避免 PheWAS 包把整数结局当连续变量做线性回归。NA 保留，文件未涵盖的人不补成对照。原始 P 值保留，P=0 只在画图时作下限处理。BH FDR 和 Bonferroni 阈值按暴露分别计算。
- **中介分析**：protein/metabolite 分层先做 X→M（BH FDR<0.05），再做调整 X 后的 M→Y（P<0.05），只对通过两道筛选的候选计算线性×Cox 系数乘积及 bootstrap，不对全部蛋白/代谢物做中介分析。`FDR_b` 会保存供参考，第二道筛选实际使用未校正的 `p_b<0.05`。连续暴露和 mediator 在完整病例样本内标准化；二元暴露保持 0/1。X→M 不再按每个结局重复拟合。

中介结果为探索性关联分解。筛选后的 `FDR_indirect` 不校正先筛选再检验的影响；`p_boot` 是 bootstrap 符号尾部比例，不应当作验证性因果检验。默认 100 次的最小双侧 `p_boot` 为 2/101。需要更稳定的探索性区间时可提高 `--med-boot`；少于 50 次有效重复，或有效率不足 80%，不输出 bootstrap 区间/P 值。`--med-boot 0` 仅输出系数及 Sobel 近似，不生成 bootstrap 区间。

## 输出与恢复

- `raw/cohort_qc.csv`、`raw/outcome_qc.csv`：队列及结局核查。
- `raw/border.cohort.rds`：分析用队列与输入来源元数据。
- `raw/main_associations.csv`：所有请求组合，包括成功、样本量不足或拟合失败的状态与原因。
- `raw/phewas.csv`：PheWAS 结果；`phewas_<exposure>.rds` 为带输入/代码指纹的缓存。
- `raw/med_<layer>_<exposure>_Xa.csv`：完整 X→M 筛查，即使没有候选也保留。
- `raw/med_<layer>_<exposure>_<outcome>_Mb.csv`：M→Y 筛查。
- `raw/mediation.csv`、`raw/mediation_methods.txt`：候选中介结果与解释。
- `fig/Fig1_phewas_border12.png`：两孩家庭次子女 vs 长子女的 PheWAS。
- `fig/Fig1_phewas_border1n.png`：最小 vs 最大的 PheWAS。
- `fig/Fig1_phewas_border_rank.png`：出生排行的 PheWAS。
- `fig/Fig1_phewas_border_relative.png`：相对出生顺序的 PheWAS。
- `fig/Fig2_main_associations.png`、`fig/Fig3_mediation.png`：主疾病关联和中介分析图。

PheWAS 图参考 `D:\analysis\maha\ukb\Fig1.phewas.png`：每个暴露单独输出常规单面板 Manhattan 图，保持四图分类配色和横坐标一致。阈值以上的正 beta 为透明空心上三角，负 beta 为实心下三角；阈值以下仍为圆点。标题二元暴露显示 cases（暴露=1）vs controls（暴露=0），连续暴露显示 N，出生排行显示 N 和 ordinal 范围；N 为暴露、协变量完整且在 PheWAS 文件中的样本数，单个疾病的实际 N 可因缺失/对照排除不同；显示值在 -log10(P)=50 处截顶，原始 P 值不变。分类汇总保存为 `raw/phewas_category_summary.csv`；标题人数保存为 `raw/phewas_plot_samples.csv`。图中不显示底部形状图例和说明文字，这些说明由 figure note 提供。共享函数为 `/mnt/d/scripts/0f/plotting.R` 的 `plot_phewas_results()`，border、UKB 的分析/重绘入口、MAHA 均使用它。

四张 PheWAS 图在 `phewas` 阶段即时生成；`./border.sh plot` 仅从保存结果重画四张 PheWAS 图，`final` 也会重绘。若手动缩小 `--exposures`，仅生成所选暴露的图；某个暴露没有任何可估计结果时，会明确提示并跳过空图。

PheWAS 缓存会校验队列、疾病文件、模型代码和包版本；匹配时复用，`--replace` 可强制重算。其他阶段每次重算。不要对同一个输出目录并发启动多个运行。图表重绘前清除同名旧图，避免本次没有候选时留下旧中介图。

依赖：R、dplyr、tidyr、tibble、survival、PheWAS、ggplot2、ggrepel、digest、patchwork。缺失的 PheWAS/ggrepel 及必要依赖已安装在本机 R 用户库，不在脚本运行时自动安装。PheWAS 来源为 [PheWAS/PheWAS](https://github.com/PheWAS/PheWAS)。

## 开发验证

按[总目录约定](../README.md)，独立回归测试、验证脚本及模拟数据临时放在 `/tmp`，不保留在正式代码目录。正式运行前可用 `./border.sh check` 检查输入和依赖。
