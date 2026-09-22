# Panome

Panome把人的分子组合、邻域和混合状态作为分析对象，检验它们能否在常规组学预测之外解释个体差异。本次下载包未附带其原先引用的 PROJECT_PROPOSAL.md、DESIGN.md 和 VALIDATION.md；本地核查与迁移记录见 `LOCAL_REVIEW.md`。当前版本为3.0.0；它是可运行的研究实现，不包含真实UKB结果。

## 文件与运行方式

沿用LE8的组织习惯：根目录保留 `panome.sh` 入口，计算代码放在 `f/`，结果放在独立analysis目录。Python负责建模，`f/export_rds.R`仅负责读取RDS并导出所需字段，不执行LE8的表型构建或遗传分析。

`f/panome.py`调度六个阶段；`io_data.py`读取数据；`preprocess.py`处理结局、分组和训练集预处理；`representation.py`学习AE与分子模块Transformer；`graph.py`建立训练参照图；`prediction.py`、`neural.py`和`evaluation.py`拟合并比较预测；`pgs.py`提供可选遗传评分注释；`interpretation.py`解释个体；`final.py`生成图及源数据；`project.py`投影新个体；`common.py`记录配置、缓存与锁。

建议Python 3.11或3.12。在本目录创建环境，或使用已有环境并设置 `PANOME_PYTHON`：

```bash
python3 -m venv .venv
.venv/bin/python -m pip install -r requirements.txt
chmod +x panome.sh
./panome.sh --help
```

RDS优先使用可用的 `Rscript`，否则使用pyreadr。含特殊R对象的文件建议使用Rscript；RDS必须是一个data.frame。CSV、TSV和gzip文本也可直接读取。Parquet需要额外安装pyarrow。CUDA须安装与机器匹配的PyTorch；`--device auto`自动选用可用GPU，`--device cpu`强制CPU。不会自动安装R包或修改输入数据。

## 数据位置

默认 `UKB_PHE=/mnt/d/data/ukb/phe`，也可通过 `--ukb-phe` 指定。以下位置来自核对后的LE8/UKB代码：

- 表型：`/mnt/d/data/ukb/phe/Rdata/all.rds`。
- 原始蛋白：`/mnt/d/data/ukb/phe/rap/raw/prot.tab.gz`。
- 原始代谢物：`/mnt/d/data/ukb/phe/rap/met.tab.gz`。
- 代谢物映射：`/mnt/d/data/ukb/phe/common/met.lst`，支持无表头或当前 data_field/met_name 表头，前两列为原始字段或算术表达式、目标特征名。保留baseline `_i0`，去掉后续访视字段，只允许字段之间的加减乘除。
- 可选已清洗矩阵：`Rdata/prot.rds`、`Rdata/met.rds`，用 `--input-source cleaned` 明确选择。上游全样本过滤或插补无法在此撤销，主分析建议原始数据。已变换、含负值的代谢物需指定 `--transform none`。
- 可选PGS：优先 `Rdata/prot.pgs.rds` 或 `Rdata/met.pgs.rds`，不存在时读取 `all.rds` 中的 `*.pgs`；`--pgs-file`覆盖这一选择。

矩阵中第一列应为唯一 `eid`，其余列只放分子测量。结局和协变量放在表型文件中；已知结局、随访或协变量混入分子矩阵会报错。蛋白名转为大写，代谢物保留映射后的名称。不固定蛋白数量。输入均要求一个人一行；本版不自行合并重复访视。

默认临床协变量为 `age,sex,tdi,PC1,PC2,center`，其中年龄采用训练拟合的样条。分子残差化默认 `age,sex,prot.plate` 或 `age,sex,met.plate`。这些字段必须实际存在；可用 `--covariates`、`--residualize`和 `--categorical`明确修改。临床变量可以加入已核验的baseline风险因子或LE8变量。本版默认协变量集合不是SCORE2、QRISK等经过验证的临床风险计算器。

默认不限制祖源，也不假定所有人无其他疾病。正式研究请准备符合纳入标准的表型文件；需要额外排除baseline疾病时，使用 `--healthy-date-cols` 提供一组首诊日期列。`--group-col family_component`使用完整亲缘连通分量分组，不能只给某个家庭成员的零散亲属编号。

## CAD与height

先查看路径、再检查文件和依赖：

```bash
./panome.sh --Y cvd_cad --biom prot --dry-run
./panome.sh --Y cvd_cad --biom prot --preflight
./panome.sh --Y cvd_cad --biom prot --run-name cad_primary
```

CAD默认读取 `fod_icd10_cvd_cad`、`date_attend`、`date_death`和 `date_lost`。末次行政随访默认2023-04-01，来自现有UKB表型代码；正式运行应按自己的数据覆盖日期指定 `--end-date`。baseline当天及之前确诊者排除，事件不得晚于死亡、失访或行政截止时间。若另有确诊标记而缺少首诊日期，可通过 `--disease-evidence-col`排除此类时间不明者。仅有日期列时，程序无法识别未编码的既往病例。

```bash
./panome.sh --Y height --outcome-type quantitative --target-col height \
  --biom prot --run-name height_primary
./panome.sh --Y cvd_cad --biom prot,met --run-name cad_layers
```

`height`须替换成all.rds中的实际列名；本环境未读取你的UKB数据来核对该列。`prot,met`执行两个独立图谱，**不是**多视图联合训练。当前一次接受一个Y；分析多个疾病可在shell中循环调用。

启用已有PGS的可选注释：

```bash
./panome.sh --Y cvd_cad --biom prot --run-name cad_pgs \
  --pgs --pgs-file /mnt/d/data/ukb/phe/Rdata/all.rds \
  --pgs-source 'describe GWAS and scoring version' --pgs-overlap unknown
```

只有 `feature.pgs` 能匹配该feature，蛋白匹配忽略大小写。PGS模块比较临床项与临床项+PGS对成年实测分子的重构，输出独立测试R²增量和个体偏离；它不参与主图谱或主疾病预测。GWAS发现样本的重叠不能通过本项目的60/20/20划分消除，需记录实际来源。PGS不是出生时的蛋白浓度，偏离不能直接解释为环境效应、亚临床疾病或因果成分。

## 阶段与结果

默认输出 `/mnt/d/analysis/panome/<Y>/<biom>/<run-name>/`；`--analysis-root`可修改根目录。默认run-name为 `v3`，不会导入旧版缓存。

- `s1_prepare`：保留原始缺失的分子矩阵、表型与输入审计。
- `s2_preprocess`：结局、60/20/20划分、仅训练拟合的处理器及QC审计。
- `s3_representation`：AE、PCA和可选模块Transformer的冻结模型与坐标。
- `s4_graph`：训练样本图、Leiden发现标签、所有人的投影坐标、状态权重和稳定性。
- `s5_predict`：所有可用模型、测试预测、校准、配对bootstrap、landmark与PGS结果。
- `s6_report`：个体分子状态、状态关联、AE归因、近邻注意力、相近风险配对及文字报告。
- `publication`：Fig1图谱、Fig2预测、Fig3个体；每组均有PNG、PDF和XLSX源数据。

本版增加状态权重调节的变系数模型、邻域汇总模型和对训练参照的cross-attention，并保留临床、PWAS score、PCA、AE、elastic net、树模型、个体内部Transformer及同结构随机邻域等对照。最多16个模型；无合格多状态解、事件不足或模型失败时会明确记录，不伪造输出。务必查看 `s5_predict/model_status.csv` 与 `metric_limitations.json`，不要只看总流程DONE。

部分运行和重绘示例：

```bash
./panome.sh --Y cvd_cad --biom prot --run-name cad_primary --to s4_graph
./panome.sh --Y cvd_cad --biom prot --run-name cad_primary --from s5_predict
./panome.sh final --Y cvd_cad --biom prot --run-name cad_primary
```

恢复运行须提供与原运行一致的分析参数。manifest记录源代码、配置、依赖版本和输入指纹；已完成阶段还核对输出指纹。大文件默认检查路径、大小和mtime；`--full-input-hash`会额外完整哈希输入。代码移动、参数或输入变化时使用新的run-name，或者明确 `--replace` 重建该run的全部阶段。`--replace` 会删除这个run下六个阶段与publication目录，保留其他run。中断的阶段会重做；本版不从半轮神经网络训练恢复。若硬中断留下锁，核实无运行进程后删除该run的 `.lock`。`final`从完成的结果重绘，不需要原始UKB文件，也不重新训练。

## 新个体投影

投影只需分子测量和baseline协变量，不需要Y。使用训练时相同的特征定义、单位、采样及平台处理。缺少训练保留的assay、样本缺失率过高或重复ID会报错；新出现的分类变量水平会记录审计，不能据此认为平台差异已经解决。

```bash
.venv/bin/python f/project.py \
  --run-dir /mnt/d/analysis/panome/cvd_cad/prot/cad_primary \
  --phenotype /path/new_baseline.csv --omics /path/new_proteins.csv \
  --output /path/new_people_panome.csv
```

若训练使用亲缘组，投影表型也需完整组编号。输入可以是CSV、TSV或RDS；原始代谢物另加 `--met-input raw`，并保留训练配置中的映射文件。近邻仅来自冻结训练参照，自身和同亲缘组排除。改变测试批量大小或顺序不会改变个体的参考人群。生存输出是死亡删失下的net risk，不是竞争风险累计发生率。

## 验证及研究边界

原始发布说明称其合成数据已经覆盖生存和连续结局、全部16个模型、新个体无Y投影、批量与顺序不变性、亲缘划分及部分RDS/代谢物接口。上述属于原始发布说明；本机现已完成真实 PPP 抽样读取、预处理与原生 Rscript 桥接核验，详见 LOCAL_REVIEW.md。仍未完成真实全量训练、GPU及全量资源验证。

正式分析应补齐强临床基准、plate-only和缺失模式对照、多次完整重训、预先定义的疾病排除、PPP选择性抽样评估、外部队列与平台校准。并行线程可用 `--cores` 调整；它不是总内存上限。可以先用 `--max-samples 5000 --run-name pilot`检查数据，但pilot不能作为全队列结果。

基于同一批分子建图后再发现组间分子差异，属于描述而非独立机制证据。注意力和integrated gradients解释模型，不识别因果。单次baseline数据不建立真实病程；本版没有执行MR、共定位、Dandelion、通路富集或动态agent模拟，也不声称优于现有模型。

重写依据：`jielab/scripts`，核对快照 `8c49af6a148e2cfd2f077176862f9d40c1ee72e7` 中的panome、le8、UKB表型函数，以及提供的DESIGN.md。迁移时建议把原panome目录另存后，使用本包的完整目录，避免旧辅助脚本混用。
