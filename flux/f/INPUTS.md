# UKB 输入、时间轴与数据边界

本适配器直接读取原始 TSV/TSV.gz，按参与者和真实记录日期整理到 `data/cohort.sqlite`。它不读取已经在全人群上插补、标准化的 `Rdata/*.rds`。训练集拟合的词表、数值均值/标准差、特征清单、拆分信息和输入文件路径记录在 `data/meta.json`。标准化参数可随模型保存和复用。

默认 UKB 根目录沿用现有仓库的布局：本地 `/mnt/d/phe/ukb`，HPC `/work/sph-huangj/phe/ukb`。可以通过 `--dir-ukb` 或 `paths.ukb` 覆盖；所有 `data.sources` 相对路径都以该根目录解析。

## 1. 原始文件自动发现

| 配置键 `data.sources.*` | 自动发现路径 | 处理方式 |
|---|---|---|
| `participants` | `participants.tsv.gz`、`participants.tsv` | 如存在，优先采用下文规范参与者表 |
| `vip` | `rap/vip.tab.gz`、`rap/vip.tab` | 招募日期、出生年/月、性别、随访退出及可用访视日期 |
| `death` | `rap/death.tab.gz`、`rap/death.tab` | 字段 p40000；采用最早有效死亡日期 |
| `pca` | `rap/pca.tab.gz`、`rap/pca.tab` | PC1–PC40 等遗传主成分；作为静态背景 |
| `events` | `events.tsv.gz`、`events.tsv` | 如存在，优先采用规范有日期的 ICD-10 事件表 |
| `hesin` | `rap/raw/`、`raw/`、`rap/`、根目录下的 `hesin.tsv.gz`、`.tsv`、`.tab.gz`、`.tab` | 住院 episode 表 |
| `hesin_diag` | 同上，文件名 `hesin_diag` | 与 episode 表按 `(eid, ins_index)` 精确连接 |
| `icd10_codes` | `rap/icd10.code.tab.gz`、`rap/icd10.code.tab` | p41270 诊断数组，仅在显式允许后作为后备输入 |
| `icd10_dates` | `rap/icd10.date.tab.gz`、`rap/icd10.date.tab` | p41280 首次出现日期，与 p41270 按 eid 和数组下标匹配 |
| `observations` | `observations.tsv.gz`、`observations.tsv` | 如存在，优先采用规范有日期的连续测量长表 |
| `bbc` | `rap/bbc.tab.gz` | 原始血液生化/计数宽表 |
| `prot` | `rap/raw/prot.tab.gz` | 原始 NPX 宽表，保留合法负值 |
| `met` | `rap/met.tab.gz` | 原始代谢组宽表 |
| `img` | `rap/img.tab.gz` | 数值型影像衍生表型 IDP；不读取原始 MRI 影像 |
| `img_extra` | `rap/img.extra.tab.gz` | 补充影像访视日期和主表未包含的数值列 |

配置为 `auto` 或省略该键时执行上述查找。配置具体路径时文件必须存在；可使用绝对路径。`null`、`none`、`disabled` 可关闭该候选源，例如关闭 `prot`。关闭某个临床事件候选源会继续检查后备源，不能让整个临床事件输入为空。

`data.modalities` 控制原始宽表中的 `clinical`、`bbc`、`prot`、`met`、`img` 自动加载。规范 `observations.tsv` 已显式指定模态，因此其中的模态均会加载。可选宽表不存在时会跳过，并写入 metadata 的说明。若缺少参与者或有效临床事件，准备阶段会停止，不会生成虚构记录。

PCA 表只接受 `PC1`、`PC2` 等明确命名列，或 UKB 字段 p22009 的数组列（如 `p22009_a0` → PC1）。其中其它数值列不会被擅自改名为遗传主成分。原始 VIP 中只有性别和明确 PC 列进入静态背景，评估年龄、测量值和访视日期不会成为永久背景。规范参与者表的 `context_*` 是使用者显式声明的静态输入，不能用它存放访视后才形成的测量或结局。

## 2. 临床事件：完整住院记录与首次出现摘要

优先级为：规范 `events` → 完整 `hesin` + `hesin_diag` → 明确开启的 p41270/p41280 后备源。

完整住院记录使用 `hesin_diag.diag_icd10`，通过 `(eid, ins_index)` 找到住院 episode。日期优先使用有效的 `epistart`，缺失时使用 `admidate`。本程序不会用出院日期或 episode 结束日期冒充诊断开始日期。缺失连接、缺失有效开始日期的记录会被跳过并计数。同一参与者在不同日期的重复 ICD-10 代码会保留；同一日、同一 ICD-10 三位码去重。

**p41270/p41280 是诊断首次出现的摘要，并不是完整的重复就诊流水。** 当 `data.allow_first_occurrence: true` 时可用它们运行有限数据版本，但 `meta.json.event_source_mode` 会明确记录 `first_occurrence_only`，终端和报告也会注明。仅凭这两个文件不能声称复现了完整 EHR 轨迹。数组可以由多个 `_aN` 列组成；代码列也支持 `|` 或 `;` 分隔的串。日期通过 eid 与同一数组位置匹配，不依赖两个文件的行顺序。

当前模型统一使用 ICD-10 三位码，例如 `I25.1` → `I25`。代码和 endpoint 配置均按此分辨率解释。原始 ICD-9、GP 的 Read/SNOMED、OPCS、处方、原始影像像素没有被自动转换。若需要纳入经核实的历史 ICD-9/GP 映射，请先生成规范的 ICD-10 事件表，并保留映射版本及日期依据；不要把未核对的代码系统直接混入 ICD-10。

死亡来自参与者表或原始 p40000 导出，并作为 `DEATH` 事件加入。死亡后的临床和测量记录会被排除。

## 3. 日期和数值测量

支持 RAP 列名 `p53_i0`、`p4080_i0_a1`、`p41280_a0`，也支持旧导出格式 `53-0.0`、`f.53.0.0`。日期接受 ISO 格式 `YYYY-MM-DD` 及以该格式开头的时间戳；不猜测 `01/02/2024` 之类的区域格式。

无效日期和 UKB 哨兵日期 `1900-01-01`、`1901-01-01`、`1902-02-02`、`1903-03-03`、`1909-09-09`、`2037-07-07` 均按缺失处理。出生日前、行政随访结束后、死亡后的记录不进入模型。原始出生年月用当月 **15 日**生成近似出生日期；这是原仓库采用的近似，metadata 记录其人数，模型年龄因此有相应日期精度限制。

同一测量列中的 `_iN` 与 p53 的第 N 次真实评估日期匹配。没有相应访视日期的测量不参与训练；计数存入 `summary`。影像的 i2/i3 日期可从 `img` 或 `img_extra` 补充，不能倒填到基线。对于独立的日期列，可设置 `data.date_columns: {prot: sample_date}` 等显式映射；也自动识别 `date`、`sample_date`、`measurement_date`。

当前用户仓库的未标注日期蛋白宽表需要显式配置：

```yaml
data:
  prot_date: baseline
```

这只适用于确认全部为基线采样的 NPX 文件。程序会在 metadata 中记录该假设。若包含复测样本，必须提供真实采样日期或改用规范长表。NPX 的负数是合法值，不能作为 UKB 缺失编码清除。临床、血液和代谢宽表中的负数缺失编码会过滤；其它模态不统一删除负数。

临床宽表默认选择以下 UKB 字段，均按所属评估 instance 定位日期：

```yaml
data:
  clinical_fields:
    - p21001
    - p21002
    - p48
    - p49
    - p50
    - p46
    - p47
    - p4080
    - p4079
    - p102
    - p874
    - p884
    - p904
```

可以替换此列表，但应选择真正的连续型变量。程序不会自动把生活方式类别编码解释为连续量，也不会推断分类词典。数值特征的固定训练清单、均值、标准差和所属模态均由训练参与者决定；验证/测试中独有的特征不会新增模型输入。`min_feature_count` 约束训练测量数；`max_missing_fraction` 的分母是**训练集中有该模态的参与者数**，不是全 UKB 人数。缺失保持为掩码，供模型识别，未在全队列上插补。

原始 BBC 和代谢宽表若存在本地 `common/bbc.lst`、`common/met.lst`，默认使用其第一列的独立 UKB 字段 ID 作为白名单。这与现有 `phe.R` 的无表头、制表符分隔、前两列分别为 ID 与名称的字典格式一致。第二列明确标记为 date、assay、batch、QC 等的行被排除；`p23400/p23401` 等计算表达式不会当作原始字段，也不会在这里自动生成比值。字典作为输入源记录路径、大小和修改时间。若字典存在但没有可识别的独立字段，程序报错；不会猜测其格式。若字典不存在，metadata 会提醒检查可能被数字字段 ID 隐藏的实验技术变量。

可用以下配置明确控制**原始宽表**的特征：

```yaml
data:
  feature_include:
    bbc: [p30600, p30700]
    prot: [GDF15, APOE]
  feature_exclude:
    bbc: [p30700_i2]
```

选择器采用区分大小写的精确匹配，支持原始列名、标准化后的特征名、UKB 基础字段 ID，或 `模态:特征`。例如 `p30700` 匹配这个字段的所有 instance/array，`p30700_i2` 只匹配这一原始列。它们不是正则表达式或 glob。某模态配置了 `feature_include` 后，以显式列表取代本地字典白名单；空列表表示该模态不取任何数值特征。`feature_exclude` 最后生效并优先于 include。日期、ID 等已知元数据列仍不会成为数值特征。规范长表已由使用者明确提供特征和值，不再套用原始宽表白名单。

特征统计会使用训练参与者的所有有效日期测量，供共享的自监督模型拟合。风险预测时另行按 landmark 截断每位参与者的输入；验证和测试参与者的任何测量都不参与词表/标准化/缺失率筛选的拟合。元数据的 `numeric_scaler` 和 `context_scaler` 可用于审查这一边界。

## 4. 规范 TSV 输入接口

如果原始导出结构与自动发现路径不一致，可显式指定三张表：

```yaml
data:
  sources:
    participants: /absolute/path/participants.tsv.gz
    events: /absolute/path/events.tsv.gz
    observations: /absolute/path/observations.tsv.gz
```

所有文件均为制表符分隔，UTF-8，第一行是列名；可 gzip 压缩。缺失值可用空串、`NA`、`NaN`、`null`。不会把 `eid` 转为浮点数。未知 eid 的事件和测量会跳过并计数。

### participants.tsv

每个参与者一行，eid 不得重复。

| 列 | 必需 | 含义 |
|---|---|---|
| `eid` | 是 | 稳定的参与者标识；字符串 |
| `birth_date` | 是 | 出生日期或已明确声明的出生日期代理 |
| `baseline_date` | 是 | 招募/基线评估日期 |
| `exit_date` | Optional | Verified earlier individual exit or source-coverage limit; cannot extend the canonical administrative limit |
| `death_date` | 否 | 已知死亡日期，不用它替代行政随访结束 |
| `sex` | 否 | UKB 0/1；其它值按缺失 |
| `PC1`、`PC2` 等 | 否 | 数值遗传主成分 |
| `context_*` | 否 | 经确认不包含未来结局的静态数值背景 |
| `split` | 否 | 默认忽略；仅 `data.use_input_split: true` 时使用 `train`、`validation`、`test` |

默认按 eid 和 seed 稳定哈希拆分，比例由 `data.split: [0.7, 0.15, 0.15]` 指定。各数据模态共用同一拆分，不会让同一人进入多个集合。出生日期也可由 `birth_year` + `birth_month` 提供，使用上述月中代理。

The administrative maximum is resolved at runtime from `0f/phenotype.R::date_follow_end`; YAML and CLI date overrides are not supported. Verified earlier individual exit, loss to follow-up or source-coverage dates still shorten follow-up. The last diagnosis or measurement date is not a censoring date. Existing death-event and competing-event rules are unchanged.

### events.tsv

| 列 | 必需 | 含义 |
|---|---|---|
| `eid` | 是 | 参与者标识 |
| `date` | 是 | 实际事件/诊断记录日期 |
| `code` | 是 | ICD-10 代码或 `DEATH`；统一到三位码 |

表可以不排序。每行是一个有日期的代码，同日多个不同代码可占多行。事件来源优先级是互斥选择；不会自动把 canonical 与 HES 再合并，避免重复和未声明的来源改变。

### observations.tsv

| 列 | 必需 | 含义 |
|---|---|---|
| `eid` | 是 | 参与者标识 |
| `date` | 是 | 实际采样/测量日期 |
| `modality` | 是 | 如 `clinical`、`bbc`、`prot`、`met`、`img` |
| `feature` | 是 | 同一模态内稳定的特征名 |
| `value` | 是 | 有限数值；缺失测量可省略该行 |

任意顺序均可；SQLite 外部排序后归并同人、同日、同模态的测量。同一 `(eid,date,modality,feature)` 重复时采用后写值，因此应在上游先明确技术复测的合并规则。全队列不构建“参与者 × 时间 × 全特征”的巨型稠密数组；只在读取单人窗口时构造该窗口的数值张量。

## 5. landmark、历史前缀与事件标签

`data.landmark: midpoint` 是面向论文比较的默认设置。先取该人所有非死亡临床事件的**不同日期**并排序；以零起始下标 `min(n // 2, n - 2)` 选日期作为 landmark。`data.min_clinical_days: 3` 要求至少三个不同临床日期；最后一次临床日期不能作为中点。选择使用的是整个已观测轨迹，这是回顾性基准的取样定义，应与前瞻性部署时点区分。

所有同日事件和测量都纳入截至 landmark 的历史前缀。风险推断保留完整前缀，不受预训练的 `max_nodes: 96` 窗口上限限制。模型还接收明确的 `query_age`，从最近观测状态推进到 landmark；不能简单拿最近一次测量年龄替代预测年龄。风险输入的 `next_event`、`next_dt`、`label_mask` 全部置零，结局只在独立的评估字段中保存。

其它可选取样方案：

| `landmark` | 额外参数 | 含义 |
|---|---|---|
| `baseline` | 无 | 招募日期；用于易解释的基线预测复现变体 |
| `age` | `landmark_age: 60` | 固定达到年龄 |
| `date` | `landmark_date: YYYY-MM-DD` | 固定日历日期 |
| `last_event` | 无 | 基线当日或以前的最后一个临床事件日；没有事件则使用基线 |

`min_clinical_days` 对这些风险取样方案同样生效；若研究目标是纯基线预测而不要求 EHR 天数，可显式设为 0。landmark 必须早于行政随访结束，且参与者在该日存活。

预训练从连续时间轴抽取至多 96 个节点；训练窗口的抽样可由 seed、epoch 和 eid 重现。验证/测试使用固定末尾窗口。一天内先写入该日数值测量节点，再按字典顺序写入 ICD-10 节点；此排序不代表已知临床先后。目标是**严格晚于当前日**的下一个临床日期的首个字典序代码，不预测同一天内的伪造时间顺序。这是当前独立实现的明确建模选择，并不声称恢复了原论文未公开的内部编码细节。

每个 endpoint 根据完整存储历史寻找第一次日期。landmark 当日或以前已有该病者不进入该病的 incident-risk 分析，即使该疾病发生在很久以前也不会因窗口截短而被当成“新发”。`endpoints: all` 使用训练词表中的三位疾病代码以及 `DEATH`；也可指定 `{t2d: [E11], cad: [I20, I21, I22, I23, I24, I25]}`。死亡作为非死亡疾病的竞争事件，行政随访与死亡日期分开保存；评估阶段应依据其评估口径处理不足五年的未发病者，不得一律当作有完整五年随访的阴性对照。

## 6. 合成演示与本地结果审查

`demo` 生成的是程序合成的 320 人 TSV，并显式标记为 `SYNTHETIC.json`。其中有先后发生的多模态测量、同日多代码、既往病、死亡和不完整随访，并人工保证小样本拆分中有死亡阳性，以检查执行路径。其 AUROC、AUPRC 等仅证明代码能运行，不能与论文性能比较。

相同 seed、人数和生成器版本会复用已有合成输入，保持文件修改时间稳定；参数或生成器版本改变时重新生成。真实数据准备后可审查 `meta.json` 中的 `sources`、`event_source_mode`、`summary`、`warnings`、`split_counts`、`preprocessing_fit_split`、特征统计和行政随访日期。输入目录、模型和输出中包含受控健康数据，应按已有 UKB/机构环境要求保留在获批本地或 HPC 位置。
