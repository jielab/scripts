# LE8

## 2026-10-06 最新版与完整重跑

当前默认是 `reference / selective_attention / cuda`，版本 `7.1.0-cuda-parallel-qc-20261006`。Transformer 训练及预测必须使用 CUDA；默认检索也在 CUDA，失败直接报错。当前实现和测试以 [FINAL_20261006.md](FINAL_20261006.md) 为准。下文旧日期的修复记录是历史记录；其中提及的临时结果已按用户要求清理，不代表当前还有这些分析结果。

可直接复制（先预检，成功后完整重跑 C1–C5）：

```bash
cd /mnt/d/scripts/le8
./le8.sh --Y cvd_cad --biom prot,met --preflight && \
./le8.sh --Y cvd_cad --biom prot,met --replace TRUE
```

全部分析完成后：

```bash
cd /mnt/d/scripts/le8
./le8.sh final --Y cvd_cad --biom prot,met
./le8.sh shiny --Y cvd_cad --biom prot,met --no-reindex
```

`--replace TRUE` 强制重新拟合。中断恢复时使用同样命令但去掉该参数。每个模块成功后立即发布已有结果；后续模块失败不会阻止前面已完成结果发布。定时更新已停用，代码不会再次被定时覆盖。本轮没有运行新的 UKB 全量分析，也没有恢复用户删除的正式结果。

统一入口为 `le8.sh`。`f/` 中公用文件以 `0.` 开头，分析文件以 `c1.` 至 `c5.` 标明归属。全项目整理约定统一见 [根目录 README](../README.md)。

| 模块 | 文件 |
|---|---|
| 公用、调度、索引、验证和导出缓存 | `0.common.R`、`0.common.py`、`0.engine.sh` |
| C1 相关与 ABM | `c1.correlate.R`、`c1.abm.py` |
| C2 因果分析、MR-link-2、cis 索引 | `c2.cause.R`、`c2.cause.sh`、`c2.mr_link2.py` |
| C3 共定位 | `c3.coloc.R`、`c3.coloc_GPU.py`、`c3.coloc_GPU.sh` |
| C4 连接分析 | `c4.connect.R` |
| C5 细胞分析 | `c5.cellulation.py` |
| 最终汇总、表格、科学问题报告 | `final.R`、`final.py` |

C1 使用 ABM（agent-based modeling）命名。模型运行、注释及图形代码集中在 `c1.abm.py`，注释和重绘分别使用 `annotations`、`figures` 子命令。新模型统一使用 `c1_abm` 序列化身份；不注册已经合并掉的旧 Python 模块名。

`final.py` 合并原先分散的最终报告、表格和科学问题处理。Shiny 入口为 `shiny/app.R`。正式执行统一使用 `le8.sh`；辅助文件是运行目录内的内部步骤。

## 2026-10-05 方法更新

- C1 保留常规关联及 PGS 分析，默认 reference ABM 使用 `selective_attention`：真实 Transformer 与个体检索模型在 CUDA 上训练和预测，开发集内选择、校准和冻结门控，验证集评估覆盖率与风险增益。默认 C1–C5 流程包含该 ABM，无需额外参数；结果独立保存在 `c1_correlate/abm_selective_attention`。TabICLv2 仍可通过 `--abm-backend tabicl` 或 `both` 选择。已有匹配结果按缓存规则复用，旧 CPU selective 结果保留在原目录。
- C2 将交叉拟合的遗传预测部分 G 与残差 R 同时放入调整模型，直接检验二者差异；残差不解释为纯环境作用。MR 使用边际效应、独立工具和坐标/等位基因身份；冲突重复变异排除并记录质量控制。MR-link-2 按完整预定检验家族校正。DANDELION 原生结果、全局多重校正、逐位点剔除结果分别展示；真实 WES 与 GWAS 基因适配证据分开。冻结的 Yang 状态投影及近远期风险属于补充分析。
- C3 固定同一组、有序且方向一致的变异供 CPU/GPU 比较，稳定计算 H3，并保留先验敏感性与输入条件。未通过敏感性检查的 H4 不标记为稳健共定位。SuSiE 的可信集与各 MR 工具逐一对应，只有部分工具信号获得支持时标记 partial_signal_support，不升级整个 MR 汇总结果。SuSiE 可选分支要求有符号 LD、参考样本量、版本和祖源元数据；缺失条件时给出不可用状态。
- C4 分开基本调整的总连接与条件特异性，允许一项分子连接多个 LE8 域。YS 按无疾病结局参与的连接强度和冗余排序；固定预算下比较 YS、NS、YSplus、原始分子及嵌套交叉拟合的概念表示。保存临床与分子线性预测值贡献、惩罚路径和非零系数，避免把零分子贡献解释为模型优势。路径乘积及其有符号比例只作描述；bootstrap 不足时不输出推断。模块稳定性是固定候选分子集合下的重聚类稳定性。
- C5 以全部已测分子为背景，保留未知注释，使用等预算且排他基因的细胞对比。原生 CIGMA、外部导入结果和不可用状态严格区分；只有符合输入条件时才运行对应分析，不从缺失证据推导衰老或细胞谱系结论。
- Final/Shiny 增加上述方法、有效性状态、验证比较及下载表。缺失分析明确显示缺失，汇总报告不把不可用结果补成阳性证据。

有家系或相关个体时，可用 `LE8_GROUP_FILE` 指定含 `eid,group` 的文件，或用 `LE8_GROUP_COLUMN` 指定个体表中的分组列；交叉验证、开发/验证划分和相关 bootstrap 按组进行。未提供时以个体为组，程序不会猜测家系。疾病 PRS 是可选临床基线输入，需显式提供 `C4_DISEASE_PRS_COLUMN` 或 `C4_DISEASE_PRS_FILE`，同时设置 `C4_DISEASE_PRS_TRAIT` 与 `--Y` 一致；不会拿单分子 PGS 替代疾病 PRS。

分析缓存检查代码、输入文件及相关设置；输入或方法变化时重新计算。同一结果根目录下，不同 `Y/biom` 的命令可以同时计算；相同 `Y/biom` 的命令排队，避免重写同一组结果。读取结果快照和发布结果时使用短期共享目录锁，每个分析任务只发布自己的疾病/数据层。`final,shiny` 等汇总任务等待正在运行的分析结束，再合并已有疾病/数据层；已有 Shiny 服务时不重复启动。测试请使用独立的 `--analysis-root`，不要将缩小特征数或 bootstrap 次数的测试结果发布到正式目录。

默认分析包含 reference ABM；新环境使用 `./install.sh --abm` 安装其依赖，单独 `./install.sh` 只安装报告依赖。Conda 配方现已列入 `pyreadr`、LightGBM 和 sklearn 等默认 ABM 依赖。只补本次 RDS 读取依赖可运行 `/home/huangj/anaconda3/envs/le8/bin/python3 -m pip install 'pyreadr>=0.5,<1'`；使用其他环境时换成实际的 `ABM_PYTHON`。

选择 ABM 训练时，调度器会在原生 C1 扫描前，使用实际 ABM Python、R 路径及后端参数执行预检；读取器和树模型包会实际导入，并执行 CUDA 矩阵前向/反向计算。默认 `--device cuda`，缺包、CUDA 不可用或 GPU 内核失败立即报错，不自动回退 CPU。日志显示实际 Python、PyTorch/CUDA 版本及 GPU 名称。独立预检可运行 `./le8.sh c1_abm --Y cvd_cad --biom prot --preflight`，不会读取整个人群或拟合模型。显式 CPU 测试用 `--abm-args '--device cpu'`；旧 CPU 基线用 `--abm-args '--abm-design selective --device cpu'`。数据预处理、sklearn/LightGBM 基线和统计汇总仍使用 CPU；保存的神经模型权重转到 CPU 以便跨设备加载，运行中的 Transformer 预测保持请求的设备。

本次 `pyreadr` 失败留下的代谢组 C1 已恢复到 `/mnt/d/analysis/le8/cvd_cad/met/c1_correlate`，MWAS、PGS 和汇总 RDS 校验值均未改变，20 张 PNG 已配套工作簿。原诊断目录 `/tmp/le8-run-2k8coq7z` 保留。可继续用原命令 `./le8.sh --Y cvd_cad --biom met`；保持默认 `--replace FALSE`，程序按数值方法及输入签名复用扫描缓存，必要时重建图表。这次修复没有改动 R 分析方法或重新拟合 C1。 新增真实 RDS 读取、缺依赖提前退出及调度顺序回归检查，当前 Python/C5/dispatcher 共 65 项通过。

真实代谢组 ABM 在训练和外部预测的输入阶段统一将负值设为缺失（包括 `Non_HDL_C`），再执行 `log1p` 等变换；不把负值当成零，也不在全体数据上插补。默认输入是 `rap/met.tab.gz`，仅修改 `Rdata/met.rds` 不会影响这一入口。蛋白组 Olink NPX 为有符号的 log₂ 尺度，保留负值，默认 `--transform none`。reference selective、attention 与 TabICL 共用输入规则；质控阈值、插补和标准化参数仍只在训练数据内拟合。`0f/phenotype.R` 的 `clean_biom(nonnegative = TRUE)` 用于 met、bbc，清理发生于缺失率筛选之前，原负值位置在 KNN 插补后仍保留 NA；prot 和影像保持其原有尺度。

原始代谢物映射改为一次构建 DataFrame，避免逐列插入造成的碎片化；RDS 读取将 R `NA_real_` 的特殊 NaN 位模式统一为标准缺失值，保留有限数值（包括用作类别编码的极小值）。输入日志记录变换定义域诊断，selective 日志显示读取、样本质控和每折开始/完成；依赖预检通过并不表示真实数据训练已完成。 公共入口的 ABM 签名使用正式结果根目录作为身份，临时工作区变化不会使同一分析失效；数据、代码或科学参数变化仍会阻止错误复用。

## 2026-10-05 实现审计修复（A01–A08）

本轮按 `LE8_implementation_audit_20261005.md` 修复公共配置、比较定义和证据链。没有删除旧结果，也没有启动真实 UKB 全量重跑。改变这些定义后，应让对应模块检查并更新受影响的缓存；不要把旧结果仅改标签后当成新版结果。

- **A01：共同配置与家庭隔离。** `--group-file` / `LE8_GROUP_FILE` 接收 `eid,group`；`--group-col` / `LE8_GROUP_COLUMN` 指定表型家系列。CSV/TSV 的 ID 与组号按字符串读取，缺失、重复及文件/列冲突报错。ABM reference、attention 和 TabICL 共用配置解析；截止日期为 `--end-date` / `DATE_FOLLOW_END`，诊断日期为 `--Y-date` / `LE8_Y_DATE`。优先级是显式公共 CLI > 环境 > 默认值；覆盖环境时输出提示并记录来源。`--abm-args` 的显式后端参数仍可覆盖公共参数，实际值写入模型 manifest。默认 ABM 协变量含 center，而 R basic 默认不含；需要一致时显式传 `--vars.adj`，实际设计仍以各模型输出为准。
- **共同外层人员。** `--outer-roster` / `LE8_OUTER_ROSTER` 接收 `eid,role`，角色仅为 `training,test`，连接后的分析人员必须完整覆盖，家系不得跨角色。C4/Final 内部把 test 标为 validation；ABM 继续在 training 家庭内部划分 build/tune/calibration。模型保存实际家系来源、覆盖、日期、协变量、角色数、人员/结局/组学 hash；这些信息随现有 workbook/RDS 保存，不新增散落的来源文件。不同模块的排除规则仍可能改变实际测试人员，跨模块配对前需核对最终人员及结局，不能仅凭同 seed 宣称配对。
- **A02：真正的 LE8 replacement。** 基本背景、实测 LE8 和疾病 PRS 分开定义。replacement 从最终风险设计中排除被替代领域的原始测量、积分及已登记的确定性派生项，并在拟合前和实际系数列上断言。默认替代全部八领域；`C4_REPLACE_COMPONENTS=bmi.pts,bp.pts` 可指定部分替代，并记录保留/移除领域。自定义表型派生列通过 `C4_LE8_MEASURE_MAP` 的 `component,variable` 映射补充。原 additive 分支的背景定义保持一致。
- **A03：两种重建评价分开。** `posthoc_panel_reconstruction` 是另外拟合的 panel OLS；`deployed_concept_fidelity` 直接评价进入风险模型的冻结 `cz` 测试预测，输出 R²、RMSE、校准截距/斜率及家庭 bootstrap 区间。R² 使用 Yin 训练均值作参照，两类表均记录目标尺度、队列和缺失处理。新增 fidelity 图及对应工作簿内容，Final 也读取该结果。区间不包含重新选择/拟合模型的不确定性；LP 分量加和检查保留。
- **A04–A05：统一主分析与证据政策。** 正向 MR 主类固定 protein cis / metabolite local；trans/distal 单列，失败/未检验保留，不按跨类别最小 P 选择主结果。C3、共享 helper、C4 候选集及 Final 使用同一政策：`C3_H4`（默认 0.70）、`C3_MR_FDR`（默认 0.05）、完整先验/覆盖以及 `LE8_EVIDENCE_LEVEL=region_or_signal`（默认）或 `signal_only`。政策及 hash 进入结果和缓存；已保存的新格式 C3 结果与当前政策不一致时拒绝混合汇总。逐位点/信号证据保留，区域支持不标为已证明因果。
- **A07：统一聚类。** 观察数据和每次家庭 bootstrap 均调用 `fit_le8_modules()`，复用 discovery/replication profile、方向规则、候选 K、最小模块大小和无可行 K 处理。使用 `cluster::silhouette()`，singleton 为 0。少量特征、单模块、无可行 K 和 B=0 有独立状态；稳定性仍是给定 assay roster 的条件稳定性。

### A06：变异身份与参考规范化

保留 orientation matching key，同时输出 `BUILD,REF,ALT,variant_id,normalization_status,normalization_proof_hash`。不同 build 和同坐标不同 ALT 不混并；冲突重复记录不会按最强 P 保留一条。普通 SNV 仍可按效应方向匹配，但只有验证过 REF/ALT 的条目才标为真实参考规范化身份。

有真实上游规范化证明时，在 `LE8_GWAS_MANIFEST` 对应行指定 `normalization_proof` JSON 路径。证明需包含 `file_sha256,reference_sha256,tool,tool_version,build,reference_checked,left_aligned,multiallelic_split`，后三项为 true；程序验证数据文件 hash、build 和效应等位基因一致性。证明应来自实际的参考校验、左对齐与多等位拆分过程，不能从排序的 EA/NEA 反推。

没有证明时，可在 manifest 提供 `reference_fasta`，或设置 `LE8_REFERENCE_FASTA_37` / `LE8_REFERENCE_FASTA_38`，并提供 FASTA 的 `.fai`。C3 对实际读取的区域运行 `bcftools norm --check-ref e -f ... -m -any`；只处理具有明确且一致 REF/ALT 的标量关联，不猜参考等位基因。无法确认的 INDEL 标为 `indel_normalization_unverified` 并从共定位对齐中排除；不会要求无差别重跑全部 GWAS。规范化后的工具变量使用同一身份匹配。

### A08：C5 公共入口

以下参数由公共入口验证，并且只转发给 C5：`--contrasts`、`--matched-draws`、`--cigma-cells`、`--allow-untested-cigma`、`--no-plots`。C5 使用公共 `--seed`（或 `SEED`），来源路径在写入时统一为正式路径，避免临时工作区改写破坏输出 hash 和缓存复用。`--matched-draws` 为 0 或至少 100；`--cigma-cells` 配合 `--cigma-results`；原生 manifest 与导入结果互斥。`--preflight` 检查显式输入路径与 contrasts 表结构；`--dry-run` 显示实际命令。

```sh
./le8.sh c5_cellulation --Y cvd_cad --biom prot \
  --universe /path/universe.csv --atlas /path/atlas.csv --panels /path/panels.csv \
  --contrasts /path/contrasts.csv --matched-draws 1000 --seed 2026
```

contrasts 表必须有 `model,reference`，对比双方的实际 assay 数及 fold 必须一致，模型名与解析后的 panel 名一致。未提供 contrasts 时明确记录 `not_requested`，不解释成检验阴性。原生 CIGMA 仍需其真实输入；注释完成不代表 CIGMA 已运行。

### 可复现验收

按本次审计“将前次测试套件纳入当前代码版本”的要求，长期回归源码保存在 `validation/`；合成输入、拟合、图片、日志、备份和缓存仍只放 `/tmp`。前次交付中仅适用于安装器/旧宿主打补丁的检查已移除，科学检查改为导入当前整合代码；通过数按本次实际运行计数。

在已经安装 LE8 依赖的环境执行：

```sh
PYTHONDONTWRITEBYTECODE=1 python validation/check_c1.py
PYTHONDONTWRITEBYTECODE=1 python -m pytest -q -p no:cacheprovider validation/test_c5.py validation/test_audit.py
Rscript validation/audit_acceptance.R
```

本轮结果：C1 17 项、Python/C5/dispatcher 58 项、R 18 组验收通过。R 包括实际嵌套 concept 拟合、replacement 风险模型系数、真实 bcftools INDEL 规范化、规范化证明失效、COJO 身份继承、prepared pair 与 `coloc.abf` 数值一致，以及 20 次家庭聚类 bootstrap。

另通过公共 `le8.sh` 完成单结局/蛋白层的隔离 smoke：C1 使用 2,400 人合成数据、家庭文件及自定义截止日期（质控后 2,383 人）；C5 使用显式同预算比较，均完成正常结果打包；实际公共入口与直接调用的配置签名及全部输出 hash 一致。新增 C4 保真度图和 Final/Shiny 索引也完成隔离导出。Python/R 语法与 shell 解析通过。这些是当前代码的合成/入口验证，不是全量 UKB、真实 GPU、TabICL 权重或原生 CIGMA 的验证。

## C1 selective attention（2026-10-05）

按 `C1_Codex_Implementation_20261005` 接入独立的 `--abm-design selective_attention`。默认仍是原 S6 `selective`；原 `attention` 和 TabICL 仍是独立实验。新模式只接受 `--abm-backend reference`，结果放在 `c1_correlate/abm_selective_attention/`，不覆盖 `abm_reference/`。逐项实现、测试与限制见 [C1.md](C1.md)。

新模式保留连续分子的 bulk、上下 tail、missing 通道及二元特征；线性/树和神经模型使用同一原始 assay 与临床信息。复用 RowEncoder/ContextStack，实际训练 Q/K 检索注意力。全流程按家庭分为 build、tune_model、tune_gate、calibration_fit、calibration_audit、test；训练权重来自 recipient 家庭完全不参与的 pilot/gate。释放 gate 比较最终已校准候选与实际 fallback，独立审计证据不足时输出 fallback，不强制释放固定人数。

显式选择新模式的命令如下。家庭文件、固定 outer roster 需替换为实际核实的路径；本轮没有启动这条 UKB 正式命令。没有家庭数据时可不传 `--group-file`，但此时只按个人隔离，不能声称控制了亲缘。

```sh
./le8.sh c1_abm --Y cvd_cad --biom prot \
  --analysis-root /mnt/d/analysis/le8_s7 \
  --group-file /actual/path/families.csv \
  --outer-roster /actual/path/outer.csv \
  --abm-backend reference --seed 2026 \
  --abm-args '--abm-design selective_attention --device cuda --cores 8 --selective-coverage 0.60 --s7-primary attention_retrieval_weighted --bootstrap 500'
```

家庭文件为字符串 `eid,group`；outer 文件为 `eid,role`，角色仅 `training,test`。默认训练家庭命名空间为 `ukb`，其他队列应显式使用 `--group-namespace`。外部投射必须提供当前外部家庭文件或列及核实过的命名空间；不会继承训练时的临时文件，也不会要求外部 CSV 自带 `.le8_family`。开发集所有角色的 ID 均禁止进入外部验证，同命名空间的开发家庭也禁止；改变命名空间不能绕过 ID 重叠检查。

```sh
./le8.sh c1_abm --Y cvd_cad --biom prot \
  --analysis-root /mnt/d/analysis/le8_s7 --abm-backend reference \
  --abm-args 'project --run-dir /mnt/d/analysis/le8_s7/cvd_cad/prot/c1_correlate/abm_selective_attention --phe-file /actual/path/external_baseline.csv --omics-file /actual/path/external_prot.csv --projection-group-file /actual/path/external_families.csv --projection-group-namespace ukb --output /tmp/external_prediction.csv --device cpu --cores 4'
```

这里的 baseline 输入仅需 ID 和训练时所用的基线临床/技术字段，不需结局。正式私有预测可将 `--output` 指向分析目录中的 CSV 交换路径，由公共事务封入同名 RDS。已有旧 bundle 如果没有完整 development roster，外部验证会明确拒绝；需要重拟合才能补齐，不能把仅 context 名单当作完整名单。

新模式的主指标为固定 horizon 的 IPCW 风险指标与真正的 Uno C；家庭 bootstrap 对各模型使用相同抽样。年龄三分位/性别的条件删失敏感性仅在 build 拟合，作为另表诊断，不替换主策略；该敏感性不声称估计条件删失 Uno C。Final/Shiny 分开展示新旧模式、实际架构、候选/fallback、研究/释放覆盖和审计状态。

依赖沿用 `environment.yml` / `requirements.txt` 的 PyTorch、scikit-learn、SciPy、pandas、joblib、threadpoolctl；默认树为 LightGBM，CPU smoke 显式选 `--tree hist`。RDS 仍需 pyreadr；打包需现有 R 环境。新模式不下载 TabICL 权重。以下测试均使用合成数据，输入与结果留在 `/tmp`；结构模拟不要求神经模型胜出。

```sh
PYTHONDONTWRITEBYTECODE=1 python -m pytest -q -p no:cacheprovider validation/test_c1_s7_kernels.py validation/test_c1_s7_production.py
PYTHONDONTWRITEBYTECODE=1 python -m pytest -q -p no:cacheprovider validation/test_c1_s7_simulations.py
PYTHONDONTWRITEBYTECODE=1 python -m pytest -q -p no:cacheprovider validation/test_c1_s7_devices.py
PYTHONDONTWRITEBYTECODE=1 python -m pytest -q -p no:cacheprovider validation/test_c1_s7_public.py
```

## 分开分析，最后汇总

可以在不同终端运行以下命令，使用同一个默认结果目录：

```sh
./le8.sh --Y cvd_cad --biom prot
./le8.sh --Y cvd_cad --biom met
./le8.sh --Y ra --biom prot
```

全部分析完成后，统一汇总多个疾病和数据层：

```sh
./le8.sh final,shiny --Y cvd_cad,ra --biom prot,met
```

默认流程为 `c1_correlate,c1_abm,c2_cause,c3_coloc,c4_connect,c4_panel_validation,c5_cellulation`。`final,shiny` 只汇总已有结果，不自动训练 ABM；显式模块列表只运行所选模块及其必要前置分析。需要跳过默认 ABM 时使用 `--skip-abm`，需要在其他模块组合中加入 ABM 时仍可使用 `--run-abm`。

C1 PGS 按特征原子保存检查点，终端显示已完成数、复用数和估计剩余时间。相同输入、模型和参数下重新执行原命令即可续跑；输入或方法变化不会混用检查点，`--replace TRUE` 会重新计算。每个特征仍计算完整的六类模型，完成整个特征家族后统一计算 FDR。`--cores N` 现在也控制 PGS 特征并行；默认仍为 1，以控制完整 PGS 矩阵与模型数据的内存开销，多终端并发时需合计各任务的内存。已完成的 PWAS/MWAS 扫描同样保留在 `/tmp/le8-cache/` 供重启复用，检查原始输入、协变量、数值方法及相关设置；修改调度或报告代码不会单独触发这些扫描重算。清理 `/tmp` 或系统重启后，未发布的检查点可能丢失。

## 结果文件

- 每张 PNG 配同名 XLSX，只放对应图的分析结果；同类 panel 用 panel 列区分，尽量合并为少量 worksheet。例如 `c3.Fig2.regional_top_loci.png` 对应 `c3.Fig2.regional_top_loci.xlsx`。重要的未绘图结果按主题另外保存。
- 未绘图结果按分析内容组织，不按表数量分卷：C1 使用 `c1.association.xlsx`，PGS 配对比较、分解、时间分析分别为 `c1.pgs.comparison.xlsx`、`c1.pgs.decomposition.xlsx`、`c1.pgs.temporal.xlsx`；C4 验证按表现、面板、模型和遗传证据区分。Final 和 Shiny 同样按主题组织，不再生成 `.2.xlsx`、`.3.xlsx`。
- 工作簿内部同时保存原始汇总导出及校验信息，供重绘、缓存复用和 Shiny 读取；不再另存通用汇总 RDS。完整区域／SNP 结果留在可复用分析 RDS，图表只导出实际展示的范围，不重复嵌入旧工作簿。程序可以在 `/tmp` 完整恢复这些导出。
- 工作簿是可复用结果，需保留。Excel 等编辑器重新保存时可能移除内部来源记录；需要编辑展示版时，另存到其他目录。
- 带个体 ID 的表按具体数据名称保存为 `.rds`，例如 `test_individuals.rds`、`individual_explanations.rds`、`c4.focus.roles.rds`，不进入汇总工作簿或 Shiny。嵌套的个体解释和 attention 数组也保存在 RDS 内，不另留 JSONL、NPZ 文件。
- 不单独输出分析参数、代码版本或输入路径清单，不生成 `analysis_options`、`revision_manifest`、`analysis_manifest` 文件。缓存复用所需的最少一致性信息保存在相应分析结果对象内部。
- 可复用的 R 模型和分阶段结果保留为 RDS。名称简短且固定，例如 `c2.dandelion.rds`、`c2.instruments.rds`、`c2.mr.rds`、`c2.reverse_mr.rds`；版本信息放在对象内部。
- `le8.sh` 在 `/tmp/le8-run-*` 中展开交换表、执行分析和重绘。成功后发布结果工作簿、PNG 和个体/模型 RDS；运行失败时保留 `/tmp` 中的诊断现场。

## 目录与阅读顺序

| 目录 | 内容 |
|---|---|
| `<疾病>/<prot或met>/c1_correlate/` | C1 工作簿、PNG 和可复用分析对象；`abm_reference/`、`abm_tabicl/` 各保存一个 ABM 方法的结果 |
| `<疾病>/<数据层>/c2_cause/` | 按图配对的 XLSX、PNG 和命名明确的 RDS |
| `<疾病>/<数据层>/c3_coloc/` | 共定位汇总、PNG 和可复用结果；GPU 格式转换和逐区域计算缓存在 `/tmp` |
| `<疾病>/<数据层>/c4_connect/` | LE8 连接、代理、交互、非线性及固定预算验证结果；每张 PNG 配同名 XLSX，验证拟合保存在 `c4.validation.rds`；Yin、YinYang 区分代理发现人群 |
| `<疾病>/<数据层>/c5_cellulation/` | 细胞注释汇总和 PNG；自动构建的输入表在 `/tmp` |
| `final/` | Fig1–8、补充图及各自的 XLSX、`report.html`、`index.html` 和图注；研究问题综合结果也在本层，不另建 `overview/`、`tables/`、`panels/` 或结果 README |
| `shiny/` | 按主题拆分的索引工作簿，含查看器索引和综合视图；表格详情、图像仍通过相对路径引用各模块结果 |

先看 `final/report.html` 或 Shiny 总览，再按疾病、数据层和模块查阅结果。工作簿使用简短的结果页名，不保留空表、重复表和无用的管理页。已有 PNG 时不生成同名 PDF；中间 panel 图、PGS 扫描缓存、MR-link-2 工作文件、GPU 转换文件和日志均留在 `/tmp`。ABM 的输入副本、逐折神经网络检查点也放在 `/tmp/le8-cache/`；正式目录保留自包含的最终 `model_bundle.joblib` 和必要的个体 RDS。ABM 的冻结状态、校准信息等复用所需信息收进方法工作簿内部，不再散落为 JSON 文件。不因整理重新拟合模型。

## 分享 Shiny

只发送 `shiny/` 和 `final/` 不足以打开全部详情与图片。生成完整的汇总查看包：

```sh
./le8.sh share --share-out /mnt/d/analysis/le8-share.zip
```

发送这个 ZIP 即可。包内只有两份查看代码（`shiny/app.R`、`f/0.common.R`）、说明和 `results/` 中所需的汇总工作簿与 PNG；不包含个体 RDS、训练权重或原始 UKB 数据。创建时会核对来源校验值、检查个体 ID 列，并在独立目录验证 Shiny 读取。

接收方安装 R，以及 `shiny`、`DT`、`data.table`、`ggplot2`、`digest`、`openxlsx`、`jsonlite` 包。解压后在包的根目录运行：

```sh
Rscript shiny/app.R
```

浏览器打开 `http://127.0.0.1:3839`。查看包不需要 Python、WSL 或原作者的绝对路径；请保留包内相对目录结构。Windows 的 RStudio 也可切换工作目录到包内 `shiny/`，然后运行 `source('app.R')`。

自己的原结果仍可这样打开：

```sh
./le8.sh shiny --no-reindex --analysis-root /mnt/d/analysis/le8
```

运行 `./le8.sh` 默认只执行 C1–C5（包含 reference ABM 和 C4 面板验证），完成后退出。汇总与查看单独使用 `./le8.sh final,shiny`；只生成汇总而不启动服务，可加 `--prepare-only`。临时验证输出、备份和 Python 缓存均放在 `/tmp`；本次审计要求纳入版本的回归源码见 `validation/`。

当前测试范围与实测通过数见上方“可复现验收”。SuSiE 的真实分析仍需提供与各 GWAS 的 `ancestry` 一致的 LD 元数据；未知效应尺度、缺失强信号或未完成的分析不会升级为完整 MR–coloc 支持。Final 默认展示 10-assay 预算，若未运行该预算则展示最近的已配置预算，选择不依据预测结果。

可选扩展的边界：本版未启用组织 eQTL/sQTL 扩展或 SuSiE 信号 LBF 的 GPU 批处理；GPU 主线比较的是相同区域边际 BF。风险尺度因果中介和年龄断点搜索也未启用，现有输出分别为描述性路径乘积和基线年龄平滑曲线。原生 CIGMA 仍需匹配的单细胞表达和基因型输入。
