# LE8 5C 多组学分析流水线

LE8 是一个在 WSL/Linux 环境中运行的 Bash、R 和 Python 流水线，用于对 UK Biobank 蛋白组（`prot`）和代谢组（`met`）数据执行 C1–C5 主分析，以及 S1–S2 补充分析。统一入口为 `le8.sh`。

## 快速开始

代码位于 Windows 的 `D:\scripts\le8`；在 WSL 中对应 `/mnt/d/scripts/le8`。

```bash
cd /mnt/d/scripts/le8

# 检查脚本、数据路径和软件依赖
./le8.sh -Y cvd_cad,masld --biom prot,met --preflight

# 只打印执行计划，不启动分析
./le8.sh -Y cvd_cad,masld --biom prot,met --dry-run

# 正式运行
./le8.sh -Y cvd_cad,masld --biom prot,met
```

默认会复用有效缓存。只有需要忽略模块缓存并重新计算时才加入 `--replace TRUE`。

## 分析模块

模块按下列顺序执行：

| 模块 | 功能 |
|---|---|
| `c1_correlate` | 观察性 ProtWAS/MWAS、诊断前时间模式及协变量敏感性分析 |
| `c2_cause` | cis/local 与 trans/distal MR、反向 MR、MR-link-2、DANDELION 和个体遗传分解 |
| `c3_coloc` | 位点共定位、精细定位及可信集诊断 |
| `c4_connect` | LE8 监督的代理特征发现、聚类、网络和中介分析 |
| `c5_consolidate` | 多来源证据整合、预测面板和可选的嵌套交叉验证 |
| `s1_interact` | 性别分层及 LE8 两两交互 |
| `s2_nonlin` | 非线性关联和交叉验证惩罚模型 |

依赖会自动补齐：运行 C2、C3、C4 或 S1 时会加入 C1；运行 C5 时会加入 C1–C4。

## 安装环境

建议使用仓库内的 Conda 环境文件：

```bash
cd /mnt/d/scripts/le8

# 首次创建
conda env create -f environment.yml

# 已有 le8 环境时更新
conda env update -n le8 -f environment.yml --prune
```

`environment.yml` 包含主要运行依赖，包括 R、Python、PLINK/PLINK2、MR-link-2 所需包和 GPU-coloc。`le8.sh` 会在非交互式 WSL 会话中尝试定位用户目录下的 Conda，并通过 `conda run -n le8` 调用相应程序。

## 默认数据路径

`le8.sh` 当前使用以下默认路径，均可用命令行参数覆盖：

| 内容 | 默认路径 |
|---|---|
| 分析输出根目录 | `/mnt/d/analysis/le8` |
| 结局 GWAS 项目根目录 | `/mnt/d/data.BIG/gwas/main` |
| pQTL 项目根目录 | `/mnt/d/data.BIG/gwas/prot` |
| mQTL 项目根目录 | `/mnt/d/data.BIG/gwas/met` |
| 1000 Genomes 参考目录 | `/mnt/d/data.BIG/refGen/1kg` |
| UKB phenotype 目录 | `/mnt/d/data/ukb/phe` |
| 共享 phenotype Shell 库 | `/mnt/d/scripts/0f/0phe.f.sh` |
| 匹配过程临时目录 | `/mnt/d/tmp` |
| 蛋白遗传力文件 | `/mnt/d/files/ppp.h2.csv` |
| 蛋白遗传力来源说明 | `/mnt/d/files/ppp.h2.txt` |

UKB phenotype 目录通常包含 `Rdata/all.rds`、`Rdata/prot.rds`、`Rdata/met.rds`、`Rdata/prot.pgs.rds` 和 `Rdata/met.pgs.rds`。蛋白注释 BED 会优先从 pQTL 项目目录自动定位 `ppp_3k.b38.bed`，并可回退到 `/mnt/d/files/ppp_3k.38.bed`。

### GWAS 目录约定

按 trait 组织的基本结构为：

```text
<project>/common/<trait>/gwas/<trait>.gz
<project>/common/<trait>/gwas/<trait>.cis.gz
<project>/common/<trait>/gwas/<trait>.jma.cojo
<project>/common/<trait>/gwas/<trait>.ldr.cojo
<main-project>/common/<trait>/magma/<trait>.genes.out
```

`--grch auto` 会按 trait 自动识别 GRCh37/GRCh38。单 trait 诊断运行也可以用 `--outcome-gwas FILE` 直接指定结局 GWAS。

## 已内置的本地配置

以下配置已经写入 `le8.sh`，当前运行无需额外 `export`：

```bash
C1_LE4_COVARS="diet.pts,pa.pts,smoke.pts,sleep.pts"
C1_FULL_LE8_SENSITIVITY="TRUE"
C1_TREATMENT_VARS="drug.lipid"

C2_PROT_HERITABILITY_FILE="/mnt/d/files/ppp.h2.csv"
C2_MET_HERITABILITY_FILE=""
C2_HERITABILITY_FILE=""
```

`drug.lipid` 是由 UKB 6153 和 6177 字段各次访视中的降脂药类别生成的分析用二元变量。C1 会把它作为治疗协变量，并把基础协变量加 LE4 作为主调整集；完整 LE8 调整保留为敏感性输出。

蛋白和代谢物遗传力输入分别配置，避免把蛋白遗传力误用于代谢物。遗传力文件支持 CSV 或 RDS，至少应包含：

| 字段 | 必需 | 含义 |
|---|---:|---|
| `feature` | 是 | 与组学特征匹配的标识 |
| `snp_h2` | 是 | SNP 遗传力估计 |
| `snp_h2_se` | 否 | SNP 遗传力标准误 |

这些值仍可在启动命令前通过同名环境变量覆盖。

## 命令行

查看程序内置、最权威的参数说明：

```bash
./le8.sh --help
```

常用参数：

| 参数 | 说明 |
|---|---|
| `-Y, --trait TRAITS` | 一个或多个 trait，多个值用逗号分隔 |
| `-b, --biom LAYERS` | `prot`、`met` 或 `prot,met` |
| `-s, --steps JOBS` | 指定模块，多个值用逗号分隔 |
| `--from JOB` / `--to JOB` | 按固定模块顺序选择起止范围 |
| `--replace TRUE\|FALSE` | 是否忽略模块缓存；默认 `FALSE` |
| `--preflight` | 检查环境、路径和必需文件后退出 |
| `--dry-run` | 打印执行计划后退出 |
| `--analysis-root DIR` | 覆盖输出根目录 |
| `--cores N` | R worker 数；默认 4 |
| `--memory-limit-gb N` | 整个 LE8 进程树的 RAM 硬上限；默认 32 GiB，设为 `0` 可禁用 |
| `--memory-swap-gb N` | cgroup 可额外使用的 swap；默认 4 GiB |
| `--seed N` | 随机种子；默认 2026 |
| `--grch auto\|37\|38` | 基因组版本；默认自动识别 |

数据路径参数：

- `--gwas-dir DIR`、`--pqtl-dir DIR`、`--mqtl-dir DIR`
- `--prot-bed FILE`、`--refgen-root DIR`、`--ukb-phe DIR`
- `--outcome-gwas FILE`、`--shared-shell FILE`、`--r-bin FILE`

可选计算参数：

- `--run-mrlink2 Top|All|None`：MR-link-2，默认 `Top`。
- `--run-dandelion Top|All|None`：蛋白 DANDELION，默认 `Top`。
- `--dandelion-gene-p FILE`：基因层面的结局 P 值。
- `--dandelion-allow-magma TRUE|FALSE`：未指定基因 P 值时是否允许使用 MAGMA 输入。
- `--dandelion-gene-annotation FILE`、`--dandelion-snp-file FILE`、`--dandelion-snp-gene-map FILE`、`--dandelion-fdr FLOAT`。
- `--run-gpu-coloc TRUE|FALSE`：是否运行 GPU-coloc，默认 `TRUE`。
- `--nested-cv TRUE|FALSE`：是否运行 C5 嵌套交叉验证，默认 `FALSE`。
- `--inc_prot FILE`、`--inc_met FILE`：限制 C5 使用的候选特征。
- `--direction-anchors CSV`：覆盖 C1/C2/C5 的方向锚点。
- `--c4-module-boot N`、`--c4-module-stability X`：C4 模块稳定性参数。
- `--c5-topn-max N`、`--c5-topn-boot N`、`--c5-lead-min-cases N`、`--c5-mechanism-n N`：C5 面板参数。

MR-link-2 参考数据可通过 `--mrlink2-ref-bed`、`--mrlink2-ref-pfile-dir`、`--mrlink2-ref-pop`、`--mrlink2-ref-id-dir` 和 `--mrlink2-ref-samples` 配置。

## 常用运行方式

仅运行 C1 和 C2：

```bash
./le8.sh -Y cvd_cad --biom prot --steps c1_correlate,c2_cause
```

从 C4 运行到最后：

```bash
./le8.sh -Y cvd_cad --biom met --from c4_connect
```

运行多个 trait 和两层组学：

```bash
./le8.sh -Y cvd_cad,masld --biom prot,met
```

临时关闭可选的 MR-link-2、DANDELION 和 GPU-coloc：

```bash
./le8.sh -Y cvd_cad --biom prot \
  --run-mrlink2 None \
  --run-dandelion None \
  --run-gpu-coloc FALSE
```

强制重算所选模块：

```bash
./le8.sh -Y cvd_cad --biom prot,met --replace TRUE
```

## 缓存与续跑

默认 `--replace FALSE`。每个模块会检查已有产物和缓存签名，并尽量从有效结果继续运行。C1/C2 的最终缓存包含代码和关键配置签名；部分耗时阶段还使用可独立复用的阶段缓存。

已有的 `/mnt/d/analysis/le8` 可以保留用于续跑。修改输入、关键配置或代码后，如果签名不再匹配，相关结果会重新生成；需要完全忽略模块缓存时使用 `--replace TRUE`。

## 输出

默认输出根目录：

```text
/mnt/d/analysis/le8/
└── <trait>/
    ├── prot/
    │   └── <module outputs>
    ├── met/
    │   └── <module outputs>
    └── logs/
```

各模块按适用情况写出 CSV、RDS、XLSX 和图形文件，并生成输入状态、匹配率、方法范围、缓存状态或恢复任务等审计文件。运行日志按 trait 和时间记录在 `logs` 下。

## 资源控制

SNP/基因组版本匹配默认采用磁盘外排，减少大文件在内存中的复制：

| 参数或变量 | 默认值 |
|---|---:|
| `--match-memory-mb` / `LE8_MATCH_MEMORY_MB` | 4096 MB |
| `--match-sort-memory-mb` / `LE8_MATCH_SORT_MEMORY_MB` | 512 MB |
| `--match-tmp-dir` / `LE8_MATCH_TMP_DIR` | `/mnt/d/tmp` |
| `--cores` / `N_CORES` | 4 |
| `--memory-limit-gb` / `LE8_MEMORY_LIMIT_GB` | 32 GiB（整个任务进程树） |
| `--memory-swap-gb` / `LE8_MEMORY_SWAP_GB` | 4 GiB |

包装器默认通过 WSL 的 cgroup v2 为完整任务树设置 `32 GiB RAM + 4 GiB swap` 的硬上限，并把常见 BLAS/OpenMP 线程数设为 1。达到上限时，systemd 只终止 LE8 任务组，不会把调用它的 terminal 放进同一个受限组。没有可用的 user systemd 时，脚本会明确警告并退回到单进程 `ulimit`。运行前应确保临时目录和输出目录有足够磁盘空间。

例如，把总 RAM 上限改为 24 GiB、禁止使用 swap：

```bash
./le8.sh --memory-limit-gb 24 --memory-swap-gb 0 [其他参数]
```

## 故障排查

1. 先运行与正式命令相同参数的 `--preflight`。
2. 再运行 `--dry-run`，确认 trait、组学层、模块顺序和解析后的路径。
3. 数据未找到时，用相应的路径参数覆盖默认值。
4. 依赖缺失时，用 `environment.yml` 更新 `le8` Conda 环境。
5. 可选工具暂不可用时，可把 MR-link-2、DANDELION 或 GPU-coloc 设为 `None`/`FALSE`。
6. 大型匹配任务失败时，检查 `/mnt/d/tmp` 空间，并按需要调整匹配内存或临时目录。
7. 从旧输出续跑异常时，查看对应 trait 的日志和缓存审计；确认需要全量重算后再使用 `--replace TRUE`。

## 项目结构

```text
le8/
├── le8.sh                 # 统一命令行入口
├── environment.yml        # Conda 环境
├── README.md              # 本说明
└── f/
    ├── c1_correlate.R
    ├── c2_cause.R
    ├── c2_genetic_decompose.R
    ├── c2_mr_link2.py
    ├── c3_coloc.R
    ├── c3_coloc_GPU.py
    ├── c4_connect.R
    ├── c5_consolidate.R
    ├── s1_interact.R
    ├── s2_nonlin.R
    └── comm.f.R
```
