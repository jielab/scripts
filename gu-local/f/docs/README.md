# GU FINAL 2026：从新代码、新结果开始运行

版本日期：2026-10-31。入口为 `gu.sh`、`arg.sh` 和 `install.sh`；HPC 包另有 LSF 调度入口 `gu-hpc.sh`。所有科学分析从本版新建的结果目录开始，结果数据库版本为 **2026**。

IBDmix 只有一个模块名 **`ibdmix`**：`--GPU TRUE` 调用 CUDA 分段实现，`--GPU FALSE` 或省略此参数使用原生 CPU 实现。两者共用 `ibdmix/` 结果布局。默认参考在调用、检查、命令生成和汇总中均为：

```text
Altai,Chagyr,Denisova,Denisova25,Vindija
```

汇总只显示实际完成的参考；没有运行的参考不会自动插入。`Altai.2013`、`Denisova.2013` 只接受显式 `--archaic-gen` 选择，必须提供对应版本的真实 VCF 与 mask，不把当前参考文件重命名使用。

## 1. 部署目录与配置

ZIP 包顶层包含 `gu/`、`0f/`、`grid/f/`、`RELEASE.json` 和 `SHA256SUMS`。三个代码目录必须保持同级；共享辅助文件已随包提供。

| 内容 | 本机 | HPC |
| --- | --- | --- |
| 代码根目录 | `/mnt/d/scripts/` | `/work/sph-huangj/scripts.HPC/` |
| GU 入口目录 | `/mnt/d/scripts/gu/` | `/work/sph-huangj/scripts.HPC/` |
| 数据盘根目录 | `/mnt/f/` | `/data/sph-huangj/` |
| 基因组参考根目录 | `/mnt/f/gen/` | `/data/sph-huangj/gen/` |
| 结果默认根目录 | `/mnt/f/analysis/gu/` | `/data/sph-huangj/analysis/gu/` |
| 软件安装目录，可覆盖 | `/mnt/d/software/gu/` | `/work/sph-huangj/software/gu/` |

将 ZIP 内的三个代码目录放到对应代码根目录下。**不要形成 `scripts/gu-FINAL-local-2026/gu/` 这样的额外层级作为预定安装位置。** 解压后可在包顶层先运行 `sha256sum -c SHA256SUMS` 核验文件，再复制代码。

每个包的 `gu/f/runtime.env` 已设置对应路径，公开脚本自动加载，无需手工替换源码中的盘符。显式导出的环境变量优先；要使用另一份配置，可设置 `GU_RUNTIME_ENV=/absolute/path/runtime.env`。交互式 shell 中若要使用下面的 `$GU_REF_ROOT` 等变量，先加载一次配置：

```bash
cd /mnt/d/scripts/gu
source ./f/runtime.env
./gu.sh --help
```

本机 conda 默认 `/mnt/d/software/anaconda3`；HPC 默认 `/share/apps/anaconda3/2020.7`。这两个软件路径可按实际安装通过 `GU_CONDA_BASE` 覆盖。激活 conda 时也可从 `CONDA_EXE` 识别本机安装。额外软件设置示例在 `local.runtime.env.example`，HPC 详见 `HPC.README.md` 与 `hpc.runtime.env.example`。

F: 盘镜像只对应 `/data/sph-huangj/` 下相同相对路径。以下原本位于 D: 的数据需要额外复制或另行指定：

- `CHECK_GRCH_SNP_LIST`：本机 `/mnt/d/data/ukb/common/snp.lst`；HPC `/data/sph-huangj/ukb/common/snp.lst`。
- `GU_COJO_CHAIN`：本机 `/mnt/d/files/liftOver/hg38ToHg19.over.chain.gz`；HPC `/data/sph-huangj/files/liftOver/hg38ToHg19.over.chain.gz`。仅需要坐标转换的 PhyML 输入使用。

ZIP 包提供代码，不含基因型、古基因组、FASTA、AXT、mask、遗传图谱、SNP 校验表和平台二进制。

## 2. 安装

```bash
# CPU 环境及固定版本的上游程序。
bash install.sh --GPU FALSE --dir-software "$GU_SOFT"
bash install.sh --check
bash install.sh --check-references

# 本机有 NVIDIA GPU、可用驱动与 nvcc 时，额外编译 IBDmix CUDA。
bash install.sh --GPU TRUE --dir-software "$GU_SOFT"
```

安装入口固定 IBDmix `6a9d4bc67dcaf602e65cac0fa6191e1e7a1291ef` 和 TRACE `70cab6c6687900826eb0b2ded69c39554b1b4146`。环境求解后保存 `gu-conda-explicit.txt` 与 `gu-pip-freeze.txt`。`--check-references` 检查五个默认参考；按子集运行时，使用对应 `gu.sh ibdmix check` 请求检查实际输入。

主要资源路径为 `GU_REF_ROOT`、`GU_ARCHAIC37_ROOT`、`IBDMIX_MASK_ROOT`、`IBDMIX_REFERENCE_FASTA`、`GU_ANNOT_ROOT`、`IBDMIX_AXT_DIR`。后两者默认分别为数据盘的 `annot/` 与 `annot/axt/37/`，用于 CpG 排除区域计算。输入 VCF/BCF 的 `.tbi/.csi` 及 FASTA 索引应预先准备好。

ARG 的 Threads 使用独立 Python ≥3.12 环境；主 GU 环境为 Python 3.11。需要时安装：

```bash
export GU_CONDA_EXE="$GU_CONDA_BASE/bin/conda"
"$GU_CONDA_EXE" create -n gu-arg-bootstrap python=3.12 -y
export GU_THREADS_PYTHON="$GU_CONDA_BASE/envs/gu-arg-bootstrap/bin/python"
export GU_THREADS_ENV="$GU_SOFT/venvs/gu-threads"
export GU_SINGER_ENV="$GU_SOFT/venvs/gu-singer"
bash install.sh --arg
```

`--arg` 安装 Threads 与区域 SINGER 的独立环境；`tsinfer + tsdate` 在主 GU 环境中。Needle 使用独立安装、遗传图谱和其 ASMC 后端，需要自行提供真实路径。各后端的科学模型和适用样本规模不同。

## 3. 输入与公开参数

| 参数或输入 | 约定 |
| --- | --- |
| `--dir-target PREFIX` | 染色体文件前缀，例如 `/mnt/f/gen/1kg/37/pfile/chr` 对应 `chr1.pgen/.pvar/.psam` 等；也支持带索引的多样本 VCF/BCF |
| `--target NAME` | 与目标前缀一起指定数据集名称，例如 `1kg`、`ukb` |
| `--sample-panel FILE` | 样本元数据包含 `sample`、`pop`、`super_pop`、`sex`；男性为 1；ID 应唯一且与基因型一致 |
| `--archaic-path DIR` | 古基因组 VCF 根目录；每个参考有自己的子目录，每条染色体必须唯一匹配 |
| `--archaic-gen LIST` | IBDmix 的参考选择；默认五个当前版本；逗号分隔 |
| `--dir-out DIR --postfix ID` | 结果根为 `DIR` 直接拼接 `ID`；如 `ukb` + `001` → `ukb001` |
| `--chr 1-22,X` | 全染色体扫描；支持逗号、空格、区间；X 也接受 23 |
| `--loci FILE` | 受限区域或 GWAS lead 输入；普通 BED 为 0-based、half-open |
| `--keep-psam FILE,SIZE,INDEX` | 先按 PSAM 数据行切片，批次从 1 开始，再应用 `--keep-males` |
| `--jobs N --memory-cap 32G` | 单次 GU 请求的工作进程和内存预算 |
| `--dir-arg DIR` | TRACE 所需树、样本映射与 `ARG_TRACE_BUILD.tsv` 的直接所在目录 |

VCF POS 是 1-based；统一片段表为 BED 风格区间。原生调用文件保留供复核。GRCh37 是内置 IBDmix 排除 mask 的版本；其他 build 需要 `IBDMIX_PROFILE=custom` 与明确的 `IBDMIX_MASK_DIR`。默认 `IBDMIX_PROFILE=multi_reference`。

IBDmix 从本次选定现代人群估计频率；多个独立单人 VCF 的分别运行不等于一份人群 VCF 的联合扫描。TRACE 使用现代人群 ARG，不依赖古基因组参考，但需要符合其约定的样本身份、祖先状态、突变与时间信息。

## 4. IBDmix：1KG 与 UKB

### 1KG 全基因组 CPU

```bash
./gu.sh ibdmix --GPU FALSE \
  --grch 37 --target 1kg --chr 1-22,X \
  --dir-target "$GU_REF_ROOT/1kg/37/pfile/chr" \
  --sample-panel "$GU_REF_ROOT/1kg/38/samples.txt" \
  --archaic-path "$GU_REF_ROOT/archaic/37/vcf" \
  --dir-out "$GU_ANALYSIS_BASE/1kg2026" --jobs 4
```

这里 `38/samples.txt` 提供的是样本元数据，不改变 37 版基因型坐标；若样本表存放在其他位置，显式指定即可。不写 `--archaic-gen` 时就是五参考全部运行。

### 本机 UKB GPU，每批 10,000 行 PSAM

```bash
cd /mnt/d/scripts/gu
./gu.sh ibdmix --GPU TRUE \
  --grch 37 --target ukb \
  --dir-target /mnt/f/gen/ukb/37/imp/chr \
  --dir-out /mnt/f/analysis/gu/ukb --postfix 001 \
  --keep-psam /mnt/f/gen/ukb/37/imp/chr1.psam,10000,1 \
  --keep-males
```

结果根为 `/mnt/f/analysis/gu/ukb001`。下一批同时改 `--postfix 002` 与 PSAM 参数的末尾 `,2`。先取 10,000 行再筛男性，因此每批实际男性数通常小于 10,000。批次标签与实际参与分析的 roster 写入运行记录。

每批都是独立的频率估计人群，拼接几十批不保证与全 UKB 一次共同频率模型相同。比较 CPU/GPU 必须使用相同样本、参考、mask、阈值与输入版本，并写到不同结果根。CUDA 路径仍在 CPU 准备基因型和计算频率/LOD；设备上执行分段调用，不承诺整条链路全部在 GPU。

原生阈值默认 LOD ≥4、最短片段 50 kb。UKB 没有 AFR5 对照时关闭 African Denisovan subtraction，并将其状态写入 provenance；不要将这一设置全局套用于包含 AFR 对照的 1KG 请求。单选参考时，相关后处理依赖由预检查明确验证。

可用 `--archaic-gen Altai` 只分析一个当前参考。需要额外参考时显式传入完整集合；不会自动扩充参考。可选 DAF 排除 BED 通过 `GU_IBDMIX_DAF_MASK` 指定，必须配有精确参考、build、坐标与来源信息的 `.json` sidecar。缺失时明确记录未应用；不从现有统计结果反推一个 mask。示例 sidecar：

```json
{"genome_build":"GRCh37","reference":"Altai","percentile":99.9,"coordinate_system":"0-based-half-open","source":"实际文件来源与对应参考版本"}
```

### 命令计划和资源

在分析命令末尾加 `--run-cmd FALSE`，会检查请求并保存 `.cmd`，不执行推断。`--auto-final` 默认 FALSE；所有工作完成后显式运行 `final`。

本机 CPU 默认内存预算 32G，CUDA 默认 16G。`--resource-guard auto` 在实际 LSF 作业内使用 scheduler，本机优先使用 user systemd。CPU 的 `ulimit` 回退只限制单进程虚拟内存，不是整个进程树 RSS 上限。CUDA 与 tsinfer 的大虚拟地址空间需要 systemd 或真实调度器资源控制。`--jobs` 控制准备并发，GPU caller 在同一设备串行。

## 5. chrX 与汇总指标

男性 non-PAR X 使用显式的倍性检查、两条输入流一致的非 PAR 筛选和分块调用。优先提供 `chrX.male` 文件与实际男性元数据。新 X 单元目录含 `.nonpar-v2.mac1`，表示当前 non-PAR 处理及 MAC 设置。PAR locus 与整条 X 分析是不同请求。

当前 IBDmix 男性 X 将单倍型 GT 变成伪二倍体输入，仍用原始二倍体似然。一个男性替代等位基因贡献两个伪拷贝，因此 MAC 不等于男性携带人数；X 与常染色体的检出率仍需真实数据与模型校准。代码不为了让比例接近而调整结果。

令每人的片段坐标并集长度为 `U`、实际纳入汇总的染色体物理跨度为 `L`：

| 输出 | 定义 |
| --- | --- |
| `mean_mb` | `mean(U) / 1e6`，单位 Mb/人 |
| `coordinate_coverage_pct` | `100 × U/L`，X 与常染色体同为物理坐标覆盖率 |
| `coverage_pct` | 常染色体 `100 × U/(2L)`；男性 non-PAR X `100 × U/L` |
| 同男性 X/常染色体比较 | 同一参考、同一男性集合，要求全部 22 条常染色体和 X 完成 |

未分相 IBDmix 片段并集不能区分同一位置一份或两份 archaic 拷贝，这些统计不能直接当成已校准的 ancestry dosage。物理跨度也不是 callable mask 长度。成功完成且零片段记为零；未检测或未完成保留 missing。新 Shiny 表格与下载文件明确区分单位、实际样本数和参考，不把不同分母的量混为同一个百分比。

## 6. ARG → TRACE

```bash
./arg.sh build --method threads --format trace \
  --dir-gen "$GU_REF_ROOT/1kg/37" \
  --dir-map "$GU_REF_ROOT/hm3/37/gmap" \
  --sample-panel "$GU_REF_ROOT/1kg/38/samples.txt" \
  --dir-arg "$GU_REF_ROOT/1kg/37/arg.threads" \
  --dir-gen4arg "$GU_REF_ROOT/1kg/37/gen4arg.threads" \
  --grch 37 --chr 1-22,X --threads 8 --jobs 1 --memory-limit-gb 32

./arg.sh check --method threads --format trace \
  --dir-gen "$GU_REF_ROOT/1kg/37" \
  --dir-map "$GU_REF_ROOT/hm3/37/gmap" \
  --sample-panel "$GU_REF_ROOT/1kg/38/samples.txt" \
  --dir-arg "$GU_REF_ROOT/1kg/37/arg.threads" \
  --dir-gen4arg "$GU_REF_ROOT/1kg/37/gen4arg.threads" \
  --grch 37 --chr 1-22,X

./gu.sh trace --grch 37 --target 1kg --chr 1-22,X \
  --dir-target "$GU_REF_ROOT/1kg/37/pfile/chr" \
  --sample-panel "$GU_REF_ROOT/1kg/38/samples.txt" \
  --dir-arg "$GU_REF_ROOT/1kg/37/arg.threads/trace/threads" \
  --dir-out "$GU_ANALYSIS_BASE/1kg-trace2026" --memory-cap 64G
```

实际 TRACE 输入目录以 ARG 完成日志为准，必须直接包含树、sample map 与 build 记录。显式样本表会用于 ARG 的群体信息；需要另一份 ancestry 文件时使用 `--ancestry-file`。ARG scratch 放原生 Linux `/tmp` 或指定高速 scratch，主结果仍为持久目录。

所有选定常染色体进入一个联合 TRACE HMM 请求；男性 X 根据自己的样本映射与倍性单独建模。不要把 IBDmix 的每批 10,000 人规则直接套到既有全人群 ARG 上。TRACE 汇总的单一对外片段表为 `final/trace_haplotype_segments.tsv.gz`。

`TRACE_GENETIC_MAP_DIR` 若显式设置，每条所需染色体必须有恰好一个合法地图；未设置则记录默认换算 1 cM/Mb。ARG 图谱与 TRACE 的此参数分开配置。输入树、sample map、地图、mask、参数和实际 TRACE 代码均参与本版运行的缓存校验。汇总使用 NPZ 副本，保留 inference 输入完整性。

`phyml-GPU`、`trace-GPU` 是独立实验后端；本版 `--GPU` 参数只用于 IBDmix。ARG 无通用 GPU 入口。没有经过 GPU 硬件测试的结果应先在相同真实小规模数据上与 CPU 核对。

## 7. PhyML 和结果查看

```bash
./gu.sh phyml \
  --loci "$GU_DATA_ROOT/gwas/trait.jma.cojo" --loci-format cojo --grch 38 \
  --target 1kg --dir-target "$GU_REF_ROOT/1kg/37/pfile/chr" \
  --sample-panel "$GU_REF_ROOT/1kg/38/samples.txt" \
  --dir-out "$GU_ANALYSIS_BASE/trait-phyml2026" --jobs 1 --memory-cap 32G

./gu.sh final --dir-out "$GU_ANALYSIS_BASE/ukb001"
./gu.sh shiny --dir-out "$GU_ANALYSIS_BASE/ukb001"
```

COJO 的 `--grch` 描述 lead 输入坐标；风险位点序列分析使用 GRCh37，38→37 时需要 chain 与 liftOver。PhyML 默认用五个当前参考；IBDmix 的 `--archaic-gen` 不是 PhyML 参数。树证据、IBDmix 片段与 TRACE 推断分别保存，不能只凭一棵树的接近关系宣称全基因组渗入。

结果根中保留原生片段/压缩归档、`run.meta.tsv`、样本 roster 和 provenance；`final/` 保存 SQLite、汇总工作簿与当前 PhyML 报告。SQLite 从当前原始结果原子重建，schema 固定 2026，Shiny 只读本版结构。归档恢复与新运行内部的断点续跑保留，只有完整 contract、实际输入和完成收据匹配才复用；没有从未知文件自动认领完成状态的路径。

需要在本版运行中更改科学参数或输入时，应选新的输出根，或明确使用对应 `--replace-ibdmix TRUE`、`--replace-trace TRUE`、`--replace-phyml TRUE`。CPU 和 GPU 的结果验证应各写一个根目录，不同时写同一单元。

## 8. HPC 与验证

HPC 包使用相同科学核心，增加 LSF 的计划、冻结输入、数组提交、精确依赖、状态、失败续跑与汇总。部署和命令详见 `HPC.README.md`。工作流为 `plan` → 核对生成命令 → `submit`；单个 CPU 任务强制 `--GPU FALSE`。队列采用现有 `0HPC/ld/` 的方式配置，实际权限和配额以服务器查询为准。

验证清单与复现命令见 [tests/README.md](tests/README.md)，本次整理说明见 [CHANGES-2026.md](CHANGES-2026.md)，精确测试数量在 ZIP 顶层 `RELEASE.json`。验证包括官方源码算法对照、真实 TRACE 合成 ARG 流程、新数据库、公开命令生成、LSF 模拟和语法检查。没有重跑完整 1KG/UKB，没有实测 GPU 设备或真实 LSF 提交，也没有在完整 R 环境中交互测试 Shiny。
