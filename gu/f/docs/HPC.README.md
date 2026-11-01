# 本机部署补充（2026-10-31）

实际入口为 /work/sph-huangj/scripts.HPC/gu/gu-hpc.sh，共享代码全部位于 f/。
本部署不需要外部 0f/、grid/f/ 或额外 gu/ 层。原共享 grid common 文件改名为 grid.common.py / grid.common.sh，避免覆盖 GU 同名文件。
默认队列 short；每个阶段至少申请 40 slots，实际科学线程仍由 --threads 控制。
代码身份检查仅覆盖 GU 入口和 f/、shiny/、tests/，不会冻结旁边的 ld/、metag/、ukb/。
已配置集群现有 py311 环境作为调度与基础运行环境；完整科学依赖仍须安装，见 DEPLOYMENT.txt。
数据尚待复制，不能把代码验证当作真实数据计算成功。文件清单见 COPY-FILES.tsv / COPY-GUIDE.txt。

以下为原包使用说明（路径与队列已更新）；原 ZIP 的三目录解压步骤不适用于这个已合并部署目录。

# GU FINAL 2026：南方科技大学 LSF / CPU 版

交付文件为 `gu-FINAL-HPC-2026.zip`，从新的代码和新的结果目录开始。`gu.sh`、`arg.sh` 执行科学分析，`gu-hpc.sh` 负责 LSF 调度。IBDmix 使用同一套 CPU 科学实现；本版不会请求 GPU。`gu.sh ibdmix --GPU FALSE` 可在一个已经分配的计算节点内直接运行。

## 1. 部署与必须提供的数据

将 ZIP 内的 `gu/`、`0f/` 和 `grid/f/` 直接部署在 `/work/sph-huangj/scripts.HPC/` 下，保留这三个目录的相对位置。可执行入口是 `/work/sph-huangj/scripts.HPC/gu/gu-hpc.sh`。不要只复制 `gu.sh`。

用户已经将 F: 盘内容对应复制到 `/data/sph-huangj/`，因此输入和输出采用以下固定映射：

| 内容 | HPC 路径 |
| --- | --- |
| 代码及同级辅助程序 | `/work/sph-huangj/scripts.HPC/` |
| F: 的 `gen/` | `/data/sph-huangj/gen/` |
| F: 的 `analysis/` | `/data/sph-huangj/analysis/` |
| 软件与独立环境 | `/work/sph-huangj/software/gu/` |

所有示例都使用这些部署路径。每次正式分析建立一个新的 `--run-dir`，本次运行内部可用 completion receipt 可靠续跑。

```bash
cd /work/sph-huangj/scripts.HPC
cp hpc.runtime.env.example hpc.runtime.env
```

编辑 `hpc.runtime.env` 一次，将示例路径改成服务器上的真实路径。该文件不是分析参数的猜测器，尤其需要确认以下输入。

| 项目 | 作用 |
| --- | --- |
| `GU_CONDA_BASE` / `GU_CONDA_EXE` | HPC 上的 conda 安装 |
| `GU_SOFT` / `IBDMIX_RUNTIME` | 已安装程序；`IBDMIX_RUNTIME/build/src/ibdmix` 应存在 |
| `GU_REF_ROOT` | `/data/sph-huangj/gen`，基因组参考根目录 |
| `GU_FILES_ROOT` / `GU_COJO_CHAIN` | D: 来源的补充文件和 hg38→hg19 chain，须额外复制 |
| `IBDMIX_MASK_ROOT` | 与古基因组、GRCh 版本相匹配的原始排除区域 mask |
| `GU_ANNOT_ROOT` / `IBDMIX_AXT_DIR` | `/data/sph-huangj/annot/axt/37` 的物种比对，供 CpG 排除区域计算 |
| `IBDMIX_REFERENCE_FASTA` | 正确版本的参考 FASTA 及其索引 |
| `UKB_REF_FASTA` / `GU_REF_FASTA` | UKB 准备及共享流程所用参考 FASTA，37/38 必须相符 |
| `CHECK_GRCH_SNP_LIST` | `/data/sph-huangj/ukb/common/snp.lst`，跨 build 坐标校验 SNP 表；这是数据，不在代码 ZIP 内 |
| `PATH` | `bcftools`、`samtools`、`plink2` 等软件所在目录 |
| `GU_THREADS_ENV` | 使用 Threads ARG 时的独立环境 |

包内的 `f/runtime.env` 自动提供本部署的路径默认值，并允许显式环境变量覆盖；作业专用 `hpc.runtime.env` 负责激活实际软件环境。安装脚本、`gu.sh` 和 `arg.sh` 使用同一套默认值。

包内不含 UKB/1KG 基因型、古基因组、FASTA、遗传图谱、mask 或 SNP 校验表。F: 盘的数据镜像提供主要输入。**SNP 校验表之前位于 D: 盘的 `data/ukb/common/snp.lst`，不属于已声明的 F: 镜像；若尚未额外复制，请将真实文件放到 `/data/sph-huangj/ukb/common/snp.lst`。** 若使用需要 hg38→hg19 转换的 COJO 输入，还需把 D: 的 `files/liftOver/hg38ToHg19.over.chain.gz` 额外复制到 `/data/sph-huangj/files/liftOver/hg38ToHg19.over.chain.gz`。不要用空文件绕过检查。输入 VCF/BCF 必须事先建立 `.tbi`/`.csi` 索引，以免数百个作业同时修改同一个索引。

软件尚未安装时，可先设置安装路径，然后运行 CPU 安装；安装需在集群允许联网/安装软件的位置进行。不要在运行中的 array 内同时安装。

```bash
export GU_CONDA_BASE=/share/apps/anaconda3/2020.7
export GU_CONDA_EXE="$GU_CONDA_BASE/bin/conda"
export GU_SOFT=/work/sph-huangj/software/gu
bash install.sh --GPU FALSE --dir-software "$GU_SOFT"

# 只有需要重新构建 ARG 时才安装独立 ARG 环境。
export GU_THREADS_ENV="$GU_SOFT/venvs/gu-threads"
export GU_SINGER_ENV="$GU_SOFT/venvs/gu-singer"
"$GU_CONDA_EXE" create -n gu-arg-bootstrap python=3.12 -y
export GU_THREADS_PYTHON="$GU_CONDA_BASE/envs/gu-arg-bootstrap/bin/python"
bash install.sh --arg
```

已有环境时直接使用。Threads 安装需要 Python ≥3.12，主 GU 环境的 Python 3.11 不能作为其安装解释器；上面的 `gu-arg-bootstrap` 只用于建立独立 Threads 环境。最终以本包主 README 中的安装与参考数据要求为准。

## 2. UKB：每 10,000 行 PSAM 一批，再筛选男性

一条计划覆盖所有样本批次及指定染色体，无需手写几十条命令。`--batch-size 10000` 按 PSAM 数据行分批，`--keep-males` 在分批后筛选男性。

```bash
./gu-hpc.sh plan ibdmix \
  --run-dir /data/sph-huangj/analysis/gu/ukb-FINAL-2026-001 \
  --env-file ./hpc.runtime.env \
  --queue short \
  --threads 4 --memory-gb 32 --walltime 168:00 \
  --max-running 64 --array-size 500 \
  --grch 37 --target ukb \
  --dir-target /data/sph-huangj/gen/ukb/37/imp/chr \
  --psam /data/sph-huangj/gen/ukb/37/imp/chr1.psam \
  --batch-size 10000 --keep-males \
  --chr 1-22,X \
  --archaic-path /data/sph-huangj/gen/archaic/37/vcf \
  --archaic-gen Altai,Chagyr,Denisova,Denisova25,Vindija

# 查看生成的完整命令；这一步尚未提交或运行任何分析。
less /data/sph-huangj/analysis/gu/ukb-FINAL-2026-001/jobs/commands.tsv

# 明确提交后，才调用 bsub。
./gu-hpc.sh submit /data/sph-huangj/analysis/gu/ukb-FINAL-2026-001
```

`ser` 是你已有 `0HPC/ld/` 示例使用的队列名，不代表所有账户当前都有权限；提交前在服务器用 `bqueues -l short` 查看实际队列、最大运行时间、内存和核数限制。需要其他队列时在 `plan` 中替换 `--queue`。

分批规则是先按原始 PSAM 数据行取第 1–10,000、10,001–20,000……行，再在每批中筛选 `SEX=1`。第一批不会被悄悄改成“前 10,000 个男性”。每个选定 IID 只属于一个批次。冻结后的 `inputs/batch00001.keep` 与原 PSAM 哈希保存在计划中；之后更改样本表会拒绝执行，需建立新计划。完全没有男性的批次有明确排除记录，不作为失败作业，也不伪装成有样本的零检出。

**每 10,000 行构成一个独立的 IBDmix 频率估计人群。** 这种并行方式不保证与把全 UKB 放在一个共同频率模型中的结果逐位相同。后台没有用每批的百分比去推算全队列比例。CPU/GPU 结果比较也应固定同一批人、同一参考、同一阈值和同一输入。

每个 `(批次, 染色体)` 都有独立的完整结果目录，例如：

```text
ukb-FINAL-2026-001/tasks/batch00001/chr1/
ukb-FINAL-2026-001/tasks/batch00001/chrX/
ukb-FINAL-2026-001/tasks/batch00002/chr1/
```

每个目录仍含正常 GU 的 `ibdmix/`、原始结果/归档和 provenance。不会让上千个作业并发改写同一个 cohort metadata、命令文件或临时导出目录。作业内部固定 `--GPU FALSE --foreground TRUE --jobs 1 --auto-final FALSE`。

每个科学任务还显式传入 `GU_BATCH_ID=batch00001` 等批次标签，供原生 `run.meta` 和后续 SQLite 保留批次身份；同一批的常染色体和 X 使用相同标签。汇总 TSV 的 `hpc_batch=1` 与 `batch00001` 对应。

## 3. 1KG：完整人群、五个当前参考

1KG 全体样本无需分批；**不提供 `--psam`** 即可。下面每条染色体是一个任务，每个任务保持完整人口结构及 AFR 对照。默认参考固定为 `Altai,Chagyr,Denisova,Denisova25,Vindija`。

```bash
./gu-hpc.sh plan ibdmix \
  --run-dir /data/sph-huangj/analysis/gu/1kg-FINAL-2026 \
  --env-file ./hpc.runtime.env --queue short \
  --threads 4 --memory-gb 32 --max-running 23 \
  --target 1kg --grch 37 --chr 1-22,X \
  --dir-target /data/sph-huangj/gen/1kg/37/pfile/chr \
  --sample-panel /data/sph-huangj/gen/1kg/38/samples.txt \
  --archaic-path /data/sph-huangj/gen/archaic/37/vcf \
  --archaic-gen Altai,Chagyr,Denisova,Denisova25,Vindija

./gu-hpc.sh submit /data/sph-huangj/analysis/gu/1kg-FINAL-2026
```

可选的 `Altai.2013`、`Denisova.2013` 仅在用户显式传入 `--archaic-gen` 并提供对应原始数据时启用，不自动加入默认分析。UKB 没有 AFR5 控制人群时关闭的 African Denisovan subtraction，不应被写成用于所有数据集的全局环境配置。

chrX 仍由 GU 的 nonPAR、男性、实际调用区域与样本数规则处理；HPC 不自行把男性 X 复制成两条染色体，也不人为调高 X 结果以逼近常染色体。

## 4. TRACE：已经有 ARG

TRACE 的常染色体要联合训练 HMM。这里整个所选常染色体集合进入**同一个 TRACE 请求**，程序内部可按样本/节点组并行。不会因为 HPC 有很多节点，就把每条常染色体分开拟合后说成相同的方法。X 的样本/倍性差异由本包 TRACE 驱动分组处理。

```bash
./gu-hpc.sh plan trace \
  --run-dir /data/sph-huangj/analysis/gu/1kg-TRACE-FINAL-2026 \
  --env-file ./hpc.runtime.env --queue short \
  --threads 16 --memory-gb 128 --walltime 168:00 \
  --target 1kg --grch 37 --chr 1-22,X \
  --dir-target /data/sph-huangj/gen/1kg/37/pfile/chr \
  --sample-panel /data/sph-huangj/gen/1kg/38/samples.txt \
  --dir-arg /data/sph-huangj/gen/1kg/37/arg.threads/trace/threads

./gu-hpc.sh submit /data/sph-huangj/analysis/gu/1kg-TRACE-FINAL-2026
```

ARG 目录需要直接含 `chrN.trees` / sample map / `ARG_TRACE_BUILD.tsv`，不是任意的树文件父目录。TRACE 计划拒绝 `--psam` 和 `--keep-males`；不能把针对 IBDmix 的 10,000 人批处理直接套到全人群 ARG/HMM 训练上。

## 5. 先并行构建 ARG，再联合运行 TRACE

```bash
./gu-hpc.sh plan arg-trace \
  --run-dir /data/sph-huangj/analysis/gu/1kg-ARG-TRACE-FINAL-2026 \
  --env-file ./hpc.runtime.env --queue short \
  --threads 8 --memory-gb 96 --max-running 23 \
  --target 1kg --grch 37 --chr 1-22,X \
  --dir-gen /data/sph-huangj/gen/1kg/37 \
  --dir-map /data/sph-huangj/gen/hm3/37/gmap \
  --dir-target /data/sph-huangj/gen/1kg/37/pfile/chr \
  --sample-panel /data/sph-huangj/gen/1kg/38/samples.txt \
  --arg-method threads

./gu-hpc.sh submit /data/sph-huangj/analysis/gu/1kg-ARG-TRACE-FINAL-2026
```

顺序为：准备检查 → 按染色体的 ARG array → 核对所有 ARG receipt 并组装 `arg-assembled/` → 一个联合 TRACE 请求 → 汇总。

`--dir-map` 显式指定与 GRCh 版本相符的 ARG 遗传图谱目录；未指定时默认 `/data/sph-huangj/gen/hm3/GRCH/gmap`，其中 GRCH 取 `--grch`。TRACE 自身的遗传图谱设置遵循主 README，与 ARG 图谱参数分开。

1KG 示例复用 `38/samples.txt` 的人群和性别元数据；这张表不包含基因组坐标，因此可用于 GRCh37 数据。冻结后的同一张 panel 同时提供给 Threads/Needle 的样本选择和 TRACE，避免两阶段读取不同人群表。

每条染色体使用自己的 `arg/` 和 `prepared/`，不会抢占同一个 ARG 目录锁。组装目录仅链接经核对的输出，保留 method、build、sample map 与每条染色体的 QA 文件。TRACE 继续核对树的样本身份、坐标、突变与覆盖范围。

`HPC_ARG_ASSEMBLY.json` 记录原始计划及所有组装文件的身份。HPC 驱动核对各染色体完成记录后才将它交给 TRACE，避免普通单根 ARG 检查错误地把组装目录的上两级当作原始 ARG 工作目录。之后再次使用这个组装目录时仍会核对原计划与输出；普通外部 ARG 目录不会跳过其常规检查。

只要 ARG 时把 `plan arg-trace` 改为 `plan arg`；成功后的 `arg-assembled/` 可用于之后的 TRACE 计划。`--arg-method` 可选 `threads`、`needle`、`tsinfer`，各自的程序、输入和地图要求遵循主 README；示例采用 Threads。区域 SINGER 小测试不在此处伪装成全染色体 ARG。

## 6. PhyML：给定 GWAS loci

```bash
./gu-hpc.sh plan phyml \
  --run-dir /data/sph-huangj/analysis/gu/gwas-PhyML-FINAL-2026 \
  --env-file ./hpc.runtime.env --queue short \
  --threads 8 --memory-gb 64 \
  --target 1kg --grch 37 --chr 3,12 \
  --dir-target /data/sph-huangj/gen/1kg/37/pfile/chr \
  --sample-panel /data/sph-huangj/gen/1kg/38/samples.txt \
  --archaic-path /data/sph-huangj/gen/archaic/37/vcf \
  --loci /data/sph-huangj/gwas/trait.jma.cojo
./gu-hpc.sh submit /data/sph-huangj/analysis/gu/gwas-PhyML-FINAL-2026
```

`--chr` 应列出 loci 实际涉及的染色体，用于冻结输入检查；PhyML 科学请求由原始 COJO leads 决定。该 HPC 入口把整份固定 loci 请求放进一个 job，沿用程序的缓存与树推断；不会对几十万个样本随意分批后把不同树混为一棵树。

## 7. 监控、失败续跑、输出

```bash
./gu-hpc.sh status /data/sph-huangj/analysis/gu/ukb-FINAL-2026-001
bjobs -l JOBID
bpeek -f JOBID

# 等前一轮相关作业结束，只重跑明确失败的任务编号。
./gu-hpc.sh retry /data/sph-huangj/analysis/gu/ukb-FINAL-2026-001 --indices 7,12-15

# 科学任务都成功，仅 TRACE 后续阶段或汇总需要重跑。
./gu-hpc.sh retry /data/sph-huangj/analysis/gu/ukb-FINAL-2026-001 --indices none
```

`status` 显示的是最后一次 worker receipt；进程被 SIGKILL 后可能仍显示 running，此时用 `bjobs` 看调度器状态。`retry` 会查询前次 JOBID，发现仍有运行/排队作业即拒绝重复提交。初次 `submit` 不可重复执行；部分提交失败的 JOBID 也会即时落盘，不会假装整次提交未发生。

任务依赖用**本次提交返回的确切 JOBID**，不使用 `done("gu.*")` 这样的作业名通配符。使用 `ended(JOBID)` 让失败后的下游任务也能启动并报告“上游未完成”，其科学动作前会严格核对完成 receipt。这避免失败 array 留下一串永远 PEND 的合并任务；它不意味着允许汇总失败输出。

| 路径 | 内容 |
| --- | --- |
| `plan.json` | 冻结参数、代码哈希、输入 metadata、批次和完整任务表 |
| `inputs/` | 冻结环境、样本批次、panel/COJO 文件 |
| `jobs/commands.tsv` | 可审核的实际科学命令 |
| `jobs/*.cmd` | 需要 LSF 分配后才能执行的脚本 |
| `submissions.jsonl` | 每次提交状态、JOBID、原始 bsub 参数 |
| `logs/*.out` / `*.err` | LSF 日志，含 `%J/%I`，不会覆盖别的任务 |
| `state/task.N.json` | running/failed/complete 与最终输出哈希 |
| `tasks/` | 每个隔离的 GU/ARG 原生结果 |
| `summary/segments.tsv.gz` | 保留 `hpc_task`、`hpc_batch` 的合并片段表 |
| `summary/task_summary.tsv` | 各任务/批次的片段数、行长度之和、结果位置 |
| `summary/COMPLETE.json` | 只有所有必需任务和输出内容核对通过才产生 |

**有效表头但没有片段，是成功的零检出；没有输出、退出非零、receipt 不匹配或输出内容被修改，是失败。** 后者不会以零替代。存在原生归档和解压副本时，合并只读取一份。

`batch_cohort_samples` 是该批次进入染色体筛选之前的人数，X 的实际调用样本数可能不同；这个字段不是百分比的分母。`sum_segment_bp` 是结果行长度相加，同一人的不同古参考可能重叠，不能读成 genome union coverage，也不能读成 archaic allele dosage。原生 GU 结果里的具体统计与科学解释仍以主 README 为准。

汇总不会输出未经分母校准的“全 UKB introgression 百分比”，不会平均不同批次的 population union 百分比。若查看单个任务的原生 Shiny/报告，使用该任务的完整结果根；本 HPC 合并 TSV 是跨任务数据交付，不是把多批不同频率模型伪装为同一个 Shiny cohort。

## 8. CPU、内存、任务规模和缓存

`--threads N` 为一个 worker 的程序线程预算，`bsub -n` 申请 `N+1` 个 slots，预留监督进程。`span[hosts=1]` 保证一个 task 的共享内存程序运行在同一节点。GU 染色体 job 并行固定为 1，避免 array × 本地多 job 的隐式放大；TRACE 内部按配置使用多个 worker，每个数值库限定单线程。

`--memory-gb M` 同时传入 GU/ARG 内存预算并生成显式单位的 `bsub -M MGB`。LSF 对 `-M` 是每进程或整作业的强制方式取决于管理员配置，不应只凭 `-M` 推断硬限制已开启。`rusage[mem]` 也可能按 slot 或 job 计算：默认 `--rusage-memory-mode per-job` 直接申请 M GiB，对按 slot 计费的系统会保守地多预留；确认站点按 slot 解释后可改为 `--rusage-memory-mode per-slot`，总预留量至少 M GiB。不要在不清楚站点设置时把申请量除以线程数。

`--array-size 500` 控制每个 array 的最大 index；按实际 LSF 限制调整。超过时分成多个 array，`--max-running` 的并发上限分摊到它们，总和不会超过设置。例如 1,150 个任务、每 array 500、总并发 64，形成 3 个 array，运行额度分为 22/21/21。

计划创建后不应编辑代码、冻结环境或样本选择；worker 检出变更会停止，避免同一运行混杂不同科学实现。科学失败重试仍会进入原程序的输入/cache 校验，不仅凭“某文件存在”跳过。大基因型输入的计划身份记录 size/mtime，小型 PSAM/panel 等有 SHA256；不声称在登录节点重新计算全部 TB 级 PGEN 的哈希。

共享参考的 native prepare 编译在单独的准备 job 完成；各 worker 不创建全局安装。ARG 暂存按任务隔离，长期结果不要写入集群有自动清理政策的 `/scratch`。本包保留任务暂存以便失败恢复；确认成功并保留原生最终结果后，可按主 README 使用 cleanup 工具，不要在其他任务仍运行时做全根目录清理。

## 9. 验证范围与依据

包内 `tests/test_hpc.py` 使用合成小表、压缩结果和 mock LSF 验证计划、批次、资源参数、精确 JOBID 依赖、拒绝重复/活动任务重试、部分失败阻止汇总、零检出、归档去重和输入变更。可运行：

```bash
python3 -m unittest discover -s tests -p 'test_hpc.py' -v
```

这些检查没有向南科大 HPC 提交实际作业，也不等于在服务器上重跑了 1KG/UKB。第一次完整运行前，建议先用单条常染色体和一批样本验证真实路径、软件和资源，再建立正式运行计划。

参考实现是 `jielab/scripts/0HPC/ld/` 的 `.cmd + bsub + merge dependency` 组织方式；本版增加独立命名、数组、完成 receipt 与严格聚合。官方说明截至 2026-10-31 可检索到的版本包括：

- 南科大 [Quick Start for Supercomputers，2025.9](https://hpc.sustech.edu.cn/ref/QuickStartforClusterPlatformofCCSE.pdf)：使用 LSF、避免登录节点长期计算、存储与自动清理政策。
- 南科大 [超级计算机使用手册，2025.3](https://hpc.sustech.edu.cn/ref/cluster_User_Manual.pdf)：`bqueues`、`bsub`、`bjobs`、`bpeek`。
- IBM [bsub -M 官方说明](https://www.ibm.com/docs/en/spectrum-lsf/10.1.0?topic=o-m)：显式 KB/MB/GB 单位及内存限制的站点依赖行为。

本次直接打开南科大站点时返回校园网访问/维护页面，因此没有把队列名、配额或最大 array 规模宣称为已实时核实的配置。它们均保留为可配置参数。
