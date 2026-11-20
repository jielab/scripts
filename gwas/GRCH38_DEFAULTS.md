# GRCh38 defaults — 2026-11-16

所有目标基因型/主分析默认 GRCh38。入口示例显式写 `--grch 38`。
HPC/LOCAL 的根目录仍来自相邻的 `0f/0os.sh`，没有写死机器根目录。

GRID: `${DIR_GEN}/ukb/38/imp` → `${ANALYSIS_ROOT}/grid/38/`。
PCA: `${DIR_GEN}/ukb/38/pca_proj/`。
Consensus: `${DIR_LD}/38/consensus/`，已复制验证过的 GRCh38 BED 目录。

```bash
./1.pca.sh --grch 38
./2.csx.sh --grch 38
./3.disco.sh --grch 38
./4.prsformer.sh --grch 38
./grid.sh --grch 38
```

依次等待上一个独立入口的任务成功后再提交下一个。`grid.sh all --grch 38` 与最后一条等价；GRID 不自动运行前四个入口。

allele3 新任务默认输出 `${ANALYSIS_ROOT}/allele3/38/`。
古代 PCA 数据产品仍在 `${DIR_GEN}/ancient/38/pca3d/`，UKB 投影在 `ukb/38/allele3/`。

```bash
./1.ancient_pca.sh --grch 38
./2.herv_scan.sh --grch 38
./3.allele3_scan.sh --grch 38
./4.avi_scan.sh --grch 38
./allele3.sh --grch 38
```

现有扫描作业继续写入原路径，没有停止或重新提交。新版本不会把未完成的旧结果认作已完成的新版输入。
需要从原始古代面板准备 GRCh38 时，单独使用 `./0.ancient_split.sh --grch 38`；这是 37→38 转换工具，故必须读取真正的 37 源数据并保留 37 中间结果。

GWAS 默认读取 `${DIR_GEN}/ukb/38/{typ,imp}`，分析输出 `${ANALYSIS_ROOT}/gwas/38/<module>/`。

```bash
./gwas.sh regenie --grch 38 --pheno-name height --type qt
```

## 不可伪装成 GRCh38 的原始参考

- PRS-CS/CSx 官方 LD 资源 `${DIR_LD}/csx` 原始坐标为 GRCh37；按 rsID 和等位基因与 GRCh38 目标对齐，不把原始 BP 当作 BP38。
- allele3 的古代 HQ、原始 UKB MFI、FASTA37 及 liftOver chain 保留真实来源版本；古代时间分析显式转换坐标，MFI 比较使用转换后的 37 目标键。
- `hg38` 注释、文件名 `*.38.bed` 也是 GRCh38，不需要仅为了目录名移动原始资源。
- 表型、样本标签、亲缘图、冻结 TEST 名单不含基因组坐标。亲缘图从 `37/typ/ukb_rel.dat` 复制到 `38/typ/ukb_rel.dat`，校验和一致，不涉及 liftOver。保留原 TEST 名单以保持可比性。

显式自定义路径仍由调用者负责真实内容；路径检查不能证明一个任意命名文件的实际坐标版本。禁止把 37 文件简单改名为 38。
此次未提交任何生产分析任务。
