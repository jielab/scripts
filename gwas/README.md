# GWAS

三个根目录入口与辅助代码一一对应。默认参数及命令示例见各入口的 `--help`。

| 模块 | Shell | Python | R |
|---|---|---|---|
| GWAS | `gwas.sh`、`f/gwas.f.sh` | `f/gwas.py` | `f/gwas.inputs.R` |
| 比较、LDSC | `compare.sh`、`f/compare.f.sh` | `f/compare.py` | `f/compare.R` |
| 格式化、后处理 | `format.sh`、`f/format.f.sh` | `f/format.py` | `f/format.R` |

`gwas.py merge` 合并 PLINK/REGENIE/SAIGE 结果，其余子命令负责计划、执行和缓存。`compare.py` 包含 compare/ldsc/shiny 流程，以及 ldsc-run、ldsc-summary、ld-blocks 工具。`format.py` 包含 liftover、h2、magma-ids、yap2018-mpb。`compare.R` 使用 compare/ldsc 模式；`format.R` 使用 thin/mplot 模式。

Shell 格式化流程、任务生成、执行及 BGZF/tabix 索引函数集中在 `f/format.f.sh`。该文件供主入口和生成的任务加载；`--base`、`--worker` 是内部加载模式。LE8 通过 `source .../format.f.sh --index` 只加载索引函数，不启动格式化流程。

`format` 默认 `--rsid TRUE --rsid-unmatched drop`：已有 rsID 保留；其他 ID 按声明的 GRCh 版本、位置和两条等位基因匹配 dbSNP，不改变效应方向。无法唯一匹配的完整标准化行写入 `common/<trait>/qc/<trait>.rsid.unresolved.tsv.gz`，统计写入 `<trait>.rsid.tsv`。`--rsid-unmatched keep` 保留未匹配的原 ID；`--rsid FALSE` 完全跳过转换。更改筛选策略时可用 `--replace TRUE` 从原始数据重新生成。

`--hm3 TRUE` 先保留 HM3 rsID/坐标或 `P < --p-hm3` 的候选，再转换 ID，并以转换后的 HM3 rsID 或 P 值确认最终保留项。默认阈值为 `1e-3`，不是全基因组显著阈值 `5e-8`。dbSNP 的 BGZF/tabix 副本和按位置查询的结果共用 `<project>/.project/rsid/`，不同 GWAS 只补查新坐标，包含未匹配结果的缓存；首次使用每个版本需要准备参考索引。`--dbsnp FILE`、`--rsid-cache DIR` 可指定其他参考和共享缓存位置。

列识别函数是 `0f/phenotype.sh` 的 `phe_header_names()`，由 `gwas_clean_header_names()` 调用，支持 CKB 的 `EA_FREQ → EAF`、`AA/A2 → NEA`。`--sample-info FILE` 接受 `phenocode,num_samples,num_cases,num_controls` 或 `GWAS,N` 的 TSV，逐性状补充缺失 N，已有有效 N 保留。可选 `grch` 列用于检查显式版本，或在 `--grch auto` 时提供版本；样本量和版本冲突会在启动任务前报错。

CKB 已按原始版本拆成 `/mnt/f/gwas/ckb_37`、`/mnt/f/gwas/ckb_38` 两个项目，各自包含 `raw0/`、`raw/`、`ckb_decrypt.sh` 和本版本的 `phenotypes.tsv`。解密脚本接受空格或 Tab 分隔，输出到各项目自己的 `raw/[gwas].gz`。自动扫描 GWAS 时跳过 `decryption_report_*.tsv` 解密报告（含 `.gz`、`.bgz` 压缩版本），保留普通 TSV 数据输入。分析分别指定 `--dir-raw /mnt/f/gwas/ckb_37/raw --dir-out /mnt/f/gwas/ckb_37 --label ckb_37 --grch 37 --sample-info /mnt/f/gwas/ckb_37/phenotypes.tsv`，或对应的 `ckb_38`/`38` 参数；模块及执行参数为 `format,magma,lead,mplot --hm3 TRUE --run-cmd TRUE --foreground TRUE --jobs 4`，不含 PGS。输出分别位于两个项目的 `common/`、`mplot/`。

三个 AWK 文件按归属命名为 `gwas.extract.awk`、`format.hm3.awk`、`format.thin.awk`。Shiny 仍在 `shiny/`；临时验证和测试放在 `/tmp`。没有旧文件名转发脚本。
