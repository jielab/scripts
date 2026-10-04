# GWAS

三个根目录入口与辅助代码一一对应。默认参数及命令示例见各入口的 `--help`。

| 模块 | Shell | Python | R |
|---|---|---|---|
| GWAS | `gwas.sh`、`f/gwas.f.sh` | `f/gwas.py` | `f/gwas.inputs.R` |
| 比较、LDSC | `compare.sh`、`f/compare.f.sh` | `f/compare.py` | `f/compare.R` |
| 格式化、后处理 | `format.sh`、`f/format.f.sh` | `f/format.py` | `f/format.R` |

`gwas.py merge` 合并 PLINK/REGENIE/SAIGE 结果，其余子命令负责计划、执行和缓存。`compare.py` 包含 compare/ldsc/shiny 流程，以及 ldsc-run、ldsc-summary、ld-blocks 工具。`format.py` 包含 liftover、h2、magma-ids、yap2018-mpb。`compare.R` 使用 compare/ldsc 模式；`format.R` 使用 thin/mplot 模式。

Shell 格式化流程、任务生成、执行及 BGZF/tabix 索引函数集中在 `f/format.f.sh`。该文件供主入口和生成的任务加载；`--base`、`--worker` 是内部加载模式。LE8 通过 `source .../format.f.sh --index` 只加载索引函数，不启动格式化流程。

三个 AWK 文件按归属命名为 `gwas.extract.awk`、`format.hm3.awk`、`format.thin.awk`。Shiny 仍在 `shiny/`；临时验证和测试放在 `/tmp`。没有旧文件名转发脚本。
