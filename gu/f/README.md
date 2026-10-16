# GU 内部模块

日常分析从上一级的 `gu.sh`、`arg.sh` 启动。此目录按功能组织内部实现。

| 文件 | 用途 |
| --- | --- |
| `ukb_prepare.py` | UKB PGEN 检查、样本选择、固定队列、VCF 导出及质控、古人类评分输入准备 |
| `0.common.py progress` | 旁路观察 UKB 子进程与文件读取进度，报告当前读取的预计剩余时间，不改动准备缓存身份 |
| `ukb_prepare.sh` | UKB 原始 BGEN/芯片输入检查、VCF 导出和祖先等位基因注释；PGEN 操作交给同名 Python 模块 |
| `gwas_loci_prepare.py` | COJO 信号归并、基因组坐标转换、GRCh37 区间及审计表生成 |
| `introgression_review.py` | 渗入证据复核、原始 lead 与单倍型标签比较、HTML/TSV 报告 |
| `results_normalize.py` | 整合各方法结果到 SQLite、参考结果缓存和渗入密度汇总 |
| `arg.sh` / `arg.py` | ARG 内部工作流及数据转换、缓存工具；入口为上一级 `arg.sh` |
| `arg.py singer` / `arg.py threads` | SINGER、Threads 各自的 ARG 后端 |
| `ibdmix.sh` / `ibdmix.py` | IBDmix 调度、掩膜、群体样本、质控和片段汇总 |
| `ibdmix.py loci` | Shiny 的群体窗口支持峰汇总，保留人数分母与去重规则 |
| `ibdmix.py resume` | UKB 整染色体续跑：预处理前核验原始输入、样本内容、调用参数与结果完整性；临时样本文件重建不使已完成结果失效 |
| `ibdmix.prepare.cc` / `ibdmix.sh build-prepare` | 共享 HTSlib 预处理器及编译入口，一次读取现代基因型供多个古人类参考使用 |
| `ibdmix.gpu.cu` / `ibdmix.sh build-gpu` | IBDmix CUDA 片段调用器及编译入口 |
| `phyml.py` / `phyml.R` | PhyML 分析、结果检查及绘图 |
| `trace.sh` / `trace.py` | TRACE 分析及片段处理 |
| `as3.sh` / `as3.py` | ArchaicSeeker 3 分析及输入、结果处理 |
| `0.common.py` / `0.common.sh` / `0.common.R` | 跨方法共享的坐标、样本、路径、缓存和结果归档工具 |

UKB 的 Shell 和 Python 模块分别承担原始 BGEN/芯片与 PGEN 准备，保持独立以避免把两套输入流程混在一起。独立的位点准备、证据复核和结果整合也保留各自入口；重复的空命令分发层和模块加载已去除。

IBDmix 原生辅助文件直接平铺在 `f/`，使用 `ibdmix.` 前缀。构建与运行说明见 [ibdmix.md](ibdmix.md)，许可证见 [ibdmix.LICENSE](ibdmix.LICENSE)。已安装的可执行文件仍位于软件目录。

IBDmix 默认 `--replace-ibdmix FALSE` 时，UKB 整染色体任务先检查历史结果。完整且兼容的结果打印 `[GU CMD] SKIP`；部分完成的结果继续补跑。历史归档中的原始输入记录会保存为结果目录内的 `.ibdmix.resume.json`，用于避免重新预处理覆盖旧质控记录后丢失续跑依据。样本选择按实际队列及顺序校验，原始 PGEN、参考基因组、古人类 VCF、掩膜和计算参数仍参与检查；不兼容时在预处理前报明原因。该检查不改写旧调用或最终片段。

IBDmix 掩膜的旧文件 SHA 保留为审计证据；只有此次整合前的确切版本、且掩膜函数与常量的语法树一致时，续跑才接受旧记录。输入、样本、参数和片段校验仍按原规则执行。新掩膜记录当前代码的真实 SHA。
