# GU Shiny

启动入口仍是项目根目录的 `gu.sh`：

```bash
cd /mnt/d/scripts/gu
./gu.sh shiny
```

默认地址为 `http://127.0.0.1:3838`；主脚本继续管理数据路径、环境设置和密度缓存准备。可用 `GU_SHINY_HOST`、`GU_SHINY_PORT` 修改监听地址。

- `app.R`：启动应用，载入各模块，注册参考序列请求。
- `ui.R`、`server.R`：页面布局与服务器交互。
- `data.R`：共享数据访问与资源初始化。
- `methods.R`、`density.R`、`summary.R`：各分析方法和汇总面板。
- `phyml.R`、`dual_lead.R`：PhyML 报告与双 lead 面板。
- `ref.R`：IGV 本地参考序列服务。
- `review.R`：已有静态 review 页面的小型 Shiny 宿主。
- `www/`：JavaScript、CSS、地图数据和 IGV 静态资源。

五参考分析默认使用 Altai、Chagyr、Vindija、Denisova、Denisova25。
PhyML 在同一棵树上分别检验 Neanderthal（三参考）和 Denisovan（两参考），
报告、携带者验证和图像均按谱系分开。加入参考后需要重新建树；旧的三参考树不能充当 Denisovan 结果。

IBDmix 默认 `IBDMIX_PROFILE=multi_reference`，对五个参考分别在各人群运行。
Neanderthal 保留非洲 Denisova 对照过滤；两个 Denisovan 参考输出独立的原生匹配片段，
不减去由 Denisova 自身产生的对照区间。Overview 按个体合并同谱系参考的重叠片段后计算覆盖率，
不将两个参考的百分比相加。原生匹配本身不等于已确认的渗入来源。
世界地图下方的 Population summary 横跨整行，分别展示五参考的 `Mb / ind.`：
各参考独立合并每个个体在已确认完成的常染色体上的重叠片段，再按个体等权平均，不跨参考相加。
另有 `Altai + Denisova`（不含 Denisova25）和 `All Five` 两列：先按个体合并所选参考的片段，
重叠位置只计一次，再按个体等权平均；仅使用组合内所有参考均已确认完成的常染色体。
缺少任一成员参考的检测结果时，不将已检测参考的数值冒充完整组合结果。
地图、Coverage、右侧比例和 Cell 2020 对照均仅使用 Altai；旧表中的 `23.2` 等数值也是 Altai 单参考。
未检测的参考显示 `N/A`，已检测但无片段为 `0.0`；人数不同或染色体不完整时另行标注。
升级后运行 `./gu.sh shiny` 会自动生成五参考及两种合并统计的汇总缓存，不需要重新运行 IBDmix。
各参考使用其发布者的独立质量掩码，统一放在 `E:/refGen/archaic/37/mask/<参考名>/`。
Altai 和 Denisova 的 minimal filters 也在各自目录，保留原始长文件名，不能以 `chrN_mask.bed.gz` 替代。
`mask/common/` 保存 1KG strict accessibility 和 genomicSuperDups 原始注释。
参考资源目录中的 `mask/` 只保存永久输入文件，不写入运行时派生文件。
分析输出目录的 `mask/derived/cpg/` 保存 CpG 中间计算结果；`mask/derived/combined/` 保存依赖现代样本、参考和参数生成的排除区域。
例如默认 1KG chr1 分析写入 `/mnt/d/analysis/gu/ibdmix/1kg/chr1/mask/derived/`。
这些派生文件可在分析不再使用时清理，原始输入齐全时可重新生成；运行中的任务和指向它们的分析链接仍依赖这些文件。
原始 mask 表示可纳入区域，派生的 `*.exclude.bed` 表示排除区域，不能互换。
分析目录使用单数 `mask/<分析单元>/`，链接到同一分析输出目录 `mask/derived/` 中的派生结果。
`E:/annot/axt/37/vs*` 仅保存跨物种比对原始文件及其校验、下载文件。
`IBDMIX_AXT_DIR`（默认 `${GU_ANNOT_ROOT:-/mnt/e/annot}/axt/37`）指定 AXT 根目录。
`IBDMIX_MASK_ROOT` 仅覆盖永久 mask 输入目录，默认是古人类 VCF 目录同级的 `mask`；派生结果的位置随分析输出目录确定。
`--replace-ibdmix FALSE` 按运行参数、输入文件记录、程序版本和实际排除 BED 的 SHA-256 判断能否复用。
已记录的 `/mnt/i/refGen/` 到 `/mnt/e/refGen/` 迁移，以及仅改变存放路径的 manifest、工作流脚本时间戳，不使结果失效。
`f/ibdmix.sh` 的 `pipeline_version` 是工作流算法兼容版本；修改基因型处理、调用、过滤或结果语义时必须更新它，不能只改脚本而沿用旧版本。
真实输入、mask 内容、原生调用程序或算法兼容版本变化时，仍需显式替换旧结果；报错会列出不兼容记录。
`IBDMIX_MASK_CACHE` 已停用，设置它会报错并提示更换。`IBDMIX_MASK_DIR` 仍专门用于 custom profile 的排除区域 BED。
运行前统一检查本地 mask，缺失、空文件或带 `.aria2` 标记时退出，不自动下载 mask。
显式设置 `IBDMIX_PROFILE=cell2020` 可保留原来的 Altai + Denisova 对照流程。

质量掩码来源：[Chagyr](https://ftp.eva.mpg.de/neandertal/Chagyrskaya/FilterBed/)、
[Vindija](https://ftp.eva.mpg.de/neandertal/Vindija/FilterBed/Vindija33.19/)、
[Denisova25](https://ftp.eva.mpg.de/denisova/Den25/FilterBed/)。

批处理也使用的 `introgression_density.py`、`phyml_panel_b.R` 等公共代码继续位于 `../f/`。

单独展示已有的 review 页面：

```bash
Rscript --vanilla shiny/app.R --review \
  --data /mnt/d/analysis/gu/final/review --port 3839 --host 127.0.0.1
```

AXT 及其发布者 MD5 清单也只从本地读取，缺失或校验失败会报错；分析代码已移除自动下载实现。
