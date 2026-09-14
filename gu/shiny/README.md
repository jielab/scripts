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
各参考使用其发布者的独立质量掩码；原始资源缓存统一放在 `I:/refGen/ibdmix/37/masks/raw`，供所有分析模式共用。
`IBDMIX_MASK_CACHE` 可覆盖缓存根目录（其下使用 `raw` 保存原始资源）。
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
