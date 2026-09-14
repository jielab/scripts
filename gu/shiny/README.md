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

批处理也使用的 `introgression_density.py`、`phyml_panel_b.R` 等公共代码继续位于 `../f/`。

单独展示已有的 review 页面：

```bash
Rscript --vanilla shiny/app.R --review \
  --data /mnt/d/analysis/gu/final/review --port 3839 --host 127.0.0.1
```
