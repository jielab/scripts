# EMS120 运行环境与数据准备

本文件针对新位置 `scripts/ems120`。调度入口和参数以 `0.README.md`、`bash ems120.sh --help` 为准。程序不会在提交任务时自动创建 Conda 环境或下载依赖；应先在启明上准备可用的环境，再提交 GPU 检查任务。

## 1. 两层 Python 环境

规划和提交作业只需要系统 Bash、Python 3.9+ 标准库，以及启明上的 `bsub`。生成计划不会导入 PyTorch、读取 EMS Excel、启动 R 或激活 Conda。计算节点上的实际分析使用 `--conda-env` 指定的环境，默认名称是 `ai`。

分析同时使用 R 和 Python。R/reticulate 必须绑定到该 Conda 环境的 Python；只在终端执行一次 `import torch` 成功，不能证明嵌入 R 后 CUDA 仍可用。入口的 GPU 检查因此在实际 R/reticulate 进程内执行。

### Python 依赖

| 组件 | 用途 / 要求 |
|---|---|
| PyTorch 的 CUDA 构建 | `ml_dx`、`ml_geo` 的 Transformer 训练和推理；CUDA 运行时应与节点驱动兼容 |
| `transformers`、`accelerate`、`safetensors` | Hugging Face 模型、Trainer、模型权重 |
| `numpy`、`pandas`、`scikit-learn` | 表格、指标、电话评分候选模型和交叉验证 |
| `openpyxl` | Python 读取训练 Excel |
| `joblib` | 电话评分模型缓存 |
| `xgboost >= 2.0`，支持 CUDA 的构建 | 电话评分 XGBoost 候选的 GPU 训练/推理；默认要求 GPU 分支可用 |
| `lightgbm` | 原程序中可选的电话评分候选；已有环境使用过它时，迁移环境也应保留 |
| `Pillow` | 共享 `png_output.py` 压缩和检查输出 PNG |
| `keplergl`（可选） | 启用 Fig2 交互地图导出时需要；当前提交用成图默认不生成 3D HTML |

这里没有把所有包强行升级到最新版本。迁移已验证的本地环境时，应记录本地 Python、R、PyTorch、CUDA、Transformers、XGBoost 和 scikit-learn 版本，在 HPC 建立兼容环境，并以小数据检查结果。已有 `joblib` 模型与生成它的库版本有关；迁移时保留版本记录和原始模型备份。

GPU 构建与设备选择的官方资料：

- PyTorch 安装选择器：<https://pytorch.org/get-started/locally/>
- XGBoost GPU：<https://xgboost.readthedocs.io/en/stable/gpu/>
- XGBoost Python API：<https://xgboost.readthedocs.io/en/stable/python/python_api.html>
- reticulate 的 Python 版本选择：<https://rstudio.github.io/reticulate/articles/versions.html>

### R 依赖

主脚本及直接使用的共享绘图函数需要：

```text
readxl writexl data.table tidyverse scales RColorBrewer reticulate
lubridate patchwork zoo broom forcats circlize jpeg png grid
jsonlite base64enc cowplot lattice
```

其中 `grid` 属于 R 自带组件，`lattice` 通常随推荐包安装；`tidyverse` 包含 dplyr、tidyr、purrr、tibble、stringr、ggplot2 等。不要把某个空的 `R_lib` 路径当成环境已经准备好的证明。缺少依赖时程序应给出包名并停止，由用户在计划之外准备依赖。

## 2. 共享脚本依赖

新版项目使用仓库同级的 `scripts/0f`，不再从 `/mnt/d/scripts/0f` 固定取文件，也不再调用强制 local 的 `local_paths.sh`。

需要保留的共享文件包括：

```text
scripts/0f/0os.sh
scripts/0f/phenotype.R
scripts/0f/plotting.R
scripts/0f/png_output.R
scripts/0f/png_output.py
```

`phenotype.R` 和 `plotting.R` 会从同一目录加载 `png_output.R`，后者调用 `png_output.py`。仅复制前两个 R 文件到一个临时目录不足以运行。请把 `ems120` 安装到完整 `scripts` 仓库的顶层。

## 3. 数据、模型、缓存

原始 EMS 输入目录需要包含：

```text
ems120/
  clean/
    2013.xlsx
    ...
    2024.xlsx
  szu_xia/
    2013.geo.xlsx
    ...
    2024.geo.xlsx
```

当前年度列表是 2013–2024，趋势图使用 2017–2024；权威定义及输入列名仍在 `f/ems120.R`。此次迁移不改变研究年份、疾病映射、评分规则或五折交叉验证定义。缺少某个年度的 `szu_xia/YYYY.geo.xlsx` 时，原逻辑会警告，后续需要地理关联变量的图可能停止；应核对输入是否齐全。

预训练基础模型目录应是已经下载完整的本地 Hugging Face 模型目录，包含模型配置、权重和 tokenizer 文件。默认使用 `--dir-big` 对应根目录下的 `AI_model/hfl`；可以用 `--dir-base-model` 明确指定。节点无外网时，不应把不存在的本地目录写成期望运行时自动下载的仓库名。

训练样本 `dat/2019.train_dx.xlsx` 仍位于分析目录，由原 R 数据整理流程创建/复用。默认模型目录是分析目录下的 `ml`，其下保留原来的 `hfl`、`hfl_geo`、`phone_sco3` 结构。年度 RDS、交叉验证输出及成图仍位于分析目录。

HPC 默认分析目录使用持久的 `/work` 路径，而临时目录由共享 `--dir-tmp` 配置控制。启明原共享配置注明带日期的 scratch 最长保留 8 天，所以模型和长期分析缓存不应默认放进这类目录。若显式把分析目录改到 scratch，应自行在清理前转存重要输出。

迁移旧 local 缓存时，应先复制完整分析目录并核对源数据/模型版本，再设置 `--dir-analysis` 指向副本。不要让本地和 HPC 两个正在运行的流程同时写同一组分析缓存或模型目录。

原分析仍按已有缓存规则复用 RDS 和模型。本次路径迁移不代表旧缓存自动适用于另一份数据；输入数据、标签、基础模型或分析定义改变时，使用新的分析目录并按需要重新训练/生成缓存。

## 4. GPU 的实际使用边界

`ml_dx` 和 `ml_geo` 使用 PyTorch Transformer，训练和推理可使用 GPU。它们仍在同一任务中顺序执行，地理模型在疾病模型阶段之前释放，避免同时占用显存。

电话模块包含规则评分、特征生成和多个候选回归器。旧版并没有给这些回归器传入 CUDA 参数；新版为支持的 XGBoost 分支启用 CUDA，并限制 CPU 候选的线程数。Elastic Net、scikit-learn Random Forest 和保留的 LightGBM 比较分支仍使用 CPU；最终模型按原交叉验证准则选择，不能因为有 GPU 就强制改变模型赢家。有效的既有 CPU 模型缓存也不会仅因迁移就被偷偷重训。

默认电话 GPU 策略为 `EMS120_PHONE_GPU=required`。需要训练 XGBoost 时，代码检查版本、构建与实际 CUDA 能力；仅仅 `torch.cuda.is_available()` 为真不能证明 XGBoost 的二进制也支持 CUDA。GPU 不可用时应清楚失败，不能把错误吞掉并把所有电话分数写成 NaN。

`--phone-gpu auto` 会在 GPU 不可用时显式报告并允许该候选使用 CPU；`--phone-gpu off` 明确关闭该候选的 GPU 路径。完整流水线的 PyTorch GPU 要求仍独立存在。只有显式 `--submit local --gpu FALSE` 才把整个本地运行切为 CPU，并把电话 GPU 策略同步设为 `off`。这用于无 GPU 的本地调试，不会成为 HPC 提交失败后的自动回退路径。

一张被调度器分配的 GPU 足够运行当前单进程设计。多申请 GPU 并不会自动把 Trainer、交叉验证或多个年度变成分布式计算。实际脚本应保留调度器的 `CUDA_VISIBLE_DEVICES`，使用分配后可见的逻辑设备，不能直接选择节点上的物理卡号。

## 5. 地图和无显示器绘图

当前 Fig2 默认保留普通二维坐标图，`draw_satellite_map` 和 `draw_3D_map` 均为 `FALSE`，此次迁移保留这些设置。只有启用可选的卫星图截图或 3D HTML 后，才需要相应浏览器/keplergl/网络依赖。`geo.map.use` 的值为 `keplergl_3d`，但它本身并不会覆盖上述关闭开关。

启用浏览器截图时，HPC 需要可执行的 Chromium / Google Chrome / Edge；可通过 `CHROME_BIN` 指定。程序使用无界面模式，不需要桌面会话。浏览器截图的 `--disable-gpu` 是浏览器自身的设置，不会关闭 PyTorch 的 CUDA。

ArcGIS 卫星瓦片、地图样式和部分 HTML 的 JavaScript 来自外部服务。无外网的计算节点可能无法首次下载这些资源；应按原代码准备已有瓦片缓存/允许访问的环境，并查看地图阶段日志。不要把缺失瓦片或未渲染的网页当作已验证的科研成图。

若 Fig2 Python 另有单独环境，可显式设置 `FIG2D_PYTHON`；否则使用所选分析环境的 Python。图中阴阳图素材从 `--dir-files` 指定目录或原工作目录候选位置读取，缺少素材会按原流程给出警告。

## 6. 验证边界

打包前可完成 CLI、路径、调度指令、退出码、并发锁、语法和轻量模型契约检查。真实启明 GPU、LSF 队列权限、R/reticulate/Conda ABI 以及未提供的 EMS 原始数据，需要在目标节点验证。先执行 GPU 检查任务，然后运行小范围分析并核对缓存/输出，再启动完整流程。
