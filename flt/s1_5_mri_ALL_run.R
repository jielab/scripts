# MRI 总入口：分析 -> R 热图 -> Python 解剖图/森林图；任一步失败即停止。
dir0 <- "D:"
project <- file.path(dir0, "projects/flt")
input <- file.path(project, "analysis/00data")

run_analysis <- TRUE    # 已有结果、只重画图时改为 FALSE
run_heatmap <- TRUE
run_python <- TRUE
skip_anatomy <- FALSE   # TRUE：Python 只画森林图，不画解剖图
python <- Sys.getenv("FLT_PYTHON", unset = "C:/Users/93677/.cache/codex-runtimes/codex-primary-runtime/dependencies/python/python.exe")
assets_project <- file.path(dir0, "projects/dose") # Python 使用的图谱和绘图工具

local({
  scripts <- file.path(project, "scripts")
  outdir <- file.path(project, "analysis/01lm/mri")
  result_file <- file.path(outdir, "mri_lm_results.csv")
  # 临时统一子脚本的项目路径，结束后恢复原设置。
  previous <- Sys.getenv("FLT_PROJECT", unset = NA_character_)
  on.exit(if (is.na(previous)) Sys.unsetenv("FLT_PROJECT") else Sys.setenv(FLT_PROJECT = previous))
  Sys.setenv(FLT_PROJECT = project)

  if (run_python && !file.exists(python)) stop("Python 不存在，请修改 python 路径：", python)
  if (run_python && !dir.exists(file.path(assets_project, ".tools"))) stop("缺少 MRI 绘图工具：", assets_project)
  run_r <- function(file) source(file.path(scripts, file), local = new.env(parent = globalenv()), encoding = "UTF-8")
  if (run_analysis) {
    message("1/3：MRI 回归分析")
    run_r("s1_5_mri_lm.R")
  }
  if ((run_heatmap || run_python) && !file.exists(result_file)) stop("缺少分析结果，请先设置 run_analysis <- TRUE：", result_file)
  if (run_heatmap) {
    message("2/3：MRI 热图（PDF）")
    run_r("s1_5_mri_lm.plot.R")
  }
  if (run_python) {
    message("3/3：MRI 解剖图与森林图")
    args <- c(shQuote(file.path(scripts, "s1_5_mri_lm.plot.py")),
              "--assets-project", shQuote(assets_project),
              "--results", shQuote(result_file),
              "--out", shQuote(file.path(outdir, "figures/anatomy")),
              if (skip_anatomy) "--skip-anatomy")
    status <- system2(python, args = args)
    if (status != 0) stop("Python 绘图失败（退出码 ", status, "），请检查控制台及 figures/anatomy 下的 plot_status.txt。")
  }
  message("MRI 流程完成：", outdir)
})
