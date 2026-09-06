# MAHA 分析

在 WSL / Ubuntu 终端运行全部队列并整合投稿图表：

```bash
cd /mnt/d/scripts/maha
bash maha.sh all --replace-steps all --jobs 4
```

也可以分别运行 `bash maha.sh ukb`、`bash maha.sh nhanes`、`bash maha.sh chns`，最后执行 `bash maha.sh final`。更多选项见 `bash maha.sh --help`。

输出目录为 `/mnt/d/analysis/maha`，程序自动创建。队列日志位于对应子目录的 `launcher.log`；运行状态见 `run_status.txt`。正式分析不要使用缩减诊断次数的 `--null-max` 测试设置。

`maha.sh` 是分析入口。`f/MAHA_ukb.R`、`f/MAHA_nhanes.R`、`f/MAHA_chns.R` 分别分析三个队列，`f/MAHA_publication.R` 整合投稿图表。`f/comm.f.R` 为共用流程函数；`f/ukb_followup.R` 和 `f/nhanes_meat_sensitivity.R` 为正式分析依赖，不能删除。

`f/ukb_diet.R` 是上游 UKB 评分构建代码。当前分析使用既有评分数据，无需为本次随访修正单独重新运行它。

UKB 从完整 WebQ 评估窗口结束后开始随访；PheWAS 保留原始 P 值，仅绘图时封顶。完整组成及红肉代理敏感性结果保存在对应队列目录的 `Reviewer.*` 文件中。这些是正式分析导出文件，与已清理的临时 `reviewer_revision` 目录不同。

分析命令不会自动修改 NatureH 下的投稿 Word。应在正式结果完整生成并核对后同步更新文稿。
