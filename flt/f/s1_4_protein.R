source(file.path(Sys.getenv("FLT_CODE"), "runtime.R"))
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 三组食物偏好/摄入与蛋白组
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# limma: 极喜欢(恰好9分) vs 极不喜欢(恰好1分)，摄入上/下四分位；
# 两项各自 BH<0.05 且同向为探索性候选，不等于联合FDR或中介证据。
# linear: 各食物交集蛋白 ~ 连续偏好 + 协变量，普通OLS标准误。
# Cox: 第2步线性回归FDR<0.05蛋白 × 对应结局；RCS: 第3步Cox FDR<0.05蛋白。

suppressPackageStartupMessages({
  library(tidyverse)
  library(data.table)
  library(survival)
  library(readxl)
  library(ggplot2)
  library(limma)
})

run_modules <- c("limma", "linear", "cox", "rcs")
resume_models <- FALSE

project <- Sys.getenv("FLT_OUT")
input <- file.path(project, "00data")
indir <- Sys.getenv("FLT_UKB")
outdir <- file.path(project, "01protein")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
config <- readRDS(file.path(input, "analysis_config.rds"))
base_file <- file.path(input, config$primary_data_file)
diet_file <- file.path(indir, "Rdata/diet.fudan.rds")
prot_file <- file.path(indir, "Rdata/prot.rds")

liking_selected <- c("Butter on bread" = "p20628", "Total salad vegetable" = "F-Salad_vegetables", "Chicken" = "p20638")
covs.all <- c("age", "sex", "bmi", "edu", "inc", "tdi", "smoke_status", "alcohol_status", "pa_met_tot")
food_pairs <- data.frame(
  pair_id=c("G717", "G824", "G732"),
  disease_id=c("sleep_disorder_r13", "sleep_disorder_r13", "mdd_2025"),
  outcome=c("sleep_disorder", "sleep_disorder", "depression"),
  food=names(liking_selected), liking=unname(liking_selected),
  intake=c("w.p26063", "w.p26123", "w.p26121"),
  intake_field_id=c(26063L,26123L,26121L),
  proxy=c("Butter", "Salad/raw vegetables", "Poultry"), intake_unit="g/day")

read_data <- function(path) {
  x <- as.data.frame(readRDS(path))
  if (!"eid" %in% names(x) || anyDuplicated(names(x))) stop("eid缺失或列名重复：", path)
  x$eid <- trimws(as.character(x$eid))
  if (anyNA(x$eid) || any(!nzchar(x$eid)) || anyDuplicated(x$eid)) stop("eid缺失或重复：", path)
  x
}
file_stamp <- function(path) {
  info <- file.info(path)
  list(path = normalizePath(path), size = info$size, mtime = as.numeric(info$mtime))
}

# ==================== LIMMA ====================
if ("limma" %in% run_modules) local({
outdir <- file.path(project, "01protein/limma")

intake_probs <- c(.25, .75)       # 低/高摄入切点；相同摄入值不拆组
positive_high_prob <- .75        # 零值过多时，高组取正摄入者的第75百分位及以上
covariates <- covs.all
protein_subset <- NULL          # NULL=prot中除eid外全部蛋白；或填写指定蛋白名
min_per_group <- 3L
max_missing <- .20              # 每个蛋白在每组的最大缺失比例，自定QC规则
alpha <- .05
pairs <- food_pairs

base <- read_data(base_file)
diet <- read_data(diet_file)
prot <- read_data(prot_file)
required <- unique(c(pairs$liking, covariates))
proteins <- if (is.null(protein_subset)) setdiff(names(prot), "eid") else protein_subset

dat <- base[base$eid %in% prot$eid, c("eid", required), drop = FALSE]
dat[pairs$intake] <- diet[match(dat$eid, diet$eid), pairs$intake, drop = FALSE]
protein_matrix <- t(as.matrix(prot[match(dat$eid, prot$eid), proteins, drop = FALSE]))
protein_matrix[!is.finite(protein_matrix)] <- NA_real_
observed <- colSums(!is.na(protein_matrix)) > 0
dat <- dat[observed, , drop = FALSE]
protein_matrix <- protein_matrix[, observed, drop = FALSE]
if (!nrow(dat)) stop("没有匹配且有蛋白测量的参与者")
rm(base, diet, prot)
invisible(gc())

cov_complete <- rep(TRUE, nrow(dat))
if (length(covariates)) cov_complete <- complete.cases(dat[covariates])
for (v in covariates) {if (is.numeric(dat[[v]])) cov_complete <- cov_complete & is.finite(dat[[v]])}
for (v in c(pairs$liking, pairs$intake)) {if (!is.numeric(dat[[v]])) stop("评分/摄入列必须为数值：", v)}
for (v in pairs$liking) {if (any(is.finite(dat[[v]]) & (dat[[v]] < 1 | dat[[v]] > 9))) stop("评分超出1-9范围：", v)}
for (v in pairs$intake) {if (any(is.finite(dat[[v]]) & dat[[v]] < 0)) stop("摄入含负值：", v)}

# 3. 同一个函数完成两组比较 -------------------------------------------
empty_result <- data.frame(food = character(), analysis = character(), protein = character(), grouping_rule = character(), effect = numeric(), P.Value = numeric(), adj.P.Val = numeric(), n_low = integer(), n_high = integer(), direction = character())
run_limma <- function(food, kind, field) {
  values <- dat[[field]]
  available <- is.finite(values) & cov_complete
  cuts <- if (kind == "liking") c(1, 9) else if (any(available))
    quantile(values[available], intake_probs, names = FALSE) else c(NA_real_, NA_real_)
  original_cuts <- cuts
  rule <- if (kind == "liking") "score_9_vs_1" else "intake_quantiles"
  # 零摄入不是缺失。仅在两分位点均为0时启用替代规则，不强拆相同值。
  positive <- values[available & values > 0]
  if (kind == "intake" && all(is.finite(cuts)) && all(cuts == 0) && length(positive)) {
    cuts <- c(0, quantile(positive, positive_high_prob, names = FALSE))
    rule <- "zero_vs_high_positive_intake"
  }
  group <- rep(NA_character_, nrow(dat))
  valid_cuts <- all(is.finite(cuts)) && cuts[1] < cuts[2]
  if (valid_cuts) {
    if (kind == "liking") {
      group[available & values == 1] <- "low"
      group[available & values == 9] <- "high"
    } else {
      group[available & values <= cuts[1]] <- "low"
      group[available & values >= cuts[2]] <- "high"
    }
  }
  keep <- !is.na(group)
  counts <- table(factor(group[keep], levels = c("low", "high")))
  pair <- pairs[match(food, pairs$food), ]
  summary <- data.frame(food = food, analysis = kind, input_field = field,
                        grouping_rule = rule,
                        field_id = if (kind == "intake") pair$intake_field_id else NA_integer_,
                        input_name = if (kind == "intake") pair$proxy else paste(food, "liking"),
                        unit = if (kind == "intake") pair$intake_unit else "score (1-9)",
                        n_available = sum(available),
                        n_zero = if (kind == "intake") sum(values[available] == 0) else NA_integer_,
                        n_positive = if (kind == "intake") length(positive) else NA_integer_,
                        quantile_low = if (kind == "intake") original_cuts[1] else NA_real_,
                        quantile_high = if (kind == "intake") original_cuts[2] else NA_real_,
                        positive_high_prob = if (rule == "zero_vs_high_positive_intake")
                          positive_high_prob else NA_real_,
                        low_cut = cuts[1], high_cut = cuts[2],
                        n_low = as.integer(counts[1]), n_high = as.integer(counts[2]),
                        n_tested = 0L, status = "ok")
  result <- empty_result
  if (!valid_cuts) {
    summary$status <- if (kind == "intake" && any(available) && !length(positive))
      "no_positive_intake" else "missing_or_tied_cutoffs"
  } else if (min(counts) < min_per_group) {
    summary$status <- "insufficient_group_size"
  } else {
    model_data <- dat[keep, covariates, drop = FALSE]
    model_data$group <- factor(group[keep], levels = c("low", "high"))
    design <- model.matrix(reformulate(c("group", covariates)), model_data)
    if (qr(design)$rank < ncol(design) || nrow(design) <= ncol(design)) stop("协变量共线或模型自由度不足：", food, " / ", kind)

    y <- protein_matrix[, keep, drop = FALSE]
    low <- model_data$group == "low"
    n_low <- rowSums(!is.na(y[, low, drop = FALSE]))
    n_high <- rowSums(!is.na(y[, !low, drop = FALSE]))
    variance <- apply(y, 1, var, na.rm = TRUE)
    pass <- n_low >= min_per_group & n_high >= min_per_group & n_low / counts[1] >= 1 - max_missing & n_high / counts[2] >= 1 - max_missing & is.finite(variance) & variance > 0
    if (sum(pass) >= 2L) {
      fit <- limma::lmFit(y[pass, , drop = FALSE], design)
      estimable <- fit$df.residual > 0 & is.finite(fit$sigma) & rowSums(!is.finite(fit$coefficients)) == 0 & rowSums(!is.finite(fit$stdev.unscaled)) == 0
      summary$n_tested <- sum(estimable)
      if (sum(estimable) >= 2L) {
        fit <- limma::eBayes(fit[estimable, ])
        tt <- limma::topTable(fit, coef = "grouphigh", number = Inf, sort.by = "none", adjust.method = "BH")
        # effect=高组减低组，在输入单位下解释；不自动称为倍数变化。
        result <- data.frame(food = food, analysis = kind, protein = rownames(tt), grouping_rule = rule,
                             effect = tt$logFC, P.Value = tt$P.Value, adj.P.Val = tt$adj.P.Val, n_low = unname(n_low[rownames(tt)]), n_high = unname(n_high[rownames(tt)]),
                             direction = ifelse(tt$adj.P.Val < alpha & tt$logFC != 0, ifelse(tt$logFC > 0, "up", "down"), "not_significant"))
      }
    }
    if (!nrow(result)) summary$status <- "fewer_than_two_testable_proteins"
  }
  list(result = result, summary = summary)
}

fits <- list()
for (i in seq_len(nrow(pairs))) {
  for (kind in c("liking", "intake")) {
    key <- paste(pairs$food[i], kind, sep = "_")
    fits[[key]] <- run_limma(pairs$food[i], kind, pairs[[kind]][i])
  }
}
results <- bind_rows(lapply(fits, `[[`, "result"))
groups <- bind_rows(lapply(fits, `[[`, "summary"))

# 4. 按食物和蛋白取交集：两项BH-FDR<0.05，且效应同向 -------------------
# BH分别在每个食物、每次比较内校正；交集筛选不等于联合FDR检验。
cols <- c("food", "protein", "grouping_rule", "effect", "P.Value", "adj.P.Val")
paired <- merge(results[results$analysis == "liking", cols], results[results$analysis == "intake", cols], by = c("food", "protein"), suffixes = c("_liking", "_intake"))
shared <- paired %>% filter(is.finite(adj.P.Val_liking), is.finite(adj.P.Val_intake), adj.P.Val_liking < alpha, adj.P.Val_intake < alpha, effect_liking * effect_intake > 0)

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
csv <- function(x, name) write.csv(x, file.path(outdir, name), row.names = FALSE, na = "")
csv(pairs, "food_pair_mapping.csv")
csv(groups, "group_counts_and_status.csv")
csv(results, "all_limma_results.csv")
csv(paired, "paired_protein_results.csv")
csv(shared, "shared_same_direction_proteins.csv")
saveRDS(unique(shared$protein), file.path(outdir, "candidate_proteins_for_disease.rds"))
capture.output(sessionInfo(), file = file.path(outdir, "sessionInfo.txt"))
if (!nrow(results)) stop("没有成功模型，请检查group_counts_and_status.csv")
if (any(groups$status != "ok")) warning("部分比较跳过，原因见group_counts_and_status.csv")
message("完成：", outdir, "；同向显著食物-蛋白配对数：", nrow(shared), "；去重蛋白数：", length(unique(shared$protein)))

# 5. 参考图式双层效应散点图（不是传统火山图） -----------------------
# 可单独运行本段：直接读取limma输出，无需重新拟合模型。
# y轴使用输入尺度差异；只有确认输入为log2尺度后才可标为log2FC。
draw_effect_strip <- function(limma_dir, protein_scale = "input") {
  # 确认prot.rds为未做z标准化/秩逆正态变换的log2 NPX后，才设置为"log2_NPX"。
  stopifnot(protein_scale %in% c("input", "log2_NPX"))
  axis_label <- if(protein_scale == "log2_NPX") expression("Adjusted " * log[2] * "FC") else
    "Adjusted protein difference (high - low; input scale)"
  d <- read.csv(file.path(limma_dir, 'all_limma_results.csv'))
  m <- read.csv(file.path(limma_dir, 'food_pair_mapping.csv'))
  shared <- read.csv(file.path(limma_dir, 'shared_same_direction_proteins.csv'))
  d <- d[is.finite(d$effect) & is.finite(d$adj.P.Val), ]
  if (!nrow(d)) { warning('无有效limma结果，跳过效应散点图'); return(invisible(NULL)) }
  d$food <- factor(d$food, levels = m$food)
  d$panel <- factor(d$analysis, levels = c('liking', 'intake'), labels = c('Liking', 'Intake'))
  d$direction <- ifelse(d$adj.P.Val < .05,
    ifelse(d$effect > 0, 'Up', 'Down'), 'Not significant')
  # 固定抖动种子，横坐标仅用于分散点，没有统计含义。
  set.seed(42)
  d$x <- as.integer(d$food) + runif(nrow(d), -.40, .40)
  key <- function(food, protein) paste(food, protein, sep = '\r')
  d$is_shared <- key(d$food, d$protein) %in% key(shared$food, shared$protein)
  food_labels <- m$food
  food_labels[food_labels == "Total salad vegetable"] <- "Salad vegetables"
  # 每个面板使用自身的效应范围。
  panel_names <- levels(d$panel)
  panel_extent <- setNames(vapply(panel_names, function(panel) {
    max(c(.01, abs(d$effect[d$panel == panel])))
  }, numeric(1)), panel_names)
  d$panel_extent <- unname(panel_extent[as.character(d$panel)])
  bg <- expand.grid(food = m$food, panel = panel_names, stringsAsFactors = FALSE)
  bg$x <- match(bg$food, m$food)
  bg$half <- unname(panel_extent[bg$panel]) * .085
  bg$ymin <- vapply(seq_len(nrow(bg)), function(i) {
    z <- d$effect[d$food == bg$food[i] & d$panel == bg$panel[i]]
    min(c(z, -bg$half[i]))
  }, numeric(1))
  bg$ymax <- vapply(seq_len(nrow(bg)), function(i) {
    z <- d$effect[d$food == bg$food[i] & d$panel == bg$panel[i]]
    max(c(z, bg$half[i]))
  }, numeric(1))
  food_colors <- setNames(c('#3366A5', '#16854B', '#C66A27'), m$food)
  labs <- d[d$is_shared, ]
  # 固定双列分层标注，避免窄食物区域内ggrepel标签互相挤压。
  labs$label_x <- labs$label_y <- NA_real_
  for (food_i in seq_len(nrow(m))) for (panel_i in levels(d$panel)) for (sgn in c(-1,1)) {
    ii <- which(labs$food==m$food[food_i] & labs$panel==panel_i & sign(labs$effect)==sgn)
    if(!length(ii)) next
    ii <- ii[order(labs$x[ii])]
    ex <- max(panel_extent[[panel_i]], .01)
    start <- max(abs(labs$effect[ii])) + .08*ex
    n_left <- ceiling(length(ii)/2)
    labs$label_x[ii] <- food_i + if(length(ii)==1) 0 else c(rep(-.23,n_left),rep(.23,length(ii)-n_left))
    tiers <- c(seq_len(n_left), if(length(ii)>n_left) seq_len(length(ii)-n_left) else integer())
    labs$label_y[ii] <- sgn*(start + (tiers-1 + ifelse(labs$label_x[ii]>food_i,.50,0))*.115*ex)
  }
  panels <- list()
  for (panel_i in panel_names) {
    z <- d[d$panel == panel_i, ]
    lab <- labs[labs$panel == panel_i, ]
    bk <- bg[bg$panel == panel_i, ]
    ex <- panel_extent[[panel_i]]
    # 色带只占用显示空间。正负点分开画，真实零值保留在色带中心。
    gap <- .13 * ex
    display_y <- function(value) value + sign(value) * gap
    z$y_plot <- display_y(z$effect)
    lab$y_plot <- display_y(lab$effect)
    lab$label_plot <- display_y(lab$label_y)
    raw_range <- range(c(z$effect, lab$label_y, bk$ymin, bk$ymax, 0))
    ticks <- pretty(raw_range, n = 6)
    ticks <- sort(unique(c(0, ticks[ticks >= raw_range[1] & ticks <= raw_range[2]])))
    tick_labels <- ifelse(ticks == 0, '0', format(ticks, trim = TRUE))
    yrange <- range(c(display_y(raw_range), -gap, gap)) + c(-1, 1) * ex * .08
    bk$label <- food_labels[bk$x]
    bk$label[bk$label=='Butter on bread'] <- 'Butter on\nbread'
    bk$label[bk$label=='Salad vegetables'] <- 'Salad\nvegetables'
    bk$ymin <- display_y(bk$ymin)
    bk$ymax <- display_y(bk$ymax)
    panels[[panel_i]] <- ggplot2::ggplot() +
      ggplot2::geom_rect(data=bk,ggplot2::aes(xmin=x-.45,xmax=x+.45,ymin=ymin,ymax=ymax),fill='#F0F1F3') +
      ggplot2::geom_rect(data=bk,ggplot2::aes(xmin=x-.5,xmax=x+.5,fill=food),
        ymin=-gap,ymax=gap) +
      ggplot2::geom_text(data=bk,ggplot2::aes(x,0,label=label),color='white',fontface='bold',size=4,lineheight=.95) +
      ggplot2::geom_point(data=z[z$direction=='Not significant', ],ggplot2::aes(x,y_plot,color=direction),size=.85,alpha=.65) +
      ggplot2::geom_point(data=z[z$direction!='Not significant', ],ggplot2::aes(x,y_plot,color=direction),size=1.25,alpha=.95) +
      ggplot2::geom_segment(data=lab,ggplot2::aes(x=x,y=y_plot,xend=label_x,yend=label_plot),color='#7A8591',linewidth=.32) +
      ggplot2::geom_label(data=lab,ggplot2::aes(label_x,label_plot,label=protein),fill='white',color='#26313E',
        fontface='bold',size=3.3,label.padding=grid::unit(.12,'lines'),linewidth=.25) +
      ggplot2::scale_x_continuous(limits=c(.35,nrow(m)+.65),breaks=NULL,
        expand=ggplot2::expansion(mult=0)) +
      ggplot2::scale_y_continuous(breaks=display_y(ticks),labels=tick_labels,
        limits=yrange,expand=ggplot2::expansion(mult=0)) +
      ggplot2::scale_fill_manual(values=food_colors,guide='none') +
      ggplot2::scale_color_manual(values=c(Up='#D73027',Down='#3366A5','Not significant'='#929292'),
        breaks=c('Up','Down','Not significant'),labels=c('Up: FDR < 0.05','Down: FDR < 0.05','FDR >= 0.05')) +
      ggplot2::guides(color=ggplot2::guide_legend(nrow=2,byrow=TRUE)) +
      ggplot2::labs(title=panel_i,x=NULL,y=axis_label,color=NULL) +
      ggplot2::theme_classic(base_size=13) +
      ggplot2::theme(axis.line.x=ggplot2::element_blank(),
        axis.line.y=ggplot2::element_line(color="black",linewidth=.42),axis.ticks.x=ggplot2::element_blank(),
        plot.title=ggplot2::element_text(size=20,face='bold',hjust=.5),
        legend.position=if(panel_i=='Intake') 'bottom' else 'none',
        legend.text=ggplot2::element_text(size=11),
        axis.text.y=ggplot2::element_text(size=12),axis.title.y=ggplot2::element_text(size=14),
        plot.margin=ggplot2::margin(10,18,10,14))
  }
  # 竖版尺寸（英寸）：统一控制PNG/PDF，方便后续调整。
  figure_width <- 7
  figure_height <- 12
  dest <- file.path(limma_dir, 'effect_strip')
  dir.create(dest, recursive=TRUE,showWarnings=FALSE)
  render <- function() {
    grid::grid.newpage()
    grid::pushViewport(grid::viewport(layout=grid::grid.layout(3,1,heights=grid::unit(c(1,1.10,.20),'null'))))
    for(i in seq_along(panels)) print(panels[[i]],vp=grid::viewport(layout.pos.row=i,layout.pos.col=1))
    grid::grid.text('Liking: score 9 vs 1. Intake: high vs low.\nLabels: concordant shared proteins (both FDR < .05).\nThe colour band marks a display gap at zero.\nPanels use different y-axis ranges.',
      x=.03,hjust=0,gp=grid::gpar(fontsize=10),vp=grid::viewport(layout.pos.row=3,layout.pos.col=1))
    grid::popViewport()
  }
  grDevices::png(file.path(dest,'liking_intake_effect_strip.png'),width=figure_width,height=figure_height,units='in',res=300)
  tryCatch(render(),finally=grDevices::dev.off())
  grDevices::pdf(file.path(dest,'liking_intake_effect_strip.pdf'),width=figure_width,height=figure_height)
  tryCatch(render(),finally=grDevices::dev.off())
  write.csv(labs, file.path(dest, 'annotated_shared_proteins.csv'), row.names = FALSE)
  message('双层效应散点图已保存：', dest)
  invisible(panels)
}
draw_effect_strip(file.path(project, '01protein/limma'), protein_scale = "log2_NPX")

})

# 后续严格按食物取交集，不将三组候选混用。
read_candidates <- function() {
  path <- file.path(project, "01protein/limma/shared_same_direction_proteins.csv")
  if (!file.exists(path)) stop("请先运行limma模块：", path)
  x <- read.csv(path, stringsAsFactors = FALSE)
  required <- c("food", "protein", "adj.P.Val_liking", "adj.P.Val_intake", "effect_liking", "effect_intake")
  if (!all(required %in% names(x))) stop("交集文件缺列，请重跑limma")
  x <- x[with(x, is.finite(adj.P.Val_liking) & is.finite(adj.P.Val_intake) &
    adj.P.Val_liking < .05 & adj.P.Val_intake < .05 & effect_liking * effect_intake > 0), ]
  merge(food_pairs, unique(x[c("food", "protein")]), by = "food", sort = FALSE)
}

# 按上一阶段显著结果逐级筛选；先限定当前limma候选，防止跨食物混用。
read_stage_targets <- function(pair, stage = c("linear", "cox")) {
  stage <- match.arg(stage)
  targets <- read_candidates()
  targets <- targets[targets$pair_id == pair$pair_id, , drop = FALSE]
  if (!nrow(targets)) return(targets)
  path <- file.path(project, "01protein/linear/FLT_protein_linear_all.csv")
  if (!file.exists(path)) stop("请先运行linear模块：", path)
  x <- read.csv(path, stringsAsFactors = FALSE)
  if (!all(c("food", "protein", "status", "FDR_food") %in% names(x)))
    stop("线性回归结果缺列，请重跑linear")
  x <- x[x$food == pair$food & x$status == "ok" &
    is.finite(x$FDR_food) & x$FDR_food < .05, , drop = FALSE]
  targets <- targets[targets$protein %in% x$protein, , drop = FALSE]
  if (stage == "linear" || !nrow(targets)) return(targets)
  path <- file.path(project, "01protein/cox", pair$pair_id, "protein_cox_all.csv")
  if (!file.exists(path)) stop("请先运行cox模块：", path)
  x <- read.csv(path, stringsAsFactors = FALSE)
  if (!all(c("pair_id", "outcome", "protein", "status", "FDR_outcome") %in% names(x)))
    stop("Cox结果缺列，请重跑cox")
  x <- x[x$pair_id == pair$pair_id & x$outcome == pair$outcome & x$status == "ok" &
    is.finite(x$FDR_outcome) & x$FDR_outcome < .05, , drop = FALSE]
  targets[targets$protein %in% x$protein, , drop = FALSE]
}

# 无候选时归档旧表和图，防止旧的显著结果被误认为本次结果。
record_empty_stage <- function(outdir, reason) {
  old <- list.files(outdir, pattern = "\\.(csv|png|pdf|txt)$", full.names = TRUE)
  old <- old[!basename(old) %in% c("selected_targets.csv", "stage_status.csv")]
  if (length(old)) {
    archive <- tempfile("previous_results_", tmpdir = outdir)
    dir.create(archive)
    if (!all(file.rename(old, file.path(archive, basename(old)))))
      stop("归档旧结果失败：", outdir)
  }
  write.csv(data.frame(status = "skipped", reason = reason, n_selected = 0L),
    file.path(outdir, "stage_status.csv"), row.names = FALSE)
  message(reason)
}

# ==================== LINEAR ====================
if ("linear" %in% run_modules) local({
outdir <- file.path(project, "01protein/linear")
linear_dir <- outdir
linear_covariates <- covs.all
alpha <- .05
overview_top_n <- 10L       # 每食物概览最多10个；按本次回归FDR排序，仅限制展示
rows_per_page <- 40L        # 完整PDF分页展示所有配对，不按本次P值删减

stopifnot(overview_top_n >= 1, rows_per_page >= 1)
dir.create(linear_dir, recursive = TRUE, showWarnings = FALSE)
write_linear_csv <- function(x, name) {write.csv(x, file.path(linear_dir, name), row.names = FALSE, na = '')}

# 2. 读取原始评分、协变量和全部蛋白 -----------------------------------
mapping <- food_pairs

base <- read_data(base_file)
prot <- read_data(prot_file)
targets <- read_candidates()
if (!nrow(targets)) {
  write_linear_csv(targets, 'selected_targets.csv')
  record_empty_stage(outdir, '无交集蛋白，跳过线性回归'); return(invisible(NULL))
}
proteins <- unique(targets$protein)
if (!all(proteins %in% names(prot))) stop('交集蛋白不在prot中，请重跑limma')
write_linear_csv(targets, 'selected_targets.csv')
targets <- targets[order(match(targets$food, mapping$food), targets$protein), ]
required <- unique(c(mapping$liking, linear_covariates))
if (!all(required %in% names(base))) stop("base缺少评分或协变量")
if (!length(proteins)) stop("无蛋白列")
numeric_vars <- unique(c(mapping$liking, setdiff(linear_covariates, c("sex","edu","inc","smoke_status","alcohol_status"))))

# 依eid一对一匹配；沿用base人群，不额外按疾病或摄入情况筛选。
dat <- base[base$eid %in% prot$eid, c('eid', required), drop = FALSE]
protein_data <- prot[match(dat$eid, prot$eid), proteins, drop = FALSE]
for(v in c("sex", "edu", "inc", "smoke_status", "alcohol_status")) dat[[v]] <- factor(dat[[v]])
for(v in names(dat)) if(is.numeric(dat[[v]])) dat[[v]][!is.finite(dat[[v]])] <- NA_real_
for(v in proteins) {
  if(!is.numeric(protein_data[[v]])) stop("非数值蛋白：", v)
  protein_data[[v]][!is.finite(protein_data[[v]])] <- NA_real_
}
rm(base, prot)
invisible(gc())

# 3. 每个“食物—蛋白”配对分别做普通线性回归 ---------------------------
# β：偏好每增加1分，蛋白平均改变多少输入单位（若为NPX，即NPX单位）。
# beta_protein_SD：同一β除以该模型样本中的蛋白SD，便于跨蛋白比较。
# 每个模型只剔除自身评分、蛋白、上述9项协变量缺失者。
model_formula <- reformulate(c('liking', linear_covariates), response = 'protein')
run_linear <- function(i) {
  z <- targets[i, , drop = FALSE]
  field <- mapping$liking[match(z$food, mapping$food)]
  d <- dat[linear_covariates]
  d$liking <- dat[[field]]
  d$protein <- protein_data[[z$protein]]
  d <- droplevels(d[complete.cases(d), , drop = FALSE])
  z$liking_field <- field
  z$n <- nrow(d)
  z$n_liking_levels <- length(unique(d$liking))
  z$protein_SD <- if (nrow(d) > 1L) sd(d$protein) else NA_real_
  for (v in c('beta', 'SE', 'CI_low', 'CI_high', 'P', 'beta_protein_SD', 'CI_low_protein_SD', 'CI_high_protein_SD')) z[[v]] <- NA_real_
  z$status <- 'ok'
  z$reason <- ''
  tryCatch({
    if (nrow(d) < 3L) stop('完整观测不足')
    if (z$n_liking_levels < 2L) stop('偏好评分无变异')
    if (!is.finite(z$protein_SD) || z$protein_SD <= 0) stop('蛋白无变异')
    x <- model.matrix(model_formula, d)
    if (qr(x)$rank < ncol(x)) stop('设计矩阵不满秩，请检查协变量')
    if (nrow(x) <= ncol(x)) stop('模型残差自由度不足')
    fit <- lm(model_formula, data = d)
    a <- coef(summary(fit))['liking', ]
    if (!all(is.finite(a))) stop('回归估计非有限值')
    margin <- qt(.975, df.residual(fit)) * a['Std. Error']
    z$beta <- unname(a['Estimate'])
    z$SE <- unname(a['Std. Error'])
    z$CI_low <- z$beta - margin
    z$CI_high <- z$beta + margin
    z$P <- unname(a['Pr(>|t|)'])
    z$beta_protein_SD <- z$beta / z$protein_SD
    z$CI_low_protein_SD <- z$CI_low / z$protein_SD
    z$CI_high_protein_SD <- z$CI_high / z$protein_SD
    z
  }, error = function(e) {
    z$status <- 'failed'
    z$reason <- conditionMessage(e)
    z
  })
}
linear_results <- do.call(rbind, lapply(mapping$food, function(f) {
  message('Linear regression: ', f, ' x ', length(proteins), ' proteins')
  do.call(rbind, lapply(which(targets$food == f), run_linear))
}))
# BH按每种食物内的交集蛋白校正；另附三种食物全部配对一起校正。
# 失败模型留在计划检验数中，不用减少检验数来降低FDR。
linear_results$FDR_food <- rep(NA_real_, nrow(linear_results))
for (f in mapping$food) {
  ii <- which(linear_results$food == f & is.finite(linear_results$P))
  if (length(ii)) linear_results$FDR_food[ii] <- p.adjust(linear_results$P[ii], 'BH', n = sum(linear_results$food == f))
}
linear_results$FDR_all <- rep(NA_real_, nrow(linear_results))
ii <- which(is.finite(linear_results$P))
if (length(ii)) linear_results$FDR_all[ii] <- p.adjust(linear_results$P[ii], 'BH', n = nrow(linear_results))
linear_results$significant <- !is.na(linear_results$FDR_food) & linear_results$FDR_food < alpha
linear_results$stars <- with(linear_results, ifelse(is.na(FDR_food), '',
  ifelse(FDR_food < .001, '***', ifelse(FDR_food < .01, '**', ifelse(FDR_food < .05, '*', '')))))
write_linear_csv(linear_results, 'FLT_protein_linear_all.csv')
write_linear_csv(linear_results[linear_results$significant, ], 'FLT_protein_linear_significant.csv')
linear_counts <- do.call(rbind, lapply(mapping$food, function(f) {
  z <- linear_results[linear_results$food == f, ]
  data.frame(food = f, n_proteins = nrow(z), n_fitted = sum(z$status == 'ok'),
    n_failed = sum(z$status != 'ok'), n_FDR_significant = sum(z$significant),
    n_positive = sum(z$significant & z$beta > 0, na.rm = TRUE),
    n_negative = sum(z$significant & z$beta < 0, na.rm = TRUE))
}))
write_linear_csv(linear_counts, 'FLT_protein_linear_counts.csv')

# 4. Fig.4b风格：颜色表示β，星号表示FDR，旁边打印β数值 ----------------
# 主图使用未标准化β；原文图注未说明是否标准化，故不强行匹配其数值尺度。
finite_beta <- linear_results$beta[is.finite(linear_results$beta)]
plot_limit <- if (length(finite_beta)) max(abs(finite_beta)) else 1
if (plot_limit == 0) plot_limit <- 1
plot_names <- setNames(gsub('_', ' ', mapping$food), mapping$food)
plot_note <- paste0('Adjusted for all nine project covariates. Protein input units per liking point.\n',
  '* BH-FDR < .05; ** < .01; *** < .001 (within food). NA: failed model.\n',
  'Only food-specific same-direction limma intersection proteins analyzed.')
draw_linear_heatmap <- function(z, subtitle) {
  if (!nrow(z)) return(ggplot() + annotate('text', x = 0, y = 0,
    label = 'No proteins available') + theme_void() +
    theme(plot.background = element_rect(fill = 'white', color = NA)) +
    labs(title = 'Fig.4b', subtitle = subtitle))
  z$row_id <- factor(seq_len(nrow(z)), levels = rev(seq_len(nrow(z))))
  row_labels <- setNames(z$protein, as.character(seq_len(nrow(z))))
  z$food <- factor(z$food, levels = mapping$food)
  z$beta_label <- ifelse(is.finite(z$beta), sprintf('%.3f', z$beta), 'NA')
  ggplot(z, aes(y = row_id)) +
    geom_tile(aes(x = 1, fill = beta), width = 1, height = 1, color = 'grey60', linewidth = .25) +
    geom_tile(aes(x = 2), width = 1, height = 1, fill = 'white', color = 'grey60', linewidth = .25) +
    geom_text(aes(x = 1, label = stars), size = 3.2) +
    geom_text(aes(x = 2, label = beta_label), size = 3.1) +
    facet_grid(food ~ ., scales = 'free_y', space = 'free_y', labeller = labeller(food = plot_names)) +
    scale_fill_gradient2(low = '#425D8B', mid = 'white', high = '#C84D4D',
      midpoint = 0, limits = c(-plot_limit, plot_limit), na.value = 'grey90', name = 'Beta') +
    scale_x_continuous(breaks = c(1, 2), labels = c('Beta / FDR', 'Beta'),
      limits = c(.5, 2.5), expand = expansion(mult = 0), position = 'top') +
    scale_y_discrete(labels = row_labels, expand = expansion(add = .5)) +
    labs(x = NULL, y = NULL, title = 'FLT - circulating protein associations',
      subtitle = subtitle, caption = plot_note) +
    theme_minimal(base_size = 11) +
    theme(plot.background = element_rect(fill = 'white', color = NA),
      panel.grid = element_blank(), axis.ticks = element_blank(),
      axis.text.y = element_text(face = 'bold', color = 'black', size = 9),
      strip.text.y = element_text(angle = 0, face = 'bold'),
      legend.position = 'bottom', plot.caption = element_text(size = 7, hjust = 0),
      panel.spacing.y = grid::unit(.5, 'lines'))
}
# 分析完成并校正后，才选取概览；完整CSV/PDF仍包含所有配对和失败原因。
linear_results <- linear_results[order(match(linear_results$food, mapping$food), linear_results$FDR_food, linear_results$P, linear_results$protein, na.last = TRUE), ]
overview_rows <- unlist(lapply(mapping$food, function(f)
  head(which(linear_results$food == f), overview_top_n)), use.names = FALSE)
overview <- linear_results[overview_rows, , drop = FALSE]
write_linear_csv(overview, 'Fig4b_overview_rows.csv')
subtitle <- paste0('Top ', overview_top_n, ' per food by regression BH-FDR; ',
  nrow(targets), ' selected food-protein pairs analyzed')
p <- draw_linear_heatmap(overview, subtitle)
h <- max(5, 3 + nrow(overview) * .23)
ggsave(file.path(linear_dir, 'Fig4b_overview.png'), p, width = 6.5, height = h, dpi = 220, limitsize = FALSE)
open_linear_pdf <- function(path, height) {
  if (capabilities('cairo')) cairo_pdf(path, width = 6.5, height = height, onefile = TRUE)
  else pdf(path, width = 6.5, height = height, onefile = TRUE)
}
open_linear_pdf(file.path(linear_dir, 'Fig4b_overview.pdf'), h)
tryCatch(print(p), finally = dev.off())
open_linear_pdf(file.path(linear_dir, 'Fig4b_all_proteins.pdf'), 3 + rows_per_page * .23)
tryCatch({
  if (!nrow(linear_results)) print(draw_linear_heatmap(linear_results, 'No proteins'))
  else for (f in mapping$food) {
    z <- linear_results[linear_results$food == f, ]
    if (!nrow(z)) next
    pages <- split(seq_len(nrow(z)), ceiling(seq_len(nrow(z)) / rows_per_page))
    for (k in seq_along(pages)) print(draw_linear_heatmap(z[pages[[k]], ],
      paste0(plot_names[[f]], ': ', nrow(z), ' proteins | page ', k, '/', length(pages))))
  }
}, finally = dev.off())
writeLines(c(paste('Formula:', paste(deparse(model_formula), collapse = ' ')),
  paste('Base data:', base_file), paste('Protein data:', prot_file),
  'Continuous liking: all valid scores 1-9; no intake-data requirement.',
  'Beta/CI: protein input units per 1-point higher liking; ordinary lm, two-sided t tests.',
  'FDR_food: BH across selected intersection proteins within each food; FDR_all: BH across all three foods.',
  'Only food-specific same-direction limma intersection proteins are tested.',
  paste('Input proteins per food:', length(proteins)),
  'Overview is ranked by regression FDR and limits display only; complete CSV/PDF retain all pairs.'),
  file.path(linear_dir, 'analysis_notes.txt'))
capture.output(sessionInfo(), file = file.path(linear_dir, 'sessionInfo.txt'))
print(linear_counts, row.names = FALSE)
message('Fig.4b完成：', linear_dir)

})

# ==================== COX ====================
if ("cox" %in% run_modules) for (pair_index in seq_len(nrow(food_pairs))) local({
pair_current <- food_pairs[pair_index, , drop = FALSE]
selected_targets <- read_stage_targets(pair_current, "linear")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 🚩 各组交集蛋白与对应疾病结局的Cox关联和森林图
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Ys <- pair_current$outcome
covariates <- covs.all
race_var <- 'ethnic.c3'     # 仅用于审计；主模型与s1_1一致，不调整族群
resume <- resume_models
checkpoint_every <- 100L
alpha <- .05
min_events <- 10L          # 自定计算门槛，不表示模型稳定性已得到保证

# 默认HR为蛋白每增加1 SD；可改成'input'表示每增加1个输入单位。
cox_scale <- 'SD'
run_PH <- TRUE              # 比例风险假设诊断；单独记录，不按其P值删蛋白
forest_top_n <- 10L         # 概览每结局前10个；完整PDF保留全部有效模型
forest_rows_per_page <- 30L
outdir <- file.path(project, "01protein/cox", pair_current$pair_id)
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
write.csv(selected_targets, file.path(outdir, "selected_targets.csv"), row.names = FALSE)
if (!nrow(selected_targets)) {
  record_empty_stage(outdir, paste0(pair_current$pair_id, ": 无第2步线性回归FDR显著蛋白，跳过cox")); return(invisible(NULL))
}
write.csv(data.frame(status = "selected", n_selected = nrow(selected_targets)),
  file.path(outdir, "stage_status.csv"), row.names = FALSE)
model_version <- 'flt_s1_4_cox_sequential_v3'
model_options <- list(cox_scale = cox_scale, run_PH = run_PH)

# 2. 数据、风险集和共用函数 -------------------------------------------
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
checkpoint_dir <- file.path(outdir, 'checkpoints')
dir.create(checkpoint_dir, showWarnings = FALSE)
write_csv <- function(x, name) write.csv(x, file.path(outdir, name), row.names = FALSE, na = '')
base <- read_data(base_file)
prot <- read_data(prot_file)
proteins <- unique(selected_targets$protein)
if (!all(proteins %in% names(prot))) stop('交集蛋白不在prot中，请重跑limma')
race_included <- FALSE  # 与s1_1一致，不额外调整族群
if (!all(covariates %in% names(base))) stop("缺少项目协变量")
if (!length(proteins) || !all(vapply(prot[proteins], is.numeric, logical(1)))) stop("蛋白列须为数值")
fields <- as.vector(outer(Ys, c('.eligible', '.t2e', '.Yt2e'), paste0))
if (!all(fields %in% names(base))) stop('base缺少结局字段：', paste(setdiff(fields, names(base)), collapse = ', '))
dat <- base[base$eid %in% prot$eid, unique(c('eid', covariates, fields)), drop = FALSE]
if (!nrow(dat)) stop('base和prot没有匹配的eid')
protein_data <- prot[match(dat$eid, prot$eid), proteins, drop = FALSE]
for (v in proteins) protein_data[[v]][!is.finite(protein_data[[v]])] <- NA_real_
factor_vars <- c("sex", "edu", "inc", "smoke_status", "alcohol_status")
for (v in factor_vars) dat[[v]] <- factor(dat[[v]])
for (v in setdiff(covariates, factor_vars)) {
  if (!is.numeric(dat[[v]])) stop('连续协变量必须为数值：', v)
  dat[[v]][!is.finite(dat[[v]])] <- NA_real_
}
for (Y in Ys) {
  eligible <- dat[[paste0(Y, '.eligible')]]
  time <- dat[[paste0(Y, '.t2e')]]; event <- dat[[paste0(Y, '.Yt2e')]]
  if (!is.logical(eligible) || anyNA(eligible) || !is.numeric(time) || !is.numeric(event))
    stop('无效eligible、随访时间或事件类型：', Y)
  if (any(!is.finite(time[eligible]) | time[eligible] <= 0) || anyNA(event[eligible]) || any(!event[eligible] %in% 0:1)) stop('无效随访/事件：', Y)
}
rm(base, prot); invisible(gc())
write_csv(data.frame(variable = unique(c(covariates, if (is.null(race_var)) 'race' else race_var)),
  included = unique(c(covariates, if (is.null(race_var)) 'race' else race_var)) %in% covariates),
  'covariate_audit.csv')
# 恢复时比较输入文件路径/大小/修改时间与模型设置，防止混合不同分析。
settings <- list(version = model_version, base = file_stamp(base_file), prot = file_stamp(prot_file),
  pair = pair_current, Ys = Ys, proteins = proteins, covariates = covariates, race_included = race_included,
  min_events = min_events, model_options = model_options)
saveRDS(settings, file.path(outdir, 'analysis_settings.rds'))
writeLines(c(paste('Covariates:', paste(covariates, collapse = ', ')),
  paste('Race included:', race_included),
  'Only step-2 linear FDR_food < 0.05 proteins from the food-specific limma intersection.',
  'Outcome-specific eligible risk sets; complete cases for each protein and model covariates.',
  'No new imputation; continuous covariates remain continuous.',
  'Resume checkpoints compare file path/size/mtime and model options; failed models are also retained.',
  'Set resume_models <- FALSE after changing inputs/options or to retry failed models.'),
  file.path(outdir, 'analysis_notes.txt'))
model_data <- function(Y, gene) {
  keep <- dat[[paste0(Y, '.eligible')]]
  d <- dat[keep, covariates, drop = FALSE]
  d$.time <- dat[[paste0(Y, '.t2e')]][keep]
  d$.event <- dat[[paste0(Y, '.Yt2e')]][keep]
  d$.protein <- protein_data[[gene]][keep]
  droplevels(d[complete.cases(d), , drop = FALSE])
}
fit_cox <- function(formula, d) {
  fit <- withCallingHandlers(coxph(formula, data = d, ties = 'efron', na.action = na.fail,
    singular.ok = FALSE, x = TRUE, control = coxph.control(iter.max = 50)),
    warning = function(w) stop(conditionMessage(w), call. = FALSE))
  if (any(!is.finite(coef(fit))) || any(!is.finite(vcov(fit)))) stop('invalid_cox_estimate')
  fit
}
validate_model_data <- function(d) {
  if (nrow(d) < 3L || sum(d$.event) < min_events) stop('insufficient_samples_or_events')
  if (!is.finite(sd(d$.protein)) || sd(d$.protein) <= 0) stop('constant_protein')
}
load_checkpoint <- function(Y) {
  path <- file.path(checkpoint_dir, paste0(Y, '.rds'))
  if (resume && file.exists(path)) {
    state <- readRDS(path)
    if (!identical(state$settings, settings))
      stop('断点与当前输入/设置不符：', path, '；请设resume_models <- FALSE重算')
    message(Y, ': resume ', length(state$rows), '/', length(proteins))
    state
  } else list(settings = settings, rows = list(), curves = list())
}
save_checkpoint <- function(state, Y) {
  path <- file.path(checkpoint_dir, paste0(Y, '.rds'))
  tmp <- paste0(path, '.tmp')
  saveRDS(state, tmp, compress = FALSE)
  if (!file.rename(tmp, path)) {
    if (!file.copy(tmp, path, overwrite = TRUE)) stop('断点保存失败：', path)
    unlink(tmp)
  }
}
bh <- function(p, n) {
  answer <- rep(NA_real_, length(p)); ok <- is.finite(p)
  if (any(ok)) answer[ok] <- p.adjust(p[ok], 'BH', n = n)
  answer
}
open_pdf <- function(path, width = 10, height = 9) {
  if (capabilities('cairo')) cairo_pdf(path, width = width, height = height, onefile = TRUE)
  else pdf(path, width = width, height = height, onefile = TRUE)
}
format_p <- function(p) ifelse(is.finite(p), formatC(p, digits = 2, format = 'g'), 'NA')

# 3. 全蛋白逐个Cox模型；每100个保存断点 -------------------------------
cox_formula <- reformulate(c('.z', covariates), response = 'Surv(.time, .event)')
run_one <- function(Y, gene) {
  d <- model_data(Y, gene)
  a <- data.frame(pair_id = pair_current$pair_id, disease_id = pair_current$disease_id, food = pair_current$food, outcome = Y, protein = gene, n = nrow(d), events = sum(d$.event),
    person_years = sum(d$.time), mean_protein = if (nrow(d)) mean(d$.protein) else NA_real_,
    sd_protein = if (nrow(d) > 1L) sd(d$.protein) else NA_real_,
    unit = if (cox_scale == 'SD') 'per 1 SD' else 'per input unit',
    beta = NA_real_, SE = NA_real_, HR = NA_real_, CI_low = NA_real_, CI_high = NA_real_,
    P = NA_real_, PH_P = NA_real_, PH_global_P = NA_real_, PH_status = 'not_run',
    status = 'ok', reason = '', covariates = paste(covariates, collapse = ';'),
    race_included = race_included)
  tryCatch({
    validate_model_data(d)
    d$.z <- if (cox_scale == 'SD') (d$.protein - a$mean_protein) / a$sd_protein else d$.protein
    fit <- fit_cox(cox_formula, d)
    s <- summary(fit)
    a$beta <- s$coefficients['.z', 'coef']; a$SE <- s$coefficients['.z', 'se(coef)']
    a$P <- s$coefficients['.z', 'Pr(>|z|)']
    a$HR <- s$conf.int['.z', 'exp(coef)']
    a$CI_low <- s$conf.int['.z', 'lower .95']; a$CI_high <- s$conf.int['.z', 'upper .95']
    if (!all(is.finite(unlist(a[c('beta','SE','P','HR','CI_low','CI_high')]))) ||
        a$SE <= 0 || a$CI_low <= 0) stop('invalid_estimate_or_CI')
    if (run_PH) {
      ph <- tryCatch(cox.zph(fit), warning = function(w) w, error = function(e) e)
      if (inherits(ph, 'condition')) a$PH_status <- conditionMessage(ph)
      else {
        a$PH_P <- ph$table['.z', 'p']; a$PH_global_P <- ph$table['GLOBAL', 'p']
        a$PH_status <- 'ok'
      }
    }
    a
  }, error = function(e) {
    a$status <- 'failed'; a$reason <- conditionMessage(e); a
  })
}
all_rows <- list()
for (Y in Ys) {
  state <- load_checkpoint(Y)
  todo <- setdiff(proteins, names(state$rows))
  for (i in seq_along(todo)) {
    gene <- todo[i]
    state$rows[[gene]] <- run_one(Y, gene)
    if (i %% checkpoint_every == 0L || i == length(todo)) {
      save_checkpoint(state, Y)
      message('Cox ', Y, ': ', length(state$rows), '/', length(proteins))
    }
  }
  all_rows[[Y]] <- do.call(rbind, state$rows[proteins])
}
results <- do.call(rbind, all_rows)
results$FDR_outcome <- NA_real_
for (Y in Ys) {
  j <- results$outcome == Y
  results$FDR_outcome[j] <- bh(ifelse(results$status[j] == 'ok', results$P[j], NA_real_), length(proteins))
}
results$FDR_all <- bh(ifelse(results$status == 'ok', results$P, NA_real_), nrow(results))
results$significant <- !is.na(results$FDR_outcome) & results$FDR_outcome < alpha
write_csv(results, 'protein_cox_all.csv')
write_csv(results[results$significant, ], 'protein_cox_significant.csv')
write_csv(results[results$status != 'ok', ], 'protein_cox_failed.csv')
write_csv(do.call(rbind, lapply(Ys, function(Y) {
  z <- results[results$outcome == Y, ]
  data.frame(outcome = Y, n_proteins = nrow(z), n_fitted = sum(z$status == 'ok'),
    n_failed = sum(z$status != 'ok'), n_significant = sum(z$significant),
    n_higher_risk = sum(z$significant & z$HR > 1, na.rm = TRUE),
    n_lower_risk = sum(z$significant & z$HR < 1, na.rm = TRUE))
})), 'protein_cox_counts.csv')

# 4. 自动生成森林图；使用数值HR/CI，无需人工整理Excel ------------------
forest_plot <- function(z, title) {
  if (!nrow(z)) return(ggplot() + annotate('text', x = 0, y = 0, label = 'No valid Cox models') + theme_void() + labs(title = title))
  z$protein <- factor(z$protein, levels = rev(z$protein))
  z$label <- sprintf('%.2f (%.2f, %.2f) | %s', z$HR, z$CI_low, z$CI_high, format_p(z$FDR_outcome))
  z$direction <- ifelse(!z$significant, 'Other', ifelse(z$HR > 1, 'Higher risk', 'Lower risk'))
  z$case_label <- paste0(z$protein, '  [', z$events, '/', z$n, ']')
  labs_y <- setNames(z$case_label, as.character(z$protein))
  ggplot(z, aes(HR, protein, color = direction)) +
    geom_vline(xintercept = 1, linetype = 2, color = 'grey40') +
    geom_segment(aes(x = CI_low, xend = CI_high, yend = protein), linewidth = .45) +
    geom_point(size = 1.6) +
    scale_y_discrete(labels = labs_y) + scale_x_log10(expand = expansion(mult = c(.05, .05))) +
    coord_cartesian(clip = 'off') +
    geom_text(aes(x = max(CI_high) * 1.05, label = label), hjust = 0, color = 'black', size = 2.6) +
    scale_color_manual(values = c('Higher risk' = '#C84D4D', 'Lower risk' = '#425D8B', Other = 'grey60')) +
    labs(title = title, subtitle = 'Protein [cases/N] | right: HR (95% CI) | BH-FDR',
      x = paste('Hazard ratio', if (cox_scale == 'SD') 'per 1 SD higher protein' else 'per input unit higher protein'),
      y = NULL, color = NULL,
      caption = paste0('Dashed line: HR = 1. BH-FDR within each outcome. Race adjusted: ', race_included)) +
    theme_classic(base_size = 10) + theme(legend.position = 'bottom',
      plot.margin = margin(10, 140, 10, 10),
      plot.background = element_rect(fill = 'white', color = NA))
}
overviews <- list()
for (Y in Ys) {
  z <- results[results$outcome == Y & results$status == 'ok', ]
  z <- z[order(z$FDR_outcome, z$P, z$protein), ]
  top <- head(z, forest_top_n); overviews[[Y]] <- top
  p <- forest_plot(top, paste0(Y, ' | Top ', forest_top_n, ' by BH-FDR'))
  ggsave(file.path(outdir, paste0('forest_', Y, '_overview.png')), p, width = 11,
    height = max(4, 2 + nrow(top) * .3), dpi = 220, bg = 'white')
  open_pdf(file.path(outdir, paste0('forest_', Y, '_all.pdf')), 11, 10)
  tryCatch({
    if (!nrow(z)) print(forest_plot(z, Y))
    else for (j in seq(1L, nrow(z), by = forest_rows_per_page))
      print(forest_plot(z[j:min(j + forest_rows_per_page - 1L, nrow(z)), ], Y))
  }, finally = dev.off())
}
write_csv(do.call(rbind, overviews), 'forest_overview_rows.csv')
capture.output(sessionInfo(), file = file.path(outdir, 'sessionInfo.txt'))
message('Cox及森林图完成：', outdir)

})

# ==================== RCS ====================
if ("rcs" %in% run_modules) for (pair_index in seq_len(nrow(food_pairs))) local({
pair_current <- food_pairs[pair_index, , drop = FALSE]
selected_targets <- read_stage_targets(pair_current, "cox")
# Fig.4d扩展：各组交集蛋白与对应疾病结局的5结点RCS；本文件是R脚本
# 按三组配对分别分析；依赖第2步线性回归及第3步Cox显著结果。
# 沿用00data中每个结局自身的eligible、随访时间和事件；不要求偏好或摄入数据。
# 连续协变量保留连续形式；分类变量沿用上游整理值；不在此插补或重复log2。

# 1. 路径和参数：两个脚本的协变量、输入路径应保持一致 -----------------
Ys <- pair_current$outcome
covariates <- covs.all
race_var <- 'ethnic.c3'     # 仅用于审计；主模型与s1_1一致，不调整族群
resume <- resume_models
checkpoint_every <- 100L
alpha <- .05              # RCS使用原始P值，不做多重检验校正
min_events <- 10L          # 自定计算门槛，不表示模型稳定性已得到保证

rcs_probs <- c(.05, .25, .50, .75, .95)  # 论文方法的5个结点；不自动减少结点
rcs_grid_n <- 151L
rcs_overview_n <- 3L       # 每结局按RCS总体P选前3个，仅用于概览
outdir <- file.path(project, "01protein/rcs", pair_current$pair_id)
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
write.csv(selected_targets, file.path(outdir, "selected_targets.csv"), row.names = FALSE)
if (!nrow(selected_targets)) {
  record_empty_stage(outdir, paste0(pair_current$pair_id, ": 无第3步Cox FDR显著蛋白，跳过rcs")); return(invisible(NULL))
}
stopifnot(length(rcs_probs) == 5L, all(diff(rcs_probs) > 0),
          rcs_grid_n >= 20L, checkpoint_every >= 1L, rcs_overview_n >= 1L)
write.csv(data.frame(status = "selected", n_selected = nrow(selected_targets)),
  file.path(outdir, "stage_status.csv"), row.names = FALSE)
model_version <- 'flt_s1_4_rcs_sequential_v3'
model_options <- list(rcs_probs = rcs_probs, rcs_grid_n = rcs_grid_n)

# 2. 数据、风险集和共用函数 -------------------------------------------
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
checkpoint_dir <- file.path(outdir, 'checkpoints')
dir.create(checkpoint_dir, showWarnings = FALSE)
write_csv <- function(x, name) write.csv(x, file.path(outdir, name), row.names = FALSE, na = '')
base <- read_data(base_file)
prot <- read_data(prot_file)
proteins <- unique(selected_targets$protein)
if (!all(proteins %in% names(prot))) stop('交集蛋白不在prot中，请重跑limma')
if (!length(proteins) || !all(vapply(prot[proteins], is.numeric, logical(1))))
  stop('prot中除eid外必须全部为数值蛋白列，且至少有一个蛋白')
if (!all(covariates %in% names(base)))
  stop('base缺少协变量：', paste(setdiff(covariates, names(base)), collapse = ', '))
race_included <- FALSE  # 与s1_1一致，不额外调整族群
if (!all(covariates %in% names(base))) stop("缺少项目协变量")
if (!length(proteins) || !all(vapply(prot[proteins], is.numeric, logical(1)))) stop("蛋白列须为数值")
fields <- as.vector(outer(Ys, c('.eligible', '.t2e', '.Yt2e'), paste0))
if (!all(fields %in% names(base)))
  stop('base缺少结局字段：', paste(setdiff(fields, names(base)), collapse = ', '))
dat <- base[base$eid %in% prot$eid, unique(c('eid', covariates, fields)), drop = FALSE]
if (!nrow(dat)) stop('base和prot没有匹配的eid')
protein_data <- prot[match(dat$eid, prot$eid), proteins, drop = FALSE]
for (v in proteins) protein_data[[v]][!is.finite(protein_data[[v]])] <- NA_real_
factor_vars <- c("sex", "edu", "inc", "smoke_status", "alcohol_status")
for (v in factor_vars) dat[[v]] <- factor(dat[[v]])
for (v in setdiff(covariates, factor_vars)) {
  if (!is.numeric(dat[[v]])) stop('连续协变量必须为数值：', v)
  dat[[v]][!is.finite(dat[[v]])] <- NA_real_
}
for (Y in Ys) {
  eligible <- dat[[paste0(Y, '.eligible')]]
  time <- dat[[paste0(Y, '.t2e')]]; event <- dat[[paste0(Y, '.Yt2e')]]
  if (!is.logical(eligible) || anyNA(eligible) || !is.numeric(time) || !is.numeric(event))
    stop('无效eligible、随访时间或事件类型：', Y)
  if (any(!is.finite(time[eligible]) | time[eligible] <= 0) ||
      anyNA(event[eligible]) || any(!event[eligible] %in% 0:1)) stop('无效随访/事件：', Y)
}
rm(base, prot); invisible(gc())
write_csv(data.frame(variable = unique(c(covariates, if (is.null(race_var)) 'race' else race_var)),
  included = unique(c(covariates, if (is.null(race_var)) 'race' else race_var)) %in% covariates),
  'covariate_audit.csv')
# 恢复时比较输入文件路径/大小/修改时间与模型设置，防止混合不同分析。
settings <- list(version = model_version, base = file_stamp(base_file), prot = file_stamp(prot_file),
  pair = pair_current, Ys = Ys, proteins = proteins, covariates = covariates, race_included = race_included,
  min_events = min_events, model_options = model_options)
saveRDS(settings, file.path(outdir, 'analysis_settings.rds'))
writeLines(c(paste('Covariates:', paste(covariates, collapse = ', ')),
  paste('Race included:', race_included),
  'Only step-3 Cox FDR_outcome < 0.05 proteins; all valid curves plotted regardless of RCS P.',
  'RCS significance: raw overall/nonlinear LRT P < 0.05; no FDR adjustment. Overview ranked by raw overall P.',
  'Outcome-specific eligible risk sets; complete cases for each protein and model covariates.',
  'No new imputation; continuous covariates remain continuous.',
  'Resume checkpoints compare file path/size/mtime and model options; failed models are also retained.',
  'Set resume_models <- FALSE after changing inputs/options or to retry failed models.'),
  file.path(outdir, 'analysis_notes.txt'))
model_data <- function(Y, gene) {
  keep <- dat[[paste0(Y, '.eligible')]]
  d <- dat[keep, covariates, drop = FALSE]
  d$.time <- dat[[paste0(Y, '.t2e')]][keep]
  d$.event <- dat[[paste0(Y, '.Yt2e')]][keep]
  d$.protein <- protein_data[[gene]][keep]
  droplevels(d[complete.cases(d), , drop = FALSE])
}
fit_cox <- function(formula, d) {
  fit <- withCallingHandlers(coxph(formula, data = d, ties = 'efron', na.action = na.fail,
    singular.ok = FALSE, x = TRUE, control = coxph.control(iter.max = 50)),
    warning = function(w) stop(conditionMessage(w), call. = FALSE))
  if (any(!is.finite(coef(fit))) || any(!is.finite(vcov(fit)))) stop('invalid_cox_estimate')
  fit
}
validate_model_data <- function(d) {
  if (nrow(d) < 3L || sum(d$.event) < min_events) stop('insufficient_samples_or_events')
  if (!is.finite(sd(d$.protein)) || sd(d$.protein) <= 0) stop('constant_protein')
}
load_checkpoint <- function(Y) {
  path <- file.path(checkpoint_dir, paste0(Y, '.rds'))
  if (resume && file.exists(path)) {
    state <- readRDS(path)
    if (!identical(state$settings, settings))
      stop('断点与当前输入/设置不符：', path, '；请设resume_models <- FALSE重算')
    message(Y, ': resume ', length(state$rows), '/', length(proteins))
    state
  } else list(settings = settings, rows = list(), curves = list())
}
save_checkpoint <- function(state, Y) {
  path <- file.path(checkpoint_dir, paste0(Y, '.rds'))
  tmp <- paste0(path, '.tmp')
  saveRDS(state, tmp, compress = FALSE)
  if (!file.rename(tmp, path)) {
    if (!file.copy(tmp, path, overwrite = TRUE)) stop('断点保存失败：', path)
    unlink(tmp)
  }
}
open_pdf <- function(path, width = 10, height = 9) {
  if (capabilities('cairo')) cairo_pdf(path, width = width, height = height, onefile = TRUE)
  else pdf(path, width = width, height = height, onefile = TRUE)
}
format_p <- function(p) ifelse(is.finite(p), formatC(p, digits = 2, format = 'g'), 'NA')

# 3. 5结点RCS：同一完整样本拟合协变量模型、线性模型和样条模型 ----------
# ns的三个内部结点+两个边界结点，与5结点RCS为同一自然三次样条函数空间。
# 总体检验：协变量模型 vs 样条模型（4 df）；非线性：线性 vs 样条（3 df）。
# 另附线性模型 vs 协变量模型的1 df检验；均为部分似然比检验（ANOVA）。
null_formula <- reformulate(covariates, response = 'Surv(.time, .event)')
linear_formula <- reformulate(c('.z', covariates), response = 'Surv(.time, .event)')
lrt <- function(small, large) {
  df <- length(coef(large)) - length(coef(small))
  if (df <= 0) stop('nonpositive_LRT_df')
  statistic <- 2 * (tail(large$loglik, 1L) - tail(small$loglik, 1L))
  if (statistic < -1e-5) stop('nested_model_loglik_inconsistent')
  c(P = pchisq(max(0, statistic), df, lower.tail = FALSE), df = df)
}
run_rcs <- function(Y, gene) {
  d <- model_data(Y, gene)
  a <- data.frame(pair_id = pair_current$pair_id, disease_id = pair_current$disease_id, food = pair_current$food, outcome = Y, protein = gene, n = nrow(d), events = sum(d$.event),
    reference = NA_real_, knot05 = NA_real_, knot25 = NA_real_, knot50 = NA_real_,
    knot75 = NA_real_, knot95 = NA_real_, P_linear = NA_real_, P_overall = NA_real_,
    P_nonlinear = NA_real_, df_linear = NA_integer_, df_overall = NA_integer_,
    df_nonlinear = NA_integer_, status = 'ok', reason = '',
    covariates = paste(covariates, collapse = ';'), race_included = race_included)
  curve <- NULL
  tryCatch({
    validate_model_data(d)
    knots <- as.numeric(quantile(d$.protein, rcs_probs, names = FALSE))
    if (any(diff(knots) <= 0)) stop('tied_spline_knots')
    a[1L, c('knot05','knot25','knot50','knot75','knot95')] <- as.list(knots)
    a$reference <- median(d$.protein)
    d$.z <- as.numeric(scale(d$.protein))
    basis <- splines::ns(d$.protein, knots = knots[2:4], Boundary.knots = knots[c(1,5)])
    spline_names <- paste0('.s', seq_len(ncol(basis)))
    d[spline_names] <- as.data.frame(basis)
    null <- fit_cox(null_formula, d)
    linear <- fit_cox(linear_formula, d)
    full <- fit_cox(reformulate(c(spline_names, covariates), response = 'Surv(.time, .event)'), d)
    test_linear <- lrt(null, linear); test_overall <- lrt(null, full); test_nonlinear <- lrt(linear, full)
    a$P_linear <- test_linear['P']; a$df_linear <- test_linear['df']
    a$P_overall <- test_overall['P']; a$df_overall <- test_overall['df']
    a$P_nonlinear <- test_nonlinear['P']; a$df_nonlinear <- test_nonlinear['df']
    # 模型使用全部观测；仅将绘图范围限定为第5至95百分位。中位数为HR=1参考。
    xx <- sort(unique(c(seq(knots[1], knots[5], length.out = rcs_grid_n), a$reference)))
    contrast <- sweep(predict(basis, xx), 2L, as.numeric(predict(basis, a$reference)), '-')
    b <- coef(full)[spline_names]; V <- vcov(full)[spline_names, spline_names, drop = FALSE]
    eta <- as.numeric(contrast %*% b)
    SE <- sqrt(pmax(0, rowSums((contrast %*% V) * contrast)))
    curve <- data.frame(pair_id = pair_current$pair_id, disease_id = pair_current$disease_id, food = pair_current$food, outcome = Y, protein = gene, x = xx,
      HR = exp(eta), CI_low = exp(eta - 1.96 * SE), CI_high = exp(eta + 1.96 * SE))
    if (any(!is.finite(as.matrix(curve[c('HR','CI_low','CI_high')])))) stop('invalid_curve')
    list(row = a, curve = curve)
  }, error = function(e) {
    a$status <- 'failed'; a$reason <- conditionMessage(e)
    list(row = a, curve = NULL)
  })
}
all_tests <- list()
for (Y in Ys) {
  state <- load_checkpoint(Y)
  todo <- setdiff(proteins, names(state$rows))
  for (i in seq_along(todo)) {
    gene <- todo[i]; ans <- run_rcs(Y, gene)
    state$rows[[gene]] <- ans$row
    if (!is.null(ans$curve)) state$curves[[gene]] <- ans$curve
    if (i %% checkpoint_every == 0L || i == length(todo)) {
      save_checkpoint(state, Y)
      message('RCS ', Y, ': ', length(state$rows), '/', length(proteins))
    }
  }
  tests <- do.call(rbind, state$rows[proteins])
  all_tests[[Y]] <- tests
  write_csv(tests, paste0('protein_rcs_tests_', Y, '.csv'))
  curves <- if (length(state$curves)) do.call(rbind, state$curves) else
    data.frame(outcome = character(), protein = character(), x = numeric(),
               HR = numeric(), CI_low = numeric(), CI_high = numeric())
  write_csv(curves, paste0('protein_rcs_curves_', Y, '.csv'))
  # 4. 当前结局所有曲线分页PDF；另存前3个概览PNG/PDF --------------------
  draw_curves <- function(genes, title) {
    d <- curves[curves$protein %in% genes, ]
    if (!nrow(d)) return(ggplot() + annotate('text', x = 0, y = 0, label = 'No valid RCS curves') +
      theme_void() + labs(title = title))
    stat <- tests[match(genes, tests$protein), ]
    panel_labels <- paste0(stat$protein, '\nP overall=', format_p(stat$P_overall),
      '; P nonlinear=', format_p(stat$P_nonlinear))
    d$panel <- factor(d$protein, levels = genes, labels = panel_labels)
    ggplot(d, aes(x, HR)) +
      geom_hline(yintercept = 1, linetype = 2, color = 'grey40') +
      geom_ribbon(aes(ymin = CI_low, ymax = CI_high), fill = '#F1948A', alpha = .4) +
      geom_line(color = '#D55E40', linewidth = .7) +
      facet_wrap(~panel, scales = 'free', ncol = 3) +
      labs(title = title, x = 'Protein level (input scale)', y = 'HR relative to median protein',
        caption = paste0('5 knots: 5th, 25th, 50th, 75th, 95th percentiles. Pointwise 95% CI.\n',
          'Plots show 5th-95th percentiles; all observations fitted. Labels: nominal LRT P. Race adjusted: ', race_included)) +
      theme_classic(base_size = 10) + theme(plot.background = element_rect(fill = 'white', color = NA),
        strip.background = element_blank(), strip.text = element_text(size = 9))
  }
  valid <- tests[tests$status == 'ok', ]
  valid <- valid[order(valid$P_overall, valid$protein), ]
  overview <- head(valid, rcs_overview_n)
  write_csv(overview, paste0('RCS_', Y, '_overview_rows.csv'))
  p <- draw_curves(overview$protein, paste0(Y, ' | Top ', rcs_overview_n, ' RCS associations by raw overall P'))
  ggsave(file.path(outdir, paste0('RCS_', Y, '_overview.png')), p, width = 13, height = 5,
    dpi = 220, bg = 'white')
  open_pdf(file.path(outdir, paste0('RCS_', Y, '_overview.pdf')), 13, 5)
  tryCatch(print(p), finally = dev.off())
  open_pdf(file.path(outdir, paste0('RCS_', Y, '_all.pdf')), 13, 8)
  tryCatch({
    if (!nrow(valid)) print(draw_curves(character(), Y))
    else for (j in seq(1L, nrow(valid), by = 6L))
      print(draw_curves(valid$protein[j:min(j + 5L, nrow(valid))], Y))
  }, finally = dev.off())
  rm(state, curves); invisible(gc())
}
results <- do.call(rbind, all_tests)
results$significant_overall <- results$status == 'ok' & is.finite(results$P_overall) & results$P_overall < alpha
results$significant_nonlinear <- results$status == 'ok' & is.finite(results$P_nonlinear) & results$P_nonlinear < alpha
write_csv(results, 'protein_rcs_tests_all.csv')
write_csv(results[results$significant_overall, ], 'protein_rcs_significant_overall.csv')
write_csv(results[results$status != 'ok', ], 'protein_rcs_failed.csv')
write_csv(do.call(rbind, lapply(Ys, function(Y) {
  z <- results[results$outcome == Y, ]
  data.frame(outcome = Y, n_proteins = nrow(z), n_fitted = sum(z$status == 'ok'),
    n_failed = sum(z$status != 'ok'), n_P_overall = sum(z$significant_overall),
    n_P_nonlinear = sum(z$significant_nonlinear))
})), 'protein_rcs_counts.csv')
capture.output(sessionInfo(), file = file.path(outdir, 'sessionInfo.txt'))
message('交集蛋白RCS完成：', outdir)

})
