library(dplyr)

dir0 <- "D:"
project <- file.path(dir0, "projects/flt")
input <- file.path(project, "analysis/00data")
outdir <- file.path(project, "analysis/01lm")
dir.create(outdir, recursive=TRUE, showWarnings=FALSE)

cfg <- readRDS(file.path(input, "analysis_config.rds"))
liking_selected <- c(`Butter on bread`="p20628",`Total salad vegetable`="F-Salad_vegetables", Chicken="p20638")
flt_vars <- unname(liking_selected)
covs.all <- c("age", "sex", "bmi", "edu", "inc", "tdi", "smoke_status", "alcohol_status", "pa_met_tot")

dat <- as.data.frame(readRDS(file.path(input, cfg$primary_data_file)))
lifestyle <- as.data.frame(readRDS(file.path(dir0, "data/phe/ukb/phe/Rdata/lifestyle_36vars.rds")))
brain <- readRDS(file.path(dir0, "data/phe/ukb/phe/Rdata/behavioral.rds"))
brain <- as.data.frame(brain)

required <- c("eid", flt_vars, covs.all)
if (!all(required %in% names(dat))) stop("缺少分析变量：", paste(setdiff(required, names(dat)), collapse=", "))
dat <- dat %>% select(all_of(required))
body_vars <- c("fat_pct", "whole_body_fat_mass", "bmi",
  paste0("leg_", c("fat_pct", "fat_mass", "fat_free_mass", "predicted_mass"), "_right"),
  paste0("leg_", c("fat_pct", "fat_mass", "fat_free_mass", "predicted_mass"), "_left"),
  paste0("arm_", c("fat_pct", "fat_mass", "fat_free_mass", "predicted_mass"), "_right"),
  paste0("arm_", c("fat_pct", "fat_mass", "fat_free_mass", "predicted_mass"), "_left"))
lifestyle_vars <- c("smoke_status", "alcohol_status", "alcohol_freq", "sleep_duration",
  "phone_call_time", "tv_time", "computer_time", "pa_met_tot", "summed_mvpa", "walking_duration",
  "walk_time", "liking_gym", "liking_stairs", "liking_sweat", "education_degree", "tdi", "inc")
brain_continuous <- paste0(c("fluid_intelligence", "reaction_time", "pairs_matching", "numeric_memory", "life_satisfaction", "subjective_wellbeing", "phq4", "phq9", "cidi", "gad7"), ".z")
brain_binary <- c("prospective_memory.b", "general_mental_health.b")
traits <- data.frame(trait=c(body_vars, lifestyle_vars, brain_continuous, brain_binary),
  source=c(rep("lifestyle",36),rep("brain",12)),
  category=c(rep("Anthropometry",19),rep("Lifestyle and social factors",17), rep("Brain continuous",10),rep("Brain binary",2)))
stopifnot(all(traits$trait[traits$source=="lifestyle"] %in% names(lifestyle)), all(traits$trait[traits$source=="brain"] %in% names(brain)))
ordinal <- c("smoke_status", "alcohol_status", "alcohol_freq", "phone_call_time", "walk_time", "inc")
traits$model <- ifelse(traits$trait %in% c("education_degree",brain_binary), "Linear probability (HC1)", ifelse(traits$trait %in% ordinal, "Numeric category score (HC1)", "Linear regression (HC1)"))
# 自身或同一构念的替代测量不再作调整项；其余保持配置中的协变量。
# 身体组成仍调整 BMI；因此解释为相同 BMI 下的关联。
exclude_covs <- list(bmi="bmi", smoke_status="smoke_status", alcohol_status="alcohol_status",
  alcohol_freq="alcohol_status", pa_met_tot="pa_met_tot", summed_mvpa="pa_met_tot",
  walking_duration="pa_met_tot", walk_time="pa_met_tot", education_degree="edu", tdi="tdi", inc="inc")
for (nm in c("dat", "lifestyle", "brain")) {
  z <- get(nm)
  z$eid <- as.character(z$eid)
  stopifnot(!anyNA(z$eid), all(nzchar(z$eid)), !anyDuplicated(z$eid))
  assign(nm, z)
}
# 按 eid 映射，避免同名结局和协变量 join 后产生 .x/.y。
idx_lifestyle <- match(dat$eid, lifestyle$eid)
idx_brain <- match(dat$eid, brain$eid)
cat("分析人群：", nrow(dat), "; FLT：", length(flt_vars), "; 表型：", nrow(traits), "\n")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Adjusted linear models
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# beta：FLT 每增加 1 分，对应结局原单位的平均差；.z 为上游秩-logit标准化单位。
# 二分类 beta 为概率差；等级变量假设相邻编码等距，不解释为真实剂量。
# 所有模型使用 HC1 异方差稳健标准误和残差自由度的双侧 t 检验。
# 不插补 FLT/结局；逐对使用完整数据。在线问卷时间不同，不作因果或前瞻性解释。
fit_pair <- function(base, x, covs) {
  d <- base
  d$.flt <- x
  for (nm in names(d)) if (is.numeric(d[[nm]])) d[[nm]][!is.finite(d[[nm]])] <- NA_real_
  d <- droplevels(d[complete.cases(d), , drop=FALSE])
  ans <- data.frame(n=nrow(d), beta=NA_real_, se=NA_real_, CI_low=NA_real_, CI_high=NA_real_, p=NA_real_, status="ok", message="")
  warnings <- character()
  fit <- tryCatch(withCallingHandlers({
    if (nrow(d)<3L || length(unique(d$.flt))<2L || length(unique(d$.outcome))<2L) stop("样本不足或暴露/结局无变化")
    constant <- covs[vapply(d[covs], function(z) length(unique(z))<2L, logical(1))]
    if (length(constant)) stop("协变量无变化：",paste(constant,collapse=","))
    f <- reformulate(c(".flt", covs), response=".outcome")
    fit <- lm(f, data=d, na.action=na.fail, x=TRUE, y=FALSE)
    if (fit$rank<ncol(fit$x) || fit$df.residual<=0L) stop("设计矩阵秩不足或残差自由度不足")
    # HC1 的 FLT 方差：n/(n-p) * sum(e^2 * (X (X'X)^-1[,flt])^2)。
    bread <- summary(fit)$cov.unscaled
    influence <- drop(fit$x[,rownames(bread),drop=FALSE] %*% bread[,".flt"])
    se <- sqrt(nrow(d)/fit$df.residual * sum((residuals(fit)*influence)^2))
    beta <- unname(coef(fit)[".flt"])
    if (!is.finite(se) || se<=0 || !is.finite(beta)) stop("无效估计")
    critical <- qt(0.975, fit$df.residual)
    c(beta=beta, se=se, CI_low=beta-critical*se, CI_high=beta+critical*se, p=2*pt(-abs(beta/se), df=fit$df.residual))
  }, warning=function(w) {warnings <<- c(warnings, conditionMessage(w)); invokeRestart("muffleWarning")}),
  error=function(e) e)
  if (inherits(fit,"error")) {ans$status <- "error"; ans$message <- conditionMessage(fit)
  } else {
    ans[names(fit)] <- as.list(fit)
    if (length(warnings)) {ans$status <- "warning"; ans$message <- paste(unique(warnings),collapse="; ")}
  }
  ans
}
cat("Adjustment covariates: ", paste(covs.all, collapse=", "), "\n")
traits$covariates <- vapply(traits$trait, function(y) paste(setdiff(covs.all, exclude_covs[[y]]), collapse=";"), character(1))
traits$excluded_covariates <- vapply(traits$trait, function(y) paste(intersect(covs.all, exclude_covs[[y]]), collapse=";"), character(1))
rows <- vector("list", nrow(traits)*length(flt_vars))
k <- 0L
for (i in seq_len(nrow(traits))) {
  Y <- traits$trait[i]
  y <- if (traits$source[i]=="lifestyle") lifestyle[[Y]][idx_lifestyle] else brain[[Y]][idx_brain]
  if (!is.numeric(y)) stop("结局不是数值型：",Y)
  cat(sprintf("\n[%d/%d] 表型：%s | 表型有效人数：%s / %s\n", i, nrow(traits), Y, format(sum(is.finite(y)), big.mark=","), format(nrow(dat), big.mark=",")))
  covs <- setdiff(covs.all, exclude_covs[[Y]])
  base <- dat[covs]
  base$.outcome <- y
  for (v in flt_vars) {
    if (!is.numeric(dat[[v]])) stop("FLT 不是数值型：",v)
    k <- k+1L
    rows[[k]] <- cbind(traits[i,,drop=FALSE], data.frame(variable=v, food_liking=names(liking_selected)[match(v, liking_selected)]), fit_pair(base, dat[[v]], covs))
    cat(sprintf("  [%d/%d] %s | FLT：%s | 实际分析 N=%s | 状态：%s\n", k, length(rows), Y, v, format(rows[[k]]$n, big.mark=","), rows[[k]]$status))
  }
}
results <- bind_rows(rows)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 3. Multiple testing and save
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 主要 FDR 在每个结局的 3 个指定 FLT 内校正。
# 同时提供全部 3 × 48 = 144 项检验的全局 FDR；失败模型仍计入计划检验数。
results <- results %>% group_by(trait) %>% mutate(FDR = p.adjust(ifelse(status == "ok", p, NA_real_), method = "BH", n = length(flt_vars))) %>% ungroup() %>%
  mutate(FDR_all = p.adjust(ifelse(status == "ok", p, NA_real_), method = "BH"), significant = !is.na(FDR) & FDR < 0.05)
write.csv(results, file.path(outdir,"lm_results.csv"), row.names=FALSE, na="", fileEncoding="UTF-8")
print(aggregate(significant ~ trait, results, sum))
print(table(results$status))
cat("已保存：",outdir,"\n")
