source(file.path(Sys.getenv("FLT_CODE"), "runtime.R"))
suppressPackageStartupMessages({
  library(tidyverse)
  library(data.table)
  library(survival)
  library(readxl)
})

project <- Sys.getenv("FLT_OUT")
input <- file.path(project, "00data")
flt_file <-  Sys.getenv("FLT_DICTIONARY")
outdir <- file.path(project, "01cox")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

# m1: "age", "sex", "bmi", "edu", "inc", "tdi"
# m2: "age", "sex", "bmi", "edu", "inc", "tdi", "smoke_status", "alcohol_status", "pa_met_tot"

covs.all <- c("age", "sex", "bmi", "edu", "inc", "tdi", "smoke_status", "alcohol_status", "pa_met_tot") 
Ys <- c("acd", "stroke", "sleep_disorder", "depression", "anxiety") # 
out_map <- setNames(c("All-cause dementia", "Stroke", "Sleep disorders", "Depression", "Anxiety disorders"), Ys)
meta <- as.data.frame(readxl::read_excel(flt_file, sheet = "item182"))
meta$Food_liking <- meta$item_raw
flt_vars <- meta$name

cfg <- readRDS(file.path(input, "analysis_config.rds"))
dat <- as.data.frame(readRDS(file.path(input, cfg$primary_data_file)))
required <- c(covs.all, flt_vars, paste0(rep(Ys, each = 3), c(".eligible", ".t2e", ".Yt2e")))

model_formula <- as.formula(paste("Surv(.time, .event) ~ .flt +", paste(covs.all, collapse = " + ")))
# Continuous: original score per point. Grouped: low (1-3) is reference.
# Invalid/noninteger scores are excluded from grouped analysis only.
rows <- list()
for (Y in Ys) {
  cat("Outcome:", Y, "\n")
  keep <- dat[[paste0(Y, ".eligible")]]
  stopifnot(is.logical(keep), !anyNA(keep))
  base <- dat[keep, , drop = FALSE]
  base$.time <- base[[paste0(Y, ".t2e")]]; base$.event <- base[[paste0(Y, ".Yt2e")]]
  stopifnot(all(is.finite(base$.time) & base$.time > 0), all(base$.event %in% 0:1))
  for (v in flt_vars) {
    for (mode in c("continuous", "grouped")) {
      d <- base[c(covs.all, ".time", ".event")]
      x <- base[[v]]
      x[!is.finite(x)] <- NA_real_
      if (mode == "grouped") {
        x[!is.na(x) & !(x %in% 1:9)] <- NA_real_
        d$.flt <- cut(x, breaks = c(0, 3, 6, 9), labels = c("low", "medium", "high"))
        contrasts(d$.flt) <- contr.treatment(c("low", "medium", "high"), base = 1)
      } else d$.flt <- x
      d <- d[complete.cases(d), , drop = FALSE]
      types <- if (mode == "continuous") "continuous" else c("medium_vs_low", "high_vs_low")
      terms <- if (mode == "continuous") ".flt" else c(".fltmedium", ".flthigh")
      fit <- tryCatch({
        if (!nrow(d) || sum(d$.event) == 0 || length(unique(d$.flt)) < 2) stop("No events or no exposure variation")
        if (mode == "grouped" && any(table(d$.flt) == 0)) stop("One or more preference groups absent")
        coxph(model_formula, data = d, ties = "efron", na.action = na.fail, singular.ok = FALSE, control = coxph.control(iter.max = 50))
      }, warning = function(w) w, error = function(e) e)
      for (k in seq_along(types)) {
        ans <- data.frame(outcome = Y, name = v, analysis = types[k],
          n = nrow(d), events = sum(d$.event), person_years = sum(d$.time),
          n_low = if (mode == "grouped") sum(d$.flt == "low") else NA_integer_,
          n_medium = if (mode == "grouped") sum(d$.flt == "medium") else NA_integer_,
          n_high = if (mode == "grouped") sum(d$.flt == "high") else NA_integer_,
          beta = NA_real_, se = NA_real_, HR = NA_real_, CI_low = NA_real_, CI_high = NA_real_,
          p = NA_real_, status = "ok", message = "")
        if (inherits(fit, "condition")) {
          ans$status <- if (inherits(fit, "warning")) "warning" else "error"
          ans$message <- conditionMessage(fit)
        } else {
          s <- summary(fit)
          ans[c("beta", "se", "p")] <- as.list(s$coefficients[terms[k], c("coef", "se(coef)", "Pr(>|z|)")])
          ans[c("HR", "CI_low", "CI_high")] <- as.list(s$conf.int[terms[k], c("exp(coef)", "lower .95", "upper .95")])
          if (!all(is.finite(unlist(ans[c("beta", "se", "HR", "CI_low", "CI_high", "p")]))) || ans$se <= 0 || ans$CI_low <= 0) ans$status <- "invalid_estimate"
        }
        rows[[length(rows) + 1L]] <- ans
      }
    }
  }
}
results <- do.call(rbind, rows)
# Each disease AND contrast is a separate 182-test family; failed fits retain NA.
results$FDR <- NA_real_
for (Y in Ys) for (a in unique(results$analysis)) {
  j <- results$outcome == Y & results$analysis == a
  results$FDR[j] <- p.adjust(ifelse(results$status[j] == "ok", results$p[j], NA_real_), method = "BH", n = 182L)
}
i <- match(results$name, meta$name)
cox_table <- data.frame(name = results$name, "Food liking" = meta$Food_liking[i],
  Category = meta$category[i], Disease = results$outcome,
  "Disease label" = unname(out_map[results$outcome]), Analysis = results$analysis,
  HR = results$HR, Conf.low = results$CI_low, Conf.high = results$CI_high,
  "P-value" = results$p, FDR = results$FDR,
  Direction = ifelse(results$status != "ok", NA_character_, ifelse(results$HR > 1, "harmful", ifelse(results$HR < 1, "protective", "neutral"))),
  Beta = results$beta, SE = results$se, N = results$n, Events = results$events,
  "Person-years" = results$person_years, N_low = results$n_low, N_medium = results$n_medium,
  N_high = results$n_high, Status = results$status, Message = results$message, check.names = FALSE)
cox_table <- cox_table[order(i, match(cox_table$Disease, Ys), match(cox_table$Analysis, unique(results$analysis))), ]
write.csv(cox_table, file.path(outdir, "cox_results_table.csv"), row.names = FALSE, na = "", fileEncoding = "UTF-8")
# Snapshot metadata used for analysis so plotting does not depend on subsequent Excel edits.
saveRDS(list(table = cox_table, meta = meta, Ys = Ys, out_map = out_map, covariates = covs.all, config = cfg), file.path(outdir, "cox_results.rds"))
cat("Saved:", file.path(outdir, "cox_results_table.csv"), "\n")

# Demographic characteristics -----------------------------------------------
# 描述 s0_prepare_data.R 保存的未插补队列（至少一个结局符合随访条件）。
# 不是某个疾病/食物偏好模型的完整病例样本；缺失值保留为单独类别。
demo <- readRDS(file.path(project, "00data/dat.base.rds"))
stopifnot(nrow(demo) > 0, !anyDuplicated(demo$eid))
demo_n <- nrow(demo)
demo_cut <- function(x, breaks, labels) {
  x[!is.finite(x)] <- NA_real_
  cut(x, breaks, labels = labels, right = FALSE)
}
# pa_met_tot 仅按 MET-min/week 数值分组，不等同于完整 IPAQ 分类算法。
# age、bmi、tdi、pa_met_tot 在 Cox 中仍为连续变量；以下切点仅用于展示。
demo_vars <- list(
  Age = demo_cut(demo$age, c(-Inf, 60, Inf), c("< 60", ">= 60")),
  Sex = demo$sex,
  BMI = demo_cut(demo$bmi, c(-Inf, 18.5, 25, 30, Inf), c("< 18.5", "18.5 to < 25", "25 to < 30", ">= 30")),
  Education = demo$edu,
  Income = demo$inc,
  TDI = demo_cut(demo$tdi, c(-Inf, -2.5, 0, 2.5, Inf), c("< -2.5", "-2.5 to < 0", "0 to < 2.5", ">= 2.5")),
  `Current smoking` = demo$smoke_status,
  `Current drinking` = demo$alcohol_status,
  `Physical activity` = demo_cut(demo$pa_met_tot, c(0, 600, 3000, Inf), c("< 600", "600 to < 3000", ">= 3000"))
)
demo_counts <- lapply(seq_along(demo_vars), function(i) {
  x <- demo_vars[[i]]
  stopifnot(length(x) == demo_n, is.factor(x))
  lev <- levels(x)
  v <- as.character(x)
  if (anyNA(v)) { lev <- c(lev, "Missing"); v[is.na(v)] <- "Missing" }
  n <- as.integer(table(factor(v, levels = lev)))
  data.frame(characteristic = names(demo_vars)[i], category = lev, n = n, percent = 100 * n / demo_n)
})
demo_table <- do.call(rbind, demo_counts)
stopifnot(all(vapply(demo_counts, function(x) sum(x$n), numeric(1)) == demo_n))
write.csv(demo_table[c("characteristic", "category", "n", "percent")], file.path(outdir, "demographic_characteristics.csv"), row.names = FALSE)
