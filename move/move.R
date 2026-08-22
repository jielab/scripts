
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 🚩 Move v6 pipeline: multi-trait TDI transition, SES, disease risk, and omics
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

get_arg <- function(key, default = NA_character_) {
  args <- commandArgs(trailingOnly = TRUE)
  i <- match(key, args)
  if (!is.na(i) && i < length(args)) args[i + 1] else default
}

step_order <- c("data_prep", "qc", "ses_disease", "omics_signature", "deep_dive", "consolidate")
start_step <- Sys.getenv("START_STEP", unset = get_arg("--start_step", "all"))
if (is.na(start_step) || start_step %in% c("", "all", "ALL")) start_step <- step_order[1]
if (start_step %in% c("help", "--help", "-h")) {
  cat("Step order:\n")
  cat(paste0(seq_along(step_order), ". ", step_order, collapse = "\n"), "\n")
  quit(save = "no")
}
if (!start_step %in% step_order) stop("Unknown START_STEP: ", start_step, "; available: ", paste(step_order, collapse = ", "))
run_step <- function(step) match(step, step_order) >= match(start_step, step_order)
run_if <- function(step, expr) {
  if (run_step(step)) {
    cat("\n========== RUN STEP:", step, "==========\n")
    eval(substitute(expr), envir = parent.frame())
    cat("========== DONE STEP:", step, "==========\n")
  } else {
    cat("Skip step:", step, "\n")
  }
}

Y <- Sys.getenv("Y", unset = get_arg("--trait", "mi"))
main_transition <- Sys.getenv("MOVE_TRANSITION", unset = get_arg("--transition", "birth_current"))
if (!main_transition %in% c("birth_current", "home_future")) stop("MOVE_TRANSITION must be birth_current or home_future.")

if (!requireNamespace("pacman", quietly = TRUE)) {
  install.packages("pacman", repos = "https://cloud.r-project.org")
}
suppressPackageStartupMessages(pacman::p_load(
  data.table, ggplot2, patchwork, writexl, matrixStats, lubridate, scales,
  RANN, survival, broom, ggrepel, stringr, R.utils, forcats
))

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Paths, logging, and global settings
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

dir0 <- Sys.getenv("DIR0", unset = "")
if (!nzchar(dir0)) dir0 <- ifelse(Sys.info()[["sysname"]] == "Windows", "D:", "/mnt/d")
indir <- Sys.getenv("UKB_PHE", unset = file.path(dir0, "data", "ukb", "phe"))
outroot <- Sys.getenv("MOVE_OUTDIR", unset = file.path(dir0, "analysis", "move"))
Y_safe <- gsub("[^A-Za-z0-9_.-]+", "_", Y)
trait_outdir <- Sys.getenv("MOVE_TRAIT_OUTDIR", unset = file.path(outroot, Y_safe))
outdir <- outroot
rawdir <- file.path(outroot, "raw")
dir.create(outroot, recursive = TRUE, showWarnings = FALSE)
dir.create(rawdir, recursive = TRUE, showWarnings = FALSE)
setwd(outroot)

set_output_scope <- function(scope = c("root", "trait")) {
  scope <- match.arg(scope)
  if (scope == "root") {
    outdir <<- outroot
    rawdir <<- file.path(outroot, "raw")
  } else {
    outdir <<- trait_outdir
    rawdir <<- file.path(trait_outdir, "raw")
  }
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  dir.create(rawdir, recursive = TRUE, showWarnings = FALSE)
  setwd(outdir)
  invisible(NULL)
}
with_output_scope <- function(scope = c("root", "trait"), expr) {
  old_outdir <- outdir
  old_rawdir <- rawdir
  old_wd <- getwd()
  on.exit({
    outdir <<- old_outdir
    rawdir <<- old_rawdir
    setwd(old_wd)
  }, add = TRUE)
  set_output_scope(scope)
  eval(substitute(expr), envir = parent.frame())
}
fig_file <- function(n, stem, ext = ".png", trait = FALSE) {
  stem <- gsub("^\\.+|\\.+$", "", stem)
  parts <- c(paste0("Fig", n), if (trait) Y_safe else NULL, stem)
  paste0(paste(parts[nzchar(parts)], collapse = "."), ext)
}
read_core <- function() {
  f_trait <- file.path(trait_outdir, "raw", "move_core.rds")
  f_root_y <- file.path(outroot, "raw", paste0("move_core.", Y_safe, ".rds"))
  f_root  <- file.path(outroot, "raw", "move_core.rds")
  f_root_old <- file.path(outroot, "move_core.rds")
  f <- if (file.exists(f_trait)) f_trait else if (file.exists(f_root_y)) f_root_y else if (file.exists(f_root)) f_root else f_root_old
  if (!file.exists(f)) stop("Cannot find move_core.rds. Run START_STEP=data_prep first. Tried: ", paste(c(f_trait, f_root_y, f_root, f_root_old), collapse = "; "), call. = FALSE)
  readRDS(f)
}
read_home_long <- function() {
  f_root <- file.path(outroot, "raw", "move_home_long.rds")
  f_old <- file.path(outroot, "move_home_long.rds")
  f <- if (file.exists(f_root)) f_root else f_old
  if (!file.exists(f)) stop("Cannot find move_home_long.rds. Run START_STEP=data_prep first.", call. = FALSE)
  readRDS(f)
}

log_file <- file.path(outroot, "raw", paste0("run.", Y_safe, ".log"))
log_con <- file(log_file, open = if (file.exists(log_file) && start_step != step_order[1]) "at" else "wt")
sink(log_con, split = TRUE)
sink(log_con, type = "message")
on.exit({
  while (sink.number(type = "message") > 0) sink(type = "message")
  while (sink.number() > 0) sink()
  try(close(log_con), silent = TRUE)
}, add = TRUE)

cat("Move v6 pipeline\n")
cat("Time:", as.character(Sys.time()), "\n")
cat("Trait Y:", Y, "\n")
cat("Main transition:", main_transition, "\n")
cat("dir0:", dir0, "\n")
cat("indir:", indir, "\n")
cat("outroot:", outroot, "\n")
cat("trait_outdir:", trait_outdir, "\n")
cat("START_STEP:", start_step, "\n")

set.seed(as.integer(Sys.getenv("SEED", unset = "2026")))
move_date_follow_end <- as.Date(Sys.getenv("MOVE_DATE_FOLLOW_END", unset = "2023-04-01"))
rich_poor_q <- as.numeric(Sys.getenv("MOVE_RICH_POOR_Q", unset = "0.25"))
main_nn_km <- as.numeric(Sys.getenv("MOVE_TDI_NN_KM", unset = "5"))
min_n_model <- as.integer(Sys.getenv("MOVE_MIN_N", unset = "500"))
min_case_model <- as.integer(Sys.getenv("MOVE_MIN_CASE", unset = "30"))
max_prot_features <- as.integer(Sys.getenv("MOVE_MAX_PROT_FEATURES", unset = "0"))
max_met_features <- as.integer(Sys.getenv("MOVE_MAX_MET_FEATURES", unset = "0"))
prot_do <- toupper(Sys.getenv("PROT_DO", unset = "TRUE")) %in% c("TRUE", "T", "1", "YES", "Y")
met_do <- toupper(Sys.getenv("MET_DO", unset = "TRUE")) %in% c("TRUE", "T", "1", "YES", "Y")
resume_omics <- toupper(Sys.getenv("MOVE_RESUME_OMICS", unset = "FALSE")) %in% c("TRUE", "T", "1", "YES", "Y")
outcome_scan_env <- Sys.getenv("MOVE_OUTCOME_SCAN", unset = Sys.getenv("MOVE_CAD_SCAN", unset = "TRUE"))
cad_scan_do <- toupper(outcome_scan_env) %in% c("TRUE", "T", "1", "YES", "Y")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Helper functions
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

`%||%` <- function(a, b) if (!is.null(a)) a else b
bt <- function(x) paste0("`", gsub("`", "", x), "`")

asDate2 <- function(x) {
  if (inherits(x, "Date")) return(x)
  if (inherits(x, "IDate")) return(as.Date(x))
  if (inherits(x, "POSIXct") || inherits(x, "POSIXt")) return(as.Date(x))
  if (is.character(x)) return(suppressWarnings(as.Date(x)))
  y <- suppressWarnings(as.numeric(unclass(x)))
  y[!is.finite(y)] <- NA_real_
  as.Date(y, origin = "1970-01-01")
}
safe_num <- function(x) suppressWarnings(as.numeric(as.character(x)))
safe_mean <- function(x) { x <- safe_num(x); x <- x[is.finite(x)]; if (length(x)) mean(x) else NA_real_ }
safe_median <- function(x) { x <- safe_num(x); x <- x[is.finite(x)]; if (length(x)) median(x) else NA_real_ }
safe_sd <- function(x) { x <- safe_num(x); x <- x[is.finite(x)]; if (length(x) > 1) sd(x) else NA_real_ }
min_date2 <- function(x) { x <- asDate2(x); x <- x[!is.na(x)]; if (length(x)) min(x) else as.Date(NA) }
max_date2 <- function(x) { x <- asDate2(x); x <- x[!is.na(x)]; if (length(x)) max(x) else as.Date(NA) }
first_non_na <- function(x) { x <- x[!is.na(x)]; if (length(x)) x[1] else NA }
last_non_na <- function(x) { x <- x[!is.na(x)]; if (length(x)) x[length(x)] else NA }
zscore <- function(x) { x <- safe_num(x); s <- sd(x, na.rm = TRUE); m <- mean(x, na.rm = TRUE); if (!is.finite(s) || s == 0) return(rep(NA_real_, length(x))); (x - m) / s }
inormal2 <- function(x) {
  x <- safe_num(x)
  ok <- is.finite(x)
  out <- rep(NA_real_, length(x))
  if (sum(ok) >= 3) out[ok] <- qnorm((rank(x[ok], ties.method = "average") - 0.5) / sum(ok))
  out
}
label_n_pct <- function(n, total) paste0(format(n, big.mark = ","), " (", scales::percent(n / total, accuracy = 0.1), ")")
fig_theme <- function(base_size = 10) {
  theme_classic(base_size = base_size) +
    theme(plot.title = element_text(face = "bold", hjust = 0.5),
          axis.title = element_text(face = "bold"),
          legend.title = element_text(face = "bold"))
}

# Paper-style plotting helpers, adapted from the 5C scripts.
BLUE <- "#2C7FB8"
GOLD <- "#D8A03D"
TEAL <- "#00A6B2"
ORANGE <- "#E07A3F"
PURPLE <- "#7B61A8"
GREEN <- "#4D9221"
DARK <- "grey25"
LIGHT <- "grey88"
cols_transition <- c(
  "poor -> poor" = "#7F8C8D",
  "poor -> rich" = GOLD,
  "rich -> poor" = TEAL,
  "rich -> rich" = BLUE
)
cols_ses <- c("Low SES" = ORANGE, "Middle SES" = "#BDBDBD", "High SES" = BLUE)
cols_candidate <- c(
  high_priority_resilience_driver = GOLD,
  medium_priority_driver = TEAL,
  exploratory = "grey70",
  stress_or_consequence = "#B35806"
)
cols_signature <- c(
  `FDR<0.05, positive` = GOLD,
  `FDR<0.05, inverse` = BLUE,
  `P<0.05 only` = TEAL,
  `Not significant` = "grey76"
)
theme_5c <- function(base_size = 12) {
  theme_bw(base_size = base_size) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0, size = base_size + 1),
      plot.subtitle = element_text(color = "grey30", size = base_size - 1),
      axis.title = element_text(face = "bold"),
      axis.text = element_text(color = "grey15"),
      strip.background = element_rect(fill = "grey92", color = NA),
      strip.text = element_text(face = "bold"),
      panel.grid.minor = element_blank(),
      legend.position = "bottom",
      legend.title = element_text(face = "bold")
    )
}
fig_theme <- function(base_size = 10) theme_5c(base_size)
short_transition <- function(x) {
  x <- as.character(x)
  dplyr::recode(x,
    "poor -> poor" = "Poor → poor",
    "poor -> rich" = "Poor → rich",
    "rich -> poor" = "Rich → poor",
    "rich -> rich" = "Rich → rich",
    .default = x
  )
}
clean_hypothesis_label <- function(x) {
  x <- as.character(x)
  dplyr::recode(x,
    poor_to_rich_lowSES = "Poor → rich + low SES",
    poor_to_rich_highSES = "Poor → rich + high SES",
    rich_to_rich_lowSES = "Rich → rich + low SES",
    rich_to_rich_highSES = "Rich → rich + high SES",
    poor_to_poor_lowSES = "Poor → poor + low SES",
    poor_to_poor_highSES = "Poor → poor + high SES",
    rich_to_poor_lowSES = "Rich → poor + low SES",
    rich_to_poor_highSES = "Rich → poor + high SES",
    .default = x
  )
}
clean_cox_label <- function(term) {
  out <- as.character(term)
  out <- gsub("^tdi_transition", "", out)
  out <- gsub("^hypothesis_group", "", out)
  out <- gsub("^ses_extreme", "", out)
  out <- gsub("low_sesTRUE", "Low SES", out)
  out <- gsub(":", " × ", out)
  out <- gsub("poor -> poor", "Poor → poor", out, fixed = TRUE)
  out <- gsub("poor -> rich", "Poor → rich", out, fixed = TRUE)
  out <- gsub("rich -> poor", "Rich → poor", out, fixed = TRUE)
  out <- gsub("rich -> rich", "Rich → rich", out, fixed = TRUE)
  out <- clean_hypothesis_label(out)
  out <- gsub("_", " ", out)
  trimws(out)
}

wrap_vs_label <- function(x) {
  x <- as.character(x)
  x <- gsub(" vs ", "\nvs\n", x, fixed = TRUE)
  x <- gsub(" \\+ ", "\n+ ", x)
  x
}
wrap_group_label <- function(x, width = 24) {
  x <- as.character(x)
  x <- gsub(" \\+ ", "\n+ ", x)
  stringr::str_wrap(x, width = width)
}
score_signed_logp <- function(beta, p, positive_good = TRUE, cap = 8) {
  lp <- -log10(pmax(as.numeric(p), 1e-300))
  lp[!is.finite(lp)] <- 0
  lp <- pmin(lp, cap)
  sgn <- sign(as.numeric(beta))
  sgn[!is.finite(sgn)] <- 0
  out <- sgn * lp
  if (!positive_good) out <- -out
  out
}
score_pos <- function(x) {
  x <- as.numeric(x)
  x[!is.finite(x)] <- 0
  pmax(x, 0)
}
score_abs <- function(x) {
  x <- as.numeric(x)
  x[!is.finite(x)] <- 0
  abs(x)
}

plot_forest_clean <- function(d, title = NULL, ref_label = NULL, base_size = 11, xlim = NULL) {
  d <- as.data.table(d)
  if (!nrow(d) || !all(c("HR", "lo", "hi", "label") %in% names(d))) {
    return(ggplot() + labs(title = title %||% "Forest plot unavailable") + theme_5c(base_size))
  }
  d <- d[is.finite(HR) & is.finite(lo) & is.finite(hi) & HR > 0 & lo > 0 & hi > 0]
  if (!nrow(d)) return(ggplot() + labs(title = title %||% "Forest plot unavailable") + theme_5c(base_size))
  d[, sig := is.finite(p) & p < 0.05]
  d[, label := factor(label, levels = rev(unique(label)))]
  p <- ggplot(d, aes(HR, label)) +
    geom_vline(xintercept = 1, linetype = "dashed", color = "grey55") +
    geom_errorbar(aes(xmin = lo, xmax = hi, color = sig), orientation = "y", width = .18, linewidth = .85) +
    geom_point(aes(fill = sig), shape = 21, size = 2.8, color = "grey20") +
    scale_x_log10() +
    scale_color_manual(values = c(`TRUE` = GOLD, `FALSE` = "grey45"), guide = "none") +
    scale_fill_manual(values = c(`TRUE` = GOLD, `FALSE` = "white"), guide = "none") +
    labs(title = title, subtitle = ref_label, x = "Hazard ratio", y = NULL) +
    theme_5c(base_size)
  if (!is.null(xlim)) p <- p + coord_cartesian(xlim = xlim)
  p
}
plot_adjusted_10y_risk <- function(dat, covs, t0 = 10, title = NULL) {
  d <- as.data.table(dat)
  d <- d[!is.na(y_t2e) & !is.na(y_event) & y_t2e > 0 & y_prevalent != 1 & !is.na(tdi_transition) & !is.na(ses_extreme)]
  covs <- covs[covs %in% names(d)]
  covs <- covs[vapply(covs, function(v) length(unique(na.omit(d[[v]]))) > 1, logical(1))]
  if (nrow(d) < min_n_model || sum(d$y_event == 1, na.rm = TRUE) < min_case_model) {
    return(list(plot = ggplot() + labs(title = title %||% "Adjusted risk unavailable") + theme_5c(), table = data.table()))
  }
  d[, tdi_transition := factor(as.character(tdi_transition), levels = names(cols_transition))]
  d[, ses_extreme := factor(as.character(ses_extreme), levels = c("Low SES", "High SES"))]
  f <- as.formula(paste0("survival::Surv(y_t2e, y_event) ~ tdi_transition * ses_extreme",
                         if (length(covs)) paste0(" + ", paste(bt(covs), collapse = " + ")) else ""))
  fit <- tryCatch(survival::coxph(f, data = d), error = function(e) NULL)
  if (is.null(fit)) return(list(plot = ggplot() + labs(title = title %||% "Adjusted risk unavailable") + theme_5c(), table = data.table()))
  nd <- expand.grid(tdi_transition = factor(names(cols_transition), levels = names(cols_transition)),
                    ses_extreme = factor(c("Low SES", "High SES"), levels = c("Low SES", "High SES")))
  for (v in covs) nd[[v]] <- if (is.numeric(d[[v]])) mean(d[[v]], na.rm = TRUE) else names(which.max(table(d[[v]])))
  ss <- summary(survival::survfit(fit, newdata = nd), times = t0, extend = TRUE)
  nd$risk <- 1 - as.numeric(ss$surv)
  nd$lo <- 1 - as.numeric(ss$upper)
  nd$hi <- 1 - as.numeric(ss$lower)
  nd$label <- sprintf("%.1f%%", 100 * nd$risk)
  nd$transition_label <- factor(short_transition(nd$tdi_transition), levels = short_transition(names(cols_transition)))
  p <- ggplot(nd, aes(transition_label, risk, fill = ses_extreme)) +
    geom_col(position = position_dodge(.76), width = .68, color = "white") +
    geom_errorbar(aes(ymin = lo, ymax = hi), width = .16, position = position_dodge(.76), linewidth = .7) +
    geom_text(aes(label = label, y = hi), position = position_dodge(.76), vjust = -.45, size = 3.2, fontface = "bold") +
    scale_fill_manual(values = cols_ses[c("Low SES", "High SES")]) +
    scale_y_continuous(labels = scales::label_percent(accuracy = 1), expand = expansion(mult = c(0, .22))) +
    labs(title = title, subtitle = paste0("Cox-adjusted ", t0, "-year risk; prevalent cases excluded"), x = "TDI transition", y = paste0(t0, "-year risk"), fill = NULL) +
    theme_5c(11)

  list(plot = p, table = as.data.table(nd))
}
plot_crude_origin_tile <- function(dat, ses_level, title = NULL) {
  d <- as.data.table(dat)
  d <- d[!is.na(source_area) & !is.na(dest_area) & !is.na(ses_extreme) & ses_extreme == ses_level & y_prevalent != 1]
  if (!nrow(d)) return(ggplot() + labs(title = title %||% paste("No data", ses_level)) + theme_5c())
  tab <- d[, .(N = .N, events = sum(y_event == 1, na.rm = TRUE), risk = mean(y_event == 1, na.rm = TRUE)), by = .(source_area, dest_area)]
  tab[, source_area := factor(source_area, levels = c("poor", "rich"), labels = c("Source poor", "Source rich"))]
  tab[, dest_area := factor(dest_area, levels = c("poor", "rich"), labels = c("Current poor", "Current rich"))]
  tab[, label := paste0(percent(risk, accuracy = .1), "\n", events, " events\nN=", comma(N))]
  ggplot(tab, aes(dest_area, source_area, fill = risk)) +
    geom_tile(color = "white", linewidth = .8) +
    geom_text(aes(label = label), fontface = "bold", size = 3.4) +
    scale_fill_gradient(low = "#F7FBFF", high = ORANGE, labels = percent, name = "Incident\nrisk") +
    labs(title = title, x = "Current community", y = "Origin community") + theme_5c(11) +
    theme(panel.grid = element_blank(), axis.text = element_text(face = "bold"))
}

plot_crude_origin_tile_combined <- function(dat, title = NULL) {
  d <- as.data.table(dat)
  d <- d[!is.na(source_area) & !is.na(dest_area) & !is.na(ses_extreme) & y_prevalent != 1]
  if (!nrow(d)) return(ggplot() + labs(title = title %||% "Crude incidence unavailable") + theme_5c())
  tab <- d[, .(N = .N, events = sum(y_event == 1, na.rm = TRUE), risk = mean(y_event == 1, na.rm = TRUE)),
           by = .(ses_extreme, source_area, dest_area)]
  tab[, source_area := factor(source_area, levels = c("poor", "rich"), labels = c("Source poor", "Source rich"))]
  tab[, dest_area := factor(dest_area, levels = c("poor", "rich"), labels = c("Current poor", "Current rich"))]
  tab[, ses_extreme := factor(as.character(ses_extreme), levels = c("Low SES", "High SES"))]
  tab[, label := paste0(percent(risk, accuracy = .1), "\n", events, " events\nN=", comma(N))]
  ggplot(tab, aes(dest_area, source_area, fill = risk)) +
    geom_tile(color = "white", linewidth = .8) +
    geom_text(aes(label = label), fontface = "bold", size = 3.2) +
    facet_wrap(~ses_extreme, nrow = 1) +
    scale_fill_gradient(low = "#F7FBFF", high = ORANGE, labels = percent, name = "Incident\nrisk") +
    labs(title = title, x = "Current community", y = "Origin community") + theme_5c(10) +
    theme(panel.grid = element_blank(), axis.text = element_text(face = "bold"))
}
plot_od_forest_combined <- function(d, title = NULL, base_size = 10) {
  dd <- as.data.table(d)
  if (!nrow(dd) || !all(c("HR", "lo", "hi", "label", "contrast_family") %in% names(dd))) {
    return(ggplot() + labs(title = title %||% "Origin/destination effects unavailable") + theme_5c(base_size))
  }
  dd <- dd[is.finite(HR) & is.finite(lo) & is.finite(hi) & HR > 0 & lo > 0 & hi > 0]
  if (!nrow(dd)) return(ggplot() + labs(title = title %||% "Origin/destination effects unavailable") + theme_5c(base_size))
  dd[, sig := is.finite(p) & p < 0.05]
  dd[, family_lab := fifelse(grepl("^Origin", contrast_family),
                             "Origin effect\nreference = rich origin",
                             "Destination effect\nreference = current poor")]
  dd[, label_clean := gsub("\\n", " ", label)]
  dd[, label_clean := stringr::str_wrap(label_clean, width = 35)]
  dd[, label_clean := factor(label_clean, levels = rev(unique(label_clean)))]
  ggplot(dd, aes(HR, label_clean)) +
    geom_vline(xintercept = 1, linetype = "dashed", color = "grey55") +
    geom_errorbar(aes(xmin = lo, xmax = hi, color = sig), orientation = "y", width = .18, linewidth = .80) +
    geom_point(aes(fill = sig), shape = 21, size = 2.6, color = "grey20") +
    facet_grid(family_lab ~ ., scales = "free_y", space = "free_y") +
    scale_x_log10() +
    scale_color_manual(values = c(`TRUE` = GOLD, `FALSE` = "grey45"), guide = "none") +
    scale_fill_manual(values = c(`TRUE` = GOLD, `FALSE` = "white"), guide = "none") +
    labs(title = title, x = "Hazard ratio", y = NULL) + theme_5c(base_size) +
    theme(axis.text.y = element_text(face = "bold", size = base_size - 1), strip.text.y = element_text(face = "bold", angle = 0))
}

cox_binary_contrast <- function(dat, subset_expr, var, case_level, ctrl_level, label, covs_extra = c("source_tdi", "dest_tdi")) {
  d <- as.data.table(dat)
  d <- d[eval(parse(text = subset_expr))]
  d <- d[!is.na(y_t2e) & !is.na(y_event) & y_t2e > 0 & y_prevalent != 1]
  if (!nrow(d) || !var %in% names(d)) return(data.table())
  d <- d[get(var) %in% c(case_level, ctrl_level)]
  d[, contrast01 := fifelse(get(var) == case_level, 1L, fifelse(get(var) == ctrl_level, 0L, NA_integer_))]
  d <- d[!is.na(contrast01)]
  if (nrow(d) < min_n_model || sum(d$y_event == 1, na.rm = TRUE) < min_case_model || length(unique(d$contrast01)) < 2) return(data.table())
  covs <- model_covs(d, extra = covs_extra)
  covs <- setdiff(covs, c("source_area", "dest_area", "tdi_transition", "individual_ses_z"))
  f <- as.formula(paste0("survival::Surv(y_t2e, y_event) ~ contrast01", if (length(covs)) paste0(" + ", paste(bt(covs), collapse = " + ")) else ""))
  fit <- tryCatch(survival::coxph(f, data = d), error = function(e) NULL)
  if (is.null(fit)) return(data.table())
  out <- tidy_cox_term(fit)[term == "contrast01"]
  if (!nrow(out)) return(data.table())
  out[, `:=`(label = label, N = nrow(d), events = sum(d$y_event == 1, na.rm = TRUE), case_level = case_level, ctrl_level = ctrl_level)]
  out
}
make_origin_destination_contrasts <- function(dat) {
  specs1 <- rbindlist(lapply(c("poor", "rich"), function(dst) {
    rbindlist(lapply(c("Low SES", "High SES"), function(ses) {
      cox_binary_contrast(dat,
        subset_expr = sprintf('dest_area == "%s" & ses_extreme == "%s"', dst, ses),
        var = "source_area", case_level = "poor", ctrl_level = "rich",
        label = paste0("Poor origin vs rich origin\ncurrent ", dst, " + ", ses))
    }), fill = TRUE)
  }), fill = TRUE)
  if (nrow(specs1)) specs1[, contrast_family := "Origin effect within same current community and SES"]
  specs2 <- rbindlist(lapply(c("poor", "rich"), function(src) {
    rbindlist(lapply(c("Low SES", "High SES"), function(ses) {
      cox_binary_contrast(dat,
        subset_expr = sprintf('source_area == "%s" & ses_extreme == "%s"', src, ses),
        var = "dest_area", case_level = "rich", ctrl_level = "poor",
        label = paste0("Current rich vs current poor\norigin ", src, " + ", ses))
    }), fill = TRUE)
  }), fill = TRUE)
  if (nrow(specs2)) specs2[, contrast_family := "Destination effect within same origin and SES"]
  rbindlist(list(specs1, specs2), fill = TRUE)
}
plot_origin_gap <- function(dat) {
  d <- as.data.table(dat)
  d <- d[!is.na(dest_area) & !is.na(source_area) & !is.na(ses_extreme) & y_prevalent != 1]
  tab <- d[, .(N = .N, events = sum(y_event == 1, na.rm = TRUE), risk = mean(y_event == 1, na.rm = TRUE)), by = .(dest_area, ses_extreme, source_area)]
  w <- dcast(tab, dest_area + ses_extreme ~ source_area, value.var = "risk")
  if (!all(c("poor", "rich") %in% names(w))) return(ggplot() + labs(title = "f. Origin risk gap unavailable") + theme_5c())
  w[, gap := poor - rich]
  w[, label := paste0("Current ", dest_area, "\n", ses_extreme)]
  w[, label := factor(label, levels = label[order(gap)])]
  ggplot(w, aes(gap, label, fill = gap > 0)) +
    geom_vline(xintercept = 0, color = "grey55", linetype = "dashed") +
    geom_col(width = .68, color = "white") +
    geom_text(aes(label = paste0(ifelse(gap > 0, "+", ""), percent(gap, accuracy = .1))), hjust = ifelse(w$gap > 0, -0.15, 1.15), fontface = "bold") +
    scale_fill_manual(values = c(`TRUE` = ORANGE, `FALSE` = BLUE), guide = "none") +
    scale_x_continuous(labels = percent, expand = expansion(mult = c(.18, .25))) +
    labs(title = "f. Excess crude risk from poor origin", x = "Poor-origin risk minus rich-origin risk", y = NULL) + theme_5c(11) +
    theme(axis.text.y = element_text(face = "bold"))
}
save_plot <- function(p, fn, width = 12, height = 8, dpi = 300) {
  ggsave(file.path(outdir, fn), p, width = width, height = height, dpi = dpi, bg = "white", limitsize = FALSE)
  cat("Wrote:", file.path(outdir, fn), "\n")
}
write_book <- function(x, fn) {
  x <- x[vapply(x, function(z) is.data.frame(z) && nrow(z) > 0, logical(1))]
  if (!length(x)) x <- list(empty = data.frame(note = "No rows"))
  names(x) <- substr(make.unique(gsub("[^A-Za-z0-9_]+", "_", names(x))), 1, 31)
  writexl::write_xlsx(x, file.path(outdir, fn))
  cat("Wrote:", file.path(outdir, fn), "\n")
}
write_raw_book <- function(x, fn) {
  x <- x[vapply(x, function(z) is.data.frame(z) && nrow(z) > 0, logical(1))]
  if (!length(x)) x <- list(empty = data.frame(note = "No rows"))
  names(x) <- substr(make.unique(gsub("[^A-Za-z0-9_]+", "_", names(x))), 1, 31)
  f <- file.path(rawdir, fn)
  dir.create(dirname(f), recursive = TRUE, showWarnings = FALSE)
  writexl::write_xlsx(x, f)
  cat("Wrote:", f, "\n")
}
write_raw <- function(x, fn) {
  f <- file.path(rawdir, fn)
  data.table::fwrite(as.data.table(x), f)
  cat("Wrote:", f, "\n")
}
save_rds_out <- function(obj, fn, compress = "xz") {
  f1 <- file.path(rawdir, fn)
  dir.create(dirname(f1), recursive = TRUE, showWarnings = FALSE)
  saveRDS(obj, f1, compress = compress)
  cat("Wrote:", f1, "\n")
}
sync_root_rds_to_raw <- function() invisible(NULL)
pick_col <- function(DT, candidates) {
  hit <- candidates[candidates %in% names(DT)][1]
  if (length(hit) && !is.na(hit)) hit else NA_character_
}
slot_id <- function(x, base) as.integer(sub(paste0("^", base, "_a"), "", x))
mk_slot_mat <- function(DT, prefix, idx, min_val = -Inf, max_val = Inf) {
  if (!length(idx)) return(matrix(NA_real_, nrow(DT), 0))
  m <- matrix(NA_real_, nrow(DT), length(idx)); colnames(m) <- paste0(prefix, "_a", idx)
  for (j in seq_along(idx)) {
    cc <- paste0(prefix, "_a", idx[j])
    if (cc %in% names(DT)) {
      x <- safe_num(DT[[cc]]); x[x < min_val | x > max_val] <- NA_real_; m[, j] <- x
    }
  }
  m
}
row_min_na <- function(m) { if (!ncol(m)) return(rep(NA_real_, nrow(m))); x <- matrixStats::rowMins(m, na.rm = TRUE); x[!is.finite(x)] <- NA_real_; x }
row_max_na <- function(m) { if (!ncol(m)) return(rep(NA_real_, nrow(m))); x <- matrixStats::rowMaxs(m, na.rm = TRUE); x[!is.finite(x)] <- NA_real_; x }
row_adj_change <- function(m) {
  if (ncol(m) < 2) return(rep(0L, nrow(m)))
  out <- matrix(FALSE, nrow(m), ncol(m) - 1)
  for (j in 2:ncol(m)) out[, j - 1] <- !is.na(m[, j - 1]) & !is.na(m[, j]) & m[, j - 1] != m[, j]
  rowSums(out, na.rm = TRUE)
}
min_date_cols <- function(DT, cols) {
  cols <- intersect(cols, names(DT))
  if (!length(cols)) return(as.Date(rep(NA, nrow(DT))))
  m <- do.call(cbind, lapply(cols, function(cc) as.numeric(asDate2(DT[[cc]]))))
  x <- matrixStats::rowMins(m, na.rm = TRUE); x[!is.finite(x)] <- NA_real_
  as.Date(x, origin = "1970-01-01")
}
km_between <- function(e1, n1, e2, n2) {
  out <- sqrt((safe_num(e2) - safe_num(e1))^2 + (safe_num(n2) - safe_num(n1))^2) / 1000
  out[!is.finite(out)] <- NA_real_
  out
}

find_input <- function(name) {
  f <- file.path(indir, "Rdata", name)
  if (!file.exists(f)) stop("Missing input file: ", f, call. = FALSE)
  f
}

read_main_data <- function() {
  cat("\n# Read UKB data\n")
  dat0 <- readRDS(find_input("all.rds"))
  DT <- as.data.table(dat0)
  if (!"eid" %in% names(DT)) stop("all.rds must contain eid.", call. = FALSE)
  DT[, eid := as.character(eid)]
  if ("date_attend" %in% names(DT)) DT[, date_attend := asDate2(date_attend)] else stop("all.rds must contain date_attend.", call. = FALSE)
  DT[, attend_year := as.integer(format(date_attend, "%Y"))]
  cat("nrow(all.rds)=", nrow(DT), " ncol=", ncol(DT), "\n", sep = "")
  DT
}

read_layer <- function(file, required = FALSE) {
  f <- file.path(indir, "Rdata", file)
  if (!file.exists(f)) {
    if (required) stop("Missing input file: ", f, call. = FALSE)
    return(NULL)
  }
  x <- readRDS(f)
  x <- as.data.table(x)
  if (!"eid" %in% names(x)) stop(file, " must contain eid.", call. = FALSE)
  x[, eid := as.character(eid)]
  x
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Derivation functions
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

make_area_class <- function(tdi, q25, q75) {
  out <- rep(NA_character_, length(tdi))
  out[is.finite(tdi) & tdi <= q25] <- "rich"
  out[is.finite(tdi) & tdi >= q75] <- "poor"
  out[is.finite(tdi) & tdi > q25 & tdi < q75] <- "middle"
  out
}
make_transition <- function(a, b) {
  out <- ifelse(is.na(a) | is.na(b), NA_character_, paste(a, "->", b))
  factor(out, levels = c("poor -> poor", "poor -> rich", "rich -> poor", "rich -> rich",
                         "middle -> rich", "middle -> poor", "rich -> middle", "poor -> middle",
                         "middle -> middle"))
}
make_transition_rp <- function(a, b) {
  out <- ifelse(a %in% c("poor", "rich") & b %in% c("poor", "rich"), paste(a, "->", b), NA_character_)
  factor(out, levels = c("poor -> poor", "poor -> rich", "rich -> poor", "rich -> rich"))
}

map_tdi_to_coords <- function(query, ref, max_km = 5, prefix = "tdi") {
  query <- as.data.table(query)
  ref <- as.data.table(ref)
  out <- copy(query)
  map_col <- paste0(prefix, "_mapped")
  dist_col <- paste0(prefix, "_nn_km")
  out[, (map_col) := NA_real_]
  out[, (dist_col) := NA_real_]
  okq <- is.finite(safe_num(query$east)) & is.finite(safe_num(query$north))
  okr <- is.finite(safe_num(ref$east)) & is.finite(safe_num(ref$north)) & is.finite(safe_num(ref$tdi))
  if (!any(okq) || !any(okr)) return(out)
  r0 <- ref[okr, .(east = safe_num(east), north = safe_num(north), tdi = safe_num(tdi))]
  # Collapse duplicated 1-km grid coordinates.
  r0 <- r0[, .(tdi = median(tdi, na.rm = TRUE)), by = .(east, north)]
  q0 <- query[okq, .(east = safe_num(east), north = safe_num(north))]
  nn <- RANN::nn2(data = as.matrix(r0[, .(east, north)]), query = as.matrix(q0), k = 1)
  val <- r0$tdi[nn$nn.idx[, 1]]
  dst <- nn$nn.dists[, 1] / 1000
  val[!is.finite(dst) | dst > max_km] <- NA_real_
  idx <- which(okq)
  data.table::set(out, i = idx, j = map_col, value = val)
  data.table::set(out, i = idx, j = dist_col, value = dst)
  out
}

make_home_long <- function(DT) {
  home_date_cols <- grep("^home_date_hist_a[0-9]+$", names(DT), value = TRUE)
  home_east_cols <- grep("^home_east_hist_a[0-9]+$", names(DT), value = TRUE)
  home_north_cols <- grep("^home_north_hist_a[0-9]+$", names(DT), value = TRUE)
  home_idx <- sort(unique(c(slot_id(home_date_cols, "home_date_hist"),
                            slot_id(home_east_cols, "home_east_hist"),
                            slot_id(home_north_cols, "home_north_hist"))))
  home_idx <- home_idx[!is.na(home_idx)]
  if (!length(home_idx)) return(data.table(eid = character()))
  home_long <- rbindlist(lapply(home_idx, function(i) data.table(
    eid = DT$eid,
    slot = as.integer(i),
    date_attend = DT$date_attend,
    home_date = if (paste0("home_date_hist_a", i) %in% names(DT)) asDate2(DT[[paste0("home_date_hist_a", i)]]) else as.Date(NA),
    hist_east = if (paste0("home_east_hist_a", i) %in% names(DT)) safe_num(DT[[paste0("home_east_hist_a", i)]]) else NA_real_,
    hist_north = if (paste0("home_north_hist_a", i) %in% names(DT)) safe_num(DT[[paste0("home_north_hist_a", i)]]) else NA_real_
  )), use.names = TRUE, fill = TRUE)
  home_long <- home_long[!is.na(home_date) | !is.na(hist_east) | !is.na(hist_north)]
  if (!nrow(home_long)) return(home_long)
  home_long[, coord_ok := is.finite(hist_east) & is.finite(hist_north)]
  home_long[, rel_year := as.numeric(home_date - date_attend) / 365.25]
  home_long[, home_vs_baseline := fcase(
    is.na(home_date) | is.na(date_attend), "Missing date",
    home_date < date_attend, "Before baseline",
    home_date == date_attend, "Same day",
    home_date > date_attend, "After baseline"
  )]
  setorder(home_long, eid, home_date, slot)
  home_long[, prev_east := shift(hist_east), by = eid]
  home_long[, prev_north := shift(hist_north), by = eid]
  home_long[, prev_ok := shift(coord_ok), by = eid]
  home_long[, coord_changed_from_prev := coord_ok & prev_ok & (hist_east != prev_east | hist_north != prev_north)]
  home_long[, hist_step_km := fifelse(coord_changed_from_prev, km_between(prev_east, prev_north, hist_east, hist_north), NA_real_)]
  home_long
}

make_ses <- function(DT) {
  out <- data.table(eid = DT$eid)
  camsis_col <- pick_col(DT, c("camsis", "job_camsis", "CAMSIS"))
  inc_col <- pick_col(DT, c("inc.c", "inc", "income", "household_income"))
  edu_col <- pick_col(DT, c("edu.c", "edu", "education"))
  age_edu_col <- pick_col(DT, c("age_edu"))

  out[, camsis := if (!is.na(camsis_col)) safe_num(DT[[camsis_col]]) else NA_real_]
  ordinal_from_any <- function(x) {
    if (is.null(x)) return(rep(NA_real_, nrow(DT)))
    y <- safe_num(x)
    if (sum(is.finite(y)) > 0.5 * length(y)) return(y)
    z <- as.character(x)
    z[z %in% c("", "NA", "Prefer not to answer", "Do not know")] <- NA_character_
    suppressWarnings(as.numeric(factor(z, levels = sort(unique(z[!is.na(z)])), ordered = TRUE)))
  }
  out[, income_ord := if (!is.na(inc_col)) ordinal_from_any(DT[[inc_col]]) else NA_real_]
  out[, edu_ord := if (!is.na(edu_col)) ordinal_from_any(DT[[edu_col]]) else NA_real_]
  out[, age_edu := if (!is.na(age_edu_col)) safe_num(DT[[age_edu_col]]) else NA_real_]
  # Higher age at completing education generally indicates longer education; keep positive direction.
  zmat <- cbind(zscore(out$camsis), zscore(out$income_ord), zscore(out$edu_ord), zscore(out$age_edu))
  out[, ses_n_component := rowSums(is.finite(zmat))]
  out[, individual_ses_z := rowMeans(zmat, na.rm = TRUE)]
  out[ses_n_component == 0, individual_ses_z := NA_real_]
  q <- quantile(out$individual_ses_z, c(.25, .75), na.rm = TRUE)
  out[, person_ses_class := fifelse(is.na(individual_ses_z), NA_character_,
                                    fifelse(individual_ses_z <= q[1], "low",
                                            fifelse(individual_ses_z >= q[2], "high", "middle")))]
  out[, person_ses_class := factor(person_ses_class, levels = c("low", "middle", "high"))]
  attr(out, "ses_sources") <- data.frame(
    component = c("camsis", "income_ord", "edu_ord", "age_edu"),
    source_column = c(camsis_col, inc_col, edu_col, age_edu_col),
    stringsAsFactors = FALSE
  )
  out
}

make_y_outcome <- function(DT, Y) {
  out <- data.table(eid = DT$eid, date_attend = asDate2(DT$date_attend))
  tvar0 <- paste0(Y, ".t2e")
  evar0 <- paste0(Y, ".Yt2e")
  rvar0 <- paste0(Y, ".Yr2e")
  if (all(c(tvar0, evar0) %in% names(DT))) {
    out[, y_t2e := safe_num(DT[[tvar0]])]
    out[, y_event := as.integer(safe_num(DT[[evar0]]) == 1)]
    out[, y_prevalent := if (rvar0 %in% names(DT)) as.integer(safe_num(DT[[rvar0]]) == 1) else NA_integer_]
    out[, y_date := as.Date(NA)]
    out[, y_source := paste(tvar0, evar0)]
    return(out)
  }
  y_patterns <- c(
    paste0("^fod_.*_", Y, "$"),
    paste0("^date_", Y, "$"),
    paste0("^", Y, "_date$"),
    paste0("^", Y, "\\.date$")
  )
  date_cols <- unique(unlist(lapply(y_patterns, function(p) grep(p, names(DT), value = TRUE))))
  date_cols <- date_cols[!grepl("death$|lost$", date_cols)]
  y_date <- min_date_cols(DT, date_cols)

  date_lost_col <- pick_col(DT, c("date_lost"))
  date_death_col <- pick_col(DT, c("date_death"))
  date_lost <- if (!is.na(date_lost_col)) asDate2(DT[[date_lost_col]]) else as.Date(NA)
  date_death <- if (!is.na(date_death_col)) asDate2(DT[[date_death_col]]) else as.Date(NA)
  censor <- do.call(pmin, c(list(date_lost, date_death, rep(move_date_follow_end, nrow(DT))), na.rm = TRUE))
  censor[is.na(censor)] <- move_date_follow_end

  out[, y_date := y_date]
  out[, y_prevalent := as.integer(!is.na(y_date) & !is.na(date_attend) & y_date <= date_attend)]
  out[, y_event := as.integer(!is.na(y_date) & !is.na(date_attend) & y_date > date_attend)]
  out[, y_censor_date := censor]
  y_end_date <- censor
  idx_event <- which(out$y_event == 1 & !is.na(y_date))
  y_end_date[idx_event] <- y_date[idx_event]
  out[, y_end_date := y_end_date]
  out[, y_t2e := as.numeric(y_end_date - date_attend) / 365.25]
  out[!is.finite(y_t2e) | y_t2e <= 0, y_t2e := NA_real_]
  out[, y_source := if (length(date_cols)) paste(date_cols, collapse = ";") else "no_date_column_found"]
  out
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 🚩 Step 1: data preparation
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

run_if("data_prep", {
  set_output_scope("root")
  DT <- read_main_data()
  prot <- if (prot_do) read_layer("prot.rds", required = FALSE) else NULL
  met  <- if (met_do)  read_layer("met.rds", required = FALSE) else NULL
  prot_eid <- if (!is.null(prot)) unique(prot$eid) else character()
  met_eid <- if (!is.null(met)) unique(met$eid) else character()

  if (!"tdi" %in% names(DT)) stop("all.rds must contain baseline/current tdi.", call. = FALSE)
  tdi_q <- quantile(safe_num(DT$tdi), c(rich_poor_q, 1 - rich_poor_q), na.rm = TRUE)
  names(tdi_q) <- c("rich_cut_low_tdi", "poor_cut_high_tdi")
  cat("TDI rich/poor cuts: rich <=", tdi_q[1], "; poor >=", tdi_q[2], "\n")

  home_east_col <- pick_col(DT, c("home_east", "home_easting", "east", "home_x"))
  home_north_col <- pick_col(DT, c("home_north", "home_northing", "north", "home_y"))
  birth_east_col <- pick_col(DT, c("birth_east", "place_birth_east", "birth_easting"))
  birth_north_col <- pick_col(DT, c("birth_north", "place_birth_north", "birth_northing"))

  core <- data.table(
    eid = DT$eid,
    date_attend = DT$date_attend,
    age = if ("age" %in% names(DT)) safe_num(DT$age) else NA_real_,
    sex = if ("sex" %in% names(DT)) as.factor(DT$sex) else NA,
    center = if ("center" %in% names(DT)) as.factor(DT$center) else NA,
    PC1 = if ("PC1" %in% names(DT)) safe_num(DT$PC1) else NA_real_,
    PC2 = if ("PC2" %in% names(DT)) safe_num(DT$PC2) else NA_real_,
    current_tdi = safe_num(DT$tdi),
    current_east = if (!is.na(home_east_col)) safe_num(DT[[home_east_col]]) else NA_real_,
    current_north = if (!is.na(home_north_col)) safe_num(DT[[home_north_col]]) else NA_real_,
    birth_east = if (!is.na(birth_east_col)) safe_num(DT[[birth_east_col]]) else NA_real_,
    birth_north = if (!is.na(birth_north_col)) safe_num(DT[[birth_north_col]]) else NA_real_
  )
  core[, current_area := make_area_class(current_tdi, tdi_q[1], tdi_q[2])]

  # Birth-place TDI, using baseline home grid as a TDI reference.
  ref_tdi <- data.table(
    east = core$current_east,
    north = core$current_north,
    tdi = core$current_tdi
  )
  birth_query <- data.table(eid = core$eid, east = core$birth_east, north = core$birth_north)
  birth_map <- map_tdi_to_coords(birth_query, ref_tdi, max_km = main_nn_km, prefix = "birth_tdi")
  core <- merge(core, birth_map[, .(eid, birth_tdi = birth_tdi_mapped, birth_tdi_nn_km)], by = "eid", all.x = TRUE, sort = FALSE)
  core[, birth_area := make_area_class(birth_tdi, tdi_q[1], tdi_q[2])]
  core[, birth_current_transition_all := make_transition(birth_area, current_area)]
  core[, birth_current_transition := make_transition_rp(birth_area, current_area)]
  core[, birth_current_tdi_delta := current_tdi - birth_tdi]
  core[, birth_current_tdi_improve := birth_tdi - current_tdi]  # positive = moved toward lower-TDI/richer area
  core[, birth_current_distance_km := km_between(birth_east, birth_north, current_east, current_north)]

  # Home-history baseline-to-future transition.
  home_long <- make_home_long(DT)
  if (nrow(home_long)) {
    q_coord <- unique(home_long[coord_ok == TRUE, .(eid, east = hist_east, north = hist_north)])
    q_map <- map_tdi_to_coords(q_coord[, .(eid, east, north)], ref_tdi, max_km = main_nn_km, prefix = "hist_tdi")
    home_long <- merge(home_long, q_map[, .(eid, hist_east = east, hist_north = north, hist_tdi = hist_tdi_mapped, hist_tdi_nn_km)],
                       by = c("eid", "hist_east", "hist_north"), all.x = TRUE, sort = FALSE)
    home_long[, hist_area := make_area_class(hist_tdi, tdi_q[1], tdi_q[2])]

    base_home <- home_long[coord_ok == TRUE & !is.na(home_date) & !is.na(date_attend) & home_date <= date_attend]
    setorder(base_home, eid, home_date, slot)
    base_home <- base_home[, .SD[.N], by = eid]
    base_home <- base_home[, .(eid, base_home_date = home_date, base_home_east = hist_east, base_home_north = hist_north,
                               base_home_tdi = hist_tdi, base_home_area = hist_area, base_home_tdi_nn_km = hist_tdi_nn_km)]
    post_home <- home_long[coord_ok == TRUE & !is.na(home_date) & !is.na(date_attend) & home_date > date_attend]
    setorder(post_home, eid, home_date, slot)
    post_count <- post_home[, .(n_post_home_records = .N), by = eid]
    future_home <- post_home[, .SD[.N], by = eid]
    future_home <- future_home[, .(eid, future_home_date = home_date, future_home_east = hist_east, future_home_north = hist_north,
                                   future_home_tdi = hist_tdi, future_home_area = hist_area, future_home_tdi_nn_km = hist_tdi_nn_km)]
    core <- merge(core, base_home, by = "eid", all.x = TRUE, sort = FALSE)
    core <- merge(core, post_count, by = "eid", all.x = TRUE, sort = FALSE)
    core <- merge(core, future_home, by = "eid", all.x = TRUE, sort = FALSE)
    core[is.na(n_post_home_records), n_post_home_records := 0L]
    core[, home_future_transition_all := make_transition(base_home_area, future_home_area)]
    core[, home_future_transition := make_transition_rp(base_home_area, future_home_area)]
    core[, home_future_tdi_delta := future_home_tdi - base_home_tdi]
    core[, home_future_tdi_improve := base_home_tdi - future_home_tdi]
    core[, home_future_distance_km := km_between(base_home_east, base_home_north, future_home_east, future_home_north)]
  } else {
    home_long <- data.table(eid = character())
    core[, `:=`(
      base_home_date = as.Date(NA), base_home_east = NA_real_, base_home_north = NA_real_,
      base_home_tdi = NA_real_, base_home_area = NA_character_, base_home_tdi_nn_km = NA_real_,
      n_post_home_records = NA_integer_,
      future_home_date = as.Date(NA), future_home_east = NA_real_, future_home_north = NA_real_,
      future_home_tdi = NA_real_, future_home_area = NA_character_, future_home_tdi_nn_km = NA_real_,
      home_future_transition_all = factor(NA), home_future_transition = factor(NA),
      home_future_tdi_delta = NA_real_, home_future_tdi_improve = NA_real_, home_future_distance_km = NA_real_
    )]
  }

  # Choose primary transition for downstream analyses.
  if (main_transition == "birth_current") {
    core[, source_tdi := birth_tdi]
    core[, dest_tdi := current_tdi]
    core[, source_area := birth_area]
    core[, dest_area := current_area]
    core[, tdi_transition := birth_current_transition]
    core[, tdi_transition_all := birth_current_transition_all]
    core[, tdi_delta := birth_current_tdi_delta]
    core[, tdi_improve := birth_current_tdi_improve]
    core[, move_distance_km := birth_current_distance_km]
  } else {
    core[, source_tdi := base_home_tdi]
    core[, dest_tdi := future_home_tdi]
    core[, source_area := base_home_area]
    core[, dest_area := future_home_area]
    core[, tdi_transition := home_future_transition]
    core[, tdi_transition_all := home_future_transition_all]
    core[, tdi_delta := home_future_tdi_delta]
    core[, tdi_improve := home_future_tdi_improve]
    core[, move_distance_km := home_future_distance_km]
  }

  ses <- make_ses(DT)
  ses_sources <- attr(ses, "ses_sources")
  core <- merge(core, ses, by = "eid", all.x = TRUE, sort = FALSE)

  ydat <- make_y_outcome(DT, Y)
  core <- merge(core, ydat[, .(eid, y_date, y_prevalent, y_event, y_t2e, y_source)], by = "eid", all.x = TRUE, sort = FALSE)
  core[, has_prot := eid %in% prot_eid]
  core[, has_met := eid %in% met_eid]

  # Analysis groups that directly encode the hypothesis.
  core[, low_ses := person_ses_class == "low"]
  core[, high_ses := person_ses_class == "high"]
  core[, contrast_up_vs_stable_poor := fifelse(tdi_transition == "poor -> rich", 1L, fifelse(tdi_transition == "poor -> poor", 0L, NA_integer_))]
  core[, contrast_up_vs_down := fifelse(tdi_transition == "poor -> rich", 1L, fifelse(tdi_transition == "rich -> poor", 0L, NA_integer_))]
  core[, contrast_up_lowSES_vs_rr_lowSES := fifelse(tdi_transition == "poor -> rich" & low_ses == TRUE, 1L,
                                                    fifelse(tdi_transition == "rich -> rich" & low_ses == TRUE, 0L, NA_integer_))]
  core[, hypothesis_group := fcase(
    tdi_transition == "poor -> rich" & low_ses == TRUE, "poor_to_rich_lowSES",
    tdi_transition == "poor -> rich" & high_ses == TRUE, "poor_to_rich_highSES",
    tdi_transition == "rich -> rich" & low_ses == TRUE, "rich_to_rich_lowSES",
    tdi_transition == "rich -> rich" & high_ses == TRUE, "rich_to_rich_highSES",
    tdi_transition == "poor -> poor" & low_ses == TRUE, "poor_to_poor_lowSES",
    tdi_transition == "poor -> poor" & high_ses == TRUE, "poor_to_poor_highSES",
    tdi_transition == "rich -> poor" & low_ses == TRUE, "rich_to_poor_lowSES",
    tdi_transition == "rich -> poor" & high_ses == TRUE, "rich_to_poor_highSES",
    default = NA_character_
  )]
  core[, hypothesis_group_label := clean_hypothesis_label(hypothesis_group)]

  meta <- list(
    Y = Y,
    main_transition = main_transition,
    move_date_follow_end = move_date_follow_end,
    rich_poor_q = rich_poor_q,
    tdi_cuts = data.frame(cut = names(tdi_q), value = as.numeric(tdi_q)),
    coordinate_columns = data.frame(role = c("home_east", "home_north", "birth_east", "birth_north"),
                                    column = c(home_east_col, home_north_col, birth_east_col, birth_north_col)),
    ses_sources = ses_sources,
    y_source = unique(core$y_source)
  )

  save_rds_out(core, "move_core.rds", compress = "xz")
  save_rds_out(core, paste0("move_core.", Y_safe, ".rds"), compress = "xz")
  save_rds_out(home_long, "move_home_long.rds", compress = "xz")
  save_rds_out(meta, "move_meta.rds", compress = "xz")
  save_rds_out(meta, paste0("move_meta.", Y_safe, ".rds"), compress = "xz")
  write_raw(core, "move_core.tsv.gz")
  write_raw(meta$tdi_cuts, "tdi_rich_poor_cuts.tsv")
  write_raw(meta$coordinate_columns, "coordinate_columns.tsv")
  write_raw(meta$ses_sources, "ses_sources.tsv")

  prep_summary <- rbindlist(list(
    data.table(item = "all_rds_N", value = nrow(core)),
    data.table(item = "has_prot_N", value = sum(core$has_prot, na.rm = TRUE)),
    data.table(item = "has_met_N", value = sum(core$has_met, na.rm = TRUE)),
    data.table(item = "main_transition_nonmissing_N", value = sum(!is.na(core$tdi_transition)))
  ), fill = TRUE)
  write_book(list(prep_summary = prep_summary,
                  transition_counts = core[, .N, by = tdi_transition][order(tdi_transition)],
                  transition_all_counts = core[, .N, by = tdi_transition_all][order(tdi_transition_all)],
                  tdi_cuts = meta$tdi_cuts,
                  coordinate_columns = meta$coordinate_columns,
                  ses_sources = meta$ses_sources),
             "Fig0_data_prep.out.xlsx")
  print(prep_summary)
})

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 🚩 Step 2: QC plots and tables
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


run_if("qc", {
  set_output_scope("root")
  DT <- read_main_data()
  core <- read_core()
  core <- as.data.table(core)
  home_long <- read_home_long()

  # Job history QC, compressed to three useful panels.
  cat("\n# Job history QC\n")
  job_start_cols <- grep("^job_year0_a[0-9]+$", names(DT), value = TRUE)
  job_end_cols   <- grep("^job_year1_a[0-9]+$", names(DT), value = TRUE)
  job_code_cols  <- grep("^job_code_a[0-9]+$",  names(DT), value = TRUE)
  job_soc_cols   <- grep("^job_soc_a[0-9]+$",   names(DT), value = TRUE)
  job_idx <- sort(unique(c(slot_id(job_start_cols, "job_year0"), slot_id(job_end_cols, "job_year1"),
                           slot_id(job_code_cols, "job_code"), slot_id(job_soc_cols, "job_soc"))))
  job_idx <- job_idx[!is.na(job_idx)]
  job_start <- mk_slot_mat(DT, "job_year0", job_idx, 1900, 2037)
  job_end   <- mk_slot_mat(DT, "job_year1", job_idx, 1900, 2037)
  job_code  <- mk_slot_mat(DT, "job_code",  job_idx)
  job_soc   <- mk_slot_mat(DT, "job_soc",   job_idx)
  job_valid <- (!is.na(job_start)) | (!is.na(job_end)) | (!is.na(job_code)) | (!is.na(job_soc))
  attend_year <- DT$attend_year
  job_sum <- data.table(
    eid = DT$eid,
    n_job_valid = as.integer(rowSums(job_valid, na.rm = TRUE)),
    first_job_start_year = row_min_na(job_start), last_job_start_year = row_max_na(job_start),
    first_job_end_year = row_min_na(job_end), last_job_end_year = row_max_na(job_end),
    n_job_start_post_baseline = as.integer(rowSums(!is.na(job_start) & job_start > attend_year, na.rm = TRUE)),
    n_job_spanning_baseline = as.integer(rowSums(!is.na(job_start) & job_start <= attend_year & (is.na(job_end) | job_end >= attend_year), na.rm = TRUE)),
    n_job_code_changes = as.integer(row_adj_change(job_code)),
    n_job_soc_changes = as.integer(row_adj_change(job_soc))
  )
  job_sum[, job_change_n_derived := as.integer(pmax(n_job_valid - 1L, 0L))]
  job_sum[, job_timing_class := fifelse(n_job_valid == 0, "No job history",
                                        fifelse(n_job_start_post_baseline > 0, "Job start after baseline",
                                                fifelse(n_job_spanning_baseline > 0, "Job spans baseline", "Other/unclear")))]
  job_sum[, job_change_cat := fifelse(job_change_n_derived == 0, "0",
                                      fifelse(job_change_n_derived == 1, "1",
                                              fifelse(job_change_n_derived == 2, "2",
                                                      fifelse(job_change_n_derived <= 4, "3-4", "5+"))))]
  job_timing <- job_sum[, .(n = .N, pct = .N / nrow(job_sum), median_job_records = safe_median(n_job_valid),
                            median_job_changes = safe_median(job_change_n_derived),
                            median_spanning_jobs = safe_median(n_job_spanning_baseline)), by = job_timing_class][order(-n)]
  job_change_tab <- job_sum[, .N, by = .(job_timing_class, job_change_cat)][order(job_timing_class, job_change_cat)]
  job_slot <- data.table(slot = job_idx)
  if (length(job_idx)) {
    job_slot[, start_n := colSums(!is.na(job_start))]
    job_slot[, end_n := colSums(!is.na(job_end))]
    job_slot[, median_start_year := apply(job_start, 2, safe_median)]
    job_slot[, median_end_year := apply(job_end, 2, safe_median)]
    job_slot[, pct_start_after_baseline := sapply(seq_len(ncol(job_start)), function(j) { ok <- !is.na(job_start[, j]) & !is.na(attend_year); if (any(ok)) mean(job_start[ok, j] > attend_year[ok]) else NA_real_ })]
  }

  # Home history QC.
  cat("\n# Home history QC\n")
  if (nrow(home_long)) {
    home_sum <- home_long[, .(n_home_records = .N, n_home_dates = sum(!is.na(home_date)),
                              n_home_pre = sum(home_vs_baseline == "Before baseline", na.rm = TRUE),
                              n_home_same = sum(home_vs_baseline == "Same day", na.rm = TRUE),
                              n_home_post = sum(home_vs_baseline == "After baseline", na.rm = TRUE),
                              first_home_date = min_date2(home_date), last_home_date = max_date2(home_date),
                              n_coord_records = sum(coord_ok),
                              n_coord_changes = sum(coord_changed_from_prev, na.rm = TRUE),
                              total_hist_move_km = sum(hist_step_km, na.rm = TRUE)), by = eid]
    home_sum <- merge(data.table(eid = DT$eid), home_sum, by = "eid", all.x = TRUE, sort = FALSE)
    home_sum[, home_timing_class := fifelse(is.na(n_home_records), "No home history",
                                            fifelse(n_home_post > 0 & n_home_pre == 0, "Only post-baseline",
                                                    fifelse(n_home_pre > 0 & n_home_post == 0, "Only pre-baseline",
                                                            fifelse(n_home_pre > 0 & n_home_post > 0, "Both pre/post", "Same-day/missing"))))]
    home_timing <- home_sum[, .(n = .N, pct = .N / nrow(home_sum), median_records = safe_median(n_home_records),
                                median_coord_changes = safe_median(n_coord_changes), median_total_move_km = safe_median(total_hist_move_km)), by = home_timing_class][order(-n)]
    home_relative <- home_long[!is.na(home_date) & !is.na(date_attend), .(n = .N), by = home_vs_baseline][, pct := n / sum(n)][order(-n)]
    home_slot <- home_long[!is.na(home_date), .(n = .N, median_rel_year = safe_median(rel_year),
                                                p10_rel_year = quantile(rel_year, .10, na.rm = TRUE),
                                                p90_rel_year = quantile(rel_year, .90, na.rm = TRUE),
                                                min_date = min_date2(home_date), max_date = max_date2(home_date)), by = slot][order(slot)]
  } else {
    home_sum <- data.table(eid = DT$eid, home_timing_class = "No home history")
    home_timing <- home_sum[, .(n = .N), by = home_timing_class]
    home_relative <- data.table(home_vs_baseline = character(), n = integer(), pct = numeric())
    home_slot <- data.table(slot = integer())
  }

  # Main move/transition QC.
  transition_tab <- core[!is.na(tdi_transition), .N, by = tdi_transition][order(tdi_transition)]
  transition_all_tab <- core[!is.na(tdi_transition_all), .N, by = tdi_transition_all][order(tdi_transition_all)]
  omics_transition <- core[!is.na(tdi_transition), .(N = .N, proteomics = sum(has_prot), metabolomics = sum(has_met), both = sum(has_prot & has_met),
                                                    pct_prot = mean(has_prot), pct_met = mean(has_met)), by = tdi_transition][order(tdi_transition)]
  omics_long <- melt(omics_transition, id.vars = "tdi_transition", measure.vars = intersect(c("proteomics", "metabolomics", "both"), names(omics_transition)), variable.name = "omics", value.name = "N")

  # 3 x 3 summary figure: job, home, and TDI.
  p1 <- ggplot(job_timing, aes(reorder(job_timing_class, n), n, fill = job_timing_class)) +
    geom_col(width = .72, color = "white") +
    geom_text(aes(label = label_n_pct(n, nrow(job_sum))), hjust = -0.03, size = 3.2, fontface = "bold") +
    coord_flip(clip = "off") + scale_fill_manual(values = c("No job history" = "grey75", "Job spans baseline" = BLUE, "Job start after baseline" = GOLD, "Other/unclear" = TEAL), guide = "none") +
    scale_y_continuous(labels = comma, expand = expansion(mult = c(0, .30))) +
    labs(title = "a. Job-history timing", x = NULL, y = "Participants") + theme_5c(10)
  p2 <- ggplot(job_sum[, .N, by = job_change_cat], aes(factor(job_change_cat, c("0", "1", "2", "3-4", "5+")), N, fill = job_change_cat)) +
    geom_col(width = .70, color = "white") + scale_fill_manual(values = c("0" = "grey75", "1" = BLUE, "2" = TEAL, "3-4" = GOLD, "5+" = ORANGE), guide = "none") +
    scale_y_continuous(labels = comma) + labs(title = "b. Derived job changes", x = "Number of job changes", y = "Participants") + theme_5c(10)
  p3 <- ggplot(job_slot[slot <= 16], aes(slot)) +
    geom_hline(yintercept = median(DT$attend_year, na.rm = TRUE), linetype = "dotted", color = "grey40") +
    geom_line(aes(y = median_start_year, color = "Start"), linewidth = 1) + geom_point(aes(y = median_start_year, color = "Start"), size = 1.8) +
    geom_line(aes(y = median_end_year, color = "End"), linewidth = 1, linetype = "dashed") + geom_point(aes(y = median_end_year, color = "End"), size = 1.8) +
    scale_color_manual(values = c(Start = BLUE, End = GOLD)) + labs(title = "c. Median job years by slot", x = "Employment-history slot", y = "Calendar year", color = NULL) + theme_5c(10)
  p4 <- if (nrow(home_long)) ggplot(home_long[!is.na(rel_year) & rel_year >= -5 & rel_year <= 15], aes(rel_year)) +
    geom_histogram(binwidth = .25, boundary = 0, fill = BLUE, color = "white", linewidth = .05) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey25") +
    labs(title = "d. Home-history records", x = "Years from baseline", y = "Records") + scale_y_continuous(labels = comma) + theme_5c(10) else ggplot() + labs(title = "d. No home history") + theme_5c(10)
  p5 <- ggplot(home_timing, aes(reorder(home_timing_class, n), n, fill = home_timing_class)) +
    geom_col(width = .72, color = "white") + geom_text(aes(label = label_n_pct(n, nrow(home_sum))), hjust = -0.03, size = 3.1, fontface = "bold") +
    coord_flip(clip = "off") + scale_fill_manual(values = c("Only pre-baseline" = BLUE, "Both pre/post" = GOLD, "No home history" = "grey70", "Same-day/missing" = TEAL, "Only post-baseline" = ORANGE), guide = "none") +
    scale_y_continuous(labels = comma, expand = expansion(mult = c(0, .32))) + labs(title = "e. Person-level home timing", x = NULL, y = "Participants") + theme_5c(10)
  p6 <- if (nrow(home_long)) ggplot(home_long[!is.na(rel_year) & slot <= 10 & rel_year >= -5 & rel_year <= 15], aes(factor(slot), rel_year)) +
    geom_boxplot(aes(fill = factor(slot)), outlier.size = .12, color = "grey35", show.legend = FALSE) + geom_hline(yintercept = 0, linetype = "dashed", color = "grey35") +
    labs(title = "f. Home-date timing by slot", x = "Home-history slot", y = "Years from baseline") + theme_5c(10) else ggplot() + labs(title = "f. No slot timing") + theme_5c(10)
  p7 <- ggplot(transition_tab, aes(tdi_transition, N, fill = tdi_transition)) + geom_col(width = .70, color = "white") +
    geom_text(aes(label = comma(N)), vjust = -0.35, size = 3.2, fontface = "bold") +
    scale_fill_manual(values = cols_transition, guide = "none") + scale_y_continuous(labels = comma, expand = expansion(mult = c(0, .15))) +
    scale_x_discrete(labels = short_transition) + labs(title = paste0("g. Rich/poor TDI transition: ", main_transition), x = "Source → destination", y = "Participants") + theme_5c(10)
  p8 <- ggplot(core[!is.na(tdi_improve) & !is.na(tdi_transition)], aes(tdi_improve, fill = tdi_transition)) +
    geom_histogram(bins = 55, alpha = .78, color = "white", linewidth = .05) + geom_vline(xintercept = 0, linetype = "dashed", color = "grey25") +
    scale_fill_manual(values = cols_transition, name = NULL, labels = short_transition) + labs(title = "h. TDI improvement", x = "Source TDI − destination TDI", y = "Participants") + theme_5c(10)
  p9 <- ggplot(omics_long, aes(tdi_transition, N, fill = omics)) + geom_col(position = position_dodge(.74), width = .65, color = "white") +
    scale_fill_manual(values = c(proteomics = GOLD, metabolomics = BLUE, both = TEAL), labels = c(proteomics = "Proteomics", metabolomics = "Metabolomics", both = "Both")) +
    scale_x_discrete(labels = short_transition) + scale_y_continuous(labels = comma) + labs(title = "i. Omics sample size", x = "Source → destination", y = "N", fill = NULL) + theme_5c(10)
  Fig1 <- (p1 | p2 | p3) / (p4 | p5 | p6) / (p7 | p8 | p9)
  save_plot(Fig1, "Fig1.job_home_TDI_summary.png", 17, 14)

  write_book(list(
    job_timing_summary = job_timing,
    job_change_table = job_change_tab,
    job_slot_summary = job_slot,
    home_timing_summary = home_timing,
    record_relative_timing = home_relative,
    home_slot_summary = home_slot,
    main_transition = transition_tab,
    all_transition = transition_all_tab,
    omics_by_transition = omics_transition,
    note = data.frame(note = c(
      "Fig1 compresses the old Fig1-Fig3 into one 3 x 3 summary figure.",
      "Job spans baseline means at least one recorded job interval has start_year <= baseline year and end_year missing or >= baseline year; it is not proof of post-baseline job change.",
      "The peak around baseline in home history likely reflects UKB-maintained address records around invitation/recruitment rather than lifetime residential history."
    ))
  ), "Fig1.out.xlsx")
})

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 🚩 Step 3: SES and disease-risk analyses
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

model_covs <- function(dat, extra = character()) {
  covs <- c("age", "sex", "PC1", "PC2", extra)
  covs <- covs[covs %in% names(dat)]
  covs <- covs[vapply(covs, function(v) sum(!is.na(dat[[v]])) > 0 && length(unique(na.omit(dat[[v]]))) > 1, logical(1))]
  unique(covs)
}

tidy_lm_term <- function(fit) {
  sm <- summary(fit)$coefficients
  data.table(term = rownames(sm), beta = sm[, "Estimate"], se = sm[, "Std. Error"], statistic = sm[, "t value"], p = sm[, "Pr(>|t|)"])
}
tidy_glm_term <- function(fit) {
  sm <- summary(fit)$coefficients
  data.table(term = rownames(sm), beta = sm[, "Estimate"], se = sm[, "Std. Error"], statistic = sm[, "z value"], p = sm[, "Pr(>|z|)"], OR = exp(sm[, "Estimate"]), lo = exp(sm[, "Estimate"] - 1.96 * sm[, "Std. Error"]), hi = exp(sm[, "Estimate"] + 1.96 * sm[, "Std. Error"]))
}
tidy_cox_term <- function(fit) {
  sm <- summary(fit)
  cf <- sm$coefficients
  ci <- sm$conf.int
  data.table(term = rownames(cf), beta = cf[, "coef"], se = cf[, "se(coef)"], z = cf[, "z"], p = cf[, "Pr(>|z|)"], HR = ci[, "exp(coef)"], lo = ci[, "lower .95"], hi = ci[, "upper .95"])
}


run_if("ses_disease", {
  set_output_scope("trait")
  core <- read_core()
  core <- as.data.table(core)
  core[, tdi_transition := factor(as.character(tdi_transition), levels = names(cols_transition))]
  core[, person_ses_class := factor(as.character(person_ses_class), levels = c("low", "middle", "high"))]
  core[, ses_extreme := fcase(person_ses_class == "low", "Low SES", person_ses_class == "high", "High SES", default = NA_character_)]
  core[, ses_extreme := factor(ses_extreme, levels = c("Low SES", "High SES"))]
  core[, hypothesis_group_label := clean_hypothesis_label(hypothesis_group)]

  dat_rp <- core[!is.na(tdi_transition)]
  if (!nrow(dat_rp)) stop("No rich/poor transition records available for ses_disease step.", call. = FALSE)

  desc_transition <- dat_rp[, .(
    N = .N,
    age_mean = safe_mean(age),
    female_pct = mean(as.character(sex) %in% c("Female", "0", "2"), na.rm = TRUE),
    source_tdi_mean = safe_mean(source_tdi),
    dest_tdi_mean = safe_mean(dest_tdi),
    tdi_improve_mean = safe_mean(tdi_improve),
    move_km_median = safe_median(move_distance_km),
    individual_ses_mean = safe_mean(individual_ses_z),
    low_ses_pct = mean(person_ses_class == "low", na.rm = TRUE),
    high_ses_pct = mean(person_ses_class == "high", na.rm = TRUE),
    prevalent_Y_pct = mean(y_prevalent == 1, na.rm = TRUE),
    incident_Y_pct = mean(y_event == 1, na.rm = TRUE),
    incident_Y_events = sum(y_event == 1, na.rm = TRUE),
    proteomics_N = sum(has_prot), metabolomics_N = sum(has_met)
  ), by = tdi_transition][order(tdi_transition)]
  desc_transition[, transition_label := short_transition(tdi_transition)]

  desc_ses <- dat_rp[!is.na(ses_extreme), .(
    N = .N,
    incident_Y_pct = mean(y_event == 1, na.rm = TRUE),
    incident_Y_events = sum(y_event == 1, na.rm = TRUE),
    individual_ses_mean = safe_mean(individual_ses_z),
    proteomics_N = sum(has_prot), metabolomics_N = sum(has_met)
  ), by = .(tdi_transition, ses_extreme)][order(tdi_transition, ses_extreme)]
  desc_ses[, transition_label := short_transition(tdi_transition)]

  desc_hypothesis <- dat_rp[!is.na(hypothesis_group), .(
    N = .N,
    individual_ses_mean = safe_mean(individual_ses_z),
    prevalent_Y_pct = mean(y_prevalent == 1, na.rm = TRUE),
    incident_Y_pct = mean(y_event == 1, na.rm = TRUE),
    incident_Y_events = sum(y_event == 1, na.rm = TRUE),
    move_km_median = safe_median(move_distance_km),
    proteomics_N = sum(has_prot), metabolomics_N = sum(has_met)
  ), by = .(hypothesis_group, hypothesis_group_label)][order(hypothesis_group)]

  # Individual SES models.
  covs_ses <- model_covs(dat_rp, extra = c("source_tdi", "dest_tdi"))
  ses_lm_tab <- data.table()
  ses_model_dat <- dat_rp[, unique(c("individual_ses_z", "tdi_transition", covs_ses)), with = FALSE]
  ses_model_dat <- ses_model_dat[complete.cases(ses_model_dat)]
  if (nrow(ses_model_dat) >= min_n_model) {
    ses_formula <- as.formula(paste("individual_ses_z ~ tdi_transition", if (length(covs_ses)) paste("+", paste(bt(covs_ses), collapse = " + ")) else ""))
    ses_lm <- tryCatch(lm(ses_formula, data = ses_model_dat), error = function(e) NULL)
    if (!is.null(ses_lm)) ses_lm_tab <- tidy_lm_term(ses_lm)
  }

  low_glm_tab <- data.table()
  low_model_dat <- dat_rp[, unique(c("low_ses", "tdi_transition", covs_ses)), with = FALSE]
  low_model_dat <- low_model_dat[!is.na(low_ses) & complete.cases(low_model_dat)]
  if (nrow(low_model_dat) >= min_n_model && length(unique(low_model_dat$low_ses)) == 2) {
    low_glm <- tryCatch(glm(as.formula(paste("low_ses ~ tdi_transition", if (length(covs_ses)) paste("+", paste(bt(covs_ses), collapse = " + ")) else "")), data = low_model_dat, family = binomial()), error = function(e) NULL)
    if (!is.null(low_glm)) low_glm_tab <- tidy_glm_term(low_glm)
  }

  # Cox models for disease risk.
  surv_base <- dat_rp[!is.na(y_t2e) & !is.na(y_event) & y_t2e > 0 & y_prevalent != 1]
  covs_base <- model_covs(surv_base, extra = c("source_tdi", "dest_tdi", "individual_ses_z"))
  cox_transition <- cox_group <- cox_interaction <- data.table()
  risk10 <- data.table()
  if (nrow(surv_base) >= min_n_model && sum(surv_base$y_event == 1, na.rm = TRUE) >= min_case_model) {
    surv_base[, tdi_transition := relevel(tdi_transition, ref = "poor -> poor")]
    f1 <- as.formula(paste0("survival::Surv(y_t2e, y_event) ~ tdi_transition",
                            if (length(covs_base)) paste0(" + ", paste(bt(covs_base), collapse = " + ")) else ""))
    fit1 <- tryCatch(survival::coxph(f1, data = surv_base), error = function(e) { message("transition Cox failed: ", conditionMessage(e)); NULL })
    if (!is.null(fit1)) {
      cox_transition <- tidy_cox_term(fit1)[grepl("^tdi_transition", term)]
      cox_transition[, `:=`(model = "TDI transition", label = paste0(clean_cox_label(term), " vs Poor → poor"))]
    }

    # Extreme SES interaction, excluding middle SES for interpretability.
    d_int <- surv_base[!is.na(ses_extreme)]
    if (nrow(d_int) >= min_n_model && sum(d_int$y_event == 1, na.rm = TRUE) >= min_case_model) {
      covs_int <- setdiff(model_covs(d_int, extra = c("source_tdi", "dest_tdi")), c("individual_ses_z"))
      f2 <- as.formula(paste0("survival::Surv(y_t2e, y_event) ~ tdi_transition * ses_extreme",
                              if (length(covs_int)) paste0(" + ", paste(bt(covs_int), collapse = " + ")) else ""))
      fit2 <- tryCatch(survival::coxph(f2, data = d_int), error = function(e) { message("transition x SES Cox failed: ", conditionMessage(e)); NULL })
      if (!is.null(fit2)) {
        cox_interaction <- tidy_cox_term(fit2)[grepl("tdi_transition|ses_extreme", term)]
        cox_interaction[, `:=`(model = "Transition × SES", label = clean_cox_label(term))]
      }
      rr <- plot_adjusted_10y_risk(d_int, covs = covs_int, t0 = as.numeric(Sys.getenv("MOVE_RISK_T0", unset = "10")), title = "b. Adjusted risk by transition and SES")
      p_risk <- rr$plot
      risk10 <- rr$table
    } else {
      p_risk <- ggplot() + labs(title = "b. Adjusted risk unavailable") + theme_5c()
    }

    # Hypothesis-group model with clean labels. Reference is Rich->rich + low SES.
    hyp <- surv_base[!is.na(hypothesis_group)]
    if (nrow(hyp) >= min_n_model && sum(hyp$y_event == 1, na.rm = TRUE) >= min_case_model) {
      hyp[, hypothesis_group := factor(hypothesis_group)]
      ref <- "rich_to_rich_lowSES"
      if (ref %in% hyp$hypothesis_group) hyp[, hypothesis_group := relevel(hypothesis_group, ref = ref)]
      covs_hyp <- model_covs(hyp, extra = c("source_tdi", "dest_tdi"))
      f3 <- as.formula(paste0("survival::Surv(y_t2e, y_event) ~ hypothesis_group",
                              if (length(covs_hyp)) paste0(" + ", paste(bt(covs_hyp), collapse = " + ")) else ""))
      fit3 <- tryCatch(survival::coxph(f3, data = hyp), error = function(e) { message("hypothesis-group Cox failed: ", conditionMessage(e)); NULL })
      if (!is.null(fit3)) {
        cox_group <- tidy_cox_term(fit3)[grepl("^hypothesis_group", term)]
        cox_group[, `:=`(model = "Hypothesis groups", label = paste0(clean_cox_label(term), " vs Rich → rich + low SES"))]
      }
    }
  } else {
    p_risk <- ggplot() + labs(title = "b. Adjusted risk unavailable") + theme_5c()
  }

  # Final disease-risk figure: four high-information panels.
  inc_base <- dat_rp[!is.na(y_event) & y_prevalent != 1]
  pA <- plot_crude_origin_tile_combined(inc_base, title = "a. Crude incidence by origin, destination and SES")

  od_contrasts <- make_origin_destination_contrasts(dat_rp[!is.na(ses_extreme)])
  od_origin <- od_contrasts[contrast_family == "Origin effect within same current community and SES"]
  od_dest <- od_contrasts[contrast_family == "Destination effect within same origin and SES"]
  pC <- plot_od_forest_combined(od_contrasts, title = "c. Origin and destination effects within matched strata", base_size = 10)
  pD <- plot_origin_gap(inc_base) + labs(title = "d. Excess crude risk from poor origin")

  Fig3 <- (pA | p_risk) / (pC | pD)
  save_plot(Fig3, fig_file(3, "SES_risk", trait = TRUE), 18, 10.5)

  # New paradox decomposition: does risk track origin, destination, improvement, or individual SES?
  cont_risk <- data.table()
  cont_base <- copy(surv_base)
  if (nrow(cont_base) >= min_n_model && sum(cont_base$y_event == 1, na.rm = TRUE) >= min_case_model) {
    cont_base[, source_tdi_z := zscore(source_tdi)]
    cont_base[, dest_tdi_z := zscore(dest_tdi)]
    cont_base[, tdi_improve_z := zscore(tdi_improve)]
    cont_base[, move_distance_z := zscore(move_distance_km)]
    cont_base[, individual_ses_z2 := zscore(individual_ses_z)]
    covs_cont <- setdiff(model_covs(cont_base, extra = c()), c("source_tdi", "dest_tdi", "tdi_improve", "move_distance_km", "individual_ses_z"))
    specs_cont <- list(
      `Origin + destination` = c("source_tdi_z", "dest_tdi_z", "individual_ses_z2"),
      `Improvement + distance` = c("tdi_improve_z", "move_distance_z", "individual_ses_z2")
    )
    cont_risk <- rbindlist(lapply(names(specs_cont), function(mn) {
      rhs <- c(specs_cont[[mn]], covs_cont)
      dd <- cont_base[, unique(c("y_t2e", "y_event", rhs)), with = FALSE]
      dd <- dd[complete.cases(dd)]
      if (nrow(dd) < min_n_model || sum(dd$y_event == 1, na.rm = TRUE) < min_case_model) return(data.table())
      fit <- tryCatch(survival::coxph(as.formula(paste0("survival::Surv(y_t2e, y_event) ~ ", paste(bt(rhs), collapse = " + "))), data = dd), error = function(e) NULL)
      if (is.null(fit)) return(data.table())
      out <- tidy_cox_term(fit)[term %in% specs_cont[[mn]]]
      out[, `:=`(model = mn, N = nrow(dd), events = sum(dd$y_event == 1, na.rm = TRUE))]
      out
    }), fill = TRUE)
    if (nrow(cont_risk)) {
      cont_risk[, label := dplyr::recode(term,
        source_tdi_z = "Poorer origin TDI",
        dest_tdi_z = "Poorer destination TDI",
        tdi_improve_z = "TDI improvement toward rich area",
        move_distance_z = "Move distance",
        individual_ses_z2 = "Higher individual SES",
        .default = term
      )]
    }
  }
  p6A <- if (nrow(cont_risk)) {
    dd <- copy(cont_risk); dd[, label := factor(label, levels = rev(unique(label)))]
    ggplot(dd, aes(HR, label)) +
      geom_vline(xintercept = 1, linetype = "dashed", color = "grey55") +
      geom_errorbar(aes(xmin = lo, xmax = hi), orientation = "y", width = .18, linewidth = .75, color = "grey45") +
      geom_point(aes(fill = p < 0.05), shape = 21, size = 2.8, color = "grey20") +
      facet_wrap(~model, scales = "free_y") + scale_x_log10() +
      scale_fill_manual(values = c(`TRUE` = GOLD, `FALSE` = "white"), guide = "none") +
      labs(title = "a. Continuous origin/destination decomposition", x = "Hazard ratio per 1 s.d.", y = NULL) + theme_5c(10)
  } else ggplot() + labs(title = "a. Continuous decomposition unavailable") + theme_5c(10)
  p6B <- ggplot(dat_rp[!is.na(tdi_transition) & is.finite(individual_ses_z)], aes(short_transition(tdi_transition), individual_ses_z, fill = tdi_transition)) +
    geom_violin(alpha = .35, color = NA, trim = TRUE) + geom_boxplot(width = .14, outlier.size = .12, alpha = .90) +
    scale_fill_manual(values = cols_transition, guide = "none") +
    labs(title = "b. Individual SES by TDI transition", x = "TDI transition", y = "Individual SES z-score") + theme_5c(10)
  risk_ses_tab <- dat_rp[!is.na(tdi_transition) & !is.na(ses_extreme) & y_prevalent != 1,
                         .(N = .N, events = sum(y_event == 1, na.rm = TRUE), risk = mean(y_event == 1, na.rm = TRUE)),
                         by = .(tdi_transition, ses_extreme)]
  p6C <- if (nrow(risk_ses_tab)) ggplot(risk_ses_tab, aes(short_transition(tdi_transition), risk, fill = ses_extreme)) +
    geom_col(position = position_dodge(.72), width = .65, color = "white") +
    geom_text(aes(label = percent(risk, accuracy = .1)), position = position_dodge(.72), vjust = -.35, size = 3.0, fontface = "bold") +
    scale_fill_manual(values = cols_ses[c("Low SES", "High SES")], name = NULL) + scale_y_continuous(labels = percent, expand = expansion(mult = c(0, .18))) +
    labs(title = "c. Crude risk by transition and SES", x = "TDI transition", y = "Incident risk") + theme_5c(10) else ggplot() + labs(title = "c. Crude risk unavailable") + theme_5c(10)
  improve_tab <- cont_base[is.finite(tdi_improve) & y_prevalent != 1]
  p6D <- if (nrow(improve_tab) >= min_n_model) {
    qs <- unique(quantile(improve_tab$tdi_improve, probs = seq(0, 1, 0.25), na.rm = TRUE))
    if (length(qs) >= 3) improve_tab[, improve_q := cut(tdi_improve, breaks = qs, include.lowest = TRUE, labels = paste0("Q", seq_len(length(qs) - 1)))] else improve_tab[, improve_q := factor(NA_character_)]
    tt <- improve_tab[!is.na(improve_q), .(N = .N, events = sum(y_event == 1, na.rm = TRUE), risk = mean(y_event == 1, na.rm = TRUE), mean_improve = mean(tdi_improve, na.rm = TRUE)), by = improve_q]
    ggplot(tt, aes(improve_q, risk, fill = mean_improve)) + geom_col(width = .68, color = "white") +
      geom_text(aes(label = percent(risk, accuracy = .1)), vjust = -.35, size = 3.0, fontface = "bold") +
      scale_fill_gradient2(low = BLUE, mid = "white", high = GOLD, midpoint = 0, name = "Mean\nTDI improvement") +
      scale_y_continuous(labels = percent, expand = expansion(mult = c(0, .18))) +
      labs(title = "d. Risk across TDI-improvement quartiles", x = "TDI-improvement quartile", y = "Incident risk") + theme_5c(10)
  } else ggplot() + labs(title = "d. TDI-improvement quartiles unavailable") + theme_5c(10)
  Fig7 <- (p6A | p6B) / (p6C | p6D)
  save_plot(Fig7, fig_file(7, "paradox", trait = TRUE), 16, 10)
  write_book(list(continuous_decomposition = cont_risk, risk_by_transition_SES = risk_ses_tab), fig_file(7, "paradox", ".out.xlsx", trait = TRUE))

  cox_all <- rbindlist(list(cox_transition, cox_interaction, cox_group), fill = TRUE)
  write_book(list(
    desc_by_transition = desc_transition,
    desc_by_transition_SES = desc_ses,
    desc_hypothesis_group = desc_hypothesis,
    adjusted_10y_risk = risk10,
    lm_individual_SES = ses_lm_tab,
    glm_low_SES = low_glm_tab,
    cox_transition = cox_transition,
    cox_transition_x_SES = cox_interaction,
    cox_hypothesis_group = cox_group,
    origin_destination_contrasts = od_contrasts,
    continuous_decomposition = cont_risk,
    risk_by_transition_SES = risk_ses_tab,
    model_notes = data.frame(note = c(
      "Fig3 now has four panels: combined crude origin/destination incidence, adjusted 10-year risk, combined origin/destination HRs, and poor-origin crude risk gap.",
      "highSES means the top individual-SES quartile; lowSES means the bottom individual-SES quartile; middle SES is not folded into highSES.",
      "rich = lowest TDI quartile; poor = highest TDI quartile; middle-TDI areas are excluded from main transition models.",
      "Cox models exclude prevalent disease and use incident disease after date_attend."
    ))
  ), fig_file(3, "SES_risk", ".out.xlsx", trait = TRUE))
  save_rds_out(list(desc_transition = desc_transition, desc_ses = desc_ses, desc_hypothesis = desc_hypothesis,
               ses_lm = ses_lm_tab, low_glm = low_glm_tab, cox = cox_all, risk10 = risk10,
               origin_destination_contrasts = od_contrasts, continuous_decomposition = cont_risk,
               risk_by_transition_SES = risk_ses_tab),
          "ses_disease.res.rds", compress = "xz")
})

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 🚩 Step 4: omics signatures
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

feature_list <- function(biom, max_features = 0) {
  xs <- setdiff(names(biom), "eid")
  if (max_features > 0) xs <- head(xs, max_features)
  xs
}

scan_lm_features <- function(dat, features, outcome, covs, layer, contrast_label, min_n = 500, min_case = 30) {
  covs <- covs[covs %in% names(dat)]
  base_cols <- unique(c("eid", outcome, covs))
  features <- intersect(features, names(dat))
  if (!length(features) || !outcome %in% names(dat)) return(data.table())
  y0 <- dat[[outcome]]
  ok_y <- !is.na(y0)
  n_case <- sum(y0 == 1, na.rm = TRUE)
  n_ctrl <- sum(y0 == 0, na.rm = TRUE)
  if (sum(ok_y) < min_n || n_case < min_case || n_ctrl < min_case) {
    message("Skip ", layer, " contrast ", contrast_label, ": N/cases too small (N=", sum(ok_y), ", cases=", n_case, ", controls=", n_ctrl, ").")
    return(data.table())
  }
  progress <- file.path(rawdir, paste0("omics_progress_", layer, "_", contrast_label, ".tsv"))
  writeLines("time\tindex\ttotal\tfeature\tstatus", progress)
  out <- vector("list", length(features))
  for (i in seq_along(features)) {
    x <- features[i]
    dd <- dat[, c(base_cols, x), with = FALSE]
    setnames(dd, x, "feature_value")
    dd[, feature_value := inormal2(feature_value)]
    dd <- dd[complete.cases(dd)]
    if (nrow(dd) < min_n || sum(dd[[outcome]] == 1) < min_case || sum(dd[[outcome]] == 0) < min_case || length(unique(dd$feature_value)) < 3) {
      out[[i]] <- data.table(layer = layer, contrast = contrast_label, feature = x, beta = NA_real_, se = NA_real_, statistic = NA_real_, p = NA_real_, N = nrow(dd), n_case = sum(dd[[outcome]] == 1), n_ctrl = sum(dd[[outcome]] == 0), model = "linear_probability")
    } else {
      form <- as.formula(paste(outcome, "~ feature_value", if (length(covs)) paste("+", paste(bt(covs), collapse = " + ")) else ""))
      fit <- tryCatch(lm(form, data = dd), error = function(e) NULL)
      if (is.null(fit)) {
        out[[i]] <- data.table(layer = layer, contrast = contrast_label, feature = x, beta = NA_real_, se = NA_real_, statistic = NA_real_, p = NA_real_, N = nrow(dd), n_case = sum(dd[[outcome]] == 1), n_ctrl = sum(dd[[outcome]] == 0), model = "linear_probability")
      } else {
        sm <- summary(fit)$coefficients
        if ("feature_value" %in% rownames(sm)) {
          out[[i]] <- data.table(layer = layer, contrast = contrast_label, feature = x, beta = sm["feature_value", "Estimate"], se = sm["feature_value", "Std. Error"], statistic = sm["feature_value", "t value"], p = sm["feature_value", "Pr(>|t|)"], N = nrow(dd), n_case = sum(dd[[outcome]] == 1), n_ctrl = sum(dd[[outcome]] == 0), model = "linear_probability")
        } else {
          out[[i]] <- data.table(layer = layer, contrast = contrast_label, feature = x, beta = NA_real_, se = NA_real_, statistic = NA_real_, p = NA_real_, N = nrow(dd), n_case = sum(dd[[outcome]] == 1), n_ctrl = sum(dd[[outcome]] == 0), model = "linear_probability")
        }
      }
    }
    if (i %% 100 == 0 || i == length(features)) {
      cat(format(Sys.time(), "%F %T"), i, length(features), x, "done", sep = "\t", file = progress, append = TRUE); cat("\n", file = progress, append = TRUE)
      message(layer, " ", contrast_label, ": ", i, "/", length(features))
    }
  }
  res <- rbindlist(out, fill = TRUE)
  res[, FDR := p.adjust(p, method = "BH")]
  res[order(p)]
}

scan_cox_features <- function(dat, features, covs, layer, min_n = 500, min_event = 30) {
  covs <- covs[covs %in% names(dat)]
  features <- intersect(features, names(dat))
  if (!length(features) || !all(c("y_t2e", "y_event") %in% names(dat))) return(data.table())
  d0 <- dat[!is.na(y_t2e) & !is.na(y_event) & y_t2e > 0 & y_prevalent != 1]
  if (nrow(d0) < min_n || sum(d0$y_event == 1, na.rm = TRUE) < min_event) return(data.table())
  progress <- file.path(rawdir, paste0("omics_progress_", layer, "_outcome_cox.tsv"))
  writeLines("time\tindex\ttotal\tfeature\tstatus", progress)
  base_cols <- unique(c("eid", "y_t2e", "y_event", covs))
  out <- vector("list", length(features))
  for (i in seq_along(features)) {
    x <- features[i]
    dd <- d0[, c(base_cols, x), with = FALSE]
    setnames(dd, x, "feature_value")
    dd[, feature_value := inormal2(feature_value)]
    dd <- dd[complete.cases(dd)]
    if (nrow(dd) < min_n || sum(dd$y_event == 1) < min_event || length(unique(dd$feature_value)) < 3) {
      out[[i]] <- data.table(layer = layer, feature = x, beta_cad = NA_real_, se_cad = NA_real_, p_cad = NA_real_, HR_cad = NA_real_, lo_cad = NA_real_, hi_cad = NA_real_, N_cad = nrow(dd), events_cad = sum(dd$y_event == 1))
    } else {
      form <- as.formula(paste0("survival::Surv(y_t2e, y_event) ~ feature_value", if (length(covs)) paste0(" + ", paste(bt(covs), collapse = " + ")) else ""))
      fit <- tryCatch(survival::coxph(form, data = dd), error = function(e) NULL)
      if (is.null(fit)) {
        out[[i]] <- data.table(layer = layer, feature = x, beta_cad = NA_real_, se_cad = NA_real_, p_cad = NA_real_, HR_cad = NA_real_, lo_cad = NA_real_, hi_cad = NA_real_, N_cad = nrow(dd), events_cad = sum(dd$y_event == 1))
      } else {
        tb <- tidy_cox_term(fit)
        tb <- tb[term == "feature_value"]
        if (nrow(tb)) {
          out[[i]] <- data.table(layer = layer, feature = x, beta_cad = tb$beta, se_cad = tb$se, p_cad = tb$p, HR_cad = tb$HR, lo_cad = tb$lo, hi_cad = tb$hi, N_cad = nrow(dd), events_cad = sum(dd$y_event == 1))
        } else {
          out[[i]] <- data.table(layer = layer, feature = x, beta_cad = NA_real_, se_cad = NA_real_, p_cad = NA_real_, HR_cad = NA_real_, lo_cad = NA_real_, hi_cad = NA_real_, N_cad = nrow(dd), events_cad = sum(dd$y_event == 1))
        }
      }
    }
    if (i %% 100 == 0 || i == length(features)) {
      cat(format(Sys.time(), "%F %T"), i, length(features), x, "done", sep = "\t", file = progress, append = TRUE); cat("\n", file = progress, append = TRUE)
      message(layer, " outcome Cox: ", i, "/", length(features))
    }
  }
  res <- rbindlist(out, fill = TRUE)
  res[, FDR_cad := p.adjust(p_cad, method = "BH")]
  res[order(p_cad)]
}


plot_volcano <- function(res, beta_col = "beta", p_col = "p", fdr_col = "FDR", title = "", top_n = 20, xlab = "Effect per 1 s.d. higher biomarker") {
  if (!nrow(res) || !all(c(beta_col, p_col) %in% names(res))) return(ggplot() + labs(title = paste(title, "unavailable")) + theme_5c())
  d <- as.data.table(res)
  d <- d[is.finite(get(beta_col)) & is.finite(get(p_col))]
  if (!nrow(d)) return(ggplot() + labs(title = paste(title, "unavailable")) + theme_5c())
  d[, logp := -log10(pmax(get(p_col), 1e-300))]
  if (fdr_col %in% names(d)) d[, sig := is.finite(get(fdr_col)) & get(fdr_col) < 0.05] else d[, sig := FALSE]
  d[, direction := fifelse(sig & get(beta_col) > 0, "FDR<0.05, positive",
                           fifelse(sig & get(beta_col) < 0, "FDR<0.05, inverse",
                                   fifelse(get(p_col) < 0.05, "P<0.05 only", "Not significant")))]
  lab <- data.table::rbindlist(list(
    d[sig == TRUE][order(get(p_col))][seq_len(min(top_n, .N))],
    d[order(get(p_col))][seq_len(min(ceiling(top_n / 2), .N))]
  ), fill = TRUE)
  if (nrow(lab)) lab <- unique(lab, by = "feature")
  ggplot(d, aes(x = .data[[beta_col]], y = logp)) +
    geom_hline(yintercept = -log10(0.05 / max(1, nrow(d))), linetype = "dashed", color = "grey45") +
    geom_vline(xintercept = 0, color = "grey70") +
    geom_point(aes(color = direction, size = sig), alpha = .82) +
    ggrepel::geom_text_repel(data = lab, aes(label = feature), size = 2.9, max.overlaps = top_n, seed = 11, min.segment.length = 0) +
    scale_color_manual(values = cols_signature, breaks = names(cols_signature)) +
    scale_size_manual(values = c(`TRUE` = 2.0, `FALSE` = 1.25), guide = "none") +
    labs(title = title, x = xlab, y = expression(-log[10](P)), color = NULL) + theme_5c(11)
}
plot_driver_heatmap <- function(driver, layer, top_n = 18) {
  d <- as.data.table(driver)
  if (!nrow(d)) return(ggplot() + labs(title = "d. Driver evidence unavailable") + theme_5c())
  d <- d[order(-priority_score_cont, p_up, p_up_lowSES, p_cad)]
  d <- d[seq_len(min(top_n, .N))]
  d[, feature := factor(feature, levels = rev(feature))]
  long <- melt(d[, .(feature, `Upward move` = score_up, `Low-SES
resilience` = score_lowSES,
                    `Opposite-move
contrast` = score_opposite, `Lower
outcome risk` = score_cad_protect)],
               id.vars = "feature", variable.name = "component", value.name = "score")
  ggplot(long, aes(component, feature, fill = score)) +
    geom_tile(color = "white") +
    geom_text(aes(label = sprintf("%.1f", score)), size = 2.7, fontface = "bold") +
    scale_fill_gradient2(low = BLUE, mid = "white", high = GOLD, midpoint = 0, name = "Favorable
signed
evidence") +
    labs(title = paste0("d. ", layer, " prioritized driver evidence"), x = NULL, y = NULL) +
    theme_5c(10) + theme(axis.text.y = element_text(face = "bold", size = 9))
}
plot_top_candidates <- function(d, score_col = "priority_score_cont", title = "", top_n = 20, fill_values = cols_candidate) {
  dd <- copy(as.data.table(d))
  if (!nrow(dd) || !score_col %in% names(dd)) return(ggplot() + labs(title = title) + theme_5c())
  dd[, candidate_type2 := fifelse(stress_consequence_flag %in% TRUE, "stress_or_consequence", as.character(candidate_type))]
  aux_p_up <- if ("p_up" %in% names(dd)) dd$p_up else rep(NA_real_, nrow(dd))
  aux_p_low <- if ("p_up_lowSES" %in% names(dd)) dd$p_up_lowSES else rep(NA_real_, nrow(dd))
  aux_p_cad <- if ("p_cad" %in% names(dd)) dd$p_cad else rep(NA_real_, nrow(dd))
  ord <- order(-dd[[score_col]], aux_p_up, aux_p_low, aux_p_cad, na.last = TRUE)
  dd <- dd[ord]
  dd <- dd[seq_len(min(top_n, .N))]
  dd[, label := factor(paste(layer, feature, sep = ":"), levels = rev(paste(layer, feature, sep = ":")))]
  ggplot(dd, aes(.data[[score_col]], label)) +
    geom_segment(aes(x = 0, xend = .data[[score_col]], y = label, yend = label, color = candidate_type2), linewidth = .95) +
    geom_point(aes(fill = candidate_type2), shape = 21, size = 3.0, color = "grey20") +
    scale_color_manual(values = fill_values, guide = "none") +
    scale_fill_manual(values = fill_values, guide = "none") +
    labs(title = title, x = "Continuous evidence score", y = NULL) + theme_5c(9)
}

clean_feature_label <- function(layer, feature = NULL) {
  if (is.null(feature)) {
    x <- as.character(layer)
    x <- sub("^(protein|metabolite):", "", x)
    return(x)
  }
  paste0(ifelse(as.character(layer) == "protein", "P: ", "M: "), as.character(feature))
}
plot_top_candidates_combined <- function(d, top_n = 12, title = "") {
  dd <- copy(as.data.table(d))
  if (!nrow(dd)) return(ggplot() + labs(title = title) + theme_5c())
  dd[, candidate_type2 := fifelse(stress_consequence_flag %in% TRUE, "stress_or_consequence", as.character(candidate_type))]
  res <- dd[stress_consequence_flag %in% FALSE | is.na(stress_consequence_flag)]
  if (nrow(res)) {
    res <- res[order(-priority_score_cont, p_up, p_up_lowSES, p_cad, na.last = TRUE)][seq_len(min(top_n, .N))]
    res[, `:=`(axis = "Resilience / fuel", score = priority_score_cont)]
  }
  str <- dd[stress_score_cont > 0]
  if (nrow(str)) {
    str <- str[order(-stress_score_cont, p_up, p_cad, na.last = TRUE)][seq_len(min(top_n, .N))]
    str[, `:=`(axis = "Stress / consequence", score = stress_score_cont)]
  }
  z <- rbindlist(list(res, str), fill = TRUE)
  if (!nrow(z)) return(ggplot() + labs(title = title) + theme_5c())
  z[, feature_label := clean_feature_label(layer, feature)]
  z[, axis := factor(axis, levels = c("Resilience / fuel", "Stress / consequence"))]
  z[, label_id := paste(axis, feature_label, sep = "|||")]
  z[, label_id := factor(label_id, levels = rev(unique(label_id[order(axis, score)])))]
  ggplot(z, aes(score, label_id, fill = candidate_type2)) +
    geom_col(width = .70, color = "white") +
    facet_wrap(~axis, ncol = 1, scales = "free_y") +
    scale_y_discrete(labels = function(x) sub("^.*\\|\\|\\|", "", x)) +
    scale_fill_manual(values = cols_candidate, guide = "none") +
    labs(title = title, x = "Continuous evidence score", y = NULL) + theme_5c(9) +
    theme(axis.text.y = element_text(face = "bold", size = 7.5))
}

plot_evidence_heatmap <- function(d, top_n = 18, title = "") {
  dd <- copy(as.data.table(d))
  if (!nrow(dd)) return(ggplot() + labs(title = title) + theme_5c())
  dd <- dd[order(-priority_score_cont, p_up, p_up_lowSES, p_cad)]
  dd <- dd[seq_len(min(top_n, .N))]
  dd[, feat_lab := clean_feature_label(layer, feature)]
  dd[, feat_lab := factor(feat_lab, levels = rev(feat_lab))]
  long <- melt(dd[, .(feat_lab, `Upward move` = score_up, `Low-SES
resilience` = score_lowSES,
                     `Opposite move` = score_opposite, `Lower outcome risk` = score_cad_protect,
                     `Higher outcome risk` = score_cad_harm)],
               id.vars = "feat_lab", variable.name = "component", value.name = "score")
  ggplot(long, aes(component, feat_lab, fill = score)) +
    geom_tile(color = "white") +
    scale_fill_gradient2(low = BLUE, mid = "white", high = GOLD, midpoint = 0, name = "Signed
evidence") +
    labs(title = title, x = NULL, y = NULL) + theme_5c(9) + theme(axis.text.y = element_text(size = 8, face = "bold"))
}

save_root_omics_transition_figure <- function(res_list) {
  layers <- c(protein = "Protein", metabolite = "Metabolite")
  plots <- list()
  sheets <- list()
  panel_i <- 0L
  for (layer in names(layers)) {
    obj <- res_list[[layer]]
    if (is.null(obj) || is.null(obj$contrast) || !nrow(obj$contrast)) next
    up <- obj$contrast[contrast == "up_vs_stable_poor"]
    low <- obj$contrast[contrast == "up_lowSES_vs_rr_lowSES"]
    if (!nrow(up) && !nrow(low)) next
    if (nrow(up)) {
      panel_i <- panel_i + 1L
      plots[[length(plots) + 1L]] <- plot_volcano(
        up,
        title = paste0(letters[panel_i], ". ", layers[[layer]], " upward-mobility signature"),
        xlab = "Association with Poor -> rich vs Poor -> poor"
      )
      sheets[[paste0(layer, "_upward")]] <- up[order(p)][seq_len(min(1000, nrow(up)))]
    }
    if (nrow(low)) {
      panel_i <- panel_i + 1L
      plots[[length(plots) + 1L]] <- plot_volcano(
        low,
        title = paste0(letters[panel_i], ". ", layers[[layer]], " low-SES resilience contrast"),
        xlab = "Association with Poor -> rich + low SES vs Rich -> rich + low SES"
      )
      sheets[[paste0(layer, "_lowSES")]] <- low[order(p)][seq_len(min(1000, nrow(low)))]
    }
  }
  if (!length(plots)) return(invisible(NULL))
  fig <- if (length(plots) >= 4L) {
    (plots[[1]] | plots[[2]]) / (plots[[3]] | plots[[4]])
  } else if (length(plots) == 3L) {
    (plots[[1]] | plots[[2]]) / plots[[3]]
  } else if (length(plots) == 2L) {
    plots[[1]] | plots[[2]]
  } else {
    plots[[1]]
  }
  with_output_scope("root", {
    save_plot(fig, fig_file(2, "omics.signature"), width = if (length(plots) >= 4L) 16 else 14, height = if (length(plots) >= 4L) 12 else 6.5)
    write_book(sheets, fig_file(2, "omics.signature", ".out.xlsx"))
  })
  invisible(NULL)
}

run_omics_layer <- function(layer, biom, max_features = 0) {
  if (is.null(biom)) return(NULL)
  core <- read_core()
  core <- as.data.table(core)
  biom <- as.data.table(biom)
  biom[, eid := as.character(eid)]
  xs <- feature_list(biom, max_features)
  if (!length(xs)) return(NULL)

  message("Running ", layer, " omics scans; features=", length(xs), "; N rows=", nrow(biom))
  dat <- merge(core, biom[, c("eid", xs), with = FALSE], by = "eid", all = FALSE, sort = FALSE)
  if (layer == "protein") dat <- dat[has_prot == TRUE]
  if (layer == "metabolite") dat <- dat[has_met == TRUE]

  # Convert all features to numeric once.
  for (x in xs) dat[[x]] <- safe_num(dat[[x]])
  miss <- sapply(xs, function(x) mean(is.na(dat[[x]])))
  keep <- names(miss)[miss <= as.numeric(Sys.getenv("MOVE_OMICS_MAX_MISSING", unset = "0.30"))]
  xs <- keep
  if (!length(xs)) return(NULL)
  message(layer, ": features after missingness filter=", length(xs))

  cov_move <- model_covs(dat, extra = c("source_tdi", "dest_tdi", "individual_ses_z", "move_distance_km", "y_prevalent"))
  cov_move <- setdiff(cov_move, c("center"))
  contrasts <- list(
    up_vs_stable_poor = "contrast_up_vs_stable_poor",
    up_vs_down = "contrast_up_vs_down",
    up_lowSES_vs_rr_lowSES = "contrast_up_lowSES_vs_rr_lowSES"
  )
  contrast_res <- rbindlist(lapply(names(contrasts), function(nm) {
    scan_lm_features(dat, xs, contrasts[[nm]], cov_move, layer, nm, min_n = min_n_model, min_case = min_case_model)
  }), fill = TRUE)
  if (!nrow(contrast_res)) {
    contrast_res <- data.table(layer = character(), contrast = character(), feature = character(), beta = numeric(), se = numeric(), statistic = numeric(), p = numeric(), N = integer(), n_case = integer(), n_ctrl = integer(), model = character(), FDR = numeric())
  }
  if (nrow(contrast_res)) {
    write_raw(contrast_res, paste0(layer, "_transition_signature.tsv.gz"))
  }

  cov_cad <- model_covs(dat, extra = c("source_tdi", "dest_tdi", "individual_ses_z", "tdi_transition"))
  cad_res <- if (cad_scan_do) scan_cox_features(dat, xs, cov_cad, layer, min_n = min_n_model, min_event = min_case_model) else data.table()
  if (nrow(cad_res)) write_raw(cad_res, paste0(layer, "_outcome_signature.tsv.gz"))

  up <- contrast_res[contrast == "up_vs_stable_poor", .(feature, beta_up = beta, p_up = p, FDR_up = FDR, N_up = N, n_case_up = n_case)]
  low <- contrast_res[contrast == "up_lowSES_vs_rr_lowSES", .(feature, beta_up_lowSES = beta, p_up_lowSES = p, FDR_up_lowSES = FDR, N_up_lowSES = N, n_case_up_lowSES = n_case)]
  down <- contrast_res[contrast == "up_vs_down", .(feature, beta_up_vs_down = beta, p_up_vs_down = p, FDR_up_vs_down = FDR)]
  if (!nrow(up)) up <- data.table(feature = xs, beta_up = NA_real_, p_up = NA_real_, FDR_up = NA_real_, N_up = NA_integer_, n_case_up = NA_integer_)
  if (!nrow(low)) low <- data.table(feature = xs, beta_up_lowSES = NA_real_, p_up_lowSES = NA_real_, FDR_up_lowSES = NA_real_, N_up_lowSES = NA_integer_, n_case_up_lowSES = NA_integer_)
  if (!nrow(down)) down <- data.table(feature = xs, beta_up_vs_down = NA_real_, p_up_vs_down = NA_real_, FDR_up_vs_down = NA_real_)
  if (!nrow(cad_res)) cad_res <- data.table(feature = xs, beta_cad = NA_real_, se_cad = NA_real_, p_cad = NA_real_, HR_cad = NA_real_, lo_cad = NA_real_, hi_cad = NA_real_, N_cad = NA_integer_, events_cad = NA_integer_, FDR_cad = NA_real_)
  driver <- Reduce(function(a, b) merge(a, b, by = "feature", all = TRUE, sort = FALSE), list(up, low, down, cad_res))
  if (!nrow(driver)) driver <- data.table(feature = xs)
  driver[, layer := layer]
  driver[, driver_score := 0L]
  driver[is.finite(FDR_up) & FDR_up < 0.05 & beta_up > 0, driver_score := driver_score + 1L]
  driver[is.finite(FDR_up_lowSES) & FDR_up_lowSES < 0.05 & beta_up_lowSES > 0, driver_score := driver_score + 1L]
  driver[is.finite(p_up_vs_down) & p_up_vs_down < 0.05 & beta_up_vs_down > 0, driver_score := driver_score + 1L]
  driver[is.finite(p_cad) & p_cad < 0.05 & beta_cad < 0, driver_score := driver_score + 1L]
  driver[, candidate_type := fifelse(driver_score >= 3, "high_priority_resilience_driver",
                                     fifelse(driver_score == 2, "medium_priority_driver", "exploratory"))]
  driver[, score_up := score_signed_logp(beta_up, p_up, positive_good = TRUE)]
  driver[, score_lowSES := score_signed_logp(beta_up_lowSES, p_up_lowSES, positive_good = TRUE)]
  driver[, score_opposite := score_signed_logp(beta_up_vs_down, p_up_vs_down, positive_good = TRUE)]
  driver[, score_cad_protect := score_signed_logp(beta_cad, p_cad, positive_good = FALSE)]
  driver[, score_cad_harm := score_signed_logp(beta_cad, p_cad, positive_good = TRUE)]
  driver[, priority_score_cont := score_pos(score_up) + score_pos(score_lowSES) + score_pos(score_opposite) + score_pos(score_cad_protect)]
  driver[, stress_score_cont := score_pos(score_up) + score_pos(score_cad_harm)]
  # A separate stress/consequence flag: upward mover signature but higher disease risk.
  driver[, stress_consequence_flag := is.finite(FDR_up) & FDR_up < 0.05 & beta_up > 0 & is.finite(p_cad) & p_cad < 0.05 & beta_cad > 0]
  setorder(driver, -driver_score, p_up, p_up_lowSES, p_cad)
  write_raw(driver, paste0(layer, "_driver_candidates.tsv.gz"))
  save_rds_out(list(contrast = contrast_res, cad = cad_res, driver = driver), paste0(layer, "_omics_signature.res.rds"), compress = "xz")

  # Transition panels are root-level; disease folders keep outcome and driver panels only.
  pA <- plot_volcano(contrast_res[contrast == "up_vs_stable_poor"],
                     title = "a. Upward-mobility signature",
                     xlab = "Association with Poor → rich vs Poor → poor")
  pB <- plot_volcano(contrast_res[contrast == "up_lowSES_vs_rr_lowSES"],
                     title = "b. Low-SES resilience contrast",
                     xlab = "Association with Poor → rich + low SES vs Rich → rich + low SES")
  pC <- if (nrow(cad_res)) plot_volcano(cad_res[, .(feature, beta = beta_cad, p = p_cad, FDR = FDR_cad)],
                                        title = "a. Incident outcome signature",
                                        xlab = "Cox log-HR for outcome") else ggplot() + labs(title = "a. Outcome Cox scan skipped") + theme_5c()
  pD <- plot_driver_heatmap(driver, layer, top_n = 18) + labs(title = paste0("b. ", layer, " prioritized driver evidence"))
  Fig <- pC | pD
  layer_tag <- if (layer == "protein") "prot" else if (layer == "metabolite") "met" else layer
  fig_name <- fig_file(4, paste(layer_tag, "signature", sep = "."), trait = TRUE)
  out_xlsx <- fig_file(4, paste(layer_tag, "signature", sep = "."), ".out.xlsx", trait = TRUE)
  save_plot(Fig, fig_name, 16, 6.5)

  top_driver <- driver[seq_len(min(200, nrow(driver)))]
  write_book(list(
    driver_candidates_top200 = top_driver,
    transition_signature_top = contrast_res[order(p)][seq_len(min(1000, nrow(contrast_res)))],
    outcome_signature_top = if (nrow(cad_res)) cad_res[order(p_cad)][seq_len(min(1000, nrow(cad_res)))] else data.table(note = "Outcome scan skipped"),
    notes = data.frame(note = c(
      "Y-independent upward-mobility and low-SES resilience panels are written once at the root as Fig2.omics.signature.png.",
      "Panel b shows component-level evidence instead of a score ending at 1 or 2, so each candidate can be interpreted across move, low-SES resilience, opposite-move, and outcome-risk axes.",
      "These are discovery signatures, not proof of causality. For causal evidence, use pQTL/mQTL instruments and an independent move GWAS in follow-up MR/coloc."
    ))
  ), out_xlsx)
  list(contrast = contrast_res, cad = cad_res, driver = driver)
}

read_existing_omics_layer <- function(layer) {
  rds_file <- file.path(rawdir, paste0(layer, "_omics_signature.res.rds"))
  if (!file.exists(rds_file)) return(NULL)
  res <- tryCatch(readRDS(rds_file), error = function(e) {
    message("Could not read existing ", layer, " omics result: ", conditionMessage(e))
    NULL
  })
  if (is.null(res) || !is.list(res) || !"driver" %in% names(res)) return(NULL)
  message("Reusing existing ", layer, " omics result: ", rds_file)
  res
}

run_if("omics_signature", {
  set_output_scope("trait")
  res_list <- list()
  if (prot_do) {
    res_list$protein <- if (resume_omics) read_existing_omics_layer("protein") else NULL
    if (is.null(res_list$protein)) {
      prot <- read_layer("prot.rds", required = FALSE)
      if (!is.null(prot)) res_list$protein <- run_omics_layer("protein", prot, max_features = max_prot_features)
    }
  }
  if (met_do) {
    res_list$metabolite <- if (resume_omics) read_existing_omics_layer("metabolite") else NULL
    if (is.null(res_list$metabolite)) {
      met <- read_layer("met.rds", required = FALSE)
      if (!is.null(met)) res_list$metabolite <- run_omics_layer("metabolite", met, max_features = max_met_features)
    }
  }
  save_root_omics_transition_figure(res_list)
  save_rds_out(res_list, "omics_signature.res.rds", compress = "xz")

  # Combined cross-layer integration figure.
  drivers <- rbindlist(lapply(res_list, function(z) if (!is.null(z)) z$driver else data.table()), fill = TRUE)
  if (nrow(drivers)) {
    drivers[, label := paste(layer, feature, sep = ":")]
    drivers[, candidate_type2 := fifelse(stress_consequence_flag %in% TRUE, "stress_or_consequence", as.character(candidate_type))]
    drivers[, layer := factor(layer, levels = c("protein", "metabolite"))]
    score_long <- melt(drivers[is.finite(priority_score_cont) | is.finite(stress_score_cont)],
                       id.vars = c("layer", "feature"),
                       measure.vars = c("priority_score_cont", "stress_score_cont"),
                       variable.name = "score_type", value.name = "score")
    score_long <- score_long[is.finite(score)]
    score_long[, score_type := fifelse(as.character(score_type) == "priority_score_cont", "Priority / resilience", "Stress / consequence")]
    score_long[, score_type := factor(score_type, levels = c("Priority / resilience", "Stress / consequence"))]
    p4a <- ggplot(score_long, aes(layer, score, fill = layer)) +
      geom_violin(alpha = .30, color = NA) + geom_boxplot(width = .16, outlier.size = .25, alpha = .90) +
      facet_wrap(~score_type, nrow = 1, scales = "free_y") +
      scale_fill_manual(values = c(protein = GOLD, metabolite = BLUE), guide = "none") +
      labs(title = "a. Priority and stress score spectra", x = NULL, y = "Score") + theme_5c(10)

    lab2 <- drivers[is.finite(beta_up) & is.finite(beta_cad)][order(-priority_score_cont)][seq_len(min(18, .N))]
    p4b <- ggplot(drivers[is.finite(beta_up) & is.finite(beta_cad)], aes(beta_up, beta_cad)) +
      geom_hline(yintercept = 0, color = "grey70") + geom_vline(xintercept = 0, color = "grey70") +
      geom_point(aes(color = layer, shape = candidate_type2), alpha = .72, size = 1.6) +
      ggrepel::geom_text_repel(data = lab2, aes(label = feature), size = 2.6, max.overlaps = 20, seed = 12, min.segment.length = 0) +
      scale_color_manual(values = c(protein = GOLD, metabolite = BLUE), name = "Layer") +
      labs(title = "b. Upward mobility vs outcome risk", x = "Association with Poor → rich", y = "Association with outcome risk", shape = NULL) + theme_5c(10)

    lab3 <- drivers[is.finite(beta_up_lowSES) & is.finite(beta_cad)][order(-priority_score_cont)][seq_len(min(18, .N))]
    p4c <- ggplot(drivers[is.finite(beta_up_lowSES) & is.finite(beta_cad)], aes(beta_up_lowSES, beta_cad)) +
      geom_hline(yintercept = 0, color = "grey70") + geom_vline(xintercept = 0, color = "grey70") +
      geom_point(aes(color = layer), alpha = .72, size = 1.6) +
      ggrepel::geom_text_repel(data = lab3, aes(label = feature), size = 2.6, max.overlaps = 18, seed = 13, min.segment.length = 0) +
      scale_color_manual(values = c(protein = GOLD, metabolite = BLUE), guide = "none") +
      labs(title = "c. Low-SES resilience vs outcome risk", x = "Association with Poor → rich + low SES", y = "Association with outcome risk") + theme_5c(10)

    lab4 <- drivers[is.finite(beta_up) & is.finite(beta_up_lowSES)][order(-priority_score_cont)][seq_len(min(18, .N))]
    p4d <- ggplot(drivers[is.finite(beta_up) & is.finite(beta_up_lowSES)], aes(beta_up, beta_up_lowSES)) +
      geom_hline(yintercept = 0, color = "grey70") + geom_vline(xintercept = 0, color = "grey70") +
      geom_point(aes(color = layer), alpha = .72, size = 1.6) +
      ggrepel::geom_text_repel(data = lab4, aes(label = feature), size = 2.6, max.overlaps = 18, seed = 14, min.segment.length = 0) +
      scale_color_manual(values = c(protein = GOLD, metabolite = BLUE), guide = "none") +
      labs(title = "d. Upward mobility vs resilience contrast", x = "Association with Poor → rich", y = "Association with Poor → rich + low SES") + theme_5c(10)

    p4e <- plot_top_candidates_combined(drivers, top_n = 14, title = "e. Top resilience and stress candidates")
    p4f <- plot_evidence_heatmap(drivers, top_n = 16, title = "f. Component-level evidence map")

    Fig6 <- (p4a | p4b) / (p4c | p4d) / (p4e | p4f)
    save_plot(Fig6, fig_file(6, "cross_layer_driver", trait = TRUE), 18, 15)
    write_book(list(cross_layer_driver_candidates = drivers[seq_len(min(1000, nrow(drivers)))],
                    top_resilience_candidates = drivers[order(-priority_score_cont)][seq_len(min(200, nrow(drivers)))],
                    top_stress_candidates = drivers[order(-stress_score_cont)][seq_len(min(200, nrow(drivers)))]) , fig_file(6, "cross_layer_driver", ".out.xlsx", trait = TRUE))
  }
})


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 🚩 Step 5: biomarker-score and attenuation analysis
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

read_raw_if_exists <- function(fn) {
  f <- file.path(rawdir, fn)
  if (file.exists(f) && file.size(f) > 0) data.table::fread(f) else data.table()
}
metab_class <- function(x) {
  x <- as.character(x)
  fifelse(grepl("Pyruvate|Lactate|Citrate|Acetate", x, ignore.case = TRUE), "Energy / organic acids",
  fifelse(grepl("HDL|LDL|VLDL|IDL|TG|CE|PL|MUFA|PUFA|LA|Apo", x, ignore.case = TRUE), "Lipoprotein / fatty acids",
  fifelse(grepl("Gly|Ala|His|Ile|Leu|Phe|Tyr|Val|Gln|Glu|Creatinine", x, ignore.case = TRUE), "Amino acid / renal",
         "Other metabolites")))
}
make_signature_score <- function(layer, biom, top_n = 20) {
  sig <- read_raw_if_exists(paste0(layer, "_transition_signature.tsv.gz"))
  sig <- sig[contrast == "up_vs_stable_poor" & is.finite(beta) & is.finite(p)]
  if (!nrow(sig) || is.null(biom)) return(NULL)
  sig[, pass := is.finite(FDR) & FDR < 0.05]
  sig <- sig[order(!pass, p)][seq_len(min(top_n, .N))]
  biom <- as.data.table(biom); biom[, eid := as.character(eid)]
  vars <- intersect(sig$feature, names(biom))
  if (!length(vars)) return(NULL)
  sig <- sig[feature %in% vars]
  signs <- setNames(sign(sig$beta), sig$feature)
  dat <- biom[, c("eid", vars), with = FALSE]
  for (v in vars) dat[[v]] <- inormal2(dat[[v]]) * signs[[v]]
  dat[, score := rowMeans(.SD, na.rm = TRUE), .SDcols = vars]
  dat[, n_score_features := rowSums(!is.na(.SD)), .SDcols = vars]
  dat[n_score_features < max(3, ceiling(length(vars) * 0.5)), score := NA_real_]
  list(score = dat[, .(eid, score, n_score_features)], signature = sig, vars = vars)
}
score_cox_models <- function(dat, score_vars) {
  out <- list(); k <- 0L
  base <- dat[!is.na(y_t2e) & !is.na(y_event) & y_t2e > 0 & y_prevalent != 1]
  for (sv in score_vars) {
    if (!sv %in% names(base)) next
    dd <- base[is.finite(get(sv))]
    if (nrow(dd) < min_n_model || sum(dd$y_event == 1, na.rm = TRUE) < min_case_model) next
    dd[, score_z := zscore(get(sv))]
    model_list <- list(
      base = model_covs(dd, extra = c("source_tdi", "dest_tdi", "individual_ses_z")),
      transition_SES = model_covs(dd, extra = c("source_tdi", "dest_tdi", "individual_ses_z", "tdi_transition", "ses_extreme"))
    )
    for (mn in names(model_list)) {
      covs <- model_list[[mn]]
      rhs <- c("score_z", covs)
      f <- as.formula(paste0("survival::Surv(y_t2e, y_event) ~ ", paste(bt(rhs), collapse = " + ")))
      fit <- tryCatch(survival::coxph(f, data = dd), error = function(e) NULL)
      if (is.null(fit)) next
      z <- tidy_cox_term(fit)[term == "score_z"]
      if (nrow(z)) { k <- k + 1L; z[, `:=`(score = sv, model = mn, N = nrow(dd), events = sum(dd$y_event == 1, na.rm = TRUE), label = paste0(sv, "\n", mn))]; out[[k]] <- z }
    }
  }
  rbindlist(out, fill = TRUE)
}
attenuation_by_score <- function(dat, score_var) {
  if (!score_var %in% names(dat)) return(data.table())
  d <- dat[!is.na(y_t2e) & !is.na(y_event) & y_t2e > 0 & y_prevalent != 1 & !is.na(tdi_transition) & is.finite(get(score_var))]
  if (nrow(d) < min_n_model || sum(d$y_event == 1, na.rm = TRUE) < min_case_model) return(data.table())
  d[, tdi_transition := relevel(factor(as.character(tdi_transition), levels = names(cols_transition)), ref = "poor -> poor")]
  d[, score_z := zscore(get(score_var))]
  covs <- model_covs(d, extra = c("source_tdi", "dest_tdi", "individual_ses_z"))
  f0 <- as.formula(paste0("survival::Surv(y_t2e, y_event) ~ tdi_transition", if (length(covs)) paste0(" + ", paste(bt(covs), collapse = " + ")) else ""))
  f1 <- as.formula(paste0("survival::Surv(y_t2e, y_event) ~ tdi_transition + score_z", if (length(covs)) paste0(" + ", paste(bt(covs), collapse = " + ")) else ""))
  fit0 <- tryCatch(survival::coxph(f0, data = d), error = function(e) NULL)
  fit1 <- tryCatch(survival::coxph(f1, data = d), error = function(e) NULL)
  if (is.null(fit0) || is.null(fit1)) return(data.table())
  b0 <- tidy_cox_term(fit0)[grepl("^tdi_transition", term), .(term, beta0 = beta, HR0 = HR, p0 = p)]
  b1 <- tidy_cox_term(fit1)[grepl("^tdi_transition", term), .(term, beta1 = beta, HR1 = HR, p1 = p)]
  z <- merge(b0, b1, by = "term")
  z[, attenuation_pct := 100 * (beta0 - beta1) / beta0]
  z[, `:=`(score = score_var, label = clean_cox_label(term))]
  z
}
plot_score_by_transition <- function(dat, score_var, title) {
  d <- dat[!is.na(tdi_transition) & is.finite(get(score_var))]
  if (!nrow(d)) return(ggplot() + labs(title = title) + theme_5c())
  d[, score_z := zscore(get(score_var))]
  ggplot(d, aes(short_transition(tdi_transition), score_z, fill = tdi_transition)) +
    geom_violin(alpha = .35, color = NA, trim = TRUE) +
    geom_boxplot(width = .14, outlier.size = .12, alpha = .9) +
    scale_fill_manual(values = cols_transition, guide = "none") +
    labs(title = title, x = "TDI transition", y = "Signature score\nz-score") + theme_5c(10)
}
plot_score_heatmap <- function(dat, sig, layer, top_n = 16) {
  if (is.null(sig) || !nrow(sig)) return(ggplot() + labs(title = "i. Feature heatmap unavailable") + theme_5c())
  vars <- sig[order(p)][seq_len(min(top_n, .N))]$feature
  vars <- intersect(vars, names(dat))
  if (!length(vars)) return(ggplot() + labs(title = "i. Feature heatmap unavailable") + theme_5c())
  dd <- dat[!is.na(tdi_transition), c("tdi_transition", vars), with = FALSE]
  for (v in vars) dd[[v]] <- zscore(dd[[v]])
  long <- melt(dd, id.vars = "tdi_transition", measure.vars = vars, variable.name = "feature", value.name = "z")
  hm <- long[, .(mean_z = mean(z, na.rm = TRUE)), by = .(tdi_transition, feature)]
  hm[, tdi_transition := factor(short_transition(tdi_transition), levels = short_transition(names(cols_transition)))]
  hm[, feature := factor(feature, levels = rev(vars))]
  ggplot(hm, aes(tdi_transition, feature, fill = mean_z)) +
    geom_tile(color = "white") +
    scale_fill_gradient2(low = BLUE, mid = "white", high = GOLD, midpoint = 0, name = "Mean z") +
    labs(title = paste0("i. Top ", layer, " features by transition"), x = NULL, y = NULL) + theme_5c(9) + theme(axis.text.y = element_text(size = 7, face = "bold"))
}

run_if("deep_dive", {
  set_output_scope("trait")
  core <- as.data.table(read_core())
  core[, eid := as.character(eid)]
  met_obj <- if (met_do) make_signature_score("metabolite", read_layer("met.rds", required = FALSE), top_n = as.integer(Sys.getenv("MOVE_SCORE_TOP_N", unset = "20"))) else NULL
  prot_obj <- if (prot_do) make_signature_score("protein", read_layer("prot.rds", required = FALSE), top_n = as.integer(Sys.getenv("MOVE_SCORE_TOP_N", unset = "20"))) else NULL
  score_dat <- copy(core)
  if (!is.null(met_obj)) { setnames(met_obj$score, "score", "met_mobility_score"); score_dat <- merge(score_dat, met_obj$score, by = "eid", all.x = TRUE) }
  if (!is.null(prot_obj)) { setnames(prot_obj$score, "score", "prot_mobility_score"); score_dat <- merge(score_dat, prot_obj$score, by = "eid", all.x = TRUE) }
  score_vars <- intersect(c("met_mobility_score", "prot_mobility_score"), names(score_dat))
  score_cox <- score_cox_models(score_dat, score_vars)
  atten <- rbindlist(lapply(score_vars, function(sv) attenuation_by_score(score_dat, sv)), fill = TRUE)
  write_raw(score_cox, "mobility_score_outcome_cox.tsv.gz")
  write_raw(atten, "mobility_score_transition_attenuation.tsv.gz")
  save_rds_out(list(score_cox = score_cox, attenuation = atten,
                    met_signature = if (!is.null(met_obj)) met_obj$signature else data.table(),
                    prot_signature = if (!is.null(prot_obj)) prot_obj$signature else data.table()),
               "biomarker_score.res.rds", compress = "xz")

  pA <- if (!is.null(met_obj)) {
    d <- copy(met_obj$signature); d[, class := metab_class(feature)]; d <- d[order(p)][seq_len(min(25, .N))]
    d[, feature := factor(feature, levels = rev(feature))]
    ggplot(d, aes(beta, feature, fill = class)) + geom_col(width = .72, color = "white") +
      geom_vline(xintercept = 0, color = "grey60") + labs(title = "a. Top upward-mobility metabolites", x = "Association with Poor → rich", y = NULL, fill = "Class") + theme_5c(9)
  } else ggplot() + labs(title = "a. No metabolite score") + theme_5c()
  pB <- if (!is.null(prot_obj)) {
    d <- copy(prot_obj$signature); d <- d[order(p)][seq_len(min(25, .N))]; d[, feature := factor(feature, levels = rev(feature))]
    ggplot(d, aes(beta, feature, fill = beta > 0)) + geom_col(width = .72, color = "white") +
      geom_vline(xintercept = 0, color = "grey60") + scale_fill_manual(values = c(`TRUE` = GOLD, `FALSE` = BLUE), guide = "none") +
      labs(title = "b. Top upward-mobility proteins", x = "Association with Poor → rich", y = NULL) + theme_5c(9)
  } else ggplot() + labs(title = "b. No protein score") + theme_5c()
  pC <- if (!is.null(met_obj)) plot_score_by_transition(score_dat, "met_mobility_score", "c. Metabolite mobility score by transition") else ggplot() + labs(title = "c. No metabolite score") + theme_5c()
  pD <- if (!is.null(prot_obj)) plot_score_by_transition(score_dat, "prot_mobility_score", "d. Protein mobility score by transition") else ggplot() + labs(title = "d. No protein score") + theme_5c()
  pE <- if (nrow(score_cox)) plot_forest_clean(score_cox[, .(HR, lo, hi, p, label)], title = "e. Mobility-score association with outcome risk", ref_label = "Per 1 s.d. higher score", base_size = 10) + theme(axis.text.y = element_text(face = "bold", size = 8.5)) else ggplot() + labs(title = "e. Score-outcome model unavailable") + theme_5c()
  pF <- if (nrow(atten)) {
    a <- atten[is.finite(attenuation_pct)]; a[, score := gsub("_mobility_score", "", score)]; a[, label2 := paste0(score, "\n", clean_cox_label(term))]
    ggplot(a, aes(attenuation_pct / 100, reorder(label2, attenuation_pct), fill = score)) +
      geom_vline(xintercept = 0, color = "grey55", linetype = "dashed") + geom_col(width = .72, color = "white") +
      scale_x_continuous(labels = percent) + scale_fill_manual(values = c(met = BLUE, prot = GOLD, metabolite = BLUE, protein = GOLD), guide = "none") +
      labs(title = "f. Transition-HR attenuation after adding score", x = "Beta attenuation", y = NULL) + theme_5c(9)
  } else ggplot() + labs(title = "f. Attenuation unavailable") + theme_5c()
  pG <- if (!is.null(met_obj)) {
    sig <- copy(met_obj$signature); sig[, class := metab_class(feature)]; sig[, direction := ifelse(beta > 0, "Higher in Poor→rich", "Lower in Poor→rich")]
    tab <- sig[, .N, by = .(class, direction)]
    ggplot(tab, aes(N, forcats::fct_reorder(class, N, sum), fill = direction)) + geom_col(color = "white") +
      scale_fill_manual(values = c(`Higher in Poor→rich` = GOLD, `Lower in Poor→rich` = BLUE)) +
      labs(title = "g. Metabolite class composition", x = "Top features", y = NULL, fill = NULL) + theme_5c(9)
  } else ggplot() + labs(title = "g. Class composition unavailable") + theme_5c()
  pH <- if (all(c("met_mobility_score", "tdi_improve") %in% names(score_dat))) {
    d <- score_dat[is.finite(met_mobility_score) & is.finite(tdi_improve)]
    qs <- unique(quantile(d$tdi_improve, probs = seq(0, 1, 0.25), na.rm = TRUE))
    if (length(qs) < 3) d[, improve_q := factor(NA_character_)] else d[, improve_q := cut(tdi_improve, breaks = qs, include.lowest = TRUE, labels = paste0("Q", seq_len(length(qs) - 1)))]
    d[, met_score_z := zscore(met_mobility_score)]
    ggplot(d[!is.na(improve_q)], aes(improve_q, met_score_z, fill = improve_q)) + geom_boxplot(outlier.size = .12, show.legend = FALSE) +
      labs(title = "h. Metabolite score by TDI-improvement quartile", x = "TDI-improvement quartile", y = "Score z") + theme_5c(9)
  } else ggplot() + labs(title = "h. Score vs TDI improvement unavailable") + theme_5c()
  # Heatmap needs biomarker columns; merge only top met features for memory.
  pI <- if (!is.null(met_obj)) {
    met <- as.data.table(read_layer("met.rds", required = FALSE)); met[, eid := as.character(eid)]
    vars <- intersect(met_obj$vars, names(met)); tmp <- merge(score_dat[, .(eid, tdi_transition)], met[, c("eid", vars), with = FALSE], by = "eid", all = FALSE)
    plot_score_heatmap(tmp, met_obj$signature, "metabolite", top_n = 16)
  } else ggplot() + labs(title = "i. Heatmap unavailable") + theme_5c()
  Fig5 <- (pA | pB | pC) / (pD | pE | pF) / (pG | pH | pI)
  save_plot(Fig5, fig_file(5, "biomarker_score", trait = TRUE), 18, 16)
  write_book(list(score_outcome_cox = score_cox, attenuation = atten,
                  met_score_signature = if (!is.null(met_obj)) met_obj$signature else data.table(),
                  prot_score_signature = if (!is.null(prot_obj)) prot_obj$signature else data.table()), fig_file(5, "biomarker_score", ".out.xlsx", trait = TRUE))
})

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 🚩 Step 6: consolidate outputs
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

run_if("consolidate", {
  set_output_scope("trait")
  core <- read_core()
  files <- data.table(file = list.files(outdir, recursive = TRUE, full.names = FALSE))
  files[, size_MB := round(file.info(file.path(outdir, file))$size / 1024^2, 3)]
  files[, type := fifelse(grepl("\\.png$", file), "figure", fifelse(grepl("\\.xlsx$", file), "excel", fifelse(grepl("\\.rds$", file), "rds", "other")))]

  transition_summary <- as.data.table(core)[, .(
    N = .N,
    has_prot = sum(has_prot),
    has_met = sum(has_met),
    low_ses_pct = mean(person_ses_class == "low", na.rm = TRUE),
    incident_events = sum(y_event == 1, na.rm = TRUE),
    incident_pct = mean(y_event == 1, na.rm = TRUE)
  ), by = tdi_transition][order(tdi_transition)]

  driver_files <- list.files(rawdir, pattern = "_driver_candidates.tsv.gz$", full.names = TRUE)
  drivers <- rbindlist(lapply(driver_files, data.table::fread), fill = TRUE)
  if (nrow(drivers)) {
    setorder(drivers, -driver_score, p_up, p_up_lowSES, p_cad)
  } else {
    drivers <- data.table(note = "No omics driver files found. Run START_STEP=omics_signature.")
  }

  run_notes <- data.frame(
    item = c("Main hypothesis", "Main transition", "Rich/poor definition", "Disease outcome", "Causality caution", "Next causal step"),
    value = c(
      "TDI-stratified social mobility has SES, outcome-risk, and omics signatures; poor -> rich movers with low current SES may show resilience/driver biomarkers.",
      main_transition,
      paste0("rich = TDI lowest ", rich_poor_q * 100, "%; poor = TDI highest ", rich_poor_q * 100, "%"),
      Y,
      "Current UKB proteomics/metabolomics are baseline biomarkers. Birth/current transition is retrospective; home_future transition is prospective but sparse. Observational signatures are not causal proof.",
      "Use top driver candidates to build pQTL/mQTL instrument manifests and test biomarker -> move / biomarker -> outcome with MR and coloc."
    ), stringsAsFactors = FALSE
  )

  write_raw_book(list(
    run_notes = run_notes,
    transition_summary = transition_summary,
    top_driver_candidates = if (nrow(drivers) && "driver_score" %in% names(drivers)) drivers[seq_len(min(500, nrow(drivers)))] else drivers,
    output_index = files[order(type, file)]
  ), "move_v6_project_summary.xlsx")
  save_rds_out(list(notes = run_notes, transition_summary = transition_summary, drivers = drivers, files = files), "move_v6_project_summary.rds", compress = "xz")
  sync_root_rds_to_raw()

  cat("\nKey output files:\n")
  cat("  raw/move_core.rds\n")
  cat("  Fig1.job_home_TDI_summary.png\n")
  cat("  ", file.path(outroot, fig_file(2, "omics.signature")), "\n", sep = "")
  cat("  ", fig_file(3, "SES_risk", trait = TRUE), "\n", sep = "")
  cat("  ", fig_file(4, "prot.signature", trait = TRUE), " / ", fig_file(4, "met.signature", trait = TRUE), "\n", sep = "")
  cat("  ", fig_file(5, "biomarker_score", trait = TRUE), "\n", sep = "")
  cat("  ", fig_file(6, "cross_layer_driver", trait = TRUE), "\n", sep = "")
  cat("  ", fig_file(7, "paradox", trait = TRUE), "\n", sep = "")
  cat("  raw/move_v6_project_summary.xlsx\n")
  cat("\nFinished successfully.\n")
})
