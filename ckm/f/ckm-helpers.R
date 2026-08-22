#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Generic helpers
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

as_date2 <- function(x) {
  if (inherits(x, "Date")) return(x)
  if (inherits(x, "IDate")) return(as.Date(x))
  if (inherits(x, c("POSIXct", "POSIXt"))) return(as.Date(x))
  if (is.character(x)) return(suppressWarnings(as.Date(x)))
  y <- suppressWarnings(as.numeric(unclass(x)))
  y[!is.finite(y)] <- NA_real_
  as.Date(y, origin = "1970-01-01")
}

pmin_date2 <- function(...) {
  x <- list(...)
  if (length(x) == 1 && is.data.frame(x[[1]])) x <- as.list(x[[1]])
  x <- Filter(Negate(is.null), x)
  if (!length(x)) return(as.Date(character()))
  x <- lapply(x, as_date2)
  out <- do.call(pmin, c(x, list(na.rm = TRUE)))
  out <- as_date2(out)
  out[!is.finite(as.numeric(out))] <- as.Date(NA)
  out
}

first_present <- function(dat, candidates) {
  z <- intersect(candidates, names(dat))
  if (length(z)) z[1] else NA_character_
}

safe_num <- function(x) suppressWarnings(as.numeric(as.character(x)))
bt <- function(x) paste0("`", gsub("`", "", x), "`")

safe_covars <- function(dat, covars) {
  covars <- unique(intersect(covars, names(dat)))
  covars[vapply(covars, function(v) {
    z <- dat[[v]]
    sum(!is.na(z)) > 20 && length(unique(z[!is.na(z)])) > 1
  }, logical(1))]
}

get_demog_covars <- function(dat) {
  x <- ckm_vars_basic
  miss <- setdiff(x, names(dat))
  if (length(miss)) stop("Variables listed in vars.basic are missing: ", paste(miss, collapse = ", "), call. = FALSE)
  safe_covars(dat, x)
}

get_le8_covars <- function(dat) {
  miss <- setdiff(ckm_vars_le8, names(dat))
  if (length(miss)) stop("Variables listed in vars.le8 are missing: ", paste(miss, collapse = ", "), call. = FALSE)
  safe_covars(dat, ckm_vars_le8)
}

get_base_covars <- function(dat, mode = pred_base_mode) {
  demog <- get_demog_covars(dat)
  behavior <- safe_covars(dat, clock_behavior_vars)
  le8 <- get_le8_covars(dat)
  genetic <- grep("prs|polygenic|apoe|rare|lof|plof", names(dat), value = TRUE, ignore.case = TRUE)
  # Primary prediction uses demographics + behavioral LE8 only.  Full LE8 is a
  # sensitivity model because BMI/BP/HbA1c/non-HDL partly define CKM stage.
  out <- switch(
    mode,
    demog = demog,
    behavior = c(demog, behavior),
    le8 = c(demog, le8),
    full = c(demog, le8, genetic),
    c(demog, behavior)
  )
  safe_covars(dat, out)
}

get_clock_covariate_sets <- function(dat) {
  demog <- get_demog_covars(dat)
  behavior <- safe_covars(dat, clock_behavior_vars)
  list(
    crude = character(),
    demog = demog,
    le8_behavior = unique(c(demog, behavior)),
    # le8_full is retained as an explicitly overadjusted sensitivity model
    # because BMI/BP/HbA1c/non-HDL help define CKM stages 1-3.
    le8_full = unique(c(demog, get_le8_covars(dat)))
  )
}

get_clock_covars <- function(dat, mode = clock_adjust_mode) {
  sets <- get_clock_covariate_sets(dat)
  if (!mode %in% names(sets)) mode <- "le8_behavior"
  safe_covars(dat, sets[[mode]])
}

winsorize <- function(x, p = c(0.001, 0.999)) {
  x <- safe_num(x)
  ok <- is.finite(x)
  if (sum(ok) < 20) return(x)
  q <- quantile(x[ok], p, na.rm = TRUE, names = FALSE, type = 8)
  pmin(pmax(x, q[1]), q[2])
}

standardize_ref <- function(x, ref) {
  x <- winsorize(x)
  ref <- ref %in% TRUE & is.finite(x)
  if (sum(ref) < 100) ref <- is.finite(x)
  mu <- mean(x[ref], na.rm = TRUE)
  s <- sd(x[ref], na.rm = TRUE)
  if (!is.finite(s) || s == 0) return(rep(NA_real_, length(x)))
  (x - mu) / s
}

standardize_ref_train <- function(x, ref, train_scope, p = c(0.001, 0.999)) {
  x <- safe_num(x)
  train_scope <- train_scope %in% TRUE & is.finite(x)
  if (sum(train_scope) >= 20) {
    q <- quantile(x[train_scope], p, na.rm = TRUE, names = FALSE, type = 8)
    x <- pmin(pmax(x, q[1]), q[2])
  }
  ref <- ref %in% TRUE & train_scope & is.finite(x)
  if (sum(ref) < 100) ref <- train_scope
  mu <- mean(x[ref], na.rm = TRUE)
  ss <- sd(x[ref], na.rm = TRUE)
  if (!is.finite(ss) || ss == 0) return(rep(NA_real_, length(x)))
  (x - mu) / ss
}

write_tab <- function(x, filename, dir = rawdir) {
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  path <- file.path(dir, filename)
  data.table::fwrite(as.data.table(x), path)
  progress_log("Wrote: ", path)
  invisible(path)
}

write_book <- function(x, filename, dir = outroot) {
  x <- x[vapply(x, function(z) is.data.frame(z) || data.table::is.data.table(z), logical(1))]
  if (!length(x)) x <- list(empty = data.frame(note = "No rows"))
  names(x) <- substr(make.unique(gsub("[^A-Za-z0-9_]+", "_", names(x))), 1, 31)
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  path <- file.path(dir, filename)
  writexl::write_xlsx(lapply(x, as.data.frame), path)
  progress_log("Wrote: ", path)
  invisible(path)
}

save_plot <- function(p, filename, width = 12, height = 8, dir = outroot) {
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  path <- file.path(dir, filename)
  ggplot2::ggsave(path, p, width = width, height = height, dpi = 300, bg = "white", limitsize = FALSE)
  progress_log("Wrote: ", path)
  invisible(path)
}

biom_fig_filename <- function(filename) {
  if (!grepl("^Fig[^.]+\\.", filename)) return(filename)
  label <- gsub("[^A-Za-z0-9]+", "_", biom)
  label <- gsub("^_+|_+$", "", label)
  if (!nzchar(label)) label <- biom
  if (grepl(paste0("^Fig[^.]+\\.", label, "\\."), filename)) return(filename)
  sub("^([^.]+)\\.", paste0("\\1.", label, "."), filename)
}

save_rds_safe <- function(object, path, compress = "gzip") {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  t0 <- Sys.time()
  progress_log("Saving lightweight RDS: ", path, "; compression=", as.character(compress))
  saveRDS(object, path, compress = compress)
  progress_log("Wrote: ", path, " in ", elapsed_since(t0))
  invisible(path)
}

theme_ckm <- function(base_size = 11) {
  theme_classic(base_size = base_size) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0, margin = margin(b = 4)),
      plot.subtitle = element_text(color = "grey30", lineheight = 1.05, margin = margin(b = 6)),
      axis.title = element_text(face = "bold"),
      axis.text = element_text(color = "black"),
      strip.background = element_rect(fill = "grey92", color = NA),
      strip.text = element_text(face = "bold"),
      legend.title = element_text(face = "bold"),
      legend.text = element_text(color = "black")
    )
}

theme_minimal_ckm <- function(base_size = 11) {
  theme_minimal(base_size = base_size) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0, margin = margin(b = 4)),
      plot.subtitle = element_text(color = "grey30", lineheight = 1.05, margin = margin(b = 6)),
      axis.title = element_text(face = "bold"),
      axis.text = element_text(color = "black"),
      strip.background = element_rect(fill = "grey92", color = NA),
      strip.text = element_text(face = "bold"),
      legend.title = element_text(face = "bold"),
      panel.grid.minor = element_blank()
    )
}

stage_palette <- c(
  `0` = "#2E8B57",  # s0 green
  `1` = "#2C7FB8",  # s1 blue
  `2` = "#F39C12",  # s2 orange
  `3` = "#8E63CE",  # s3 purple
  `4` = "#D73027",  # s4 red
  `5` = "#000000"   # death black
)
stage_short_labels <- c(`0` = "s0", `1` = "s1", `2` = "s2", `3` = "s3", `4` = "s4", `5` = "death")
stage_full_labels <- c(`0` = "Stage 0", `1` = "Stage 1", `2` = "Stage 2", `3` = "Stage 3", `4` = "Stage 4", `5` = "Death")
adjust_palette <- c(crude = "#F26C64", demog = "#7CAE00", le8_behavior = "#00BFC4", le8 = "#00BFC4", le8_full = "#C77CFF", full = "#C77CFF")
risk_horizon_palette <- c(`5-year risk` = "#F26C64", `10-year risk` = "#00BFC4")
model_palette <- c(
  M0_base = "#7F8C8D",
  M1_categorical_stage = "#2C7FB8",
  M2_risk_age_stage = "#F39C12",
  M2b_explainable_Yin_score = "#4D4D4D",
  M3_risk_age_stage_YY_score = "#D73027",
  M3a_persistent_YY = "#1B9E77",
  M3b_stage_context_YY = "#7570B3",
  M3c_yang_supported_YY = "#E6AB02",
  M4_Yin_glmnet_benchmark = "#8E63CE"
)
biom_refine_palette <- c(
  Dnull_always_stage2 = "#BDBDBD",
  D0_demographic = "#7F8C8D",
  D1_explainable_biomarker_stage = "#D73027"
)

ckm_stage6_factor <- function(stage, death = NULL, short = TRUE) {
  s <- safe_num(stage)
  if (!is.null(death)) s[death %in% TRUE] <- 5
  lab <- if (short) stage_short_labels else stage_full_labels
  factor(s, levels = 0:5, labels = unname(lab[as.character(0:5)]))
}

scale_stage_color <- function(...) scale_color_manual(values = stage_palette, ...)
scale_stage_fill <- function(...) scale_fill_manual(values = stage_palette, ...)


candidate_date_cols <- function(dat, y, source = "all") {
  if (identical(y, "death")) return(intersect(c("date_death", "fod_icd10_death"), names(dat)))
  # Strict ICD-10 mode means the phenotype generated by phe.R/get_icd_fast():
  # fod_icd10_[trait]. It does not silently fall back to aggregate or other-source dates.
  icd10 <- intersect(paste0("fod_icd10_", y), names(dat))
  all <- intersect(c(
    # Prefer an aggregate first-occurrence date created by the phenotype/t2e
    # pipeline when it exists, then retain each source-specific date for audit.
    paste0("fod_", y), paste0("fod_all_", y), paste0("fod_combined_", y),
    paste0("fod_icd10_", y), paste0("fod_icd9_", y), paste0("fod_srd_", y),
    paste0("fod_gp_", y), paste0("fod_ref_", y), paste0("date_", y), paste0(y, "_date")
  ), names(dat))
  if (y %in% c("cvd_cad", "ihd", "mi")) all <- unique(c(all, intersect("date_mi", names(dat))))
  if (y %in% c("cvd_stroke", "stroke")) all <- unique(c(all, intersect("date_stroke", names(dat))))
  if (source == "icd10") icd10 else all
}

trait_prevalent_simple <- function(dat, y, baseline = "date_attend") {
  cols <- candidate_date_cols(dat, y, "all")
  if (length(cols)) {
    d <- pmin_date2(dat[, cols, drop = FALSE])
    flag <- !is.na(d) & !is.na(dat[[baseline]]) & d <= as_date2(dat[[baseline]])
  } else if (y %in% names(dat)) {
    # Use an undated phenotype flag only when no dated source exists.  Otherwise
    # an "ever diagnosed" flag could incorrectly turn post-baseline disease into
    # a baseline CKM criterion.
    flag <- safe_num(dat[[y]]) == 1
  } else {
    flag <- rep(FALSE, nrow(dat))
  }
  flag[is.na(flag)] <- FALSE
  flag
}

mode_value <- function(x) {
  z <- x[!is.na(x)]
  if (!length(z)) return(NA)
  names(sort(table(z), decreasing = TRUE))[1]
}
