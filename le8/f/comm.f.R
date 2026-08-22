# Shared functions for the LE8-supervised 5C omics pipeline.
# The functions in this file deliberately keep protein and metabolite branches separate.

.required_pkgs <- c("data.table", "dplyr", "tidyr", "purrr", "stringr", "tibble",
                    "ggplot2", "ggrepel", "patchwork", "survival", "scales", "openxlsx", "forcats")
if (!requireNamespace("pacman", quietly = TRUE)) install.packages("pacman", repos = "https://cloud.r-project.org")
suppressPackageStartupMessages(
  pacman::p_load(char = .required_pkgs)
)

`%||%` <- function(x, y) if (is.null(x) || length(x) == 0) y else x
bt <- function(x) paste0("`", x, "`")

# Reproducible outcome-stratified folds shared by C5 and S2. Prediction-model
# functions themselves live in 0f/pred.f.R and C5-only extensions in c5_pred.R.
make_folds <- function(dat, event_var, k = 5, seed = 2026) {
  set.seed(seed)
  fold <- integer(nrow(dat)); event <- dat[[event_var]]
  for (value in unique(event[!is.na(event)])) {
    id <- which(event == value)
    fold[id] <- sample(rep(seq_len(k), length.out = length(id)))
  }
  fold
}
cap <- function(x, limit) pmax(pmin(x, limit), -limit)
safe_log <- function(x) ifelse(is.finite(x) & x > 0, log(x), NA_real_)
truthy <- function(x) toupper(as.character(x)) %in% c("TRUE", "T", "1", "YES", "Y")
std_num <- function(x) {
  x <- suppressWarnings(as.numeric(x)); s <- stats::sd(x, na.rm = TRUE)
  if (!is.finite(s) || s <= 0) return(rep(NA_real_, length(x)))
  as.numeric(scale(x))
}
inormal2 <- function(x) {
  x <- suppressWarnings(as.numeric(x)); n <- sum(is.finite(x))
  if (n < 5 || !is.finite(sd(x, na.rm = TRUE)) || sd(x, na.rm = TRUE) == 0) return(rep(NA_real_, length(x)))
  qnorm((rank(x, na.last = "keep", ties.method = "average") - 0.5) / n)
}

# -----------------------------------------------------------------------------
# Paths and global settings
# -----------------------------------------------------------------------------
if (!exists("dir0")) dir0 <- ifelse(Sys.info()[["sysname"]] == "Windows", "D:", "/mnt/d")
find_first_existing <- function(paths, default = paths[[1]]) {
  paths <- unique(paths[!is.na(paths) & nzchar(paths)])
  hit <- paths[file.exists(paths) | dir.exists(paths)]
  if (length(hit)) normalizePath(hit[[1]], winslash = "/", mustWork = FALSE) else default
}
Y <- if (exists("Y")) Y else Sys.getenv("Y", unset = "cvd_cad")
BIOM <- if (exists("BIOM")) BIOM else Sys.getenv("BIOM", unset = "prot,met")
LE8_JOB <- if (exists("LE8_JOB")) LE8_JOB else Sys.getenv("LE8_JOB", unset = "")
prot_DO <- truthy(Sys.getenv("PROT_DO", unset = "TRUE"))
met_DO <- truthy(Sys.getenv("MET_DO", unset = "TRUE"))
LE8_FORCE <- truthy(Sys.getenv("LE8_FORCE", unset = "FALSE"))
N_CORES <- max(1L, suppressWarnings(as.integer(Sys.getenv("N_CORES", unset = "4"))))
if (!is.finite(N_CORES)) N_CORES <- 1L
SEED <- as.integer(Sys.getenv("SEED", unset = "2026"))
set.seed(SEED)

indir <- find_first_existing(c(Sys.getenv("UKB_PHE", unset = ""), file.path(dir0, "data/ukb/phe"),
                               "/mnt/ddata/ukb/phe", "D:/data/ukb/phe"), file.path(dir0, "data/ukb/phe"))
common_dir <- find_first_existing(c(Sys.getenv("UKB_COMMON", unset = ""), file.path(indir, "common"),
                                    "/mnt/ddata/ukb/phe/common", "D:/data/ukb/phe/common"), file.path(indir, "common"))
analysis_root <- Sys.getenv("LE8_ANALYSIS_ROOT", unset = file.path(dir0, "analysis/le8"))
out.base <- file.path(analysis_root, Y)
out.prot <- file.path(out.base, "prot")
out.met <- file.path(out.base, "met")
invisible(lapply(c(out.base, out.prot, out.met), dir.create, recursive = TRUE, showWarnings = FALSE))

dir.Y <- Sys.getenv("LE8_GWAS_DIR", unset = file.path(dir0, "data.BIG/gwas/main/clean"))
dir.X <- Sys.getenv("LE8_PQTL_IV_DIR", unset = file.path(dir0, "data.BIG/gwas/ppp/clean"))
dir.met.gwas <- Sys.getenv("LE8_MQTL_IV_DIR", unset = file.path(dir0, "data.BIG/gwas/met/clean"))
ppp_bed_file <- find_first_existing(c(Sys.getenv("LE8_PPP_BED", unset = ""),
                                      file.path(dir0, "data.BIG/gwas/ppp/ppp_3k.b38.bed"),
                                      file.path(dir0, "data.BIG/gwas/ppp/ppp_3k_b38.bed"),
                                      "D:/data.BIG/gwas/ppp/ppp_3k.b38.bed"),
                                    file.path(dir0, "data.BIG/gwas/ppp/ppp_3k.b38.bed"))
met_list_file <- find_first_existing(c(Sys.getenv("LE8_MET_LIST", unset = ""), file.path(common_dir, "met.lst"),
                                       "D:/data/ukb/phe/common/met.lst"), file.path(common_dir, "met.lst"))
ukb_bgen_dir <- Sys.getenv("UKB_BGEN_DIR", unset = "/mnt/d/data/ukb/gen/imp")

# Source existing project helpers when present. The new code does not require assoc.f.R/mr.f.R,
# but uses the user's existing t2e() definition from 0phe.f.R when available.
helper_names <- c("0phe.f.R", "assoc.f.R", "mr.f.R", "plot.f.R", "pred.f.R")
helper_dirs <- file.path(dir0, "scripts/0f")
for (f0 in unique(unlist(lapply(helper_dirs, function(d0) file.path(d0, helper_names))))) {
  if (file.exists(f0)) try(source(f0), silent = TRUE)
}
if (!exists("date_follow_end")) date_follow_end <- as.Date(Sys.getenv("DATE_FOLLOW_END", unset = "2023-04-01"))
if (!exists("vars.basic")) stop("vars.basic was not loaded from scripts/0f/0phe.f.R.", call. = FALSE)
if (!exists("vars.le8")) stop("vars.le8 was not loaded from scripts/0f/0phe.f.R.", call. = FALSE)
if (!exists("names.le8")) stop("names.le8 was not loaded from scripts/0f/0phe.f.R.", call. = FALSE)
if (!exists("vars.adj2")) vars.adj2 <- unique(c(vars.basic, vars.le8))
# Analysis-wide covariate choice.  Keep vars.basic as the default; replace the
# next line with `covs_use <- vars.adj2` when downstream modules should consume
# the comprehensively adjusted C1 scan.
if (!exists("covs_use")) covs_use <- vars.adj2
covs_use_name <- if (identical(unique(covs_use), unique(vars.adj2))) "adj2" else "basic"
LE8_LABS <- c(diet = "Diet", pa = "Physical activity", smoke = "Smoking", bmi = "BMI",
              nonhdl = "Non-HDL-C", hba1c = "HbA1c", bp = "Blood pressure", sleep = "Sleep")
cols_le8 <- c(diet = "#E68613", pa = "#5B9BD5", smoke = "#CC79A7", bmi = "#E76F51",
              nonhdl = "#00A6B2", hba1c = "#009E73", bp = "#B79F00", sleep = "#4063D8")
cols_evidence <- c(Observational = "#4C78A8", MR = "#F58518", Colocalization = "#54A24B",
                   Connection = "#B279A2", Prediction = "#E45756")

# -----------------------------------------------------------------------------
# IO and output management
# -----------------------------------------------------------------------------
required_file <- function(path, label = "file") {
  if (is.na(path) || !nzchar(path) || !file.exists(path) || isTRUE(file.size(path) == 0))
    stop("Missing required ", label, ": ", path, call. = FALSE)
  normalizePath(path, winslash = "/", mustWork = FALSE)
}
read_all <- function(select_vars = NULL) {
  x <- readRDS(required_file(file.path(indir, "Rdata/all.rds"), "UKB phenotype all.rds"))
  if (!is.null(select_vars)) x <- x[, intersect(unique(select_vars), names(x)), drop = FALSE]
  x
}
read_prot <- function(required = TRUE) {
  f <- file.path(indir, "Rdata/prot.rds")
  if (!file.exists(f) || file.size(f) == 0) {
    if (required) stop("Missing proteomics RDS: ", f, call. = FALSE)
    return(NULL)
  }
  readRDS(f)
}
read_met <- function(required = TRUE) {
  f <- file.path(indir, "Rdata/met.rds")
  if (!file.exists(f) || file.size(f) == 0) {
    if (required) stop("Missing metabolomics RDS: ", f, call. = FALSE)
    return(NULL)
  }
  readRDS(f)
}
setwd2 <- function(path) { dir.create(path, recursive = TRUE, showWarnings = FALSE); setwd(path); invisible(path) }
le8_job_dir <- function(outdir = getwd(), job = LE8_JOB) {
  if (is.null(job) || !nzchar(job)) job <- "raw"
  bn <- basename(normalizePath(outdir, winslash = "/", mustWork = FALSE))
  if (bn %in% c("prot", "met")) file.path(outdir, job) else outdir
}
safe_sheet <- function(x) {
  x <- gsub("[\\[\\]:*?/\\\\]", "_", x)
  x <- substr(x, 1, 31); make.unique(x, sep = "_")
}
write_raw_csv <- function(x, file, rawdir = le8_job_dir()) {
  dir.create(rawdir, recursive = TRUE, showWarnings = FALSE)
  data.table::fwrite(as.data.frame(x), file.path(rawdir, file), sep = ",", na = "NA")
  invisible(file.path(rawdir, file))
}
write_raw_tsv <- function(x, file, rawdir = le8_job_dir()) {
  dir.create(rawdir, recursive = TRUE, showWarnings = FALSE)
  data.table::fwrite(as.data.frame(x), file.path(rawdir, file), sep = "\t", na = "")
  invisible(file.path(rawdir, file))
}
write_xlsx2 <- function(x, file) {
  if (is.data.frame(x)) x <- list(data = x)
  x <- x[!vapply(x, is.null, logical(1))]
  if (!length(x)) x <- list(note = data.frame(note = "No table was produced."))
  names(x) <- safe_sheet(names(x))
  x <- lapply(x, function(z) {
    if (is.matrix(z)) z <- as.data.frame(z)
    if (!is.data.frame(z)) z <- data.frame(value = I(list(z)))
    # openxlsx::write.xlsx(asTable = TRUE) cannot serialize a data frame with
    # zero columns and fails internally with `df[[i]]: subscript out of bounds`.
    # Keep the worksheet, but make the absence of rows explicit.
    if (ncol(z) == 0L) z <- data.frame(note = "No table was produced.")
    z
  })
  openxlsx::write.xlsx(x, file, overwrite = TRUE, asTable = TRUE, freezePane = TRUE, autoFilter = TRUE)
  invisible(file)
}
module_cache <- function(outdir, module = LE8_JOB) file.path(le8_job_dir(outdir, module), paste0(module, ".res.rds"))
cache_valid <- function(path) !LE8_FORCE && file.exists(path) && file.size(path) > 0
cache_message <- function(label, path) {
  message(label, " already exists, skip ", label, " step: ", path,
          ". Delete this existing file (or its step folder) to re-run.")
}
parallel_map <- function(x, fun) {
  if (.Platform$OS.type != "windows" && N_CORES > 1L && length(x) > 1L) {
    message("Parallel scan: ", min(N_CORES, length(x)), " workers for ", length(x), " tasks")
    parallel::mclapply(x, fun, mc.cores = min(N_CORES, length(x)), mc.preschedule = TRUE)
  } else lapply(x, fun)
}
module_meta <- function(layer, module = LE8_JOB, extra = list()) {
  c(list(module = module, layer = layer, trait = Y, generated = format(Sys.time(), "%F %T %z"),
         seed = SEED, R = R.version.string, biom = BIOM), extra)
}
finalize_outputs <- function(module, outdir = getwd()) {
  rawdir <- le8_job_dir(outdir, module); dir.create(rawdir, recursive = TRUE, showWarnings = FALSE)
  png <- list.files(outdir, pattern = "\\.png$", full.names = FALSE)
  xlsx <- list.files(outdir, pattern = "\\.xlsx$", full.names = FALSE)
  index <- data.frame(module = module, trait = Y, directory = normalizePath(outdir, winslash = "/", mustWork = FALSE),
                      png = paste(png, collapse = "; "), workbook = paste(xlsx, collapse = "; "),
                      generated = format(Sys.time(), "%F %T %z"))
  write_raw_csv(index, paste0(module, ".output_index.csv"), rawdir)
  message("Completed ", module, ": ", outdir)
  invisible(index)
}

# -----------------------------------------------------------------------------
# Plotting
# -----------------------------------------------------------------------------
theme_5c <- function(base_size = 12) {
  theme_classic(base_size = base_size) +
    theme(plot.title = element_text(face = "bold", hjust = 0),
          plot.subtitle = element_text(color = "grey35"),
          axis.text = element_text(color = "black"),
          legend.title = element_text(face = "bold"),
          panel.grid.major.y = element_line(color = "grey91", linewidth = 0.25),
          strip.background = element_blank(), strip.text = element_text(face = "bold"),
          plot.margin = margin(8, 12, 8, 12))
}
save_plot <- function(p, file, w = 8, h = 6, dpi = 320, outdir = getwd()) {
  ggplot2::ggsave(file.path(outdir, file), p, width = w, height = h, dpi = dpi, bg = "white", limitsize = FALSE)
  invisible(file.path(outdir, file))
}
blank_plot <- function(title, subtitle = "Required data were not available") {
  ggplot() + theme_void(base_size = 13) +
    annotate("text", x = 0, y = .08, label = title, fontface = "bold", size = 5) +
    annotate("text", x = 0, y = -.06, label = subtitle, color = "grey35", size = 3.7)
}
forest_theme <- function(base_size = 10) theme_5c(base_size) + theme(panel.grid.major.y = element_blank())

# -----------------------------------------------------------------------------
# Phenotype and association models
# -----------------------------------------------------------------------------
make_outcome <- function(dat, outcome = Y) {
  tvar <- paste0(outcome, ".t2e"); evar <- paste0(outcome, ".Yt2e")
  if (all(c(tvar, evar) %in% names(dat))) return(dat)
  ydate <- paste0("fod_icd10_", outcome)
  need <- c("birth_date", "date_attend", "date_lost", "date_death", ydate)
  miss <- setdiff(need, names(dat))
  if (length(miss)) stop("Cannot construct outcome ", outcome, "; missing: ", paste(miss, collapse = ", "), call. = FALSE)
  if (!exists("t2e", mode = "function")) stop("t2e() was not found; keep /mnt/d/scripts/0f/0phe.f.R available.", call. = FALSE)
  old <- grep(paste0("^", outcome, "\\.(Y)?(t2e|r2e|b2e)$"), names(dat), value = TRUE)
  if (length(old)) dat[old] <- NULL
  # The user's current t2e() calls use NA as the second argument. A fallback is retained for older versions.
  first_error <- NULL
  ans <- tryCatch(t2e(dat, NA, ydate, "birth_date", "date_attend", "date_lost", "date_death",
                      date_follow_end, outcome, "year"), error = function(e) { first_error <<- conditionMessage(e); NULL })
  if (is.null(ans)) ans <- tryCatch(t2e(dat, "cvd", ydate, "birth_date", "date_attend", "date_lost", "date_death",
                                       date_follow_end, outcome, "year"), error = function(e) {
                                         stop("t2e() failed: ", conditionMessage(e),
                                              if (!is.null(first_error)) paste0("; first attempt: ", first_error) else "",
                                              call. = FALSE)
                                       })
  if (!all(c(tvar, evar) %in% names(ans))) stop("t2e() did not create ", tvar, " and ", evar,
                                                "; returned columns include: ", paste(grep("t2e|r2e|b2e", names(ans), value=TRUE), collapse=", "),
                                                call. = FALSE)
  ans
}
filter_analysis_cohort <- function(dat) {
  if ("ethnic.c" %in% names(dat) && truthy(Sys.getenv("LE8_WHITE_ONLY", unset = "TRUE"))) {
    keep <- as.character(dat$ethnic.c) %in% c("White", "1")
    dat <- dat[keep, , drop = FALSE]
  }
  dat
}
cox_scan <- function(dat, xs, covars, outcome = Y, scale_x = TRUE, min_n = 500, min_event = 20,
                     time_var = paste0(outcome, ".t2e"), event_var = paste0(outcome, ".Yt2e")) {
  tvar <- time_var; evar <- event_var
  xs <- intersect(xs, names(dat)); covars <- intersect(covars, names(dat))
  bind_rows(parallel_map(xs, function(x) {
    need <- unique(c(tvar, evar, x, covars)); d <- dat[, need, drop = FALSE]
    d <- d[stats::complete.cases(d), , drop = FALSE]
    ne <- sum(d[[evar]] == 1, na.rm = TRUE)
    empty <- tibble(term = x, estimate = NA_real_, beta = NA_real_, std.error = NA_real_, conf.low = NA_real_,
                    conf.high = NA_real_, statistic = NA_real_, p.value = NA_real_, N_total = nrow(d), N_event = ne)
    if (nrow(d) < min_n || ne < min_event || length(unique(d[[evar]])) < 2) return(empty)
    d[[x]] <- suppressWarnings(as.numeric(d[[x]])); sx <- sd(d[[x]], na.rm = TRUE)
    if (!is.finite(sx) || sx == 0) return(empty)
    if (scale_x) d[[x]] <- as.numeric(scale(d[[x]]))
    fit <- tryCatch(coxph(as.formula(paste0("Surv(", bt(tvar), ",", bt(evar), ") ~ ",
                                             paste(bt(c(x, covars)), collapse = " + "))), d, ties = "efron"), error = function(e) NULL)
    if (is.null(fit)) return(empty)
    sm <- coef(summary(fit)); if (!x %in% rownames(sm)) return(empty)
    b <- sm[x, "coef"]; se <- sm[x, "se(coef)"]
    tibble(term = x, estimate = exp(b), beta = b, std.error = se, conf.low = exp(b - 1.96 * se),
           conf.high = exp(b + 1.96 * se), statistic = b / se,
           p.value = 2 * pnorm(abs(b / se), lower.tail = FALSE), N_total = nrow(d), N_event = ne)
  })) |> mutate(FDR = p.adjust(p.value, "BH")) |> arrange(p.value)
}
lm_scan <- function(dat, xs, y, covars, scale_x = TRUE, min_n = 500) {
  xs <- intersect(xs, names(dat)); covars <- intersect(covars, names(dat))
  bind_rows(parallel_map(xs, function(x) {
    d <- dat[, unique(c(y, x, covars)), drop = FALSE]; d <- d[complete.cases(d), , drop = FALSE]
    empty <- tibble(term = x, beta = NA_real_, std.error = NA_real_, statistic = NA_real_, p.value = NA_real_, N_total = nrow(d))
    if (nrow(d) < min_n) return(empty)
    d[[x]] <- suppressWarnings(as.numeric(d[[x]])); d[[y]] <- suppressWarnings(as.numeric(d[[y]]))
    if (!is.finite(sd(d[[x]], na.rm = TRUE)) || sd(d[[x]], na.rm = TRUE) == 0) return(empty)
    if (scale_x) d[[x]] <- std_num(d[[x]])
    fit <- tryCatch(lm(as.formula(paste0(bt(y), " ~ ", paste(bt(c(x, covars)), collapse = " + "))), d), error = function(e) NULL)
    if (is.null(fit)) return(empty)
    sm <- coef(summary(fit)); if (!x %in% rownames(sm)) return(empty)
    tibble(term = x, beta = sm[x, "Estimate"], std.error = sm[x, "Std. Error"], statistic = sm[x, "t value"],
           p.value = sm[x, "Pr(>|t|)"], N_total = nrow(d))
  })) |> mutate(FDR = p.adjust(p.value, "BH")) |> arrange(p.value)
}
partial_r2 <- function(full, reduced) {
  r2f <- summary(full)$r.squared; r2r <- summary(reduced)$r.squared
  pmax(0, pmin(1, r2f - r2r))
}

# -----------------------------------------------------------------------------
# Protein and metabolite annotations
# -----------------------------------------------------------------------------
read_ppp_bed <- function(proteins = NULL) {
  f <- required_file(ppp_bed_file, "PPP bed annotation")
  x <- data.table::fread(f, header = FALSE, fill = TRUE, showProgress = FALSE)
  if (ncol(x) < 5) stop("PPP bed must have at least five columns (chr,start,end,protein,group): ", f, call. = FALSE)
  ans <- tibble(chr = as.character(x[[1]]), start = as.numeric(x[[2]]), end = as.numeric(x[[3]]),
                protein = as.character(x[[4]]), group = as.character(x[[5]])) |>
    mutate(chr = str_remove(chr, "^chr"), pos = floor((start + end) / 2),
           group = coalesce(na_if(group, ""), "Other")) |> distinct(protein, .keep_all = TRUE)
  if (!is.null(proteins)) tibble(protein = proteins) |> left_join(ans, by = "protein") else ans
}
met_super_group <- function(group, subgroup = NA_character_) {
  txt <- str_to_lower(paste(coalesce(group, ""), coalesce(subgroup, "")))
  case_when(str_detect(txt, "amino") ~ "Amino acids",
            str_detect(txt, "apolipoprotein|cholesterol|cholesteryl") ~ "Cholesterol and apolipoproteins",
            str_detect(txt, "fatty acid|omega|pufa|mufa|sfa") ~ "Fatty acids",
            str_detect(txt, "triglyceride|lipoprotein|hdl|ldl|vldl|lipid|phospholipid|sphingomyelin|ceramide") ~ "Lipoprotein lipids",
            str_detect(txt, "glycolysis|glucose|ketone|citrate|lactate") ~ "Energy metabolism",
            TRUE ~ "Other metabolites")
}
read_met_annotation <- function(vars.met = NULL) {
  f <- required_file(met_list_file, "met.lst")
  # Let fread detect whether the file has a header; several historical met.lst versions are headerless.
  x <- data.table::fread(f, fill = TRUE, check.names = FALSE, showProgress = FALSE)
  nm0 <- names(x); nml <- tolower(nm0)
  pick <- function(pattern, fallback) { i <- grep(pattern, nml)[1]; if (is.na(i)) i <- fallback; nm0[[i]] }
  var_col <- pick("^var_name$|^variable$|^trait$|^field$|^name$", 1)
  if (!is.null(vars.met) && length(vars.met)) {
    vm <- unique(c(vars.met, str_remove(vars.met, "^met_"), paste0("met_", str_remove(vars.met, "^met_"))))
    overlap <- vapply(x, function(z) sum(unique(as.character(z)) %in% vm, na.rm = TRUE), numeric(1))
    if (length(overlap) && max(overlap, na.rm = TRUE) > 0) var_col <- nm0[[which.max(overlap)]]
  }
  label_col <- pick("^full_name$|biomarker|description|label|long_name", min(2, ncol(x)))
  # Prefer the explicit header in current met.lst files.  Only use the
  # historical penultimate-column convention when no group header exists.
  group_i <- grep("^group$|^class$|^category$", nml)[1]
  group_col <- if (!is.na(group_i)) nm0[[group_i]] else nm0[[max(1, ncol(x) - 1)]]
  subgroup_col <- if (any(grepl("subgroup|sub_group|class", nml))) nm0[[grep("subgroup|sub_group|class", nml)[1]]] else nm0[[ncol(x)]]
  ans <- x |> transmute(trait0 = as.character(.data[[var_col]]), label = as.character(.data[[label_col]]),
                       group = as.character(.data[[group_col]]), subgroup = as.character(.data[[subgroup_col]])) |>
    mutate(trait_prefixed = ifelse(str_detect(trait0, "^met_"), trait0, paste0("met_", trait0)),
           trait_stripped = str_remove(trait0, "^met_"),
           trait = case_when(!is.null(vars.met) & trait0 %in% vars.met ~ trait0,
                             !is.null(vars.met) & trait_prefixed %in% vars.met ~ trait_prefixed,
                             !is.null(vars.met) & trait_stripped %in% vars.met ~ trait_stripped,
                             TRUE ~ trait_prefixed),
           label = coalesce(na_if(label, ""), trait), group = coalesce(na_if(group, ""), "Other"),
           subgroup = coalesce(na_if(subgroup, ""), group), super_group = met_super_group(group, subgroup)) |>
    select(-trait_prefixed, -trait_stripped) |> distinct(trait, .keep_all = TRUE)
  if (!is.null(vars.met)) tibble(trait = vars.met) |> left_join(ans, by = "trait") |>
      mutate(label = coalesce(label, trait), group = coalesce(group, "Other"), subgroup = coalesce(subgroup, group),
             super_group = coalesce(super_group, "Other metabolites")) else ans
}
layer_annotation <- function(layer, features) {
  if (layer == "protein") read_ppp_bed(features) |> transmute(feature = protein, label = protein, group, chr, start, end, pos)
  else read_met_annotation(features) |> transmute(feature = trait, label, group, subgroup, super_group)
}
qtl_name_candidates <- function(feature, layer) {
  z <- unique(c(feature, str_remove(feature, "^met_"), str_replace_all(str_remove(feature, "^met_"), "[^A-Za-z0-9_.-]", "_")))
  z[nzchar(z)]
}
find_qtl_directory <- function(feature, base_dir, layer) {
  cand <- qtl_name_candidates(feature, layer)
  hit <- cand[dir.exists(file.path(base_dir, cand))][1]
  if (length(hit) && !is.na(hit)) list(name = hit, dir = file.path(base_dir, hit)) else list(name = cand[[1]], dir = file.path(base_dir, cand[[1]]))
}
find_qtl_files <- function(feature, base_dir, layer = c("protein", "metabolite")) {
  layer <- match.arg(layer); loc <- find_qtl_directory(feature, base_dir, layer); nm <- loc$name; d <- loc$dir
  first <- function(x) { y <- x[file.exists(x) & file.size(x) > 0]; if (length(y)) normalizePath(y[[1]], winslash = "/", mustWork = FALSE) else NA_character_ }
  list(feature = feature, qtl_name = nm, dir = d,
       joint = first(c(file.path(d, paste0(nm, ".jma.cojo")), file.path(d, paste0(feature, ".jma.cojo")))),
       full = first(c(file.path(d, paste0(nm, ".gz")), file.path(d, paste0(feature, ".gz")), file.path(d, paste0(nm, ".4gcta")))),
       cis = first(c(file.path(d, paste0(nm, ".cis.gz")), file.path(d, paste0(feature, ".cis.gz")))),
       trans = first(c(file.path(d, paste0(nm, ".trans.gz")), file.path(d, paste0(nm, ".clump.assoc")), file.path(d, paste0(feature, ".clump.assoc")))) )
}

# -----------------------------------------------------------------------------
# Standardised summary statistics, harmonisation, MR and R2
# -----------------------------------------------------------------------------
.norm_name <- function(x) toupper(gsub("[^A-Za-z0-9]", "", x))
.pick_col <- function(nms, candidates, prefer = NULL) {
  nn <- .norm_name(nms); cc <- .norm_name(candidates)
  if (!is.null(prefer)) {
    pp <- .norm_name(prefer); i <- match(pp, nn); i <- i[!is.na(i)]; if (length(i)) return(nms[[i[[1]]]])
  }
  i <- match(cc, nn); i <- i[!is.na(i)]; if (length(i)) nms[[i[[1]]]] else NA_character_
}
read_table_auto <- function(file, nrows = -1) {
  if (is.na(file) || !file.exists(file) || file.size(file) == 0) return(data.table())
  data.table::fread(file, nrows = nrows, fill = TRUE, showProgress = FALSE, check.names = FALSE)
}
standardize_sumstat <- function(d, N_default = as.numeric(Sys.getenv("C2_SUMSTAT_N", unset = "100000")),
                                joint = FALSE, source_file = NA_character_) {
  d <- as.data.frame(d, check.names = FALSE); nms <- names(d)
  col <- list(
    SNP = .pick_col(nms, c("SNP","RSID","RS_NUMBER","VARIANT_ID","ID","MARKERNAME","POS_NAME")),
    CHR = .pick_col(nms, c("CHR","CHROM","CHROMOSOME")),
    POS = .pick_col(nms, c("POS","BP","POSITION","BASE_PAIR_LOCATION")),
    EA = .pick_col(nms, c("EA","A1","ALT","ALLELE1","EFFECT_ALLELE","EFFECTALLELE")),
    NEA = .pick_col(nms, c("NEA","A2","REF","ALLELE0","OTHER_ALLELE","REFERENCE_ALLELE","NONEFFECTALLELE")),
    REF_A = .pick_col(nms, c("REFA")),
    EAF = .pick_col(nms, c("EAF","FREQ","FREQ_GENO","A1FREQ","EFFECT_ALLELE_FREQUENCY")),
    BETA = if (joint) .pick_col(nms, c("BJ","BETAJ","BETA","B","EFFECT"), prefer = c("BJ","BETAJ")) else .pick_col(nms, c("BETA","B","EFFECT","LOG_ODDS","LOGOR","BJ")),
    SE = if (joint) .pick_col(nms, c("BJ_SE","SEJ","SE","STDERR","SEBETA"), prefer = c("BJ_SE","SEJ")) else .pick_col(nms, c("SE","STDERR","SEBETA","BJ_SE")),
    P = if (joint) .pick_col(nms, c("PJ","P_J","P","PVAL","PVALUE"), prefer = c("PJ","P_J")) else .pick_col(nms, c("P","PVAL","P_VALUE","PVALUE","PJ")),
    N = .pick_col(nms, c("N","SAMPLESIZE","N_TOTAL","NEFF","N_IIDS"))
  )
  getv <- function(nm, default = NA) if (!is.na(col[[nm]]) && col[[nm]] %in% names(d)) d[[col[[nm]]]] else rep(default, nrow(d))
  ea <- toupper(as.character(getv("EA"))); nea <- toupper(as.character(getv("NEA"))); refa <- toupper(as.character(getv("REF_A")))
  if (joint) {
    # GCTA .jma.cojo supplies refA but not both alleles. Keep refA as a provisional effect allele;
    # read_qtl_instruments() replaces it using the corresponding full QTL file.
    ea[is.na(ea) | ea == ""] <- refa[is.na(ea) | ea == ""]
  }
  ans <- tibble(SNP = as.character(getv("SNP")), CHR = str_remove(as.character(getv("CHR")), regex("^chr", ignore_case = TRUE)),
                POS = suppressWarnings(as.numeric(getv("POS"))), EA = ea, NEA = nea,
                EAF = suppressWarnings(as.numeric(getv("EAF"))), BETA = suppressWarnings(as.numeric(getv("BETA"))),
                SE = suppressWarnings(as.numeric(getv("SE"))), P = suppressWarnings(as.numeric(getv("P"))),
                N = suppressWarnings(as.numeric(getv("N"))), source_file = source_file, joint = joint)
  ans$N[!is.finite(ans$N)] <- N_default
  ans$P[!is.finite(ans$P) | ans$P <= 0 | ans$P > 1] <- 2 * pnorm(abs(ans$BETA[!is.finite(ans$P) | ans$P <= 0 | ans$P > 1] /
                                                                          ans$SE[!is.finite(ans$P) | ans$P <= 0 | ans$P > 1]), lower.tail = FALSE)
  ans |> filter(!is.na(SNP), SNP != "", is.finite(BETA), is.finite(SE), SE > 0) |> distinct(SNP, .keep_all = TRUE)
}
read_sumstat <- function(file, N_default = as.numeric(Sys.getenv("C2_SUMSTAT_N", unset = "100000")), joint = grepl("jma\\.cojo$", file %||% "")) {
  if (is.na(file) || !file.exists(file) || file.size(file) == 0) return(tibble())
  standardize_sumstat(read_table_auto(file), N_default, joint = joint, source_file = file)
}
read_sumstat_region <- function(file, chr, start, end, N_default = as.numeric(Sys.getenv("C2_SUMSTAT_N", unset = "100000"))) {
  if (is.na(file) || !file.exists(file) || file.size(file) == 0) return(tibble())
  # Try a streaming awk extraction for large gzipped files. Fall back to fread on the full file.
  hdr <- tryCatch(readLines(if (grepl("\\.gz$", file)) gzfile(file, "rt") else base::file(file, "rt"), n = 1), error = function(e) "")
  d <- NULL
  if (nzchar(hdr) && Sys.info()[["sysname"]] != "Windows" && nzchar(Sys.which("awk"))) {
    nms <- strsplit(hdr, "\t", fixed = TRUE)[[1]]
    cchr <- .pick_col(nms, c("CHR","CHROM","CHROMOSOME")); cpos <- .pick_col(nms, c("POS","BP","POSITION","BASE_PAIR_LOCATION"))
    ichr <- match(cchr, nms); ipos <- match(cpos, nms)
    if (is.finite(ichr) && is.finite(ipos)) {
      dec <- if (grepl("\\.gz$", file)) paste("gzip -cd", shQuote(file)) else paste("cat", shQuote(file))
      awk <- sprintf("awk -F '\\t' 'NR==1 || (($%d==\"%s\" || $%d==\"chr%s\") && $%d>=%d && $%d<=%d)'",
                     ichr, chr, ichr, chr, ipos, floor(start), ipos, ceiling(end))
      d <- tryCatch(data.table::fread(cmd = paste(dec, "|", awk), fill = TRUE, showProgress = FALSE, check.names = FALSE), error = function(e) NULL)
    }
  }
  if (is.null(d)) d <- read_table_auto(file)
  z <- standardize_sumstat(d, N_default, joint = grepl("jma\\.cojo$", file), source_file = file)
  z |> filter(CHR == as.character(chr), POS >= start, POS <= end)
}

read_sumstat_snps <- function(file, snps, N_default = as.numeric(Sys.getenv("C2_SUMSTAT_N", unset = "100000"))) {
  snps <- unique(as.character(snps)); snps <- snps[!is.na(snps) & nzchar(snps)]
  if (!length(snps) || is.na(file) || !file.exists(file)) return(tibble())
  hdr <- tryCatch(readLines(if (grepl("\\.gz$", file)) gzfile(file, "rt") else base::file(file, "rt"), n = 1), error = function(e) "")
  d <- NULL
  if (nzchar(hdr) && Sys.info()[["sysname"]] != "Windows" && nzchar(Sys.which("awk"))) {
    nms <- strsplit(hdr, "\t", fixed = TRUE)[[1]]
    csnp <- .pick_col(nms, c("SNP","RSID","RS_NUMBER","VARIANT_ID","ID","MARKERNAME","POS_NAME")); isnp <- match(csnp, nms)
    if (is.finite(isnp)) {
      tf <- tempfile("le8_snps_"); writeLines(snps, tf)
      dec <- if (grepl("\\.gz$", file)) paste("gzip -cd", shQuote(file)) else paste("cat", shQuote(file))
      awk <- sprintf("awk -F '\\t' 'NR==FNR {a[$1]=1; next} FNR==1 || ($%d in a)' %s -", isnp, shQuote(tf))
      d <- tryCatch(data.table::fread(cmd = paste(dec, "|", awk), fill = TRUE, showProgress = FALSE, check.names = FALSE), error = function(e) NULL)
      unlink(tf)
    }
  }
  if (is.null(d)) d <- read_table_auto(file)
  standardize_sumstat(d, N_default, joint = grepl("jma\\.cojo$", file), source_file = file) |> filter(SNP %in% snps)
}

get_y_gwas_file <- function(trait = Y, required = TRUE) {
  p <- Sys.getenv("Y_GWAS", unset = Sys.getenv("CAD_GWAS", unset = ""))
  if (!nzchar(p)) p <- file.path(dir.Y, trait, paste0(trait, ".gz"))
  if (file.exists(p) && file.size(p) > 0) return(normalizePath(p, winslash = "/", mustWork = FALSE))
  if (required) stop("Missing outcome GWAS: ", p, call. = FALSE)
  NA_character_
}
get_case_fraction <- function(outcome = Y) {
  s <- suppressWarnings(as.numeric(Sys.getenv("COLOC_CASE_FRAC", unset = "NA")))
  if (is.finite(s) && s > 0 && s < 1) return(s)
  tryCatch({
    d <- read_all(c("eid", "birth_date", "date_attend", "date_lost", "date_death", paste0("fod_icd10_", outcome), paste0(outcome, ".Yt2e"), paste0(outcome, ".t2e")))
    d <- make_outcome(d, outcome); mean(d[[paste0(outcome, ".Yt2e")]] == 1, na.rm = TRUE)
  }, error = function(e) NA_real_)
}
recover_qtl_alleles <- function(iv, full_file) {
  if (!nrow(iv) || is.na(full_file) || !file.exists(full_file)) return(iv)
  full <- read_sumstat_snps(full_file, iv$SNP)
  if (!nrow(full)) return(iv)
  a <- full |> transmute(SNP, EA_full = toupper(EA), NEA_full = toupper(NEA),
                         EAF_full = EAF, CHR_full = CHR, POS_full = POS)
  iv |> mutate(EA_joint = toupper(EA), NEA_joint = toupper(NEA)) |> left_join(a, by = "SNP") |>
    mutate(same_to_full = EA_joint == EA_full & (is.na(NEA_joint) | NEA_joint == "" | NEA_joint == NEA_full),
           reverse_to_full = EA_joint == NEA_full & (is.na(NEA_joint) | NEA_joint == "" | NEA_joint == EA_full),
           allele_ok = same_to_full | reverse_to_full) |>
    filter(allele_ok) |>
    mutate(BETA = ifelse(reverse_to_full, -BETA, BETA),
           EAF = ifelse(reverse_to_full & is.finite(EAF), 1 - EAF, EAF),
           EA = EA_full, NEA = NEA_full, EAF = coalesce(EAF, EAF_full),
           CHR = coalesce(na_if(CHR, "NA"), CHR_full), POS = coalesce(POS, POS_full)) |>
    select(-ends_with("_full"), -EA_joint, -NEA_joint, -same_to_full, -reverse_to_full, -allele_ok)
}
read_qtl_instruments <- function(feature, base_dir, layer = c("protein", "metabolite"), annotation = NULL) {
  layer <- match.arg(layer); fs <- find_qtl_files(feature, base_dir, layer)
  # MR instruments must be approximately independent. The clean pipeline's .jma.cojo file is therefore
  # the default and required source. A cis/full regional file contains correlated variants and must not
  # silently become an IV set. The explicit fallback is provided only for deliberate diagnostic runs.
  allow_unclumped <- truthy(Sys.getenv("C2_ALLOW_UNCLUMPED_FALLBACK", unset = "FALSE"))
  iv <- if (!is.na(fs$joint)) {
    read_sumstat(fs$joint, joint = TRUE)
  } else if (allow_unclumped && !is.na(fs$cis)) {
    warning("C2 diagnostic fallback: using an unclumped cis QTL file for ", feature,
            ". Do not use this result as the primary MR estimate.")
    read_sumstat(fs$cis)
  } else if (allow_unclumped && !is.na(fs$full)) {
    warning("C2 diagnostic fallback: using an unclumped full QTL file for ", feature,
            ". Do not use this result as the primary MR estimate.")
    read_sumstat(fs$full)
  } else tibble()
  if (nrow(iv) && !is.na(fs$full)) iv <- recover_qtl_alleles(iv, fs$full)
  if (!nrow(iv)) return(list(instruments = tibble(), files = fs))
  iv <- iv |> filter(!is.na(EA), EA != "", !is.na(NEA), NEA != "")
  if (layer == "protein" && !is.null(annotation) && nrow(annotation)) {
    a <- annotation[annotation$feature == feature, , drop = FALSE] |> slice(1)
    if (nrow(a) && is.finite(a$start) && is.finite(a$end) && !is.na(a$chr)) {
      pad <- as.numeric(Sys.getenv("C2_CIS_WINDOW_BP", unset = "1000000"))
      iv <- iv |> mutate(analysis = ifelse(CHR == a$chr & POS >= a$start - pad & POS <= a$end + pad, "cis", "trans"))
    } else iv <- iv |> mutate(analysis = ifelse(row_number() == which.min(P), "cis", "trans"))
  } else {
    # A metabolite has no gene-defined cis region. "local" is the lead locus and "distal" the other independent loci.
    lead <- iv |> arrange(P) |> slice(1)
    iv <- iv |> mutate(analysis = ifelse(CHR == lead$CHR & abs(POS - lead$POS) <= as.numeric(Sys.getenv("C2_LOCAL_WINDOW_BP", unset = "1000000")), "local", "distal"))
  }
  list(instruments = iv, files = fs)
}
is_palindromic <- function(a1, a2) paste0(a1, a2) %in% c("AT","TA","CG","GC")
allele_complement <- function(a) chartr("ATCG", "TAGC", toupper(a))
harmonize_sumstats <- function(x, y, drop_palindromic = TRUE) {
  d <- inner_join(x |> select(SNP, CHR_x = CHR, POS_x = POS, EA_x = EA, NEA_x = NEA, EAF_x = EAF, BETA_x = BETA, SE_x = SE, P_x = P, N_x = N, everything()),
                  y |> select(SNP, CHR_y = CHR, POS_y = POS, EA_y = EA, NEA_y = NEA, EAF_y = EAF, BETA_y = BETA, SE_y = SE, P_y = P, N_y = N), by = "SNP", suffix = c("", ".dup"))
  if (!nrow(d)) return(d)
  d <- d |> mutate(EA_x = toupper(EA_x), NEA_x = toupper(NEA_x), EA_y = toupper(EA_y), NEA_y = toupper(NEA_y),
                    EA_y_comp = allele_complement(EA_y), NEA_y_comp = allele_complement(NEA_y),
                    same = EA_x == EA_y & NEA_x == NEA_y,
                    flip = EA_x == NEA_y & NEA_x == EA_y,
                    strand_same = EA_x == EA_y_comp & NEA_x == NEA_y_comp,
                    strand_flip = EA_x == NEA_y_comp & NEA_x == EA_y_comp,
                    pal = is_palindromic(EA_x, NEA_x),
                    pal_same_diff = abs(EAF_x - EAF_y),
                    pal_flip_diff = abs(EAF_x - (1 - EAF_y)),
                    pal_freq_ok = pal & is.finite(EAF_x) & is.finite(EAF_y) &
                      pmin(pal_same_diff, pal_flip_diff) <= .10 & abs(pal_same_diff - pal_flip_diff) >= .05,
                    pal_reverse = pal_freq_ok & pal_flip_diff < pal_same_diff,
                    nonpal_keep = !pal & (same | flip | strand_same | strand_flip),
                    keep = if (drop_palindromic) nonpal_keep | pal_freq_ok else (same | flip | strand_same | strand_flip),
                    reverse = ifelse(pal, pal_reverse, flip | strand_flip),
                    BETA_y = ifelse(reverse, -BETA_y, BETA_y), EAF_y = ifelse(reverse, 1 - EAF_y, EAF_y),
                    harmonization = case_when(pal & pal_freq_ok & !pal_reverse ~ "palindrome_frequency_same",
                                              pal & pal_freq_ok & pal_reverse ~ "palindrome_frequency_swapped",
                                              !pal & same ~ "same", !pal & flip ~ "swapped",
                                              !pal & strand_same ~ "strand", !pal & strand_flip ~ "strand_swapped",
                                              TRUE ~ "mismatch")) |>
    filter(keep) |> distinct(SNP, .keep_all = TRUE) |>
    select(-EA_y_comp, -NEA_y_comp, -same, -flip, -strand_same, -strand_flip, -pal,
           -pal_same_diff, -pal_flip_diff, -pal_freq_ok, -pal_reverse, -nonpal_keep, -keep, -reverse)
  d
}
calc_iv_metrics <- function(iv) {
  if (!nrow(iv)) return(tibble(n_iv = 0L, r2_sum = 0, r2_bounded = 0, mean_F = NA_real_, min_F = NA_real_, max_F = NA_real_))
  # Partial R2 derived from each SNP's t statistic is bounded and does not depend on phenotype scaling.
  t2 <- (iv$BETA / iv$SE)^2; df <- pmax(iv$N - 2, 1)
  r2_i <- pmin(pmax(t2 / (t2 + df), 0), .999999)
  F_i <- t2
  tibble(n_iv = nrow(iv), r2_sum = sum(r2_i, na.rm = TRUE),
         r2_bounded = 1 - prod(1 - r2_i[is.finite(r2_i)]),
         mean_F = mean(F_i, na.rm = TRUE), min_F = min(F_i, na.rm = TRUE), max_F = max(F_i, na.rm = TRUE))
}
run_mr <- function(iv, ygwas, exposure, analysis) {
  empty <- tibble(exposure = exposure, analysis = analysis, method = NA_character_, n_IV = 0L,
                  b = NA_real_, se = NA_real_, pval = NA_real_, Q = NA_real_, Q_p = NA_real_,
                  egger_intercept = NA_real_, egger_intercept_p = NA_real_)
  if (!nrow(iv)) return(bind_cols(empty, calc_iv_metrics(iv)[, -1, drop = FALSE]))
  d <- harmonize_sumstats(iv, ygwas)
  if (!nrow(d)) return(bind_cols(empty, calc_iv_metrics(iv)[, -1, drop = FALSE]))
  d <- d |> filter(is.finite(BETA_x), BETA_x != 0, is.finite(BETA_y), is.finite(SE_y), SE_y > 0) |>
    mutate(ratio = BETA_y / BETA_x, ratio_se = abs(SE_y / BETA_x), w = 1 / ratio_se^2)
  if (!nrow(d)) return(bind_cols(empty, calc_iv_metrics(iv)[, -1, drop = FALSE]))
  if (nrow(d) == 1) {
    b <- d$ratio[[1]]; se <- d$ratio_se[[1]]; method <- "Wald ratio"; Q <- Qp <- NA_real_
  } else {
    b <- sum(d$w * d$ratio) / sum(d$w); se <- sqrt(1 / sum(d$w)); method <- "IVW fixed effect"
    Q <- sum(d$w * (d$ratio - b)^2); Qp <- pchisq(Q, df = nrow(d) - 1, lower.tail = FALSE)
  }
  ei <- eip <- NA_real_
  if (nrow(d) >= 3) {
    eg <- tryCatch(lm(BETA_y ~ BETA_x, weights = 1 / SE_y^2, data = d), error = function(e) NULL)
    if (!is.null(eg)) { sm <- coef(summary(eg)); ei <- sm["(Intercept)", "Estimate"]; eip <- sm["(Intercept)", "Pr(>|t|)"] }
  }
  metrics <- calc_iv_metrics(iv |> filter(SNP %in% d$SNP))
  bind_cols(tibble(exposure = exposure, analysis = analysis, method = method, n_IV = nrow(d), b = b, se = se,
                   pval = 2 * pnorm(abs(b / se), lower.tail = FALSE), Q = Q, Q_p = Qp,
                   egger_intercept = ei, egger_intercept_p = eip), metrics |> select(-n_iv))
}

# -----------------------------------------------------------------------------
# Reusable plots
# -----------------------------------------------------------------------------
plot_pwas_manhattan <- function(res, proteins, title = "Proteome-wide association") {
  bed <- read_ppp_bed(proteins)
  d <- res |> mutate(beta = coalesce(beta, safe_log(estimate))) |>
    transmute(protein = term, beta, p.value) |> left_join(bed, by = "protein") |>
    filter(!is.na(chr), is.finite(pos), is.finite(p.value)) |>
    mutate(chr_num = suppressWarnings(as.numeric(chr))) |> arrange(chr_num, pos) |> mutate(index = row_number(), sig = p.value < .05 / n(), y = -log10(pmax(p.value, 1e-300)))
  labs <- d |> arrange(p.value) |> slice_head(n = 25)
  ggplot(d, aes(index, y)) + geom_hline(yintercept = -log10(.05 / max(1, nrow(d))), linetype = 2, color = "#B2182B") +
    geom_point(aes(color = sig), alpha = .85, size = 1.5) +
    ggrepel::geom_text_repel(data = labs, aes(label = protein), size = 2.6, seed = 1, max.overlaps = Inf) +
    scale_color_manual(values = c(`TRUE` = "#2C7FB8", `FALSE` = "grey78"), guide = "none") +
    labs(title = title, x = "Protein genomic order", y = expression(-log[10](P))) + theme_5c(12)
}
plot_volcano <- function(res, xvar = "beta", label_col = "term", title = "Volcano", pth = NULL, top_n = 25) {
  if (is.null(pth)) pth <- .05 / max(1, nrow(res))
  d <- res |> mutate(x = .data[[xvar]], y = -log10(pmax(p.value, 1e-300)),
                     direction = case_when(p.value < pth & x > 0 ~ "Positive", p.value < pth & x < 0 ~ "Inverse", TRUE ~ "NS"),
                     label = ifelse(min_rank(p.value) <= top_n, .data[[label_col]], NA_character_))
  ggplot(d, aes(x, y)) + geom_hline(yintercept = -log10(pth), linetype = 2, color = "grey45") +
    geom_vline(xintercept = 0, color = "grey65") + geom_point(aes(color = direction), alpha = .82, size = 1.8) +
    ggrepel::geom_text_repel(aes(label = label), size = 2.7, seed = 2, max.overlaps = Inf, na.rm = TRUE) +
    scale_color_manual(values = c(Positive = "#D7301F", Inverse = "#2C7FB8", NS = "grey82")) +
    labs(title = title, x = "Effect estimate", y = expression(-log[10](P)), color = NULL) + theme_5c(12) + theme(legend.position = "bottom")
}
plot_met_circle <- function(d, title, label_n = 60, fdr_cap = NULL) {
  d <- d |> filter(is.finite(beta), is.finite(p.value)) |> mutate(FDR = coalesce(FDR, p.adjust(p.value, "BH")),
      group = coalesce(na_if(group, ""), "Other"), label = coalesce(label, term),
      group_label = str_replace_all(group, "_", " ")) |> arrange(group, FDR)
  if (!nrow(d)) return(blank_plot(title))
  # Compute BH-adjusted significance on the log scale. The ordinary P values
  # underflow to zero for very large Wald statistics, which previously forced
  # many metabolite bars to the identical 1e-300/cap radius.
  logp <- log(2) + pnorm(abs(d$statistic), lower.tail = FALSE, log.p = TRUE)
  ord <- order(logp); m <- length(logp); cand <- logp[ord] + log(m / seq_len(m))
  logq_ord <- rev(cummin(rev(cand))); logq <- numeric(m); logq[ord] <- pmin(logq_ord, 0)
  d$neglog10_FDR <- -logq / log(10)
  if (is.null(fdr_cap)) fdr_cap <- as.numeric(quantile(d$neglog10_FDR, .99, na.rm = TRUE, names = FALSE))
  if (!is.finite(fdr_cap) || fdr_cap <= 0) fdr_cap <- max(d$neglog10_FDR, na.rm = TRUE)
  gap <- d[1, , drop = FALSE]; gap[,] <- NA; gap$gap <- TRUE
  z <- d |> mutate(gap = FALSE) |> group_split(group) |> map_dfr(~bind_rows(.x, gap)) |>
    mutate(id = row_number(), angle0 = 90 - 360 * (id - .5) / n(), hjust = ifelse(angle0 < -90, 1, 0),
           angle = ifelse(angle0 < -90, angle0 + 180, angle0),
           height = pmin(neglog10_FDR, fdr_cap),
           fillv = beta, lab = ifelse(!gap & min_rank(FDR) <= label_n, str_trunc(label, 24), NA_character_))
  bands <- z |> filter(!gap) |> group_by(group, group_label) |>
    summarise(x1 = min(id)-.4, x2 = max(id)+.4, x=(x1+x2)/2, .groups="drop") |>
    mutate(angle0 = 90 - 360 * (x - .5) / nrow(z), hjust = ifelse(angle0 < -90, 1, 0),
           angle = ifelse(angle0 < -90, angle0 + 180, angle0))
  ymax <- max(z$height, na.rm = TRUE)
  beta_lim <- as.numeric(quantile(abs(z$fillv), .98, na.rm = TRUE))
  if (!is.finite(beta_lim) || beta_lim <= 0) beta_lim <- max(abs(z$fillv), na.rm = TRUE)
  ggplot(z, aes(id, height)) + geom_col(aes(fill = fillv), width = .9, na.rm = TRUE) +
    geom_segment(data = bands, aes(x=x1,xend=x2,y=-ymax*.10,yend=-ymax*.10), inherit.aes=FALSE, linewidth=.8) +
    geom_text(data=bands, aes(x=x,y=-ymax*.46,label=str_trunc(group_label, 24),angle=angle,hjust=hjust),
              inherit.aes=FALSE, size=1.65, fontface="bold") +
    geom_text(aes(y=height+ymax*.035,label=lab,angle=angle,hjust=hjust), size=1.75, na.rm=TRUE) +
    coord_polar(clip="off") +
    scale_fill_gradient2(low="#0878BF",mid="white",high="#FF5A52",midpoint=0,
                         limits=c(-beta_lim,beta_lim),oob=scales::squish,name="ln(HR)") +
    ylim(-ymax*.82, ymax*1.27) + labs(title=title, subtitle="Radius: -log10(FDR), winsorized at 99th percentile; colour: ln(HR)") + theme_void(11) +
    theme(plot.title=element_text(face="bold",hjust=.5),legend.position="bottom",plot.margin=margin(25,50,25,50))
}

# Find PRS columns without reading a second file. Preference is CAD, then BMI, then other score_sum variables.
find_prs_vars <- function(dat, outcome = Y, max_n = 4) {
  z <- grep("score_sum$|\\.prs$|_prs$", names(dat), value = TRUE, ignore.case = TRUE)
  if (!length(z)) return(character())
  key <- c(outcome, str_remove(outcome, "^cvd_"), "cad", "bmi")
  ord <- order(!vapply(z, function(v) any(str_detect(tolower(v), fixed(tolower(key)))), logical(1)), z)
  head(z[ord], max_n)
}
