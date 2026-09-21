# Rscript tests/test_c1_cohort_restore.R
suppressPackageStartupMessages(library(dplyr))
script <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))[1]
fdir <- file.path(dirname(normalizePath(script)), "..", "f")
# Load definitions without starting the analysis pipeline or reading UKB inputs.
load_function <- function(file, name) {
  for (expr in parse(file)) {
    if (is.call(expr) && identical(expr[[1]], as.name("<-")) &&
        identical(expr[[2]], as.name(name))) eval(expr, .GlobalEnv)
  }
}
`%||%` <- function(x, y) if (is.null(x)) y else x
load_function(file.path(fdir, "comm.f.R"), "add_attained_age_time")
load_function(file.path(fdir, "comm.f.R"), "write_raw_csv")
load_function(file.path(fdir, "c1_correlate.R"), "make_prevalent_status")
required_file <- function(path, ...) { stopifnot(file.exists(path)); path }
source(file.path(fdir, "le8_figures.R"))
vars.le8 <- c("diet.pts", "bmi.pts")
dat <- tibble(cad.Yt2e = c(1, 0, NA), cad.t2e = c(5, 10, NA),
  cad.bi2e = c(65, 70, 55), date_attend = as.Date(rep("2010-01-01", 3)),
  fod_icd10_cad = as.Date(c("2015-01-01", NA, "2009-01-01")),
  diet.pts = 1, bmi.pts = 2)
obj <- list(meta = list(trait = "cad", N = 3L, events = 1L,
  covs_use = "adj_le4", le4_covariates = "diet.pts",
  le8_covariates_removed = "bmi.pts", treatment_covariates = character(),
  pgs_matched = 2L, pgs_signature = "cached-signature"),
  input_feature_annotation_audit = tibble(feature = c("A", "B", "C", "D"),
    annotation_matched = c(TRUE, TRUE, FALSE, NA)))
rawdir <- tempfile(); dir.create(rawdir)
tryCatch({
  for (layer in c("protein", "metabolite")) {
    unlink(file.path(rawdir, "c1.cohort.csv"))
    cohort <- le8_figure_c1_cohort(obj, dat, layer, rawdir)
    stopifnot(cohort$layer == layer, cohort$N_omics == 3,
      cohort$incident_events == 1, cohort$prevalent_cases == 1,
      cohort$features == 4, cohort$annotation_matched == 2,
      cohort$annotation_unmatched == 1, cohort$attained_age_N == 2,
      cohort$bi2e_available == 3, cohort$LE8_components == 2,
      cohort$LE4_components == 1, cohort$PGS_file_signature == "cached-signature",
      file.exists(file.path(rawdir, "c1.cohort.csv")))
    # Existing CSV and embedded tables must work without individual-level data.
    from_csv <- le8_figure_c1_cohort(obj, NULL, layer, rawdir)
    stopifnot(from_csv$N_omics == 3, from_csv$features == 4)
    saved <- obj; saved$cohort <- cohort
    unlink(file.path(rawdir, "c1.cohort.csv"))
    stopifnot(identical(le8_figure_c1_cohort(saved, NULL, layer, rawdir), cohort))
  }
  changed <- obj; changed$meta$N <- 4L
  err <- tryCatch(le8_figure_c1_cohort(changed, dat, "protein", rawdir), error = identity)
  stopifnot(inherits(err, "error"), grepl("no longer match", conditionMessage(err)))
}, finally = unlink(rawdir, recursive = TRUE))
cat("C1 cohort restoration checks passed\n")
