#!/usr/bin/env Rscript
# Deliberate setup only. No analysis script sources this installer.
# Existing requested R packages are retained; there is no blanket update step.
args <- commandArgs(trailingOnly=TRUE)
allowed <- c("--check", "--with-optional")
if (length(setdiff(args, allowed))) stop("Unknown installer option: ", paste(setdiff(args, allowed), collapse=", "))
optional <- "--with-optional" %in% args
if (getRversion() < "4.5.0") stop("Use R >= 4.5; env.yml selects R 4.5.")

# Direct imports, dynamic package lists, and packages needed by the shared 0f
# functions actually used by the MAHA/DOSE workflows. Base/recommended packages
# are checked but not installed from CRAN here.
recommended <- c("MASS", "Matrix", "lattice", "survival")
cran <- c(
    "tidyverse", "dplyr", "tidyr", "readr", "purrr", "tibble", "stringr", "forcats", "ggplot2",
    "data.table", "R.utils", "readxl", "writexl", "openxlsx", "haven", "jsonlite", "digest",
    "survey", "broom", "rms", "pROC", "irr", "matrixStats", "Rcpp",
    "patchwork", "cowplot", "ggrepel", "scales", "flextable", "gtsummary",
    "ggnewscale", "corrplot", "pheatmap", "svglite",
    "mediation", "lavaan", "igraph", "sandwich",
    "sf", "raster", "rosm", "ggspatial", "rnaturalearth", "rnaturalearthdata"
)
bioc <- "impute"  # Shared phenotype.R clean_biom()/drop_10pct() when imputation is selected.
if (optional) {
    # Optional SEM diagram and other lazy-loaded 0f plotting/prediction helpers.
    cran <- c(cran, "semPlot", "circlize", "glmnet", "randomForestSRC", "xgboost")
    bioc <- c(bioc, "ComplexHeatmap")
}
github <- c(
    PheWAS="PheWAS/PheWAS@eae96fa5ed3f42e6dc5d22fad09fd7734f0df48f",
    dietaryindex="jamesjiadazhan/dietaryindex@9830a68bcfdc16204b7c2a85668548cb0b5f22d8"
)
expected_version <- c(PheWAS="0.99.6-2", dietaryindex="2.1")
cran <- unique(cran)
available <- function(x) vapply(x, requireNamespace, logical(1), quietly=TRUE)

required <- unique(c(recommended, cran, bioc, names(github)))
absent <- required[!available(required)]
cat(sprintf("R %s: %d/%d required packages loadable (%s).\n", getRversion(),
    length(required)-length(absent), length(required), if (optional) "core + optional" else "core"))
if (length(absent)) cat("Missing/unloadable: ", paste(absent, collapse=", "), "\n", sep="")
for (pkg in setdiff(names(github), absent)) {
    desc <- utils::packageDescription(pkg)
    actual <- desc$RemoteSha
    expected <- sub("^.*@", "", github[[pkg]])
    status <- if (is.null(actual)) "revision metadata unavailable" else if (identical(actual, expected)) "pinned revision" else "different revision retained"
    cat(sprintf("  %s %s: %s (reference version %s).\n", pkg, desc$Version, status, expected_version[[pkg]]))
    if (is.null(actual) || !identical(actual, expected)) cat("    Reference: ", github[[pkg]], "\n", sep="")
}
if (length(absent)) {
    cat("Prepare missing packages separately before running.\n")
    quit(save="no", status=1L)
}
