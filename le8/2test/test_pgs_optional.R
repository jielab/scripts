# Base-R regression checks: disabling PGS must not inspect its raw input files.
args <- commandArgs(trailingOnly = TRUE)
root <- if (length(args)) normalizePath(args[1]) else normalizePath(".")
env <- new.env(parent = baseenv())
wanted <- c("find_c1_pgs_file", "c1_pgs_signature", "le8_c1_pgs_mode", "le8_c1_additions")
load_definitions <- function(expr) {
    if (missing(expr) || !is.call(expr)) return(invisible(NULL))
    if (identical(expr[[1]], as.name("<-")) && length(expr) == 3L &&
        is.symbol(expr[[2]]) && as.character(expr[[2]]) %in% wanted) {
        eval(expr, envir = env)
        return(invisible(NULL))
    }
    for (part in as.list(expr)[-1L]) load_definitions(part)
}
for (file in c("f/0.common.R", "f/c1.correlate.R")) {
    for (expr in parse(file.path(root, file))) load_definitions(expr)
}
stopifnot(all(vapply(wanted, exists, logical(1), envir = env, inherits = FALSE)))
scratch <- tempfile("le8-pgs-optional-")
dir.create(scratch)
Sys.setenv(DIR_DATA = scratch)
env$os_path <- function(...) scratch
for (mode in c("FALSE", "SKIP")) {
    Sys.setenv(C1_PGS_MODE = mode)
    for (layer in c("protein", "metabolite")) {
        stopifnot(identical(env$c1_pgs_signature(layer), "skipped_by_request"))
    }
}
Sys.setenv(C1_PGS_MODE = "TRUE")
for (layer in c("protein", "metabolite")) {
    err <- tryCatch(env$c1_pgs_signature(layer), error = identity)
    stopifnot(inherits(err, "error"), grepl("Required canonical PGS file", conditionMessage(err)))
}
dir.create(file.path(scratch, "gwas/prot/pgs"), recursive = TRUE)
file <- file.path(scratch, "gwas/prot/pgs/prot.pgs.gz")
writeLines("present", file)
stopifnot(grepl(normalizePath(file), env$c1_pgs_signature("protein"), fixed = TRUE))

# Exercise the paired-review branch with raw PGS access forbidden.
env$find_c1_pgs_file <- function(...) stop("Disabled PGS attempted raw file access")
env$le8_csv_env <- function(...) character()
env$le8_time_heterogeneity <- function(...) data.frame()
env$le8_job_dir <- function(...) scratch
env$tibble <- function(...) data.frame(...)
env$c1_saved_pgs_result <- function(...) list()
env$c1_saved_pgs_table <- function(...) NULL
env$`%||%` <- function(x, y) if (is.null(x)) y else x
env$write_raw_csv <- function(...) invisible(NULL)
env$le8_plot_c1_temporal <- function(...) invisible(NULL)
env$out.prot <- env$out.met <- scratch
Sys.setenv(C1_PGS_MODE = "FALSE")
for (layer in c("protein", "metabolite")) {
    result <- env$le8_c1_additions(data.frame(), layer, character(), "time", "event")
    stopifnot(nrow(result$paired_associations) == 0L,
              grepl("PGS disabled", result$paired_status$status, fixed = TRUE))
}
unlink(scratch, recursive = TRUE)
cat("PGS optional-input regression checks passed\n")
