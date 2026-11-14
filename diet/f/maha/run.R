#!/usr/bin/env Rscript
# Internal cohort engine. Public entry points are 1.disco.sh and diet.sh.
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2L) stop("Usage: Rscript f/maha/run.R COHORT disco|phe|assoc|bayes|final")
maha_cohort <- tolower(args[[1L]])
maha_phase <- tolower(args[[2L]])
if (!maha_cohort %in% c("ukb", "nhanes", "chns", "all")) stop("Unsupported MAHA cohort: ", maha_cohort)
if (!maha_phase %in% c("disco", "phe", "assoc", "bayes", "final", "all")) stop("Unsupported MAHA phase: ", maha_phase)
if (maha_cohort == "all" && maha_phase != "final") stop("Select one cohort for analysis phases")
entry <- sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[[1L]])
maha_runner <- normalizePath(entry, winslash = "/", mustWork = TRUE)
maha_code <- dirname(maha_runner)
helper <- Sys.getenv("DIET_HELPERS", Sys.getenv("MAHA_HELPER_DIR", ""))
if (!nzchar(helper)) {
    directory <- maha_code
    repeat {
        helper <- file.path(directory, "0f")
        if (file.exists(file.path(helper, "0os.R"))) break
        parent <- dirname(directory)
        if (identical(parent, directory)) stop("Set DIET_HELPERS to the scripts/0f directory")
        directory <- parent
    }
}
source(file.path(helper, "0os.R"))
source(file.path(dirname(maha_code), "common.R"))
os_paths_init(default_mode = "local")
maha_output <- Sys.getenv("DIET_OUT", Sys.getenv("DIET_OUTDIR", file.path(os_path("ANALYSIS_ROOT"), "diet")))
maha_output <- normalizePath(maha_output, winslash = "/", mustWork = FALSE)
maha_cache <- Sys.getenv("DIET_CACHE", file.path(maha_output, ".cache"))
maha_state_dir <- Sys.getenv("DIET_MAHA_STATE", file.path(maha_cache, "maha", maha_cohort))
maha_scores_dir <- Sys.getenv("DIET_SCORES", file.path(maha_output, "disco"))
maha_score_file <- Sys.getenv("DIET_DISCO_FILE", file.path(maha_scores_dir, paste0(maha_cohort, ".disco.rds")))
maha_part_file <- Sys.getenv("DIET_DISCO_PART", file.path(maha_cache, "disco", paste0(maha_cohort, ".maha.rds")))
dir.create(maha_state_dir, recursive = TRUE, showWarnings = FALSE)
Sys.setenv(MAHA_HELPER_DIR = helper,
           MAHA_OUTDIR = Sys.getenv("MAHA_OUTDIR", file.path(maha_output, "maha")),
           MAHA_AUXDIR = Sys.getenv("MAHA_AUXDIR", file.path(maha_output, "maha", maha_cohort)),
           MAHA_JOBS = Sys.getenv("DIET_THREADS", Sys.getenv("MAHA_JOBS", "1")))
# Scoring diagnostics are cache artifacts; public MAHA directories belong to
# downstream phenotype/association/publication phases only.
if (identical(maha_phase, "disco")) {
    scoring_aux <- file.path(maha_state_dir, "diagnostics")
    Sys.setenv(MAHA_OUTDIR = scoring_aux, MAHA_AUXDIR = scoring_aux)
}
if (tolower(Sys.getenv("DIET_REPLACE", "TRUE")) %in% c("true", "1", "yes")) Sys.setenv(MAHA_REPLACE_STEPS = "all")
if (!nzchar(Sys.getenv("DIET_UKB_RDATA"))) Sys.setenv(DIET_UKB_RDATA = file.path(maha_cache, "ukb", "maha", "Rdata"))

maha_load_packages <- function(...) {
    packages <- vapply(as.list(substitute(list(...)))[-1L], function(x) {
        if (is.character(x)) x else as.character(x)
    }, character(1))
    missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
    if (length(missing)) stop("Missing R packages: ", paste(missing, collapse = ", "),
        ". Run diet/install.sh before analysis.")
    invisible(lapply(packages, function(x) suppressPackageStartupMessages(library(x, character.only = TRUE))))
}

maha_write_rds <- function(x, path) {
    dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
    tmp <- tempfile(pattern = ".diet-", tmpdir = dirname(path), fileext = ".rds")
    on.exit(unlink(tmp), add = TRUE)
    saveRDS(x, tmp, version = 3)
    if (file.exists(path) && .Platform$OS.type == "windows") unlink(path)
    if (!file.rename(tmp, path)) stop("Cannot replace ", path)
    invisible(path)
}
maha_state_file <- function(stage) file.path(maha_state_dir, paste0(stage, ".rds"))
maha_sweep_child <- function() identical(maha_cohort, "chns") &&
    identical(Sys.getenv("CHNS_MORTALITY_SWEEP_CHILD"), "1")
maha_disco_stamp <- function() {
    if (!file.exists(maha_score_file)) stop("Missing canonical diet scores: ", maha_score_file,
        ". Run 1.disco.sh ", maha_cohort, " first.")
    stamp <- unname(tools::md5sum(maha_score_file))
    if (length(stamp) != 1L || is.na(stamp)) stop("Cannot read canonical diet scores: ", maha_score_file)
    stamp
}
maha_save_state <- function(stage, vars, envir = .GlobalEnv) {
    vars <- intersect(unique(vars), ls(envir = envir, all.names = TRUE))
    values <- mget(vars, envir = envir, inherits = FALSE)
    values <- values[!vapply(values, function(x) is.function(x) || is.environment(x) || inherits(x, "connection"), logical(1))]
    # Clinical exclusions/covariates belong to this exact published score file.
    # CHNS sweep children use isolated wave-native scores instead of the parent
    # cohort's baseline canonical file and therefore do not carry this stamp.
    if (stage == "phe" && !maha_sweep_child()) attr(values, "disco_md5") <- maha_disco_stamp()
    maha_write_rds(values, maha_state_file(stage))
}
maha_save_required_state <- function(stage, files) {
    # Keep the values actually referenced by later phase files. Functions and paths
    # are reloaded from code; raw modules not used downstream are not duplicated.
    symbols <- unique(unlist(lapply(files, function(f) all.names(parse(file.path(maha_code, f)), functions = FALSE))))
    symbols <- unique(c("dat", "baseline_wave", "followup_end", "common_baseline", symbols))
    reserved <- c("maha_output", "maha_state_dir", "maha_cache", "maha_code", "maha_runner", "maha_part_file", "maha_score_file", "maha_cohort", "maha_phase", "maha_auxdir", "outdir", "dir0", "project_dir", "raw_dir", "log_con", "maha_jobs", "null_max", "run_mortality_sweep", "mortality_sweep_child", "mortality_sweep_wave_req", "mortality_model_age", "mortality_event_target")
    maha_save_state(stage, setdiff(symbols, reserved))
}
maha_load_state <- function(stage) {
    path <- maha_state_file(stage)
    if (!file.exists(path)) stop("Missing MAHA ", stage, " data: ", path,
        ". Run ", if (stage == "score") "1.disco.sh" else "diet.sh", " ", maha_cohort, " first.")
    values <- readRDS(path)
    if (!is.list(values) || is.data.frame(values)) stop("Invalid phase data in ", path)
    if (stage == "phe" && !maha_sweep_child()) {
        expected <- attr(values, "disco_md5", exact = TRUE)
        instruction <- paste0("Rerun diet.sh ", maha_cohort, " --study maha before association or Bayes analysis.")
        if (is.null(expected)) stop("MAHA phenotype cache predates the canonical score link. ", instruction)
        if (!identical(expected, maha_disco_stamp()))
            stop("Canonical ", maha_cohort, " diet scores changed after phenotype preparation. ", instruction)
    }
    list2env(values, envir = .GlobalEnv)
    invisible(NULL)
}
maha_score_frame <- function(x) {
    if (is.data.frame(x)) return(x)
    if (is.list(x) && is.data.frame(x$scores)) return(x$scores)
    if (is.list(x) && is.data.frame(x$data)) return(x$data)
    stop("Diet score file must contain a data.frame")
}
maha_export_scores <- function(dat, id, wave = NULL) {
    keep <- unique(c(id, if (maha_cohort == "nhanes") "cycle", grep("^diet\\..*\\.sum$|^(maha|dash|medi|mind)(_.*)?_c_|\\.maha(_[a-z]+)?$", names(dat), value = TRUE)))
    keep <- intersect(keep, names(dat))
    out <- as.data.frame(dat[, keep, drop = FALSE])
    out[[id]] <- as.character(out[[id]])
    if (anyNA(out[[id]]) || anyDuplicated(out[[id]])) stop("Missing or duplicate ", id, " in MAHA score output")
    attr(out, "cohort") <- maha_cohort
    attr(out, "profile") <- "maha"
    attr(out, "source_commit") <- "971a61b44b84500ab46e2f5592478505a0eae50e"
    attr(out, "baseline_wave") <- wave
    attr(out, "method_note") <- "Cohort-specific source scoring and eligibility; raw scales are not interchangeable across cohorts."
    if (maha_cohort == "nhanes") {
        attr(out, "main_cycles") <- main_cycles
        attr(out, "mort_cycles") <- mort_cycles
        attr(out, "record_unit") <- "Survey-release record (SEQN); 2017-2018 and 2017-March2020 overlap with reissued identifiers. Do not pool both releases in one analysis."
    }
    maha_write_rds(out, maha_part_file)
    message("Wrote MAHA score part: ", maha_part_file, " (", nrow(out),
        if (maha_cohort == "nhanes") " survey-release records)" else " participants)")
    invisible(out)
}
maha_join_scores <- function(dat, id) {
    # Resolve the recorded MAHA columns and membership even when another
    # profile has an identically named score. A DOSE-only file is not MAHA data.
    scores <- diet_profile(maha_score_frame(readRDS(maha_score_file)), "maha", maha_cohort)
    if (!id %in% names(scores)) {
        alternate <- intersect(c("eid", "SEQN", "IDIND", "id"), names(scores))
        if (length(alternate) != 1L) stop("Cannot resolve identifier in ", maha_score_file)
        names(scores)[names(scores) == alternate] <- id
    }
    if (anyDuplicated(as.character(scores[[id]]))) stop("Duplicate identifiers in ", maha_score_file)
    idx <- match(as.character(dat[[id]]), as.character(scores[[id]]))
    vars <- setdiff(names(scores), id)
    vars <- vars[grepl("^diet\\.|^(maha|dash|medi|mind)(_.*)?_c_|\\.maha(_[a-z]+)?$", vars)]
    for (v in vars) dat[[v]] <- scores[[v]][idx]
    dat
}
maha_write_phenotypes <- function(dat, id, mortality = NULL) {
    if (identical(Sys.getenv("CHNS_MORTALITY_SWEEP_CHILD"), "1")) return(invisible(NULL))
    path <- file.path(maha_output, "phe", paste0(maha_cohort, ".maha.rds"))
    maha_write_rds(list(data = dat, mortality = mortality, id = id, profile = "maha"), path)
}
maha_write_assoc <- function(tbl, cohort) {
    if (identical(Sys.getenv("CHNS_MORTALITY_SWEEP_CHILD"), "1")) return(invisible(NULL))
    path <- file.path(maha_output, "assoc", paste0(cohort, ".maha.tsv"))
    dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
    utils::write.table(as.data.frame(tbl), path, sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
}
maha_write_bayes <- function(tbl, cohort) {
    if (identical(Sys.getenv("CHNS_MORTALITY_SWEEP_CHILD"), "1")) return(invisible(NULL))
    path <- file.path(maha_output, "bayes", paste0(cohort, ".maha.tsv"))
    dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
    utils::write.table(as.data.frame(tbl), path, sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
}
maha_ukb_score <- function() {
    root <- Sys.getenv("DIET_UKB_RDATA")
    dat <- readRDS(file.path(root, "diet.rds"))
    comp <- readRDS(file.path(root, "diet.maha.rds"))
    if (!"eid" %in% names(dat) || !"eid" %in% names(comp)) stop("UKB score input has no eid")
    idx <- match(as.character(dat$eid), as.character(comp$eid))
    for (nm in setdiff(names(comp), names(dat))) dat[[nm]] <- comp[[nm]][idx]
    maha_export_scores(dat, "eid")
}

log_path <- file.path(maha_state_dir, paste0(maha_phase, ".log"))
log_con <- file(log_path, "wt")
sink(log_con, split = TRUE)
status <- tryCatch({
    if (maha_phase == "final") {
        sys.source(file.path(maha_code, "publication.R"), envir = .GlobalEnv)
    } else if (maha_cohort == "ukb" && maha_phase == "disco") {
        maha_ukb_score()
    } else {
        sys.source(file.path(maha_code, paste0(maha_cohort, "_init.R")), envir = .GlobalEnv)
        .maha_selected_steps <- character()
        phases <- if (maha_phase == "all") c("disco", "phe", "assoc", "bayes") else maha_phase
        for (phase in phases) {
            if (phase == "disco" && maha_cohort == "ukb") maha_ukb_score() else {
                suffix <- if (phase == "disco") "score" else phase
                sys.source(file.path(maha_code, paste0(maha_cohort, "_", suffix, ".R")), envir = .GlobalEnv)
            }
        }
    }
    0L
}, error = function(e) { message(conditionMessage(e)); 1L })
while (sink.number() > 0L) sink()
close(log_con)
if (status) quit(save = "no", status = status)
message("MAHA ", maha_cohort, " ", maha_phase, " complete")
