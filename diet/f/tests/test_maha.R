# Base-R regression checks for MAHA timing, score profiles and publication phases.
entry <- sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[[1L]])
maha_dir <- file.path(dirname(dirname(normalizePath(entry))), "maha")
get_functions <- function(path) {
    out <- list()
    for (e in parse(path)) if (is.call(e) && identical(e[[1L]], as.name("<-")) &&
        is.symbol(e[[2L]]) && is.call(e[[3L]]) && identical(e[[3L]][[1L]], as.name("function"))) {
        out[[as.character(e[[2L]])]] <- e[[3L]]
    }
    out
}
functions <- get_functions(file.path(maha_dir, "ukb_init.R"))
env <- new.env(parent = baseenv())
for (nm in c("maha_date", "maha_webq_window", "maha_set_landmark", "maha_outcome_followup")) assign(nm, eval(functions[[nm]], env), env)
raw <- data.frame(eid = as.character(1:4), p105010_i0 = as.Date(c("2010-01-01", "2010-01-01", NA, "2010-01-01")))
for (k in 1:4) raw[[paste0("p105010_i", k)]] <- as.Date(NA)
raw$p105010_i1[1L] <- as.Date("2012-01-01")
raw$p26002_i1 <- c(2000, 2000, NA, NA)
window <- env$maha_webq_window(raw)
stopifnot(window$webq_window_end[1L] == as.Date("2012-01-01"), window$webq_undated[2L], is.na(window$webq_window_end[3L]))
dat <- data.frame(eid = as.character(1:4), diet.maha.s100 = 50, age = 50,
    date_attend = as.Date("2009-01-01"), date_lost = as.Date(NA), date_death = as.Date(NA))
dat <- env$maha_set_landmark(dat, window, "diet.maha.s100")
dat$fod_icd10_cad <- as.Date(c("2012-01-01", "2013-01-01", NA, "2013-01-01"))
dat$cnt_icd10_cad <- c(1, 1, 0, 1)
out <- env$maha_outcome_followup(dat, "cad", "2015-01-01")
stopifnot(is.na(out$cad.Yt2e[1L]), is.na(out$cad.Yt2e[2L]), is.na(out$cad.Yt2e[3L]), out$cad.Yt2e[4L] == 1L)
dat$fod_icd10_cad[4L] <- as.Date(NA)
out <- env$maha_outcome_followup(dat, "cad", "2015-01-01")
stopifnot(is.na(out$cad.Yt2e[4L]))
functions <- get_functions(file.path(maha_dir, "common.R"))
bf01 <- eval(functions$ap_bf01_one)
stopifnot(bf01(0, .03, 1.10) > 1, bf01(.2, .03, 1.10) < 1,
    abs(bf01(.2, .03, 1.10) - bf01(-.2, .03, 1.10)) < 1e-12,
    is.na(bf01(0, 0, 1.10)))

check_phase_contract <- function() {
    root <- tempfile("maha-contract-")
    dir.create(root)
    on.exit(unlink(root, recursive = TRUE), add = TRUE)
    e <- new.env(parent = globalenv())
    source(file.path(dirname(maha_dir), "common.R"), local = e)
    a <- data.frame(eid = c("001", "003"), diet.dash.sum = c(30, 34), diet.maha.sum = c(4, 6))
    b <- data.frame(eid = c("002", "003"), diet.dash.sum = c(28, 31), dose9half = c(4.5, 5.5))
    saveRDS(a, file.path(root, "maha.rds")); saveRDS(b, file.path(root, "dose.rds"))
    e$diet_publish_disco("ukb", root, c(maha = file.path(root, "maha.rds"), dose = file.path(root, "dose.rds")))
    e$maha_cohort <- "ukb"
    e$maha_score_file <- file.path(root, "disco", "ukb.disco.rds")
    functions <- get_functions(file.path(maha_dir, "run.R"))
    for (nm in c("maha_score_frame", "maha_join_scores")) assign(nm, eval(functions[[nm]], e), e)
    joined <- e$maha_join_scores(data.frame(eid = c("003", "002", "001")), "eid")
    stopifnot(identical(joined$eid, c("003", "002", "001")),
        identical(joined$diet.dash.sum, c(34, NA_real_, 30)), !"dose9half" %in% names(joined))
    e$diet_publish_disco("ukb", root, c(dose = file.path(root, "dose.rds")))
    error <- tryCatch(e$maha_join_scores(data.frame(eid = "003"), "eid"), error = conditionMessage)
    stopifnot(is.character(error), grepl("do not contain profile maha", error, fixed = TRUE))

    # A new score file requires a new phenotype snapshot; old caches fail clearly.
    for (nm in c("maha_state_file", "maha_sweep_child", "maha_disco_stamp", "maha_write_rds",
        "maha_save_state", "maha_load_state")) assign(nm, eval(functions[[nm]], e), e)
    e$maha_state_dir <- file.path(root, "state")
    e$maha_save_state("phe", character(), envir = e)
    e$maha_load_state("phe")
    saveRDS(data.frame(eid = "003", diet.maha.sum = 5), e$maha_score_file)
    error <- tryCatch(e$maha_load_state("phe"), error = conditionMessage)
    stopifnot(grepl("changed after phenotype preparation", error, fixed = TRUE))
    saveRDS(list(), e$maha_state_file("phe"))
    error <- tryCatch(e$maha_load_state("phe"), error = conditionMessage)
    stopifnot(grepl("predates the canonical score link", error, fixed = TRUE))
    previous_child <- Sys.getenv("CHNS_MORTALITY_SWEEP_CHILD", unset = NA_character_)
    on.exit(if (is.na(previous_child)) Sys.unsetenv("CHNS_MORTALITY_SWEEP_CHILD") else
        Sys.setenv(CHNS_MORTALITY_SWEEP_CHILD = previous_child), add = TRUE)
    Sys.setenv(CHNS_MORTALITY_SWEEP_CHILD = "1")
    e$maha_cohort <- "chns"; e$maha_score_file <- file.path(root, "absent-parent-baseline.rds")
    e$maha_save_state("phe", character(), envir = e)
    e$maha_load_state("phe")
    stopifnot(is.null(attr(readRDS(e$maha_state_file("phe")), "disco_md5")))

    functions <- get_functions(file.path(maha_dir, "publication.R"))
    for (nm in c("clear_publication_outputs", "render_one", "merge_publication_status")) assign(nm, eval(functions[[nm]], e), e)
    e$outdir <- file.path(root, "final"); dir.create(e$outdir)
    owned <- c("Fig3.nhanes.mortality.png", "Fig3.nhanes.mortality.out.xlsx")
    retained <- c("Fig2.ukb.precision.png", "nhanes__Fig3.mortality.png", "notes.nhanes.txt")
    for (name in c(owned, retained)) writeLines("fixture", file.path(e$outdir, name))
    e$publication_status <- list()
    evalq(render_one("nhanes", stop("expected missing input")), e)
    stopifnot(!any(file.exists(file.path(e$outdir, owned))),
        all(file.exists(file.path(e$outdir, retained))),
        e$publication_status$nhanes$status == "pending_or_partial")
    previous <- data.frame(item = c("ukb", "nhanes"), status = "complete", reason = "")
    current <- data.frame(item = "nhanes", status = "pending_or_partial", reason = "Missing input")
    statuses <- e$merge_publication_status(current, previous)
    stopifnot(nrow(statuses) == 2L, statuses$status[statuses$item == "ukb"] == "complete",
        statuses$status[statuses$item == "nhanes"] == "pending_or_partial")
}
check_phase_contract()
check_score_refresh <- function() {
    if (!requireNamespace("dplyr", quietly = TRUE)) stop("The canonical-scale regression check requires dplyr")
    e <- new.env(parent = globalenv())
    functions <- get_functions(file.path(maha_dir, "nhanes_init.R"))
    for (nm in c("std_10_90", "score_0_100", "score_q5", "f3c", "mk_hml", "nhanes_refresh_scales"))
        assign(nm, eval(functions[[nm]], e), e)
    e$case_when <- dplyr::case_when
    e$diet.lst <- c(maha = "MAHA"); e$group_pct_th <- .4
    original <- e$nhanes_refresh_scales(data.frame(SEQN = 1:10, diet.maha.sum = 1:10))
    updated <- original; updated$diet.maha.sum[1L] <- 20
    updated <- e$nhanes_refresh_scales(updated)
    stopifnot(updated$diet.maha.pts[1L] > original$diet.maha.pts[1L],
        identical(updated$diet.maha.pts, e$std_10_90(updated$diet.maha.sum)),
        original$diet.maha.3c[1L] == "low", updated$diet.maha.3c[1L] == "high",
        updated$diet.maha.hml[1L] == "high", updated$diet.maha.q5[1L] == 100,
        identical(e$nhanes_refresh_scales(original), original))
    functions <- get_functions(file.path(maha_dir, "chns_init.R"))
    for (nm in c("num", "score100", "f3c", "chns_refresh_scales")) assign(nm, eval(functions[[nm]], e), e)
    original <- e$chns_refresh_scales(data.frame(IDIND = 1:10, diet.maha.sum = 1:10), "diet.maha.sum")
    updated <- original; updated$diet.maha.sum[1L] <- 20
    updated <- e$chns_refresh_scales(updated, "diet.maha.sum")
    stopifnot(original$diet.maha.s100[1L] == 0, updated$diet.maha.s100[1L] == 100,
        updated$diet.maha.3c[1L] == "high", identical(e$chns_refresh_scales(original, "diet.maha.sum"), original))
}
check_score_refresh()
cat("MAHA tests passed: recall landmark, prevalence boundary, BF01 direction, score profile isolation, derived-scale refresh, phase freshness, sweep isolation and scoped publication refresh.\n")
