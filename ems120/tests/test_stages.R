# Synthetic stage/cache checks: no patient data, CUDA or third-party packages.
test_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[[1]])
project <- dirname(dirname(normalizePath(test_file, mustWork = TRUE)))
Sys.setenv(EMS120_SCRIPT_DIR = project, EMS120_PHASE = "all")
source(file.path(project, "f", "runtime.R"))
source(file.path(project, "f", "stages.R"))

run_tests <- function() {
	root <- tempfile("ems120-stages-")
	dir.create(root)
	on.exit(unlink(root, recursive = TRUE), add = TRUE)
	Sys.setenv(EMS120_HANDOFF_DIR = file.path(root, "handoff"))
	cache <- file.path(root, "dat.list.rds")
	d <- list(`2019` = data.frame(score = c(1, 3, 7), dx.pipeline.version = "fixture-v1", transformations = 1L))
	saveRDS(d, cache)
	before <- tools::md5sum(cache)

	set.seed(71)
	summaries <- list(phone_score_summary = data.frame(score = 1:3, n = c(5L, 2L, 8L)))
	ems120_stage_publish("gpu_ml", cache, context = list(final_cached = FALSE, report_objects = summaries))
	expected_rng <- runif(5)
	set.seed(99)
	context <- ems120_stage_restore("gpu_ml")
	stopifnot(identical(runif(5), expected_rng), identical(context$report_objects, summaries))
	cat("PASS: cache handoff retains exact summaries and random stream\n")

	# Extract the real validation function without loading any scientific package.
	parsed <- parse(file.path(project, "f", "ems120.R"))
	assignments <- Filter(function(x) is.call(x) && identical(x[[1]], as.name("<-")) && is.symbol(x[[2]]), as.list(parsed))
	get_function <- function(name) {
		x <- Filter(function(x) identical(as.character(x[[2]]), name), assignments)
		stopifnot(length(x) == 1L)
		eval(x[[1]], envir = .GlobalEnv)
	}
	get_function("cache_complete")
	get_function("ensure_ml_python")
	assign("apply_phone_luck_use", function(x) { x$transformations <- x$transformations + 1L; x }, envir = .GlobalEnv)
	on.exit(rm(list = c("apply_phone_luck_use", "ensure_ml_python", "cache_complete"), envir = .GlobalEnv), add = TRUE)
	cols <- c("score", "dx.pipeline.version")
	for (cached in c(FALSE, TRUE)) {
		original <- readRDS(cache)
		if (cached) original <- lapply(original, apply_phone_luck_use)
		original <- lapply(original, apply_phone_luck_use)
		split <- ems120_report_cache(cache, cols, 2019, "fixture-v1", apply_cached_phone = cached)
		split <- lapply(split, apply_phone_luck_use)
		stopifnot(identical(original, split))
	}
	stopifnot(identical(before, tools::md5sum(cache)))
	cat("PASS: fresh and cached report paths apply exactly the original transforms; cache unchanged\n")

	for (args in list(list(cache, cols, 2019, "old-version"), list(cache, c(cols, "missing-column"), 2019, "fixture-v1"),
		list(file.path(root, "missing.rds"), cols, 2019, "fixture-v1"))) {
		stopifnot(inherits(try(do.call(ems120_report_cache, args), silent = TRUE), "try-error"))
	}
	stopifnot(identical(before, tools::md5sum(cache)))
	for (phase in c("cpu_prepare", "cpu_report")) {
		assign("ems120_phase", phase, envir = .GlobalEnv)
		failed <- try(ensure_ml_python(), silent = TRUE)
		stopifnot(inherits(failed, "try-error"), grepl("forbidden in CPU-only stages", as.character(failed), fixed = TRUE))
	}
	cat("PASS: missing/stale caches and ML calls on CPU fail without training\n")

	d[[1]]$score[[1]] <- 2
	saveRDS(d, cache)
	failed <- try(ems120_stage_restore("gpu_ml"), silent = TRUE)
	stopifnot(inherits(failed, "try-error"), grepl("changed after", as.character(failed), fixed = TRUE))
	Sys.setenv(EMS120_HANDOFF_DIR = file.path(root, "another-plan"))
	stopifnot(inherits(try(ems120_stage_restore("gpu_ml"), silent = TRUE), "try-error"))
	cat("PASS: changed caches and handoffs from another plan are rejected\n")
}
run_tests()
