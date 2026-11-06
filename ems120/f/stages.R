# Scheduling boundaries only; statistical and model definitions stay in ems120.R/py.
ems120_phase <- Sys.getenv("EMS120_PHASE", unset = "all")
if (!ems120_phase %in% c("all", "cpu_prepare", "gpu_ml", "cpu_report")) {
	stop("Unknown EMS120_PHASE: ", ems120_phase, call. = FALSE)
}

ems120_stage_file <- function(stage) {
	root <- Sys.getenv("EMS120_HANDOFF_DIR", unset = "")
	if (!nzchar(root)) stop("Missing plan-specific EMS120_HANDOFF_DIR", call. = FALSE)
	file.path(root, paste0(stage, ".rds"))
}

ems120_stage_checksums <- function(files) {
	files <- sort(unique(normalizePath(files, winslash = "/", mustWork = TRUE)))
	if (any(file.info(files)$size == 0)) stop("Empty stage input", call. = FALSE)
	sums <- tools::md5sum(files)
	if (anyNA(sums)) stop("Cannot checksum stage inputs", call. = FALSE)
	sums
}

ems120_stage_publish <- function(stage, files, context = list()) {
	# Bind a handoff to this plan's exact cache contents and code. Preserve R's
	# random stream across processes without introducing a new scientific seed.
	code <- file.path(Sys.getenv("EMS120_SCRIPT_DIR"), "f", c("ems120.R", "ems120.py", "runtime.R", "stages.R"))
	state <- list(stage = stage, checksums = ems120_stage_checksums(c(files, code)),
		context = context, rng_kind = RNGkind(),
		rng_seed = if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) get(".Random.seed", envir = .GlobalEnv) else NULL)
	ems120_save_rds_atomic(state, ems120_stage_file(stage))
	cat("EMS120 stage complete:", stage, "\n")
	invisible(state)
}

ems120_stage_restore <- function(stage) {
	file <- ems120_stage_file(stage)
	if (!file.exists(file)) stop("Missing successful stage handoff: ", stage, call. = FALSE)
	state <- readRDS(file)
	if (!identical(state$stage, stage) || !length(state$checksums) ||
		!identical(state$checksums, ems120_stage_checksums(names(state$checksums)))) {
		stop("Stage inputs or code changed after ", stage, "; refusing mixed-run results.", call. = FALSE)
	}
	do.call(RNGkind, as.list(state$rng_kind))
	if (!is.null(state$rng_seed)) {
		assign(".Random.seed", state$rng_seed, envir = .GlobalEnv)
	} else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
		rm(".Random.seed", envir = .GlobalEnv)
	}
	state$context
}

ems120_report_cache <- function(file, required_cols, year_set, version, apply_cached_phone = TRUE) {
	if (!file.exists(file)) stop("CPU report requires dat/dat.list.rds; run the split workflow first.", call. = FALSE)
	dat <- readRDS(file)
	# Match the original cached branch exactly. A freshly trained cache has
	# already had this transform applied before serialization.
	if (apply_cached_phone) dat <- lapply(dat, apply_phone_luck_use)
	valid <- cache_complete(dat, required_cols, year_set) && all(vapply(dat, function(d) {
		is.null(d) || !nrow(d) || ("dx.pipeline.version" %in% names(d) && all(as.character(d$dx.pipeline.version) == version, na.rm = TRUE))
	}, logical(1)))
	if (!isTRUE(valid)) stop("CPU report cache is stale/incomplete; rerun the split workflow. CPU reporting never trains models.", call. = FALSE)
	dat
}
