# Load the shared path table relative to this script (also supports source()).
local({
    files <- Filter(Negate(is.null), lapply(sys.frames(), function(x) x$ofile))
    entry <- if (length(files)) tail(files, 1L)[[1L]] else sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE))[1L]
    if (is.na(entry)) stop("Cannot locate this script; use source() or Rscript")
    directory <- dirname(normalizePath(entry, winslash = "/"))
    repeat {
        adapter <- file.path(directory, "0f", "0os.R")
        if (file.exists(adapter)) break
        parent <- dirname(directory)
        if (parent == directory) stop("Cannot locate 0f/0os.R")
        directory <- parent
    }
    source(adapter, local = parent.env(environment()))
    os_paths_init(default_mode = "hpc")
})

# EMS120 runtime helpers. This file uses base R only when sourced, so path and
# dependency checks can run without loading the scientific analysis.
ems120_env_path <- function(name, fallback = "") {
	value <- Sys.getenv(name, unset = "")
	if (!nzchar(value)) value <- fallback
	if (!nzchar(value)) return("")
	normalizePath(path.expand(value), winslash = "/", mustWork = FALSE)
}

ems120_init_paths <- function(main_file) {
	os_paths_init()
	mode <- tolower(Sys.getenv("submit", unset = Sys.getenv("SUBMIT", unset = "hpc")))
	local <- identical(mode, "local")
	windows <- identical(.Platform$OS.type, "windows")
	dat_root <- ems120_env_path("DIR_DAT", if (local) if (windows) os_path("DIR_DAT") else os_path("DIR_DAT") else os_path("DIR_WORK"))
	work_root <- ems120_env_path("DIR_WORK", if (local) dat_root else os_path("DIR_WORK"))
	script_dir <- ems120_env_path("EMS120_SCRIPT_DIR", dirname(dirname(main_file)))
	shared_dir <- ems120_env_path("EMS120_SHARED_DIR", file.path(script_dir, "..", "0f"))
	analysis_root <- ems120_env_path("ANALYSIS_ROOT", file.path(dat_root, "analysis"))
	list(
		script = script_dir,
		shared = shared_dir,
		data = ems120_env_path("EMS120_DATA_DIR", file.path(dat_root, "ems120")),
		analysis = ems120_env_path("EMS120_ANALYSIS_DIR", file.path(analysis_root, "ems120")),
		files = ems120_env_path("EMS120_FILES_DIR", ems120_env_path("DIR_FILES", file.path(dat_root, "files")))
	)
}

ems120_require_packages <- function(packages, attach = FALSE) {
	missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
	if (length(missing)) {
		stop(
			"Missing R packages: ", paste(missing, collapse = ", "),
			". Install them in the selected R/Conda environment before submitting EMS120; the worker does not install packages.",
			call. = FALSE
		)
	}
	if (isTRUE(attach)) invisible(lapply(packages, function(pkg) suppressPackageStartupMessages(library(pkg, character.only = TRUE))))
	invisible(TRUE)
}

ems120_configure_python <- function() {
	ems120_require_packages("reticulate")
	py <- Sys.getenv("RETICULATE_PYTHON", unset = "")
	if (!nzchar(py) || !file.exists(py)) {
		stop("RETICULATE_PYTHON must identify the selected environment's Python. Run ems120.sh, or set this variable before starting R.", call. = FALSE)
	}
	py <- normalizePath(py, winslash = "/", mustWork = TRUE)
	Sys.setenv(RETICULATE_PYTHON = py, RETICULATE_AUTOCONFIGURE = "FALSE")
	reticulate::use_python(py, required = TRUE)
	modules <- c("torch", "numpy", "pandas", "transformers", "sklearn", "openpyxl", "joblib", "accelerate", "PIL")
	missing <- modules[!vapply(modules, reticulate::py_module_available, logical(1))]
	if (length(missing)) {
		stop("Missing Python modules in ", py, ": ", paste(missing, collapse = ", "), ". Install them before submitting the job.", call. = FALSE)
	}
	torch <- reticulate::import("torch", delay_load = FALSE)
	available <- isTRUE(torch$cuda$is_available())
	require_gpu <- !tolower(Sys.getenv("EMS120_REQUIRE_GPU", unset = "1")) %in% c("0", "false", "no", "off")
	if (require_gpu && !available) {
		stop("CUDA is unavailable inside R/reticulate. EMS120 requires its allocated GPU; CPU fallback is disabled.", call. = FALSE)
	}
	if (require_gpu && as.integer(torch$cuda$device_count()) != 1L) {
		stop("EMS120 is configured for one GPU but multiple devices are visible. Check the LSF GPU allocation/masking, or set CUDA_VISIBLE_DEVICES to one device for local execution. Refusing unallocated multi-GPU training.", call. = FALSE)
	}
	message("Python: ", reticulate::py_config()$python)
	message("CUDA available inside R/reticulate: ", available)
	try(torch$set_float32_matmul_precision("high"), silent = TRUE)
	invisible(TRUE)
}

ems120_save_rds_atomic <- function(object, file, ...) {
	dir.create(dirname(file), recursive = TRUE, showWarnings = FALSE)
	tmp <- tempfile(paste0(".", basename(file), "."), tmpdir = dirname(file))
	on.exit(if (file.exists(tmp)) unlink(tmp), add = TRUE)
	saveRDS(object, file = tmp, ...)
	# A same-directory rename publishes the complete RDS atomically on QiMing
	# Linux and WSL. Never remove an existing good cache before serialization.
	if (suppressWarnings(file.rename(tmp, file))) return(invisible(file))
	# Windows does not replace existing files with rename(). Preserve a rollback
	# copy there; all EMS120 writers are serialized by the shell's analysis lock.
	if (identical(.Platform$OS.type, "windows") && file.exists(file)) {
		backup <- tempfile(paste0(".", basename(file), ".previous."), tmpdir = dirname(file))
		if (!file.rename(file, backup)) stop("Cannot preserve existing RDS cache: ", file, call. = FALSE)
		if (file.rename(tmp, file)) {
			unlink(backup)
			return(invisible(file))
		}
		if (!file.rename(backup, file)) stop("RDS publication failed. The previous cache remains at: ", backup, call. = FALSE)
	}
	stop("Could not publish RDS cache; previous file was preserved: ", file, call. = FALSE)
}
