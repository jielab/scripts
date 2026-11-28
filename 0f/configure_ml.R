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

# 🚩 Packages, paths, and shared inputs
sys <- Sys.info()[["sysname"]]
resolve_dir0 <- function() {
	normalize_dir0 <- function(path) {
		path <- gsub("\\\\", "/", path)
		if (sys == "Windows" && grepl("^[A-Za-z]:/?$", path)) return(sub("/$", "", path))
		normalizePath(path, winslash = "/", mustWork = FALSE)
	}
	env_dir0 <- Sys.getenv("EMS120_DIR0", unset = "")
	if (!nzchar(env_dir0)) env_dir0 <- Sys.getenv("DIR0", unset = "")
	if (nzchar(env_dir0)) return(normalize_dir0(env_dir0))
	candidates <- if (sys == "Windows") c(os_path("DIR_DAT")) else c(os_path("DIR_DAT"))
	hit <- candidates[file.exists(candidates)]
	if (length(hit)) normalize_dir0(hit[1]) else normalize_dir0(candidates[1])
}
dir0 <- resolve_dir0()
Sys.setenv(EMS120_DIR0 = dir0)

r_lib <- file.path(dir0, "R_lib", ifelse(sys == "Windows", "Windows", "Linux"))
dir.create(r_lib, recursive = TRUE, showWarnings = FALSE)
.libPaths(unique(c(r_lib, .libPaths())))

pick_python <- function(paths) {
	paths <- paths[nzchar(paths) & file.exists(paths)]
	if (length(paths)) normalizePath(paths[1], winslash = "/", mustWork = TRUE) else NA_character_
}

py_ai <- switch(
	sys,
	Windows = pick_python(c(
		Sys.getenv("RETICULATE_PYTHON", ""),
		file.path(Sys.getenv("USERPROFILE"), "anaconda3", "envs", "ai", "python.exe"),
		file.path(Sys.getenv("USERPROFILE"), "miniconda3", "envs", "ai", "python.exe"),
		file.path(os_path("DIR_WINDOWS_PROGRAMDATA"), "anaconda3/envs/ai/python.exe"),
		file.path(os_path("DIR_WINDOWS_PROGRAMDATA"), "miniconda3/envs/ai/python.exe"),
		file.path(os_path("DIR_DAT"), "anaconda3/envs/ai/python.exe"),
		file.path(os_path("DIR_DAT"), "miniconda3/envs/ai/python.exe")
	)),
	Linux = pick_python(c(
		Sys.getenv("RETICULATE_PYTHON", ""),
		file.path(Sys.getenv("HOME"), "anaconda3", "envs", "ai", "bin", "python"),
		file.path(os_path("DIR_HOME"), "anaconda3/envs/ai/bin/python"),
		file.path(os_path("DIR_HOME"), "anaconda3/envs/ai/bin/python")
	)),
	Darwin = pick_python(c(
		Sys.getenv("RETICULATE_PYTHON", ""),
		file.path(Sys.getenv("HOME"), "anaconda3", "envs", "ai", "bin", "python"),
		file.path(Sys.getenv("HOME"), "miniconda3", "envs", "ai", "bin", "python")
	)),
	NA_character_
)

if (is.na(py_ai)) {
	stop(
		"Cannot find conda env 'ai' Python. Please check conda path.\n",
		"Linux expected example: /home/huangj/anaconda3/envs/ai/bin/python\n",
		"Windows expected example: C:/Users/<USER>/anaconda3/envs/ai/python.exe",
		call. = FALSE
	)
}

Sys.setenv(RETICULATE_PYTHON = py_ai)
Sys.setenv(RETICULATE_AUTOCONFIGURE = "FALSE")

if (!requireNamespace("reticulate", quietly = TRUE)) {
	install.packages("reticulate", repos = "https://cloud.r-project.org")
}

library(reticulate)
reticulate::use_python(py_ai, required = TRUE)

required_py <- c("torch", "numpy", "pandas", "transformers", "sklearn", "openpyxl", "joblib")
missing_py <- required_py[!vapply(required_py, reticulate::py_module_available, logical(1))]
if (length(missing_py)) {
	stop(
		"Missing Python modules in conda env 'ai': ", paste(missing_py, collapse = ", "), "\n",
		"Run: conda activate ai && pip install ", paste(missing_py, collapse = " "),
		call. = FALSE
	)
}

torch <- reticulate::import("torch", delay_load = FALSE)
message(" Python: ", reticulate::py_config()$python)
message(" CUDA Is Available: ", torch$cuda$is_available())
try(torch$set_float32_matmul_precision("high"), silent = TRUE)
