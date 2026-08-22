#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 🚩 Packages, paths, and shared inputs
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
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
	candidates <- if (sys == "Windows") c("D:") else c("/mnt/d", "/d")
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
		"C:/ProgramData/anaconda3/envs/ai/python.exe",
		"C:/ProgramData/miniconda3/envs/ai/python.exe",
		"D:/anaconda3/envs/ai/python.exe",
		"D:/miniconda3/envs/ai/python.exe"
	)),
	Linux = pick_python(c(
		Sys.getenv("RETICULATE_PYTHON", ""),
		file.path(Sys.getenv("HOME"), "anaconda3", "envs", "ai", "bin", "python"),
		"/home/huangj/anaconda3/envs/ai/bin/python",
		"/home/jiehuang001/anaconda3/envs/ai/bin/python"
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
message("✅ Python: ", reticulate::py_config()$python)
message("✅ CUDA Is Available: ", torch$cuda$is_available())
try(torch$set_float32_matmul_precision("high"), silent = TRUE)
