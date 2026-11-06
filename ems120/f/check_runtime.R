# Worker-only environment check. Launch through ems120.sh --check-gpu so HPC
# checks execute inside a GPU allocation rather than on the login node.
script_dir <- Sys.getenv("EMS120_SCRIPT_DIR", unset = "")
if (!nzchar(script_dir)) {
	arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
	if (!length(arg)) stop("Set EMS120_SCRIPT_DIR or run this file with Rscript.", call. = FALSE)
	script_dir <- dirname(dirname(normalizePath(sub("^--file=", "", arg[[1]]), winslash = "/", mustWork = TRUE)))
}
source(file.path(script_dir, "f", "runtime.R"))
ems120_paths <- ems120_init_paths(file.path(script_dir, "f", "ems120.R"))
ems120_require_packages(c("readxl", "writexl", "data.table", "tidyverse", "scales", "RColorBrewer", "reticulate", "lubridate", "patchwork", "zoo", "broom", "forcats", "circlize", "jpeg", "png", "grid", "lattice", "cowplot", "jsonlite", "base64enc"))
ems120_configure_python()

reticulate::py_run_string(paste(c(
	"import sys, torch, PIL",
	"print('EMS120 runtime Python:', sys.executable, flush=True)",
	"print('PyTorch:', torch.__version__, 'CUDA runtime:', torch.version.cuda, 'Pillow:', PIL.__version__, flush=True)",
	"if torch.cuda.is_available():",
	"    print('Visible CUDA devices:', torch.cuda.device_count(), flush=True)",
	"    print('Allocated logical device 0:', torch.cuda.get_device_name(0), flush=True)",
	"    _ems120_probe_x = torch.ones((16, 16), device='cuda:0')",
	"    _ems120_probe_y = _ems120_probe_x @ _ems120_probe_x",
	"    torch.cuda.synchronize()",
	"    assert bool(torch.all(_ems120_probe_y == 16).item()), 'CUDA arithmetic check failed'",
	"    del _ems120_probe_x, _ems120_probe_y",
	"    torch.cuda.empty_cache()",
	"    print('R/reticulate CUDA allocation and arithmetic check: PASS', flush=True)",
	"else:",
	"    print('R/reticulate CPU mode explicitly enabled; CUDA check skipped.', flush=True)"
), collapse = "\n"))

# Explicit --check-gpu verifies both neural-network CUDA and the phone
# regression candidate. Ordinary runs defer the XGBoost probe until training
# is needed, so a valid cached CPU-selected phone estimator remains reusable.
check_gpu_only <- tolower(Sys.getenv("EMS120_CHECK_GPU_ONLY", unset = "0")) %in% c("1", "true", "yes", "on")
phone_gpu_policy <- tolower(Sys.getenv("EMS120_PHONE_GPU", unset = "required"))
if (check_gpu_only && !identical(phone_gpu_policy, "off")) {
	phone_probe_env <- new.env(parent = globalenv())
	reticulate::source_python(file.path(script_dir, "f", "ems120.py"), envir = phone_probe_env)
	phone_cuda <- isTRUE(phone_probe_env$check_phone_gpu())
	cat(if (phone_cuda) "Phone XGBoost CUDA training probe: PASS\n" else "Phone XGBoost CUDA unavailable; optional auto policy permits CPU comparison.\n")
}
cat("EMS120 runtime check completed; no patient data or model files were opened.\n")
