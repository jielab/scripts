# R session setup, atomic caches, stage handoffs and explicit processor checks.
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
	analysis_root <- if (local) ems120_env_path("ANALYSIS_ROOT", file.path(dat_root, "analysis")) else ems120_env_path("PERSISTENT_ANALYSIS_ROOT", file.path(work_root, "analysis"))
	default_data <- file.path(dat_root, "phe", "ems120")
	list(
		script = script_dir,
		shared = shared_dir,
		data = ems120_env_path("EMS120_DATA_DIR", default_data),
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
	require_gpu <- identical(Sys.getenv("EMS120_DEVICE", unset = "cpu"), "gpu")
	if (require_gpu && !available) {
		stop("CUDA is unavailable inside R/reticulate. This EMS120 stage explicitly requested GPU; CUDA is required.", call. = FALSE)
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


# Scheduling boundaries only; statistical and model definitions stay in ems120.R/py.
ems120_phase <- Sys.getenv("EMS120_PHASE", unset = "prepare")
ems120_replace <- identical(toupper(Sys.getenv("EMS120_REPLACE", unset = "TRUE")), "TRUE")
if (!ems120_phase %in% c("prepare", "phone", "geo", "dx", "report")) {
	stop("Unknown EMS120_PHASE: ", ems120_phase, call. = FALSE)
}

ems120_stage_file <- function(stage) {
	root <- Sys.getenv("EMS120_HANDOFF_DIR", unset = "")
	if (!nzchar(root)) stop("Missing plan-specific EMS120_HANDOFF_DIR", call. = FALSE)
	file.path(root, paste0(stage, ".rds"))
}

ems120_stage_inputs <- function(files) {
    files <- sort(unique(normalizePath(files, winslash = "/", mustWork = TRUE)))
    info <- file.info(files)
    if (anyNA(info$size) || any(info$size <= 0)) stop("Missing/empty stage input", call. = FALSE)
    data.frame(path = files, size = info$size, modified = as.numeric(info$mtime), stringsAsFactors = FALSE)
}

ems120_stage_publish <- function(stage, files, context = list()) {
    # Bind each handoff to the actual files read by its successor. The R RNG
    # stream and report objects are carried across the process boundary.
    state <- list(stage = stage, inputs = ems120_stage_inputs(files), context = context,
        rng_kind = RNGkind(), rng_seed = if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) get(".Random.seed", envir = .GlobalEnv) else NULL)
    ems120_save_rds_atomic(state, ems120_stage_file(stage))
    cat("EMS120 stage complete:", stage, "\n")
    invisible(state)
}

ems120_stage_restore <- function(stage) {
    file <- ems120_stage_file(stage)
    if (!file.exists(file)) stop("Missing successful stage handoff: ", stage, call. = FALSE)
    state <- readRDS(file)
    if (!identical(state$stage, stage) || !is.data.frame(state$inputs) || !nrow(state$inputs) ||
        !identical(state$inputs, ems120_stage_inputs(state$inputs$path))) {
        stop("Stage inputs changed after ", stage, "; refusing mixed-run results.", call. = FALSE)
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
	if (!file.exists(file)) stop("CPU report requires dat/dat.list.rds; run the full workflow first.", call. = FALSE)
	if (!ems120_cache_call("check",file)) stop("CPU report cache provenance is stale/unverified; rerun upstream stages. Reporting never trains models.")
	dat <- ems120_read_rds(file)
	# Match the original cached branch exactly. A freshly trained cache has
	# already had this transform applied before serialization.
	if (apply_cached_phone) dat <- lapply(dat, apply_phone_luck_use)
	valid <- cache_complete(dat, required_cols, year_set) && all(vapply(dat, function(d) {
		is.null(d) || !nrow(d) || ("dx.pipeline.version" %in% names(d) && all(!is.na(d$dx.pipeline.version) & as.character(d$dx.pipeline.version) == version))
	}, logical(1)))
	if (!isTRUE(valid)) stop("CPU report cache is stale/incomplete; rerun the full workflow. CPU reporting never trains models.", call. = FALSE)
	dat
}


ems120_check_session <- function() {
	# Worker-only environment check. The planner selects the processor and checks run inside the allocated stage.
	script_dir <- Sys.getenv("EMS120_SCRIPT_DIR", unset = "")
	if (!nzchar(script_dir)) {
		arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
		if (!length(arg)) stop("Set EMS120_SCRIPT_DIR or run this file with Rscript.", call. = FALSE)
		script_dir <- dirname(dirname(normalizePath(sub("^--file=", "", arg[[1]]), winslash = "/", mustWork = TRUE)))
	}
	ems120_paths <- ems120_init_paths(file.path(script_dir, "f", "ems120.R"))
	ems120_require_packages(c("readxl", "writexl", "data.table", "tidyverse", "scales", "RColorBrewer", "reticulate", "lubridate", "patchwork", "zoo", "broom", "forcats", "circlize", "jpeg", "png", "grid", "lattice", "cowplot", "jsonlite", "base64enc"))
	ems120_configure_python()

	reticulate::py_run_string(paste(c(
		"import sys, torch, PIL",
		"print('EMS120 runtime Python:', sys.executable, flush=True)",
		"print('PyTorch:', torch.__version__, 'CUDA runtime:', torch.version.cuda, 'Pillow:', PIL.__version__, flush=True)",
		"if __import__('os').environ.get('EMS120_DEVICE', 'cpu') == 'gpu':",
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

    # An explicitly requested phone GPU must execute a CUDA XGBoost probe even
    # when the selected cached winner is a CPU comparator.
    if (identical(ems120_phase, "phone") && identical(Sys.getenv("EMS120_DEVICE"), "gpu")) {
        phone_probe_env <- new.env(parent = globalenv())
        reticulate::source_python(file.path(script_dir, "f", "ems120.py"), envir = phone_probe_env)
        stopifnot(isTRUE(phone_probe_env$check_phone_gpu()))
        cat("Phone XGBoost CUDA training probe: PASS\n")
    }

	cat("EMS120 runtime check completed; no patient data or model files were opened.\n")

	invisible(TRUE)
}

# Sourcing defines helpers only. The worker explicitly requests the runtime check.
if (sys.nframe() == 0L) {
    args <- commandArgs(trailingOnly = TRUE)
    if (!identical(args, "check")) stop("Usage: Rscript session.R check", call. = FALSE)
    ems120_check_session()
}


# Cross-run cache receipts; the existing within-plan handoff checks remain intact.
ems120_cache_spec <- function(stage, inputs, parameters=list(), dependencies=character()) {
    list(stage=stage,inputs=as.list(inputs),parameters=parameters,dependencies=as.list(dependencies))
}
.ems120_rds_cache <- new.env(parent=emptyenv())
ems120_read_rds <- function(file) {
    key <- normalizePath(file,winslash="/",mustWork=FALSE)
    if (exists(key,envir=.ems120_rds_cache,inherits=FALSE)) {
        value <- get(key,envir=.ems120_rds_cache,inherits=FALSE)
        rm(list=key,envir=.ems120_rds_cache)
        return(value)
    }
    readRDS(file)
}

ems120_cached_rds_valid <- function(file,spec) {
    # Read a large aggregate once, then transfer the validated object to its
    # caller instead of decoding the same multi-year RDS a second time.
    tryCatch({
        value <- readRDS(file)
        if (!is.null(spec)) {
            stage <- sub("/.*$","",spec$stage)
            required <- switch(stage,raw=c("ID","phone",vars),phone=c("ID","phone",phone.required),
                               geo=c("ID","phone",phone.required,"geo.type1"),
                               dx=final.required,dx_predictions=final.required,character())
            valid_year <- function(d) is.null(d) || (is.data.frame(d) && all(required %in% names(d)))
            valid <- if (!is.null(spec$parameters$years)) {
                is.list(value) && all(as.character(spec$parameters$years) %in% names(value)) &&
                    all(vapply(value[as.character(spec$parameters$years)],valid_year,logical(1)))
            } else valid_year(value)
            if (!isTRUE(valid)) return(FALSE)
        }
        assign(normalizePath(file,winslash="/",mustWork=FALSE),value,envir=.ems120_rds_cache)
        TRUE
    },error=function(e) FALSE)
}

ems120_cache_call <- function(action, output, spec=NULL) {
    if (action == "check") {
        if (!is.null(spec) && ems120_replace) return(FALSE)
        if (!file.exists(output) || is.na(file.info(output)$size) || file.info(output)$size <= 0) return(FALSE)
    }
    receipt <- paste0(output,".dependencies.json")
    args <- c(shQuote(file.path(script.dir,"f","cache_manifest.py")),action,"--receipt",shQuote(receipt))
    if (action == "check" && identical(ems120_phase, "report")) args <- c(args,"--outputs-only")
    if (!is.null(spec)) {
        tmp <- tempfile(fileext=".json"); on.exit(unlink(tmp),add=TRUE)
        jsonlite::write_json(spec,tmp,auto_unbox=TRUE,null="null",digits=NA)
        args <- c(args,"--spec",shQuote(tmp))
    }
    if (action=="seal") {
        outputs <- output
        if (!is.null(spec$parameters$years) && spec$stage %in% c("raw","phone","geo","dx_predictions")) {
            year_stage <- if (identical(spec$stage,"dx_predictions")) "dx" else spec$stage
            outputs <- c(outputs,vapply(spec$parameters$years,function(y) ems120_year_file(year_stage,y),character(1)))
        }
        args <- c(args,as.vector(rbind("--output",shQuote(outputs))))
    }
    py <- Sys.getenv("RETICULATE_PYTHON",unset="python3")
    status <- system2(py,args,stdout=FALSE,stderr="")
    if (action=="seal" && status!=0L) stop("Cannot publish EMS120 cache receipt: ",output)
    if (action=="check" && identical(status,0L) && grepl("\\.rds$",output)) return(ems120_cached_rds_valid(output,spec))
    identical(status,0L)
}

ems120_year_file <- function(stage,year) file.path(dir.cache,"years",stage,paste0(year,".rds"))
ems120_year_spec <- function(spec,stage,year,source) {
    spec$stage <- paste(stage,year,sep="/")
    spec$parameters$years <- NULL
    spec$parameters$year <- year
    if (stage=="raw") spec$parameters$phone <- vars.phone[as.character(year)]
    # Aggregate inputs are replaced by the one upstream year; model receipts
    # remain dependencies, so label changes invalidate all dependent years.
    keep <- spec$inputs[vapply(spec$inputs,function(x) grepl("dependencies.json$",x),logical(1))]
    spec$inputs <- c(list(source),keep)
    deps <- unlist(spec$dependencies,use.names=FALSE)
    deps <- deps[!grepl("\\.rds.dependencies.json$",deps)]
    if (stage!="raw") deps <- c(deps,paste0(source,".dependencies.json"))
    if (stage=="geo") spec$inputs <- c(spec$inputs,list(file.path(dir.dat,"szu_xia",paste0(year,".geo.xlsx"))))
    spec$dependencies <- as.list(deps)
    spec
}


# Scope-level report completion. The R script calls finish on successful end-step
# exits as well as the full run; a failed report never leaves a complete receipt.
ems120_report_snapshot <- function() {
    paths <- list.files(dir.analysis, recursive=TRUE, full.names=TRUE, all.files=FALSE)
    relative <- substring(paths, nchar(dir.analysis) + 2L)
    keep <- !grepl("^(ml|logs|tmp)/|^raw/hfl_cv5/|\\.rds$|\\.dependencies\\.json$|^dat/2019\\.train_dx\\.xlsx$", relative)
    paths <- paths[keep & !dir.exists(paths)]
    info <- file.info(paths)
    data.frame(path=paths,size=info$size,modified=as.numeric(info$mtime),stringsAsFactors=FALSE)
}

ems120_report_begin <- function() {
    stopifnot(identical(ems120_phase, "report"))
    scope <- paste(ems120_start_step_requested, if (nzchar(ems120_end_step)) ems120_end_step else "end", sep="-")
    .ems120_report_state <<- list(
        receipt=file.path(dir.analysis,paste0(".report-",scope,".complete.json")),
        spec=tempfile("ems120-report-spec-",fileext=".json"),
        before=ems120_report_snapshot())
    spec <- ems120_cache_spec(paste0("report/",scope), c(file.dat.rds,file.dat.phone.rds,list.files(dir_dx_cv,full.names=TRUE,pattern="\\.csv$")),
        parameters=list(start=ems120_start_step_requested,end=ems120_end_step,
                        sensitivity_permutations=Sys.getenv("EMS120_SENS_N_PERM",unset="")))
    jsonlite::write_json(spec,.ems120_report_state$spec,auto_unbox=TRUE,null="null",digits=NA)
    args <- c(shQuote(file.path(script.dir,"f","cache_manifest.py")),"check",
              "--receipt",shQuote(.ems120_report_state$receipt),"--spec",shQuote(.ems120_report_state$spec))
    py <- Sys.getenv("RETICULATE_PYTHON",unset="python3")
    reusable <- !ems120_replace && identical(system2(py,args,stdout=FALSE,stderr=""),0L)
    if (reusable) {
        cat("Complete report scope reused:",scope,"\n")
        unlink(.ems120_report_state$spec)
        return(TRUE)
    }
    unlink(.ems120_report_state$receipt)
    FALSE
}

ems120_report_finish <- function() {
    if (!exists(".ems120_report_state", envir=.GlobalEnv, inherits=FALSE)) return(invisible(FALSE))
    state <- .ems120_report_state
    after <- ems120_report_snapshot()
    idx <- match(after$path,state$before$path)
    changed <- is.na(idx) | after$size != state$before$size[idx] | after$modified != state$before$modified[idx]
    outputs <- after$path[which(changed & !is.na(after$size) & after$size > 0)]
    if (!length(outputs)) stop("Report completed without any recorded output files")
    args <- c(shQuote(file.path(script.dir,"f","cache_manifest.py")),"seal",
              "--receipt",shQuote(state$receipt),"--spec",shQuote(state$spec),
              as.vector(rbind("--output",shQuote(outputs))))
    py <- Sys.getenv("RETICULATE_PYTHON",unset="python3")
    status <- system2(py,args,stdout=FALSE,stderr="")
    if (!identical(status,0L)) stop("Cannot publish complete report receipt")
    unlink(state$spec)
    invisible(TRUE)
}
