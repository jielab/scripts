# Execution-only adapter: glmnet retains the solver, lambda path and CV aggregation.
# Fresh R workers receive matrices, not the analysis process's full cohort heap.
le8_glmnet_memory <- function() {
    status <- if (file.exists('/proc/self/status')) readLines('/proc/self/status') else character()
    kb <- function(name) {
        line <- status[startsWith(status, paste0(name, ':'))]
        if (length(line)) as.numeric(gsub('[^0-9]', '', line[1])) / 1024^2 else 0
    }
    limits <- as.numeric(c(Sys.getenv('LE8_PHASE_MEMORY_GIB', '32'),
        Sys.getenv('LE8_MEMORY_LIMIT_GB', '32'), Sys.getenv('LE8_RESOURCE_MEMORY_GIB', '32')))
    budget <- min(c(32, limits[is.finite(limits) & limits > 0]))
    parent <- kb('VmRSS') + kb('VmSwap')
    # Respect tighter cgroup limits and account for every process in that scope.
    if (file.exists('/proc/self/cgroup')) {
        line <- readLines('/proc/self/cgroup')
        line <- line[startsWith(line, '0::')]
        if (length(line)) {
            dir <- paste0('/sys/fs/cgroup', sub('^0::', '', line[1]))
            while (startsWith(dir, '/sys/fs/cgroup')) {
                read_number <- function(name) suppressWarnings(tryCatch(
                    as.numeric(readLines(file.path(dir, name), warn = FALSE)[1]) / 1024^3,
                    error = function(e) NA_real_))
                limit <- read_number('memory.max')
                if (is.finite(limit) && limit <= budget) {
                    budget <- limit
                    used <- read_number('memory.current')
                    if (is.finite(used)) parent <- max(parent, used)
                }
                if (dir == '/sys/fs/cgroup') break
                dir <- dirname(dir)
            }
        }
    }
    list(parent = parent, budget = budget)
}

le8_glmnet_worker_plan <- function(x, y, folds, parent_gib, budget_gib,
                                  cores = 8L, requested = cores) {
    if (any(!is.finite(c(parent_gib, budget_gib, cores, requested, folds))) ||
        parent_gib < 0 || budget_gib <= 0 || cores < 1 || requested < 1 || folds < 3)
        stop('Invalid glmnet resource budget')
    matrix_gib <- as.numeric(object.size(x)) / 1024^3
    # Exported matrix, fold subset, solver copy and serialization/GC headroom.
    # Include the dense Gram workspace used by Gaussian fits and R/package RSS.
    worker_gib <- .75 + 4 * matrix_gib + 8 * ncol(x)^2 / 1024^3 +
        4 * as.numeric(object.size(y)) / 1024^3
    reserve <- 2
    available <- budget_gib - parent_gib - reserve
    workers <- max(1L, min(8L, cores, requested, folds, floor(available / worker_gib)))
    if (available < max(.25, 3 * matrix_gib))
        stop('Insufficient RAM for glmnet; free memory or reduce the parent data footprint')
    list(workers = as.integer(workers), worker_gib = worker_gib, reserve_gib = reserve,
         parent_gib = parent_gib, budget_gib = budget_gib)
}

le8_glmnet_rng <- function() {
    if (exists('.Random.seed', .GlobalEnv, inherits = FALSE)) get('.Random.seed', .GlobalEnv) else NULL
}
le8_glmnet_restore_rng <- function(seed) {
    if (is.null(seed)) {
        if (exists('.Random.seed', .GlobalEnv, inherits = FALSE)) rm('.Random.seed', envir = .GlobalEnv)
    } else assign('.Random.seed', seed, .GlobalEnv)
}

# Independent R sessions use files/pipes rather than localhost TCP: delayed
# PSOCK replies can stall on WSL mirrored networking even for a two-second job.
# doParallel builds/evaluates exactly the same glmnet foreach expression.
le8_glmnet_fold_worker <- function(payload_path, job) {
    payload <- readRDS(payload_path)
    Sys.setenv(LE8_IN_WORKER = '1')
    if (requireNamespace('data.table', quietly = TRUE)) data.table::setDTthreads(1L)
    if (requireNamespace('RhpcBLASctl', quietly = TRUE)) {
        RhpcBLASctl::blas_set_num_threads(1L); RhpcBLASctl::omp_set_num_threads(1L)
    }
    do.call(glmnet::glmnet.control, payload$control)
    problem <- doParallel:::workerInit(payload$expr, payload$exports,
        payload$pkgname, payload$packages)
    if (!is.null(problem)) stop(problem)
    warnings <- character()
    result <- withCallingHandlers(doParallel:::evalWrapper(job), warning = function(w) {
        warnings <<- c(warnings, conditionMessage(w))
        invokeRestart('muffleWarning')
    })
    if (inherits(result, 'condition')) stop(result)
    list(value = result, warnings = unique(warnings))
}

le8_glmnet_jobs <- function(payload, jobs, workers, timeout) {
    scratch <- tempfile('le8-glmnet-', tmpdir = '/tmp')
    dir.create(scratch, mode = '0700')
    running <- list()
    on.exit({
        for (task in running) try(task$process$kill(), silent = TRUE)
        unlink(scratch, recursive = TRUE)
    }, add = TRUE)
    path <- file.path(scratch, 'inputs.rds')
    saveRDS(payload, path, compress = FALSE)
    results <- vector('list', length(jobs))
    next_job <- 1L
    completed <- 0L
    seen_warnings <- character()
    while (completed < length(jobs)) {
        while (next_job <= length(jobs) && length(running) < workers) {
            j <- next_job
            logfile <- file.path(scratch, paste0('fold-', j, '.log'))
            process <- callr::r_bg(le8_glmnet_fold_worker, list(path, jobs[[j]]),
                libpath = .libPaths(), stdout = logfile, stderr = logfile,
                user_profile = FALSE, system_profile = FALSE, supervise = TRUE,
                env = c(callr::rcmd_safe_env(), OMP_NUM_THREADS = '1',
                    OPENBLAS_NUM_THREADS = '1', MKL_NUM_THREADS = '1', BLIS_NUM_THREADS = '1'))
            running[[as.character(j)]] <- list(process = process, log = logfile,
                started = proc.time()[['elapsed']])
            next_job <- j + 1L
        }
        for (name in names(running)) {
            task <- running[[name]]
            elapsed <- proc.time()[['elapsed']] - task$started
            if (task$process$is_alive()) {
                if (elapsed >= timeout) stop('glmnet CV fold ', name, ' timed out after ',
                    timeout, ' seconds; no partial model was saved', call. = FALSE)
                next
            }
            result <- tryCatch(task$process$get_result(), error = function(e) {
                detail <- if (file.exists(task$log)) tail(readLines(task$log, warn = FALSE), 12) else character()
                stop('glmnet CV fold ', name, ' failed: ', conditionMessage(e),
                    if (length(detail)) paste0('\n', paste(detail, collapse = '\n')) else '', call. = FALSE)
            })
            for (w in setdiff(result$warnings, seen_warnings)) warning(w, call. = FALSE)
            seen_warnings <- union(seen_warnings, result$warnings)
            results[as.integer(name)] <- list(result$value)
            running[[name]] <- NULL
            completed <- completed + 1L
            message('[LE8] PROGRESS glmnet CV: fold ', name, ' finished (', completed, '/', length(jobs),
                '); elapsed=', round(elapsed, 1), ' s')
        }
        if (length(running)) Sys.sleep(.1)
    }
    results
}

le8_glmnet_backend <- function(obj, expr, envir, data) {
    message('[LE8] PROGRESS glmnet CV: full-data fit finished; dispatching folds to ', data$workers, ' workers')
    # Mirror doParallel's export and dots handling without exporting the live fit
    # frame. The expression and its variables still come directly from glmnet.
    exports <- environment(do.call(doParallel:::makeDotsEnv, eval(quote(list(...)), envir)))
    packages <- foreach::getexports(expr, exports, envir, bad = union(obj$noexport, obj$argnames))
    for (name in setdiff(obj$export, ls(exports))) assign(name, get(name, envir), exports)
    iterator <- iterators::iter(obj)
    accumulator <- foreach::makeAccum(iterator)
    payload <- list(expr = expr, exports = exports, pkgname = utils::packageName(envir),
        packages = unique(c(packages, obj$packages)), control = data$control)
    results <- le8_glmnet_jobs(payload, as.list(iterator), data$workers, data$timeout)
    accumulator(results, seq_along(results))
    problem <- foreach::getErrorValue(iterator)
    if (!is.null(problem)) stop(problem)
    message('[LE8] PROGRESS glmnet CV: fold fits finished; aggregating validation predictions')
    foreach::getResult(iterator)
}
le8_glmnet_backend_info <- function(data, item) switch(item,
    name = 'le8-callr', version = as.character(utils::packageVersion('callr')), workers = data$workers)

le8_cv_glmnet <- function(x, y, ..., parallel = NULL,
                          checkpoint_dir = getOption('le8.glmnet_checkpoint_dir', NULL)) {
    if (!requireNamespace('glmnet', quietly = TRUE)) stop('Prediction requires glmnet')
    args <- list(...)
    # This adapter owns execution only. All scientific arguments pass unchanged.
    args$parallel <- NULL
    folds <- if (is.null(args$foldid)) if (is.null(args$nfolds)) 10L else args$nfolds else max(args$foldid)
    before <- le8_glmnet_rng()
    control <- glmnet::glmnet.control()
    key <- path <- NULL
    if (!is.null(checkpoint_dir)) {
        # Match exact data, arguments, library versions and RNG for generated folds.
        key <- digest::digest(list(format = 1L, x = x, y = y, args = args,
            glmnet = as.character(utils::packageVersion('glmnet')),
            survival = as.character(utils::packageVersion('survival')),
            R = R.version.string, control = control,
            rng_kind = RNGkind(), rng = if (is.null(args$foldid)) before else NULL), algo = 'sha256')
        dir.create(checkpoint_dir, recursive = TRUE, showWarnings = FALSE, mode = '0700')
        path <- file.path(checkpoint_dir, paste0(key, '.rds'))
        replace <- isTRUE(get0('LE8_REPLACE', ifnotfound = FALSE))
        saved <- if (!replace && file.exists(path)) tryCatch(readRDS(path), error = function(e) NULL) else NULL
        if (is.list(saved) && identical(saved$key, key) && inherits(saved$fit, 'cv.glmnet')) {
            if (is.null(args$foldid)) le8_glmnet_restore_rng(saved$rng_after)
            for (w in saved$warnings) warning(w, call. = FALSE)
            message('[LE8] SKIP glmnet CV: validated model checkpoint ', basename(path))
            return(saved$fit)
        }
    }
    invisible(gc())
    memory <- le8_glmnet_memory()
    cores <- min(as.integer(Sys.getenv('N_CORES', '8')), as.integer(Sys.getenv('LE8_PHASE_CORES', '8')))
    choice <- Sys.getenv('LE8_GLMNET_WORKERS', 'auto')
    requested <- if (choice == 'auto') cores else as.integer(choice)
    # Avoid starting eight R sessions for a fit that completes in milliseconds.
    if (choice == 'auto' && length(x) < 50000L) requested <- 1L
    if (identical(Sys.getenv('LE8_IN_WORKER'), '1')) requested <- 1L
    plan <- le8_glmnet_worker_plan(x, y, folds, memory$parent, memory$budget, cores, requested)
    workers <- plan$workers
    if (workers > 1L && (!requireNamespace('doParallel', quietly = TRUE) || !requireNamespace('callr', quietly = TRUE)))
        stop('Parallel glmnet requires doParallel and callr in the selected R library')
    message('[LE8] START glmnet CV: N=', nrow(x), ' features=', ncol(x), ' folds=', folds,
        ' workers=', workers, ' parent_GiB=', round(plan$parent_gib, 2),
        ' worker_GiB=', round(plan$worker_gib, 2), ' reserve_GiB=', plan$reserve_gib,
        ' cap_GiB=', plan$budget_gib, '; full-data fit first, then independent folds')
    started <- proc.time()[['elapsed']]
    if (workers > 1L) {
        old_backend <- foreach:::getDoPar()
        globals <- get('.foreachGlobals', envir = asNamespace('foreach'))
        # Restore the actual info function, not a closure capturing this large fit.
        old_info <- if (exists('info', globals, inherits = FALSE)) get('info', globals) else NULL
        thread_names <- c('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'BLIS_NUM_THREADS')
        old_threads <- Sys.getenv(thread_names, unset = NA_character_)
        on.exit(for (nm in thread_names) {
            if (is.na(old_threads[[nm]])) Sys.unsetenv(nm)
            else do.call(Sys.setenv, setNames(list(old_threads[[nm]]), nm))
        }, add = TRUE)
        do.call(Sys.setenv, setNames(as.list(rep('1', length(thread_names))), thread_names))
        fold_timeout <- as.numeric(Sys.getenv('LE8_GLMNET_FOLD_TIMEOUT', '3600'))
        if (length(fold_timeout) != 1L || !is.finite(fold_timeout) || fold_timeout <= 0)
            stop('LE8_GLMNET_FOLD_TIMEOUT must be positive seconds')
        on.exit({
            if (is.null(old_info)) foreach::registerDoSEQ()
            else foreach::setDoPar(old_backend$fun, old_backend$data, old_info)
        }, add = TRUE)
        foreach::setDoPar(le8_glmnet_backend,
            list(workers = workers, control = control, timeout = fold_timeout), le8_glmnet_backend_info)
    }
    # Worker setup must not affect glmnet's original random fold generation.
    le8_glmnet_restore_rng(before)
    warnings <- character()
    fit <- tryCatch(withCallingHandlers(glmnet::cv.glmnet(x, y, ..., parallel = workers > 1L),
        warning = function(w) warnings <<- c(warnings, conditionMessage(w))), error = function(e) {
            message('[LE8] ERROR glmnet CV: ', conditionMessage(e))
            stop(e)
        })
    after <- le8_glmnet_rng()
    if (!is.null(path)) {
        tmp <- tempfile('.glmnet-', tmpdir = checkpoint_dir)
        on.exit(unlink(tmp), add = TRUE)
        saveRDS(list(key = key, fit = fit, rng_after = after, warnings = unique(warnings)), tmp, compress = FALSE)
        if (!file.rename(tmp, path)) stop('Cannot publish glmnet checkpoint: ', path)
        Sys.chmod(path, '0600')
    }
    message('[LE8] DONE glmnet CV: folds=', folds, ' workers=', workers,
        ' elapsed=', round((proc.time()[['elapsed']] - started) / 60, 2), ' min')
    fit
}

# Narrow compatibility receipt for reviewed execution-only changes.
# Require exact reviewed code versions AND unchanged inputs/settings. Saved
# scientific values and their original provenance are never rewritten.
le8_prediction_screen_cache <- function(path) {
    # A compatible execution change is a cache hit, not a recomputation warning.
    value <- suppressMessages(read_stage_cache(path))
    if (!is.null(value) || !cache_valid(path)) return(value)
    receipt <- file.path(dirname(Sys.getenv('LE8_FDIR')), 'config', 'c1_execution_compatibility.rds')
    if (!file.exists(receipt)) return(NULL)
    pair <- tryCatch(readRDS(receipt), error = function(e) NULL)
    if (!is.list(pair) || !identical(pair$after, .le8_loaded_code)) return(NULL)
    old <- tryCatch(readRDS(path), error = function(e) NULL)
    if (!is.list(old) || is.null(old$data) || !identical(old$analysis_options, le8_analysis_options())) return(NULL)
    fingerprint <- le8_stage_fingerprint
    scope <- new.env(parent = environment(fingerprint))
    environment(fingerprint) <- scope
    compatible <- FALSE
    for (code in c(list(pair$before), pair$also_before)) {
        scope$.le8_loaded_code <- code
        if (identical(old$source_signature, fingerprint())) { compatible <- TRUE; break }
    }
    if (!compatible) return(NULL)
    message('[LE8] SKIP prediction screen: unchanged inputs/settings; reviewed execution-only code change; ', basename(path))
    old$data
}
