# Shared adapters for data produced by FLT analysis stages.
# Add inside R as the HPC R launcher resets R_LIBS_USER.
flt_library <- c(Sys.getenv("FLT_R_LIBRARY"), Sys.getenv("DIET_R_LIBRARY"))
flt_library <- flt_library[nzchar(flt_library) & dir.exists(flt_library)]
if (length(flt_library)) .libPaths(c(flt_library, .libPaths()))
flt_phewas_cached <- function() {
    tryCatch({
        leads <- read.delim(file.path(Sys.getenv("FLT_FUMA_DIR"), "leadSNPs.txt"))$rsID
        if (!length(leads)) return(FALSE)
        roots <- unique(c(Sys.getenv("FLT_PHEWAS_CACHE"), Sys.getenv("FLT_PHEWAS_INPUT")))
        all(vapply(leads, function(snp) any(vapply(roots[nzchar(roots)], function(root) {
            tryCatch({
                file <- file.path(root, paste0("OpenGWAS_", snp, ".rds"))
                if (!file.exists(file)) return(FALSE)
                z <- readRDS(file)
                cols <- c("query_snp", "source", "study_id", "trait", "p", "beta", "n", "ea", "nea", "record_id")
                isTRUE(z$complete) && identical(z$p_limit, 1e-5) && is.data.frame(z$data) &&
                    all(cols %in% names(z$data)) && all(z$data$query_snp == snp) &&
                    all(is.finite(z$data$p) & z$data$p >= 0 & z$data$p < 1e-5) &&
                    !anyNA(z$data$trait) && all(nzchar(z$data$trait))
            }, error = function(e) FALSE)
        }, logical(1))), logical(1)))
    }, error = function(e) FALSE)
}
flt_root <- Sys.getenv("FLT_OUT")
flt_placo_file <- function(suffix) {
    root <- Sys.getenv("FLT_PLACO_DIR", file.path(Sys.getenv("FLT_OUT"), "03placo"))
    file.path(root, paste0(Sys.getenv("FLT_PLACO_PREFIX", "G732-mdd_2025_PLACO"), suffix))
}
flt_labels <- c(dementia_r13="All-cause dementia", anystroke_eur="Stroke",
                sleep_disorder_r13="Sleep disorders", mdd_2025="Depression", anxiety_2026="Anxiety disorders")
flt_label <- function(x) {
    x <- sub("^.*/", "", x)
    ifelse(x %in% names(flt_labels), unname(flt_labels[x]), x)
}
flt_tab <- function(name) {
    file <- file.path(Sys.getenv("FLT_YGB"), "rap", paste0(name, ".tab.gz"))
    if (file.exists(file)) return(file)
    path <- file.path(Sys.getenv("FLT_PREP"), paste0(name, ".tab.gz"))
    status <- system2(Sys.getenv("FLT_PYTHON"), shQuote(c(file.path(Sys.getenv("FLT_CODE"), "extract_tab.py"), name, path)))
    if (status != 0) stop("Field extraction failed for ", name)
    path
}
flt_panel <- function() {
    d <- read.delim(file.path(Sys.getenv("FLT_CONFIG"), "panels.tsv"), stringsAsFactors = FALSE)
    if (Sys.getenv("FLT_PANEL") == "extended") d else d[d$panel != "extended", ]
}
flt_cores <- function() max(1L, as.integer(Sys.getenv("WORKFLOW_CORES", "1")))
flt_shard <- function(n) {
    shard <- as.integer(Sys.getenv("FLT_SHARD", "1"))
    total <- as.integer(Sys.getenv("FLT_SHARDS", "1"))
    which((seq_len(n) - 1L) %% total + 1L == shard)
}
flt_run <- function(type) file.path(flt_root, paste0("02", type), Sys.getenv("FLT_PANEL", "primary"))
flt_workbooks <- function(tables, path) {
    if (!length(tables)) stop("No tables for ", path)
    dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
    staging <- paste0(path, ".tmp.xlsx")
    on.exit(unlink(staging))
    openxlsx::write.xlsx(tables, staging, overwrite = TRUE, asTable = FALSE)
    if (!file.rename(staging, path)) stop("Unable to publish ", path)
}
flt_png <- function(file, width, height, res = 300, pointsize = 12) {
    dir.create(dirname(file), recursive = TRUE, showWarnings = FALSE)
    grDevices::png(file, width = width, height = height, units = "in",
                   res = res, pointsize = pointsize, type = "cairo", bg = "white")
}
flt_cox <- function() {
    d <- read.csv(file.path(flt_root, "01cox/cox_results_table.csv"), check.names=FALSE)
    d <- d[d$Analysis == "continuous",]
    data.frame(name=d$name, label=d[["Food liking"]], category=d$Category, disease=d[["Disease label"]],
               N=d$N, events=d$Events, HR=d$HR, CI_low=d$Conf.low, CI_high=d$Conf.high, P=d[["P-value"]], FDR=d$FDR)
}
flt_uv <- function(which_direction, primary_only = TRUE) {
    path <- flt_run("uvmr")
    main <- read.csv(file.path(path, "uvmr_main_table.csv"), check.names=FALSE)
    main <- main[main$direction == which_direction, ]
    if (primary_only) {
        id <- if (which_direction == "forward") main$outcome else main$exposure
        main <- main[sub("^brain/", "", id) %in% names(flt_labels), ]
    }
    raw <- read.csv(file.path(path, "uvmr_results.csv"), check.names=FALSE)
    if (!nrow(raw)) {
        for (v in c("OR","b","OR_low","OR_high","beta_low","beta_high","pval")) raw[[v]] <- numeric(0)
    }
    at <- match(paste(main$pair_id, main$method_id), paste(raw$pair_id, raw$method_id))
    odds <- main$effect_type == "OR"
    point <- ifelse(odds, raw$OR[at], raw$b[at])
    low <- ifelse(odds, raw$OR_low[at], raw$beta_low[at])
    high <- ifelse(odds, raw$OR_high[at], raw$beta_high[at])
    p <- raw$pval[at]
    valid <- main$method_id == "mr_ivw" & is.finite(p) & main$status == "ok"
    family <- if (which_direction == "forward") main$outcome else main$exposure
    q <- rep(NA_real_, length(p))
    for (idx in split(which(valid), family[valid])) q[idx] <- p.adjust(p[idx], "BH")
    ci <- ifelse(is.finite(point), sprintf("%.12g (%.12g, %.12g)", point, low, high), NA_character_)
    data.frame(Exposure=flt_label(main$exposure), Outcome=flt_label(main$outcome), Method=main$method,
               effect_type=main$effect_type, estimate=point, lower=low, upper=high,
               "OR (95% CI)"=ifelse(odds, ci, NA_character_), "Beta (95% CI)"=ifelse(odds, NA_character_, ci),
               "P-value"=p, "P-FDR"=q, planned_FDR=main$FDR_method, status=main$status, check.names=FALSE)
}
flt_mv <- function() {
    path <- file.path(flt_run("mvmr"), "mvmr_all_results.csv")
    if (!file.exists(path)) return(data.frame(Exposure=character(), Outcome=character(),
        Method=character(), estimate=numeric(), lower=numeric(), upper=numeric(),
        "OR (95%CI)"=character(), "P-value"=numeric(), check.names=FALSE))
    d <- read.csv(path)
    d <- d[d$exposure == d$food, ]
    data.frame(Exposure=flt_label(d$food), Outcome=flt_label(d$outcome), Method=d$method, estimate=d$OR,
        lower=d$OR_low, upper=d$OR_high, "OR (95%CI)"=sprintf("%.12g (%.12g, %.12g)", d$OR, d$OR_low, d$OR_high),
        "P-value"=d$pval, check.names=FALSE)
}
