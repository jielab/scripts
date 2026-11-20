# Shared adapters for data produced by FLT analysis stages.
flt_root <- Sys.getenv("FLT_OUT")
flt_labels <- c(dementia_r13="All-cause dementia", anystroke_eur="Stroke",
                sleep_disorder_r13="Sleep disorders", mdd_2025="Depression", anxiety_2026="Anxiety disorders")
flt_label <- function(x) {
    x <- sub("^.*/", "", x)
    ifelse(x %in% names(flt_labels), unname(flt_labels[x]), x)
}
flt_tab <- function(name) {
    file <- file.path(Sys.getenv("FLT_UKB"), "rap", paste0(name, ".tab.gz"))
    if (file.exists(file)) return(file)
    path <- file.path(Sys.getenv("FLT_PREP"), paste0(name, ".tab.gz"))
    status <- system2(Sys.getenv("FLT_PYTHON"), shQuote(c(file.path(Sys.getenv("FLT_CODE"), "extract_tab.py"), name, path)))
    if (status != 0) stop("Field extraction failed for ", name)
    path
}
flt_record_run <- function(type, path) {
    target <- file.path(flt_root, paste0("02", type), "current.txt")
    tmp <- tempfile(tmpdir=dirname(target))
    writeLines(normalizePath(path), tmp)
    if (!file.rename(tmp, target)) stop("Could not record completed run in ", target)
}
flt_run <- function(type) {
    override <- Sys.getenv(paste0("FLT_", toupper(type), "_RUN"))
    if (nzchar(override)) return(override)
    target <- file.path(flt_root, paste0("02", type), "current.txt")
    if (!file.exists(target)) stop("Run ./flt.sh ", type, " first, or set FLT_", toupper(type), "_RUN.")
    readLines(target, n=1)
}
flt_cox <- function() {
    d <- read.csv(file.path(flt_root, "01cox/cox_results_table.csv"), check.names=FALSE)
    d <- d[d$Analysis == "continuous",]
    data.frame(name=d$name, label=d[["Food liking"]], category=d$Category, disease=d[["Disease label"]],
               N=d$N, events=d$Events, HR=d$HR, CI_low=d$Conf.low, CI_high=d$Conf.high, P=d[["P-value"]], FDR=d$FDR)
}
flt_uv <- function(which_direction) {
    path <- flt_run("uvmr")
    main <- read.csv(file.path(path, "uvmr_main_table.csv"), check.names=FALSE)
    main <- main[main$direction == which_direction, ]
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
    d <- read.csv(file.path(flt_run("mvmr"), "mvmr_all_results.csv"))
    d <- d[d$exposure == d$food, ]
    data.frame(Exposure=flt_label(d$food), Outcome=flt_label(d$outcome), Method=d$method, estimate=d$OR,
        lower=d$OR_low, upper=d$OR_high, "OR (95%CI)"=sprintf("%.12g (%.12g, %.12g)", d$OR, d$OR_low, d$OR_high),
        "P-value"=d$pval, check.names=FALSE)
}
