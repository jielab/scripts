source(file.path(Sys.getenv("FLT_CODE"), "runtime.R"))
library(data.table)
library(survival)

data_dir <- Sys.getenv("FLT_COHORT")
cfg <- readRDS(file.path(data_dir, "analysis_config.rds"))
base <- readRDS(file.path(data_dir, "dat.base.rds"))
imp <- readRDS(file.path(data_dir, "dat.base.imp.rds"))
legacy <- readRDS(file.path(data_dir, "recruitment-sensitivity.rds"))
meta <- as.data.table(readxl::read_excel(Sys.getenv("FLT_DICTIONARY"), sheet = "item182"))
raw <- readxl::read_excel(Sys.getenv("FLT_DICTIONARY"), sheet = "Sheet1")$name
# Only disjoint raw rating items enter the food-specific comparator. Restrict
# the focal comparisons to a fixed set of people with an available comparator.
chicken <- "p20638"
items <- intersect(grep("^p[0-9]+$", setdiff(raw, chicken), value = TRUE), names(imp))
stopifnot(length(items) > 10L, "depression" %in% cfg$Ys)
rating <- imp[items]
rating[] <- lapply(rating, function(x) replace(x, !(x %in% 1:9), NA_real_))
imp$general_liking <- rowMeans(rating, na.rm = TRUE)
imp$general_liking[rowMeans(!is.na(rating)) < .70] <- NA_real_

# Reproduce the older time origin on its own eligible cohort.
for (v in cfg$covs.all) {
    x <- legacy[[v]]
    value <- if (is.numeric(x)) mean(x, na.rm = TRUE) else names(which.max(table(x)))
    legacy[[v]][is.na(x)] <- value
}
base$age <- base$age_at_fpq
imp$age <- imp$age_at_fpq

fit <- function(data, y, food, mode, lag = 0, recruit = FALSE, comparator = FALSE,
                same_sample = FALSE, intake = FALSE, matched_intake = FALSE) {
    ending <- if (recruit) c(".recruit.eligible", ".recruit.time", ".recruit.event") else
        c(".eligible", ".t2e", ".Yt2e")
    key <- paste0(y, ending)
    use <- data[[key[1]]] & data[[key[2]]] > lag
    if (same_sample) use <- use & is.finite(data$general_liking)
    if (matched_intake) use <- use & is.finite(data$poultry_intake)
    covariates <- c(cfg$covs.all, if (comparator) "general_liking", if (intake) "poultry_intake")
    d <- data[which(use), covariates, drop = FALSE]
    d$time <- data[[key[2]]][which(use)] - lag
    d$event <- data[[key[3]]][which(use)]
    d$liking <- data[[food]][which(use)]
    d <- d[complete.cases(d), , drop = FALSE]
    row <- data.table(outcome = y, name = food, model = mode, N = nrow(d),
                      events = sum(d$event), person_years = sum(d$time),
                      beta = NA_real_, SE = NA_real_, HR = NA_real_,
                      lower = NA_real_, upper = NA_real_, P = NA_real_, PH_P = NA_real_,
                      status = "ok", detail = "")
    tryCatch({
        if (!nrow(d) || sum(d$event) < 5L || length(unique(d$liking)) < 2L)
            stop("Insufficient variation or events")
        model <- coxph(reformulate(c("liking", covariates), "Surv(time, event)"), d, singular.ok = FALSE)
        b <- coef(model)["liking"]
        se <- sqrt(vcov(model)["liking", "liking"])
        if (!is.finite(b) || !is.finite(se) || se <= 0) stop("Invalid model estimate")
        ph <- if (food == chicken && y == "depression")
            tryCatch(cox.zph(model)$table["liking", "p"], error = function(e) NA_real_) else NA_real_
        row[, `:=`(beta = b, SE = se, HR = exp(b), lower = exp(b - 1.96 * se),
                   upper = exp(b + 1.96 * se), P = 2 * pnorm(-abs(b / se)), PH_P = ph)]
    }, warning = function(w) row[, `:=`(status = "warning", detail = conditionMessage(w))],
       error = function(e) row[, `:=`(status = "unavailable", detail = conditionMessage(e))])
    row
}

panel <- c("acd", "stroke", "sleep_disorder", "depression", "anxiety")
all_pairs <- CJ(y = panel, food = meta$name)
cores <- if (.Platform$OS.type == "unix") flt_cores() else 1L
data.table::setDTthreads(1L)
comparison <- parallel::mclapply(seq_len(nrow(all_pairs)), function(i) {
    y <- all_pairs$y[i]; food <- all_pairs$food[i]
    rbind(fit(base, y, food, "Complete case"),
          fit(imp, y, food, "Two-year landmark", lag = 2),
          fit(legacy, y, food, "Recruitment time origin", recruit = TRUE))
}, mc.cores = cores, mc.preschedule = TRUE)
if (any(vapply(comparison, inherits, logical(1), "try-error"))) stop("Cox worker failed")
focal <- rbind(fit(imp, "depression", chicken, "Common sample", same_sample = TRUE),
               fit(imp, "depression", chicken, "Other-food liking adjusted",
                   comparator = TRUE, same_sample = TRUE))
intake_path <- file.path(Sys.getenv("FLT_YGB"), "Rdata", "diet.fudan.rds")
if (file.exists(intake_path)) {
    diet <- readRDS(intake_path)
    if (all(c("eid", "w.p26069", "w.p26121") %in% names(diet))) {
        stopifnot(!anyDuplicated(diet$eid))
        idx <- match(as.character(imp$eid), as.character(diet$eid))
        poultry <- diet$w.p26069[idx] + diet$w.p26121[idx]
        poultry[!is.finite(poultry) | poultry < 0] <- NA_real_
        imp$poultry_intake <- log1p(poultry)
        focal <- rbind(focal,
            fit(imp, "depression", chicken, "FPQ and WebQ common sample", same_sample = TRUE, matched_intake = TRUE),
            fit(imp, "depression", chicken, "Other liking and poultry intake adjusted",
                comparator = TRUE, same_sample = TRUE, intake = TRUE, matched_intake = TRUE))
    }
}
d <- rbindlist(c(comparison, list(focal)))
d[model %in% c("Complete case", "Two-year landmark", "Recruitment time origin"),
  FDR := p.adjust(P, "BH", n = .N), by = .(model, outcome)]
out <- file.path(flt_root, "01cox")
dir.create(out, recursive = TRUE, showWarnings = FALSE)
fwrite(d, file.path(out, "cox-sensitivity.csv"))
fwrite(data.table(item = items), file.path(out, "general-liking-items.csv"))
