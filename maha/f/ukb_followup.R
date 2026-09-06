# MAHA-specific timing contract. No outcome-dependent choice of landmark.
maha_followup_version <- "2026-09-06-webq-window-v1"
maha_date <- function(x) {
  if (inherits(x, "Date")) return(x)
  if (is.numeric(x)) return(as.Date(x, origin = "1970-01-01"))
  as.Date(as.character(x))
}

maha_webq_window <- function(raw) {
  date_cols <- paste0("p105010_i", 0:4)
  stopifnot(all(c("eid", date_cols) %in% names(raw)), !anyDuplicated(raw$eid))
  dates <- lapply(raw[date_cols], maha_date)
  last <- do.call(pmax, c(dates, list(na.rm = TRUE)))
  last[!is.finite(as.numeric(last))] <- as.Date(NA)
  # A recorded duration/energy without a date makes its temporal ordering unknown.
  undated <- rep(FALSE, nrow(raw))
  for (k in 0:4) {
    present <- rep(FALSE, nrow(raw))
    for (v in paste0(c("p20082_i", "p26002_i"), k)) {
      if (v %in% names(raw)) present <- present | !is.na(raw[[v]])
    }
    undated <- undated | (present & is.na(dates[[k + 1L]]))
  }
  data.frame(eid = as.character(raw$eid), webq_window_end = last,
             webq_undated = undated, webq_n_dates = rowSums(!is.na(as.data.frame(dates))))
}

maha_set_landmark <- function(dat, window, score_cols) {
  i <- match(as.character(dat$eid), window$eid)
  dat$webq_window_end <- window$webq_window_end[i]
  dat$webq_undated <- window$webq_undated[i]
  scored <- rowSums(!is.na(dat[score_cols])) > 0
  invalid <- scored & (is.na(dat$webq_window_end) | is.na(dat$webq_undated) | dat$webq_undated)
  dat$webq_timing_excluded <- invalid
  if (any(invalid)) message("[LANDMARK] Excluding ", sum(invalid), " scored participants with missing/undated WebQ dates from survival analyses.")
  attend <- maha_date(dat$date_attend)
  dat$diet_landmark <- pmax(attend, dat$webq_window_end)
  dat$diet_landmark[invalid] <- as.Date(NA)
  if ("start_date" %in% names(dat)) {
    inconsistent <- !is.na(dat$start_date) & !is.na(dat$diet_landmark) & maha_date(dat$start_date) > dat$diet_landmark
    if (any(inconsistent)) stop("Stored score start_date is later than the raw WebQ window.")
  }
  dat$age.landmark <- as.numeric(dat$age) + as.numeric(dat$diet_landmark - attend) / 365.25
  dat
}

maha_outcome_followup <- function(dat, outcome, end_date) {
  required <- c("diet_landmark", "date_lost", "date_death")
  event_col <- if (outcome == "death") "date_death" else paste0("fod_icd10_", outcome)
  if (!all(c(required, event_col) %in% names(dat))) stop("Missing follow-up input for ", outcome)
  start <- maha_date(dat$diet_landmark)
  event_date <- maha_date(dat[[event_col]])
  dead <- maha_date(dat$date_death)
  censor <- pmin(maha_date(dat$date_lost), dead, as.Date(end_date), na.rm = TRUE)
  prevalent <- !is.na(start) & !is.na(event_date) & event_date <= start
  unknown <- rep(FALSE, nrow(dat))
  # Use only evidence for THIS disease, never a general CVD flag for every outcome.
  if (outcome != "death") for (v in c(paste0("cnt_icd10_", outcome), paste0("icd10Ct_", outcome))) {
    if (v %in% names(dat)) unknown <- unknown | (is.na(event_date) & !is.na(dat[[v]]) & dat[[v]] > 0)
  }
  eligible <- !is.na(start) & !is.na(censor) & censor > start & !prevalent & !unknown
  incident <- eligible & !is.na(event_date) & event_date > start & event_date <= censor
  stop_date <- censor
  stop_date[incident] <- event_date[incident]
  dat[[paste0(outcome, ".t2e")]] <- ifelse(eligible, as.numeric(stop_date - start) / 365.25, NA_real_)
  dat[[paste0(outcome, ".Yt2e")]] <- ifelse(eligible, as.integer(incident), NA_integer_)
  attr(dat, "maha_outcome_audit") <- data.frame(outcome = outcome, N_input = nrow(dat),
    N_eligible = sum(eligible), events = sum(incident), prevalent_through_landmark = sum(prevalent),
    unknown_event_date = sum(unknown), no_positive_followup = sum(is.na(start) | is.na(censor) | censor <= start),
    timing_excluded = sum(dat$webq_timing_excluded %in% TRUE),
    time_origin = "max(enrolment date, last of all recorded WebQ dates)", admin_end = as.character(end_date))
  dat
}

# Raw P values are immutable; numerical floor and common display cap are separate.
maha_phewas_p <- function(res) {
  stopifnot(all(c("p", "snp") %in% names(res)))
  if (any(!is.na(res$p) & (res$p < 0 | res$p > 1))) stop("Invalid PheWAS P value")
  res$p_raw <- res$p
  res$p_underflow <- !is.na(res$p) & res$p == 0
  res$logp <- -log10(pmax(res$p, .Machine$double.xmin))
  res$logp_display <- pmin(res$logp, 50)
  res$display_capped <- !is.na(res$logp) & res$logp > 50
  res$FDR <- ave(res$p, res$snp, FUN = function(p) p.adjust(p, "BH"))
  res
}

maha_raw_phewas <- function(dat, Xs, covs, phenotype_file, cache_dir, replace = FALSE) {
  # This package's addPhecodeInfo expects its pheinfo dataset on the search path.
  suppressPackageStartupMessages(library(PheWAS))
  dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
  phenos <- readRDS(phenotype_file)
  if (!"id" %in% names(phenos)) stop("PheWAS phenotype input has no id column")
  phenos$id <- as.character(phenos$id)
  dat$eid <- as.character(dat$eid)
  # Full-content fingerprint prevents stale resumption after input or model changes.
  tf <- tempfile(fileext = ".rds"); on.exit(unlink(tf), add = TRUE)
  saveRDS(list(version = "raw-p-v1", dat = dat[, c("eid", Xs, covs)],
               phenotypes = phenos, package = as.character(packageVersion("PheWAS"))), tf, version = 2)
  signature <- unname(tools::md5sum(tf))
  cores <- suppressWarnings(as.integer(Sys.getenv("PHEWAS_CORES", "1")))
  if (is.na(cores) || cores < 1) cores <- 1L
  out <- lapply(Xs, function(x) {
    path <- file.path(cache_dir, paste0(x, ".raw.rds"))
    cached <- if (!replace && file.exists(path)) readRDS(path) else NULL
    if (is.list(cached) && identical(cached$signature, signature)) return(cached$result)
    message("[RAW PHEWAS] Fitting ", x)
    g <- dat[, c("eid", x)]; names(g)[1] <- "id"
    cv <- dat[, c("eid", covs)]; names(cv)[1] <- "id"
    r <- PheWAS::phewas(phenotypes = phenos, genotypes = g, covariates = cv,
                       cores = cores, significance.threshold = "bonferroni")
    r <- as.data.frame(PheWAS::addPhecodeInfo(r))
    names(r)[names(r) == "group"] <- "category"
    r <- maha_phewas_p(r)
    saveRDS(list(signature = signature, result = r), path)
    r
  })
  list(res = dplyr::bind_rows(out), signature = signature)
}
