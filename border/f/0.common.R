options(stringsAsFactors = FALSE, width = 180)
suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(tibble); library(survival)
})

dir0 <- ifelse(Sys.info()[["sysname"]] == "Windows", "D:", "/mnt/d")
indir <- Sys.getenv("BORDER_INDIR", file.path(dir0, "data", "ukb", "phe"))
outdir <- Sys.getenv("BORDER_OUTDIR", file.path(dir0, "analysis", "border"))
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
rawdir <- file.path(outdir, "raw"); dir.create(rawdir, recursive = TRUE, showWarnings = FALSE)
figdir <- file.path(outdir, "fig"); dir.create(figdir, recursive = TRUE, showWarnings = FALSE)

csv_env <- function(name, default = "") {
  z <- trimws(strsplit(Sys.getenv(name, default), ",", fixed = TRUE)[[1]])
  unique(z[nzchar(z)])
}
EXPOSURES <- csv_env("BORDER_EXPOSURES", "border12,border1n,border_rank,border_relative")
OUTCOMES  <- csv_env("BORDER_OUTCOMES", "cvd_cad,ra")
REPLACE   <- toupper(Sys.getenv("BORDER_REPLACE", "FALSE")) %in% c("TRUE","1","YES")

ALL_EXPOSURES <- c("border12", "border1n", "border_rank", "border_relative")
if (!length(EXPOSURES) || length(setdiff(EXPOSURES, ALL_EXPOSURES))) stop("Invalid BORDER_EXPOSURES")
if (!length(OUTCOMES) || any(!grepl("^[A-Za-z][A-Za-z0-9_]*$", OUTCOMES))) stop("Invalid BORDER_OUTCOMES")

validate_ids <- function(d, id, label) {
  if (!is.data.frame(d) || !id %in% names(d)) stop(label, ": missing ID column ", id)
  d <- as.data.frame(d)
  d[[id]] <- as.character(d[[id]])
  if (anyNA(d[[id]]) || any(!nzchar(d[[id]])) || anyDuplicated(d[[id]])) stop(label, ": missing or duplicate IDs")
  if (anyDuplicated(names(d))) stop(label, ": duplicate column names")
  d
}

as_num <- function(x) {
  if (is.numeric(x) || is.logical(x)) return(as.numeric(x))
  suppressWarnings(as.numeric(as.character(x)))
}
as_date <- function(x) {
  if (inherits(x, "Date")) return(x)
  if (inherits(x, "IDate")) return(as.Date(x))
  if (inherits(x, "POSIXt")) return(as.Date(x))
  if (is.factor(x)) x <- as.character(x)
  if (is.character(x)) return(suppressWarnings(as.Date(x, format = "%Y-%m-%d")))
  z <- suppressWarnings(as.numeric(unclass(x))); z[!is.finite(z)] <- NA_real_
  as.Date(z, origin = "1970-01-01")
}
pmin_date <- function(...) {
  z <- lapply(list(...), as_date)
  x <- do.call(pmin, c(z, list(na.rm = TRUE)))
  x <- as_date(x); x[!is.finite(as.numeric(x))] <- as.Date(NA); x
}

read_all <- function() {
  f <- file.path(indir, "Rdata", "all.rds")
  if (!file.exists(f)) stop("Missing ", f)
  validate_ids(readRDS(f), "eid", "all.rds")
}

# UKB variables used here:
# p5057 -> cnt_older   Number of older siblings
# p1873 -> cnt_brother Number of full brothers
# p1883 -> cnt_sister  Number of full sisters
# p1767 -> adopted     Adopted as a child (already generated in all.rds)
clean_count <- function(x) {
  z <- as_num(x)
  z[!is.finite(z) | z < 0 | z != floor(z)] <- NA_real_  # UKB negative special codes, including -1/-3
  z
}

make_birth_order <- function(dat) {
  req <- c("cnt_older","cnt_brother","cnt_sister","adopted")
  miss <- setdiff(req, names(dat)); if (length(miss)) stop("Missing birth-order columns in all.rds: ", paste(miss, collapse=", "))
  dat %>%
    mutate(
      cnt_older   = clean_count(cnt_older),
      cnt_brother = clean_count(cnt_brother),
      cnt_sister  = clean_count(cnt_sister),
      adopted_num = if_else(as_num(adopted) %in% c(0, 1), as_num(adopted), NA_real_),
      sibling_n = cnt_brother + cnt_sister,       # siblings excluding self
      sibship_size = sibling_n + 1,               # all children including self
      birth_rank = cnt_older + 1,
      impossible_rank = !is.na(cnt_older) & !is.na(sibling_n) & cnt_older > sibling_n,
      eligible_border = !is.na(cnt_older) & !is.na(cnt_brother) & !is.na(cnt_sister) &
                        !impossible_rank & (is.na(adopted_num) | adopted_num != 1),
      # Nature Health replication: exactly two children; first-born=0, second-born=1.
      border12 = case_when(
        !eligible_border | sibship_size != 2 ~ NA_real_,
        birth_rank == 1 ~ 0,
        birth_rank == 2 ~ 1,
        TRUE ~ NA_real_
      ),
      # All sibships >=2, retaining only the two extremes; oldest=0, youngest=1.
      border1n = case_when(
        !eligible_border | sibship_size < 2 ~ NA_real_,
        birth_rank == 1 ~ 0,
        birth_rank == sibship_size ~ 1,
        TRUE ~ NA_real_
      ),
      # Ordinal birth rank, per one-position later birth; singleton excluded.
      border_rank = if_else(eligible_border & sibship_size >= 2, birth_rank, NA_real_),
      # 0=oldest, 1=youngest; middle children lie strictly between 0 and 1.
      border_relative = if_else(eligible_border & sibship_size >= 2,
                                (birth_rank - 1)/(sibship_size - 1), NA_real_)
    )
}

prep_covars <- function(dat) {
  # Match current UKB pipeline naming. Center is categorical; age/edu.sco/tdi numeric.
  # Raw UKB edu contains multi-select codes such as "1|2"; it is not an ordinal score.
  # Require the score already derived by ukb/f/phe.R.
  req <- c("age","sex","edu.sco","tdi","center")
  miss <- setdiff(req, names(dat)); if (length(miss)) stop("Missing covariates: ", paste(miss, collapse=", "))
  dat %>% mutate(
    age = as_num(age), sex = factor(sex), edu.sco = as_num(edu.sco), tdi = as_num(tdi), center = factor(center)
  )
}

covars <- c("age","sex","edu.sco","tdi","center")

outcome_date_columns <- function(dat, y) {
  # Use the same hospital ICD10 definition as the UKB/LE8 outcome pipeline when available.
  primary <- paste0("fod_icd10_", y)
  if (primary %in% names(dat)) return(primary)
  n <- names(dat)
  n[grepl("^(fod_|date_)", n) & endsWith(n, paste0("_", y))]
}

validate_follow_end <- function(value, label = "Invalid --follow-end") {
  d <- as_date(value)
  if (!nzchar(value) || !grepl("^[0-9]{4}-[0-9]{2}-[0-9]{2}$", value) ||
      is.na(d) || format(d, "%Y-%m-%d") != value)
    stop(label, ": use a real calendar date in YYYY-MM-DD format.")
  d
}
default_follow_end_path <- file.path(dir0, "scripts", "0f", "phenotype.R")
read_default_follow_end <- function(path = default_follow_end_path) {
  if (!file.exists(path)) stop("Default follow-up source missing: ", path, "; supply --follow-end YYYY-MM-DD")
  first <- readLines(path, n = 1L, warn = FALSE, encoding = "UTF-8")
  expr <- tryCatch(parse(text = first, keep.source = FALSE), error = function(e) NULL)
  # Inspect syntax only; the shared analysis script must never execute here.
  node <- if (length(expr) == 1L) expr[[1L]] else NULL
  valid <- is.call(node) && length(node) == 3L &&
    as.character(node[[1L]]) %in% c("<-", "=") && identical(node[[2L]], as.name("date_follow_end"))
  rhs <- if (valid) node[[3L]] else NULL
  valid <- valid && is.call(rhs) && length(rhs) == 2L &&
    identical(rhs[[1L]], as.name("as.Date")) && is.character(rhs[[2L]]) && length(rhs[[2L]]) == 1L
  if (!valid) stop("Invalid default follow-up assignment on the first line of ", path,
                   "; expected date_follow_end <- as.Date(\"YYYY-MM-DD\")")
  validate_follow_end(rhs[[2L]], paste0("Invalid default follow-up date in ", path))
}
follow_end <- function(default_path = default_follow_end_path) {
  value <- Sys.getenv("BORDER_FOLLOW_END", "")
  if (nzchar(value)) validate_follow_end(value) else read_default_follow_end(default_path)
}
# Validate an explicit override before reading the large input files.
if (nzchar(Sys.getenv("BORDER_FOLLOW_END", ""))) invisible(follow_end())

ensure_outcome <- function(dat, y) {
  tv <- paste0(y, ".t2e"); ev <- paste0(y, ".Yt2e")
  if (all(c(tv, ev) %in% names(dat))) {
    dat[[tv]] <- as_num(dat[[tv]]); dat[[ev]] <- as_num(dat[[ev]])
    bad <- !is.finite(dat[[tv]]) | dat[[tv]] <= 0 | !dat[[ev]] %in% c(0, 1)
    dat[[tv]][bad] <- NA_real_; dat[[ev]][bad] <- NA_real_
    return(dat)
  }
  if (!"date_attend" %in% names(dat)) stop("Cannot reconstruct ", y, ": date_attend missing")
  date_cols <- outcome_date_columns(dat, y)
  if (!length(date_cols)) stop("No dated columns found to reconstruct outcome: ", y)
  event_date <- do.call(pmin_date, dat[date_cols])
  # UKB special dates represent missing/unknown, never actual incident disease.
  event_date[event_date %in% as.Date(c("1900-01-01", "1901-01-01", "1902-02-02", "1903-03-03", "1909-09-09", "2037-07-07"))] <- as.Date(NA)
  base <- as_date(dat$date_attend)
  admin <- follow_end()
  lost <- if ("date_lost" %in% names(dat)) as_date(dat$date_lost) else rep(as.Date(NA), nrow(dat))
  death <- if ("date_death" %in% names(dat)) as_date(dat$date_death) else rep(as.Date(NA), nrow(dat))
  censor <- pmin_date(rep(admin, nrow(dat)), lost, death)
  prevalent <- !is.na(event_date) & !is.na(base) & event_date <= base
  event <- !prevalent & !is.na(event_date) & !is.na(base) & event_date > base & event_date <= censor
  exit <- as.Date(ifelse(event, as.numeric(event_date), as.numeric(censor)), origin="1970-01-01")
  dat[[tv]] <- as.numeric(exit - base)/365.25
  dat[[ev]] <- as.integer(event)
  bad <- prevalent | !is.finite(dat[[tv]]) | dat[[tv]] <= 0
  dat[[tv]][bad] <- NA_real_; dat[[ev]][bad] <- NA_integer_
  dat
}

complete_model_data <- function(dat, vars) {
  miss <- setdiff(vars, names(dat))
  if (length(miss)) stop("Missing model columns: ", paste(miss, collapse=", "))
  d <- as.data.frame(dat[, vars, drop=FALSE])
  ok <- complete.cases(d)
  for (v in names(d)[vapply(d, is.numeric, logical(1))]) ok <- ok & is.finite(d[[v]])
  droplevels(d[ok, , drop=FALSE])
}

# Preserve errors/warnings as diagnostics. Non-converged or infinite coefficients are not reported as estimates.
model_try <- function(expr) {
  notes <- character()
  value <- tryCatch(withCallingHandlers(expr, warning=function(w) {
    notes <<- c(notes, conditionMessage(w)); invokeRestart("muffleWarning")
  }), error=function(e) { notes <<- c(notes, conditionMessage(e)); NULL })
  if (any(grepl("converg|infinite|singular", notes, ignore.case=TRUE))) value <- NULL
  list(fit=value, note=paste(unique(notes), collapse="; "))
}

load_cohort <- function() {
  f <- file.path(rawdir, "border.cohort.rds")
  if (!file.exists(f)) stop("Run border.sh prep first")
  d <- validate_ids(readRDS(f), "eid", "border cohort")
  meta <- attr(d, "border_metadata")
  if (is.null(meta)) stop("Cohort is from an old pipeline; rerun prep")
  input <- file.path(indir, "Rdata", "all.rds")
  info <- file.info(input)
  if (!identical(meta$input, normalizePath(input, mustWork=TRUE)) ||
      !identical(meta$input_size, info$size) || !identical(meta$input_mtime, as.numeric(info$mtime)))
    stop("all.rds changed since prep; rerun prep")
  if (length(setdiff(OUTCOMES, meta$outcomes))) stop("Outcomes changed; rerun prep")
  requested <- Sys.getenv("BORDER_FOLLOW_END", "")
  if (nzchar(requested) && requested != meta$follow_end) stop("Follow-up end changed; rerun prep")
  d
}

fit_cox <- function(dat, exposure, y) {
  dat <- ensure_outcome(dat, y)
  tv <- paste0(y, ".t2e"); ev <- paste0(y, ".Yt2e")
  d <- complete_model_data(dat, unique(c(tv, ev, exposure, covars)))
  row <- tibble(outcome=y, exposure=exposure, N=nrow(d), events=sum(d[[ev]] == 1),
    beta=NA_real_, se=NA_real_, HR=NA_real_, conf.low=NA_real_, conf.high=NA_real_, p=NA_real_,
    status="skipped", note="", unit=switch(exposure,
      border12="second-born vs first-born (two-child families)",
      border1n="youngest vs oldest (all sibships; middle excluded)",
      border_rank="per one-position later birth rank",
      border_relative="youngest (1) vs oldest (0) on relative scale"))
  if (nrow(d) < 500 || row$events < 20 || length(unique(d[[exposure]])) < 2) {
    row$note <- "Need N>=500, events>=20 and exposure variation"; return(row)
  }
  f <- reformulate(c(exposure,covars), response=paste0("Surv(`",tv,"`,`",ev,"`)"))
  fit <- model_try(coxph(f, data=d, ties="efron"))
  row$note <- fit$note
  if (is.null(fit$fit)) { row$status <- "failed"; return(row) }
  sm <- coef(summary(fit$fit))
  if (!exposure %in% rownames(sm) || any(!is.finite(sm[exposure, ]))) {
    row$status <- "failed"; row$note <- "Exposure coefficient is not estimable"; return(row)
  }
  b <- sm[exposure,"coef"]; se <- sm[exposure,"se(coef)"]
  row$beta <- b; row$se <- se; row$HR <- exp(b)
  row$conf.low <- exp(b - qnorm(.975)*se); row$conf.high <- exp(b + qnorm(.975)*se)
  row$p <- sm[exposure,"Pr(>|z|)"]; row$status <- "ok"
  row
}

write_csv2 <- function(x, name) write.csv(x, file.path(rawdir,name), row.names=FALSE, na="")
