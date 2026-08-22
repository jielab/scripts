#!/usr/bin/env Rscript

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 🚩 Packages, paths, and shared inputs
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
`%||%` <- function(a, b) if (!is.null(a) && length(a) > 0 && !all(is.na(a))) a else b

args <- commandArgs(trailingOnly = TRUE)
MODE <- if (length(args) >= 1) args[1] else "all"

CHARLS_CLEAN <- Sys.getenv("CHARLS_CLEAN", "/mnt/d/data/charls/clean")
CHARLS_OUT <- Sys.getenv("CHARLS_OUT", "/mnt/d/analysis/charls")
CHARLS_HARMONIZED <- Sys.getenv("CHARLS_HARMONIZED", file.path(CHARLS_CLEAN, "Harmonized", "H_CHARLS_D_Data.dta"))
SEED <- as.integer(Sys.getenv("SEED", "2026"))
USE_MICE <- Sys.getenv("USE_MICE", "1") %in% c("1", "true", "TRUE", "yes", "YES")
USE_HARMONIZED <- Sys.getenv("USE_HARMONIZED", "1") %in% c("1", "true", "TRUE", "yes", "YES")
FAIL_ON_VAR_MISS <- Sys.getenv("FAIL_ON_VAR_MISS", "1") %in% c("1", "true", "TRUE", "yes", "YES")
ALLOW_HIGH_MISSING <- Sys.getenv("ALLOW_HIGH_MISSING", "0") %in% c("1", "true", "TRUE", "yes", "YES")
START_STEP <- Sys.getenv("START_STEP", "")
CHARLS_MANIFEST <- Sys.getenv("CHARLS_MANIFEST", "")
CHARLS_VERSION <- "v008_explicit_id_mice_robust_cox_extract"

DIR_DAT <- file.path(CHARLS_OUT, "dat")
DIR_FIG <- CHARLS_OUT
DIR_TAB <- CHARLS_OUT
DIR_QC  <- file.path(CHARLS_OUT, "qc")
DIR_LOG <- file.path(CHARLS_OUT, "log")
invisible(lapply(c(CHARLS_OUT, DIR_DAT, DIR_QC, DIR_LOG), dir.create, recursive = TRUE, showWarnings = FALSE))
LOG_FILE <- file.path(DIR_LOG, "charls.R.log")

msg <- function(...) {
  z <- paste0(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " | ", paste0(..., collapse = ""))
  cat(z, "\n")
  cat(z, "\n", file = LOG_FILE, append = TRUE)
}

need_pkg <- function(pkgs, required = TRUE) {
  miss <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(miss) && required) {
    stop("Missing R package(s): ", paste(miss, collapse = ", "),
         "\nInstall example:\ninstall.packages(c(",
         paste(sprintf('"%s"', miss), collapse = ", "), "))", call. = FALSE)
  }
  invisible(!length(miss))
}

need_pkg(c("data.table", "haven", "survival", "ggplot2", "openxlsx"), required = TRUE)
suppressPackageStartupMessages({
  library(data.table)
  library(haven)
  library(survival)
  library(ggplot2)
  library(openxlsx)
})

opt_pkg <- function(pkg) requireNamespace(pkg, quietly = TRUE)
set.seed(SEED)

# Global variable map and QC logs.
VAR_MAP <- data.table(wave = character(), concept = character(), var = character(), source = character(), label = character(), required = logical(), status = character())
FI_MAP <- data.table(wave = integer(), item = character(), var = character(), status = character(), missing_prop = numeric())
FILE_CHOICE <- data.table(wave = character(), file_role = character(), path = character(), archive_rel = character(), internal_file = character(), action = character(), size_bytes = numeric(), note = character())
ID_MATCH_QC <- data.table(wave = character(), file = character(), id_source = character(), id_transform = character(), overlap = integer(), raw_n = integer(), harmonized_n = integer(), note = character())
MANIFEST_DF <- NULL

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 🚩 Low-level utilities
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
qsave <- function(x, f) { dir.create(dirname(f), recursive = TRUE, showWarnings = FALSE); saveRDS(x, f); invisible(f) }
qread <- function(f) { if (!file.exists(f)) stop("Missing file: ", f, ". Run `bash charls.sh prep` first.", call. = FALSE); readRDS(f) }
write_xlsx <- function(x, f) {
  dir.create(dirname(f), recursive = TRUE, showWarnings = FALSE)
  if (inherits(x, "data.table") || is.data.frame(x)) {
    openxlsx::write.xlsx(as.data.frame(x), f, overwrite = TRUE)
  } else if (is.list(x)) {
    wb <- openxlsx::createWorkbook()
    for (nm in names(x)) {
      sh <- substr(gsub("[\\[\\]\\*\\?/\\\\:]", "_", nm), 1, 31)
      openxlsx::addWorksheet(wb, sh)
      openxlsx::writeData(wb, sh, as.data.frame(x[[nm]]))
    }
    openxlsx::saveWorkbook(wb, f, overwrite = TRUE)
  } else stop("Unsupported object for write_xlsx")
  invisible(f)
}

fmt_p <- function(p) {
  ifelse(is.na(p), "", ifelse(p < 0.001, "<0.001", sprintf("%.3f", p)))
}
fmt_num <- function(x, d = 2) {
  ifelse(is.na(x), "", sprintf(paste0("%0.", d, "f"), x))
}
fmt_hr <- function(hr, lo, hi, p = NULL) {
  out <- sprintf("%.2f (%.2f-%.2f)", hr, lo, hi)
  if (!is.null(p)) out <- paste0(out, "; P=", fmt_p(p))
  out
}

clean_names <- function(nm) {
  nm <- gsub("\\s+", "_", nm)
  nm <- gsub("[^A-Za-z0-9_]+", "_", nm)
  nm <- gsub("_+", "_", nm)
  nm
}

as_dt <- function(x) data.table::as.data.table(x)

read_dta_dt <- function(f, n_max = Inf) {
  msg("READ ", f)
  x <- haven::read_dta(f, n_max = n_max)
  data.table::as.data.table(x)
}

vlabel <- function(x) {
  z <- attr(x, "label", exact = TRUE)
  if (is.null(z) || length(z) == 0 || is.na(z)) "" else as.character(z)[1]
}
labels_dt <- function(dt) data.table(var = names(dt), label = vapply(dt, function(x) vlabel(x), character(1)))

num <- function(x) {
  z <- suppressWarnings(as.numeric(as.vector(x)))
  z[z %in% c(-99, -98, -97, -9, -8, -7, -2, -1, 97, 98, 99, 997, 998, 999, 9997, 9998, 9999)] <- NA_real_
  z
}

mode_value <- function(x) {
  x <- x[!is.na(x)]
  if (!length(x)) return(NA)
  names(sort(table(x), decreasing = TRUE))[1]
}

binary01 <- function(x, name = "") {
  z <- num(x)
  u <- sort(unique(z[!is.na(z)]))
  out <- rep(NA_real_, length(z))
  if (!length(u)) return(out)
  if (all(u %in% c(0, 1))) {
    out[z == 1] <- 1; out[z == 0] <- 0
  } else if (all(u %in% c(1, 2))) {
    out[z == 1] <- 1; out[z == 2] <- 0
  } else {
    out[z == 1] <- 1
    out[z %in% c(0, 2, 5)] <- 0
  }
  out
}

row_any01 <- function(...) {
  m <- cbind(...)
  if (!is.matrix(m)) m <- matrix(m, ncol = 1)
  apply(m, 1, function(x) {
    if (any(x == 1, na.rm = TRUE)) return(1)
    if (all(is.na(x))) return(NA_real_)
    0
  })
}

row_sum_na <- function(...) {
  m <- cbind(...)
  if (!is.matrix(m)) m <- matrix(m, ncol = 1)
  apply(m, 1, function(x) if (all(is.na(x))) NA_real_ else sum(x == 1, na.rm = TRUE))
}

make_quartile <- function(x, prefix = "Q") {
  q <- quantile(x, probs = seq(0, 1, 0.25), na.rm = TRUE, names = FALSE, type = 7)
  q <- unique(q)
  if (length(q) < 5) {
    r <- rank(x, ties.method = "first", na.last = "keep")
    g <- cut(r, breaks = quantile(r, seq(0, 1, 0.25), na.rm = TRUE), include.lowest = TRUE, labels = paste0(prefix, 1:4))
  } else {
    g <- cut(x, breaks = q, include.lowest = TRUE, labels = paste0(prefix, 1:4))
  }
  factor(g, levels = paste0(prefix, 1:4))
}

winsor <- function(x, p = c(0.01, 0.99)) {
  q <- quantile(x, p, na.rm = TRUE)
  pmin(pmax(x, q[1]), q[2])
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 🚩 Variable resolving from Harmonized/raw CHARLS
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
append_map <- function(wave, concept, var, source, label, required, status) {
  VAR_MAP <<- rbind(VAR_MAP, data.table(wave = as.character(wave), concept = concept, var = var %||% NA_character_, source = source %||% NA_character_, label = label %||% "", required = required, status = status), fill = TRUE)
}

manual_var <- function(wave, concept) {
  # Optional user override, e.g. CHARLS_VAR_W1_TG=raw2011_blood_20140429_bl_tg
  keys <- unique(c(
    paste0("CHARLS_VAR_W", wave, "_", toupper(gsub("[^A-Za-z0-9]+", "_", concept))),
    paste0("CHARLS_", toupper(gsub("[^A-Za-z0-9]+", "_", concept)), "_W", wave)
  ))
  for (k in keys) {
    v <- Sys.getenv(k, "")
    if (nzchar(v)) return(v)
  }
  ""
}

write_var_debug <- function(dt, wave, concept) {
  # When a required variable cannot be resolved, write all names/labels so the
  # next edit can be based on facts rather than guessing.
  z <- labels_dt(dt)
  z[, text := tolower(paste(var, label))]
  key <- tolower(concept)
  z[, likely := grepl(key, text) |
       grepl("tg|trig|triglycer|甘油", text) |
       grepl("hdl|high.?density|高密度", text) |
       grepl("ldl|low.?density|低密度", text) |
       grepl("chol|cholesterol|胆固醇", text) |
       grepl("glucose|glu|血糖", text) |
       grepl("hba1c|a1c|glyc", text)]
  f <- file.path(DIR_QC, paste0("debug_available_variables_wave", wave, "_", concept, ".xlsx"))
  suppressWarnings(write_xlsx(list(likely = z[likely == TRUE], all_variables = z), f))
  msg("DEBUG wrote available variable list: ", f)
}

find_var <- function(dt, wave, concept, candidates = character(), regex = NULL, exclude_regex = NULL, required = FALSE) {
  nm <- names(dt)
  mv <- manual_var(wave, concept)
  if (nzchar(mv)) {
    if (mv %in% nm) {
      append_map(wave, concept, mv, "manual_env", vlabel(dt[[mv]]), required, "ok")
      return(mv)
    } else {
      append_map(wave, concept, mv, "manual_env", "", required, "manual_var_not_found")
      if (required && FAIL_ON_VAR_MISS) {
        write_var_debug(dt, wave, concept)
        stop("Manual variable override not found for wave=", wave, ", concept=", concept, ": ", mv, call. = FALSE)
      }
    }
  }
  for (v in candidates) {
    if (!is.na(v) && nzchar(v) && v %in% nm) {
      append_map(wave, concept, v, "candidate", vlabel(dt[[v]]), required, "ok")
      return(v)
    }
  }
  lab <- labels_dt(dt)
  lab[, text := tolower(paste(var, label))]
  lab[, score := 0]
  if (!is.null(regex) && nzchar(regex)) inc <- grepl(regex, lab$text, ignore.case = TRUE, perl = TRUE) else inc <- rep(FALSE, nrow(lab))

  # Extra concept-specific rescue keywords for common CHARLS blood variables.
  ck <- toupper(concept)
  rescue_regex <- switch(ck,
    "TG"    = "tg|trig|triglycer|triacyl|甘油",
    "HDLC"  = "hdl|high[ _-]?density|high density lipoprotein|高密度",
    "LDLC"  = "ldl|low[ _-]?density|low density lipoprotein|低密度",
    "TC"    = "(^|[_ -])(tc|cho|chol)($|[_ -])|total[ _-]?chol|cholesterol|胆固醇",
    "FBG"   = "fbg|fglu|glucose|blood[ _-]?sugar|glu|血糖",
    "HBA1C" = "hba1c|hbalc|a1c|glycated|glycosylated|glyhb",
    ""
  )
  if (nzchar(rescue_regex)) inc <- inc | grepl(rescue_regex, lab$text, ignore.case = TRUE, perl = TRUE)

  if (!is.null(exclude_regex) && nzchar(exclude_regex)) inc <- inc & !grepl(exclude_regex, lab$text, ignore.case = TRUE, perl = TRUE)
  if (any(inc)) {
    lab[inc, score := score + 1]
    lab[inc & grepl(paste0("^r", wave), var, ignore.case = TRUE), score := score + 8]
    yrs <- wave_raw_years[[as.character(wave)]] %||% character()
    if (length(yrs)) lab[inc & grepl(paste0("^raw(", paste(yrs, collapse = "|"), ")_"), var, ignore.case = TRUE), score := score + 7]
    # Strongly prefer compact biomarker variable names over descriptions in unrelated labels.
    lab[inc & grepl("(^|_)(bl_)?(tg|trig|triglycer|hdl|hdlc|ldl|ldlc|cho|chol|tc|glu|glucose|hba1c|hbalc)(_|$)", var, ignore.case = TRUE), score := score + 6]
    lab[inc & grepl("blood|biomarker", var, ignore.case = TRUE), score := score + 2]
    lab[inc & grepl("(_|^)(e|s|l|a)$", var, ignore.case = TRUE), score := score + 1]
    lab[inc & grepl("flag|miss|imp|weight|wt|psu|strata|date|month|year|proxy|version|code", text, ignore.case = TRUE), score := score - 5]
    # Prefer variables with at least some nonmissing values when multiple candidates exist.
    cand_vars <- lab[inc][order(-score, nchar(var))]$var
    nn <- vapply(cand_vars, function(v) sum(!is.na(val_num(dt, v))), integer(1))
    z <- lab[var %in% cand_vars]
    z[, nonmissing_n := nn[match(var, cand_vars)]]
    z[, score2 := score + fifelse(nonmissing_n > 0, 3, 0)]
    hit <- z[order(-score2, nchar(var))][1]
    append_map(wave, concept, hit$var, "regex_rescue", hit$label, required, paste0("ok; nonmissing_n=", hit$nonmissing_n))
    return(hit$var)
  }
  append_map(wave, concept, NA_character_, NA_character_, NA_character_, required, if (required) "missing" else "not_found_optional")
  if (required && FAIL_ON_VAR_MISS) {
    write_var_debug(dt, wave, concept)
    suppressWarnings(write_xlsx(VAR_MAP, file.path(DIR_QC, "var_map.partial_on_error.xlsx")))
    stop("Required variable not found for wave=", wave, ", concept=", concept,
         ". See qc/debug_available_variables_wave", wave, "_", concept, ".xlsx and qc/var_map.partial_on_error.xlsx", call. = FALSE)
  }
  NA_character_
}

val <- function(dt, var) if (!is.na(var) && var %in% names(dt)) dt[[var]] else rep(NA_real_, nrow(dt))
val_num <- function(dt, var) num(val(dt, var))
val_bin <- function(dt, var) binary01(val(dt, var), var)

id_var <- function(dt) {
  cands <- c("ID", "id", "hhidpn", "HHIDPN", "mergeid", "MERGEID", "personid", "PERSONID", "pid", "PID", "individualID", "IndividualID")
  hit <- cands[cands %in% names(dt)][1]
  if (!is.na(hit)) return(hit)
  hit <- names(dt)[grepl("^(id|.*id$|hhidpn|merge)", names(dt), ignore.case = TRUE)][1]
  if (is.na(hit)) stop("Cannot identify participant ID variable in Harmonized CHARLS file.", call. = FALSE)
  hit
}

wave_prefix <- function(w) paste0("r", w)
wv_year <- c(`1` = 2011, `2` = 2013, `3` = 2015, `4` = 2018, `5` = 2020)

wave_raw_years <- list(`1` = "2011", `2` = "2013", `3` = "2015", `4` = "2018", `5` = "2020")

cand <- function(w, stems) {
  p <- wave_prefix(w)
  yrs <- wave_raw_years[[as.character(w)]] %||% character()
  raw <- unlist(lapply(yrs, function(y) c(paste0("raw", y, "_", stems), paste0("raw", y, "_bl_", stems))))
  unique(c(paste0(p, stems), paste0(p, "_", stems), raw, stems))
}

resolve_wave_vars <- function(dt, w) {
  p <- wave_prefix(w)
  list(
    age = find_var(dt, w, "age", cand(w, c("agey", "age", "ageyr")), paste0("^", p, ".*(agey|age).*"), required = (w == 1)),
    gender = find_var(dt, w, "gender", c("ragender", "gender", "sex", cand(w, c("gender", "sex"))), "(ragender|gender| sex$|^sex$)", required = (w == 1)),
    birth_year = find_var(dt, w, "birth_year", c("rabyear", "birthyear", "birth_year", "ba002_1"), "(birth.*year|rabyear)", required = FALSE),
    smoking = find_var(dt, w, "smoking", cand(w, c("smoken", "smokev", "smoke", "smokef", "smokec")), paste0("^", p, ".*smok"), required = FALSE),
    drinking = find_var(dt, w, "drinking", cand(w, c("drinkn", "drink", "drinkev", "drinkd")), paste0("^", p, ".*drink|alcohol"), required = FALSE),
    bmi = find_var(dt, w, "BMI", cand(w, c("mbmi", "bmi", "bmi_m")), paste0("^", p, ".*(bmi|body mass)"), required = FALSE),
    sbp = find_var(dt, w, "SBP", cand(w, c("systo", "sbp", "systolic", "systolic_mean")), paste0("^", p, ".*(systo|sbp|systolic)"), required = FALSE),
    dbp = find_var(dt, w, "DBP", cand(w, c("diasto", "dbp", "diastolic", "diastolic_mean")), paste0("^", p, ".*(diasto|dbp|diastolic)"), required = FALSE),
    fbg = find_var(dt, w, "FBG", cand(w, c("fglu", "glucose", "glu", "bloodglucose", "bl_glu", "gluc")), paste0("(^|_)(bl_)?(glu|glucose|fglu)($|_)|^", p, ".*(glucose|fglu|blood sugar|glu)"), "urine|ketone|self|diagnos", required = FALSE),
    hba1c = find_var(dt, w, "HbA1c", cand(w, c("hba1c", "hbalc", "bl_hbalc", "glyhb", "ghb", "a1c")), paste0("(^|_)(bl_)?(hba1c|hbalc|a1c)($|_)|^", p, ".*(hba1c|glycated|glycosylated|a1c|glyhb|hbalc)"), required = FALSE),
    tc = find_var(dt, w, "TC", cand(w, c("tc", "cho", "chol", "tchol", "totalchol", "cholesterol", "bl_cho", "bl_chol")), paste0("(^|_)(bl_)?(tc|cho|chol|tchol|totalchol|cholesterol)($|_)|^", p, ".*(total cholesterol|cholesterol|tchol|tc)"), "hdl|ldl|medicine|diagnos|ratio", required = FALSE),
    tg = find_var(dt, w, "TG", cand(w, c("tg", "trig", "trigly", "triglycer", "triglycerides", "triglyceride", "bl_tg", "bltg", "bl_trig", "bl_trigly", "bl_triglycerides")), paste0("(^|_)(bl_?)?(tg|trig|trigly|triglycer|triglyceride|triglycerides)($|_)|^", p, ".*(triglycer|trig| tg |甘油)"), required = (w %in% c(1, 3))),
    hdl = find_var(dt, w, "HDLC", cand(w, c("hdl", "hdlc", "hdl_c", "hdchol", "hdlchol", "bl_hdl", "blhdl", "bl_hdl_c", "bl_hdlchol")), paste0("(^|_)(bl_?)?(hdl|hdlc|hdl_c|hdchol|hdlchol)($|_)|^", p, ".*(hdl|high.density|high density|高密度)"), required = (w %in% c(1, 3))),
    ldl = find_var(dt, w, "LDLC", cand(w, c("ldl", "ldlc", "ldl_c", "ldchol", "bl_ldl")), paste0("(^|_)(bl_)?(ldl|ldlc|ldl_c|ldchol)($|_)|^", p, ".*(ldl|low.density|low density)"), required = FALSE),
    hibp_self = find_var(dt, w, "self_hypertension", cand(w, c("hibpe", "hibp", "hbp", "highblood", "hypertension")), paste0("^", p, ".*(hibp|hypertension|high blood pressure)"), "rx|med|treat|medicine|drug", required = FALSE),
    diab_self = find_var(dt, w, "self_diabetes", cand(w, c("diabe", "diab", "diabetes")), paste0("^", p, ".*(diab|diabetes|hyperglyc)"), "rx|med|treat|medicine|drug", required = FALSE),
    heart_self = find_var(dt, w, "self_heart", cand(w, c("hearte", "heart", "heartd", "heartpr")), paste0("^", p, ".*(heart|angina|coronary|cardiac)"), "rx|med|treat|medicine|drug", required = FALSE),
    stroke_self = find_var(dt, w, "self_stroke", cand(w, c("stroke", "strokee")), paste0("^", p, ".*stroke"), "rx|med|treat|medicine|drug", required = FALSE),
    dyslip_self = find_var(dt, w, "self_dyslipidemia", cand(w, c("dyslipe", "dyslip", "lipid", "lipo")), paste0("^", p, ".*(dyslip|lipid disorder|cholesterol problem|hyperlip)"), "rx|med|treat|medicine|drug", required = FALSE),
    rx_htn = find_var(dt, w, "rx_hypertension", cand(w, c("rxhibp", "hibprx", "rxhbp", "bpmed", "htnmed")), paste0("^", p, ".*(rx|med|treat|medicine|drug).*(hibp|hypertension|blood pressure)|^", p, ".*(hibp|hypertension).*(rx|med|treat|medicine|drug)"), required = FALSE),
    rx_diab = find_var(dt, w, "rx_diabetes", cand(w, c("rxdiab", "diabrx", "rxdiabe", "diabmed")), paste0("^", p, ".*(rx|med|treat|medicine|drug).*(diab|glucose)|^", p, ".*(diab|glucose).*(rx|med|treat|medicine|drug)"), required = FALSE),
    rx_lipid = find_var(dt, w, "rx_lipid", cand(w, c("rxdyslip", "rxlipid", "lipidmed", "cholmed")), paste0("^", p, ".*(rx|med|treat|medicine|drug).*(lipid|cholesterol)|^", p, ".*(lipid|cholesterol).*(rx|med|treat|medicine|drug)"), required = FALSE),
    rx_heart = find_var(dt, w, "rx_heart", cand(w, c("rxheart", "heartrx", "heartmed", "cvdmed")), paste0("^", p, ".*(rx|med|treat|medicine|drug).*(heart|cardiac|angina)|^", p, ".*(heart|cardiac|angina).*(rx|med|treat|medicine|drug)"), required = FALSE),
    rx_stroke = find_var(dt, w, "rx_stroke", cand(w, c("rxstroke", "strokemed")), paste0("^", p, ".*(rx|med|treat|medicine|drug).*stroke|^", p, ".*stroke.*(rx|med|treat|medicine|drug)"), required = FALSE)
  )
}

convert_glucose <- function(x) {
  # CHARLS blood glucose in the paper is mg/dL. If values look like mmol/L, convert.
  med <- suppressWarnings(median(x, na.rm = TRUE))
  if (is.finite(med) && med > 0 && med < 30) x * 18 else x
}
convert_lipid <- function(x) {
  # CHARLS lipid values in the paper are mg/dL. If total cholesterol/HDL/TG look like mmol/L, convert roughly.
  med <- suppressWarnings(median(x, na.rm = TRUE))
  if (is.finite(med) && med > 0 && med < 20) x * 38.67 else x
}
convert_tg <- function(x) {
  med <- suppressWarnings(median(x, na.rm = TRUE))
  if (is.finite(med) && med > 0 && med < 20) x * 88.57 else x
}

aip_from_mgdl <- function(tg_mgdl, hdl_mgdl) {
  ifelse(!is.na(tg_mgdl) & !is.na(hdl_mgdl) & tg_mgdl > 0 & hdl_mgdl > 0,
         log10((tg_mgdl / 88.57) / (hdl_mgdl / 38.67)), NA_real_)
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 🚩 Frailty-index construction
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
resolve_fi_vars <- function(dt, w) {
  p <- wave_prefix(w)
  list(
    cancer = find_var(dt, w, "FI_cancer", cand(w, c("cancre", "cancer", "cancr")), paste0("^", p, ".*(cancer|malignant|tumor)"), required = FALSE),
    arthritis = find_var(dt, w, "FI_arthritis", cand(w, c("arthre", "arthritis", "arthr")), paste0("^", p, ".*(arthritis|rheumatism)"), required = FALSE),
    lung = find_var(dt, w, "FI_lung", cand(w, c("lunge", "lung", "clunge")), paste0("^", p, ".*(lung|pulmonary|emphysema|bronchitis)"), "asthma", required = FALSE),
    asthma = find_var(dt, w, "FI_asthma", cand(w, c("asthmae", "asthma")), paste0("^", p, ".*asthma"), required = FALSE),
    psych = find_var(dt, w, "FI_psych", cand(w, c("psyche", "psych", "psychiatric", "emotional")), paste0("^", p, ".*(psychiatric|emotional|nervous|psych)"), required = FALSE),
    memory_dx = find_var(dt, w, "FI_memory_disease", cand(w, c("memrye", "memorye", "alzhe", "dementia", "memory")), paste0("^", p, ".*(memory.related|memory disease|dementia|alzheimer)"), "test|score|recall|word", required = FALSE),
    vision = find_var(dt, w, "FI_vision", cand(w, c("sight", "vision", "eyesight")), paste0("^", p, ".*(sight|vision|eyesight)"), required = FALSE),
    hearing = find_var(dt, w, "FI_hearing", cand(w, c("hear", "hearing")), paste0("^", p, ".*(hearing|hear)"), required = FALSE),
    health = find_var(dt, w, "FI_general_health", cand(w, c("shlt", "selfhealth", "health")), paste0("^", p, ".*(self.*health|general health|shlt)"), "insurance|care|history|disease", required = FALSE),
    dress = find_var(dt, w, "FI_dressing", cand(w, c("dressa", "dress", "dressing")), paste0("^", p, ".*(dress|dressing)"), required = FALSE),
    bath = find_var(dt, w, "FI_bathing", cand(w, c("batha", "bath", "bathing", "shower")), paste0("^", p, ".*(bath|shower)"), required = FALSE),
    eat = find_var(dt, w, "FI_eating", cand(w, c("eata", "eat", "eating")), paste0("^", p, ".*(eat|eating)"), required = FALSE),
    bed = find_var(dt, w, "FI_bed", cand(w, c("beda", "bed", "bedinout")), paste0("^", p, ".*(bed|get.*out.*bed)"), required = FALSE),
    toilet = find_var(dt, w, "FI_toilet", cand(w, c("toilta", "toilet", "toileting")), paste0("^", p, ".*(toilet)"), required = FALSE),
    money = find_var(dt, w, "FI_money", cand(w, c("moneya", "money", "managemoney")), paste0("^", p, ".*(money|financ)"), required = FALSE),
    meds = find_var(dt, w, "FI_medications", cand(w, c("medsa", "meds", "takemed", "medicine")), paste0("^", p, ".*(taking medication|take.*med|medicine)"), required = FALSE),
    shop = find_var(dt, w, "FI_shopping", cand(w, c("shopa", "shop", "shopping")), paste0("^", p, ".*(shop|groceries)"), required = FALSE),
    meals = find_var(dt, w, "FI_meals", cand(w, c("mealsa", "meal", "meals", "cook")), paste0("^", p, ".*(meal|cook|prepar.*food)"), required = FALSE),
    housework = find_var(dt, w, "FI_housework", cand(w, c("housewka", "housework", "housewk", "house")), paste0("^", p, ".*(housework|house work|household chores)"), required = FALSE),
    walk100 = find_var(dt, w, "FI_walk100", cand(w, c("walk100a", "walk100", "walk")), paste0("^", p, ".*(walk.*100|100.*yard|walk)"), required = FALSE),
    chair = find_var(dt, w, "FI_chair", cand(w, c("chaira", "chair")), paste0("^", p, ".*(chair|get.*up)"), required = FALSE),
    climb = find_var(dt, w, "FI_climb", cand(w, c("climsa", "climb", "stairs")), paste0("^", p, ".*(climb|stairs)"), required = FALSE),
    lift = find_var(dt, w, "FI_lift", cand(w, c("lifta", "lift", "carry")), paste0("^", p, ".*(lift|carry|10.*pound|10.*jin)"), required = FALSE),
    dime = find_var(dt, w, "FI_coin", cand(w, c("dimea", "dime", "coin", "pickup")), paste0("^", p, ".*(coin|dime|pick.*up)"), required = FALSE),
    stoop = find_var(dt, w, "FI_stoop", cand(w, c("stoopa", "stoop", "kneel", "crouch")), paste0("^", p, ".*(stoop|kneel|crouch)"), required = FALSE),
    arms = find_var(dt, w, "FI_arms", cand(w, c("armsa", "arms", "reach")), paste0("^", p, ".*(arm|reach.*shoulder|above.*shoulder)"), required = FALSE),
    cesd10 = find_var(dt, w, "FI_CESD10", cand(w, c("cesd10", "cesd", "depress", "depression")), paste0("^", p, ".*(cesd|depress)"), required = FALSE),
    cog_total = find_var(dt, w, "FI_cognition_total", cand(w, c("cogtot", "cognition", "cog", "mental", "memory_orientation")), paste0("^", p, ".*(cog.*score|cognition|mental status|memory.*orientation)"), "disease|diagnos|problem", required = FALSE),
    orient = find_var(dt, w, "FI_orientation", cand(w, c("orient", "orientation")), paste0("^", p, ".*(orientation|date naming|serial)"), required = FALSE),
    memory_score = find_var(dt, w, "FI_memory_score", cand(w, c("imrc", "dlrc", "tr20", "wordrecall", "memoryscore", "recall")), paste0("^", p, ".*(word recall|memory score|immediate recall|delayed recall|recall)"), "disease|diagnos", required = FALSE)
  )
}

ordinal_poor_fair <- function(x) {
  z <- num(x); out <- rep(NA_real_, length(z))
  out[z %in% c(4, 5)] <- 1; out[z %in% c(1, 2, 3)] <- 0
  # If already binary difficulty-like coding.
  u <- sort(unique(z[!is.na(z)]))
  if (length(u) && all(u %in% c(0, 1))) out <- z
  out
}

ordinal_problem <- function(x) {
  z <- num(x); out <- rep(NA_real_, length(z))
  u <- sort(unique(z[!is.na(z)]))
  if (!length(u)) return(out)
  if (all(u %in% c(0, 1))) return(z)
  if (all(u %in% c(1, 2))) return(ifelse(z == 1, 1, ifelse(z == 2, 0, NA_real_)))
  # For health/vision/hearing scales from excellent/good to poor, code fair/poor as deficit.
  out[z >= 4 & z <= 5] <- 1; out[z >= 1 & z <= 3] <- 0
  out
}

construct_fi <- function(dt, w, derived) {
  fv <- resolve_fi_vars(dt, w)
  n <- nrow(dt)
  items <- list(
    hypertension = derived$hypertension,
    diabetes = derived$diabetes,
    heart_disease = derived$heart_disease,
    stroke = derived$stroke,
    cancer = val_bin(dt, fv$cancer),
    arthritis = val_bin(dt, fv$arthritis),
    chronic_lung = val_bin(dt, fv$lung),
    asthma = val_bin(dt, fv$asthma),
    psychiatric = val_bin(dt, fv$psych),
    memory_disease = val_bin(dt, fv$memory_dx),
    vision_problem = ordinal_problem(val(dt, fv$vision)),
    hearing_problem = ordinal_problem(val(dt, fv$hearing)),
    poor_fair_health = ordinal_poor_fair(val(dt, fv$health)),
    dressing = val_bin(dt, fv$dress),
    bathing = val_bin(dt, fv$bath),
    eating = val_bin(dt, fv$eat),
    bed = val_bin(dt, fv$bed),
    toilet = val_bin(dt, fv$toilet),
    money = val_bin(dt, fv$money),
    medications = val_bin(dt, fv$meds),
    shopping = val_bin(dt, fv$shop),
    meals = val_bin(dt, fv$meals),
    housework = val_bin(dt, fv$housework),
    walk100 = val_bin(dt, fv$walk100),
    chair = val_bin(dt, fv$chair),
    climb_stairs = val_bin(dt, fv$climb),
    lift = val_bin(dt, fv$lift),
    coin = val_bin(dt, fv$dime),
    stoop = val_bin(dt, fv$stoop),
    arms = val_bin(dt, fv$arms)
  )
  cesd <- val_num(dt, fv$cesd10)
  items$depression <- ifelse(is.na(cesd), NA_real_, ifelse(cesd > 10, 1, 0))
  cog <- val_num(dt, fv$cog_total)
  if (all(is.na(cog))) {
    mem <- val_num(dt, fv$memory_score); ori <- val_num(dt, fv$orient)
    cog <- ifelse(is.na(mem) & is.na(ori), NA_real_, rowSums(cbind(mem, ori), na.rm = TRUE))
  }
  # Deficit is poor cognition. The paper describes cognition score / 14 as continuous 0-1;
  # frailty deficit should increase as cognition declines, so use 1 - score/14.
  items$cognition <- ifelse(is.na(cog), NA_real_, pmin(pmax(1 - cog / 14, 0), 1))

  mat <- as.data.frame(items)
  stopifnot(ncol(mat) == 32)
  n_nonmiss <- rowSums(!is.na(mat))
  fi <- ifelse(n_nonmiss >= 24, rowMeans(mat, na.rm = TRUE) * 100, NA_real_)
  tmp <- data.table(item = names(mat), var = c("derived", "derived", "derived", "derived", unlist(fv[c("cancer", "arthritis", "lung", "asthma", "psych", "memory_dx", "vision", "hearing", "health", "dress", "bath", "eat", "bed", "toilet", "money", "meds", "shop", "meals", "housework", "walk100", "chair", "climb", "lift", "dime", "stoop", "arms", "cesd10", "cog_total")])), status = "ok")
  tmp[, wave := w]
  tmp[, missing_prop := vapply(mat, function(x) mean(is.na(x)), numeric(1))]
  FI_MAP <<- rbind(FI_MAP, tmp[, .(wave, item, var, status, missing_prop)], fill = TRUE)
  list(fi = fi, n_nonmiss = n_nonmiss, items = mat)
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 🚩 Build per-wave analytic variables
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
normal_gender <- function(x) {
  z <- num(x); out <- rep(NA_character_, length(z))
  # Harmonized and raw CHARLS commonly use 1=male, 2=female.
  out[z == 1] <- "Male"; out[z == 2] <- "Female"
  # Some harmonized variables may use 0=female, 1=male.
  u <- sort(unique(z[!is.na(z)]))
  if (length(u) && all(u %in% c(0, 1))) { out[z == 1] <- "Male"; out[z == 0] <- "Female" }
  factor(out, levels = c("Female", "Male"))
}

build_wave <- function(h, w) {
  year <- unname(wv_year[as.character(w)])
  idv <- id_var(h)
  rv <- resolve_wave_vars(h, w)
  out <- data.table(id = as.character(h[[idv]]))
  age <- val_num(h, rv$age)
  by <- val_num(h, rv$birth_year)
  if (all(is.na(age)) && !all(is.na(by))) age <- year - by
  out[, age := age]
  out[, gender := normal_gender(val(h, rv$gender))]
  out[, smoking := val_bin(h, rv$smoking)]
  out[, drinking := val_bin(h, rv$drinking)]
  out[, BMI := val_num(h, rv$bmi)]
  out[, SBP := val_num(h, rv$sbp)]
  out[, DBP := val_num(h, rv$dbp)]
  out[, FBG := convert_glucose(val_num(h, rv$fbg))]
  out[, HbA1c := val_num(h, rv$hba1c)]
  out[, TC := convert_lipid(val_num(h, rv$tc))]
  out[, TG := convert_tg(val_num(h, rv$tg))]
  out[, HDLC := convert_lipid(val_num(h, rv$hdl))]
  out[, LDLC := convert_lipid(val_num(h, rv$ldl))]
  out[, rx_htn := val_bin(h, rv$rx_htn)]
  out[, rx_diab := val_bin(h, rv$rx_diab)]
  out[, rx_lipid := val_bin(h, rv$rx_lipid)]
  out[, rx_heart := val_bin(h, rv$rx_heart)]
  out[, rx_stroke := val_bin(h, rv$rx_stroke)]

  self_htn <- val_bin(h, rv$hibp_self)
  self_diab <- val_bin(h, rv$diab_self)
  self_heart <- val_bin(h, rv$heart_self)
  self_stroke <- val_bin(h, rv$stroke_self)
  self_dyslip <- val_bin(h, rv$dyslip_self)

  out[, hypertension := row_any01(self_htn, rx_htn, ifelse(!is.na(SBP) & SBP >= 140, 1, 0), ifelse(!is.na(DBP) & DBP >= 90, 1, 0))]
  out[, diabetes := row_any01(self_diab, rx_diab, ifelse(!is.na(FBG) & FBG >= 126, 1, 0), ifelse(!is.na(HbA1c) & HbA1c >= 6.5, 1, 0))]
  out[, heart_disease := row_any01(self_heart, rx_heart)]
  out[, stroke := row_any01(self_stroke, rx_stroke)]
  out[, dyslipidemia := row_any01(self_dyslip, rx_lipid, ifelse(!is.na(TG) & TG >= 150, 1, 0), ifelse(!is.na(TC) & TC >= 240, 1, 0), ifelse(!is.na(HDLC) & HDLC < 40, 1, 0), ifelse(!is.na(LDLC) & LDLC >= 160, 1, 0))]
  out[, n_cmd := row_sum_na(diabetes, heart_disease, stroke)]
  out[, CMM := ifelse(is.na(n_cmd), NA_real_, as.numeric(n_cmd >= 2))]
  out[, AIP := aip_from_mgdl(TG, HDLC)]

  fi <- construct_fi(h, w, out)
  out[, FI := fi$fi]
  out[, FI_items_nonmissing := fi$n_nonmiss]
  out[, AIPFI := AIP * FI]
  out[, wave := w]
  out[, year := year]
  out[]
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 🚩 Dataset building, missingness, and imputation
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
load_harmonized <- function() {
  if (!file.exists(CHARLS_HARMONIZED)) {
    stop("Harmonized CHARLS file not found: ", CHARLS_HARMONIZED,
         "\nThis reproduction workflow is designed to use Harmonized/H_CHARLS_D_Data.dta as the primary source.", call. = FALSE)
  }
  read_dta_dt(CHARLS_HARMONIZED)
}


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 🚩 Manifest-aware raw CHARLS fallback
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# The extraction manifest is optional at runtime. It is used only to select the
# canonical extracted files and to avoid renamed_conflict files such as
# 2013/PSU__PSU.dta. If no manifest is present, the code falls back to a small
# hard-coded list of canonical Blood/Biomarker files and then to filename search.
raw_file_regex <- "(?i)(blood|biomarker)"

manifest_candidates <- function() {
  z <- c(
    CHARLS_MANIFEST,
    file.path(CHARLS_CLEAN, "_extract_qc", "extracted_manifest.tsv"),
    file.path(CHARLS_CLEAN, "extracted_manifest.tsv"),
    file.path(dirname(CHARLS_CLEAN), "_extract_qc", "extracted_manifest.tsv")
  )
  unique(z[nzchar(z)])
}

load_manifest <- function() {
  if (!is.null(MANIFEST_DF)) return(MANIFEST_DF)
  fs <- manifest_candidates()
  fs <- fs[file.exists(fs)]
  if (!length(fs)) {
    msg("Manifest not found; use built-in canonical Blood/Biomarker file list and filesystem fallback.")
    MANIFEST_DF <<- data.table()
    return(MANIFEST_DF)
  }
  f <- fs[1]
  x <- tryCatch(fread(f, sep = "\t", encoding = "UTF-8"), error = function(e) {
    msg("WARN cannot read CHARLS manifest: ", f, "; ", conditionMessage(e))
    data.table()
  })
  need <- c("wave", "archive_rel", "internal_file", "output_file", "action", "size_bytes")
  if (!all(need %in% names(x))) {
    msg("WARN manifest lacks required columns; ignored: ", f)
    x <- data.table()
  }
  if (nrow(x)) {
    x[, output_file := gsub("\\\\", "/", output_file)]
    x[, wave := as.character(wave)]
    msg("Manifest loaded: ", f, "; rows=", nrow(x), "; renamed_conflict=", sum(x$action == "renamed_conflict", na.rm = TRUE), "; skipped_identical=", sum(x$action == "skipped_identical", na.rm = TRUE))
  }
  MANIFEST_DF <<- x
  MANIFEST_DF
}

add_file_choice <- function(wave, role, path, archive_rel = NA_character_, internal_file = NA_character_, action = NA_character_, size_bytes = NA_real_, note = "") {
  FILE_CHOICE <<- rbind(FILE_CHOICE, data.table(
    wave = as.character(wave), file_role = role, path = path,
    archive_rel = archive_rel, internal_file = internal_file,
    action = action, size_bytes = suppressWarnings(as.numeric(size_bytes)), note = note
  ), fill = TRUE)
}

# Canonical raw files inferred from the uploaded extraction manifest. These are
# the only raw files needed for blood biomarkers. None of the three manifest
# renamed_conflict files is needed by this analysis.
canonical_blood_biomarker_rel <- list(
  `1` = c("2011/Blood_20140429.dta", "2011/biomarkers.dta"),
  `2` = c("2013/Biomarker.dta"),
  `3` = c("2015/Biomarker.dta", "2015/Blood.dta"),
  `4` = character(),
  `5` = character()
)

manifest_wave_files <- function(w, role_regex = raw_file_regex) {
  m <- load_manifest()
  if (!nrow(m)) return(character())
  yrs <- wave_raw_years[[as.character(w)]] %||% character()
  z <- m[wave %in% yrs]
  if (!nrow(z)) return(character())
  z <- z[grepl("\\.dta$", output_file, ignore.case = TRUE) & grepl(role_regex, basename(output_file), perl = TRUE)]
  if (!nrow(z)) return(character())
  # Keep canonical files copied to the clean folder. skipped_identical rows point
  # to an already copied output_file; renamed_conflict rows are deliberately excluded.
  z <- z[action == "copied"]
  if (!nrow(z)) return(character())
  z <- unique(z, by = "output_file")
  z[, abs_path := file.path(CHARLS_CLEAN, output_file)]
  z <- z[file.exists(abs_path)]
  if (nrow(z)) {
    for (i in seq_len(nrow(z))) add_file_choice(w, "manifest_copied_blood_biomarker", z$abs_path[i], z$archive_rel[i], z$internal_file[i], z$action[i], z$size_bytes[i], "selected from manifest action=copied")
  }
  z$abs_path
}

builtin_wave_files <- function(w) {
  rel <- canonical_blood_biomarker_rel[[as.character(w)]] %||% character()
  fs <- file.path(CHARLS_CLEAN, rel)
  fs <- fs[file.exists(fs)]
  if (length(fs)) for (f in fs) add_file_choice(w, "builtin_canonical_blood_biomarker", f, note = "selected from built-in manifest-derived canonical list")
  fs
}

raw_wave_file_paths <- function(w) {
  fs <- manifest_wave_files(w)
  if (!length(fs)) fs <- builtin_wave_files(w)
  if (!length(fs)) {
    yrs <- wave_raw_years[[as.character(w)]] %||% character()
    fs <- unlist(lapply(yrs, function(y) {
      d <- file.path(CHARLS_CLEAN, y)
      if (!dir.exists(d)) return(character())
      list.files(d, pattern = "\\.dta$", full.names = TRUE, ignore.case = TRUE)
    }))
    fs <- fs[grepl(raw_file_regex, basename(fs), perl = TRUE)]
    fs <- fs[!grepl("__", basename(fs), fixed = TRUE)]
    if (length(fs)) for (f in fs) add_file_choice(w, "filesystem_fallback_blood_biomarker", f, note = "manifest absent; selected by Blood/Biomarker filename; files with __ are excluded")
  }
  unique(fs)
}

write_manifest_qc <- function() {
  m <- load_manifest()
  if (nrow(m)) {
    conflicts <- m[action == "renamed_conflict", .(wave, archive_rel, internal_file, output_file, action, size_bytes)]
    selected <- copy(FILE_CHOICE)
    write_xlsx(list(
      selected_files = selected,
      renamed_conflicts = conflicts,
      manifest_action_counts = m[, .N, by = action][order(action)]
    ), file.path(DIR_QC, "manifest_file_choice.xlsx"))
  } else {
    write_xlsx(FILE_CHOICE, file.path(DIR_QC, "manifest_file_choice.xlsx"))
  }
}


norm_id_exact <- function(x) {
  # Exact ID normalization only for safe comparison/join:
  # - Stata-labelled values are unlabelled.
  # - Numeric IDs are printed without scientific notation.
  # - Leading/trailing spaces and a terminal ".0" are removed.
  # No zero-padding, no leading-zero removal, no household+person concatenation,
  # and no "choose the largest overlap" heuristic are used.
  if (inherits(x, "haven_labelled")) x <- haven::zap_labels(x)
  if (is.numeric(x)) {
    y <- ifelse(is.na(x), NA_character_, sprintf("%.0f", x))
  } else {
    y <- as.character(x)
  }
  y <- trimws(y)
  y <- gsub("\\.0$", "", y)
  y[y %in% c("", "NA", "NaN", ".")] <- NA_character_
  y
}

explicit_raw_id_plan <- function() {
  # This plan is based on the user's wave-1 ID QC:
  #   raw_2011_blood       ID  matches Harmonized ID_w1 exactly: 11,847/11,847 raw IDs.
  #   raw_2011_biomarkers  ID  matches Harmonized ID_w1 exactly: 13,974/13,974 raw IDs.
  #   raw_2012_Biomarkers  ID  does NOT match Harmonized ID_w1/ID and is not used.
  #
  # Wave 2/3 raw biomarker files in the current clean folder use the stable
  # Harmonized person ID, so they are joined by H$ID <-> raw$ID. If a local
  # file fails this exact check, the pipeline stops rather than trying to guess.
  data.table(
    wave = c(1L, 1L, 1L, 2L, 3L, 3L),
    rel_path = c(
      "2011/Blood_20140429.dta",
      "2011/biomarkers.dta",
      "2012/Biomarkers.dta",
      "2013/Biomarker.dta",
      "2015/Biomarker.dta",
      "2015/Blood.dta"
    ),
    h_key = c("ID_w1", "ID_w1", NA_character_, "ID", "ID", "ID"),
    raw_key = c("ID", "ID", NA_character_, "ID", "ID", "ID"),
    use_file = c(TRUE, TRUE, FALSE, TRUE, TRUE, TRUE),
    min_exact_overlap = c(10000L, 10000L, NA_integer_, 1000L, 1000L, 1000L),
    note = c(
      "Use exact join: Harmonized ID_w1 <-> raw ID. QC exact_overlap_n=11847 and exact_overlap_pct_of_raw=1.",
      "Use exact join: Harmonized ID_w1 <-> raw ID. QC exact_overlap_n=13974 and exact_overlap_pct_of_raw=1.",
      "Do not use. QC showed exact overlap with Harmonized IDs is essentially zero; this is not a safe source for wave-1 blood lipids.",
      "Use exact join: Harmonized stable ID <-> raw ID.",
      "Use exact join: Harmonized stable ID <-> raw ID.",
      "Use exact join: Harmonized stable ID <-> raw ID."
    )
  )
}

raw_plan_abs <- function() {
  plan <- explicit_raw_id_plan()
  plan[, path := file.path(CHARLS_CLEAN, rel_path)]
  plan[, exists := file.exists(path)]
  plan[]
}

file_plan_to_qc <- function(plan) {
  for (i in seq_len(nrow(plan))) {
    add_file_choice(
      wave = plan$wave[i],
      role = ifelse(isTRUE(plan$use_file[i]), "explicit_id_plan_selected", "explicit_id_plan_skipped"),
      path = plan$path[i],
      note = paste0("h_key=", plan$h_key[i], "; raw_key=", plan$raw_key[i], "; ", plan$note[i])
    )
  }
}

merge_one_raw_explicit <- function(h, plan_row) {
  f <- plan_row$path
  if (!file.exists(f)) {
    msg("RAW explicit wave=", plan_row$wave, ": file not found and skipped: ", f)
    ID_MATCH_QC <<- rbind(ID_MATCH_QC, data.table(
      wave = as.character(plan_row$wave), file = plan_row$rel_path,
      id_source = plan_row$raw_key, id_transform = paste0("exact_to_", plan_row$h_key),
      overlap = 0L, raw_n = 0L, harmonized_n = 0L,
      note = "file_not_found"
    ), fill = TRUE)
    return(h)
  }

  if (!plan_row$h_key %in% names(h)) {
    stop("Explicit raw ID plan refers to missing Harmonized key: ", plan_row$h_key,
         " for file ", plan_row$rel_path, call. = FALSE)
  }

  x <- read_dta_dt(f)
  if (!plan_row$raw_key %in% names(x)) {
    stop("Explicit raw ID plan refers to missing raw key: ", plan_row$raw_key,
         " in file ", plan_row$rel_path, call. = FALSE)
  }

  hkey <- norm_id_exact(h[[plan_row$h_key]])
  rkey <- norm_id_exact(x[[plan_row$raw_key]])
  h_unique <- unique(hkey[!is.na(hkey) & nzchar(hkey)])
  r_unique <- unique(rkey[!is.na(rkey) & nzchar(rkey)])
  overlap <- length(intersect(h_unique, r_unique))
  raw_n <- length(r_unique)
  h_n <- length(h_unique)

  ID_MATCH_QC <<- rbind(ID_MATCH_QC, data.table(
    wave = as.character(plan_row$wave),
    file = plan_row$rel_path,
    id_source = plan_row$raw_key,
    id_transform = paste0("exact_to_H_", plan_row$h_key),
    overlap = overlap,
    raw_n = raw_n,
    harmonized_n = h_n,
    note = plan_row$note
  ), fill = TRUE)

  msg("RAW explicit wave=", plan_row$wave, ": ", plan_row$rel_path,
      "; join H$", plan_row$h_key, " <-> raw$", plan_row$raw_key,
      "; exact_overlap=", overlap, "/", raw_n, " raw IDs")

  minov <- suppressWarnings(as.integer(plan_row$min_exact_overlap))
  if (is.finite(minov) && overlap < minov) {
    suppressWarnings(write_xlsx(ID_MATCH_QC, file.path(DIR_QC, "raw_id_match_qc.failed.xlsx")))
    stop("Unsafe raw-Harmonized ID join for ", plan_row$rel_path,
         ": exact overlap is ", overlap, ", below min_exact_overlap=", minov,
         ". No automatic ID conversion will be attempted. Check qc/raw_id_match_qc.failed.xlsx.",
         call. = FALSE)
  }

  x[, .merge_key_raw := rkey]
  x <- x[!is.na(.merge_key_raw) & nzchar(.merge_key_raw)]
  dup <- x[, .N, by = .merge_key_raw][N > 1]
  if (nrow(dup)) {
    write_xlsx(dup, file.path(DIR_QC, paste0("duplicate_raw_ids_", gsub("[^A-Za-z0-9]+", "_", plan_row$rel_path), ".xlsx")))
    msg("WARN raw file has duplicate IDs; first row kept after writing duplicate QC: ", plan_row$rel_path)
    x <- x[, .SD[1], by = .merge_key_raw]
  }

  yr <- basename(dirname(f))
  base <- tolower(clean_names(tools::file_path_sans_ext(basename(f))))
  keep <- setdiff(names(x), c(".merge_key_raw", plan_row$raw_key))
  keep <- keep[!grepl("^(id|hhid|householdid|communityid|version|release)$", keep, ignore.case = TRUE)]
  if (!length(keep)) return(h)
  x <- x[, c(".merge_key_raw", keep), with = FALSE]
  setnames(x, keep, paste0("raw", yr, "_", base, "_", tolower(clean_names(keep))))

  h[, .rowid_h := .I]
  h[, .merge_key_raw := hkey]
  before <- ncol(h)
  out <- merge(h, x, by = ".merge_key_raw", all.x = TRUE, sort = FALSE)
  setorder(out, .rowid_h)
  out[, c(".rowid_h", ".merge_key_raw") := NULL]
  msg("RAW explicit wave=", plan_row$wave, ": merged ", ncol(out) - before + 2L,
      " columns from ", plan_row$rel_path)
  out[]
}

merge_raw_fallback <- function(h) {
  # Strict explicit-ID raw merge. This fixes the previous wave-1 mismatch:
  # v003/v004 used the stable Harmonized ID for all raw files, but wave-1
  # 2011 Blood/Biomarker files use the 11-character wave-specific ID.
  load_manifest()
  plan <- raw_plan_abs()
  file_plan_to_qc(plan)
  write_xlsx(plan, file.path(DIR_QC, "explicit_raw_id_plan.xlsx"))

  skipped <- plan[use_file == FALSE | exists == FALSE]
  if (nrow(skipped)) {
    for (i in seq_len(nrow(skipped))) {
      msg("RAW explicit plan skip: ", skipped$rel_path[i], "; ", skipped$note[i])
    }
  }

  use_plan <- plan[use_file == TRUE & exists == TRUE]
  out <- h
  for (i in seq_len(nrow(use_plan))) {
    out <- merge_one_raw_explicit(out, use_plan[i])
  }

  write_manifest_qc()
  suppressWarnings(write_xlsx(ID_MATCH_QC, file.path(DIR_QC, "raw_id_match_qc.xlsx")))
  suppressWarnings(write_xlsx(labels_dt(out), file.path(DIR_QC, "available_variables_after_raw_merge.xlsx")))
  out
}

first_event <- function(waves_dt, baseline_wave = 1, follow_waves = c(2, 3, 4, 5)) {
  base <- waves_dt[wave == baseline_wave]
  fu <- waves_dt[wave %in% follow_waves, .(id, wave, year, CMM, diabetes, heart_disease, stroke)]
  fu <- fu[!is.na(CMM)]
  ev <- fu[CMM == 1, .SD[which.min(year)], by = id]
  last <- fu[, .SD[which.max(year)], by = id]
  ans <- base[, .(id)]
  ans <- merge(ans, ev[, .(id, event_year = year, event_wave = wave)], by = "id", all.x = TRUE)
  ans <- merge(ans, last[, .(id, last_year = year, last_wave = wave)], by = "id", all.x = TRUE)
  byear <- unname(wv_year[as.character(baseline_wave)])
  ans[, event := as.numeric(!is.na(event_year))]
  ans[, time := fifelse(event == 1, event_year - byear, last_year - byear)]
  ans[is.na(time) | time <= 0, `:=`(event = NA_real_, time = NA_real_)]
  ans
}

simple_impute <- function(dt, vars) {
  x <- copy(dt)
  for (v in vars) {
    if (!v %in% names(x)) next
    if (is.factor(x[[v]]) || is.character(x[[v]])) {
      m <- mode_value(as.character(x[[v]])); x[is.na(get(v)), (v) := m]
    } else {
      if (length(unique(na.omit(x[[v]]))) <= 2) {
        m <- suppressWarnings(as.numeric(mode_value(x[[v]])))
      } else m <- suppressWarnings(median(x[[v]], na.rm = TRUE))
      x[is.na(get(v)), (v) := m]
    }
  }
  x
}

impute_baseline <- function(dt) {
  vars <- c("SBP", "DBP", "FBG", "HbA1c", "BMI", "smoking", "dyslipidemia", "drinking", "TC", "LDLC")
  miss <- data.table(variable = vars, missing_prop = vapply(vars, function(v) if (v %in% names(dt)) mean(is.na(dt[[v]])) else NA_real_, numeric(1)))
  write_xlsx(miss, file.path(DIR_QC, "missingness_for_imputation.xlsx"))
  bad <- miss[!is.na(missing_prop) & missing_prop > 0.20]
  if (nrow(bad) && !ALLOW_HIGH_MISSING) {
    stop("Some key variables have >20% missingness. Check qc/var_map.xlsx and extracted variables, or rerun with ALLOW_HIGH_MISSING=1. Variables: ", paste(bad$variable, collapse = ", "), call. = FALSE)
  }
  if (USE_MICE) {
    if (!opt_pkg("mice")) {
      stop("USE_MICE=1 but R package 'mice' is not available in this R library path. Install it in the same R used by charls.sh, or run with USE_MICE=0.
Check with: Rscript -e 'library(mice); sessionInfo()'", call. = FALSE)
    }
    msg("Imputing selected covariates using mice: PMM for quantitative variables and logistic regression for binary variables")
    suppressPackageStartupMessages(library(mice))
    imp_vars <- unique(c("event", "time", "AIPFI", vars, "age", "gender", "hypertension", "diabetes", "heart_disease", "stroke", "TG", "HDLC"))
    tmp <- as.data.frame(dt[, ..imp_vars])
    meth <- mice::make.method(tmp)
    pred <- mice::make.predictorMatrix(tmp)
    meth[] <- ""
    for (v in vars) {
      if (!v %in% names(tmp)) next
      uv <- unique(na.omit(tmp[[v]]))
      if (length(uv) <= 2) meth[v] <- "logreg" else meth[v] <- "pmm"
    }
    suppressWarnings({
      imp <- mice::mice(tmp, m = 1, maxit = 5, method = meth, predictorMatrix = pred, seed = SEED, printFlag = FALSE)
      comp <- data.table(mice::complete(imp, 1))
    })
    out <- copy(dt)
    for (v in vars) if (v %in% names(comp)) out[[v]] <- comp[[v]]
    out
  } else {
    msg("USE_MICE=0; using median/mode single imputation fallback")
    simple_impute(dt, vars)
  }
}

make_baseline_dataset <- function(waves_dt) {
  base <- waves_dt[wave == 1]
  surv <- first_event(waves_dt, baseline_wave = 1, follow_waves = c(2, 3, 4, 5))
  dat <- merge(base, surv, by = "id", all.x = TRUE)
  dat[, exclude_prevalent_CMM := CMM == 1]
  dat[, exclude_missing_baseline_AIPFI := is.na(AIPFI) | is.na(AIP) | is.na(FI)]
  dat[, exclude_missing_baseline_CMM := is.na(CMM)]
  dat[, exclude_incomplete_followup := is.na(time) | is.na(event)]
  dat <- dat[is.na(exclude_prevalent_CMM) | exclude_prevalent_CMM == FALSE]
  dat <- dat[exclude_missing_baseline_AIPFI == FALSE & exclude_missing_baseline_CMM == FALSE & exclude_incomplete_followup == FALSE]
  dat[, AIPFI_Q := make_quartile(AIPFI)]
  dat[, AIPFI_sd := as.numeric(scale(AIPFI))]
  dat[, age_group := factor(ifelse(age < 60, "<60", ">=60"), levels = c("<60", ">=60"))]
  dat <- impute_baseline(dat)
  dat[, gender := factor(gender, levels = c("Female", "Male"))]
  dat[]
}

make_longitudinal_dataset <- function(waves_dt, baseline_dat) {
  w1 <- waves_dt[wave == 1, .(id, AIPFI_2011 = AIPFI, AIP_2011 = AIP, FI_2011 = FI, CMM_2011 = CMM)]
  w3 <- waves_dt[wave == 3, .(id, AIPFI_2015 = AIPFI, AIP_2015 = AIP, FI_2015 = FI, CMM_2015 = CMM)]
  basecov <- baseline_dat[, .(id, age, gender, smoking, drinking, BMI, SBP, DBP, FBG, HbA1c, TC, TG, HDLC, LDLC, hypertension, diabetes, dyslipidemia, heart_disease, stroke, AIPFI, AIPFI_Q, AIPFI_sd, AIP, FI)]
  dat <- merge(w1, w3, by = "id")
  dat <- merge(dat, basecov, by = "id")
  surv <- first_event(waves_dt, baseline_wave = 3, follow_waves = c(4, 5))
  dat <- merge(dat, surv[, .(id, event, time, event_year, last_year)], by = "id", all.x = TRUE)
  # Dynamic cohort: complete AIPFI at 2011 and 2015, no CMM during that period, complete follow-up after 2015.
  dat <- dat[!is.na(AIPFI_2011) & !is.na(AIPFI_2015) & !is.na(CMM_2011) & !is.na(CMM_2015)]
  dat <- dat[CMM_2011 == 0 & CMM_2015 == 0]
  dat <- dat[!is.na(event) & !is.na(time)]
  dat[, cum_AIPFI := (AIPFI_2011 + AIPFI_2015) / 2 * (2015 - 2012)]
  dat[, cum_AIPFI_Q := make_quartile(cum_AIPFI)]
  dat[, cum_AIPFI_sd := as.numeric(scale(cum_AIPFI))]
  dat[, delta_AIPFI := AIPFI_2015 - AIPFI_2011]
  cl <- cluster_aipfi(dat, k = 2)
  dat <- merge(dat, cl$assign[, .(id, cluster2)], by = "id", all.x = TRUE)
  dat[, cluster2 := factor(cluster2, levels = c("Cluster 1", "Cluster 2"))]
  dat[]
}

cluster_aipfi <- function(dat, k = 2) {
  x <- as.matrix(dat[, .(AIPFI_2011, AIPFI_2015)])
  ok <- complete.cases(x)
  xs <- scale(x[ok, , drop = FALSE])
  set.seed(SEED)
  km <- kmeans(xs, centers = k, nstart = 50)
  centers_orig <- data.table::as.data.table(
    aggregate(x[ok, , drop = FALSE], list(cluster = km$cluster), mean)
  )
  centers_orig[, mean_level := rowMeans(.SD), .SDcols = c("AIPFI_2011", "AIPFI_2015")]
  ord <- order(centers_orig$mean_level)
  lab <- paste0("Cluster ", match(km$cluster, ord))
  assign <- data.table(id = dat$id[ok], cluster_raw = km$cluster, cluster = lab)
  if (k == 2) setnames(assign, "cluster", "cluster2") else setnames(assign, "cluster", paste0("cluster", k))
  list(km = km, assign = assign, centers = centers_orig)
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 🚩 Prep: read CHARLS, construct analytic datasets
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
run_prep <- function() {
  msg("🚩 prep_dat: build wave-level, baseline, and longitudinal datasets")
  msg("Code version: ", CHARLS_VERSION)
  h <- load_harmonized()
  h <- merge_raw_fallback(h)
  waves <- rbindlist(lapply(1:5, function(w) build_wave(h, w)), fill = TRUE)
  qsave(waves, file.path(DIR_DAT, "charls_waves_2011_2020.rds"))
  baseline <- make_baseline_dataset(waves)
  longitudinal <- make_longitudinal_dataset(waves, baseline)
  qsave(baseline, file.path(DIR_DAT, "baseline_2011_2020.rds"))
  qsave(longitudinal, file.path(DIR_DAT, "longitudinal_2015_2020.rds"))
  write_xlsx(VAR_MAP, file.path(DIR_QC, "var_map.xlsx"))
  write_xlsx(FI_MAP, file.path(DIR_QC, "fi_item_map.xlsx"))
  flow <- make_flow_table(waves, baseline, longitudinal)
  write_xlsx(flow, file.path(DIR_QC, "cohort_counts.xlsx"))
  msg("Baseline analytic n=", nrow(baseline), "; events=", sum(baseline$event == 1, na.rm = TRUE))
  msg("Longitudinal analytic n=", nrow(longitudinal), "; events=", sum(longitudinal$event == 1, na.rm = TRUE))
}

make_flow_table <- function(waves, baseline, longitudinal) {
  w1 <- waves[wave == 1]
  data.table(
    step = c("CHARLS baseline wave 2011/2012", "Eligible baseline analytic sample", "Longitudinal AIPFI sample"),
    n = c(uniqueN(w1$id), nrow(baseline), nrow(longitudinal)),
    events = c(NA_integer_, sum(baseline$event == 1, na.rm = TRUE), sum(longitudinal$event == 1, na.rm = TRUE)),
    note = c("All participants detected in Harmonized CHARLS baseline", "No prevalent CMM; nonmissing baseline AIPFI/CMM; follow-up CMM ascertainment", "Complete AIPFI in 2011 and 2015; no CMM through 2015; follow-up 2018/2020")
  )
}

load_dat <- function(which = c("baseline", "longitudinal", "waves")) {
  which <- match.arg(which)
  if (which == "baseline") qread(file.path(DIR_DAT, "baseline_2011_2020.rds"))
  else if (which == "longitudinal") qread(file.path(DIR_DAT, "longitudinal_2015_2020.rds"))
  else qread(file.path(DIR_DAT, "charls_waves_2011_2020.rds"))
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 🚩 Tables: descriptive and Cox models
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
is_binary_var <- function(x) {
  u <- sort(unique(na.omit(x)))
  length(u) <= 2 && all(u %in% c(0, 1))
}

desc_one <- function(x) {
  if (is.factor(x) || is.character(x) || is_binary_var(x)) {
    if (is.factor(x) || is.character(x)) {
      tab <- table(x, useNA = "no")
      paste(sprintf("%s: %d (%.2f%%)", names(tab), as.integer(tab), 100 * as.integer(tab) / sum(tab)), collapse = "; ")
    } else {
      n <- sum(x == 1, na.rm = TRUE); d <- sum(!is.na(x)); sprintf("%d (%.2f%%)", n, 100 * n / d)
    }
  } else {
    sprintf("%.2f (%.2f, %.2f)", median(x, na.rm = TRUE), quantile(x, 0.25, na.rm = TRUE), quantile(x, 0.75, na.rm = TRUE))
  }
}

p_group <- function(x, g) {
  ok <- !is.na(g) & !is.na(x)
  if (sum(ok) < 3 || length(unique(g[ok])) < 2) return(NA_real_)
  if (is.factor(x) || is.character(x) || is_binary_var(x)) {
    suppressWarnings(chisq.test(table(x[ok], g[ok]))$p.value)
  } else {
    suppressWarnings(kruskal.test(x[ok] ~ g[ok])$p.value)
  }
}

make_desc_table <- function(dt, group = NULL, vars = NULL, labels = NULL) {
  if (is.null(vars)) vars <- c("age", "gender", "diabetes", "hypertension", "heart_disease", "stroke", "dyslipidemia", "drinking", "smoking", "SBP", "DBP", "BMI", "FBG", "HbA1c", "TC", "TG", "HDLC", "LDLC", "AIP", "FI", "AIPFI")
  vars <- vars[vars %in% names(dt)]
  if (is.null(labels)) labels <- vars
  names(labels) <- vars
  if (is.null(group)) {
    data.table(variable = labels[vars], Total = vapply(vars, function(v) desc_one(dt[[v]]), character(1)))
  } else {
    g <- dt[[group]]
    lev <- levels(as.factor(g))
    out <- data.table(variable = labels[vars], Total = vapply(vars, function(v) desc_one(dt[[v]]), character(1)))
    for (lv in lev) out[[lv]] <- vapply(vars, function(v) desc_one(dt[g == lv][[v]]), character(1))
    out <- data.table::copy(out)
    out[, P := fmt_p(vapply(vars, function(v) p_group(dt[[v]], g), numeric(1)))]
    out
  }
}

covars_m1 <- character()
covars_m2 <- c("age", "gender")
covars_m3 <- c("age", "gender", "hypertension", "diabetes", "dyslipidemia", "heart_disease", "stroke", "smoking", "drinking")
covars_m4 <- c("age", "gender", "hypertension", "diabetes", "dyslipidemia", "heart_disease", "stroke", "smoking", "drinking", "BMI", "DBP", "SBP", "FBG", "HbA1c", "TC", "TG", "HDLC", "LDLC")
model_covars <- list(Model1 = covars_m1, Model2 = covars_m2, Model3 = covars_m3, Model4 = covars_m4)

cox_fit <- function(dt, exposure, covars = character(), weights = NULL, strata_vars = NULL) {
  cv <- covars[covars %in% names(dt)]
  if (!is.null(strata_vars)) {
    cv <- setdiff(cv, strata_vars)
    str <- paste0("strata(", strata_vars[strata_vars %in% names(dt)], ")")
    cv <- c(cv, str)
  }
  f <- as.formula(paste("Surv(time, event) ~", paste(c(exposure, cv), collapse = " + ")))
  if (is.null(weights)) survival::coxph(f, data = dt, x = TRUE, model = TRUE) else survival::coxph(f, data = dt, weights = weights, robust = TRUE, x = TRUE, model = TRUE)
}

cox_extract <- function(fit, exposure_label = NULL) {
  sm <- summary(fit)

  cf <- data.table(term = rownames(sm$coefficients), as.data.frame(sm$coefficients, check.names = FALSE))
  ci <- data.table(term = rownames(sm$conf.int), as.data.frame(sm$conf.int, check.names = FALSE))

  # Robust across survival/data.table versions:
  # summary.coxph() always stores conf.int columns in this order:
  #   exp(coef), exp(-coef), lower .95, upper .95
  # but their names may appear as exp(coef), exp.coef, exp.coef., etc.
  # Therefore do not depend on exact column names.
  ci_cols <- setdiff(names(ci), "term")
  if (length(ci_cols) < 4L) {
    stop("cox_extract failed: summary(fit)$conf.int has fewer than 4 non-term columns. Available columns: ",
         paste(names(ci), collapse = ", "), call. = FALSE)
  }
  ci2 <- ci[, .(
    term = term,
    HR  = as.numeric(get(ci_cols[1L])),
    LCL = as.numeric(get(ci_cols[3L])),
    UCL = as.numeric(get(ci_cols[4L]))
  )]

  # P-value column is usually Pr(>|z|), but robust Cox summaries can vary.
  # Prefer a column whose name starts with Pr; otherwise use the last numeric
  # coefficient-summary column as a fallback and record the available columns in errors.
  pcol <- grep("^Pr", names(cf), value = TRUE)[1]
  if (is.na(pcol) || !nzchar(pcol)) {
    cf_cols <- setdiff(names(cf), "term")
    num_cols <- cf_cols[vapply(cf[, ..cf_cols], is.numeric, logical(1))]
    pcol <- tail(num_cols, 1L)
  }
  if (is.na(pcol) || !nzchar(pcol)) {
    stop("cox_extract failed: P-value column not found. Available coefficient columns: ",
         paste(names(cf), collapse = ", "), call. = FALSE)
  }
  cf2 <- cf[, .(term = term, p = as.numeric(get(pcol)))]

  z <- merge(cf2, ci2, by = "term", all.x = TRUE)
  z[, HR_CI := fmt_hr(HR, LCL, UCL)]
  z[, P := fmt_p(p)]
  if (!is.null(exposure_label)) z[, exposure := exposure_label]
  z[]
}

run_cox_models <- function(dt, exposure, label = exposure, covar_list = model_covars) {
  res <- rbindlist(lapply(names(covar_list), function(m) {
    fit <- cox_fit(dt, exposure, covar_list[[m]])
    ex <- cox_extract(fit, label)
    ex[, model := m]
    ex
  }), fill = TRUE)
  res[]
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 🚩 Fig1: flow chart
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
run_fig1 <- function() {
  msg("🚩 Fig1: participant inclusion/exclusion flow diagram")
  waves <- load_dat("waves"); baseline <- load_dat("baseline"); long <- load_dat("longitudinal")
  flow <- make_flow_table(waves, baseline, long)
  write_xlsx(flow, file.path(DIR_TAB, "Fig1.out.xlsx"))
  png(file.path(DIR_FIG, "Fig1.png"), width = 1800, height = 1200, res = 180)
  par(mar = c(1, 1, 1, 1))
  plot.new(); plot.window(xlim = c(0, 1), ylim = c(0, 1))
  rect_box <- function(x, y, text, w = 0.72, h = 0.14) {
    rect(x - w/2, y - h/2, x + w/2, y + h/2, border = "black", lwd = 2)
    text(x, y, text, cex = 0.9)
  }
  rect_box(0.5, 0.82, sprintf("CHARLS baseline 2011/2012\nN = %s", format(flow$n[1], big.mark = ",")))
  arrows(0.5, 0.74, 0.5, 0.62, length = 0.08, lwd = 2)
  rect_box(0.5, 0.55, sprintf("Baseline analytic sample\nNo prevalent CMM; complete AIPFI/CMM/follow-up\nN = %s; incident CMM = %s", format(nrow(baseline), big.mark = ","), format(sum(baseline$event == 1), big.mark = ",")), h = 0.18)
  arrows(0.5, 0.45, 0.5, 0.32, length = 0.08, lwd = 2)
  rect_box(0.5, 0.24, sprintf("Longitudinal AIPFI sample\nComplete AIPFI in 2011 and 2015; CMM-free through 2015\nN = %s; incident CMM = %s", format(nrow(long), big.mark = ","), format(sum(long$event == 1), big.mark = ",")), h = 0.18)
  title("Fig. 1 Participant inclusion and exclusion criteria", cex.main = 1.1)
  dev.off()
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 🚩 Tab1: baseline characteristics by AIPFI quartile
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
run_tab1 <- function() {
  msg("🚩 Tab1: baseline characteristics by AIPFI quartiles")
  dt <- load_dat("baseline")
  tab <- make_desc_table(dt, group = "AIPFI_Q")
  write_xlsx(tab, file.path(DIR_TAB, "Tab1.out.xlsx"))
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 🚩 Tab2: baseline AIPFI Cox models
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
run_tab2 <- function() {
  msg("🚩 Tab2: Cox models for baseline AIPFI")
  dt <- load_dat("baseline")
  dt[, AIPFI_Q := relevel(AIPFI_Q, "Q1")]
  res1 <- run_cox_models(dt, "AIPFI", "AIPFI per 1 unit")
  res2 <- run_cox_models(dt, "AIPFI_sd", "AIPFI per 1 SD")
  res3 <- run_cox_models(dt, "AIPFI_Q", "AIPFI quartile")
  res <- rbindlist(list(res1, res2, res3), fill = TRUE)
  keep <- grepl("AIPFI", res$term)
  res <- res[keep, .(exposure, model, term, HR, LCL, UCL, HR_CI, p, P)]
  write_xlsx(res, file.path(DIR_TAB, "Tab2.out.xlsx"))
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 🚩 Tab3: longitudinal/cumulative AIPFI Cox models
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
run_tab3 <- function() {
  msg("🚩 Tab3: Cox models for cumulative AIPFI and trajectory clusters")
  dt <- load_dat("longitudinal")
  dt[, cum_AIPFI_Q := relevel(cum_AIPFI_Q, "Q1")]
  dt[, cluster2 := relevel(cluster2, "Cluster 1")]
  res1 <- run_cox_models(dt, "cum_AIPFI", "Cumulative AIPFI per 1 unit")
  res2 <- run_cox_models(dt, "cum_AIPFI_sd", "Cumulative AIPFI per 1 SD")
  res3 <- run_cox_models(dt, "cum_AIPFI_Q", "Cumulative AIPFI quartile")
  res4 <- run_cox_models(dt, "cluster2", "AIPFI trajectory cluster")
  res <- rbindlist(list(res1, res2, res3, res4), fill = TRUE)
  keep <- grepl("cum_AIPFI|cluster2", res$term)
  res <- res[keep, .(exposure, model, term, HR, LCL, UCL, HR_CI, p, P)]
  write_xlsx(res, file.path(DIR_TAB, "Tab3.out.xlsx"))
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 🚩 Fig2: Kaplan-Meier curves
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
km_data <- function(dt, group) {
  f <- as.formula(paste0("Surv(time, event) ~ ", group))
  fit <- survival::survfit(f, data = dt)
  s <- summary(fit)
  z <- data.table(time = s$time, surv = s$surv, lower = s$lower, upper = s$upper, strata = as.character(s$strata))
  z[, group := sub(paste0("^", group, "="), "", strata)]
  z[, cuminc := 1 - surv]
  z[, lower_ci := 1 - upper]
  z[, upper_ci := 1 - lower]
  z
}

plot_km <- function(dt, group, title) {
  kd <- km_data(dt, group)
  ggplot(kd, aes(x = time, y = cuminc, color = group)) +
    geom_step(linewidth = 0.9) +
    labs(x = "Follow-up time, years", y = "Cumulative incidence of CMM", color = NULL, title = title) +
    theme_bw(base_size = 12) + theme(legend.position = "bottom")
}

save_plot <- function(p, f, w = 7.2, h = 5.2) {
  ggplot2::ggsave(f, p, width = w, height = h, dpi = 300)
  invisible(f)
}

combine_gg <- function(plots, f, ncol = 2, w = 12, h = 6) {
  if (opt_pkg("gridExtra")) {
    png(f, width = w, height = h, units = "in", res = 300)
    gridExtra::grid.arrange(grobs = plots, ncol = ncol)
    dev.off()
  } else {
    png(f, width = w, height = h, units = "in", res = 300)
    grid::grid.newpage(); grid::grid.text("Install package 'gridExtra' to render combined panel figure.", 0.5, 0.5)
    dev.off()
  }
}

run_fig2 <- function() {
  msg("🚩 Fig2: KM curves for baseline AIPFI quartiles and trajectory clusters")
  b <- load_dat("baseline"); l <- load_dat("longitudinal")
  p1 <- plot_km(b, "AIPFI_Q", "A. Baseline AIPFI quartiles")
  p2 <- plot_km(l, "cluster2", "B. AIPFI trajectory clusters")
  combine_gg(list(p1, p2), file.path(DIR_FIG, "Fig2.png"), ncol = 2, w = 12, h = 5.5)
  out <- list(
    AIPFI_quartile_survival = km_data(b, "AIPFI_Q"),
    cluster_survival = km_data(l, "cluster2"),
    logrank = data.table(analysis = c("AIPFI_Q", "cluster2"), p = c(survdiff_p(b, "AIPFI_Q"), survdiff_p(l, "cluster2")))
  )
  write_xlsx(out, file.path(DIR_TAB, "Fig2.out.xlsx"))
}

survdiff_p <- function(dt, group) {
  f <- as.formula(paste0("Surv(time, event) ~ ", group))
  sd <- survival::survdiff(f, data = dt)
  stats::pchisq(sd$chisq, df = length(sd$n) - 1, lower.tail = FALSE)
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 🚩 Fig3: K-means AIPFI trajectories
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
run_fig3 <- function() {
  msg("🚩 Fig3: clustering of AIPFI change from 2011/2012 to 2015")
  dt <- load_dat("longitudinal")
  x <- scale(as.matrix(dt[, .(AIPFI_2011, AIPFI_2015)]))
  wss <- data.table(k = 1:6, wss = sapply(1:6, function(k) kmeans(x, centers = k, nstart = 30)$tot.withinss))
  p1 <- ggplot(wss, aes(k, wss)) + geom_line() + geom_point() + theme_bw(base_size = 12) + labs(title = "A. Elbow criterion", x = "Number of clusters", y = "Total within-cluster SS")
  p2 <- ggplot(dt, aes(AIPFI_2011, AIPFI_2015, color = cluster2)) + geom_point(alpha = 0.55, size = 1.1) + theme_bw(base_size = 12) + labs(title = "B. K-means clusters", x = "AIPFI 2011/2012", y = "AIPFI 2015", color = NULL)
  traj <- melt(dt[, .(id, cluster2, AIPFI_2011, AIPFI_2015)], id.vars = c("id", "cluster2"), variable.name = "wave", value.name = "AIPFI")
  traj[, wave := factor(wave, levels = c("AIPFI_2011", "AIPFI_2015"), labels = c("2011/2012", "2015"))]
  mean_traj <- traj[, .(mean = mean(AIPFI, na.rm = TRUE), q25 = quantile(AIPFI, 0.25, na.rm = TRUE), q75 = quantile(AIPFI, 0.75, na.rm = TRUE)), by = .(cluster2, wave)]
  p3 <- ggplot(mean_traj, aes(wave, mean, group = cluster2, color = cluster2)) + geom_line(linewidth = 1.0) + geom_point(size = 2.5) + geom_errorbar(aes(ymin = q25, ymax = q75), width = 0.08) + theme_bw(base_size = 12) + labs(title = "C. Mean AIPFI trajectories", x = NULL, y = "AIPFI", color = NULL)
  inc <- dt[, .(n = .N, events = sum(event == 1), incidence = mean(event == 1)), by = cluster2]
  p4 <- ggplot(inc, aes(cluster2, incidence)) + geom_col(width = 0.65) + geom_text(aes(label = sprintf("%d/%d\n%.1f%%", events, n, incidence * 100)), vjust = -0.25, size = 3.5) + theme_bw(base_size = 12) + labs(title = "D. CMM incidence", x = NULL, y = "Incidence proportion") + ylim(0, max(inc$incidence) * 1.25)
  combine_gg(list(p1, p2, p3, p4), file.path(DIR_FIG, "Fig3.png"), ncol = 2, w = 12, h = 9)
  out <- list(elbow = wss, cluster_summary = dt[, .(n = .N, events = sum(event == 1), incidence = mean(event == 1), AIPFI_2011_median = median(AIPFI_2011, na.rm = TRUE), AIPFI_2015_median = median(AIPFI_2015, na.rm = TRUE), cum_AIPFI_mean = mean(cum_AIPFI, na.rm = TRUE)), by = cluster2], mean_trajectory = mean_traj)
  write_xlsx(out, file.path(DIR_TAB, "Fig3.out.xlsx"))
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 🚩 Fig4: RCS / spline curves
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
reference_newdata <- function(dt, vars) {
  nd <- dt[1]
  for (v in vars) {
    if (!v %in% names(dt)) next
    if (is.factor(dt[[v]])) nd[[v]] <- factor(levels(dt[[v]])[1], levels = levels(dt[[v]]))
    else if (is.character(dt[[v]])) nd[[v]] <- mode_value(dt[[v]])
    else if (is_binary_var(dt[[v]])) nd[[v]] <- suppressWarnings(as.numeric(mode_value(dt[[v]])))
    else nd[[v]] <- median(dt[[v]], na.rm = TRUE)
  }
  nd
}

spline_curve <- function(dt, xvar, covars = covars_m4, df = 4) {
  cv <- covars[covars %in% names(dt)]
  model_vars <- unique(c("time", "event", xvar, cv))
  missing_vars <- setdiff(model_vars, names(dt))
  if (length(missing_vars)) {
    stop("spline_curve missing variable(s): ", paste(missing_vars, collapse = ", "), call. = FALSE)
  }
  z <- data.table::copy(dt[complete.cases(dt[, ..model_vars])])
  if (nrow(z) == 0L) stop("spline_curve has no complete rows for ", xvar, call. = FALSE)

  rhs <- function(terms) if (length(terms)) paste(terms, collapse = " + ") else "1"
  f_s <- as.formula(paste("Surv(time, event) ~", rhs(c(sprintf("splines::ns(%s, df = %d)", xvar, df), cv))))
  f_l <- as.formula(paste("Surv(time, event) ~", rhs(c(xvar, cv))))
  f_n <- as.formula(paste("Surv(time, event) ~", rhs(cv)))
  fit_s <- survival::coxph(f_s, data = z, x = TRUE)
  fit_l <- survival::coxph(f_l, data = z, x = TRUE)
  fit_n <- survival::coxph(f_n, data = z, x = TRUE)
  xs <- seq(quantile(z[[xvar]], 0.01, na.rm = TRUE), quantile(z[[xvar]], 0.99, na.rm = TRUE), length.out = 150)
  nd <- reference_newdata(z, unique(c(xvar, cv)))
  nd <- nd[rep(1, length(xs))]
  nd[[xvar]] <- xs
  pr <- predict(fit_s, newdata = nd, type = "lp", se.fit = TRUE)
  ref <- nd[1]; ref[[xvar]] <- median(z[[xvar]], na.rm = TRUE)
  ref_lp <- as.numeric(predict(fit_s, newdata = ref, type = "lp"))
  out <- data.table(x = xs, HR = exp(pr$fit - ref_lp), LCL = exp(pr$fit - ref_lp - 1.96 * pr$se.fit), UCL = exp(pr$fit - ref_lp + 1.96 * pr$se.fit))
  p_overall <- anova(fit_n, fit_s, test = "LRT")$`Pr(>|Chi|)`[2]
  p_nonlin <- anova(fit_l, fit_s, test = "LRT")$`Pr(>|Chi|)`[2]
  list(curve = out, p = data.table(xvar = xvar, p_overall = p_overall, p_nonlinear = p_nonlin), fit = fit_s)
}

plot_spline <- function(sc, title, xlab) {
  ptxt <- sprintf("P overall=%s; P nonlinear=%s", fmt_p(sc$p$p_overall), fmt_p(sc$p$p_nonlinear))
  ggplot(sc$curve, aes(x, HR)) + geom_ribbon(aes(ymin = LCL, ymax = UCL), alpha = 0.18) + geom_line(linewidth = 0.9) + geom_hline(yintercept = 1, linetype = 2) + theme_bw(base_size = 12) + labs(title = title, subtitle = ptxt, x = xlab, y = "Hazard ratio")
}

run_fig4 <- function() {
  msg("🚩 Fig4: spline associations for AIPFI and cumulative AIPFI")
  b <- load_dat("baseline"); l <- load_dat("longitudinal")
  s1 <- spline_curve(b, "AIPFI", covars_m4)
  s2 <- spline_curve(l, "cum_AIPFI", covars_m4)
  p1 <- plot_spline(s1, "A. Baseline AIPFI", "Baseline AIPFI")
  p2 <- plot_spline(s2, "B. Cumulative AIPFI", "Cumulative AIPFI")
  combine_gg(list(p1, p2), file.path(DIR_FIG, "Fig4.png"), ncol = 2, w = 12, h = 5.5)
  write_xlsx(list(AIPFI_curve = s1$curve, cum_AIPFI_curve = s2$curve, P_values = rbind(s1$p, s2$p)), file.path(DIR_TAB, "Fig4.out.xlsx"))
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 🚩 Fig5: subgroup analyses
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
subgroup_cox <- function(dt, exposure, subgroup_vars, covars = covars_m4) {
  rbindlist(lapply(subgroup_vars, function(sg) {
    if (!sg %in% names(dt)) return(NULL)
    z <- dt[!is.na(get(sg))]
    lev <- levels(as.factor(z[[sg]]))
    per <- rbindlist(lapply(lev, function(lv) {
      zz <- z[as.character(get(sg)) == lv]
      if (nrow(zz) < 50 || sum(zz$event == 1) < 5) return(NULL)
      fit <- try(cox_fit(zz, exposure, setdiff(covars, sg)), silent = TRUE)
      if (inherits(fit, "try-error")) return(NULL)
      ex <- cox_extract(fit)[term == exposure]
      ex[, `:=`(subgroup = sg, level = lv, n = nrow(zz), events = sum(zz$event == 1))]
      ex
    }), fill = TRUE)
    fit_int <- try(cox_fit(z, paste0(exposure, "*", sg), covars), silent = TRUE)
    pint <- NA_real_
    if (!inherits(fit_int, "try-error")) {
      cf <- summary(fit_int)$coefficients
      ints <- grep(":", rownames(cf), value = TRUE)
      if (length(ints)) pint <- min(cf[ints, ncol(cf)], na.rm = TRUE)
    }
    per[, P_interaction := fmt_p(pint)]
    per
  }), fill = TRUE)
}

plot_forest <- function(res, title) {
  res <- copy(res)
  res[, label := paste0(subgroup, ": ", level, " (", events, "/", n, ")")]
  res[, label := factor(label, levels = rev(label))]
  ggplot(res, aes(HR, label)) + geom_vline(xintercept = 1, linetype = 2) + geom_errorbarh(aes(xmin = LCL, xmax = UCL), height = 0.2) + geom_point(size = 2) + scale_x_log10() + theme_bw(base_size = 11) + labs(title = title, x = "Hazard ratio", y = NULL)
}

run_fig5 <- function() {
  msg("🚩 Fig5: subgroup analyses")
  b <- load_dat("baseline"); l <- load_dat("longitudinal")
  sg <- c("age_group", "gender", "hypertension", "diabetes", "dyslipidemia", "heart_disease", "stroke", "smoking", "drinking")
  r1 <- subgroup_cox(b, "AIPFI_sd", sg)
  r2 <- subgroup_cox(l, "cum_AIPFI_sd", sg)
  write_xlsx(list(baseline_AIPFI = r1, cumulative_AIPFI = r2), file.path(DIR_TAB, "Fig5.out.xlsx"))
  p1 <- plot_forest(r1, "A. Baseline AIPFI per 1 SD")
  p2 <- plot_forest(r2, "B. Cumulative AIPFI per 1 SD")
  combine_gg(list(p1, p2), file.path(DIR_FIG, "Fig5.png"), ncol = 2, w = 14, h = 8.5)
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 🚩 Fig6: clinical prediction model
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
auc_rank <- function(y, score) {
  ok <- !is.na(y) & !is.na(score)
  y <- y[ok]; score <- score[ok]
  if (length(unique(y)) < 2) return(NA_real_)
  r <- rank(score)
  n1 <- sum(y == 1); n0 <- sum(y == 0)
  (sum(r[y == 1]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}

predict_risk_cox <- function(fit, newdata, t0 = 8) {
  bh <- survival::basehaz(fit, centered = FALSE)
  h0 <- approx(bh$time, bh$hazard, xout = t0, rule = 2)$y
  lp <- predict(fit, newdata = newdata, type = "lp")
  1 - exp(-h0 * exp(lp))
}

roc_curve <- function(y, score) {
  ok <- !is.na(y) & !is.na(score)
  y <- y[ok]; score <- score[ok]
  thr <- sort(unique(score), decreasing = TRUE)
  data.table(threshold = thr, TPR = sapply(thr, function(t) mean(score[y == 1] >= t)), FPR = sapply(thr, function(t) mean(score[y == 0] >= t)))
}

calibration_dt <- function(dt, risk, t0 = 8, g = 10) {
  z <- copy(dt); z[, risk := risk]
  z <- z[!is.na(risk)]
  z[, decile := cut(rank(risk, ties.method = "first"), breaks = quantile(rank(risk, ties.method = "first"), seq(0, 1, length.out = g + 1)), include.lowest = TRUE, labels = paste0("D", 1:g))]
  z[, .(n = .N, predicted = mean(risk), observed = mean(event == 1 & time <= t0)), by = decile]
}

dca_dt <- function(y, risk, thresholds = seq(0.01, 0.30, by = 0.01)) {
  n <- length(y)
  rbindlist(lapply(thresholds, function(pt) {
    pred <- risk >= pt
    tp <- sum(pred & y == 1, na.rm = TRUE); fp <- sum(pred & y == 0, na.rm = TRUE)
    data.table(threshold = pt, net_benefit = tp / n - fp / n * pt / (1 - pt))
  }))
}

run_fig6 <- function() {
  msg("🚩 Fig6: prediction model, ROC, calibration, DCA, SHAP-like explanation")
  dt <- load_dat("baseline")
  set.seed(SEED)
  idx <- sample.int(nrow(dt), floor(0.70 * nrow(dt)))
  train <- dt[idx]; valid <- dt[-idx]
  pred_vars <- c("diabetes", "hypertension", "heart_disease", "stroke", "dyslipidemia", "FBG", "HbA1c", "BMI", "TG", "AIPFI")
  fit <- cox_fit(train, paste(pred_vars, collapse = " + "), character())
  risk_tr <- predict_risk_cox(fit, train, 8)
  risk_va <- predict_risk_cox(fit, valid, 8)
  y_tr <- as.numeric(train$event == 1 & train$time <= 8)
  y_va <- as.numeric(valid$event == 1 & valid$time <= 8)
  roc_tr <- roc_curve(y_tr, risk_tr); roc_va <- roc_curve(y_va, risk_va)
  aucs <- data.table(cohort = c("Training", "Validation"), n = c(nrow(train), nrow(valid)), events = c(sum(train$event == 1), sum(valid$event == 1)), AUC_96m = c(auc_rank(y_tr, risk_tr), auc_rank(y_va, risk_va)), C_index = c(summary(fit)$concordance[1], survival::concordance(Surv(time, event) ~ risk_va, data = valid)$concordance))
  cal_tr <- calibration_dt(train, risk_tr, 8); cal_tr[, cohort := "Training"]
  cal_va <- calibration_dt(valid, risk_va, 8); cal_va[, cohort := "Validation"]
  dca_tr <- dca_dt(y_tr, risk_tr); dca_tr[, cohort := "Training"]
  dca_va <- dca_dt(y_va, risk_va); dca_va[, cohort := "Validation"]
  # True SHAP via fastshap if available; otherwise Cox linear-predictor contribution approximation.
  X <- as.data.frame(train[, ..pred_vars])
  shap_note <- "linear predictor contribution beta*x; install fastshap for permutation SHAP"
  beta <- coef(fit)[pred_vars]
  contr <- sweep(as.matrix(X[, names(beta), drop = FALSE]), 2, beta, `*`)
  shap_long <- data.table(feature = rep(colnames(contr), each = nrow(contr)), value = as.vector(as.matrix(X[, colnames(contr), drop = FALSE])), shap = as.vector(contr))
  shap_imp <- shap_long[, .(mean_abs_SHAP = mean(abs(shap), na.rm = TRUE)), by = feature][order(-mean_abs_SHAP)]
  if (opt_pkg("fastshap")) {
    pred_fun <- function(object, newdata) predict_risk_cox(object, as.data.table(newdata), 8)
    set.seed(SEED)
    sh <- fastshap::explain(fit, X = X, pred_wrapper = pred_fun, nsim = 200, adjust = TRUE)
    shap_long <- data.table(feature = rep(colnames(sh), each = nrow(sh)), value = as.vector(as.matrix(X[, colnames(sh), drop = FALSE])), shap = as.vector(as.matrix(sh)))
    shap_imp <- shap_long[, .(mean_abs_SHAP = mean(abs(shap), na.rm = TRUE)), by = feature][order(-mean_abs_SHAP)]
    shap_note <- "fastshap permutation SHAP with 200 simulations"
  }
  shap_imp[, note := shap_note]

  p1 <- ggplot(rbind(cbind(roc_tr, cohort = "Training"), cbind(roc_va, cohort = "Validation")), aes(FPR, TPR, color = cohort)) + geom_line(linewidth = 0.9) + geom_abline(slope = 1, intercept = 0, linetype = 2) + theme_bw(base_size = 12) + labs(title = "A/B. ROC at 96 months", subtitle = paste(sprintf("Training AUC=%.3f", aucs$AUC_96m[1]), sprintf("Validation AUC=%.3f", aucs$AUC_96m[2]), sep = "; "), x = "1 - specificity", y = "Sensitivity", color = NULL)
  p2 <- ggplot(rbind(cal_tr, cal_va), aes(predicted, observed, color = cohort)) + geom_point() + geom_line() + geom_abline(slope = 1, intercept = 0, linetype = 2) + theme_bw(base_size = 12) + labs(title = "Calibration", x = "Mean predicted risk", y = "Observed risk", color = NULL)
  p3 <- ggplot(rbind(dca_tr, dca_va), aes(threshold, net_benefit, color = cohort)) + geom_line(linewidth = 0.9) + theme_bw(base_size = 12) + labs(title = "Decision curve analysis", x = "Threshold probability", y = "Net benefit", color = NULL)
  p4 <- ggplot(shap_imp, aes(mean_abs_SHAP, reorder(feature, mean_abs_SHAP))) + geom_col() + theme_bw(base_size = 12) + labs(title = "Mean absolute SHAP", subtitle = shap_note, x = "Mean absolute contribution", y = NULL)
  combine_gg(list(p1, p2, p3, p4), file.path(DIR_FIG, "Fig6.png"), ncol = 2, w = 12, h = 9)
  write_xlsx(list(AUC = aucs, ROC_training = roc_tr, ROC_validation = roc_va, calibration = rbind(cal_tr, cal_va), DCA = rbind(dca_tr, dca_va), SHAP_importance = shap_imp, SHAP_values = shap_long), file.path(DIR_TAB, "Fig6.out.xlsx"))
  write_xlsx(make_desc_table(train, group = "event", vars = c(pred_vars, covars_m4)), file.path(DIR_TAB, "TableS25.out.xlsx"))
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 🚩 Supplementary sensitivity analyses
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
vif_table <- function(dt, vars) {
  vars <- vars[vars %in% names(dt)]
  dd <- dt[, ..vars]
  dd <- simple_impute(dd, vars)
  rbindlist(lapply(vars, function(v) {
    others <- setdiff(vars, v)
    fit <- try(lm(as.formula(paste(v, "~", paste(others, collapse = "+"))), data = dd), silent = TRUE)
    if (inherits(fit, "try-error")) return(data.table(variable = v, VIF = NA_real_))
    r2 <- summary(fit)$r.squared
    data.table(variable = v, VIF = 1 / (1 - r2))
  }))
}

cox_model4_exposure <- function(dt, exposure, label) {
  fit <- cox_fit(dt, exposure, covars_m4)
  cox_extract(fit, label)
}

run_supp <- function() {
  msg("🚩 Supplementary analyses")
  b <- load_dat("baseline"); l <- load_dat("longitudinal")
  out <- list()
  # S3/S4/S5/S6 descriptive tables.
  out$TableS3_missing_proportions <- data.table(variable = names(b), missing_prop = vapply(b, function(x) mean(is.na(x)), numeric(1)))
  out$TableS4_baseline_by_outcome_2013_2020 <- make_desc_table(b, group = "event")
  out$TableS5_longitudinal_by_outcome_2018_2020 <- make_desc_table(l, group = "event", vars = c("age", "gender", "diabetes", "hypertension", "heart_disease", "stroke", "dyslipidemia", "drinking", "smoking", "SBP", "DBP", "BMI", "FBG", "HbA1c", "TC", "TG", "HDLC", "LDLC", "FI_2011", "FI_2015", "AIP_2011", "AIP_2015", "AIPFI_2011", "AIPFI_2015", "cum_AIPFI", "cum_AIPFI_Q"))
  out$TableS6_by_AIPFI_cluster <- make_desc_table(l, group = "cluster2")
  # S7 medication-adjusted model.
  meds <- c("rx_diab", "rx_htn", "rx_lipid", "rx_heart", "rx_stroke")
  cov_meds <- unique(c(covars_m4, meds[meds %in% names(b)]))
  out$TableS7_medication_adjusted <- run_cox_models(b, "AIPFI_Q", "AIPFI quartile", list(Model4_plus_meds = cov_meds))[grepl("AIPFI_Q", term)]
  # S8 untreated baseline subset.
  if (all(meds %in% names(b))) {
    untreated_i <- rowSums(as.data.frame(b[, ..meds]) == 1, na.rm = TRUE) == 0
    sub <- b[untreated_i]
    out$TableS8_no_treatment <- run_cox_models(sub, "AIPFI_Q", "AIPFI quartile", list(Model4 = covars_m4))[grepl("AIPFI_Q", term)]
  }
  # S9 no baseline CMD components.
  sub0 <- b[diabetes == 0 & heart_disease == 0 & stroke == 0]
  out$TableS9_no_baseline_CMD <- run_cox_models(sub0, "AIPFI_Q", "AIPFI quartile", list(Model4 = covars_m4))[grepl("AIPFI_Q", term)]
  # S10-S12 individual outcomes. Recompute first events for each component.
  waves <- load_dat("waves")
  for (component in c("diabetes", "stroke", "heart_disease")) {
    fu <- waves[wave %in% c(2,3,4,5), .(id, wave, year, y = get(component))][!is.na(y)]
    ev <- fu[y == 1, .SD[which.min(year)], by = id]
    last <- fu[, .SD[which.max(year)], by = id]
    tmp <- copy(b)
    prior_surv_cols <- intersect(c("event", "time", "event_year", "event_wave", "last_year"), names(tmp))
    tmp <- merge(tmp[, !prior_surv_cols, with = FALSE], ev[, .(id, event_year = year)], by = "id", all.x = TRUE)
    tmp <- merge(tmp, last[, .(id, last_year = year)], by = "id", all.x = TRUE)
    tmp[, event := as.numeric(!is.na(event_year))]
    tmp[, time := fifelse(event == 1, event_year - 2011, last_year - 2011)]
    tmp <- tmp[!is.na(time) & time > 0]
    nm <- paste0("TableS", switch(component, diabetes = "10", stroke = "11", heart_disease = "12"), "_", component)
    out[[nm]] <- run_cox_models(tmp, "AIPFI_Q", "AIPFI quartile", list(Model4 = covars_m4))[grepl("AIPFI_Q", term)]
  }
  # S13 2015 baseline.
  b2015 <- copy(l)
  b2015[, AIPFI_Q_2015 := make_quartile(AIPFI_2015)]
  b2015[, AIPFI_sd_2015 := as.numeric(scale(AIPFI_2015))]
  out$TableS13_wave2015_baseline <- run_cox_models(b2015, "AIPFI_Q_2015", "AIPFI 2015 quartile", list(Model4 = covars_m4))[grepl("AIPFI_Q_2015", term)]
  # S14 standardization/sign/outlier removal.
  b[, AIPFI_winsor := winsor(AIPFI)]
  b[, AIPFI_sign := factor(ifelse(AIPFI < 0, "<0", ">=0"))]
  b[, AIPFI_winsor_sd := as.numeric(scale(AIPFI_winsor))]
  out$TableS14_standardization_outlier_sign <- rbindlist(list(
    run_cox_models(b, "AIPFI_sd", "standardized AIPFI", list(Model4 = covars_m4))[term == "AIPFI_sd"],
    run_cox_models(b, "AIPFI_winsor_sd", "winsorized standardized AIPFI", list(Model4 = covars_m4))[term == "AIPFI_winsor_sd"],
    run_cox_models(b, "AIPFI_sign", "AIPFI sign", list(Model4 = covars_m4))[grepl("AIPFI_sign", term)]
  ), fill = TRUE)
  # S16 IPW high vs low.
  b[, AIPFI_high := as.numeric(AIPFI > median(AIPFI, na.rm = TRUE))]
  psfit <- glm(as.formula(paste("AIPFI_high ~", paste(covars_m4[covars_m4 %in% names(b)], collapse = "+"))), data = b, family = binomial())
  ps <- pmin(pmax(predict(psfit, type = "response"), 0.01), 0.99)
  wt <- ifelse(b$AIPFI_high == 1, mean(b$AIPFI_high) / ps, (1 - mean(b$AIPFI_high)) / (1 - ps))
  fit_ipw <- cox_fit(b, "AIPFI_high", covars_m4, weights = wt)
  out$TableS16_IPW_weighted_Cox <- cox_extract(fit_ipw, "AIPFI high vs low")[term == "AIPFI_high"]
  # S18/S19 VIF and reduced model.
  vifvars <- c("TC", "LDLC", "TG", "HDLC", "SBP", "DBP", "gender", "smoking", "FBG", "AIPFI", "HbA1c", "drinking", "hypertension", "BMI", "age", "diabetes", "dyslipidemia", "heart_disease", "stroke")
  out$TableS18_VIF <- vif_table(b, vifvars)
  red_cov <- setdiff(covars_m4, c(out$TableS18_VIF[VIF > 5]$variable, "TG"))
  out$TableS19_exclude_high_VIF_and_TG <- run_cox_models(b, "AIPFI_Q", "AIPFI quartile", list(Model4_reduced = red_cov))[grepl("AIPFI_Q", term)]
  # S20 three clusters.
  cl3 <- cluster_aipfi(l, k = 3)
  l3 <- merge(l, cl3$assign, by = "id", all.x = TRUE)
  l3[, cluster3 := factor(cluster3)]
  out$TableS20_three_clusters <- run_cox_models(l3, "cluster3", "AIPFI 3-cluster", list(Model4 = covars_m4))[grepl("cluster3", term)]
  # S21/S22 PH and stratified Cox.
  fit_m4 <- cox_fit(b, "AIPFI", covars_m4)
  zph <- survival::cox.zph(fit_m4)
  out$TableS21_PH_test <- data.table(term = rownames(zph$table), chisq = zph$table[, "chisq"], df = zph$table[, "df"], p = zph$table[, "p"])
  fit_str <- cox_fit(b, "AIPFI", covars_m4, strata_vars = c("diabetes", "hypertension", "heart_disease", "stroke"))
  out$TableS22_stratified_Cox <- cox_extract(fit_str, "AIPFI stratified Cox")
  # S23 incremental prediction using AUC difference as a pragmatic output.
  y <- as.numeric(b$event == 1 & b$time <= 8)
  base_fit <- cox_fit(b, paste(c("age", "gender", "hypertension", "dyslipidemia", "diabetes"), collapse = " + "), character())
  full_fit <- cox_fit(b, paste(c("age", "gender", "hypertension", "dyslipidemia", "diabetes", "AIPFI"), collapse = " + "), character())
  rb <- predict_risk_cox(base_fit, b, 8); rf <- predict_risk_cox(full_fit, b, 8)
  out$TableS23_incremental_prediction <- data.table(model = c("traditional", "traditional+AIPFI", "delta"), AUC = c(auc_rank(y, rb), auc_rank(y, rf), auc_rank(y, rf) - auc_rank(y, rb)), note = "NRI/IDI requires survIDINRI/riskRegression; this sheet provides AUC increment.")
  write_xlsx(out, file.path(DIR_TAB, "Supplementary_Tables.out.xlsx"))
  # Supplementary figures.
  if (opt_pkg("cmprsk")) {
    # No death variable is reliably harmonized here; Fine-Gray is not run unless user adds competing event variable.
    write_xlsx(data.table(note = "cmprsk is installed, but no validated death/competing event variable was resolved by default. Add death status before Fine-Gray reproduction."), file.path(DIR_TAB, "FigS1.out.xlsx"))
  } else {
    write_xlsx(data.table(note = "Install cmprsk and add a validated death/competing event variable to reproduce Fine-Gray Fig.S1."), file.path(DIR_TAB, "FigS1.out.xlsx"))
  }
  # Fig S2: 3-cluster version.
  p_s2 <- ggplot(l3, aes(AIPFI_2011, AIPFI_2015, color = cluster3)) + geom_point(alpha = 0.55, size = 1.1) + theme_bw(base_size = 12) + labs(title = "Fig.S2 Three-cluster AIPFI trajectories", x = "AIPFI 2011/2012", y = "AIPFI 2015", color = NULL)
  save_plot(p_s2, file.path(DIR_FIG, "FigS2.png"))
  # Fig S3 forest plot of stratified model.
  s22 <- out$TableS22_stratified_Cox[!grepl("strata", term)]
  s22[, term := factor(term, levels = rev(term))]
  p_s3 <- ggplot(s22, aes(HR, term)) + geom_vline(xintercept = 1, linetype = 2) + geom_errorbarh(aes(xmin = LCL, xmax = UCL), height = 0.2) + geom_point() + scale_x_log10() + theme_bw(base_size = 11) + labs(title = "Fig.S3 Multivariable stratified Cox model", x = "Hazard ratio", y = NULL)
  save_plot(p_s3, file.path(DIR_FIG, "FigS3.png"), w = 7.5, h = 7)
  # Fig S4 ROC for AIP, FI, AIPFI.
  roc3 <- rbindlist(lapply(c("AIP", "FI", "AIPFI"), function(v) { z <- roc_curve(y, b[[v]]); z[, marker := v]; z[, AUC := auc_rank(y, b[[v]])]; z }))
  p_s4 <- ggplot(roc3, aes(FPR, TPR, color = marker)) + geom_line(linewidth = 0.9) + geom_abline(slope = 1, intercept = 0, linetype = 2) + theme_bw(base_size = 12) + labs(title = "Fig.S4 ROC analysis of AIP, FI, and AIPFI", x = "1 - specificity", y = "Sensitivity", color = NULL)
  save_plot(p_s4, file.path(DIR_FIG, "FigS4.png"))
  write_xlsx(roc3, file.path(DIR_TAB, "FigS4.out.xlsx"))
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 🚩 QC helper
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
run_qc <- function() {
  cat("Output root: ", CHARLS_OUT, "\n", sep = "")
  cat("Key files:\n")
  cat("  ", file.path(DIR_QC, "var_map.xlsx"), "\n", sep = "")
  cat("  ", file.path(DIR_QC, "fi_item_map.xlsx"), "\n", sep = "")
  cat("  ", file.path(DIR_QC, "cohort_counts.xlsx"), "\n", sep = "")
  cat("  ", file.path(DIR_QC, "manifest_file_choice.xlsx"), "\n", sep = "")
  cat("  ", file.path(DIR_DAT, "baseline_2011_2020.rds"), "\n", sep = "")
  cat("  ", file.path(DIR_DAT, "longitudinal_2015_2020.rds"), "\n", sep = "")
  cat("  ", file.path(DIR_TAB, "Tab1.out.xlsx"), "\n", sep = "")
  cat("  ", file.path(DIR_FIG, "Fig1.png"), "\n", sep = "")
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 🚩 Workflow dispatcher
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
run_step <- function(step) {
  switch(step,
    prep = run_prep(),
    fig1 = run_fig1(),
    tab1 = run_tab1(),
    tab2 = run_tab2(),
    tab3 = run_tab3(),
    fig2 = run_fig2(),
    fig3 = run_fig3(),
    fig4 = run_fig4(),
    fig5 = run_fig5(),
    fig6 = run_fig6(),
    supp = run_supp(),
    qc = run_qc(),
    stop("Unknown step: ", step, call. = FALSE)
  )
}

all_steps <- c("prep", "fig1", "tab1", "tab2", "tab3", "fig2", "fig3", "fig4", "fig5", "fig6", "supp")
if (MODE == "all") {
  steps <- all_steps
  if (nzchar(START_STEP)) {
    i <- match(START_STEP, steps)
    if (is.na(i)) stop("Invalid START_STEP=", START_STEP, "; valid: ", paste(steps, collapse = ", "), call. = FALSE)
    steps <- steps[i:length(steps)]
  }
  for (s in steps) run_step(s)
} else {
  run_step(MODE)
}

msg("DONE ", MODE)
