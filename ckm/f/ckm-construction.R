#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# CKM 2026 data construction
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

resolve_col <- function(dat, env_name, candidates) {
  x <- Sys.getenv(env_name, unset = "")
  if (nzchar(x)) {
    if (!x %in% names(dat)) stop(env_name, " requests missing column: ", x, call. = FALSE)
    return(x)
  }
  first_present(dat, candidates)
}

ckd_risk_category <- function(egfr, uacr) {
  g <- cut(egfr, breaks = c(-Inf, 15, 30, 45, 60, 90, Inf),
           labels = c("G5", "G4", "G3b", "G3a", "G2", "G1"), right = FALSE)
  a <- cut(uacr, breaks = c(-Inf, 30, 300, Inf), labels = c("A1", "A2", "A3"), right = FALSE)
  key <- paste(g, a, sep = "_")
  map <- c(
    G1_A1 = "low", G1_A2 = "moderate", G1_A3 = "high",
    G2_A1 = "low", G2_A2 = "moderate", G2_A3 = "high",
    G3a_A1 = "moderate", G3a_A2 = "high", G3a_A3 = "very_high",
    G3b_A1 = "high", G3b_A2 = "very_high", G3b_A3 = "very_high",
    G4_A1 = "very_high", G4_A2 = "very_high", G4_A3 = "very_high",
    G5_A1 = "very_high", G5_A2 = "very_high", G5_A3 = "very_high"
  )
  out <- unname(map[key])
  # When UACR is unavailable, eGFR alone still identifies G3a/G3b/G4/G5 risk.
  out[is.na(out) & is.finite(egfr) & egfr < 60 & egfr >= 45] <- "moderate"
  out[is.na(out) & is.finite(egfr) & egfr < 45 & egfr >= 30] <- "high"
  out[is.na(out) & is.finite(egfr) & egfr < 30] <- "very_high"
  factor(out, levels = c("low", "moderate", "high", "very_high"), ordered = TRUE)
}

calculate_prevent <- function(dat, statin_col) {
  if (!requireNamespace("preventr", quietly = TRUE)) {
    progress_log("Installing missing R package: preventr")
    install_try <- try(install.packages("preventr", repos = "https://cloud.r-project.org"), silent = TRUE)
    if (inherits(install_try, "try-error") || !requireNamespace("preventr", quietly = TRUE)) {
      stop("The preventr package is required to construct CKM stage 3 but could not be installed", call. = FALSE)
    }
  }
  progress_log("PREVENT: preparing uniform base-model primary and best-available sensitivity inputs")
  use <- dat %>% transmute(
    eid,
    age = safe_num(age),
    sex = if_else(safe_num(sex) == 1, "male", "female"),
    sbp = safe_num(ckm_sbp),
    bp_tx = safe_num(drug.htn) == 1,
    total_c = safe_num(bb_TC),
    hdl_c = safe_num(bb_HDL),
    statin = safe_num(.data[[statin_col]]) == 1,
    dm = T2D_def,
    smoking = safe_num(smoke_status) == 2,
    egfr = safe_num(egfr),
    bmi = safe_num(bmi),
    hba1c = hba1c_pct,
    uacr = uacr_mg_g
  )
  eligible <- with(use,
    age >= 30 & age <= 79 & !dat$prior_cvd &
      is.finite(sbp) & is.finite(total_c) & is.finite(hdl_c) &
      is.finite(egfr) & is.finite(bmi)
  )
  out <- data.table(
    eid = use$eid,
    prevent_cvd10_base = NA_real_, prevent_model_base = NA_character_, prevent_problem_base = NA_character_,
    prevent_cvd10_best = NA_real_, prevent_model_best = NA_character_, prevent_problem_best = NA_character_
  )
  eligible_idx <- which(eligible)
  progress_log("PREVENT: eligible rows = ", fmt_int(length(eligible_idx)), " / ", fmt_int(nrow(use)))
  if (!length(eligible_idx)) return(out)

  batch_size <- suppressWarnings(as.integer(Sys.getenv("CKM_PREVENT_BATCH_SIZE", unset = "5000")))
  if (!is.finite(batch_size) || batch_size <= 0) batch_size <- 5000L
  batch_size <- max(100L, batch_size)
  chunks <- split(eligible_idx, ceiling(seq_along(eligible_idx) / batch_size))

  run_variant <- function(cols, prefix, model = NULL) {
    ans <- data.table(row_index = eligible_idx, risk = NA_real_, model = NA_character_, problem = NA_character_)
    done <- 0L
    progress_log("PREVENT ", prefix, ": estimating in ", length(chunks), " batch(es); model=", model %||% "auto")
    for (bb in seq_along(chunks)) {
      ii <- chunks[[bb]]
      t0 <- Sys.time()
      call_args <- list(
        use_dat = as.data.frame(use[ii, cols, drop = FALSE]),
        time = "10yr", chol_unit = "mmol/L", add_to_dat = FALSE,
        progress = FALSE, quiet = TRUE
      )
      if (!is.null(model)) call_args$model <- model
      rr <- try(do.call(preventr::estimate_risk, call_args), silent = TRUE)
      if (inherits(rr, "try-error")) {
        stop(
          "PREVENT ", prefix, " failed in batch ", bb,
          ". Partial PREVENT results will not be used for CKM staging. Error: ",
          as.character(rr), call. = FALSE
        )
      }
      rr <- as.data.table(rr)
      if (!"preventr_id" %in% names(rr)) rr[, preventr_id := seq_len(.N)]
      pid <- suppressWarnings(as.integer(rr$preventr_id))
      valid <- is.finite(pid) & pid >= 1L & pid <= length(ii)
      rr <- rr[valid]
      pid <- pid[valid]
      dest <- match(ii[pid], ans$row_index)
      col_or_na <- function(nm, default) if (nm %in% names(rr)) rr[[nm]] else rep(default, nrow(rr))
      ans$risk[dest] <- col_or_na("total_cvd", NA_real_)
      ans$model[dest] <- col_or_na("model", NA_character_)
      ans$problem[dest] <- col_or_na("input_problems", NA_character_)
      done <- done + length(ii)
      progress_log("PREVENT ", prefix, " batch ", bb, "/", length(chunks), " done in ", elapsed_since(t0),
                   "; completed=", fmt_int(done), "/", fmt_int(length(eligible_idx)))
    }
    ans
  }

  base_cols <- c("age", "sex", "sbp", "bp_tx", "total_c", "hdl_c", "statin", "dm", "smoking", "egfr", "bmi")
  best_cols <- c(base_cols, "hba1c", "uacr")
  # Research staging must use one equation across participants. Explicitly force
  # the base PREVENT model for the primary stage-3 definition.
  variant_cores <- suppressWarnings(as.integer(Sys.getenv("CKM_PREVENT_VARIANT_CORES", unset = "1")))
  if (!is.finite(variant_cores) || variant_cores < 1L) variant_cores <- 1L
  jobs <- list(
    function() run_variant(base_cols, "base", model = "base"),
    # Auto-selection of the best valid optional model is retained only as a
    # sensitivity analysis.
    function() run_variant(best_cols, "best_available", model = NULL)
  )
  if (variant_cores > 1L && .Platform$OS.type != "windows") {
    progress_log("PREVENT: running base and best_available variants in parallel")
    variants <- parallel::mclapply(jobs, function(f) f(), mc.cores = min(2L, variant_cores))
  } else {
    variants <- lapply(jobs, function(f) f())
  }
  rb <- variants[[1L]]
  re <- variants[[2L]]
  out$prevent_cvd10_base[rb$row_index] <- rb$risk
  out$prevent_model_base[rb$row_index] <- rb$model
  out$prevent_problem_base[rb$row_index] <- rb$problem
  out$prevent_cvd10_best[re$row_index] <- re$risk
  out$prevent_model_best[re$row_index] <- re$model
  out$prevent_problem_best[re$row_index] <- re$problem
  out
}


# Exact CKM stage-4 ICD-10 mapping supplied by the phenotype pipeline.
# No ICD9/SRD/GP/REF or derived date_mi/date_stroke variables are used here.
ckm_stage4_specification <- function(dat) {
  spec <- data.table::data.table(
    domain = c("chd", "chd", "chd", "hf", "stroke", "pad", "af"),
    trait = c("cvd_cad", "mi", "ihd", "cvd_hfail", "stroke", "cvd_pad", "cvd_afib"),
    expected_col = paste0("fod_icd10_", c("cvd_cad", "mi", "ihd", "cvd_hfail", "stroke", "cvd_pad", "cvd_afib")),
    required = c(TRUE, TRUE, FALSE, TRUE, TRUE, TRUE, TRUE),
    note = c(
      "CHD component", "CHD component", "optional: not defined in the supplied icd10.lst",
      "heart failure", "stroke", "peripheral arterial disease", "atrial fibrillation"
    )
  )
  spec[, exists_in_all := expected_col %in% names(dat)]

  spec[, `:=`(listed_in_icd10_lst = NA, icd10_pattern = NA_character_)]
  if (file.exists(icd10_list_file)) {
    defs <- try(data.table::fread(icd10_list_file, sep = "\t", header = FALSE,
                                  select = 1:2, fill = TRUE, data.table = TRUE,
                                  showProgress = FALSE), silent = TRUE)
    if (!inherits(defs, "try-error") && nrow(defs)) {
      data.table::setnames(defs, c("icd10_pattern", "trait"))
      defs[, trait := trimws(as.character(trait))]
      defs[, icd10_pattern := as.character(icd10_pattern)]
      spec[, listed_in_icd10_lst := trait %in% defs$trait]
      spec[, icd10_pattern := defs$icd10_pattern[match(trait, defs$trait)]]
    }
  }

  missing_required <- spec[required & !exists_in_all]
  if (nrow(missing_required)) {
    stop(
      "Missing required CKM stage-4 ICD-10 FOD column(s) in all.rds: ",
      paste(missing_required$expected_col, collapse = ", "),
      ". Re-run RUN_STEPS=icd,merge in phe.sh after confirming common/icd10.lst.",
      call. = FALSE
    )
  }
  spec[]
}

make_exact_ckm4_dates <- function(dat, baseline_col = "date_attend") {
  spec <- ckm_stage4_specification(dat)
  n <- nrow(dat)
  domain_order <- c("chd", "hf", "stroke", "pad", "af")
  dates <- vector("list", length(domain_order)); names(dates) <- domain_order

  for (domain in domain_order) {
    dom <- domain
    cols <- spec[spec$domain == dom & spec$exists_in_all, expected_col]
    dates[[domain]] <- if (length(cols)) {
      pmin_date2(dat[, cols, drop = FALSE])
    } else {
      as.Date(rep(NA, n))
    }
  }
  all_date <- pmin_date2(as.data.frame(dates))

  # Record which domain supplied the first event; retain tied same-day domains.
  first_domain <- rep(NA_character_, n)
  for (domain in domain_order) {
    hit <- !is.na(dates[[domain]]) & !is.na(all_date) & dates[[domain]] == all_date
    first_domain[hit & is.na(first_domain)] <- domain
    tied <- hit & !is.na(first_domain) & first_domain != domain & !grepl(paste0("(^|\\+)", domain, "($|\\+)"), first_domain)
    first_domain[tied] <- paste(first_domain[tied], domain, sep = "+")
  }

  baseline <- as_date2(dat[[baseline_col]])
  spec[, n_with_date := vapply(expected_col, function(cc) {
    if (!cc %in% names(dat)) return(0L)
    sum(!is.na(as_date2(dat[[cc]])))
  }, integer(1))]
  spec[, n_prebaseline := vapply(expected_col, function(cc) {
    if (!cc %in% names(dat)) return(0L)
    d <- as_date2(dat[[cc]])
    sum(!is.na(d) & !is.na(baseline) & d <= baseline)
  }, integer(1))]
  spec[, n_postbaseline := vapply(expected_col, function(cc) {
    if (!cc %in% names(dat)) return(0L)
    d <- as_date2(dat[[cc]])
    sum(!is.na(d) & !is.na(baseline) & d > baseline)
  }, integer(1))]
  spec[, n_first_domain := vapply(domain, function(dd) {
    sum(!is.na(first_domain) & grepl(paste0("(^|\\+)", dd, "($|\\+)"), first_domain))
  }, integer(1))]

  list(
    domain_dates = dates,
    date_ckm4_all = all_date,
    first_domain = first_domain,
    source_audit = spec
  )
}

# Broad all-source dated definition compatible with the phenotype pipeline's
# t2e/FOD variables.  It is the default primary clinical-CVD ascertainment;
# exact ICD-10 hospital dates are retained as a transparent sensitivity.
make_all_source_ckm4_dates <- function(dat, baseline_col = "date_attend") {
  domain_traits <- list(
    chd = c("cvd_cad", "mi", "ihd"),
    hf = c("cvd_hfail"),
    stroke = c("stroke", "cvd_stroke"),
    pad = c("cvd_pad"),
    af = c("cvd_afib")
  )
  n <- nrow(dat)
  baseline <- as_date2(dat[[baseline_col]])
  source_rows <- list()
  domain_dates <- list()

  for (domain in names(domain_traits)) {
    trait_map <- rbindlist(lapply(domain_traits[[domain]], function(y) {
      cols <- candidate_date_cols(dat, y, "all")
      if (!length(cols)) return(data.table())
      data.table(domain = domain, trait = y, source_col = cols)
    }), fill = TRUE)
    trait_map <- unique(trait_map[source_col %in% names(dat)])
    cols <- unique(trait_map$source_col)
    domain_dates[[domain]] <- if (length(cols)) pmin_date2(dat[, cols, drop = FALSE]) else as.Date(rep(NA, n))

    if (nrow(trait_map)) {
      source_rows[[domain]] <- trait_map[, {
        d <- as_date2(dat[[source_col]])
        .(
          n_with_date = sum(!is.na(d)),
          n_prebaseline = sum(!is.na(d) & !is.na(baseline) & d <= baseline),
          n_postbaseline = sum(!is.na(d) & !is.na(baseline) & d > baseline)
        )
      }, by = .(domain, trait, source_col)]
    }
  }

  all_date <- pmin_date2(as.data.frame(domain_dates))
  source_map <- rbindlist(source_rows, fill = TRUE)
  domain_map <- rbindlist(lapply(names(domain_dates), function(domain) {
    d <- domain_dates[[domain]]
    data.table(
      domain = domain,
      trait = "[domain composite]",
      source_col = "all dated sources",
      n_with_date = sum(!is.na(d)),
      n_prebaseline = sum(!is.na(d) & !is.na(baseline) & d <= baseline),
      n_postbaseline = sum(!is.na(d) & !is.na(baseline) & d > baseline)
    )
  }), fill = TRUE)
  map <- rbind(source_map, domain_map, fill = TRUE)
  list(domain_dates = domain_dates, date_ckm4_all = all_date, source_audit = map)
}


build_ckm_data <- function(force = FALSE) {
  t_build <- Sys.time()
  if (file.exists(ckm_file) && !force) {
    progress_log("ckm.rds already exists; build_ckm step is skipped: ", ckm_file)
    return(invisible(readRDS(ckm_file)))
  }
  # The audit directory is intentionally created only for an actual CKM build.
  # Downstream runs that reuse ckm.rds therefore do not leave an empty audit folder.
  ensure_r_directory(auditdir, "CKM audit directory")
  progress_log("build_ckm: force=", force, "; output=", ckm_file)
  dat <- timed_eval("build_ckm: reading all.rds", readRDS(all_file) %>% as.data.frame())
  progress_log("build_ckm: loaded all.rds rows=", fmt_int(nrow(dat)), "; cols=", fmt_int(ncol(dat)))
  dat$eid <- as.character(dat$eid)
  t0 <- Sys.time()
  progress_log("build_ckm: converting baseline/death date columns")
  for (v in intersect(c("date_attend", "date_lost", "date_death", "birth_date"), names(dat))) dat[[v]] <- as_date2(dat[[v]])
  progress_log("build_ckm: date conversion done in ", elapsed_since(t0))
  sbp_col <- resolve_col(dat, "CKM_SBP_COL", c("sbp", "sbp_auto"))
  dbp_col <- resolve_col(dat, "CKM_DBP_COL", c("dbp", "dbp_auto"))
  req <- c("eid", "date_attend", "age", "sex", "bmi", "waist",
           "bb_TC", "bb_HDL", "bb_TG", "bb_HBA1C", "bb_CRE",
           "smoke_status", "drug.lipid", "drug.htn", "drug.dm")
  miss <- setdiff(req, names(dat))
  if (is.na(sbp_col)) miss <- c(miss, "sbp or sbp_auto")
  if (is.na(dbp_col)) miss <- c(miss, "dbp or dbp_auto")
  if (length(miss)) stop(
    "Missing required CKM construction variables in all.rds: ",
    paste(unique(miss), collapse = ", "),
    ". Re-run the phe.R merge step after confirming the baseline phenotype and blood-biochemistry inputs.",
    call. = FALSE
  )
  dat$ckm_sbp <- safe_num(dat[[sbp_col]])
  dat$ckm_dbp <- safe_num(dat[[dbp_col]])

  # AHA stage-4 domains using the exact UKB phenotype names requested by the
  # user: CHD=cvd_cad/mi/(optional ihd), HF=cvd_hfail, stroke=stroke,
  # PAD=cvd_pad, and AF=cvd_afib. Only fod_icd10_* dates enter this endpoint.
  ckm4_exact <- timed_eval("build_ckm: deriving exact ICD-10 stage-4 dates", make_exact_ckm4_dates(dat))
  for (nm in names(ckm4_exact$domain_dates)) dat[[paste0("date_ckm4_", nm)]] <- ckm4_exact$domain_dates[[nm]]
  dat$date_ckm4_all <- ckm4_exact$date_ckm4_all
  dat$ckm4_first_domain <- ckm4_exact$first_domain
  stage4_source_audit <- ckm4_exact$source_audit
  ckm4_broad <- timed_eval("build_ckm: deriving all-source stage-4 sensitivity dates", make_all_source_ckm4_dates(dat))
  dat$date_ckm4_allsource <- ckm4_broad$date_ckm4_all
  stage4_allsource_audit <- ckm4_broad$source_audit

  # Kidney measures. Direct UACR is preferred; otherwise derive mg/g from urine
  # albumin mg/L and urine creatinine umol/L: UACR = albumin*8840/creatinine.
  uacr_col <- resolve_col(dat, "CKM_UACR_COL", c("uacr", "UACR", "bb_UACR", "urine_uacr", "uacr_mg_g"))
  ua_col <- resolve_col(dat, "CKM_URINE_ALBUMIN_COL", c("bb_UALB", "urine_albumin", "microalbumin_urine", "ualb"))
  uc_col <- resolve_col(dat, "CKM_URINE_CREATININE_COL", c("bb_UCRE", "urine_creatinine", "ucre", "urine_creatinine_umol_l"))
  if (!is.na(uacr_col)) {
    dat$uacr_mg_g <- safe_num(dat[[uacr_col]])
    uacr_source <- uacr_col
  } else if (!is.na(ua_col) && !is.na(uc_col)) {
    dat$uacr_mg_g <- safe_num(dat[[ua_col]]) * 8840 / safe_num(dat[[uc_col]])
    dat$uacr_mg_g[!is.finite(dat$uacr_mg_g) | dat$uacr_mg_g <= 0] <- NA_real_
    uacr_source <- paste0(ua_col, " + ", uc_col)
  } else {
    dat$uacr_mg_g <- NA_real_
    uacr_source <- "unavailable"
  }
  progress_log("build_ckm: UACR source=", uacr_source)
  dat$uacr_assessed <- if ("bb_UACR_assessed" %in% names(dat)) {
    dat$bb_UACR_assessed %in% TRUE
  } else {
    is.finite(dat$uacr_mg_g)
  }
  dat$uacr_category <- if ("bb_UACR_category" %in% names(dat)) {
    as.character(dat$bb_UACR_category)
  } else {
    dplyr::case_when(
      is.finite(dat$uacr_mg_g) & dat$uacr_mg_g < 30 ~ "A1",
      is.finite(dat$uacr_mg_g) & dat$uacr_mg_g < 300 ~ "A2",
      is.finite(dat$uacr_mg_g) ~ "A3",
      TRUE ~ NA_character_
    )
  }
  dat$uacr_lower_mg_g <- if ("bb_UACR_lo" %in% names(dat)) safe_num(dat$bb_UACR_lo) else dat$uacr_mg_g
  dat$uacr_upper_mg_g <- if ("bb_UACR_hi" %in% names(dat)) safe_num(dat$bb_UACR_hi) else dat$uacr_mg_g
  dat$uacr_interval_uncertain <- if ("bb_UACR_interval_uncertain" %in% names(dat)) {
    dat$bb_UACR_interval_uncertain %in% TRUE
  } else {
    rep(FALSE, nrow(dat))
  }

  # CKD-EPI 2021 creatinine equation. Do not truncate eGFR at 15, because that
  # would make kidney failure and G5 impossible to identify.
  t0 <- Sys.time()
  progress_log("build_ckm: calculating kidney/metabolic criteria")
  scr_mgdl <- safe_num(dat$bb_CRE) / 88.4
  kappa <- ifelse(safe_num(dat$sex) == 1, 0.9, 0.7)
  alpha <- ifelse(safe_num(dat$sex) == 1, -0.302, -0.241)
  dat$egfr <- 142 * pmin(scr_mgdl / kappa, 1)^alpha * pmax(scr_mgdl / kappa, 1)^(-1.2) *
    0.9938^safe_num(dat$age) * ifelse(safe_num(dat$sex) == 1, 1, 1.012)
  dat$egfr[!is.finite(dat$egfr) | dat$egfr <= 0] <- NA_real_
  dat$egfr <- pmin(dat$egfr, 200)
  dat$ckd_risk <- ckd_risk_category(dat$egfr, dat$uacr_mg_g)

  hba1c_mmol <- safe_num(dat$bb_HBA1C)
  dat$hba1c_pct <- 0.09148 * hba1c_mmol + 2.152
  glucose <- if ("bb_GLU" %in% names(dat)) safe_num(dat$bb_GLU) else rep(NA_real_, nrow(dat))
  fasting_col <- first_present(dat, c("fasting_hours", "fasting_time", "time_since_last_meal", "fasting"))
  fasting_value <- if (!is.na(fasting_col)) safe_num(dat[[fasting_col]]) else rep(NA_real_, nrow(dat))
  fasting_ok <- is.finite(fasting_value) & fasting_value >= 8
  fasting_unknown <- !is.finite(fasting_value)
  glucose_usable <- fasting_ok | (allow_unknown_fasting_glucose & fasting_unknown)
  # The 2026 AHA guideline uses different thresholds for fasting and
  # nonfasting TG.  Unknown fasting status is conservatively treated as
  # nonfasting rather than applying the lower fasting cutoff to everyone.
  tg_value <- safe_num(dat$bb_TG)
  if (allow_unknown_fasting_tg) {
    tg_threshold <- ifelse(fasting_ok, tg_fasting_cutoff, tg_nonfasting_cutoff)
    tg_usable <- is.finite(tg_value)
    tg_fasting_mode <- if (is.na(fasting_col)) {
      "pragmatic sensitivity: fasting time unavailable; unknown status uses nonfasting threshold"
    } else {
      "pragmatic sensitivity: fasting>=8h uses fasting threshold; <8h/unknown uses nonfasting threshold"
    }
  } else {
    tg_threshold <- rep(tg_fasting_cutoff, nrow(dat))
    tg_usable <- is.finite(tg_value) & fasting_ok
    tg_fasting_mode <- if (is.na(fasting_col)) {
      "strict AHA primary: fasting time unavailable; TG criterion not assigned"
    } else {
      "strict AHA primary: only fasting>=8h TG values are used"
    }
  }
  statin_col <- resolve_col(dat, "CKM_STATIN_COL", c("drug.statin", "statin", "statin_use", "drug.lipid"))
  if (is.na(statin_col)) stop("A statin-treatment column is required for PREVENT; set CKM_STATIN_COL", call. = FALSE)

  ethnicity <- if ("ethnic.c" %in% names(dat)) as.character(dat$ethnic.c) else rep("", nrow(dat))
  asian <- grepl("Asian|Chinese|Indian|Pakistani|Bangladeshi|South Asian", ethnicity, ignore.case = TRUE)
  male <- safe_num(dat$sex) == 1
  bmi_cut <- ifelse(asian, 23, 25)
  waist_cut <- ifelse(male, ifelse(asian, 90, 102), ifelse(asian, 80, 88))
  # Stage 1 permits BMI-defined overweight/obesity OR abdominal adiposity.
  # Metabolic syndrome is different: its adiposity component is elevated waist
  # circumference, not BMI. Keeping these variables separate avoids inflating
  # the stage-2 metabolic-syndrome count.
  dat$waist_high <- safe_num(dat$waist) >= waist_cut
  dat$bmi_high <- safe_num(dat$bmi) >= bmi_cut
  dat$adiposity <- dat$bmi_high | dat$waist_high
  dat$prediabetes <- (hba1c_mmol >= 39 & hba1c_mmol < 48) |
    (glucose_usable & glucose >= 5.6 & glucose < 7.0)
  prior_t2d <- trait_prevalent_simple(dat, "t2dm")
  prior_htn <- trait_prevalent_simple(dat, "cvd_htn")
  dat$T2D_def <- prior_t2d | hba1c_mmol >= 48 | (glucose_usable & glucose >= 7.0) | safe_num(dat$drug.dm) == 1
  dat$HT_def <- prior_htn | safe_num(dat$ckm_sbp) >= 130 | safe_num(dat$ckm_dbp) >= 80 | safe_num(dat$drug.htn) == 1
  dat$hyperTG <- tg_usable & safe_num(dat$bb_TG) >= tg_threshold
  dat$lowHDL <- (male & safe_num(dat$bb_HDL) < 1.036) | (!male & safe_num(dat$bb_HDL) < 1.295)
  dat$metabolic_syndrome <- rowSums(cbind(
    dat$waist_high,
    dat$hyperTG,
    dat$lowHDL,
    dat$HT_def,
    (glucose_usable & glucose >= 5.6) | dat$T2D_def
  ), na.rm = TRUE) >= 3

  existing_ckd_mh <- if ("ckdmh" %in% names(dat)) trait_prevalent_simple(dat, "ckdmh") else rep(FALSE, nrow(dat))
  existing_ckd_ex <- if ("ckdex" %in% names(dat)) trait_prevalent_simple(dat, "ckdex") else rep(FALSE, nrow(dat))
  dat$ckd_moderate_high <- as.character(dat$ckd_risk) %in% c("moderate", "high") | existing_ckd_mh
  dat$ckd_very_high <- as.character(dat$ckd_risk) == "very_high" | existing_ckd_ex
  kidney_failure_traits <- c("kidney_failure", "esrd", "ckd_esrd", "dialysis", "kidney_transplant")
  kidney_failure_flags <- lapply(kidney_failure_traits, function(y) {
    if (y %in% names(dat) || length(candidate_date_cols(dat, y, "all"))) trait_prevalent_simple(dat, y) else rep(FALSE, nrow(dat))
  })
  dat$kidney_failure <- (is.finite(dat$egfr) & dat$egfr < 15) | Reduce(`|`, kidney_failure_flags)
  progress_log("build_ckm: kidney/metabolic criteria done in ", elapsed_since(t0))

  dat$stage1_criterion <- dat$adiposity | dat$prediabetes
  dat$stage2_criterion <- dat$T2D_def | dat$HT_def | dat$hyperTG | dat$metabolic_syndrome | dat$ckd_moderate_high
  dat$ckm_risk_any <- dat$stage1_criterion | dat$stage2_criterion | dat$ckd_very_high
  # date_first_clinical_cvd is the earliest recorded CHD/HF/stroke/PAD/AF date.
  # It is not claimed to be the exact biological onset of CKM stage 4 because
  # the historical timing of the accompanying metabolic/kidney risk state is
  # generally unavailable. date_ckm4 is retained later as a compatibility alias.
  dat$prior_cvd_icd10 <- !is.na(dat$date_ckm4_all) & dat$date_ckm4_all <= dat$date_attend
  dat$prior_cvd_allsource <- !is.na(dat$date_ckm4_allsource) & dat$date_ckm4_allsource <= dat$date_attend
  dat$date_first_clinical_cvd <- if (stage4_source_primary == "all") dat$date_ckm4_allsource else dat$date_ckm4_all
  dat$prior_cvd <- if (stage4_source_primary == "all") dat$prior_cvd_allsource else dat$prior_cvd_icd10

  # Calculate PREVENT for the candidate comparison and sensitivity audit even
  # when SCORE2 is selected for the final Stage-3 assignment.
  prevent_cache <- file.path(indir, "Rdata", "ckm.prevent.rds")
  reuse_prevent <- as_bool(Sys.getenv("CKM_REUSE_PREVENT_CACHE", unset = "TRUE"), TRUE)
  prevent <- NULL
  if (reuse_prevent && file.exists(prevent_cache)) {
    cached <- try(readRDS(prevent_cache), silent = TRUE)
    required_prevent <- c("eid", "prevent_cvd10_base", "prevent_cvd10_best")
    if (!inherits(cached, "try-error") &&
        !length(setdiff(required_prevent, names(cached))) &&
        nrow(cached) == nrow(dat) &&
        identical(as.character(cached$eid), as.character(dat$eid))) {
      prevent <- cached
      progress_log("build_ckm: reusing validated PREVENT cache: ", prevent_cache)
    } else {
      progress_log("build_ckm: PREVENT cache failed validation; recalculating")
    }
  }
  if (is.null(prevent)) {
    prevent <- timed_eval("build_ckm: calculating PREVENT-CVD risk", calculate_prevent(dat, statin_col))
    saveRDS(prevent, prevent_cache)
  }
  idate_cols <- vapply(dat, inherits, logical(1), "IDate")
  for (cc in names(dat)[idate_cols]) {
    dat[[cc]] <- as.Date(as.numeric(unclass(dat[[cc]])), origin = "1970-01-01")
  }
  dat <- left_join(dat, as.data.frame(prevent), by = "eid")
  dat$prevent_cvd10 <- if (prevent_model_primary == "best_available") dat$prevent_cvd10_best else dat$prevent_cvd10_base
  dat$prevent_model <- if (prevent_model_primary == "best_available") dat$prevent_model_best else dat$prevent_model_base
  dat$prevent_problem <- if (prevent_model_primary == "best_available") dat$prevent_problem_best else dat$prevent_problem_base
  if (prevent_model_primary == "best_available") {
    miss_best <- !is.finite(dat$prevent_cvd10)
    dat$prevent_cvd10[miss_best] <- dat$prevent_cvd10_base[miss_best]
    dat$prevent_model[miss_best] <- dat$prevent_model_base[miss_best]
    dat$prevent_problem[miss_best] <- dat$prevent_problem_base[miss_best]
  }

  subascvd_col <- resolve_col(dat, "CKM_SUBCLINICAL_ASCVD_COL", c("subclinical_ascvd", "preclinical_ascvd"))
  prehf_col <- resolve_col(dat, "CKM_PREHF_COL", c("pre_hf", "prehf", "subclinical_hf"))
  cac_col <- resolve_col(dat, "CKM_CAC_COL", c("cac", "CAC", "agatston", "coronary_calcium"))
  abi_col <- resolve_col(dat, "CKM_ABI_COL", c("abi", "ABI", "ankle_brachial_index"))

  sub_ascvd <- rep(FALSE, nrow(dat))
  sub_ascvd_observed <- rep(FALSE, nrow(dat))
  pre_hf <- rep(FALSE, nrow(dat))
  pre_hf_observed <- rep(FALSE, nrow(dat))
  if (!is.na(subascvd_col)) {
    x <- safe_num(dat[[subascvd_col]])
    sub_ascvd <- sub_ascvd | x == 1
    sub_ascvd_observed <- sub_ascvd_observed | is.finite(x)
  }
  if (!is.na(cac_col)) {
    x <- safe_num(dat[[cac_col]])
    sub_ascvd <- sub_ascvd | x >= 100
    sub_ascvd_observed <- sub_ascvd_observed | is.finite(x)
  }
  if (!is.na(abi_col)) {
    x <- safe_num(dat[[abi_col]])
    sub_ascvd <- sub_ascvd | x < 0.90
    sub_ascvd_observed <- sub_ascvd_observed | is.finite(x)
  }
  if (!is.na(prehf_col)) {
    x <- safe_num(dat[[prehf_col]])
    pre_hf <- x == 1
    pre_hf_observed <- is.finite(x)
  }
  sub_ascvd[is.na(sub_ascvd)] <- FALSE
  pre_hf[is.na(pre_hf)] <- FALSE

  # Preserve assessed-negative as 0 and unassessed as NA. phe.v26.R creates
  # these tri-state GP variables from actual CAC/ABI/LVEF/structural records.
  dat$subclinical_ascvd <- ifelse(
    sub_ascvd_observed, as.integer(sub_ascvd), NA_integer_
  )
  dat$subclinical_ascvd_observed <- sub_ascvd_observed
  dat$pre_hf <- ifelse(pre_hf_observed, as.integer(pre_hf), NA_integer_)
  dat$pre_hf_observed <- pre_hf_observed
  if (!"date_subclinical_ascvd" %in% names(dat)) {
    dat$date_subclinical_ascvd <- as.Date(rep(NA, nrow(dat)))
  } else {
    dat$date_subclinical_ascvd <- as_date2(dat$date_subclinical_ascvd)
  }
  if (!"date_pre_hf" %in% names(dat)) {
    dat$date_pre_hf <- as.Date(rep(NA, nrow(dat)))
  } else {
    dat$date_pre_hf <- as_date2(dat$date_pre_hf)
  }
  if (!"subclinical_ascvd_source" %in% names(dat)) {
    dat$subclinical_ascvd_source <- rep(NA_character_, nrow(dat))
  }
  if (!"pre_hf_source" %in% names(dat)) {
    dat$pre_hf_source <- rep(NA_character_, nrow(dat))
  }

  dat$stage2_metabolic_component <- dat$T2D_def | dat$HT_def | dat$hyperTG | dat$metabolic_syndrome
  dat$stage2_kidney_component <- dat$ckd_moderate_high
  dat$ckm_stage2_subtype <- dplyr::case_when(
    dat$stage2_metabolic_component & dat$stage2_kidney_component ~ "mixed metabolic+kidney",
    dat$stage2_metabolic_component ~ "metabolic",
    dat$stage2_kidney_component ~ "kidney",
    TRUE ~ NA_character_
  )


  calculate_score2_stage3 <- function(dat) {
    if (!requireNamespace("RiskScorescvd", quietly = TRUE)) stop("RiskScorescvd required")
    cap <- function(x, lower, upper) pmin(pmax(x, lower), upper)
    d <- data.table::as.data.table(data.table::copy(dat))

    # Keep the source names age and sex throughout this calculation.
    d$age <- cap(d$age, 40, 89)
    d$sex <- ifelse(d$sex == 1, "male", "female")
    d$smoker <- as.integer(d$smoke_status == 2)
    d$diabetes <- as.integer(d$T2D_def == TRUE)
    d$systolic.bp <- cap(d$ckm_sbp, 90, 200)
    d$total.chol <- cap(d$bb_TC, 2, 10)
    d$total.hdl <- cap(d$bb_HDL, 0.5, 3)
    d$eGFR <- cap(d$egfr, 15, 120)

    s <- d[prior_cvd == FALSE, .(eid, age, sex, smoker, systolic.bp, diabetes, total.chol, total.hdl)]
    score2_one <- function(age, sex, smoker, systolic.bp, diabetes, total.chol, total.hdl) {
      values <- c(age, smoker, systolic.bp, diabetes, total.chol, total.hdl)
      if (!all(is.finite(values)) || is.na(sex)) return(NA_character_)
      tryCatch(
        as.character(RiskScorescvd::SCORE2(
          "Low", age, sex, smoker, systolic.bp, diabetes, total.chol, total.hdl, TRUE
        )),
        error = function(e) NA_character_
      )
    }
    s$SCORE2_strat <- mapply(
      score2_one, s$age, s$sex, s$smoker, s$systolic.bp, s$diabetes, s$total.chol, s$total.hdl,
      USE.NAMES = FALSE
    )
    s2 <- s[, .(eid, SCORE2_strat)]

    dm <- d[
      prior_cvd == FALSE & diabetes == 1 & age >= 40 & age <= 69,
      .(eid, age, sex, smoker, systolic.bp, total.chol, total.hdl,
        diabetes = 1L, diabetes.age = age, HbA1c = bb_HBA1C, eGFR)
    ]
    score2_dm_one <- function(age, sex, smoker, systolic.bp, total.chol, total.hdl,
                              diabetes.age, HbA1c, eGFR) {
      values <- c(age, smoker, systolic.bp, total.chol, total.hdl, diabetes.age, HbA1c, eGFR)
      if (!all(is.finite(values)) || is.na(sex)) return(NA_real_)
      tryCatch(
        as.numeric(RiskScorescvd::SCORE2_Diabetes(
          "Low", age, sex, smoker, systolic.bp, total.chol, total.hdl,
          1L, diabetes.age, HbA1c, eGFR, FALSE
        )),
        error = function(e) NA_real_
      )
    }
    sdm <- if (nrow(dm)) {
      dm$SCORE2_DM <- mapply(
        score2_dm_one, dm$age, dm$sex, dm$smoker, dm$systolic.bp,
        dm$total.chol, dm$total.hdl, dm$diabetes.age, dm$HbA1c, dm$eGFR,
        USE.NAMES = FALSE
      )
      dm[, .(eid, SCORE2_DM)]
    } else {
      data.table::data.table(eid = d$eid[0], SCORE2_DM = numeric())
    }

    ck <- d[
      prior_cvd == FALSE,
      .(eid, age, sex, smoker, systolic.bp, total.chol, total.hdl, diabetes, eGFR, ACR = NA_real_)
    ]
    score2_ckd_one <- function(age, sex, smoker, systolic.bp, diabetes, total.chol, total.hdl, eGFR) {
      values <- c(age, smoker, systolic.bp, diabetes, total.chol, total.hdl, eGFR)
      if (!all(is.finite(values)) || is.na(sex)) return(NA_character_)
      tryCatch(
        as.character(RiskScorescvd::SCORE2_CKD(
          "Low", age, sex, smoker, systolic.bp, diabetes, total.chol, total.hdl,
          eGFR, NA_real_, NA_real_, "ACR", TRUE
        )),
        error = function(e) NA_character_
      )
    }
    ck$SCORE2_CKD_strat <- mapply(
      score2_ckd_one, ck$age, ck$sex, ck$smoker, ck$systolic.bp, ck$diabetes, ck$total.chol, ck$total.hdl, ck$eGFR,
      USE.NAMES = FALSE
    )
    sck <- ck[, .(eid, SCORE2_CKD_strat)]

    out <- Reduce(
      function(x, y) merge(x, y, by = "eid", all.x = TRUE),
      list(data.table::data.table(eid = d$eid), s2, sdm, sck)
    )
    out$stage3_score2_component <-
      (!is.na(out$SCORE2_strat) & out$SCORE2_strat == "High risk") |
      (!is.na(out$SCORE2_DM) & out$SCORE2_DM / 100 >= 0.10) |
      (!is.na(out$SCORE2_CKD_strat) & out$SCORE2_CKD_strat == "High risk")
    out[, .(eid, stage3_score2_component)]
  }
  dat <- merge(dat, calculate_score2_stage3(dat), by = "eid", all.x = TRUE)
  dat$prevent_cvd10_best_effective <- ifelse(is.finite(dat$prevent_cvd10_best), dat$prevent_cvd10_best, dat$prevent_cvd10_base)
  dat$stage3_prevent_component_base <- is.finite(dat$prevent_cvd10_base) & dat$prevent_cvd10_base >= prevent_cutoff
  dat$stage3_prevent_component_best <- is.finite(dat$prevent_cvd10_best_effective) & dat$prevent_cvd10_best_effective >= prevent_cutoff
  dat$stage3_prevent_component <- if (prevent_model_primary == "best_available") dat$stage3_prevent_component_best else dat$stage3_prevent_component_base
  dat$stage3_ckd_component <- dat$ckd_very_high
  dat$stage3_subclinical_ascvd_component <- sub_ascvd & dat$ckm_risk_any
  dat$stage3_prehf_component <- pre_hf & dat$ckm_risk_any
  make_stage3_criterion <- function(risk_component) {
    dat$stage3_ckd_component | risk_component |
      dat$stage3_subclinical_ascvd_component | dat$stage3_prehf_component
  }
  make_stage3_subtype <- function(risk_component, risk_label) {
    n_components <- rowSums(cbind(
      risk_component, dat$stage3_ckd_component,
      dat$stage3_subclinical_ascvd_component, dat$stage3_prehf_component
    ), na.rm = TRUE)
    dplyr::case_when(
      n_components > 1 ~ "mixed",
      risk_component ~ risk_label,
      dat$stage3_ckd_component ~ "very-high-risk CKD",
      dat$stage3_subclinical_ascvd_component ~ "subclinical ASCVD",
      dat$stage3_prehf_component ~ "pre-HF",
      TRUE ~ NA_character_
    )
  }
  dat$stage3_criterion_SCORE2 <- make_stage3_criterion(dat$stage3_score2_component)
  dat$stage3_criterion_PREVENT <- make_stage3_criterion(dat$stage3_prevent_component)
  dat$ckm_stage3_subtype_SCORE2 <- make_stage3_subtype(dat$stage3_score2_component, "SCORE2_high_risk")
  dat$ckm_stage3_subtype_PREVENT <- make_stage3_subtype(dat$stage3_prevent_component, "PREVENT>=cutoff")

  dat$stage4_criterion_icd10 <- dat$prior_cvd_icd10 & (!stage4_require_risk | dat$ckm_risk_any)
  dat$stage4_criterion_allsource <- dat$prior_cvd_allsource & (!stage4_require_risk | dat$ckm_risk_any)
  dat$stage4_criterion <- if (stage4_source_primary == "all") dat$stage4_criterion_allsource else dat$stage4_criterion_icd10
  dat$ckm_stage4_substage <- dplyr::case_when(
    dat$stage4_criterion & dat$kidney_failure ~ "4b",
    dat$stage4_criterion ~ "4a",
    TRUE ~ NA_character_
  )

  # Higher stages take priority. Stage 0 is the available-data negative state.
  # A separate certainty field identifies participants whose missing UACR,
  # selected risk-method, subclinical-ASCVD, or pre-HF information could permit under-staging.
  t0 <- Sys.time()
  progress_log("build_ckm: assigning CKM stages")
  core_complete <- complete.cases(dat[, intersect(c(
    "bmi", "waist", "ckm_sbp", "ckm_dbp", "bb_HDL", "bb_TG", "bb_HBA1C", "bb_CRE"
  ), names(dat)), drop = FALSE]) & is.finite(safe_num(dat$bb_TG))
  assign_ckm_stage <- function(stage3_criterion) dplyr::case_when(
    dat$stage4_criterion ~ 4L,
    !dat$prior_cvd & stage3_criterion ~ 3L,
    !dat$prior_cvd & dat$stage2_criterion ~ 2L,
    !dat$prior_cvd & dat$stage1_criterion ~ 1L,
    !dat$prior_cvd & core_complete ~ 0L,
    TRUE ~ NA_integer_
  )
  assign_stage_certainty <- function(stage, risk_assessed) {
    stage3_negative_assessed <- dat$uacr_assessed & risk_assessed & sub_ascvd_observed & pre_hf_observed
    dplyr::case_when(
      stage %in% c(3L, 4L) ~ "criterion-positive",
      stage %in% c(1L, 2L) & stage3_negative_assessed ~ "criterion-positive; higher-stage assessed",
      stage %in% c(1L, 2L) ~ "criterion-positive; possible under-staging",
      stage == 0L & stage3_negative_assessed ~ "complete available-data negative",
      stage == 0L ~ "available-data negative; possible under-staging",
      TRUE ~ "insufficient baseline information"
    )
  }
  dat$ckm.stage.SCORE2 <- assign_ckm_stage(dat$stage3_criterion_SCORE2)
  dat$ckm.stage.PREVENT <- assign_ckm_stage(dat$stage3_criterion_PREVENT)
  dat$ckm_stage_certainty_SCORE2 <- assign_stage_certainty(dat$ckm.stage.SCORE2, !is.na(dat$stage3_score2_component))
  dat$ckm_stage_certainty_PREVENT <- assign_stage_certainty(dat$ckm.stage.PREVENT, is.finite(dat$prevent_cvd10))

  # The selected aliases are used only for this invocation's Fig1/audit. Both
  # complete assignments are persisted in ckm.rds for downstream selection.
  dat$ckm.stage <- dat[[paste0("ckm.stage.", stage3_method)]]
  dat$stage3_component <- dat[[paste0("stage3_", tolower(stage3_method), "_component")]]
  dat$stage3_criterion <- dat[[paste0("stage3_criterion_", stage3_method)]]
  dat$ckm_stage3_subtype <- dat[[paste0("ckm_stage3_subtype_", stage3_method)]]
  dat$ckm_stage_certainty <- dat[[paste0("ckm_stage_certainty_", stage3_method)]]
  selected_risk_assessed <- if (stage3_method == "SCORE2") {
    !is.na(dat$stage3_score2_component)
  } else {
    is.finite(dat$prevent_cvd10)
  }
  stage3_negative_assessed <- dat$uacr_assessed & selected_risk_assessed &
    sub_ascvd_observed & pre_hf_observed
  progress_log("build_ckm: stage assignment done in ", elapsed_since(t0))

  # date_ckm4 is the compact event-date field used by downstream analyses.
  dat$date_ckm4 <- dat$date_first_clinical_cvd
  dat$date_ckm_death <- ifelse(
    !is.na(dat$date_death) & !is.na(dat$date_first_clinical_cvd) &
      dat$date_death > dat$date_first_clinical_cvd,
    as.character(dat$date_death), NA_character_
  )
  dat$date_ckm_death <- as.Date(dat$date_ckm_death)

# Primary AHA CKM stage is retained unchanged.
# GP medication/primary-care evidence is an additional phenotype-certainty layer.

gp_med_file <- file.path(indir, "Rdata", "gp.ckm.med.rds")
if (file.exists(gp_med_file)) {
  gp_med <- readRDS(gp_med_file)

  dat <- dat %>%
    left_join(gp_med, by = "eid") %>%
    mutate(
      ckm_gp_evidence_score =
        as.integer(
          gp_htn_med %in% TRUE |
          gp_dm_med %in% TRUE |
          gp_lipid_med %in% TRUE
        ),

      ckm_stage_certainty_enhanced =
        case_when(
          ckm_gp_evidence_score > 0 ~ "gp_medication_supported",
          TRUE ~ ckm_stage_certainty
        )
    )
} else {
  dat$ckm_gp_evidence_score <- 0L
  dat$ckm_stage_certainty_enhanced <- dat$ckm_stage_certainty
}

  compact <- dat %>% transmute(
    eid = as.character(eid),
    ckm_definition = "2026_AHA",
    ckm.stage.SCORE2 = as.integer(ckm.stage.SCORE2),
    ckm.stage.PREVENT = as.integer(ckm.stage.PREVENT),
    date_first_clinical_cvd = as_date2(date_first_clinical_cvd),
    date_ckm4 = as_date2(date_ckm4),
    date_ckm4_icd10 = as_date2(date_ckm4_all),
    date_ckm4_allsource = as_date2(date_ckm4_allsource),
    date_ckm_death = as_date2(date_ckm_death),
    ckm_stage_certainty_SCORE2 = as.character(ckm_stage_certainty_SCORE2),
    ckm_stage_certainty_PREVENT = as.character(ckm_stage_certainty_PREVENT),
    ckm_gp_evidence_score = as.integer(ckm_gp_evidence_score),
    gp_htn_med = if ("gp_htn_med" %in% names(dat)) as.logical(gp_htn_med) else FALSE,
    gp_dm_med = if ("gp_dm_med" %in% names(dat)) as.logical(gp_dm_med) else FALSE,
    gp_lipid_med = if ("gp_lipid_med" %in% names(dat)) as.logical(gp_lipid_med) else FALSE,
    ckm_stage2_subtype = as.character(ckm_stage2_subtype),
    ckm_stage3_subtype_SCORE2 = as.character(ckm_stage3_subtype_SCORE2),
    ckm_stage3_subtype_PREVENT = as.character(ckm_stage3_subtype_PREVENT),
    ckm_stage4_substage = as.character(ckm_stage4_substage),
    egfr = safe_num(egfr),
    uacr_mg_g = safe_num(uacr_mg_g),
    uacr_category = as.character(uacr_category),
    uacr_assessed = as.logical(uacr_assessed),
    uacr_lower_mg_g = safe_num(uacr_lower_mg_g),
    uacr_upper_mg_g = safe_num(uacr_upper_mg_g),
    uacr_interval_uncertain = as.logical(uacr_interval_uncertain),
    ckd_risk = as.character(ckd_risk),
    subclinical_ascvd = as.integer(subclinical_ascvd),
    subclinical_ascvd_observed = as.logical(subclinical_ascvd_observed),
    date_subclinical_ascvd = as_date2(date_subclinical_ascvd),
    subclinical_ascvd_source = as.character(subclinical_ascvd_source),
    pre_hf = as.integer(pre_hf),
    pre_hf_observed = as.logical(pre_hf_observed),
    date_pre_hf = as_date2(date_pre_hf),
    pre_hf_source = as.character(pre_hf_source),
    stage3_prevent_component = as.logical(stage3_prevent_component),
    stage3_score2_component = as.logical(stage3_score2_component),
    stage3_ckd_component = as.logical(stage3_ckd_component),
    stage3_subclinical_ascvd_component = as.logical(stage3_subclinical_ascvd_component),
    stage3_prehf_component = as.logical(stage3_prehf_component)
  )
  dir.create(dirname(ckm_file), recursive = TRUE, showWarnings = FALSE)
  timed_eval("build_ckm: saving compact ckm.rds", saveRDS(compact, ckm_file))
  progress_log("Saved compact ckm.rds: ", ckm_file)
  print(list(
    SCORE2 = table(compact$ckm.stage.SCORE2, useNA = "always"),
    PREVENT = table(compact$ckm.stage.PREVENT, useNA = "always")
  ))

  # Imaging variables are prospective landmark/interval-censored outcomes.
  # They are saved separately and never back-filled into recruitment-baseline
  # CKM stage, protecting baseline omic analyses from future-information leakage.
  img_names <- c(
    "date_img_i2", "date_img_i3",
    "pre_hf_img_hfref_i2", "pre_hf_img_hfref_i3",
    "pre_hf_img_broad_sensitivity_i2", "pre_hf_img_broad_sensitivity_i3",
    "carotid_imt_mean_um_i2", "carotid_imt_mean_um_i3",
    "carotid_imt_assessed_i2", "carotid_imt_assessed_i3",
    "pre_hf_img_incident_i3", "date_pre_hf_img_incident_i3"
  )
  if (any(img_names %in% names(dat))) {
    n_img <- nrow(dat)
    img_num <- function(nm) if (nm %in% names(dat)) safe_num(dat[[nm]]) else rep(NA_real_, n_img)
    img_int <- function(nm) if (nm %in% names(dat)) as.integer(safe_num(dat[[nm]])) else rep(NA_integer_, n_img)
    img_log <- function(nm) if (nm %in% names(dat)) dat[[nm]] %in% TRUE else rep(FALSE, n_img)
    img_date <- function(nm) if (nm %in% names(dat)) as_date2(dat[[nm]]) else as.Date(rep(NA, n_img))
    img_compact <- data.frame(
      eid = as.character(dat$eid),
      baseline_ckm_stage = as.integer(dat$ckm.stage),
      date_attend = as_date2(dat$date_attend),
      date_img_i2 = img_date("date_img_i2"),
      date_img_i3 = img_date("date_img_i3"),
      pre_hf_img_hfref_i2 = img_int("pre_hf_img_hfref_i2"),
      pre_hf_img_hfref_i3 = img_int("pre_hf_img_hfref_i3"),
      pre_hf_img_broad_sensitivity_i2 = img_int("pre_hf_img_broad_sensitivity_i2"),
      pre_hf_img_broad_sensitivity_i3 = img_int("pre_hf_img_broad_sensitivity_i3"),
      carotid_imt_mean_um_i2 = img_num("carotid_imt_mean_um_i2"),
      carotid_imt_mean_um_i3 = img_num("carotid_imt_mean_um_i3"),
      carotid_imt_assessed_i2 = img_log("carotid_imt_assessed_i2"),
      carotid_imt_assessed_i3 = img_log("carotid_imt_assessed_i3"),
      pre_hf_img_incident_i3 = img_int("pre_hf_img_incident_i3"),
      date_pre_hf_img_incident_i3 = img_date("date_pre_hf_img_incident_i3"),
      stringsAsFactors = FALSE
    )
    img_compact$stage3_prehf_interval_i2 <-
      img_compact$baseline_ckm_stage %in% 0:2 &
      img_compact$pre_hf_img_hfref_i2 == 1L
    img_compact$interval_lower_i2 <- img_compact$date_attend
    img_compact$interval_upper_i2 <- img_compact$date_img_i2
    img_compact$interval_lower_i2[!img_compact$stage3_prehf_interval_i2] <- as.Date(NA)
    img_compact$interval_upper_i2[!img_compact$stage3_prehf_interval_i2] <- as.Date(NA)
    img_ckm_file <- file.path(dirname(ckm_file), "img.ckm.rds")
    saveRDS(img_compact, img_ckm_file)
    progress_log("Saved imaging-landmark CKM object: ", img_ckm_file)
  }

  t0 <- Sys.time()
  progress_log("build_ckm: preparing audit tables")
  audit <- dat %>% transmute(
    eid, date_attend, ckm_stage = ckm.stage,
    ckm_stage_certainty, ckm_stage2_subtype, ckm_stage3_subtype,
    date_first_clinical_cvd, date_ckm4, date_ckm4_allsource,
    ckm4_first_domain, prior_cvd, prior_cvd_icd10, prior_cvd_allsource,
    stage4_criterion, stage4_criterion_icd10, stage4_criterion_allsource, ckm_stage4_substage,
    ckm_risk_any, stage1_criterion, stage2_criterion, stage3_criterion,
    bmi_high, waist_high, adiposity, prediabetes,
    T2D_def, HT_def, hyperTG, lowHDL, metabolic_syndrome,
    stage2_metabolic_component, stage2_kidney_component,
    stage3_prevent_component, stage3_prevent_component_base, stage3_prevent_component_best, stage3_ckd_component,
    stage3_subclinical_ascvd_component, stage3_prehf_component,
    ckm_sbp, ckm_dbp,
    egfr, uacr_mg_g, uacr_category, uacr_assessed,
    uacr_lower_mg_g, uacr_upper_mg_g, uacr_interval_uncertain,
    ckd_risk, ckd_moderate_high, ckd_very_high, kidney_failure,
    prevent_cvd10, prevent_model, prevent_problem,
    prevent_cvd10_base, prevent_model_base, prevent_problem_base,
    prevent_cvd10_best, prevent_cvd10_best_effective, prevent_model_best, prevent_problem_best,
    tg_threshold,
    subclinical_ascvd, subclinical_ascvd_observed,
    date_subclinical_ascvd, subclinical_ascvd_source,
    pre_hf, pre_hf_observed, date_pre_hf, pre_hf_source,
    core_complete, stage3_negative_assessed,
    date_ckm4_chd, date_ckm4_hf, date_ckm4_stroke, date_ckm4_pad, date_ckm4_af
  )
  summary <- audit %>% group_by(ckm_stage) %>% summarise(
    N = n(),
    prior_cvd = sum(prior_cvd, na.rm = TRUE),
    future_clinical_cvd = sum(!is.na(date_first_clinical_cvd) & date_first_clinical_cvd > date_attend, na.rm = TRUE),
    possible_under_staging = sum(grepl("possible under-staging", ckm_stage_certainty), na.rm = TRUE),
    median_prevent = ifelse(all(is.na(prevent_cvd10)), NA_real_, median(prevent_cvd10, na.rm = TRUE)),
    .groups = "drop"
  )

  stage2_subtypes <- dat %>%
    transmute(
      final_stage = ckm.stage,
      component_count = rowSums(cbind(T2D_def, HT_def, hyperTG, metabolic_syndrome, ckd_moderate_high), na.rm = TRUE),
      dominant_component = case_when(
        ckd_moderate_high ~ "CKD",
        T2D_def ~ "type 2 diabetes",
        metabolic_syndrome ~ "metabolic syndrome",
        hyperTG ~ "hypertriglyceridemia",
        HT_def ~ "hypertension",
        TRUE ~ NA_character_
      )
    ) %>%
    filter(component_count > 0) %>%
    count(final_stage, component_count, dominant_component, name = "N") %>%
    arrange(final_stage, component_count, desc(N))

  stage2_overlap <- dat %>%
    transmute(
      T2D = T2D_def, hypertension = HT_def, hyperTG,
      metabolic_syndrome, CKD_moderate_high = ckd_moderate_high,
      final_stage = ckm.stage
    ) %>%
    filter(T2D | hypertension | hyperTG | metabolic_syndrome | CKD_moderate_high) %>%
    count(final_stage, T2D, hypertension, hyperTG, metabolic_syndrome, CKD_moderate_high, name = "N") %>%
    arrange(final_stage, desc(N))

  stage3_overlap <- dat %>%
    transmute(
      SCORE2_high = stage3_score2_component,
      PREVENT_high = stage3_prevent_component,
      CKD_very_high = stage3_ckd_component,
      subclinical_ASCVD = stage3_subclinical_ascvd_component,
      pre_HF = stage3_prehf_component,
      final_stage = ckm.stage
    ) %>%
    filter(SCORE2_high | PREVENT_high | CKD_very_high | subclinical_ASCVD | pre_HF) %>%
    count(final_stage, SCORE2_high, PREVENT_high, CKD_very_high, subclinical_ASCVD, pre_HF, name = "N") %>%
    arrange(final_stage, desc(N))

  certainty_summary <- dat %>%
    transmute(ckm_stage = ckm.stage, ckm_stage_certainty) %>%
    count(ckm_stage, ckm_stage_certainty, name = "N") %>%
    arrange(ckm_stage, desc(N))

  stage4_sensitivity <- data.table(
    definition = c("all available dated sources", "exact ICD-10 hospital definition"),
    is_primary = c(stage4_source_primary == "all", stage4_source_primary == "icd10"),
    any_prebaseline_clinical_cvd = c(sum(dat$prior_cvd_allsource, na.rm = TRUE), sum(dat$prior_cvd_icd10, na.rm = TRUE)),
    stage4_with_required_ckm_risk = c(sum(dat$stage4_criterion_allsource, na.rm = TRUE), sum(dat$stage4_criterion_icd10, na.rm = TRUE)),
    difference_vs_selected_primary = c(
      sum(dat$stage4_criterion_allsource, na.rm = TRUE) - sum(dat$stage4_criterion, na.rm = TRUE),
      sum(dat$stage4_criterion_icd10, na.rm = TRUE) - sum(dat$stage4_criterion, na.rm = TRUE)
    )
  )

  stage3_ascertainment <- data.table(
    component = c(
      "SCORE2", "PREVENT base", "PREVENT best available",
      "complete KDIGO eGFR+UACR", "eGFR-only G4/G5",
      "subclinical ASCVD", "pre-HF", "UACR measured"
    ),
    assessed_N = c(
      sum(!is.na(dat$stage3_score2_component)),
      sum(is.finite(dat$prevent_cvd10_base)), sum(is.finite(dat$prevent_cvd10_best)),
      sum(is.finite(dat$egfr) & dat$uacr_assessed), sum(is.finite(dat$egfr)),
      sum(sub_ascvd_observed), sum(pre_hf_observed), sum(dat$uacr_assessed)
    ),
    positive_N = c(
      sum(dat$stage3_score2_component, na.rm = TRUE),
      sum(dat$stage3_prevent_component_base, na.rm = TRUE), sum(dat$stage3_prevent_component_best, na.rm = TRUE),
      sum(dat$stage3_ckd_component & is.finite(dat$egfr) & dat$uacr_assessed, na.rm = TRUE),
      sum(is.finite(dat$egfr) & dat$egfr < 30, na.rm = TRUE),
      sum(dat$stage3_subclinical_ascvd_component, na.rm = TRUE),
      sum(dat$stage3_prehf_component, na.rm = TRUE),
      sum(dat$uacr_assessed & is.finite(dat$uacr_mg_g) & dat$uacr_mg_g >= 30)
    ),
    note = c(
      "SCORE2/SCORE2-Diabetes/SCORE2-CKD candidate", "uniform PREVENT base equation", "auto-selected optional-model sensitivity",
      "complete KDIGO heat-map ascertainment requires both eGFR and UACR",
      "eGFR alone can identify G4/G5 very-high-risk CKD but misses albuminuria-defined risk",
      "requires supplied imaging/ABI/subclinical flag",
      "requires supplied pre-HF flag",
      "missing UACR under-classifies stages 2 and 3"
    )
  )
  stage3_definition_sensitivity <- data.table(
    definition = c(
      "SCORE2 + available stage-3 components",
      "PREVENT base + available stage-3 components",
      "PREVENT best available + available stage-3 components"
    ),
    stage3_N = c(
      sum(!dat$prior_cvd & (dat$stage3_score2_component | dat$stage3_ckd_component | dat$stage3_subclinical_ascvd_component | dat$stage3_prehf_component), na.rm = TRUE),
      sum(!dat$prior_cvd & (dat$stage3_prevent_component_base | dat$stage3_ckd_component | dat$stage3_subclinical_ascvd_component | dat$stage3_prehf_component), na.rm = TRUE),
      sum(!dat$prior_cvd & (dat$stage3_prevent_component_best | dat$stage3_ckd_component | dat$stage3_subclinical_ascvd_component | dat$stage3_prehf_component), na.rm = TRUE)
    )
  )

  # Shared construction-audit figure.  This is intentionally written to the
  # CKM root rather than a prot/met subdirectory because it does not depend on
  # the omic layer.
  source_audit <- data.frame(
    item = c("AHA framework", "primary stage-4 date-source rule", "stage-4 sensitivity rule", "ICD-10 definition file",
             "blood-biochemistry source", "metabolic-syndrome adiposity component",
             "systolic BP source", "diastolic BP source",
             "PREVENT stage-3 cutoff", "PREVENT primary model", "UACR source", "statin source",
             "fasting-time source", "triglyceride variable", "TG fasting cutoff mmol/L", "TG nonfasting cutoff mmol/L", "TG fasting mode", "unknown-fasting glucose allowed", "unknown-fasting TG option retained",
             "stage-4 risk-factor sensitivity restriction", "subclinical ASCVD column", "pre-HF column", "CAC column", "ABI column",
             "stage-date semantic warning"),
    value = c("2026 AHA/ACC/ADA/ASN implementation",
              ifelse(stage4_source_primary == "all", "all available dated phenotype sources from t2e-compatible FOD variables", "exact ICD-10 hospital FODs"),
              ifelse(stage4_source_primary == "all", "exact ICD-10 hospital FODs", "all available dated phenotype sources"),
              icd10_list_file,
              "all.rds", "elevated waist circumference only; BMI remains a stage-1 adiposity criterion",
              sbp_col, dbp_col,
              prevent_cutoff, prevent_model_primary, uacr_source, statin_col,
              ifelse(!is.na(fasting_col), fasting_col, "unavailable"), "bb_TG (mmol/L)", tg_fasting_cutoff, tg_nonfasting_cutoff, tg_fasting_mode,
              allow_unknown_fasting_glucose, allow_unknown_fasting_tg,
              stage4_require_risk, subascvd_col, prehf_col, cac_col, abi_col,
              "date_first_clinical_cvd/date_ckm4 is earliest recorded CHD/HF/stroke/PAD/AF date, not proven biological stage-4 onset")
  )
  write_tab(summary, "ckm_stage_summary.tsv", auditdir)
  write_tab(audit, "ckm_build_audit.tsv.gz", auditdir)
  write_tab(source_audit, "ckm_build_sources.tsv", auditdir)
  write_tab(stage4_source_audit, "ckm_stage4_source_map.tsv", auditdir)
  write_tab(stage4_allsource_audit, "ckm_stage4_allsource_map.tsv", auditdir)
  write_tab(stage4_sensitivity, "ckm_stage4_sensitivity.tsv", auditdir)
  write_tab(stage2_overlap, "ckm_stage2_criterion_overlap.tsv", auditdir)
  write_tab(stage2_subtypes, "ckm_stage2_subtypes.tsv", auditdir)
  write_tab(stage3_overlap, "ckm_stage3_criterion_overlap.tsv", auditdir)
  write_tab(stage3_ascertainment, "ckm_stage3_ascertainment.tsv", auditdir)
  write_tab(stage3_definition_sensitivity, "ckm_stage3_definition_sensitivity.tsv", auditdir)
  write_tab(certainty_summary, "ckm_stage_certainty.tsv", auditdir)
  write_book(list(
    stage_summary = summary,
    stage_certainty = certainty_summary,
    stage2_overlap = stage2_overlap,
    stage2_subtypes = stage2_subtypes,
    stage3_overlap = stage3_overlap,
    stage3_ascertainment = stage3_ascertainment,
    stage3_definition_sensitivity = stage3_definition_sensitivity,
    stage4_sensitivity = stage4_sensitivity,
    build_sources = source_audit,
    stage4_source_map = stage4_source_audit,
    stage4_allsource_map = stage4_allsource_audit
  ), "ckm_build_audit.xlsx", auditdir)
  progress_log("build_ckm: audit outputs done in ", elapsed_since(t0))

  if (merge_ckm_to_all) {
    t0 <- Sys.time()
    progress_log("build_ckm: merging compact CKM fields into all.rds")
    old <- intersect(c(
      "ckm.stage", "ckm.stage.SCORE2", "ckm.stage.PREVENT", "ckm_definition", "date_first_clinical_cvd", "date_ckm4", "date_ckm4_icd10", "date_ckm4_allsource", "date_ckm_death",
      "ckm_stage_certainty", "ckm_stage_certainty_SCORE2", "ckm_stage_certainty_PREVENT", "ckm_stage2_subtype", "ckm_stage3_subtype", "ckm_stage3_subtype_SCORE2", "ckm_stage3_subtype_PREVENT", "ckm_stage4_substage",
    "egfr", "uacr_mg_g", "uacr_category", "uacr_assessed",
    "uacr_lower_mg_g", "uacr_upper_mg_g", "uacr_interval_uncertain", "ckd_risk",
    "subclinical_ascvd", "subclinical_ascvd_observed", "date_subclinical_ascvd", "subclinical_ascvd_source",
    "pre_hf", "pre_hf_observed", "date_pre_hf", "pre_hf_source",
    "stage3_prevent_component", "stage3_score2_component", "stage3_ckd_component",
    "stage3_subclinical_ascvd_component", "stage3_prehf_component"
    ), names(dat))
    all0 <- dplyr::select(as.data.frame(readRDS(all_file)), -dplyr::any_of(old))
    all0$eid <- as.character(all0$eid)
    all1 <- left_join(all0, compact, by = "eid")
    saveRDS(all1, all_file)
    progress_log("Merged compact CKM variables into all.rds in ", elapsed_since(t0))
  }
  progress_log("build_ckm: completed in ", elapsed_since(t_build))
  invisible(compact)
}
