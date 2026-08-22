#!/usr/bin/env python3
"""CKM build utility: create a v18 R/shell pair from a compatible v16 pair.

The patch is deliberately strict: it aborts when an expected v16 function or
code block cannot be found, rather than silently producing a partly patched
analysis script.
"""
from __future__ import annotations

import argparse
import os
import re
import stat
from pathlib import Path


class PatchError(RuntimeError):
    pass


def find_function_span(text: str, name: str) -> tuple[int, int]:
    m = re.search(rf"(?m)^{re.escape(name)}\s*<-\s*function\s*\(", text)
    if not m:
        raise PatchError(f"Cannot find R function: {name}")
    start = m.start()
    brace = text.find("{", m.end())
    if brace < 0:
        raise PatchError(f"Cannot find opening brace for R function: {name}")

    depth = 0
    quote: str | None = None
    escaped = False
    in_comment = False
    i = brace
    while i < len(text):
        ch = text[i]
        if in_comment:
            if ch == "\n":
                in_comment = False
            i += 1
            continue
        if quote is not None:
            if escaped:
                escaped = False
            elif ch == "\\":
                escaped = True
            elif ch == quote:
                quote = None
            i += 1
            continue
        if ch == "#":
            in_comment = True
        elif ch in ("'", '"', "`"):
            quote = ch
        elif ch == "{":
            depth += 1
        elif ch == "}":
            depth -= 1
            if depth == 0:
                return start, i + 1
        i += 1
    raise PatchError(f"Unbalanced braces while reading R function: {name}")


def replace_function(text: str, name: str, replacement: str) -> str:
    start, end = find_function_span(text, name)
    return text[:start] + replacement.rstrip() + "\n" + text[end:]


def replace_once(text: str, old: str, new: str, label: str) -> str:
    n = text.count(old)
    if n != 1:
        raise PatchError(f"Expected exactly one occurrence for {label}; found {n}")
    return text.replace(old, new, 1)


def regex_replace_once(text: str, pattern: str, repl: str, label: str, flags: int = 0) -> str:
    out, n = re.subn(pattern, repl, text, count=1, flags=flags)
    if n != 1:
        raise PatchError(f"Expected exactly one regex match for {label}; found {n}")
    return out


CALCULATE_PREVENT = r'''calculate_prevent <- function(dat, statin_col) {
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
  rb <- run_variant(base_cols, "base", model = "base")
  # Auto-selection of the best valid optional model is retained only as a
  # sensitivity analysis.
  re <- run_variant(best_cols, "best_available", model = NULL)
  out$prevent_cvd10_base[rb$row_index] <- rb$risk
  out$prevent_model_base[rb$row_index] <- rb$model
  out$prevent_problem_base[rb$row_index] <- rb$problem
  out$prevent_cvd10_best[re$row_index] <- re$risk
  out$prevent_model_best[re$row_index] <- re$model
  out$prevent_problem_best[re$row_index] <- re$problem
  out
}'''


ALL_SOURCE_DATES = r'''make_all_source_ckm4_dates <- function(dat, baseline_col = "date_attend") {
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
}'''


PREPARE_TARGET = r'''prepare_target_data <- function(dat, y = "ckm4") {
  if (!identical(y, "ckm4")) {
    stop("This pipeline has one fixed target: incident clinical CKM stage 4 (ckm4).", call. = FALSE)
  }
  d <- dat
  event_date <- if ("date_first_clinical_cvd" %in% names(d)) {
    as_date2(d$date_first_clinical_cvd)
  } else {
    as_date2(d$date_ckm4)
  }
  baseline <- as_date2(d$date_attend)
  lost <- if ("date_lost" %in% names(d)) as_date2(d$date_lost) else as.Date(rep(NA, nrow(d)))
  death <- if ("date_death" %in% names(d)) as_date2(d$date_death) else as.Date(rep(NA, nrow(d)))
  follow <- rep(date_follow_end_global, nrow(d))
  censor <- pmin_date2(data.frame(lost = lost, death = death, follow = follow))
  censor[is.na(censor)] <- date_follow_end_global

  # Match 0phe.f.R::t2e(): an event recorded on the assessment date is not
  # assigned a positive or negative b2e. It remains a baseline stage-4 record
  # for categorical staging, but it is excluded from event-proximity YY curves.
  same_day_stage4 <- !is.na(event_date) & !is.na(baseline) & event_date == baseline
  baseline_stage4 <- !is.na(event_date) & !is.na(baseline) & event_date <= baseline
  yang <- !is.na(event_date) & !is.na(baseline) & event_date < baseline
  event_observed <- !is.na(event_date) & !is.na(baseline) & event_date > baseline & event_date <= censor

  stop_date <- censor
  stop_date[event_observed] <- event_date[event_observed]
  time <- as.numeric(stop_date - baseline) / 365.25
  b2e <- as.numeric(event_date - baseline) / 365.25
  b2e[same_day_stage4] <- NA_real_
  yin_eligible <- !baseline_stage4 & d$ckm.stage %in% 0:3 & is.finite(time) & time > 0

  d$target_date <- event_date
  d$censor_date <- censor
  d$b2e <- b2e
  d$baseline_stage4 <- baseline_stage4
  d$same_day_stage4 <- same_day_stage4
  d$yang <- yang
  d$yin <- yin_eligible
  d$event <- as.integer(event_observed & yin_eligible)
  d$time <- time
  d$yin_control <- yin_eligible & d$event == 0
  d
}'''


MARKER_MAP = r'''make_marker_map <- function(dat, features, meta, clock_obj) {
  t_step <- Sys.time()
  position_map <- get_clock_positions(clock_obj)
  covars <- get_clock_covars(dat, "le8_behavior")
  med_covars <- safe_covars(dat, unique(c(covars, "drug.lipid", "drug.htn", "drug.dm")))
  scan_features <- if (scan_all_biom) features else resolve_display_biomarkers(features, meta, n = min(20, length(features)))
  progress_log("Marker scan features: ", fmt_int(length(scan_features)))
  res <- rbindlist(map_features(scan_features, function(p) scan_marker_one(dat, p, position_map, covars, med_covars),
                                label = "marker_map scan"), fill = TRUE)
  res[, biomarker := as.character(biomarker)]
  meta <- as.data.table(copy(meta))
  meta[, biomarker := as.character(biomarker)]
  res <- merge(res, meta, by = "biomarker", all.x = TRUE, sort = FALSE)
  res[, axis_version := "risk_age_v18"]
  res[, `:=`(
    stage_q = p.adjust(stage_p, method = "BH"),
    yin_cox_q = p.adjust(yin_cox_p, method = "BH"),
    yang_proximity_q = p.adjust(yang_proximity_p, method = "BH"),
    prev_vs_yin_q = p.adjust(prev_vs_yin_p, method = "BH"),
    yang_proximity_q_med = p.adjust(yang_proximity_p_med, method = "BH"),
    prev_vs_yin_q_med = p.adjust(prev_vs_yin_p_med, method = "BH")
  )]
  res[, medication_robust_yang :=
    is.finite(yang_proximity_q) & yang_proximity_q < .05 &
    is.finite(yang_proximity_q_med) & yang_proximity_q_med < .05 &
    sign(yang_proximity_beta10) == sign(yang_proximity_beta10_med) &
    is.finite(prev_vs_yin_q) & prev_vs_yin_q < .05 &
    is.finite(prev_vs_yin_q_med) & prev_vs_yin_q_med < .05 &
    sign(prev_vs_yin_beta) == sign(prev_vs_yin_beta_med)
  ]
  res[, marker_class := fcase(
    is.finite(yin_cox_q) & yin_cox_q < .05 & abs(yin_cox_beta) >= marker_min_cox_loghr &
      is.finite(prev_vs_yin_q) & prev_vs_yin_q < .05 & abs(prev_vs_yin_beta) >= marker_min_prevalence_beta &
      sign(yin_cox_beta) == sign(prev_vs_yin_beta), "prospective + disease-state",
    is.finite(prev_vs_yin_q) & prev_vs_yin_q < .05 & abs(prev_vs_yin_beta) >= marker_min_prevalence_beta &
      (!is.finite(yin_cox_q) | yin_cox_q >= .05 | abs(yin_cox_beta) < marker_min_cox_loghr), "prevalent marker",
    is.finite(yang_proximity_q) & yang_proximity_q < .05 & abs(yang_proximity_beta10) >= marker_min_proximity_beta &
      (!is.finite(yin_cox_q) | yin_cox_q >= .05 | abs(yin_cox_beta) < marker_min_cox_loghr), "reactive/proximity marker",
    is.finite(yin_cox_q) & yin_cox_q < .05 & abs(yin_cox_beta) >= marker_min_cox_loghr &
      (!is.finite(prev_vs_yin_q) | prev_vs_yin_q >= .05 | abs(prev_vs_yin_beta) < marker_min_prevalence_beta), "prospective-only",
    is.finite(stage_q) & stage_q < .05 & abs(stage_beta) >= marker_min_stage_beta, "stage marker",
    default = "weak/mixed"
  )]
  res[, evidence_p := pmin(stage_p, yin_cox_p, yang_proximity_p, prev_vs_yin_p, na.rm = TRUE)]
  res[!is.finite(evidence_p), evidence_p := NA_real_]
  setorder(res, evidence_p, na.last = TRUE)

  labels <- unique(c(requested_biom, res[seq_len(min(20, .N)), biomarker]))
  res[, plot_label := ifelse(biomarker %in% labels, biomarker, "")]
  class_cols <- c(
    `prospective + disease-state` = "#D95F02", `prevalent marker` = "#7570B3",
    `reactive/proximity marker` = "#E7298A", `prospective-only` = "#1B9E77",
    `stage marker` = "#E6AB02", `weak/mixed` = "grey70"
  )
  p1 <- ggplot(res, aes(stage_beta, yin_cox_beta, color = marker_class, label = plot_label)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey55") +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey55") +
    geom_point(alpha = .8, size = 1.8) +
    ggrepel::geom_text_repel(size = 3.0, max.overlaps = 35, seed = seed) +
    scale_color_manual(values = class_cols) +
    labs(title = "a. Risk-age marker vs prospective risk",
         x = "Biomarker difference per 1-SD CKM risk-age position",
         y = "Incident clinical stage-4 Cox beta", color = NULL) +
    theme_ckm(11) + guides(color = guide_legend(nrow = 2, byrow = TRUE)) + theme(legend.position = "bottom")
  p2 <- ggplot(res, aes(prev_vs_yin_beta, yang_proximity_beta10, color = marker_class, label = plot_label)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey55") +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey55") +
    geom_point(alpha = .8, size = 1.8) +
    ggrepel::geom_text_repel(size = 3.0, max.overlaps = 35, seed = seed + 1) +
    scale_color_manual(values = class_cols) +
    labs(title = "b. Prevalence contrast vs Yang event proximity",
         x = "Yang vs disease-free Yin difference", y = "Yang proximity beta per 10 years", color = NULL) +
    theme_ckm(11) + guides(color = guide_legend(nrow = 2, byrow = TRUE)) + theme(legend.position = "bottom")
  counts <- res[, .N, by = marker_class][order(-N)]
  p3 <- ggplot(counts, aes(N, fct_reorder(marker_class, N), fill = marker_class)) +
    geom_col(color = "white") + geom_text(aes(label = N), hjust = -.1, size = 3.5) +
    scale_fill_manual(values = class_cols, guide = "none") +
    scale_x_continuous(expand = expansion(mult = c(0, .18))) +
    labs(title = "c. Evidence classes", x = "Biomarkers", y = NULL) + theme_ckm(11)

  top <- res[seq_len(min(20, .N))]
  top[, biomarker_plot := factor(biomarker, levels = rev(biomarker))]
  top_long <- melt(top[, .(
    biomarker, biomarker_plot,
    Stage = sign(stage_beta) * pmin(-log10(pmax(stage_p, 1e-300)), 20),
    Prospective = sign(yin_cox_beta) * pmin(-log10(pmax(yin_cox_p, 1e-300)), 20),
    Prevalence = sign(prev_vs_yin_beta) * pmin(-log10(pmax(prev_vs_yin_p, 1e-300)), 20),
    Proximity = sign(yang_proximity_beta10) * pmin(-log10(pmax(yang_proximity_p, 1e-300)), 20)
  )], id.vars = c("biomarker", "biomarker_plot"), variable.name = "evidence_axis", value.name = "signed_logp")
  top_long[, biomarker_plot := factor(biomarker, levels = top$biomarker)]
  axis_cols <- c(Stage = "#2C7FB8", Prospective = "#D95F02", Prevalence = "#7570B3", Proximity = "#E7298A")
  p4 <- ggplot(top_long, aes(biomarker_plot, signed_logp, color = evidence_axis, group = evidence_axis)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey55") +
    geom_line(linewidth = .85) + geom_point(size = 2.0) +
    scale_color_manual(values = axis_cols) +
    labs(title = "d. Four-axis evidence profiles",
         subtitle = "Proteins are ordered by the combined evidence ranking; direction is encoded by the sign.",
         x = "Biomarker", y = "Signed -log10(P)", color = NULL) +
    theme_ckm(10) +
    guides(color = guide_legend(nrow = 1, byrow = TRUE)) +
    theme(axis.text.x = element_text(angle = 55, hjust = 1, size = 8.5), legend.position = "top")

  save_plot((p1 | p2) / (p3 | p4), biom_fig_filename("Fig3.marker_map.png"), 17, 12)
  write_tab(res, "marker_map.tsv.gz")
  marker_thresholds <- data.table(
    item = c("minimum absolute stage beta", "minimum absolute prospective log-HR", "minimum absolute Yang-vs-Yin beta", "minimum absolute Yang-proximity beta per 10 years"),
    value = c(marker_min_stage_beta, marker_min_cox_loghr, marker_min_prevalence_beta, marker_min_proximity_beta)
  )
  write_book(list(marker_thresholds = marker_thresholds, marker_map = res,
                  class_counts = counts, evidence_profiles = top_long),
             biom_fig_filename("Fig3.marker_map.out.xlsx"))
  saveRDS(res, file.path(rawdir, "marker_map.rds"), compress = "gzip")
  progress_log("marker_map: completed in ", elapsed_since(t_step))
  res
}'''


YY_TIMELINE = r'''make_yy_timeline <- function(dat, features, meta, clock_obj, marker_scan = NULL) {
  display <- resolve_display_biomarkers(features, meta, marker_scan, n = 6)
  position_map <- get_clock_positions(clock_obj)
  ref <- dat$ckm.stage == 0 & !dat$yang
  cases <- dat[(dat$yang | (dat$yin & dat$event == 1)) & is.finite(dat$b2e) & dat$b2e >= -16 & dat$b2e <= 16, , drop = FALSE]
  if (!nrow(cases)) {
    warning("No target-event cases are available in the -16 to +16 year Yin-Yang window", call. = FALSE)
    return(list(display = display, bins = data.table(), bins_baseline = data.table(), stage_profile = data.table()))
  }
  cases$event_centered_year <- -cases$b2e
  cases$baseline_centered_year <- cases$b2e
  breaks <- seq(-16, 16, 2); mids <- head(breaks, -1) + 1

  build_bin_summary <- function(xvals, phase_fun) {
    rbindlist(lapply(display, function(p) {
      zall <- standardize_ref(dat[[p]], ref)
      z <- zall[match(cases$eid, dat$eid)]
      b <- cut(xvals, breaks = breaks, include.lowest = TRUE)
      x <- data.table(bin = b, z = z)[is.finite(z), .(mean_z = mean(z), se = sd(z) / sqrt(.N), n = .N), by = bin]
      full <- data.table(bin = factor(levels(b), levels = levels(b)), mid = mids)
      x <- merge(full, x, by = "bin", all.x = TRUE, sort = FALSE)
      x[n < 20, `:=`(mean_z = NA_real_, se = NA_real_)]
      x[, phase := phase_fun(mid)]
      x[, biomarker := p]
      x
    }), fill = TRUE)
  }

  bins <- build_bin_summary(cases$event_centered_year, function(mid) ifelse(mid < 0, "Yin (pre-event)", "Yang (post-event)"))
  bins_base <- build_bin_summary(cases$baseline_centered_year, function(mid) ifelse(mid < 0, "Yang before baseline", "Yin after baseline"))
  hist0 <- data.table(mid = mids, count = hist(cases$event_centered_year, breaks = breaks, plot = FALSE)$counts)
  histb <- data.table(mid = mids, count = hist(cases$baseline_centered_year, breaks = breaks, plot = FALSE)$counts)
  hist_max <- max(c(hist0$count, histb$count), na.rm = TRUE)
  if (!is.finite(hist_max) || hist_max <= 0) hist_max <- 1

  stage_profile <- rbindlist(lapply(display, function(p) {
    z <- standardize_ref(dat[[p]], ref)
    data.table(stage = dat$ckm.stage, z = z)[stage %in% 0:4 & is.finite(z), .(
      mean_z = mean(z), se = sd(z) / sqrt(.N), n = .N
    ), by = stage][, `:=`(position = unname(position_map[as.character(stage)]), biomarker = p)]
  }), fill = TRUE)

  legend_one_line <- guides(color = guide_legend(nrow = 1, byrow = TRUE), fill = guide_legend(nrow = 1, byrow = TRUE))
  panel_legend_theme <- theme(
    legend.position = "top", legend.justification = "left",
    legend.box.just = "left", legend.text = element_text(size = 8.3),
    legend.key.width = grid::unit(.45, "cm"), legend.spacing.x = grid::unit(.10, "cm")
  )
  draw_timeline_panel <- function(hist_dat, line_dat, xlab, title, subtitle, left_fill, right_fill) {
    ggplot() +
      geom_col(data = hist_dat[mid < 0], aes(mid, count), width = 1.9, fill = left_fill, color = "white") +
      geom_col(data = hist_dat[mid > 0], aes(mid, count), width = 1.9, fill = right_fill, color = "white") +
      geom_vline(xintercept = 0, linetype = "dotted") +
      geom_line(data = line_dat, aes(mid, mean_z * hist_max / 2 + hist_max / 2, color = biomarker, group = biomarker), linewidth = .9, na.rm = TRUE) +
      geom_point(data = line_dat, aes(mid, mean_z * hist_max / 2 + hist_max / 2, color = biomarker), size = 1.9, na.rm = TRUE) +
      scale_y_continuous("Target-event cases", sec.axis = sec_axis(~ (. - hist_max / 2) * 2 / hist_max, name = "Biomarker z-score")) +
      labs(title = title, subtitle = subtitle, x = xlab, color = NULL) +
      theme_ckm(11) + legend_one_line + panel_legend_theme
  }

  p1 <- draw_timeline_panel(
    hist0, bins, "Years relative to target event",
    "a. Yin-Yang event-centered profile for CKM stage 4",
    "Bars: pre-event Yin (left) and post-event Yang (right). Connected lines compare cohorts, not repeated measures.",
    "#D9EAF7", "#F6DDD8"
  )
  p2 <- draw_timeline_panel(
    histb, bins_base, "Years relative to baseline",
    "b. Same proteins on the baseline-centered axis",
    "Negative: Yang years before baseline; positive: Yin years after baseline. Same-day records are excluded from b2e.",
    "#F6DDD8", "#D9EAF7"
  )
  p3 <- ggplot(stage_profile, aes(position, mean_z, color = biomarker, group = biomarker)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey55") +
    geom_line(linewidth = .9) + geom_point(size = 2.1) +
    geom_errorbar(aes(ymin = mean_z - 1.96 * se, ymax = mean_z + 1.96 * se), width = .08, alpha = .5) +
    scale_x_continuous(breaks = clock_obj$primary[stage %in% 0:3, position_from_stage0], labels = unname(stage_short_labels[as.character(0:3)])) +
    labs(title = "c. Proteins on the prospective CKM risk-age road",
         subtitle = "Risk-age is equivalent hazard age, not elapsed time to the next stage.",
         x = "CKM risk-age advancement vs s0 (years)", y = "Mean biomarker z-score", color = NULL) +
    theme_ckm(11) + legend_one_line + panel_legend_theme
  p4 <- ggplot(stage_profile, aes(factor(stage, levels = 0:4, labels = unname(stage_short_labels[as.character(0:4)])), mean_z, fill = biomarker)) +
    geom_col(position = "dodge") +
    labs(title = "d. Same proteins on the categorical CKM axis",
         subtitle = "s4 is the prevalent Yang disease state; it is not an observed s3→s4 transition.",
         x = "CKM stage", y = "Mean biomarker z-score", fill = NULL) +
    theme_ckm(11) + legend_one_line + panel_legend_theme

  save_plot((p1 | p2) / (p3 | p4), biom_fig_filename("Fig2.yy_timeline.png"), 17, 13.5)
  write_tab(bins, "yy_timeline.bins.tsv.gz")
  write_tab(bins_base, "yy_timeline.baseline_bins.tsv.gz")
  write_tab(stage_profile, "yy_timeline.stage_profile.tsv")
  write_book(list(timeline_bins = bins, baseline_bins = bins_base,
                  quantitative_stage_profile = stage_profile),
             biom_fig_filename("Fig2.yy_timeline.out.xlsx"))
  list(display = display, bins = bins, bins_baseline = bins_base, stage_profile = stage_profile)
}'''


STAGE_SPECIFIC = r'''make_stage_specific <- function(dat, features, meta, marker_scan) {
  t_step <- Sys.time()
  covars <- get_base_covars(dat)
  marker_scan <- as.data.table(copy(marker_scan)); marker_scan[, biomarker := as.character(biomarker)]
  meta <- as.data.table(copy(meta)); meta[, biomarker := as.character(biomarker)]
  scan_features <- if (scan_all_biom) as.character(features) else unique(as.character(c(
    resolve_display_biomarkers(features, meta, marker_scan, 12),
    marker_scan[seq_len(min(50, .N)), biomarker]
  )))
  scan_features <- intersect(scan_features, names(dat))
  if (!length(scan_features)) stop("stage_specific has no biomarkers to scan", call. = FALSE)
  progress_log("stage_specific: scan features=", fmt_int(length(scan_features)))

  probe <- normalize_stage_specific_result(stage_specific_one(dat, scan_features[1], covars), scan_features[1])
  if (!nrow(probe$summary_row)) stop("stage_specific probe returned no summary row", call. = FALSE)
  rr_raw <- map_features(scan_features, function(feature_id) stage_specific_one(dat, as.character(feature_id)[1], covars), label = "stage_specific scan")
  if (length(rr_raw) != length(scan_features)) {
    warning("stage_specific worker-result count mismatch: expected ", length(scan_features), ", obtained ", length(rr_raw), call. = FALSE)
    length(rr_raw) <- length(scan_features)
  }
  rr <- Map(function(z, feature_id) normalize_stage_specific_result(z, feature_id), rr_raw, scan_features)
  betas <- rbindlist(lapply(rr, `[[`, "beta_rows"), fill = TRUE, use.names = TRUE)
  sums <- rbindlist(lapply(rr, `[[`, "summary_row"), fill = TRUE, use.names = TRUE)
  means <- rbindlist(lapply(rr, `[[`, "mean_rows"), fill = TRUE, use.names = TRUE)
  if (nrow(sums) != length(scan_features)) stop("stage_specific internal collection failure: expected ", length(scan_features), " summary rows, obtained ", nrow(sums), call. = FALSE)

  sums[, biomarker := as.character(biomarker)]
  if (nrow(betas)) betas[, biomarker := as.character(biomarker)]
  if (nrow(means)) means[, biomarker := as.character(biomarker)]
  sums[, interaction_q := p.adjust(interaction_p, method = "BH")]
  sums[, stage_specific := analysis_ok %in% TRUE & is.finite(interaction_q) & interaction_q < .05]
  sums <- merge(sums, meta, by = "biomarker", all.x = TRUE, sort = FALSE)
  if (nrow(betas)) {
    betas[, stage_p_q := p.adjust(p, method = "BH"), by = baseline_stage]
    betas[, baseline_stage_label := factor(baseline_stage, levels = 0:3, labels = unname(stage_short_labels[as.character(0:3)]))]
  }
  audit <- sums[, .(scanned = .N, valid = sum(analysis_ok %in% TRUE), failed = sum(!(analysis_ok %in% TRUE)),
                     interaction_p_finite = sum(is.finite(interaction_p)), interaction_q_lt_005 = sum(is.finite(interaction_q) & interaction_q < .05))]
  error_counts <- sums[analysis_ok == FALSE & !is.na(error), .N, by = error][order(-N)]

  profile_counts <- sums[analysis_ok == TRUE & !is.na(profile), .N, by = profile][order(-N)]
  p1 <- if (nrow(profile_counts)) {
    ggplot(profile_counts, aes(N, fct_reorder(profile, N), fill = profile)) +
      geom_col(color = "white") + geom_text(aes(label = N), hjust = -.1, size = 3.6) +
      scale_x_continuous(expand = expansion(mult = c(0, .2))) +
      labs(title = "a. Preclinical s0-s3 protein-profile shapes", x = "Proteins", y = NULL, fill = NULL) +
      theme_ckm(12) + theme(legend.position = "none")
  } else ggplot() + annotate("text", x = 0, y = 0, label = "No valid profile summaries") + theme_void()

  ranked <- unique(c(
    sums[stage_specific == TRUE][order(interaction_q), biomarker],
    sums[analysis_ok == TRUE & is.finite(interaction_p)][order(interaction_p), biomarker],
    marker_scan[order(evidence_p), biomarker]
  ))
  ranked <- intersect(ranked, unique(c(betas$biomarker, means$biomarker)))
  top_interaction <- head(ranked, 8)
  beta_show <- betas[biomarker %in% top_interaction]
  if (nrow(beta_show)) beta_show[, biomarker := factor(biomarker, levels = top_interaction)]
  p2 <- if (nrow(beta_show)) {
    ggplot(beta_show, aes(baseline_stage_label, beta, group = 1)) +
      geom_hline(yintercept = 0, linetype = "dashed", color = "grey60") +
      geom_line(linewidth = .85, color = "#2C7FB8") + geom_point(size = 2.0, color = "#2C7FB8") +
      geom_errorbar(aes(ymin = beta - 1.96 * se, ymax = beta + 1.96 * se), width = .08, alpha = .5, color = "#2C7FB8") +
      facet_wrap(~ biomarker, ncol = 4, scales = "free_y") +
      labs(title = "b. Top stage-context-specific associations with future CKM stage 4",
           subtitle = "Two-row display. FDR hits are prioritized; otherwise the lowest interaction P values are shown.",
           x = "Baseline CKM stage", y = "Cox beta per SD") + theme_ckm(11)
  } else ggplot() + annotate("text", x = 0, y = 0, label = "No valid stage-specific Cox estimates") + theme_void()

  means2 <- merge(means, sums[, .(biomarker, profile, interaction_q)], by = "biomarker", all.x = TRUE, sort = FALSE)
  profile_order <- c("monotone increasing", "monotone decreasing", "preclinical inverted-U", "preclinical U-shaped", "non-monotone")
  representatives <- rbindlist(lapply(profile_order, function(pr) {
    z <- unique(means2[profile == pr & is.finite(interaction_q), .(biomarker, interaction_q)])[order(interaction_q)]
    head(z, 2)
  }), fill = TRUE)
  rep_ids <- unique(representatives$biomarker)
  mean_plot <- means2[biomarker %in% rep_ids]
  if (nrow(mean_plot)) {
    mean_plot[, stage_label := factor(stage, levels = 0:4, labels = unname(stage_short_labels[as.character(0:4)]))]
    mean_plot[, facet_label := paste0(biomarker, "\n", profile)]
  }
  p3 <- if (nrow(mean_plot)) {
    ggplot(mean_plot, aes(stage_label, mean_z, group = 1)) +
      geom_hline(yintercept = 0, linetype = "dashed", color = "grey60") +
      geom_line(linewidth = .85, color = "#D95F02") + geom_point(size = 2.0, color = "#D95F02") +
      geom_errorbar(aes(ymin = mean_z - 1.96 * se, ymax = mean_z + 1.96 * se), width = .08, alpha = .45, color = "#D95F02") +
      facet_wrap(~ facet_label, ncol = 5, scales = "free_y") +
      labs(title = "c. Representative preclinical protein trajectories",
           subtitle = "Two-row display. The s3→s4 point is a separate cross-sectional Yang disease-state contrast.",
           x = "CKM stage", y = "Mean protein z-score") + theme_ckm(10)
  } else ggplot() + annotate("text", x = 0, y = 0, label = "No representative stage trajectories") + theme_void()

  save_plot(p1 / p2 / p3 + plot_layout(heights = c(.55, 1.25, 1.25)), biom_fig_filename("Fig4.stage_specific.png"), 17, 17)
  write_tab(audit, "stage_specific.scan_audit.tsv")
  write_tab(error_counts, "stage_specific.error_counts.tsv")
  write_tab(betas, "stage_specific.cox.tsv.gz")
  write_tab(sums, "stage_specific.summary.tsv.gz")
  write_tab(means, "stage_specific.means.tsv.gz")
  write_book(list(scan_audit = audit, error_counts = error_counts, stage_specific_summary = sums,
                  stage_specific_cox = betas, stage_means = means),
             biom_fig_filename("Fig4.stage_specific.out.xlsx"))
  out <- list(pipeline_version = "v18", betas = betas, summary = sums, means = means, audit = audit, error_counts = error_counts)
  saveRDS(out, file.path(rawdir, "stage_specific.rds"), compress = "gzip")
  progress_log("stage_specific: completed in ", elapsed_since(t_step))
  out
}'''


CAUSAL_TRIAGE = r'''make_causal_triage <- function(dat, features, meta, marker_scan, stage_specific) {
  t_step <- Sys.time()
  covars <- get_base_covars(dat)
  marker_scan <- as.data.table(copy(marker_scan)); marker_scan[, biomarker := as.character(biomarker)]
  if (!"medication_robust_yang" %in% names(marker_scan)) marker_scan[, medication_robust_yang := FALSE]
  meta <- as.data.table(copy(meta)); meta[, biomarker := as.character(biomarker)]
  stage_summary <- if (!is.null(stage_specific$summary)) as.data.table(copy(stage_specific$summary)) else data.table()
  if (nrow(stage_summary)) stage_summary[, biomarker := as.character(biomarker)]
  scan_features <- if (scan_all_biom) as.character(features) else {
    stage_hits <- if (nrow(stage_summary) && "stage_specific" %in% names(stage_summary)) stage_summary[stage_specific == TRUE, biomarker] else character()
    unique(as.character(c(marker_scan[seq_len(min(100, .N)), biomarker], stage_hits)))
  }
  scan_features <- intersect(scan_features, names(dat))
  progress_log("causal_triage: scan features=", fmt_int(length(scan_features)), "; lag years=", paste(lag_years, collapse = ","))

  lag_res <- rbindlist(map_features(scan_features, function(feature_id) {
    feature_id <- as.character(feature_id)[1]
    cols <- unique(c("ckm.stage", "yang", "yin", "event", "time", feature_id, covars))
    d <- dat[, intersect(cols, names(dat)), drop = FALSE]
    if (!feature_id %in% names(d)) return(data.table())
    d$.z <- standardize_ref(d[[feature_id]], d$ckm.stage == 0 & !d$yang)
    rbindlist(lapply(lag_years, function(lg) {
      rr0 <- cox_one(d[d$yin, , drop = FALSE], ".z", covars, lag = lg)
      rr0[, `:=`(biomarker = feature_id, lag_year = as.numeric(lg))]
      rr0
    }), fill = TRUE)
  }, label = "causal_triage lag scan"), fill = TRUE, use.names = TRUE)
  if (!nrow(lag_res)) stop("causal_triage produced no lagged Cox results", call. = FALSE)
  lag_res[, `:=`(biomarker = as.character(biomarker), lag_year = safe_num(lag_year))]
  lag_res[, q := p.adjust(p, method = "BH"), by = lag_year]
  wide <- data.table::dcast(lag_res, biomarker ~ lag_year, value.var = c("beta", "p", "q"))
  wide[, biomarker := as.character(biomarker)]
  setnames(wide, gsub("^beta_", "lag_beta_", names(wide)))
  setnames(wide, gsub("^p_", "lag_p_", names(wide)))
  setnames(wide, gsub("^q_", "lag_q_", names(wide)))
  key_audit <- data.table(
    object = c("marker_scan", "lag_wide", "stage_summary", "metadata"),
    class = c(class(marker_scan$biomarker)[1], class(wide$biomarker)[1], if (nrow(stage_summary)) class(stage_summary$biomarker)[1] else "empty", class(meta$biomarker)[1]),
    N = c(nrow(marker_scan), nrow(wide), nrow(stage_summary), nrow(meta)),
    unique_keys = c(uniqueN(marker_scan$biomarker), uniqueN(wide$biomarker), if (nrow(stage_summary)) uniqueN(stage_summary$biomarker) else 0L, uniqueN(meta$biomarker))
  )
  write_tab(key_audit, "causal_triage.key_audit.tsv")
  triage <- merge(marker_scan, wide, by = "biomarker", all.x = TRUE, sort = FALSE)
  stage_cols <- c("biomarker", "interaction_p", "interaction_q", "profile", "preclinical_profile", "yang_shift", "yang_direction", "delta01", "delta12", "delta23", "delta34", "analysis_ok")
  if (nrow(stage_summary)) {
    stage_add <- stage_summary[, intersect(stage_cols, names(stage_summary)), with = FALSE]
    triage <- merge(triage, stage_add, by = "biomarker", all.x = TRUE, sort = FALSE)
  } else {
    triage[, `:=`(interaction_p = NA_real_, interaction_q = NA_real_, profile = NA_character_, delta01 = NA_real_, delta12 = NA_real_, delta23 = NA_real_, delta34 = NA_real_, analysis_ok = NA)]
    warning("stage_specific summary is empty; causal triage will proceed without interaction evidence", call. = FALSE)
  }
  meta_cols <- setdiff(names(meta), names(triage))
  if (length(meta_cols)) triage <- merge(triage, meta[, c("biomarker", meta_cols), with = FALSE], by = "biomarker", all.x = TRUE, sort = FALSE)
  lag5b <- if ("lag_beta_5" %in% names(triage)) triage$lag_beta_5 else rep(NA_real_, nrow(triage))
  lag5p <- if ("lag_p_5" %in% names(triage)) triage$lag_p_5 else rep(NA_real_, nrow(triage))
  lag5q <- if ("lag_q_5" %in% names(triage)) triage$lag_q_5 else rep(NA_real_, nrow(triage))
  triage[, `:=`(lag5_beta = safe_num(lag5b), lag5_p = safe_num(lag5p), lag5_q = safe_num(lag5q))]
  triage[, retention5 := abs(lag5_beta) / pmax(abs(yin_cox_beta), 1e-8)]
  triage[, attenuation5 := 1 - retention5]
  triage[, lag_persistent :=
    is.finite(yin_cox_q) & yin_cox_q < .05 & is.finite(lag5_q) & lag5_q < .05 &
    is.finite(yin_cox_beta) & abs(yin_cox_beta) >= marker_min_cox_loghr &
    is.finite(lag5_beta) & abs(lag5_beta) >= marker_min_cox_loghr &
    sign(lag5_beta) == sign(yin_cox_beta) & is.finite(retention5) & retention5 >= lag_min_retention]
  triage[, evidence_class := fcase(
    lag_persistent %in% TRUE & is.finite(interaction_q) & interaction_q < .05, "stage-specific lag-persistent association",
    lag_persistent %in% TRUE, "lag-persistent prospective association",
    is.finite(yin_cox_beta) & is.finite(prev_vs_yin_beta) & sign(yin_cox_beta) != sign(prev_vs_yin_beta) &
      is.finite(yin_cox_q) & yin_cox_q < .05 & is.finite(prev_vs_yin_q) & prev_vs_yin_q < .05, "mirror/opposing",
    medication_robust_yang %in% TRUE & is.finite(prev_vs_yin_q) & prev_vs_yin_q < .05 &
      is.finite(yang_proximity_q) & yang_proximity_q < .05 & !(lag_persistent %in% TRUE), "medication-robust reactive marker-like",
    !(medication_robust_yang %in% TRUE) & is.finite(prev_vs_yin_q) & prev_vs_yin_q < .05 &
      is.finite(yang_proximity_q) & yang_proximity_q < .05 & !(lag_persistent %in% TRUE), "treatment-sensitive Yang pattern",
    is.finite(interaction_q) & interaction_q < .05, "stage-context-specific",
    default = "mixed/uncertain")]
  triage[, triage_score := pmin(-log10(pmax(yin_cox_p, 1e-300)), 30) +
    .7 * pmin(-log10(pmax(lag5_p, 1e-300)), 20) +
    .5 * pmin(-log10(pmax(interaction_p, 1e-300)), 20) -
    .3 * pmin(-log10(pmax(yang_proximity_p, 1e-300)), 20) *
      (evidence_class %in% c("medication-robust reactive marker-like", "treatment-sensitive Yang pattern"))]
  triage[!is.finite(triage_score), triage_score := NA_real_]
  setorder(triage, -triage_score, na.last = TRUE)
  mr_candidates <- triage[evidence_class %in% c("stage-specific lag-persistent association", "lag-persistent prospective association", "stage-context-specific")]
  mr_candidates[, note := "Observational triage only; external cis-pQTL/mQTL MR, colocalization, and biological validation are required"]
  class_cols <- c(
    `stage-specific lag-persistent association` = "#D95F02", `lag-persistent prospective association` = "#1B9E77",
    `medication-robust reactive marker-like` = "#E7298A", `treatment-sensitive Yang pattern` = "#CC79A7",
    `mirror/opposing` = "#7570B3", `stage-context-specific` = "#E6AB02", `mixed/uncertain` = "grey70")

  show0 <- triage[is.finite(yin_cox_beta) & is.finite(lag5_beta)]
  if (nrow(show0) >= 10) {
    qx <- quantile(show0$yin_cox_beta, c(.02, .98), na.rm = TRUE, names = FALSE)
    qy <- quantile(show0$lag5_beta, c(.02, .98), na.rm = TRUE, names = FALSE)
    show <- show0[yin_cox_beta >= qx[1] & yin_cox_beta <= qx[2] & lag5_beta >= qy[1] & lag5_beta <= qy[2]]
  } else show <- copy(show0)
  n_trim <- nrow(show0) - nrow(show)
  show <- show[order(-abs(yin_cox_beta - lag5_beta))][seq_len(min(100, .N))]
  show[, plot_label := ifelse(rank(-abs(yin_cox_beta) - abs(lag5_beta), ties.method = "first") <= 24, biomarker, "")]
  p1 <- ggplot(show, aes(yin_cox_beta, lag5_beta, color = evidence_class, label = plot_label)) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey55") +
    geom_hline(yintercept = 0, linetype = "dotted") + geom_vline(xintercept = 0, linetype = "dotted") +
    geom_point(size = 1.8) + ggrepel::geom_text_repel(size = 3.0, max.overlaps = 32, seed = seed) +
    scale_color_manual(values = class_cols) +
    labs(title = "a. Persistence after a 5-year lag",
         subtitle = if (n_trim > 0) paste0("Robust plotting window; ", n_trim, " extreme outlier(s) omitted from the panel only.") else NULL,
         x = "Main prospective Cox beta", y = "Lag-5-year Cox beta", color = NULL) +
    theme_ckm(11) + guides(color = guide_legend(nrow = 2, byrow = TRUE)) + theme(legend.position = "bottom")

  context <- triage[is.finite(retention5) & is.finite(interaction_q)]
  context[, retention_plot := pmin(retention5, 2.5)]
  context[, interaction_evidence := -log10(pmax(interaction_q, 1e-300))]
  context[, plot_label := ifelse(rank(-triage_score, ties.method = "first") <= 22, biomarker, "")]
  p2 <- ggplot(context, aes(retention_plot, interaction_evidence, color = evidence_class, label = plot_label)) +
    geom_vline(xintercept = lag_min_retention, linetype = "dashed", color = "grey55") +
    geom_hline(yintercept = -log10(.05), linetype = "dashed", color = "grey55") +
    geom_point(size = 1.8, alpha = .85) + ggrepel::geom_text_repel(size = 2.8, max.overlaps = 28, seed = seed + 2) +
    scale_color_manual(values = class_cols) +
    labs(title = "b. Persistence and stage-context evidence",
         subtitle = "Upper-right signals combine lag retention with biomarker-by-stage interaction support.",
         x = "Lag-5/main absolute effect ratio (capped at 2.5)", y = "-log10(stage-interaction FDR)", color = NULL) +
    theme_ckm(11) + guides(color = guide_legend(nrow = 2, byrow = TRUE)) + theme(legend.position = "bottom")

  lag_rank <- lag_res[is.finite(beta), .(lag_range = diff(range(beta, na.rm = TRUE)), absmax = max(abs(beta), na.rm = TRUE)), by = biomarker][order(-lag_range, -absmax)]
  top_lag <- lag_rank[seq_len(min(16, .N)), biomarker]
  lag_show <- lag_res[biomarker %in% top_lag]
  p3 <- ggplot(lag_show, aes(lag_year, beta, color = biomarker, group = biomarker)) +
    geom_hline(yintercept = 0, linetype = "dashed") + geom_line(linewidth = .9) + geom_point(size = 1.9) +
    labs(title = "c. Lag-response trajectories",
         subtitle = "Signals with the largest change after excluding early events; shrinking curves are more compatible with reverse causation.",
         x = "Excluded early follow-up (years)", y = "Cox beta", color = NULL) +
    theme_ckm(10) + guides(color = guide_legend(nrow = 2, byrow = TRUE)) +
    theme(legend.position = "bottom", legend.text = element_text(size = 8.5))

  yy_contrast <- triage[is.finite(yin_cox_beta) & is.finite(yang_proximity_beta10)]
  yy_contrast[, plot_label := ifelse(rank(-triage_score, ties.method = "first") <= 24, biomarker, "")]
  p4 <- ggplot(yy_contrast, aes(yin_cox_beta, yang_proximity_beta10, color = evidence_class, label = plot_label)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey55") +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey55") +
    geom_point(size = 1.8, alpha = .85) + ggrepel::geom_text_repel(size = 2.8, max.overlaps = 30, seed = seed + 3) +
    scale_color_manual(values = class_cols) +
    labs(title = "d. Prospective effect vs Yang event proximity",
         subtitle = "Large proximity-only effects suggest reactive disease-state biology rather than durable prospective information.",
         x = "Prospective Cox beta", y = "Yang proximity beta per 10 years", color = NULL) +
    theme_ckm(11) + guides(color = guide_legend(nrow = 2, byrow = TRUE)) + theme(legend.position = "bottom")

  save_plot((p1 | p2) / (p3 | p4) + plot_layout(heights = c(.92, 1.25)), biom_fig_filename("Fig5.causal_triage.png"), 17, 14.5)
  triage_thresholds <- data.table(
    criterion = c("main FDR", "lag-5 FDR", "minimum absolute log-HR", "minimum lag/main effect retention", "direction"),
    threshold = c("<0.05", "<0.05", sprintf(">= %.4f [HR >= %.2f or <= %.2f]", marker_min_cox_loghr, exp(marker_min_cox_loghr), exp(-marker_min_cox_loghr)), sprintf(">= %.2f", lag_min_retention), "same")
  )
  write_tab(triage_thresholds, "causal_triage.thresholds.tsv")
  write_tab(lag_res, "causal_triage.lagged_cox.tsv.gz")
  write_tab(triage, "causal_triage.summary.tsv.gz")
  write_tab(mr_candidates, "causal_triage.mr_coloc_candidates.tsv")
  write_book(list(key_audit = key_audit, thresholds = triage_thresholds, triage = triage,
                  lagged_cox = lag_res, mr_coloc_candidates = mr_candidates),
             biom_fig_filename("Fig5.causal_triage.out.xlsx"))
  saveRDS(list(triage = triage, lagged_cox = lag_res, mr_candidates = mr_candidates), file.path(rawdir, "causal_triage.rds"), compress = "gzip")
  progress_log("causal_triage: completed in ", elapsed_since(t_step))
  list(triage = triage, lagged_cox = lag_res, mr_candidates = mr_candidates)
}'''


YY_SCHEMA_HELPER = r'''
ensure_yy_weight_schema <- function(x) {
  if (is.null(x)) return(data.table())
  x <- as.data.table(copy(x))
  if (!nrow(x)) return(x)
  logical_cols <- c("lag_support", "yang_support", "stage_support", "reactive_penalty", "selected_yy", "selected_yin")
  numeric_cols <- c("weight", "cox_beta", "normalized_weight", "yin_normalized_weight", "filter_score")
  for (v in logical_cols) {
    if (!v %in% names(x)) x[, (v) := FALSE]
    x[is.na(get(v)), (v) := FALSE]
  }
  for (v in numeric_cols) if (!v %in% names(x)) x[, (v) := NA_real_]
  if (!"reason" %in% names(x)) x[, reason := "Yin Cox"]
  if (!"biomarker" %in% names(x)) stop("YY weight table lacks biomarker column", call. = FALSE)
  x[, biomarker := as.character(biomarker)]
  x
}
'''


COMPONENT_BLOCK_NEW = r'''    component_defs <- list(
      M3a_persistent_YY = yy_weights[lag_support %in% TRUE],
      M3b_stage_context_YY = yy_weights[stage_support %in% TRUE],
      M3c_yang_supported_YY = yy_weights[yang_support %in% TRUE]
    )
    component_models <- lapply(names(component_defs), function(m) {
      subw <- component_defs[[m]]
      if (!nrow(subw)) {
        explain_rows[[length(explain_rows) + 1L]] <<- data.table(
          fold = fold, model = m, biomarkers = 0L, abs_weight = 0,
          reasons = "No biomarkers met this evidence channel in the training fold"
        )
        # Optional explanation channels must not terminate the primary model.
        return(lp2)
      }
      subtr <- score_from_weights_or_na(bm$train, subw)
      subte <- score_from_weights_or_na(bm$test, subw)
      if (all(!is.finite(subtr)) || all(!is.finite(subte)) ||
          sd(subtr, na.rm = TRUE) <= 1e-8 || sd(subte, na.rm = TRUE) <= 1e-8) {
        lp <- lp2
        reason0 <- "Evidence channel produced a constant/invalid score; risk-age model used as neutral fallback"
      } else {
        lp0 <- safe_fit_score_model(
          ytr, etr,
          cbind(covd$train, stage_quant$train),
          cbind(covd$test, stage_quant$test),
          subtr, subte, paste0(m, "_score")
        )
        lp <- if (all(!is.finite(lp0)) || sd(lp0, na.rm = TRUE) <= 1e-10) lp2 else lp0
        reason0 <- paste(sort(unique(subw$reason)), collapse = "; ")
      }
      explain_rows[[length(explain_rows) + 1L]] <<- data.table(
        fold = fold, model = m, biomarkers = nrow(subw),
        abs_weight = sum(abs(subw$normalized_weight), na.rm = TRUE),
        reasons = reason0
      )
      lp
    })
    names(component_models) <- names(component_defs)'''


def patch_r(text: str) -> str:
    required = [
        "calculate_prevent <- function", "make_all_source_ckm4_dates <- function",
        "prepare_target_data <- function", "make_marker_map <- function",
        "make_yy_timeline <- function", "make_stage_specific <- function",
        "make_causal_triage <- function", "learn_yy_weights", "component_defs <- list("
    ]
    missing = [x for x in required if x not in text]
    if missing:
        raise PatchError("Input does not look like the expected ckm.v16.R; missing: " + ", ".join(missing))

    text = replace_function(text, "calculate_prevent", CALCULATE_PREVENT)
    text = replace_function(text, "make_all_source_ckm4_dates", ALL_SOURCE_DATES)
    text = replace_function(text, "prepare_target_data", PREPARE_TARGET)
    text = replace_function(text, "make_marker_map", MARKER_MAP)
    text = replace_function(text, "make_yy_timeline", YY_TIMELINE)
    text = replace_function(text, "make_stage_specific", STAGE_SPECIFIC)
    text = replace_function(text, "make_causal_triage", CAUSAL_TRIAGE)

    # Strict AHA primary TG definition: fasting TG >=150 mg/dL. The existing
    # opt-in flag now activates the pragmatic unknown/nonfasting sensitivity.
    old_tg = '''  tg_threshold <- ifelse(fasting_ok, tg_fasting_cutoff, tg_nonfasting_cutoff)\n  tg_usable <- is.finite(safe_num(dat$bb_TG))\n  tg_fasting_mode <- if (is.na(fasting_col)) {\n    "fasting time unavailable: nonfasting threshold applied"\n  } else {\n    "fasting>=8h uses fasting threshold; <8h/unknown uses nonfasting threshold"\n  }'''
    new_tg = '''  tg_value <- safe_num(dat$bb_TG)\n  if (allow_unknown_fasting_tg) {\n    tg_threshold <- ifelse(fasting_ok, tg_fasting_cutoff, tg_nonfasting_cutoff)\n    tg_usable <- is.finite(tg_value)\n    tg_fasting_mode <- if (is.na(fasting_col)) {\n      "pragmatic sensitivity: fasting time unavailable; unknown status uses nonfasting threshold"\n    } else {\n      "pragmatic sensitivity: fasting>=8h uses fasting threshold; <8h/unknown uses nonfasting threshold"\n    }\n  } else {\n    tg_threshold <- rep(tg_fasting_cutoff, nrow(dat))\n    tg_usable <- is.finite(tg_value) & fasting_ok\n    tg_fasting_mode <- if (is.na(fasting_col)) {\n      "strict AHA primary: fasting time unavailable; TG criterion not assigned"\n    } else {\n      "strict AHA primary: only fasting>=8h TG values are used"\n    }\n  }'''
    text = replace_once(text, old_tg, new_tg, "strict/pragmatic triglyceride block")

    # Correct the Stage-3 coverage audit: eGFR availability alone is not a
    # complete KDIGO heat-map assessment when UACR is missing.
    old_stage3_asc = '''  stage3_ascertainment <- data.table(
    component = c("PREVENT base", "PREVENT best available", "very-high-risk CKD", "subclinical ASCVD", "pre-HF", "UACR measured"),
    assessed_N = c(
      sum(is.finite(dat$prevent_cvd10_base)), sum(is.finite(dat$prevent_cvd10_best)),
      sum(is.finite(dat$egfr)), sum(sub_ascvd_observed), sum(pre_hf_observed), sum(is.finite(dat$uacr_mg_g))
    ),
    positive_N = c(
      sum(dat$stage3_prevent_component_base, na.rm = TRUE), sum(dat$stage3_prevent_component_best, na.rm = TRUE),
      sum(dat$stage3_ckd_component, na.rm = TRUE), sum(dat$stage3_subclinical_ascvd_component, na.rm = TRUE),
      sum(dat$stage3_prehf_component, na.rm = TRUE), sum(is.finite(dat$uacr_mg_g) & dat$uacr_mg_g >= 30)
    ),
    note = c(
      "uniform primary PREVENT equation", "optional HbA1c/UACR-enhanced sensitivity",
      "KDIGO very-high-risk category", "requires supplied imaging/ABI/subclinical flag",
      "requires supplied pre-HF flag", "missing UACR can under-classify CKD risk"
    )
  )'''
    new_stage3_asc = '''  stage3_ascertainment <- data.table(
    component = c(
      "PREVENT base", "PREVENT best available",
      "complete KDIGO eGFR+UACR", "eGFR-only G4/G5",
      "subclinical ASCVD", "pre-HF", "UACR measured"
    ),
    assessed_N = c(
      sum(is.finite(dat$prevent_cvd10_base)), sum(is.finite(dat$prevent_cvd10_best)),
      sum(is.finite(dat$egfr) & dat$uacr_assessed), sum(is.finite(dat$egfr)),
      sum(sub_ascvd_observed), sum(pre_hf_observed), sum(dat$uacr_assessed)
    ),
    positive_N = c(
      sum(dat$stage3_prevent_component_base, na.rm = TRUE), sum(dat$stage3_prevent_component_best, na.rm = TRUE),
      sum(dat$stage3_ckd_component & is.finite(dat$egfr) & dat$uacr_assessed, na.rm = TRUE),
      sum(is.finite(dat$egfr) & dat$egfr < 30, na.rm = TRUE),
      sum(dat$stage3_subclinical_ascvd_component, na.rm = TRUE),
      sum(dat$stage3_prehf_component, na.rm = TRUE),
      sum(dat$uacr_assessed & is.finite(dat$uacr_mg_g) & dat$uacr_mg_g >= 30)
    ),
    note = c(
      "uniform primary PREVENT base equation", "auto-selected optional-model sensitivity",
      "complete KDIGO heat-map ascertainment requires both eGFR and UACR",
      "eGFR alone can identify G4/G5 very-high-risk CKD but misses albuminuria-defined risk",
      "requires supplied imaging/ABI/subclinical flag",
      "requires supplied pre-HF flag",
      "missing UACR under-classifies stages 2 and 3"
    )
  )'''
    text = replace_once(text, old_stage3_asc, new_stage3_asc, "Stage-3 ascertainment audit")

    # Integrate baseline UACR interval handling and tri-state GP Stage-3
    # components produced by phe.v26.R.
    uacr_anchor = '  progress_log("build_ckm: UACR source=", uacr_source)\n'
    uacr_insert = '''  progress_log("build_ckm: UACR source=", uacr_source)
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
'''
    text = replace_once(text, uacr_anchor, uacr_insert, "UACR interval/assessment integration")

    stage3_component_anchor = '''  sub_ascvd[is.na(sub_ascvd)] <- FALSE
  pre_hf[is.na(pre_hf)] <- FALSE
'''
    stage3_component_insert = '''  sub_ascvd[is.na(sub_ascvd)] <- FALSE
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
'''
    text = replace_once(
        text, stage3_component_anchor, stage3_component_insert,
        "tri-state subclinical ASCVD/pre-HF integration"
    )

    text = replace_once(
        text,
        '''  stage3_negative_assessed <- is.finite(dat$uacr_mg_g) & is.finite(dat$prevent_cvd10) &
    sub_ascvd_observed & pre_hf_observed''',
        '''  stage3_negative_assessed <- dat$uacr_assessed & is.finite(dat$prevent_cvd10) &
    sub_ascvd_observed & pre_hf_observed''',
        "Stage-3 negative-assessment rule"
    )

    compact_pattern = r'''  compact <- dat %>% transmute\(
    eid = as.character\(eid\),
    ckm_definition_version = "2026_AHA_v16",
.*?
  \)
  dir.create\(dirname\(ckm_file\), recursive = TRUE, showWarnings = FALSE\)'''
    compact_new = '''  compact <- dat %>% transmute(
    eid = as.character(eid),
    ckm_definition_version = "2026_AHA_v16",
    ckm.stage = as.integer(ckm.stage),
    date_first_clinical_cvd = as_date2(date_first_clinical_cvd),
    date_ckm4 = as_date2(date_ckm4),
    date_ckm4_icd10 = as_date2(date_ckm4_all),
    date_ckm4_allsource = as_date2(date_ckm4_allsource),
    date_ckm_death = as_date2(date_ckm_death),
    ckm_stage_certainty = as.character(ckm_stage_certainty),
    ckm_stage2_subtype = as.character(ckm_stage2_subtype),
    ckm_stage3_subtype = as.character(ckm_stage3_subtype),
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
    stage3_ckd_component = as.logical(stage3_ckd_component),
    stage3_subclinical_ascvd_component = as.logical(stage3_subclinical_ascvd_component),
    stage3_prehf_component = as.logical(stage3_prehf_component)
  )
  dir.create(dirname(ckm_file), recursive = TRUE, showWarnings = FALSE)'''
    text = regex_replace_once(
        text, compact_pattern, compact_new,
        "expanded compact ckm.rds", flags=re.S
    )

    imaging_anchor = '''  progress_log("Saved compact ckm.rds: ", ckm_file)
  print(table(compact$ckm.stage, useNA = "always"))
'''
    imaging_insert = '''  progress_log("Saved compact ckm.rds: ", ckm_file)
  print(table(compact$ckm.stage, useNA = "always"))

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
'''
    text = replace_once(text, imaging_anchor, imaging_insert, "imaging-landmark CKM output")

    text = replace_once(
        text,
        '''    egfr, uacr_mg_g, ckd_risk, ckd_moderate_high, ckd_very_high, kidney_failure,''',
        '''    egfr, uacr_mg_g, uacr_category, uacr_assessed,
    uacr_lower_mg_g, uacr_upper_mg_g, uacr_interval_uncertain,
    ckd_risk, ckd_moderate_high, ckd_very_high, kidney_failure,''',
        "CKM audit UACR details"
    )
    text = replace_once(
        text,
        '''    subclinical_ascvd = sub_ascvd, subclinical_ascvd_observed = sub_ascvd_observed,
    pre_hf, pre_hf_observed,''',
        '''    subclinical_ascvd, subclinical_ascvd_observed,
    date_subclinical_ascvd, subclinical_ascvd_source,
    pre_hf, pre_hf_observed, date_pre_hf, pre_hf_source,''',
        "CKM audit Stage-3 GP details"
    )

    repeated_optional = '"ckm_stage_certainty", "ckm_stage2_subtype", "ckm_stage3_subtype", "ckm_stage4_substage"'
    expanded_optional = '''"ckm_stage_certainty", "ckm_stage2_subtype", "ckm_stage3_subtype", "ckm_stage4_substage",
    "egfr", "uacr_mg_g", "uacr_category", "uacr_assessed",
    "uacr_lower_mg_g", "uacr_upper_mg_g", "uacr_interval_uncertain", "ckd_risk",
    "subclinical_ascvd", "subclinical_ascvd_observed", "date_subclinical_ascvd", "subclinical_ascvd_source",
    "pre_hf", "pre_hf_observed", "date_pre_hf", "pre_hf_source",
    "stage3_prevent_component", "stage3_ckd_component",
    "stage3_subclinical_ascvd_component", "stage3_prehf_component"'''
    if repeated_optional not in text:
        raise PatchError("Cannot find compact/optional CKM field list for expansion")
    text = text.replace(repeated_optional, expanded_optional)

    # Supplementary, rather than main-figure, stage reconstruction diagnostic.
    text = replace_once(text, 'Fig3b.stage_refinement.png', 'FigS2.stage_refinement.png', "stage-refinement figure demotion")
    text = replace_once(text, 'Fig3b.stage_refinement.out.xlsx', 'FigS2.stage_refinement.out.xlsx', "stage-refinement workbook demotion")

    # Add a mutually interpretable Stage-2 overlap audit.
    anchor = '''  stage2_overlap <- dat %>%\n    transmute('''
    insertion = '''  stage2_subtypes <- dat %>%\n    transmute(\n      final_stage = ckm.stage,\n      component_count = rowSums(cbind(T2D_def, HT_def, hyperTG, metabolic_syndrome, ckd_moderate_high), na.rm = TRUE),\n      dominant_component = case_when(\n        ckd_moderate_high ~ "CKD",\n        T2D_def ~ "type 2 diabetes",\n        metabolic_syndrome ~ "metabolic syndrome",\n        hyperTG ~ "hypertriglyceridemia",\n        HT_def ~ "hypertension",\n        TRUE ~ NA_character_\n      )\n    ) %>%\n    filter(component_count > 0) %>%\n    count(final_stage, component_count, dominant_component, name = "N") %>%\n    arrange(final_stage, component_count, desc(N))\n\n  stage2_overlap <- dat %>%\n    transmute('''
    text = replace_once(text, anchor, insertion, "stage-2 subtype audit insertion")
    text = replace_once(
        text,
        '  write_tab(stage2_overlap, "ckm_stage2_criterion_overlap.tsv", auditdir)\n',
        '  write_tab(stage2_overlap, "ckm_stage2_criterion_overlap.tsv", auditdir)\n  write_tab(stage2_subtypes, "ckm_stage2_subtypes.tsv", auditdir)\n',
        "stage-2 subtype table write"
    )
    text = replace_once(
        text,
        '    stage2_overlap = stage2_overlap,\n',
        '    stage2_overlap = stage2_overlap,\n    stage2_subtypes = stage2_subtypes,\n',
        "stage-2 subtype workbook entry"
    )

    # Prediction schema guard and safe optional-channel fallbacks.
    score_anchor = 'score_from_weights <- function(X, weights) {'
    if YY_SCHEMA_HELPER.strip() not in text:
        text = replace_once(text, score_anchor, YY_SCHEMA_HELPER + '\n' + score_anchor, "YY schema helper insertion")
    text = regex_replace_once(
        text,
        r'weights\s*<-\s*learn_yy_weights\(dat,\s*features,\s*train_rows,\s*yang_rows,\s*base_covars,\s*label\s*=\s*paste0\("prediction fold ",\s*fold\)\)',
        'weights <- ensure_yy_weight_schema(learn_yy_weights(dat, features, train_rows, yang_rows, base_covars, label = paste0("prediction fold ", fold)))',
        "YY weight schema call"
    )
    component_pattern = r'''    component_defs <- list\(\n      M3a_persistent_YY = yy_weights\[lag_support %in% TRUE\],\n      M3b_stage_context_YY = yy_weights\[stage_support %in% TRUE\],\n      M3c_yang_supported_YY = yy_weights\[yang_support %in% TRUE\]\n    \)\n    component_models <- lapply\(names\(component_defs\), function\(m\) \{.*?\n    \}\)\n    names\(component_models\) <- names\(component_defs\)'''
    text = regex_replace_once(text, component_pattern, COMPONENT_BLOCK_NEW, "optional YY component block", flags=re.S)

    # Larger default text throughout the figure suite.
    text = text.replace('theme_ckm <- function(base_size = 10)', 'theme_ckm <- function(base_size = 11)', 1)
    text = text.replace('theme_minimal_ckm <- function(base_size = 10)', 'theme_minimal_ckm <- function(base_size = 11)', 1)

    # Use the non-misleading name biom_refine for the biochemistry-only
    # refinement analysis, its pipeline step, and all output prefixes.
    text = text.replace('stage_refinement', 'biom_refine')

    # Version/stamp changes are performed last.
    text = text.replace('2026_AHA_v16', '2026_AHA_v18')
    text = text.replace('risk_age_v16', 'risk_age_v18')
    text = text.replace('stage_refinement', 'biom_refine')
    text = re.sub(r'\bv16\b', 'v18', text)
    text = text.replace('ckm.v16.sh', 'ckm.v18.sh')
    text = text.replace('ckm.v16.R', 'ckm.v18.R')
    return text


def patch_sh(text: str) -> str:
    if 'ckm.v16.sh' not in text or 'pipeline v16' not in text:
        raise PatchError('Input does not look like ckm.v16.sh')
    text = text.replace('ckm.v16.sh', 'ckm.v18.sh')
    text = text.replace('pipeline v16', 'pipeline v18')
    text = text.replace('v16 no', 'v18 no')
    text = text.replace('v16 classifies unknown\n                             fasting status with the nonfasting TG cutoff.',
                        'v18 uses fasting TG only in the strict primary definition.\n                             This flag enables the pragmatic nonfasting/unknown sensitivity.')
    text = re.sub(r'\bv16\b', 'v18', text)
    return text


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--r', default='ckm.v16.R', help='Input ckm.v16.R')
    p.add_argument('--sh', default='ckm.v16.sh', help='Input ckm.v16.sh')
    p.add_argument('--out-r', default='ckm.v18.R', help='Output R script')
    p.add_argument('--out-sh', default='ckm.v18.sh', help='Output shell script')
    args = p.parse_args()

    r_in = Path(args.r)
    sh_in = Path(args.sh)
    if not r_in.is_file():
        raise SystemExit(f'Input R file not found: {r_in}')
    if not sh_in.is_file():
        raise SystemExit(f'Input shell file not found: {sh_in}')

    try:
        r_out = patch_r(r_in.read_text(encoding='utf-8'))
        sh_out = patch_sh(sh_in.read_text(encoding='utf-8'))
    except PatchError as e:
        raise SystemExit(f'Patch aborted: {e}') from e

    out_r = Path(args.out_r)
    out_sh = Path(args.out_sh)
    out_r.write_text(r_out, encoding='utf-8', newline='\n')
    out_sh.write_text(sh_out, encoding='utf-8', newline='\n')
    out_sh.chmod(out_sh.stat().st_mode | stat.S_IXUSR | stat.S_IXGRP | stat.S_IXOTH)
    print(f'Wrote {out_r}')
    print(f'Wrote {out_sh}')
    print('Run: ./ckm.v18.sh --biom prot,met --scan-all-biom --start-step build_ckm --stop-step consolidate')


if __name__ == '__main__':
    main()
