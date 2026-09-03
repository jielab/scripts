# C1 inherited-score scan utilities.
#
# A COJO PGS is fixed at conception, but it is not the biomarker concentration
# "measured at birth".  These functions therefore keep PGS and adult measured
# omics in parallel columns and never substitute one for the other.

find_c1_pgs_file <- function(layer) {
  explicit <- Sys.getenv("C1_PGS_FILE", unset = "")
  automatic <- file.path(indir, "Rdata", if (layer == "protein") "prot.pgs.rds" else "met.pgs.rds")
  hit <- unique(c(explicit, automatic))
  hit <- hit[nzchar(hit) & file.exists(hit) & file.size(hit) > 0]
  if (length(hit)) normalizePath(hit[[1]], winslash = "/", mustWork = FALSE) else NA_character_
}

c1_pgs_signature <- function(layer) {
  f <- find_c1_pgs_file(layer)
  if (is.na(f)) return("missing")
  i <- file.info(f)
  paste(normalizePath(f, winslash = "/", mustWork = FALSE), i$size,
        format(i$mtime, "%Y-%m-%dT%H:%M:%S%z"), sep = "|")
}

read_c1_pgs <- function(file) {
  x <- if (grepl("\\.rds$", file, ignore.case = TRUE)) readRDS(file) else
    data.table::fread(file, showProgress = FALSE, check.names = FALSE)
  x <- as_tibble(x)
  if (!"eid" %in% names(x)) stop("PGS file must contain eid: ", file, call. = FALSE)
  # Keep this conversion: UKB joins otherwise fail when one input stores eid
  # as integer64 and another stores it as character.
  x$eid <- as.character(x$eid)
  x
}

map_c1_pgs_columns <- function(features, nms) {
  one <- function(f) {
    candidates <- c(paste0(f, ".pgs"), paste0(f, "_pgs"), paste0(f, ".PGS"),
      paste0(f, "_PGS"), paste0(f, "_GRS"), paste0("GRS_", f), f)
    hit <- candidates[candidates %in% nms]
    if (length(hit)) hit[[1]] else NA_character_
  }
  ans <- setNames(vapply(features, one, character(1)), features)
  ans[!is.na(ans)]
}

empty_c1_pgs_assoc <- function(features, score_columns = character()) {
  tibble(term = features, score_column = unname(score_columns[features]),
    estimate = NA_real_, beta = NA_real_, std.error = NA_real_, conf.low = NA_real_,
    conf.high = NA_real_, statistic = NA_real_, p.value = NA_real_,
    N_total = NA_integer_, N_event = NA_integer_, FDR = NA_real_)
}

run_c1_pgs_scan <- function(layer, features, covars, outcome = Y, rawdir = NULL,
  overlap_eids = NULL) {
  score_file <- find_c1_pgs_file(layer)
  overlap_eids <- as.character(overlap_eids %||% character())
  overlap_signature <- if (!length(overlap_eids)) "none" else paste(length(overlap_eids),
    min(overlap_eids), max(overlap_eids), sep = "|")
  signature <- paste(c1_pgs_signature(layer), "omic_overlap", overlap_signature, sep = "|")
  cache <- if (is.null(rawdir)) NA_character_ else file.path(rawdir, "c1.pgs_scan.rds")
  if (!is.na(cache) && cache_valid(cache)) {
    old <- tryCatch(readRDS(cache), error = function(e) NULL)
    if (!is.null(old) && identical(old$signature %||% NA_character_, signature)) {
      message("C1/", layer, ": reuse inherited-omic PGS scans")
      return(old)
    }
  }
  empty <- list(signature = signature, score_file = score_file,
    status = tibble(status = "unavailable", detail = if (is.na(score_file))
      "prot.pgs.rds/met.pgs.rds was not found" else "No assayed feature matched a PGS column"),
    incident = empty_c1_pgs_assoc(features), prevalent = empty_c1_pgs_assoc(features),
    attained_age = empty_c1_pgs_assoc(features),
    incident_same_omic = empty_c1_pgs_assoc(features),
    prevalent_same_omic = empty_c1_pgs_assoc(features),
    attained_age_same_omic = empty_c1_pgs_assoc(features), score_map = character())
  if (is.na(score_file)) {
    warning("C1/", layer, ": inherited-score file is unavailable; PGS panels will be blank.", call. = FALSE)
    return(empty)
  }
  scores <- read_c1_pgs(score_file)
  score_map <- map_c1_pgs_columns(features, names(scores))
  max_scores <- suppressWarnings(as.integer(Sys.getenv("C1_PGS_MAX", unset = "0")))
  if (is.finite(max_scores) && max_scores > 0) score_map <- head(score_map, max_scores)
  if (!length(score_map)) return(empty)

  need <- unique(c("eid", "ethnic.c", covars, "birth_date", "date_attend", "date_lost",
    "date_death", paste0("fod_icd10_", outcome)))
  ph <- read_all(need) |> filter_analysis_cohort() |> make_outcome(outcome) |>
    add_attained_age_time(outcome)
  ph$prevalent_status <- make_prevalent_status(ph, outcome)
  # Do not remove this conversion; see read_c1_pgs().
  ph$eid <- as.character(ph$eid)
  idx <- match(scores$eid, ph$eid)
  keep <- !is.na(idx)
  scores <- scores[keep, unique(c("eid", unname(score_map))), drop = FALSE]
  base <- ph[idx[keep], , drop = FALSE]
  base$.omic_overlap <- if (length(overlap_eids)) base$eid %in% overlap_eids else TRUE
  rm(ph); invisible(gc())
  tvar <- paste0(outcome, ".t2e"); evar <- paste0(outcome, ".Yt2e")
  covars <- intersect(covars, names(base))

  message("C1/", layer, ": scan ", length(score_map),
    " inherited omic scores in the genotyped cohort (sequential, memory-safe)")
  rows <- lapply(names(score_map), function(feature) {
    gcol <- score_map[[feature]]
    d <- base
    d$.pgs_score <- suppressWarnings(as.numeric(scores[[gcol]]))
    inc <- cox_scan(d, ".pgs_score", covars, outcome, time_var = tvar, event_var = evar) |>
      mutate(term = feature, score_column = gcol)
    prev <- logistic_scan(d, ".pgs_score", covars, "prevalent_status") |>
      mutate(term = feature, score_column = gcol)
    age <- cox_scan_delayed_entry(d, ".pgs_score", setdiff(covars,
      unique(c("age", grep("^age($|[._])", covars, value = TRUE, ignore.case = TRUE)))), outcome) |>
      mutate(term = feature, score_column = gcol)
    ds <- d[d$.omic_overlap %in% TRUE, , drop = FALSE]
    inc_same <- cox_scan(ds, ".pgs_score", covars, outcome, time_var = tvar, event_var = evar) |>
      mutate(term = feature, score_column = gcol)
    prev_same <- logistic_scan(ds, ".pgs_score", covars, "prevalent_status") |>
      mutate(term = feature, score_column = gcol)
    age_same <- cox_scan_delayed_entry(ds, ".pgs_score", setdiff(covars,
      unique(c("age", grep("^age($|[._])", covars, value = TRUE, ignore.case = TRUE)))), outcome) |>
      mutate(term = feature, score_column = gcol)
    list(incident = inc, prevalent = prev, attained_age = age,
      incident_same_omic = inc_same, prevalent_same_omic = prev_same,
      attained_age_same_omic = age_same)
  })
  finish <- function(kind) bind_rows(lapply(rows, `[[`, kind)) |>
    mutate(FDR = p.adjust(p.value, "BH")) |> arrange(p.value)
  ans <- list(signature = signature, score_file = score_file,
    status = tibble(status = "ok", detail = paste(length(score_map), "matched PGS columns"),
      score_file = score_file, genotype_rows = nrow(scores), phenotype_matches = nrow(base),
      same_omic_rows = sum(base$.omic_overlap),
      adjustment = "basic covariates; LE8 deliberately not included as a downstream exposure"),
    incident = finish("incident"), prevalent = finish("prevalent"),
    attained_age = finish("attained_age"),
    incident_same_omic = finish("incident_same_omic"),
    prevalent_same_omic = finish("prevalent_same_omic"),
    attained_age_same_omic = finish("attained_age_same_omic"),score_map = score_map)
  if (!is.na(cache)) saveRDS(ans, cache, compress = "xz")
  ans
}

build_pgs_actual_concordance <- function(pgs, observed) {
  analyses <- c(incident = "Incident", prevalent = "Baseline prevalent",
    attained_age = "Attained-age delayed entry")
  map_dfr(names(analyses), function(nm) {
    g <- as_tibble(pgs[[nm]] %||% tibble())
    o <- as_tibble(observed[[nm]] %||% tibble())
    if (!nrow(g) && !nrow(o)) return(tibble())
    gg <- if (nrow(g)) g |> transmute(feature = term, score_column,
      pgs_beta = beta, pgs_p = p.value, pgs_FDR = FDR, pgs_N = N_total,
      pgs_events = N_event) else tibble(feature = character())
    oo <- if (nrow(o)) o |> transmute(feature = term, observed_beta = beta,
      observed_p = p.value, observed_FDR = FDR, observed_N = N_total,
      observed_events = N_event) else tibble(feature = character())
    full_join(gg, oo, by = "feature") |> mutate(analysis = analyses[[nm]],
      pgs_supported = is.finite(pgs_FDR) & pgs_FDR < .05,
      observed_supported = is.finite(observed_FDR) & observed_FDR < .05,
      sign_concordant = is.finite(pgs_beta) & is.finite(observed_beta) &
        sign(pgs_beta) == sign(observed_beta),
      evidence_pattern = case_when(
        pgs_supported & observed_supported & sign_concordant ~ "Both, same direction",
        pgs_supported & observed_supported & !sign_concordant ~ "Both, opposite direction",
        pgs_supported & !observed_supported ~ "PGS only",
        !pgs_supported & observed_supported ~ "Observed only",
        TRUE ~ "Neither at FDR 5%"))
  })
}

plot_pgs_actual_concordance <- function(x, anchors = character()) {
  d <- as_tibble(x) |> filter(is.finite(pgs_beta), is.finite(observed_beta))
  if (!nrow(d)) return(blank_plot("Inherited PGS versus measured adult omic level",
    "No matched PGS–observed association pair was estimable"))
  lab <- d |> group_by(analysis) |>
    mutate(.discord = abs(scale(pgs_beta)[,1] - scale(observed_beta)[,1]),
      .label = feature %in% anchors | min_rank(desc(.discord)) <= 6) |>
    ungroup() |> mutate(label = ifelse(.label, feature, NA_character_))
  cols <- c("Both, same direction" = "#1B9E77", "Both, opposite direction" = "#D7301F",
    "PGS only" = "#6A3D9A", "Observed only" = "#D95F02", "Neither at FDR 5%" = "grey72")
  pa <- ggplot(lab, aes(pgs_beta, observed_beta, color = evidence_pattern)) +
    geom_hline(yintercept = 0, color = "grey82") + geom_vline(xintercept = 0, color = "grey82") +
    geom_point(alpha = .62, size = 1.7) +
    ggrepel::geom_text_repel(aes(label = label), size = 2.3, max.overlaps = 24, seed = 91) +
    facet_wrap(~analysis, scales = "free", nrow = 1) + scale_color_manual(values = cols) +
    labs(title = "a. Inherited score and adult measured level are distinct signals",
      subtitle = "Both predictors are standardized within each model; PGS is fixed at conception but is not a birth protein measurement",
      x = "Log HR/OR per 1-SD omic PGS", y = "Log HR/OR per 1-SD measured omic",
      color = NULL) + theme_5c(8) + theme(legend.position = "bottom")
  ct <- x |> count(analysis, evidence_pattern, name = "features")
  pb <- ggplot(ct, aes(features, fct_reorder(evidence_pattern, features, sum), fill = evidence_pattern)) +
    geom_col() + facet_wrap(~analysis, scales = "free_x", nrow = 1) +
    geom_text(aes(label = features), hjust = -.12, fontface = "bold", size = 2.7) +
    scale_fill_manual(values = cols, guide = "none") +
    scale_x_continuous(expand = expansion(mult = c(0, .18))) +
    labs(title = "b. FDR-supported evidence patterns", x = "Features", y = NULL) + theme_5c(8)
  pa / pb + plot_layout(heights = c(1.4, .75))
}
