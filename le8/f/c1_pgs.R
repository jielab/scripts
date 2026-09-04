# C1 inherited-score scan utilities.
#
# A COJO PGS is fixed at conception, but it is not the biomarker concentration
# "measured at birth".  These functions therefore keep PGS and adult measured
# omics in parallel columns and never substitute one for the other.

C1_PGS_SCAN_VERSION <- "2026-09-03.vldl_conditional1"

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
  signature <- paste(C1_PGS_SCAN_VERSION,c1_pgs_signature(layer), "omic_overlap",
    overlap_signature, sep = "|")
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
    attained_age_same_omic = empty_c1_pgs_assoc(features),
    vldl_conditional=tibble(),score_map = character())
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
  vldl_conditional<-tibble()
  if(layer=="metabolite"&&"L_VLDL_TG.pct"%in%names(score_map)){
    comparator_features<-intersect(c("L_VLDL_TG","Total_TG","ApoB","VLDL_size"),
      names(score_map))
    model_sets<-c(list(`Target PGS only`=character()),
      setNames(lapply(comparator_features,c),paste0("+ ",comparator_features)))
    if(length(comparator_features)>1L)
      model_sets[["+ all burden/size PGS"]]<-comparator_features
    dd<-base
    target_name<-".pgs_L_VLDL_TG_pct"
    dd[[target_name]]<-suppressWarnings(as.numeric(scores[[score_map[["L_VLDL_TG.pct"]]]]))
    zmap<-setNames(character(length(comparator_features)),comparator_features)
    for(feature in comparator_features){
      zn<-paste0(".pgs_",make.names(feature));zmap[[feature]]<-zn
      dd[[zn]]<-suppressWarnings(as.numeric(scores[[score_map[[feature]]]]))
    }
    same_n_need<-unique(c(tvar,evar,covars,target_name,unname(zmap)))
    dd<-dd[complete.cases(dd[,same_n_need,drop=FALSE]),,drop=FALSE]
    for(nm in c(target_name,unname(zmap)))dd[[nm]]<-as.numeric(scale(dd[[nm]]))
    vldl_conditional<-map_dfr(names(model_sets),function(model_name){
      adjust_features<-model_sets[[model_name]];xs<-c(target_name,unname(zmap[adjust_features]))
      empty_row<-tibble(model=model_name,target="L_VLDL_TG.pct PGS",adjusters=
        paste(adjust_features,collapse=";"),beta=NA_real_,std.error=NA_real_,conf.low=NA_real_,
        conf.high=NA_real_,p.value=NA_real_,N_total=nrow(dd),N_event=sum(dd[[evar]]==1),
        max_abs_score_correlation=NA_real_,exposure_condition_number=NA_real_)
      if(nrow(dd)<500||sum(dd[[evar]]==1)<20)return(empty_row)
      ff<-as.formula(paste0("Surv(",bt(tvar),",",bt(evar),") ~ ",
        paste(bt(c(xs,covars)),collapse=" + ")))
      fit<-tryCatch(coxph(ff,dd,ties="efron"),error=function(e)NULL)
      if(is.null(fit))return(empty_row)
      sm<-coef(summary(fit));if(!target_name%in%rownames(sm))return(empty_row)
      b<-sm[target_name,"coef"];se<-sm[target_name,"se(coef)"]
      cm<-if(length(xs)>1L)cor(dd[,xs,drop=FALSE])else matrix(1,1,1)
      max_cor<-if(length(xs)>1L)max(abs(cm[1,-1]),na.rm=TRUE)else 0
      cond<-if(length(xs)>1L)tryCatch(kappa(cm),error=function(e)NA_real_)else 1
      tibble(model=model_name,target="L_VLDL_TG.pct PGS",adjusters=paste(adjust_features,collapse=";"),
        beta=b,std.error=se,conf.low=b-1.96*se,conf.high=b+1.96*se,
        p.value=sm[target_name,"Pr(>|z|)"],N_total=nrow(dd),N_event=sum(dd[[evar]]==1),
        max_abs_score_correlation=max_cor,exposure_condition_number=cond)
    })|>mutate(FDR=p.adjust(p.value,"BH"),
      interpretation="Exploratory multivariable PGS association; not multivariable MR")
  }
  ans <- list(signature = signature, score_file = score_file,
    status = tibble(status = "ok", detail = paste(length(score_map), "matched PGS columns"),
      score_file = score_file, genotype_rows = nrow(scores), phenotype_matches = nrow(base),
      same_omic_rows = sum(base$.omic_overlap),
      adjustment = "basic covariates; LE8 deliberately not included as a downstream exposure"),
    incident = finish("incident"), prevalent = finish("prevalent"),
    attained_age = finish("attained_age"),
    incident_same_omic = finish("incident_same_omic"),
    prevalent_same_omic = finish("prevalent_same_omic"),
    attained_age_same_omic = finish("attained_age_same_omic"),
    vldl_conditional=vldl_conditional,score_map = score_map)
  if (!is.na(cache)) saveRDS(ans, cache, compress = "xz")
  ans
}

build_pgs_actual_concordance <- function(pgs_full, observed, pgs_same = NULL) {
  analyses <- c(incident = "Incident", prevalent = "Baseline prevalent",
    attained_age = "Attained-age delayed entry")
  map_dfr(names(analyses), function(nm) {
    g <- as_tibble(pgs_full[[nm]] %||% tibble())
    s <- as_tibble((pgs_same %||% list())[[nm]] %||% tibble())
    o <- as_tibble(observed[[nm]] %||% tibble())
    if (!nrow(g) && !nrow(s) && !nrow(o)) return(tibble())
    if (nrow(g) && !"score_column" %in% names(g)) g$score_column <- NA_character_
    if (nrow(s) && !"score_column" %in% names(s)) s$score_column <- NA_character_
    gg <- if (nrow(g)) g |> transmute(feature = term, score_column,
      pgs_full_beta = beta, pgs_full_se = std.error, pgs_full_p = p.value,
      pgs_full_FDR = FDR, pgs_full_N = N_total, pgs_full_events = N_event) else
      tibble(feature = character(),score_column=character(),pgs_full_beta=double(),
        pgs_full_se=double(),pgs_full_p=double(),pgs_full_FDR=double(),
        pgs_full_N=double(),pgs_full_events=double())
    ss <- if (nrow(s)) s |> transmute(feature = term,
      pgs_same_beta = beta, pgs_same_se = std.error, pgs_same_p = p.value,
      pgs_same_FDR = FDR, pgs_same_N = N_total, pgs_same_events = N_event) else
      tibble(feature = character(),pgs_same_beta=double(),pgs_same_se=double(),
        pgs_same_p=double(),pgs_same_FDR=double(),pgs_same_N=double(),
        pgs_same_events=double())
    oo <- if (nrow(o)) o |> transmute(feature = term, observed_beta = beta,
      observed_se = std.error, observed_p = p.value, observed_FDR = FDR,
      observed_N = N_total, observed_events = N_event) else
      tibble(feature=character(),observed_beta=double(),observed_se=double(),
        observed_p=double(),observed_FDR=double(),observed_N=double(),
        observed_events=double())
    full_join(gg, ss, by = "feature") |> full_join(oo, by = "feature") |>
      mutate(analysis = analyses[[nm]],
      # Compatibility aliases now refer explicitly to the full genetic cohort.
      pgs_beta = pgs_full_beta, pgs_p = pgs_full_p, pgs_FDR = pgs_full_FDR,
      pgs_N = pgs_full_N, pgs_events = pgs_full_events,
      pgs_supported = is.finite(pgs_full_FDR) & pgs_full_FDR < .05,
      observed_supported = is.finite(observed_FDR) & observed_FDR < .05,
      sign_concordant = is.finite(pgs_full_beta) & is.finite(observed_beta) &
        sign(pgs_full_beta) == sign(observed_beta),
      same_omic_sign_concordant = is.finite(pgs_same_beta) & is.finite(observed_beta) &
        sign(pgs_same_beta) == sign(observed_beta),
      beta_difference_full = observed_beta - pgs_full_beta,
      beta_difference_same_omic = observed_beta - pgs_same_beta,
      evidence_pattern = case_when(
        pgs_supported & observed_supported & sign_concordant ~ "Both, same direction",
        pgs_supported & observed_supported & !sign_concordant ~ "Both, opposite direction",
        pgs_supported & !observed_supported ~ "PGS only",
        !pgs_supported & observed_supported ~ "Observed only",
        TRUE ~ "Neither at FDR 5%"))
  })
}

plot_pgs_actual_concordance <- function(x, anchors = character()) {
  d <- as_tibble(x)
  fmt_n <- function(v) {
    n <- sort(unique(as.integer(v[is.finite(v)])))
    if (!length(n)) return("NA")
    z <- if (length(n) == 1L) n else range(n)
    paste(format(z, big.mark = ",", scientific = FALSE, trim = TRUE), collapse = "–")
  }
  one_panel <- function(xvar, nvar, panel_title, xlab) {
    z <- d |> filter(is.finite(.data[[xvar]]), is.finite(observed_beta))
    if (!nrow(z)) return(blank_plot(panel_title, "No matched PGS–observed effect pair was estimable"))
    facets <- z |> group_by(analysis) |> summarise(
      facet = paste0(first(analysis), "\nPGS N = ", fmt_n(.data[[nvar]]),
        "; measured N = ", fmt_n(observed_N)), .groups = "drop")
    z <- z |> left_join(facets, by = "analysis") |> group_by(analysis) |>
      mutate(.discord = abs(observed_beta - .data[[xvar]]),
        label = ifelse(feature %in% anchors | min_rank(desc(.discord)) <= 6, feature, NA_character_),
        direction = ifelse(sign(observed_beta) == sign(.data[[xvar]]),
          "Same direction", "Opposite direction")) |> ungroup()
    st <- z |> group_by(facet) |> summarise(
      rho = suppressWarnings(cor(.data[[xvar]], observed_beta, method = "spearman", use = "complete.obs")),
      stat_label = ifelse(is.finite(rho), sprintf("Spearman rho = %.2f", rho), "Spearman rho = NA"),
      .groups = "drop")
    ggplot(z, aes(x = .data[[xvar]], y = observed_beta, color = direction)) +
      geom_abline(slope = 1, intercept = 0, linetype = 2, color = "grey58") +
      geom_hline(yintercept = 0, color = "grey82") + geom_vline(xintercept = 0, color = "grey82") +
      geom_point(alpha = .62, size = 1.65) +
      geom_text(data = st, aes(x = -Inf, y = Inf, label = stat_label), inherit.aes = FALSE,
        hjust = -.04, vjust = 1.25, size = 2.5, fontface = "bold", color = "grey30") +
      ggrepel::geom_text_repel(aes(label = label), size = 2.25, max.overlaps = 24,
        seed = 91, show.legend = FALSE) +
      facet_wrap(~facet, scales = "free", nrow = 1) +
      scale_color_manual(values = c("Same direction" = "#1B9E77",
        "Opposite direction" = "#D7301F")) +
      labs(title = panel_title,
        subtitle = "Effects, not P values, are compared; both predictors are standardized within their own model",
        x = xlab, y = "Measured adult omic: log HR/OR per 1 SD", color = NULL) +
      theme_5c(8) + theme(legend.position = "bottom")
  }
  pa <- one_panel("pgs_full_beta", "pgs_full_N",
    "a. Full genetic cohort: inherited versus measured effect size",
    "Inherited omic PGS: log HR/OR per 1 SD")
  pb <- one_panel("pgs_same_beta", "pgs_same_N",
    "b. Same-omic-subset sensitivity",
    "Inherited omic PGS in assay subset: log HR/OR per 1 SD")
  pa / pb + plot_layout(heights = c(1, 1)) +
    plot_annotation(caption = paste0(
      "PGS is fixed at conception but is not a biomarker measured at birth. ",
      "Effect concordance is descriptive and is not an MR or mediation estimate."))
}
