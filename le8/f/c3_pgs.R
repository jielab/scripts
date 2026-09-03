# C3 feature-level triangulation of PGS, observed level and locus colocalization.
# PGS evidence is deliberately kept out of coloc.abf(): a genome-wide score is
# not a locus-specific QTL summary statistic and cannot supply PP(H4).

make_c3_pgs_integration <- function(c1, c2, coloc_summary, h4 = .70) {
  actual0 <- as_tibble(c1$association_adj2 %||% c1$association %||% tibble())
  actual <- if (nrow(actual0)) actual0 |> transmute(feature = term,
    observed_beta = beta, observed_p = p.value, observed_FDR = FDR) else
    tibble(feature = character())
  pgs0 <- as_tibble(c1$pgs_incident %||% tibble())
  pgs <- if (nrow(pgs0)) pgs0 |> transmute(feature = term, score_column,
    pgs_beta = beta, pgs_p = p.value, pgs_FDR = FDR, pgs_N = N_total) else
    tibble(feature = character())
  pgs_same0 <- as_tibble(c1$pgs_incident_same_omic %||% tibble())
  pgs_same <- if (nrow(pgs_same0)) pgs_same0 |> transmute(feature = term,
    pgs_same_omic_beta = beta, pgs_same_omic_p = p.value,
    pgs_same_omic_FDR = FDR, pgs_same_omic_N = N_total) else tibble(feature = character())
  co0 <- as_tibble(coloc_summary %||% tibble())
  if (nrow(co0) && !"PP.H4_robust_min" %in% names(co0)) co0$PP.H4_robust_min <- co0$PP.H4
  co <- if (nrow(co0) && all(c("feature", "PP.H4") %in% names(co0))) co0 |>
    filter(status == "ok") |> mutate(.pp4 = coalesce(PP.H4_robust_min, PP.H4)) |>
    group_by(feature) |> slice_max(.pp4, n = 1, with_ties = FALSE) |> ungroup() |>
    transmute(feature, locus, PP.H4, PP.H4_robust_min, robust_PP4 = .pp4,
      credible_set_n, lead_shared) else tibble(feature = character())
  de0 <- as_tibble((c2$individual_decomposition$summary %||%
    c2$individual_genetic_decomposition %||% c2$genetic_decomposition$summary %||%
    c2$PGS_decomposition$summary %||% tibble()))
  if (nrow(de0) && !"cross_fitted" %in% names(de0)) de0$cross_fitted <- FALSE
  de_defaults <- c("genetic_partial_R2", "incident_beta_genetic", "incident_p_genetic",
    "incident_beta_residual", "incident_p_residual")
  if (nrow(de0)) for (nm in de_defaults) if (!nm %in% names(de0)) de0[[nm]] <- NA_real_
  de <- if (nrow(de0) && "feature" %in% names(de0)) de0 |>
    mutate(residual_FDR = p.adjust(incident_p_residual, "BH")) |>
    transmute(feature, genetic_partial_R2, incident_beta_genetic,
      incident_p_genetic, incident_beta_residual, incident_p_residual, residual_FDR,
      decomposition_cross_fitted = coalesce(cross_fitted, FALSE)) else
    tibble(feature = character())
  z <- Reduce(function(x, y) full_join(x, y, by = "feature"), list(actual, pgs, pgs_same, co, de))
  if (!nrow(z)) return(tibble())
  defaults <- list(observed_beta = NA_real_, observed_p = NA_real_, observed_FDR = NA_real_,
    pgs_beta = NA_real_, pgs_p = NA_real_, pgs_FDR = NA_real_, robust_PP4 = NA_real_,
    pgs_same_omic_beta = NA_real_, pgs_same_omic_p = NA_real_, pgs_same_omic_FDR = NA_real_,
    residual_FDR = NA_real_, incident_beta_residual = NA_real_)
  for (nm in names(defaults)) if (!nm %in% names(z)) z[[nm]] <- defaults[[nm]]
  z |> mutate(observed_supported = is.finite(observed_FDR) & observed_FDR < .05,
    pgs_supported = is.finite(pgs_FDR) & pgs_FDR < .05,
    coloc_supported = is.finite(robust_PP4) & robust_PP4 >= h4,
    residual_supported = is.finite(residual_FDR) & residual_FDR < .05,
    pgs_observed_sign_match = is.finite(pgs_beta) & is.finite(observed_beta) &
      sign(pgs_beta) == sign(observed_beta),
    pgs_same_observed_sign_match = is.finite(pgs_same_omic_beta) & is.finite(observed_beta) &
      sign(pgs_same_omic_beta) == sign(observed_beta),
    inherited_locus_pattern = pgs_supported & coloc_supported & observed_supported &
      pgs_observed_sign_match,
    reactive_compatible_pattern = observed_supported & residual_supported &
      !pgs_supported & !coloc_supported,
    triangulation_class = case_when(
      inherited_locus_pattern ~ "PGS + observed + coloc",
      observed_supported & coloc_supported ~ "Observed + coloc; PGS weak",
      observed_supported & pgs_supported & pgs_observed_sign_match ~ "PGS + observed; no robust coloc",
      pgs_supported & !observed_supported ~ "PGS only",
      reactive_compatible_pattern ~ "Observed/residual only; reactive-compatible",
      observed_supported ~ "Observed only",
      TRUE ~ "Insufficient evidence"),
    interpretation = "Feature-level PGS is orthogonal support only; PP(H4) remains locus-specific") |>
    arrange(desc(inherited_locus_pattern), desc(coloc_supported), pgs_p, observed_p)
}

read_c3_pgs_integration <- function(layer, outdir, coloc_summary) {
  c1f <- file.path(le8_job_dir(outdir, "c1_correlate"), "c1.res.rds")
  c2f <- file.path(le8_job_dir(outdir, "c2_cause"), "c2.res.rds")
  c1 <- if (file.exists(c1f)) tryCatch(readRDS(c1f), error = function(e) list()) else list()
  c2 <- if (file.exists(c2f)) tryCatch(readRDS(c2f), error = function(e) list()) else list()
  make_c3_pgs_integration(c1, c2, coloc_summary, H4_STRONG)
}

plot_c3_pgs_integration <- function(x, outdir) {
  d <- as_tibble(x)
  if (!nrow(d)) {
    save_plot(blank_plot("PGS–observed–colocalization triangulation",
      "C1 inherited-score results or C3 locus results were unavailable"),
      "c3.Fig7.pgs_coloc_triangulation.png", 11, 7, outdir = outdir)
    return(invisible(NULL))
  }
  cols <- c("PGS + observed + coloc" = "#1B9E77",
    "Observed + coloc; PGS weak" = "#66A61E",
    "PGS + observed; no robust coloc" = "#6A3D9A", "PGS only" = "#377EB8",
    "Observed/residual only; reactive-compatible" = "#D7301F",
    "Observed only" = "#D95F02", "Insufficient evidence" = "grey72")
  xy <- d |> mutate(pgs_plot_beta = coalesce(pgs_same_omic_beta, pgs_beta)) |>
    filter(is.finite(pgs_plot_beta), is.finite(observed_beta)) |>
    mutate(label = ifelse(inherited_locus_pattern | reactive_compatible_pattern |
      min_rank(pmin(pgs_p, observed_p, na.rm = TRUE)) <= 12, feature, NA_character_))
  pa <- if (!nrow(xy)) blank_plot("a. PGS versus measured level") else
    ggplot(xy, aes(pgs_plot_beta, observed_beta, color = triangulation_class)) +
      geom_hline(yintercept = 0, color = "grey82") + geom_vline(xintercept = 0, color = "grey82") +
      geom_point(alpha = .7, size = 1.8) +
      ggrepel::geom_text_repel(aes(label = label), size = 2.4, max.overlaps = 22, seed = 93) +
      scale_color_manual(values = cols) +
      labs(title = "a. Feature-level inherited versus measured signal",
        subtitle = "Same-omic-sample PGS estimate is used when available; full genetic cohort drives PGS significance",
        x = "Log HR per 1-SD omic PGS", y = "Log HR per 1-SD measured omic", color = NULL) +
      theme_5c(8) + theme(legend.position = "bottom")
  cp <- d |> filter(is.finite(robust_PP4), is.finite(pgs_p)) |>
    mutate(pgs_logp = -log10(pmax(pgs_p, 1e-300)),
      label = ifelse(coloc_supported | pgs_supported, feature, NA_character_))
  pb <- if (!nrow(cp)) blank_plot("b. PGS evidence versus robust PP(H4)",
    "No feature had both PGS and colocalization results") else
    ggplot(cp, aes(robust_PP4, pgs_logp, color = triangulation_class)) +
      geom_vline(xintercept = H4_STRONG, linetype = 2, color = "grey55") +
      geom_hline(yintercept = -log10(.05 / max(1, nrow(cp))), linetype = 3, color = "grey65") +
      geom_point(alpha = .75, size = 2) +
      ggrepel::geom_text_repel(aes(label = label), size = 2.4, max.overlaps = 22, seed = 94) +
      scale_color_manual(values = cols) +
      labs(title = "b. Orthogonal genetic support",
        subtitle = "PGS is not entered into coloc.abf; horizontal line is a displayed-feature Bonferroni guide",
        x = "Robust minimum PP(H4)", y = expression(-log[10](P[PGS])), color = NULL) +
      theme_5c(8) + theme(legend.position = "bottom")
  ct <- d |> count(triangulation_class, name = "features")
  pc <- ggplot(ct, aes(features, fct_reorder(triangulation_class, features), fill = triangulation_class)) +
    geom_col() + geom_text(aes(label = features), hjust = -.12, fontface = "bold") +
    scale_fill_manual(values = cols, guide = "none") +
    scale_x_continuous(expand = expansion(mult = c(0, .18))) +
    labs(title = "c. Triangulation audit", x = "Features", y = NULL) + theme_5c(8)
  save_plot((pa | pb) / pc + plot_layout(heights = c(1.35, .65)),
    "c3.Fig7.pgs_coloc_triangulation.png", 17, 11, outdir = outdir)
}
