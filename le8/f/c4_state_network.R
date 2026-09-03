# State-dependent omic co-abundance remodeling.
#
# This adapts the experimental paper's "reference network -> differential
# edges -> convergent hubs" logic without claiming that plasma covariance is a
# physical protein-protein interaction. A known-PPI file can annotate edges,
# but AP-MS/variant perturbation evidence is required before calling them PPI
# rewiring.

C4_NETWORK_TOP <- as.integer(Sys.getenv("C4_NETWORK_TOP", unset = "80"))
C4_NETWORK_DELTA <- as.numeric(Sys.getenv("C4_NETWORK_DELTA", unset = "0.20"))
C4_NETWORK_MIN_N <- as.integer(Sys.getenv("C4_NETWORK_MIN_N", unset = "250"))
C4_NETWORK_ANCHORS <- unique(trimws(strsplit(Sys.getenv("C1_DIRECTION_ANCHORS",
  unset = "PCSK9,LPA,GDF15,NTPROBNP,MMP12"), ",", fixed = TRUE)[[1]]))

read_c4_ppi_reference <- function() {
  f <- Sys.getenv("C4_PPI_FILE", unset = "")
  if (!nzchar(f) || !file.exists(f) || file.size(f) <= 0) return(tibble())
  x <- tryCatch(as_tibble(data.table::fread(f, showProgress = FALSE,
    check.names = FALSE)), error = function(e) tibble())
  if (ncol(x) < 2) return(tibble())
  names(x)[1:2] <- c("feature1", "feature2")
  x |> transmute(feature1 = as.character(feature1), feature2 = as.character(feature2),
    ppi_key = ifelse(feature1 < feature2, paste(feature1, feature2, sep = "||"),
      paste(feature2, feature1, sep = "||"))) |> distinct(ppi_key, .keep_all = TRUE)
}

select_c4_network_features <- function(disease, sets, modules, available,
  max_n = C4_NETWORK_TOP) {
  ranked <- if (nrow(disease) && all(c("term", "p.value") %in% names(disease))) disease |>
    filter(is.finite(p.value)) |> arrange(p.value) |> pull(term) else character()
  supervised <- unique(c(sets$YS_strict %||% character(), sets$YS %||% character(),
    modules$YS_core %||% character(), as.character(sets$membership$feature %||% character())))
  preferred <- unique(c(C4_NETWORK_ANCHORS, supervised, ranked))
  head(preferred[preferred %in% available], max_n)
}

c4_residual_matrix <- function(d, features, covars) {
  covars <- intersect(covars, names(d)); features <- intersect(features, names(d))
  if (!length(features)) return(list(x = matrix(numeric(), 0, 0), n = 0L))
  if (length(covars)) d <- d[complete.cases(d[, covars, drop = FALSE]), , drop = FALSE]
  if (nrow(d) < C4_NETWORK_MIN_N) return(list(x = matrix(numeric(), 0, 0), n = nrow(d)))
  x <- as.matrix(data.frame(lapply(d[, features, drop = FALSE], function(z) {
    z <- suppressWarnings(as.numeric(z)); med <- median(z, na.rm = TRUE)
    if (!is.finite(med)) med <- 0
    z[!is.finite(z)] <- med; z
  }), check.names = FALSE))
  colnames(x) <- features; keep <- apply(x, 2, function(z) is.finite(sd(z)) && sd(z) > 0)
  x <- x[, keep, drop = FALSE]
  if (!ncol(x)) return(list(x = x, n = nrow(d)))
  if (length(covars)) {
    mm <- model.matrix(reformulate(covars), d)
    x <- qr.resid(qr(mm), x)
  }
  x <- scale(x); x[!is.finite(x)] <- 0
  list(x = x, n = nrow(x))
}

c4_corr_edges <- function(ref, alt, comparison, ppi = tibble()) {
  common <- intersect(colnames(ref$x), colnames(alt$x))
  if (length(common) < 2 || ref$n < 4 || alt$n < 4) return(tibble())
  r0 <- cor(ref$x[, common, drop = FALSE], use = "pairwise.complete.obs")
  r1 <- cor(alt$x[, common, drop = FALSE], use = "pairwise.complete.obs")
  ij <- which(upper.tri(r0), arr.ind = TRUE)
  ans <- tibble(feature1 = common[ij[,1]], feature2 = common[ij[,2]],
    reference_r = r0[ij], state_r = r1[ij], delta_r = r1[ij] - r0[ij],
    reference_N = ref$n, state_N = alt$n, comparison = comparison) |>
    mutate(se_delta_z = sqrt(1 / pmax(reference_N - 3, 1) + 1 / pmax(state_N - 3, 1)),
      fisher_delta = atanh(cap(state_r, .999)) - atanh(cap(reference_r, .999)),
      z_diff = fisher_delta / se_delta_z,
      p.value = 2 * pnorm(abs(z_diff), lower.tail = FALSE), FDR = p.adjust(p.value, "BH"),
      ci_low = fisher_delta - 1.96 * se_delta_z, ci_high = fisher_delta + 1.96 * se_delta_z,
      edge_key = ifelse(feature1 < feature2, paste(feature1, feature2, sep = "||"),
        paste(feature2, feature1, sep = "||")),
      remodeling_class = case_when(
        sign(reference_r) != sign(state_r) & abs(reference_r) >= .15 & abs(state_r) >= .15 ~ "sign reversal",
        abs(reference_r) < .15 & abs(state_r) >= .30 ~ "gained co-abundance",
        abs(reference_r) >= .30 & abs(state_r) < .15 ~ "lost co-abundance",
        abs(state_r) > abs(reference_r) ~ "strengthened",
        TRUE ~ "weakened"),
      remodeled = FDR < .05 & abs(delta_r) >= C4_NETWORK_DELTA)
  if (nrow(ppi)) ans <- ans |> left_join(ppi |> transmute(edge_key = ppi_key,
    known_PPI_backbone = TRUE), by = "edge_key") |>
    mutate(known_PPI_backbone = coalesce(known_PPI_backbone, FALSE))
  else ans$known_PPI_backbone <- FALSE
  ans |> arrange(FDR, desc(abs(delta_r)))
}

make_c4_state_hubs <- function(edges, sets, modules) {
  if (!nrow(edges)) return(tibble())
  hubs <- bind_rows(edges |> transmute(comparison, feature = feature1, remodeled,
      delta_strength = abs(delta_r), known_PPI_backbone),
    edges |> transmute(comparison, feature = feature2, remodeled,
      delta_strength = abs(delta_r), known_PPI_backbone)) |>
    group_by(comparison, feature) |>
    summarise(remodeled_edges = sum(remodeled, na.rm = TRUE),
      remodeling_burden = sum(delta_strength[remodeled], na.rm = TRUE),
      known_PPI_edges = sum(remodeled & known_PPI_backbone, na.rm = TRUE), .groups = "drop")
  sm <- as_tibble(sets$membership %||% tibble())
  if (nrow(sm)) sm <- sm |> select(any_of(c("feature", "set", "primary_component", "strict_YS")))
  mm <- as_tibble(modules$membership %||% tibble())
  if (nrow(mm)) mm <- mm |> select(any_of(c("feature", "module", "module_stability", "YS_core")))
  if (nrow(sm)) hubs <- hubs |> left_join(sm, by = "feature")
  if (nrow(mm)) hubs <- hubs |> left_join(mm, by = "feature")
  hubs |> arrange(comparison, desc(remodeled_edges), desc(remodeling_burden))
}

plot_c4_state_network <- function(net, layer, outdir) {
  edges <- as_tibble(net$edges %||% tibble()); hubs <- as_tibble(net$hubs %||% tibble())
  if (!nrow(edges)) {
    msg <- as.character((net$status$detail %||% "No estimable state comparison")[[1]])
    save_plot(blank_plot("State-dependent omic co-abundance remodeling", msg),
      "c4.Fig10.state_network_remodeling.png", 11, 7, outdir = outdir)
    save_plot(blank_plot("Differential co-abundance edge audit", msg),
      "c4.Fig11.state_network_edges.png", 11, 7, outdir = outdir)
    return(invisible(NULL))
  }
  top_edges <- edges |> filter(remodeled) |> group_by(comparison) |>
    slice_min(FDR, n = 70, with_ties = FALSE) |> ungroup()
  pa <- if (!nrow(top_edges)) blank_plot("a. Differential co-abundance edges",
    "No edge passed both FDR and effect-size thresholds") else
    ggplot(top_edges, aes(feature1, feature2, fill = delta_r, size = -log10(pmax(FDR, 1e-30)))) +
      geom_point(shape = 21, color = "grey45") + facet_wrap(~comparison, scales = "free") +
      scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0) +
      labs(title = "a. State-dependent co-abundance remodeling",
        subtitle = "Residual correlation changes relative to participants disease-free through year 10",
        x = NULL, y = NULL, fill = expression(Delta*r), size = expression(-log[10](FDR))) +
      theme_5c(7) + theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "bottom")
  ht <- hubs |> group_by(comparison) |> slice_max(remodeled_edges, n = 15, with_ties = FALSE) |>
    ungroup() |> filter(remodeled_edges > 0) |>
    mutate(label = paste0(feature, " (", remodeled_edges, ")"))
  pb <- if (!nrow(ht)) blank_plot("b. Remodeling hubs") else
    ggplot(ht, aes(remodeling_burden, fct_reorder(label, remodeling_burden), fill = comparison)) +
      geom_col() + facet_wrap(~comparison, scales = "free_y") +
      labs(title = "b. Convergent remodeling hubs", x = expression(Sigma*abs(Delta*r)), y = NULL,
        fill = NULL) + theme_5c(8) + theme(legend.position = "none")
  save_plot(pa | pb, "c4.Fig10.state_network_remodeling.png", 18, 10, outdir = outdir)

  forest <- edges |> filter(remodeled) |> group_by(comparison) |>
    slice_min(FDR, n = 25, with_ties = FALSE) |> ungroup() |>
    mutate(edge = paste(feature1, feature2, sep = " — "),
      edge = fct_reorder(edge, delta_r))
  pf <- if (!nrow(forest)) blank_plot("Differential co-abundance edge audit") else
    ggplot(forest, aes(fisher_delta, edge, color = remodeling_class, shape = known_PPI_backbone)) +
      geom_vline(xintercept = 0, color = "grey70") +
      geom_errorbarh(aes(xmin = ci_low, xmax = ci_high), height = .08) + geom_point(size = 2.1) +
      facet_wrap(~comparison, scales = "free_y") +
      scale_shape_manual(values = c(`TRUE` = 17, `FALSE` = 16),
        labels = c(`TRUE` = "Known-PPI annotation", `FALSE` = "No supplied PPI annotation")) +
      labs(title = "Differential residual-correlation edge audit",
        subtitle = paste0(ifelse(layer == "protein", "Protein", "Metabolite"),
          " co-abundance is not a physical interaction assay; known-PPI status is annotation only"),
        x = expression(Delta*" Fisher z-transformed residual correlation"), y = NULL, color = NULL, shape = NULL) +
      theme_5c(8) + theme(legend.position = "bottom")
  save_plot(pf, "c4.Fig11.state_network_edges.png", 15, 11, outdir = outdir)
}

run_c4_state_network <- function(dat, disease, sets, modules, features, covars,
  tvar, evar, bvar, layer, rawdir, outdir) {
  cache <- file.path(rawdir, "c4.state_network.rds")
  selected <- select_c4_network_features(disease, sets, modules, features)
  old <- read_stage_cache(cache)
  if (!is.null(old) && identical(old$code_version %||% NA_character_, C4_CODE_VERSION) &&
      identical(old$features %||% character(), selected)) {
    plot_c4_state_network(old, layer, outdir); return(old)
  }
  distal <- is.finite(dat[[tvar]]) & dat[[tvar]] >= 10
  near <- dat[[evar]] == 1 & is.finite(dat[[tvar]]) & dat[[tvar]] > 0 & dat[[tvar]] <= 2
  prevalent <- is.finite(dat[[bvar]]) & dat[[bvar]] < 0
  groups <- list(`Near-diagnosis incident vs distal` = near,
    `Baseline-prevalent vs distal` = prevalent)
  ppi <- read_c4_ppi_reference()
  ref <- c4_residual_matrix(dat[distal, , drop = FALSE], selected, covars)
  alts <- lapply(groups, function(ii) c4_residual_matrix(dat[ii %in% TRUE, , drop = FALSE], selected, covars))
  counts <- tibble(state = c("Distal reference", names(groups)),
    N = c(ref$n, vapply(alts, `[[`, integer(1), "n")))
  if (ref$n < C4_NETWORK_MIN_N || any(vapply(alts, `[[`, integer(1), "n") < C4_NETWORK_MIN_N)) {
    ans <- list(code_version = C4_CODE_VERSION, features = selected, edges = tibble(), hubs = tibble(),
      state_counts = counts, ppi_reference = ppi,
      status = tibble(status = "unavailable", detail = paste("At least one state has fewer than",
        C4_NETWORK_MIN_N, "complete-covariate participants")))
  } else {
    edges <- bind_rows(Map(function(alt, nm) c4_corr_edges(ref, alt, nm, ppi), alts, names(alts)))
    hubs <- make_c4_state_hubs(edges, sets, modules)
    ans <- list(code_version = C4_CODE_VERSION, features = selected, edges = edges, hubs = hubs,
      state_counts = counts, ppi_reference = ppi,
      status = tibble(status = "ok", detail = paste(length(selected),
        "features; residual co-abundance remodeling, not physical PPI rewiring")))
  }
  write_stage_cache(ans, cache); plot_c4_state_network(ans, layer, outdir); ans
}
