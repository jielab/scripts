# Rebuild presentation artifacts from the persisted dual-stage CKM object.

make_ckm_summary <- function(clock_source) {
  original_method <- stage3_method
  on.exit(assign("stage3_method", original_method, envir = .GlobalEnv), add = TRUE)

  clocks <- list()
  for (method in c("SCORE2", "PREVENT")) {
    progress_log("summarize_ckm: fitting stage clock for ", method)
    assign("stage3_method", method, envir = .GlobalEnv)
    d <- as.data.frame(clock_source)
    d$ckm.stage <- as.integer(d[[paste0("ckm.stage.", method)]])
    d$ckm_stage_certainty <- as.character(d[[paste0("ckm_stage_certainty_", method)]])
    d$ckm_stage3_subtype <- as.character(d[[paste0("ckm_stage3_subtype_", method)]])
    clocks[[method]] <- make_stage_clock(prepare_target_data(d, "ckm4"), write_outputs = FALSE)
    save_rds_safe(clocks[[method]], file.path(global_rawdir, paste0("stage_clock.", method, ".rds")), compress = FALSE)
  }

  method_plot <- function(p, method, panel) {
    p + labs(title = paste0(method, " — ", panel))
  }
  combined <-
    (method_plot(clocks$SCORE2$plots$profile, "SCORE2", "a. Risk-calibrated stage profile") |
       method_plot(clocks$PREVENT$plots$profile, "PREVENT", "a. Risk-calibrated stage profile")) /
    (method_plot(clocks$SCORE2$plots$rmst, "SCORE2", "b. Stage-4-free RMST profile") |
       method_plot(clocks$PREVENT$plots$rmst, "PREVENT", "b. Stage-4-free RMST profile")) /
    (method_plot(clocks$SCORE2$plots$incidence, "SCORE2", "c. Observed cumulative incidence") |
       method_plot(clocks$PREVENT$plots$incidence, "PREVENT", "c. Observed cumulative incidence")) +
    plot_layout(guides = "collect") +
    plot_annotation(title = "CKM stage-clock comparison: SCORE2 versus PREVENT") &
    theme(legend.position = "bottom")

  save_plot(combined, "Fig1.stage_clock.png", width = 18, height = 21, dir = global_outroot)

  clock_book <- list(
    SCORE2_metadata = clocks$SCORE2$metadata,
    PREVENT_metadata = clocks$PREVENT$metadata,
    SCORE2_risk_age = clocks$SCORE2$risk_age,
    PREVENT_risk_age = clocks$PREVENT$risk_age,
    SCORE2_risk_age_intervals = clocks$SCORE2$risk_age_intervals,
    PREVENT_risk_age_intervals = clocks$PREVENT$risk_age_intervals,
    SCORE2_risk_profile = clocks$SCORE2$primary_risk,
    PREVENT_risk_profile = clocks$PREVENT$primary_risk,
    SCORE2_RMST = clocks$SCORE2$primary_rmst,
    PREVENT_RMST = clocks$PREVENT$primary_rmst,
    SCORE2_observed_risk = clocks$SCORE2$observed_risk,
    PREVENT_observed_risk = clocks$PREVENT$observed_risk,
    SCORE2_stage4_death = clocks$SCORE2$stage4_death,
    PREVENT_stage4_death = clocks$PREVENT$stage4_death
  )
  write_book(clock_book, "Fig1.stage_clock.out.xlsx", dir = global_outroot)

  summarize_ckm_audit()
  invisible(clocks)
}

summarize_ckm_audit <- function() {
  paths <- list(
    stage = file.path(auditdir, "ckm_stage_summary.tsv"),
    stage2 = file.path(auditdir, "ckm_stage2_criterion_overlap.tsv"),
    stage3 = file.path(auditdir, "ckm_stage3_ascertainment.tsv"),
    stage3_sens = file.path(auditdir, "ckm_stage3_definition_sensitivity.tsv"),
    stage4 = file.path(auditdir, "ckm_stage4_sensitivity.tsv")
  )
  missing <- names(paths)[!file.exists(unlist(paths))]
  if (length(missing)) {
    progress_log("summarize_ckm: FigS1 skipped; missing audit tables: ", paste(missing, collapse = ", "))
    return(invisible(NULL))
  }

  stage <- fread(paths$stage)
  stage[, stage_label := ifelse(
    is.na(ckm_stage), "Unclassified / insufficient information",
    unname(stage_full_labels[as.character(as.integer(ckm_stage))])
  )]
  stage[, stage_label := factor(stage_label, levels = c(
    unname(stage_full_labels[as.character(0:4)]), "Unclassified / insufficient information"
  ))]
  stage[, proportion := N / sum(N)]

  overlap <- fread(paths$stage2)
  component_cols <- intersect(c("T2D", "hypertension", "hyperTG", "metabolic_syndrome", "CKD_moderate_high"), names(overlap))
  component_labels <- c(
    T2D = "Type 2 diabetes", hypertension = "Hypertension",
    hyperTG = "Hypertriglyceridemia", metabolic_syndrome = "Metabolic syndrome",
    CKD_moderate_high = "Moderate/high-risk CKD"
  )
  stage2 <- rbindlist(lapply(component_cols, function(v) {
    data.table(component = unname(component_labels[v]), positive_N = sum(overlap$N[overlap[[v]] %in% TRUE], na.rm = TRUE))
  }))
  stage2[, proportion := positive_N / sum(stage$N)]

  stage3 <- fread(paths$stage3)
  stage3[, assessed_proportion := assessed_N / sum(stage$N)]
  stage3[, coverage_label := ifelse(
    assessed_N == 0, "Not available",
    paste0(scales::percent(assessed_proportion, accuracy = .1), " assessed; ", scales::comma(positive_N), " positive")
  )]
  stage4 <- fread(paths$stage4)
  stage4_plot <- melt(
    stage4, id.vars = c("definition", "is_primary"),
    measure.vars = c("any_prebaseline_clinical_cvd", "stage4_with_required_ckm_risk"),
    variable.name = "estimand", value.name = "N"
  )

  p1 <- ggplot(stage, aes(stage_label, proportion, fill = stage_label)) +
    geom_col(width = .72, color = "white") +
    geom_text(aes(label = paste0(scales::comma(N), "\n", scales::percent(proportion, accuracy = .1))), vjust = -.15, size = 3) +
    scale_fill_manual(values = c(unname(stage_palette), "grey70"), guide = "none", drop = FALSE) +
    scale_y_continuous(labels = percent, expand = expansion(mult = c(0, .16))) +
    labs(title = "a. Baseline CKM-stage distribution", x = NULL, y = "Proportion") + theme_ckm(10)
  p2 <- ggplot(stage2, aes(proportion, fct_reorder(component, proportion))) +
    geom_col(fill = unname(stage_palette["2"])) +
    geom_text(aes(label = scales::percent(proportion, accuracy = .1)), hjust = -.1) +
    scale_x_continuous(labels = percent, expand = expansion(mult = c(0, .18))) +
    labs(title = "b. Stage-2 criteria", x = "Population prevalence", y = NULL) + theme_ckm(10)
  p3 <- ggplot(stage3, aes(assessed_proportion, fct_reorder(component, assessed_proportion))) +
    geom_col(fill = unname(stage_palette["3"])) + geom_text(aes(label = coverage_label), hjust = -.05, size = 2.7) +
    scale_x_continuous(labels = percent, limits = c(0, 1.12)) +
    labs(title = "c. Stage-3 ascertainment", x = "Proportion assessed", y = NULL) + theme_ckm(10)
  p4 <- ggplot(stage4_plot, aes(estimand, N, fill = definition, alpha = is_primary)) +
    geom_col(position = "dodge") + scale_alpha_manual(values = c(`TRUE` = 1, `FALSE` = .55), guide = "none") +
    labs(title = "d. Stage-4 source sensitivity", x = NULL, y = "Participants", fill = NULL) + theme_ckm(10)

  save_plot((p1 | p2) / (p3 | p4), "FigS1.stage_definition_audit.png", 15, 10, dir = global_outroot)
  write_book(list(
    stage_distribution = stage,
    stage2_components = stage2,
    stage3_ascertainment = stage3,
    stage3_sensitivity = fread(paths$stage3_sens),
    stage4_sensitivity = stage4
  ), "FigS1.stage_definition_audit.out.xlsx", dir = global_outroot)
  invisible(TRUE)
}
