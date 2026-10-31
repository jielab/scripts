#!/usr/bin/env Rscript
# signature: prepare, protein_cox, metabolite_cox, metabolite_plot
# Launch with dose.sh; source this file to register definitions for reuse.
if (!exists("dose_step_code", mode="function")) {
    source(file.path(Sys.getenv("DOSE_SCRIPT_DIR", if (.Platform$OS.type == "windows")
        "D:/scripts/dose/f" else "/mnt/d/scripts/dose/f"), "0.common.R"), local=TRUE)
}

# ===========================================================================
# signature / prepare
# ===========================================================================
dose_register_step("signature", "prepare", expression({
s5_load_packages()
s5_prepare_lasso_input("protein", dir_project = dose_output_dir(), dir_ukb = dose_ukb_dir())
s5_prepare_lasso_input("metabolite", dir_project = dose_output_dir(), dir_ukb = dose_ukb_dir())
}))

# ===========================================================================
# signature / protein_cox
# ===========================================================================
dose_register_step("signature", "protein_cox", expression({
# Proteomic signature -> incident IS Cox association (continuous + quartiles),
# followed by DOSE -> proteomic signature -> IS mediation using three methods.

options(stringsAsFactors = FALSE)
DIR_PROJECT <- dose_output_dir()
FILE_FUNCTIONS <- file.path(dose_script_dir(), "0.common.R")
FILE_SCORE <- file.path(
    DIR_PROJECT, "ukb/06signature/protein/protein_signature_scores.csv"
)
if (!file.exists(FILE_FUNCTIONS)) stop("Missing shared functions: ", FILE_FUNCTIONS)

# The fitted-model prediction is the statistically coherent primary score.
# Two automatic sensitivity analyses are added below:
# (1) the exact score convention in Sugar_depression/s20_Lasso.py; and
# (2) an out-of-fold Cox analysis (not present in the reference repository).
SCORE_COLUMN <- "prot_score"

result <- s5_run_signature_analysis(
    score_file = FILE_SCORE,
    score_col = SCORE_COLUMN,
    signature_label = "DOSE-related proteomic signature",
    output_prefix = "protein_signature",
    dir_project = DIR_PROJECT,
    sims = 100L,
    seed = 12345L,
    run_mediation = TRUE
)
print(result$cox$continuous)
print(result$cox$categories)
print(dplyr::bind_rows(result$mediation))

reference_result <- s5_run_signature_analysis(
    score_file = FILE_SCORE,
    score_col = "prot_score_reference",
    signature_label = "DOSE-related proteomic signature (Sugar-depression score convention)",
    output_prefix = "protein_signature_reference_style",
    dir_project = DIR_PROJECT,
    sims = 100L,
    seed = 12345L,
    run_mediation = TRUE
)

oof_result <- s5_run_signature_analysis(
    score_file = FILE_SCORE,
    score_col = "prot_score_oof",
    signature_label = "DOSE-related proteomic signature (OOF sensitivity)",
    output_prefix = "protein_signature_oof",
    dir_project = DIR_PROJECT,
    run_mediation = FALSE
)
print(reference_result$cox$continuous)
print(oof_result$cox$continuous)
}))

# ===========================================================================
# signature / metabolite_cox
# ===========================================================================
dose_register_step("signature", "metabolite_cox", expression({
# Metabolomic signature -> incident IS Cox association (continuous + quartiles),
# followed by DOSE -> metabolomic signature -> IS mediation using three methods.

options(stringsAsFactors = FALSE)
DIR_PROJECT <- dose_output_dir()
FILE_FUNCTIONS <- file.path(dose_script_dir(), "0.common.R")
FILE_SCORE <- file.path(
    DIR_PROJECT, "ukb/06signature/metabolite/metabolite_signature_scores.csv"
)
if (!file.exists(FILE_FUNCTIONS)) stop("Missing shared functions: ", FILE_FUNCTIONS)

# Primary fitted-model prediction plus automatic reference-style and OOF
# sensitivity analyses, mirroring the protein workflow.
SCORE_COLUMN <- "met_score"

result <- s5_run_signature_analysis(
    score_file = FILE_SCORE,
    score_col = SCORE_COLUMN,
    signature_label = "DOSE-related metabolomic signature",
    output_prefix = "metabolite_signature",
    dir_project = DIR_PROJECT,
    sims = 100L,
    seed = 12345L,
    run_mediation = TRUE
)
print(result$cox$continuous)
print(result$cox$categories)
print(dplyr::bind_rows(result$mediation))

reference_result <- s5_run_signature_analysis(
    score_file = FILE_SCORE,
    score_col = "met_score_reference",
    signature_label = "DOSE-related metabolomic signature (Sugar-depression score convention)",
    output_prefix = "metabolite_signature_reference_style",
    dir_project = DIR_PROJECT,
    sims = 100L,
    seed = 12345L,
    run_mediation = TRUE
)

oof_result <- s5_run_signature_analysis(
    score_file = FILE_SCORE,
    score_col = "met_score_oof",
    signature_label = "DOSE-related metabolomic signature (OOF sensitivity)",
    output_prefix = "metabolite_signature_oof",
    dir_project = DIR_PROJECT,
    run_mediation = FALSE
)
print(reference_result$cox$continuous)
print(oof_result$cox$continuous)
}))

# ===========================================================================
# signature / metabolite_plot
# ===========================================================================
dose_register_step("signature", "metabolite_plot", expression({
# Read saved metabolomic results; no model fitting and no modification of inputs.
# install.packages(c("ggplot2", "ggrepel", "jsonlite", "patchwork"))

DIR_PROJECT <- dose_output_dir()
local_lib <- file.path(DIR_PROJECT, ".r-lib")
if (dir.exists(local_lib)) .libPaths(c(local_lib, .libPaths()))
packages <- c("ggplot2", "ggrepel", "jsonlite", "patchwork")
missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) stop("Install packages: ", paste(missing, collapse = ", "))
library(ggplot2)
library(patchwork)

DIR_IN <- file.path(DIR_PROJECT, "ukb/06signature/metabolite")
DIR_ASSOC <- file.path(DIR_IN, "association_mediation")
DIR_OUT <- file.path(DIR_IN, "figures_R")
dir.create(DIR_OUT, recursive = TRUE, showWarnings = FALSE)
MODEL <- "M2"
N_LABEL_EACH_SIGN <- 4L
read_result <- function(suffix) read.csv(file.path(
    DIR_ASSOC, paste0("metabolite_signature_", suffix, ".csv")), check.names = FALSE)
one <- function(d, column, value) {
    out <- d[d[[column]] == value, , drop = FALSE]
    if (nrow(out) != 1L) stop("Expected exactly one row: ", column, " = ", value)
    out
}
star <- function(p) ifelse(is.na(p), "", ifelse(p < .001, "***",
    ifelse(p < .01, "**", ifelse(p < .05, "*", ""))))
comma <- function(x) format(x, big.mark = ",", scientific = FALSE, trim = TRUE)
summary <- jsonlite::fromJSON(file.path(DIR_IN, "metabolite_LASSO_summary.json"))
coef <- read.csv(file.path(DIR_IN, "metabolite_LASSO_selected.csv"))
cont <- one(read_result("IS_continuous"), "model", MODEL)
quart <- read_result("IS_quartiles")
quart <- quart[quart$model == MODEL, ]
overall <- one(read_result("IS_overall_tests"), "model", MODEL)
med <- one(read_result("three_method_mediation"), "method", "cox_product")
stopifnot(nrow(coef) == summary$n_selected,
          setequal(quart$category, paste0("Q", 1:4)), nrow(quart) == 4L,
          cont$N == med$N, cont$Events == med$Events,
          all(quart$N == cont$N), all(quart$Events == cont$Events),
          all(is.finite(coef$coefficient_standardized_log1p)))

ink <- "#23343D"; purple <- "#827AB8"; cyan <- "#A7DCE6"
teal <- "#379AAF"; green <- "#78C9AD"
theme_set(theme_classic(base_size = 12, base_family = "sans"))
panel_theme <- theme(plot.title = element_text(face = "bold", size = 14, color = ink),
    plot.subtitle = element_text(size = 10, color = "#596A73"),
    plot.caption = element_text(size = 9, color = "#596A73", hjust = 0),
    plot.margin = margin(14, 16, 14, 16))
blank <- function(title, subtitle = NULL, caption = NULL) {
    ggplot() + coord_cartesian(xlim = c(0, 1), ylim = c(0, 1), expand = FALSE) +
        theme_void(base_family = "sans") + panel_theme +
        labs(title = title, subtitle = subtitle, caption = caption)
}
box <- function(x, y, label, fill, border, size = 4.2) {
    annotate("label", x = x, y = y, label = label, size = size,
             fill = fill, color = border, fontface = "bold", lineheight = 1.25,
             label.padding = grid::unit(.65, "lines"),
             label.r = grid::unit(.35, "lines"))
}
arrow <- function(x, y, xend, yend) annotate("segment", x = x, y = y,
    xend = xend, yend = yend, linewidth = .65, color = ink,
    arrow = grid::arrow(length = grid::unit(2.6, "mm"), type = "closed"))

# A: Reference Figure 1C. Use actual LASSO (not elastic net) parameters.
pA <- blank("A   Construction of the metabolomic signature",
    caption = sprintf("Training R-squared = %.3f; out-of-fold correlation with DOSE = %.3f",
                      summary$full_R2, summary$oof_correlation_with_dose_z)) +
    box(.5, .88, paste0("UK Biobank\n", comma(summary$n_participants), " participants"),
        "#E8F3F7", "#477784") + arrow(.5, .76, .5, .69) +
    box(.5, .58, sprintf("DOSE-associated metabolites: %s\nBH-FDR < 0.05; retained after QC: %s",
        summary$n_fdr_features, summary$n_qc_features), "#FBE7D1", "#9B6A37", 4) +
    annotate("text", x = .5, y = .35, size = 4, lineheight = 1.5, color = ink,
        label = sprintf("Median imputation, log1p, standardization\nLASSO; %s-fold cross-validation\nalpha = %g; minimum mean CV error",
                        summary$folds, summary$alpha_min)) +
    arrow(.5, .23, .5, .17) +
    box(.5, .075, sprintf("DOSE-related metabolomic signature\n%s selected metabolites",
                         summary$n_selected), "#E3ECD9", "#617F45", 4)

# B: Reference Figure 1D. X ordering is alphabetical, not importance ranking.
coef <- coef[order(coef$metabolite), ]
coef$index <- seq_len(nrow(coef))
coef$beta <- coef$coefficient_standardized_log1p
coef$direction <- ifelse(coef$beta > 0, "Positive", "Negative")
coef$label <- sub("^met_", "", coef$metabolite)
labels <- do.call(rbind, lapply(c("Positive", "Negative"), function(s) {
    d <- coef[coef$direction == s, ]
    head(d[order(abs(d$beta), decreasing = TRUE), ], N_LABEL_EACH_SIGN)
}))
pB <- ggplot(coef, aes(index, beta)) +
    geom_hline(yintercept = 0, color = "#DCB825", linetype = "dashed", linewidth = .6) +
    geom_point(aes(color = direction), size = 2.7, alpha = .9) +
    ggrepel::geom_text_repel(data = labels, aes(label = label, color = direction),
        seed = 42, size = 3.5, box.padding = .65, point.padding = .35,
        min.segment.length = 0, max.overlaps = Inf, show.legend = FALSE,
        max.iter = 20000) +
    scale_color_manual(values = c(Positive = purple, Negative = cyan)) +
    scale_y_continuous(expand = expansion(mult = c(.18, .18))) +
    scale_x_continuous(expand = expansion(mult = .08)) +
    labs(title = "B   Metabolite coefficients in the LASSO model",
         subtitle = sprintf("%s selected metabolites; largest %s coefficients of each sign labelled",
                            nrow(coef), N_LABEL_EACH_SIGN),
         x = "Metabolites (alphabetical order)", y = "Coefficient (standardized log1p metabolite)",
         color = NULL, caption = "Coefficients predict standardized DOSE; they are not stroke effect estimates.") +
    panel_theme + theme(legend.position = "top", legend.justification = "left")

# C: Reference Figure 1E. Primary fitted-model signature, M2 only.
forest <- rbind(data.frame(group = quart$category, quart[c("HR", "LCI", "HCI", "P")]),
                data.frame(group = "Per 1 SD", cont[c("HR", "LCI", "HCI", "P")]))
forest$group <- factor(forest$group, levels = c("Per 1 SD", "Q4", "Q3", "Q2", "Q1"))
forest$y <- as.numeric(forest$group)
forest$stars <- star(forest$P)
forest$estimate <- ifelse(forest$group == "Q1", "1.000 (reference)",
    sprintf("%.3f (%.3f - %.3f)", forest$HR, forest$LCI, forest$HCI))
pC <- ggplot(forest, aes(HR, y)) +
    geom_vline(xintercept = 1, linetype = "dashed", color = "#CE8778", linewidth = .7) +
    geom_segment(aes(x = LCI, xend = HCI, yend = y), color = green, linewidth = 1) +
    geom_point(color = teal, size = 3.5) +
    geom_text(aes(label = stars), nudge_y = .24, size = 4.5) +
    geom_text(aes(x = 1.055, label = estimate), hjust = 0, size = 3.5) +
    annotate("text", x = 1.055, y = 5.63, label = "HR (95% CI)", hjust = 0,
             size = 3.6, fontface = "bold") +
    scale_y_continuous(breaks = 1:5, labels = levels(forest$group), limits = c(.5, 5.9)) +
    scale_x_continuous(breaks = seq(.6, 1, .1), limits = c(.60, 1.48)) +
    labs(title = "C   Metabolomic signature and ischemic stroke",
         subtitle = sprintf("%s; N = %s; events = %s", MODEL, comma(cont$N), comma(cont$Events)),
         x = "Hazard ratio", y = NULL,
         caption = sprintf("P for trend = %.2e. Q1 is the reference.\n* P < 0.05; ** P < 0.01; *** P < 0.001.", overall$P_trend)) +
    panel_theme + theme(axis.line.y = element_blank(), axis.ticks.y = element_blank())

# D: Reference Figure 2A. Cox product mediation, with adjusted b and direct c'.
# Never substitute b_reported_estimate: that model does not adjust for DOSE.
pD <- blank("D   Mediation through the metabolomic signature",
    subtitle = "M2; Cox product method; exposure contrast: DOSE +1 SD",
    caption = "Path b adjusts for DOSE; mediator is standardized. Proportion CI: approximate delta method.\n*** P < 0.001. Observational mediation estimates do not establish causality.") +
    box(.50, .87, "DOSE-related\nmetabolomic signature", "#E3ECD9", "#617F45") +
    box(.13, .18, "DOSE", "#FBE7D1", "#B87637") +
    box(.85, .18, "Ischemic stroke", "#DCECF6", "#548AA4", 4) +
    arrow(.16, .29, .32, .74) + arrow(.68, .74, .84, .29) +
    arrow(.26, .18, .69, .18) +
    annotate("text", x = .17, y = .53, angle = 61, size = 4,
        label = sprintf("a: beta = %.3f%s", med$a_estimate, star(med$a_P))) +
    annotate("text", x = .84, y = .53, angle = -61, size = 4,
        label = sprintf("b: HR = %.3f%s", exp(med$b_adjusted_estimate), star(med$b_adjusted_P))) +
    annotate("text", x = .5, y = .51, size = 4.5, lineheight = 1.6,
        label = sprintf("Mediation proportion\n%.2f%%%s\n95%% CI %.2f%% - %.2f%%",
                        100 * med$Proportion_mediated, star(med$P_proportion),
                        100 * med$Proportion_LCI, 100 * med$Proportion_HCI)) +
    annotate("text", x = .48, y = .105, size = 3.6,
        label = sprintf("Direct c': HR = %.3f%s", exp(med$c_direct_model_estimate), star(med$c_direct_model_P))) +
    annotate("text", x = .5, y = .015, size = 3.8,
        label = sprintf("Total c: HR = %.3f%s", exp(med$c_total_model_estimate), star(med$c_total_model_P)))

save_plot <- function(p, name, width = 8.8, height = 6.6) {
    bounded_ggsave(file.path(DIR_OUT, paste0(name, ".png")), p, width = width, height = height,
           dpi = 300, bg = "white")
    bounded_ggsave(file.path(DIR_OUT, paste0(name, ".pdf")), p, width = width, height = height,
           device = grDevices::cairo_pdf, bg = "white")
}
plots <- list(A_construction = pA, B_coefficients = pB,
              C_cox_forest_M2 = pC, D_mediation_cox_M2 = pD)
for (name in names(plots)) save_plot(plots[[name]], paste0("metabolite_", name))
combined <- (pA | pB) / (pC | pD)
save_plot(combined, "metabolite_four_panels", 17.6, 13.2)
writeLines(c(
    "A/B/C/D correspond to reference Figures 1C/1D/1E/2A.",
    "All numeric annotations come from existing primary metabolite result files.",
    "A: actual LASSO parameters; B: selected coefficients, not causal importance.",
    "C: M2 signature Cox models, without additional DOSE adjustment.",
    "D: M2 Cox product mediation; b adjusts for DOSE; c' is direct, c is total.",
    "M2: age, sex, ethnicity, education, TDI, energy, smoking, physical activity, BMI,",
    "family history of stroke, hypertension and type 2 diabetes; no separate alcohol adjustment.",
    "PH_P is NA in source outputs; proportional hazards remains unverified.",
    "The alpha optimum is at the tested lower boundary; OOF screening was not nested.",
    "Reference: https://link.springer.com/article/10.1186/s12937-026-01352-4"
), file.path(DIR_OUT, "figure_notes.txt"))
message("Saved four individual PNG/PDF figures and a combined figure: ", DIR_OUT)
}))

if (sys.nframe() == 0L) dose_run_module("signature")
