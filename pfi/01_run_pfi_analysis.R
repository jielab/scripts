# Publication Fairness Index figures and summary tables
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) stop("Usage: Rscript 01_run_pfi_analysis.R <fairness_results.csv> <plot_dir> [data_dir]")
infile <- args[1]; plot_dir <- args[2]; data_dir <- if (length(args) >= 3) args[3] else plot_dir
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE); dir.create(data_dir, recursive = TRUE, showWarnings = FALSE)
need <- c("data.table", "ggplot2")
for (p in need) if (!requireNamespace(p, quietly = TRUE)) install.packages(p, repos = "https://cloud.r-project.org")
suppressPackageStartupMessages({library(data.table); library(ggplot2)})
message("Reading: ", infile)
dat <- fread(infile)
if (!"fairness_category" %in% names(dat)) dat[, fairness_category := fifelse(fairness_z >= 2, "strongly_unfair_favored", fifelse(fairness_z >= 1, "modestly_unfair_favored", fifelse(fairness_z <= -2, "strongly_unfair_unfavored", fifelse(fairness_z <= -1, "modestly_unfair_unfavored", "fair"))))]
if (!"domain" %in% names(dat)) dat[, domain := "unknown"]
if (!"publication_type_proxy" %in% names(dat)) dat[, publication_type_proxy := "Unknown"]
cat_levels <- c("strongly_unfair_unfavored", "modestly_unfair_unfavored", "fair", "modestly_unfair_favored", "strongly_unfair_favored")
dat[, fairness_category := factor(fairness_category, levels = cat_levels)]
summary_by_year_domain <- dat[, .(n=.N, mean_actual_journal_score=mean(actual_journal_score_0_100, na.rm=TRUE), mean_expected_journal_score=mean(expected_journal_score_0_100, na.rm=TRUE), mean_predicted_paper_score=if ("predicted_paper_score" %in% names(dat)) mean(predicted_paper_score, na.rm=TRUE) else NA_real_, mean_fairness_index=mean(fairness_index, na.rm=TRUE), sd_fairness_index=sd(fairness_index, na.rm=TRUE), strong_favored_n=sum(fairness_category=="strongly_unfair_favored", na.rm=TRUE), strong_unfavored_n=sum(fairness_category=="strongly_unfair_unfavored", na.rm=TRUE)), by=.(year, domain)][order(year, domain)]
fwrite(summary_by_year_domain, file.path(data_dir, "fairness_summary_by_year_domain.csv"))
category_counts <- dat[, .N, by=.(fairness_category)][order(fairness_category)]
fwrite(category_counts, file.path(data_dir, "fairness_category_counts.csv"))
category_by_domain <- dat[, .N, by=.(domain, fairness_category)]; category_by_domain[, prop := N/sum(N), by=domain]
fwrite(category_by_domain, file.path(data_dir, "fairness_category_by_domain.csv"))
set.seed(1); plot_dat <- dat[is.finite(fairness_z)]; if (nrow(plot_dat) > 200000) plot_dat <- plot_dat[sample(.N, 200000)]
p1 <- ggplot(dat[is.finite(fairness_z)], aes(fairness_z)) + geom_histogram(bins=100) + theme_bw(base_size=14) + labs(x="Fairness z-score", y="Number of papers")
ggsave(file.path(plot_dir, "Fig3_prediction_distribution.png"), p1, width=8.2, height=4.8, dpi=300)
p2 <- ggplot(plot_dat[is.finite(actual_journal_score_0_100) & is.finite(expected_journal_score_0_100)], aes(expected_journal_score_0_100, actual_journal_score_0_100)) + geom_point(alpha=.18, size=.45) + geom_abline(slope=1, intercept=0, linetype=2) + theme_bw(base_size=14) + labs(x="Expected journal score from blinded content", y="Actual field-normalized journal score")
ggsave(file.path(plot_dir, "Fig4_actual_vs_expected_journal_score.png"), p2, width=8.2, height=6.2, dpi=300)
if ("predicted_paper_score" %in% names(dat)) {p3 <- ggplot(plot_dat[is.finite(predicted_paper_score) & is.finite(actual_journal_score_0_100)], aes(predicted_paper_score, actual_journal_score_0_100)) + geom_point(alpha=.18, size=.45) + geom_abline(slope=1, intercept=0, linetype=2) + theme_bw(base_size=14) + labs(x="Predicted paper score from blinded content", y="Actual field-normalized journal score"); ggsave(file.path(plot_dir, "actual_journal_vs_predicted_paper_score.png"), p3, width=8.2, height=6.2, dpi=300)}
p4 <- ggplot(category_by_domain, aes(domain, prop, fill=fairness_category)) + geom_col() + coord_flip() + theme_bw(base_size=13) + labs(x=NULL, y="Proportion of papers", fill="Fairness category")
ggsave(file.path(plot_dir, "Fig5_fairness_categories_by_domain.png"), p4, width=8.2, height=5.6, dpi=300)
p5 <- ggplot(summary_by_year_domain[is.finite(mean_fairness_index)], aes(year, mean_fairness_index, group=domain)) + geom_line() + geom_point(size=.7) + theme_bw(base_size=14) + labs(x="Publication year", y="Mean fairness index")
ggsave(file.path(plot_dir, "mean_fairness_by_year_domain.png"), p5, width=8.2, height=5, dpi=300)
val_file <- file.path(dirname(data_dir), "pubmedbert_validation_predictions.csv")
if (file.exists(val_file)) {
  val <- fread(val_file)
  if (all(c("expected_journal_score_0_100", "pred_expected_journal_score_0_100") %in% names(val))) {
    p6 <- ggplot(val, aes(pred_expected_journal_score_0_100, expected_journal_score_0_100)) + geom_point(alpha=.35, size=.7) + geom_abline(slope=1, intercept=0, linetype=2) + theme_bw(base_size=14) + labs(x="Validation predicted journal score", y="Validation actual journal score")
    ggsave(file.path(plot_dir, "Fig2_training_fit_journal.png"), p6, width=8.2, height=6, dpi=300)
  }
  if (all(c("paper_score", "pred_paper_score") %in% names(val))) {
    p7 <- ggplot(val, aes(pred_paper_score, paper_score)) + geom_point(alpha=.35, size=.7) + geom_abline(slope=1, intercept=0, linetype=2) + theme_bw(base_size=14) + labs(x="Validation predicted paper score", y="Validation paper label")
    ggsave(file.path(plot_dir, "Fig2_training_fit_paper.png"), p7, width=8.2, height=6, dpi=300)
  }
}
if ("famous_institution_proxy" %in% names(dat)) {
  inst <- dat[, .(n=.N, famous_institution_rate=mean(famous_institution_proxy==1, na.rm=TRUE)), by=.(fairness_category)]
  fwrite(inst, file.path(data_dir, "famous_institution_proxy_by_fairness_category.csv"))
  p8 <- ggplot(inst, aes(fairness_category, famous_institution_rate)) + geom_col() + coord_flip() + theme_bw(base_size=13) + labs(x=NULL, y="Rate of famous-institution proxy")
  ggsave(file.path(plot_dir, "Fig6_outlier_author_institution_audit.png"), p8, width=8.2, height=4.8, dpi=300)
}
setorder(dat, -fairness_z); fwrite(head(dat, 500), file.path(data_dir, "top_strongly_favored_candidates_for_manual_audit.csv"))
setorder(dat, fairness_z); fwrite(head(dat, 500), file.path(data_dir, "top_strongly_unfavored_candidates_for_manual_audit.csv"))
cat("Done. Figures written to", plot_dir, "and data tables to", data_dir, "\n")
