# MAHA manuscript engine, reorganized from jielab/scripts 971a61b44b84500ab46e2f5592478505a0eae50e.
maha_load_state("phe")
nhanes_validation <- nhanes_validate_maha(dat_mort_obs, covs.mort)
nhanes_construct <- nhanes_construct_profile(dat1_obs)
write_xlsx(nhanes_meat_sensitivity(dat_mort_obs, covs.mort), "Reviewer.meat_proxy_sensitivity.xlsx")

maha_run_step("fig3", {
step_header("Figure 3: NHANES all-cause mortality")

Y.mort.main <- "death"
fig3.all_map <- c(
	maha = "MAHA", dash = "DASH", mind = "MIND", medi = "MEDI",
	maha_bal = "MAHA-balanced", maha_strict = "MAHA-strict",
	maha_nodairy = "MAHA-no dairy", maha_noprotein = "MAHA-no protein"
)
fig3.discord_map <- c(
	dash = "DASH", medi = "MEDI", mind = "MIND",
	maha_bal = "MAHA-balanced", maha_strict = "MAHA-strict",
	maha_nodairy = "MAHA-no dairy", maha_noprotein = "MAHA-no protein"
)
fig3.shared_cols <- c(
	"MAHA" = "#F26D60", "DASH" = "#16B9C0", "MIND" = "#B97AF7", "MEDI" = "#7CAE00",
	"MAHA-balanced" = "#E76F51", "MAHA-strict" = "#C77DFF",
	"MAHA-no dairy" = "#2A9D8F", "MAHA-no protein" = "#577590"
)
fig3.legend_levels <- unname(fig3.all_map)
fig3.discord_levels <- unname(fig3.discord_map)

Fig3.all <- run_fig3_assoc_data(dat_mort_obs, Y.mort.main, covs.mort, names(fig3.all_map), fig3.all_map, family = "cox")
Fig3.disc <- run_fig3_discord_data(dat_mort_obs, Y.mort.main, covs.mort, names(fig3.discord_map), fig3.discord_map, ref_nm = "maha", ref_lab = "MAHA", family = "cox", min_n_pair = 20)

Fig3.ab <- list(
	A = plot_fig3_line_panel(
		Fig3.all$assoc, Y.mort.main, "Diet", fig3.legend_levels, fig3.shared_cols,
		"a. Dietary patterns", "Hazard Ratio (per SD)", xlim = c(0.80, 1.10),
		show_y = TRUE, show_legend = TRUE, max_shift = 0.58, y_expand = 0.08, show_hr_text = TRUE
	)
)

Fig3.cd <- make_mortality_bar_panels(dat_mort_obs, outcome = "death", t0 = 10)
legend_ab <- cowplot::get_legend(Fig3.ab$A + theme(legend.position = "bottom"))
Fig3.top <- cowplot::plot_grid(
	Fig3.ab$A + theme(legend.position = "none"),
	legend_ab,
	ncol = 1, rel_heights = c(1, 0.16)
)
Fig3.right <- cowplot::plot_grid(Fig3.cd$C, Fig3.cd$D, ncol = 1, rel_heights = c(1, 1), align = "v", axis = "lr")
Fig3 <- cowplot::plot_grid(
	Fig3.top, NULL, Fig3.right,
	nrow = 1, rel_widths = c(0.69, 0.08, 0.60),
	align = "h", axis = "tb"
)

print(Fig3.all$diagnostics)
print(Fig3.all$assoc)
print(Fig3.disc$diagnostics)
print(Fig3.disc$disc)
})

maha_run_step("figS1", {
step_header("Figure S3: NHANES diet score distributions and correlation")
maha_load_packages(tidyverse, patchwork, cowplot)
dat1 <- dat1_obs

figS1_score_map <- c(
	maha = "MAHA", dash = "DASH", mind = "MIND", medi = "MEDI",
	maha_bal = "MAHA-balanced", maha_strict = "MAHA-strict",
	maha_nodairy = "MAHA-no dairy", maha_noprotein = "MAHA-no protein"
)
figS1_cols <- c(
	"MAHA" = "#F26D60", "DASH" = "#16B9C0", "MIND" = "#B97AF7", "MEDI" = "#7CAE00",
	"MAHA-balanced" = "#E76F51", "MAHA-strict" = "#C77DFF",
	"MAHA-no dairy" = "#2A9D8F", "MAHA-no protein" = "#577590"
)
figS1_dist_axis_text_size <- 8.2
figS1_top_height <- 0.76 * 1.30
figS1_gap_height <- 0.10
figS1_bottom_height <- 1.24
figS1_plot_height <- 9.6 * (figS1_top_height + figS1_gap_height + figS1_bottom_height) / (0.76 + 1.24)

for (nm in names(figS1_score_map)) {
	sumv <- paste0("diet.", nm, ".sum"); s100v <- paste0("diet.", nm, ".s100")
	if (!s100v %in% names(dat1) && sumv %in% names(dat1)) dat1[[s100v]] <- score_0_100(dat1[[sumv]])
}

figS1_vars <- paste0("diet.", names(figS1_score_map), ".s100")
missing_figS1 <- setdiff(figS1_vars, names(dat1))
if (length(missing_figS1)) stop("Missing FigS1 score variables: ", paste(missing_figS1, collapse = ", "))
figS1_var_labs <- setNames(unname(figS1_score_map), figS1_vars)

plot_dat <- dat1 %>%
	dplyr::select(all_of(figS1_vars)) %>%
	pivot_longer(everything(), names_to = "var", values_to = "score") %>%
	mutate(score = as.numeric(score), Diet_score = factor(figS1_var_labs[var], levels = rev(unname(figS1_score_map)))) %>%
	filter(is.finite(score), !is.na(Diet_score))

plot_sum <- plot_dat %>%
	group_by(Diet_score) %>%
	summarise(
		N = n(), Mean = mean(score), SD = sd(score), Median = median(score),
		P25 = quantile(score, 0.25) %>% as.numeric(), P75 = quantile(score, 0.75) %>% as.numeric(),
		Min = min(score), Max = max(score), .groups = "drop"
	) %>%
	mutate(across(c(Mean, SD, Median, P25, P75, Min, Max), ~ round(.x, 1))) %>%
	as.data.frame()
density_dat <- plot_dat %>% mutate(Diet_score = as.character(Diet_score)) %>%
	group_by(Diet_score) %>% group_modify( ~ {
		dn <- density(.x$score, from = 0, to = 100, n = 256, na.rm = TRUE)
		tibble(score = dn$x, density = dn$y)
	}) %>% ungroup()

p_dist <- ggplot(plot_dat, aes(x = score, y = Diet_score, fill = Diet_score)) +
	geom_violin(width = 0.82, scale = "width", trim = TRUE, alpha = 0.72, color = NA) +
	geom_boxplot(width = 0.12, outlier.shape = NA, alpha = 0.92, color = "grey25") +
	scale_fill_manual(values = figS1_cols, guide = "none") +
	scale_x_continuous(limits = c(0, 100), breaks = seq(0, 100, 20), expand = expansion(mult = c(0.01, 0.01))) +
	labs(x = "Harmonized dietary score (0-100)", y = NULL, title = "b. Score distribution (NHANES)") +
	theme_classic(base_size = 10) +
	theme(
		legend.position = "none",
		plot.title = element_text(face = "bold", hjust = 0.5, size = 12),
		axis.text.y = element_text(face = "bold", size = figS1_dist_axis_text_size, color = "black"),
		axis.text.x = element_text(face = "bold", size = figS1_dist_axis_text_size),
		axis.title.x = element_text(face = "bold", size = figS1_dist_axis_text_size),
		panel.grid.major = element_blank(),
		panel.grid.minor = element_blank(),
		plot.margin = margin(5.5, 12, 5.5, 5.5)
	)

score_mat <- dat1 %>% dplyr::select(all_of(figS1_vars)) %>% mutate(across(everything(), as.numeric))
cor_mat <- cor(score_mat, use = "pairwise.complete.obs", method = "spearman")
dimnames(cor_mat) <- rep(list(unname(figS1_score_map)), 2)
cor_dat <- as.data.frame(as.table(cor_mat), stringsAsFactors = FALSE) %>%
	mutate(Var1 = factor(Var1, levels = rev(unname(figS1_score_map))), Var2 = factor(Var2, levels = unname(figS1_score_map)), lab = sprintf("%.2f", Freq))

p_cor <- ggplot(cor_dat, aes(Var2, Var1, fill = Freq)) +
	geom_tile(color = "white", linewidth = 0.35) +
	geom_text(aes(label = lab), size = 2.55, fontface = "bold") +
	scale_fill_gradient2(low = "#4575B4", mid = "white", high = "#D73027", midpoint = 0.5, limits = c(0, 1), breaks = c(0, 0.5, 1), oob = scales::squish, name = NULL) +
	coord_fixed() +
	labs(x = NULL, y = NULL, title = "d. Spearman correlation (NHANES)") +
	theme_minimal(base_size = 10) +
	theme(
		plot.title = element_text(face = "bold", hjust = 0.5, size = 12),
		axis.text.x = element_text(angle = 38, hjust = 1, vjust = 1, face = "bold", size = 8.2),
		axis.text.y = element_text(face = "bold", size = figS1_dist_axis_text_size, color = "black"),
		panel.grid = element_blank(),
		legend.position = "right",
		legend.key.height = unit(0.55, "cm"),
		plot.margin = margin(5.5, 5.5, 5.5, 5.5)
	)

FigS1_right_nhanes <- p_dist / plot_spacer() / p_cor +
	plot_layout(heights = c(figS1_top_height, figS1_gap_height, figS1_bottom_height))
save_plot(FigS1_right_nhanes, "FigS3.score_distribution_concordance.png", width = 6.9, height = figS1_plot_height, dpi = 500, bg = "white")

write_xlsx(
	list(
		score_distribution = plot_sum,
		density_curve = density_dat,
		spearman_matrix = as.data.frame(cor_mat, check.names = FALSE),
		pairwise_N = as.data.frame(crossprod(!is.na(score_mat))),
		metadata = data.frame(N = nrow(score_mat), sample = "dat1_obs: observed repeated-recall mean, age >=20, eligible dietary weights; descriptive main cycles",
			cycles = paste(main_cycles, collapse = "; "), method = "Unweighted Spearman, pairwise complete observations",
			shared_figure = "Shared FigS1 and NHANES construct figure read this exact same correlation matrix")
	),
	"FigS3.score_distribution_concordance.out.xlsx"
)
print(plot_sum)
print(round(cor_mat, 2))
})

maha_save_state("assoc", c("Fig3.all", "Fig3.disc", "Fig3.ab", "Fig3.cd", "fig3.shared_cols"))
maha_write_assoc(Fig3.all$assoc, "nhanes")
