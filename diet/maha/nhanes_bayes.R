# MAHA manuscript engine, reorganized from jielab/scripts 971a61b44b84500ab46e2f5592478505a0eae50e.
maha_load_state("phe")
maha_load_state("assoc")

step_header("Figure 4: NHANES mortality validation and precision checks")
maha_load_packages(tidyverse, patchwork, writexl, scales, cowplot)

ap_margin_main <- 1.05
ap_prior_main <- 1.10

ap_equiv_main <- ap_from_assoc(Fig3.all$assoc, "NHANES mortality")
ap_head2head <- purrr::map2_dfr(c("dash", "mind", "medi"), c("DASH", "MIND", "MEDI"), function(nm, lab) {
	ap_fit_head2head_nhanes(dat_mort_obs, Y = "death", trad_nm = nm, trad_lab = lab, covs.in = covs.mort, maha_nm = "maha", maha_lab = "MAHA")
})

ap_equiv_clean <- ap_clean_equiv(ap_equiv_main, "NHANES")
ap_head_clean <- ap_clean_head2head(ap_head2head, "NHANES")
ap_maha_levels <- c("MAHA", "MAHA-balanced", "MAHA-strict", "MAHA-no dairy", "MAHA-no protein")
ap_maha_cols <- fig3.shared_cols[ap_maha_levels]

p4a <- Fig3.ab$A +
	labs(title = "a. Dietary patterns", y = "\n\n\nAll-cause mortality") +
	guides(color = guide_legend(nrow = 2, byrow = TRUE), shape = guide_legend(nrow = 2, byrow = TRUE)) +
	theme(
		legend.position = "top",
		legend.box = "vertical",
		legend.margin = margin(0, 0, 2, 0),
		legend.box.margin = margin(0, 0, 2, 0),
		plot.title = element_text(margin = margin(b = 6)),
		axis.title.y = element_text(face = "bold", angle = 90, vjust = 1.8, margin = margin(r =  - 28)),
		plot.margin = margin(18, 18, 5.5, - 6)
	)
p4b <- Fig3.cd$C +
	labs(title = "b. DASH gradient") +
	theme(
		panel.grid.major = element_blank(),
		panel.grid.minor = element_blank(),
		axis.title.y = element_text(face = "bold", margin = margin(r =  - 5))
	)
p4c <- Fig3.cd$D +
	labs(title = "c. MAHA gradient") +
	theme(
		panel.grid.major = element_blank(),
		panel.grid.minor = element_blank(),
		axis.title.y = element_text(face = "bold", margin = margin(r =  - 4)),
		plot.margin = margin(5.5, 5.5, 5.5, 0)
	)
p4d <- ap_plot_maha_equiv(ap_equiv_clean, "d. MAHA equivalence")
p4e <- ap_plot_bf(ap_equiv_clean, "e. Bayes evidence for near-zero")
p4f <- ap_plot_power(ap_equiv_clean, "f. MAHA detectable effect")
p4g <- ap_plot_head2head(ap_head_clean, "g. Direct attenuation versus established scores")

Fig4_top <- p4a + patchwork::plot_spacer() + p4b + patchwork::plot_spacer() + p4c + plot_layout(widths = c(1.08, 0.07, 1.00, 0.07, 1.00))
Fig4_bottom <- (p4d + p4e) / patchwork::plot_spacer() / (p4f + p4g) + plot_layout(heights = c(1.00, 0.08, 1.00))
Fig4_nhanes <- Fig4_top / patchwork::plot_spacer() / Fig4_bottom + plot_layout(heights = c(0.86, 0.06, 1.44))

save_plot(Fig4_nhanes, "Fig4.mortality_validation.png", width = 17.4, height = 13.0, dpi = 320)
write_xlsx(
	list(
		Fig4A_assoc = Fig3.all$assoc,
		Fig4A_diag = Fig3.all$diagnostics,
		Fig4_disc = Fig3.disc$disc,
		Fig4_disc_diag = Fig3.disc$diagnostics,
		Fig4B_risk_DASH = Fig3.cd$risk10_DASH_within_MAHA,
		Fig4C_risk_MAHA = Fig3.cd$risk10_MAHA_within_DASH,
		Fig4D_equiv = ap_equiv_clean %>% filter(Diet %in% ap_maha_levels) %>% dplyr::select(Dataset, Outcome, Diet, N_total, N_event, HR_90CI, TOST_p, TOST_P_label, Equivalence_5pct, margin_hr, CI90_low_HR, CI90_high_HR),
		Fig4E_BF = ap_equiv_clean %>% filter(Diet %in% ap_maha_levels) %>% dplyr::select(Dataset, Outcome, Diet, N_total, N_event, BF01, BF01_label, BF_summary, prior_95_hr),
		Fig4F_power = ap_equiv_clean %>% filter(Diet %in% ap_maha_levels) %>% dplyr::select(Dataset, Outcome, Diet, N_total, N_event, MDE80_HR_label, Powered_5pct, MDE_HR_lower, MDE_HR_upper, margin_hr),
		Fig4G_head2head = ap_head_clean,
		Fig4_raw_equiv = ap_equiv_main,
		Fig4_raw_h2h = ap_head2head,
		curves_DASH_MAHA = Fig3.cd$curves_DASH_within_MAHA,
		curves_MAHA_DASH = Fig3.cd$curves_MAHA_within_DASH
	),
	"Fig4.mortality_validation.out.xlsx"
)

print(ap_equiv_clean %>% filter(Diet %in% ap_maha_levels))
print(ap_head_clean)

maha_write_bayes(ap_equiv_main, "nhanes")
