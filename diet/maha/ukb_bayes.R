# MAHA manuscript engine, reorganized from jielab/scripts 971a61b44b84500ab46e2f5592478505a0eae50e.
maha_load_state("phe")
maha_load_state("assoc")

ap_ukb_equiv_main <- ap_effect_from_assoc(
	assoc.ab %>% filter(as.character(Diet) %in% c("MAHA", "DASH", "MIND", "MEDI")),
	"UK Biobank main outcomes"
)

ap_ukb_equiv_alt <- if (exists("assoc.cd")) {
	ap_effect_from_assoc(assoc.cd, "UK Biobank alternative MAHA")
} else {
	tibble()
}

ap_ukb_head2head <- purrr::map_dfr(Y.inc, function(y) {
	purrr::map2_dfr(c("dash", "mind", "medi24"), c("DASH", "MIND", "MEDI"), function(nm, lab) {
		ap_fit_head2head_ukb(dat, y, nm, lab, covs, maha_nm = "maha", maha_lab = "MAHA")
	})
})

ap_margin_main <- 1.05
ap_prior_main <- 1.10

if (!exists("diet_levels.ab")) diet_levels.ab <- c("MAHA", "DASH", "MIND", "MEDI")
if (!exists("outcome_levels")) outcome_levels <- rev(unname(dx.lst[Y.inc]))
if (!exists("fig3_outcome_cols")) {
	fig3_outcome_cols <- c(
		"Coronary artery disease" = "#B97AF7",
		"Ischemic stroke" = "#F26D60",
		"Heart failure" = "#E68613",
		"Type 2 diabetes" = "#7CAE00",
		"Chronic kidney disease" = "#16B9C0",
		"All-cause mortality" = "#2C7FB8"
	)
}

ap_shape_guide <- guides(
	shape = guide_legend(
		nrow = 1,
		byrow = TRUE,
		keywidth = grid::unit(0.56, "cm"),
		override.aes = list(size = 3.4, stroke = 1.2)
	)
)

ap_leg_theme <- theme(
	legend.position = "bottom",
	legend.spacing.x = grid::unit(0.95, "cm"),
	legend.key.width = grid::unit(0.52, "cm"),
	legend.text = element_text(face = "bold", margin = margin(l = 4, r = 26))
)

ap_fig3_ukb <- ap_make_fig3_outputs(
	ap_ukb_equiv_main,
	ap_ukb_head2head,
	"UK Biobank",
	"Fig3.precision_attenuation",
	equiv_alt_tbl = if (exists("ap_ukb_equiv_alt")) ap_ukb_equiv_alt else NULL
)

maha_write_bayes(ap_ukb_equiv_main, "ukb")
