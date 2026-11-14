# MAHA manuscript engine, reorganized from jielab/scripts 971a61b44b84500ab46e2f5592478505a0eae50e.
maha_load_state("phe")
# 🚩 7) FigS1 anti-random validation
step_header("CHNS FigS1: anti-random validation")

# The negative control deliberately keeps the SAME components, SAME participants,
# SAME equal weighting and SAME score scale; only component directions are changed.
# This directly addresses the criticism that a newly defined diet score could be arbitrary.
dval <- dat %>% drop_na(all_of(component_cols))
X <- as.matrix(dval[, component_cols, drop = FALSE])
k <- ncol(X)
if(nrow(dval) >= 300 && k >= 5) {
	grid <- expand.grid(rep(list(c( - 1, 1)), k), KEEP.OUT.ATTRS = FALSE)
	obs_i <- which(rowSums(grid == 1) == k); nul_i <- setdiff(seq_len(nrow(grid)), obs_i)
	if(null_max > 0 && length(nul_i) > null_max) nul_i <- nul_i[unique(round(seq(1, length(nul_i), length.out = null_max)))]
	sel <- c(obs_i, nul_i); grid <- grid[sel, , drop = FALSE]

	# Convergent anchor from established/China-established pattern scores.
	az <- sapply(anchor_cols, function(v) zstd(dval[[v]]))
	if(is.null(dim(az))) az <- matrix(az, ncol = 1)
	anchor <- rowMeans(az, na.rm = TRUE); anchor[!is.finite(anchor)] <- NA_real_

	# Completely independent questionnaire preference criterion.
	pref_anchor <- zstd(dval$pref_health_index)

	# Independent observed-behavior criterion: healthier means less processed/deep-fried/
	# sugary-drink exposure and a higher plant-vs-animal fat-source ratio.
	bz <- cbind( - zstd(dval$processed_food_g), - zstd(dval$deepfried_g), - zstd(dval$sugary_drinks_week), zstd(dval$fat_proxy_ratio))
	behavior_anchor <- rowMeans(bz, na.rm = TRUE); behavior_anchor[!is.finite(behavior_anchor)] <- NA_real_

	pat <- bind_rows(par_lapply(seq_len(nrow(grid)), function(i){
		signs <- as.numeric(grid[i, ])
		raw <- rowMeans(sweep(X, 2, signs, "*"))
		m <- fit_primary_cox(dval, raw)
		z_mort <- if(nrow(m) && is.finite(m$beta[1]) && is.finite(m$se[1]) && m$se[1] > 0) - m$beta[1] / m$se[1] else NA_real_
		tibble(
			pattern = paste(ifelse(signs > 0, "+", "-"), collapse = ""),
			is_prespecified = all(signs == 1),
			rho_established = safe_cor(raw, anchor),
			rho_preference = safe_cor(raw, pref_anchor),
			rho_behavior = safe_cor(raw, behavior_anchor),
			protective_Z = z_mort,
			mortality_N = if(nrow(m))m$N[1]else NA_real_,
			mortality_events = if(nrow(m))m$Events[1]else NA_real_,
			mortality_HR = if(nrow(m))m$HR[1]else NA_real_
		)
	}))
	obs <- pat %>% filter(is_prespecified); nul <- pat %>% filter(!is_prespecified)
	emp_high <- function(x, obsx) (1 + sum(x >= obsx, na.rm = TRUE)) / (1 + sum(is.finite(x)))
	pct <- function(x, obsx) mean(x <= obsx, na.rm = TRUE)
	sumval <- tibble(
		Cohort = "CHNS", Mode = analysis_mode, N = nrow(dval), Components = k, Null_scores = nrow(nul),
		Established_anchor_rho = obs$rho_established,
		Established_anchor_percentile = pct(nul$rho_established, obs$rho_established),
		Established_anchor_empirical_P = emp_high(nul$rho_established, obs$rho_established),
		Preference_rho = obs$rho_preference,
		Preference_percentile = pct(nul$rho_preference, obs$rho_preference),
		Preference_empirical_P = emp_high(nul$rho_preference, obs$rho_preference),
		Behavior_rho = obs$rho_behavior,
		Behavior_percentile = pct(nul$rho_behavior, obs$rho_behavior),
		Behavior_empirical_P = emp_high(nul$rho_behavior, obs$rho_behavior),
		Mortality_protective_Z = obs$protective_Z,
		Mortality_percentile = if(is.finite(obs$protective_Z))pct(nul$protective_Z, obs$protective_Z)else NA_real_,
		Mortality_empirical_P = if(is.finite(obs$protective_Z))emp_high(nul$protective_Z, obs$protective_Z)else NA_real_
	)

	primary_raw <- rowMeans(X)
	loo <- map_dfr(seq_along(component_keys), function(j){
		raw <- rowMeans(X[, - j, drop = FALSE])
		mm <- fit_primary_cox(dval, raw)
		tibble(
			omitted = component_labels[[component_keys[j]]],
			rho_primary = safe_cor(raw, primary_raw),
			rho_established = safe_cor(raw, anchor),
			rho_preference = safe_cor(raw, pref_anchor),
			mortality_HR = if(nrow(mm))mm$HR[1]else NA_real_,
			mortality_CI_low = if(nrow(mm))mm$CI_low[1]else NA_real_,
			mortality_CI_high = if(nrow(mm))mm$CI_high[1]else NA_real_
		)
	})

	p1 <- ggplot(nul, aes(rho_established)) + geom_histogram(bins = 35) + geom_vline(xintercept = obs$rho_established, linewidth = .9) +
		labs(title = "a. Established-score convergence", x = "Spearman rho with DASH/MIND/MEDI/CHDI anchor", y = "Matched arbitrary scores") + theme_classic(11)
	p2 <- ggplot(nul, aes(rho_preference)) + geom_histogram(bins = 35) + geom_vline(xintercept = obs$rho_preference, linewidth = .9) +
		labs(title = "b. Independent food-preference criterion", x = "Spearman rho with healthy preference index", y = "Matched arbitrary scores") + theme_classic(11)
	if(is.finite(obs$protective_Z) && sum(is.finite(nul$protective_Z)) >= 20) {
		p3 <- ggplot(nul, aes(protective_Z)) + geom_histogram(bins = 35) + geom_vline(xintercept = obs$protective_Z, linewidth = .9) +
			labs(title = "c. Mortality criterion", x = "Protective Z (-beta/SE)", y = "Matched arbitrary scores") + theme_classic(11)
	} else {
		p3 <- ggplot(nul, aes(rho_behavior)) + geom_histogram(bins = 35) + geom_vline(xintercept = obs$rho_behavior, linewidth = .9) +
			labs(title = "c. Independent observed-behavior criterion", x = "Spearman rho with healthier behavior anchor", y = "Matched arbitrary scores") + theme_classic(11)
	}
	lp <- loo %>% mutate(omitted = factor(omitted, levels = rev(omitted)))
	p4 <- ggplot(lp, aes(rho_primary, omitted)) + geom_segment(aes(x = 0, xend = rho_primary, yend = omitted)) + geom_point() +
		coord_cartesian(xlim = c(0, 1)) + labs(title = "d. Leave-one-domain-out rank stability", x = "Spearman rho with primary score", y = NULL) + theme_classic(11)
	FigS1 <- (p1 | p2) / (p3 | p4) + plot_annotation(title = paste0("CHNS: matched negative-control validation of ", active_label))
	save_plot(FigS1, "FigS1.validate.png", 14, 10.5, 450)
	write_xlsx(list(summary = sumval, direction_null = pat, leave_one_out = loo), "FigS1.validate.out.xlsx")
} else {
	warning("Insufficient complete cases for matched anti-random validation.")
}


# 🚩 8) FigS2 construct profile
step_header("CHNS FigS2: construct profile")
dcp <- dat %>% drop_na(all_of(component_cols), all_of(active_maha))
if(nrow(dcp) >= 200) {
	XX <- as.data.frame(dcp[, component_cols]); names(XX) <- component_keys
	XX$score <- dcp[[active_maha]]; XX$decile <- ntile(XX$score, 10)
	profile <- XX %>% group_by(decile) %>% summarise(across(all_of(component_keys), mean), .groups = "drop") %>%
		pivot_longer(all_of(component_keys), names_to = "component", values_to = "mean_score") %>%
		mutate(Component = unname(component_labels[component]))
	metrics <- map_dfr(component_keys, function(kk){
		rest <- rowMeans(as.matrix(XX[, setdiff(component_keys, kk), drop = FALSE]))
		tibble(Component = component_labels[[kk]], item_rest_rho = safe_cor(XX[[kk]], rest),
					 D10_minus_D1 = mean(XX[[kk]][XX$decile == 10]) - mean(XX[[kk]][XX$decile == 1]))
	})
	conv <- tibble(Comparator = comparison_labels,
				Spearman_rho = vapply(comparison_cols, function(v)safe_cor(dcp[[active_maha]], dcp[[v]]), numeric(1)))
	p1 <- ggplot(profile, aes(decile, mean_score, color = Component, group = Component)) + geom_line() + geom_point() +
		labs(title = "a. Component profile", x = paste(active_label, "decile"), y = "Mean component score", color = NULL) + theme_classic(10)
	mp <- metrics %>% mutate(Component = factor(Component, levels = rev(Component)))
	p2 <- ggplot(mp, aes(D10_minus_D1, Component)) + geom_vline(xintercept = 0, linetype = 2) + geom_segment(aes(x = 0, xend = D10_minus_D1, yend = Component)) + geom_point() +
		labs(title = "b. Decile separation", x = "Decile 10 - decile 1", y = NULL) + theme_classic(11)
	p3 <- ggplot(mp, aes(item_rest_rho, Component)) + geom_vline(xintercept = 0, linetype = 2) + geom_point() +
		labs(title = "c. Corrected item-rest rho", subtitle = "Descriptive only; score is formative", x = "Spearman rho", y = NULL) + theme_classic(11)
	p4 <- ggplot(conv, aes(Spearman_rho, reorder(Comparator, Spearman_rho))) + geom_segment(aes(x = 0, xend = Spearman_rho, yend = reorder(Comparator, Spearman_rho))) + geom_point() +
		coord_cartesian(xlim = c( - 1, 1)) + labs(title = "d. Convergent validity", x = "Spearman rho", y = NULL) + theme_classic(11)
	save_plot((p1 | p2) / (p3 | p4) + plot_annotation(title = paste0("CHNS: ", active_label, " construct profile")),
				"FigS2.construct_profile.png", 14, 10.5, 450)
	write_xlsx(list(decile_profile = profile, component_metrics = metrics, convergent_validity = conv), "FigS2.construct_profile.out.xlsx")
}


# 🚩 9) Figprimary independent food-preference validation
step_header("CHNS Figprimary: independent food-preference validation")
pref_outcomes <- c("pref_fruit", "pref_veg", "pref_fast", "pref_salty", "pref_ssb", "pref_health_index")
pref_lab <- c(pref_fruit = "Fruit preference", pref_veg = "Vegetable preference", pref_fast = "Fast-food preference",
			pref_salty = "Salty-snack preference", pref_ssb = "Sugary-drink preference", pref_health_index = "Healthy preference index")
pref_fit <- map_dfr(pref_outcomes, function(y){
	fit_linear_pref(dat, y, dat[[active_maha]]) %>% mutate(Outcome = pref_lab[[y]], Variable = y)
})
pref_fit$Expected <- ifelse(pref_fit$Variable %in% c("pref_fruit", "pref_veg", "pref_health_index"), "Positive", "Negative")
ddp <- dat %>% filter(is.finite(.data[[active_maha]])) %>% mutate(score_decile = ntile(.data[[active_maha]], 10))
pref_profile <- ddp %>% group_by(score_decile) %>% summarise(across(all_of(pref_outcomes), ~ mean(.x, na.rm = TRUE)), .groups = "drop") %>%
	pivot_longer( - score_decile, names_to = "Variable", values_to = "Mean") %>% mutate(Outcome = unname(pref_lab[Variable]))
pref_corr <- tibble(
	Criterion = c("Healthy preference index", "Fruit preference", "Vegetable preference", "Fast-food preference", "Salty-snack preference", "Sugary-drink preference"),
	Spearman_rho = c(safe_cor(dat[[active_maha]], dat$pref_health_index), safe_cor(dat[[active_maha]], dat$pref_fruit),
				safe_cor(dat[[active_maha]], dat$pref_veg), safe_cor(dat[[active_maha]], dat$pref_fast),
				safe_cor(dat[[active_maha]], dat$pref_salty), safe_cor(dat[[active_maha]], dat$pref_ssb))
)
empty_panel <- function(title, message = "No analyzable score-preference pairs") {
	ggplot() + annotate("text", x = 0, y = 0, label = message, size = 4) +
		xlim( - 1, 1) + ylim( - 1, 1) + labs(title = title) + theme_void()
}
if(any(is.finite(pref_fit$beta))) {
	p3a <- pref_fit %>% filter(is.finite(beta), is.finite(se)) %>% mutate(Outcome = factor(Outcome, levels = rev(Outcome))) %>%
		ggplot(aes(beta, Outcome)) + geom_vline(xintercept = 0, linetype = 2) +
		geom_segment(aes(x = beta - 1.96 * se, xend = beta + 1.96 * se, yend = Outcome)) + geom_point() +
		labs(title = "a. Adjusted preference association", x = paste("Change in preference per 1 SD", active_label), y = NULL) + theme_classic(11)
} else p3a <- empty_panel("a. Adjusted preference association")
pref_profile_plot <- pref_profile %>% filter(Variable != "pref_health_index", is.finite(Mean))
if(nrow(pref_profile_plot)) {
	p3b <- pref_profile_plot %>%
		ggplot(aes(score_decile, Mean, group = Outcome)) + geom_line() + geom_point() + facet_wrap( ~ Outcome, scales = "free_y") +
		labs(title = "b. Preference gradient across score deciles", x = paste(active_label, "decile"), y = "Mean preference (1-5)") + theme_classic(9)
} else p3b <- empty_panel("b. Preference gradient across score deciles")
if(any(is.finite(pref_corr$Spearman_rho))) {
	p3c <- pref_corr %>% filter(is.finite(Spearman_rho)) %>%
		ggplot(aes(Spearman_rho, reorder(Criterion, Spearman_rho))) + geom_vline(xintercept = 0, linetype = 2) +
		geom_segment(aes(x = 0, xend = Spearman_rho, yend = reorder(Criterion, Spearman_rho))) + geom_point() +
		coord_cartesian(xlim = c( - 1, 1)) + labs(title = "c. Independent convergent validity", x = "Spearman rho", y = NULL) + theme_classic(11)
} else p3c <- empty_panel("c. Independent convergent validity")
save_plot((p3a | p3c) / p3b + plot_annotation(title = "CHNS: intake-based score validated against independent food preferences"),
			"Figprimary.food_preference_validation.png", 14, 10.5, 450)
write_xlsx(list(adjusted_associations = pref_fit, decile_profile = pref_profile, correlations = pref_corr), "Figprimary.food_preference_validation.out.xlsx")


# 🚩 8) FigS4 independent processing/preparation/seasoning/fat-source validation
step_header("CHNS FigS4: processing, preparation and CFCT subclass validation")
proc_cov <- tibble(
	Variable = c("V39C processed-food classification", "V42 preparation method", "CFCT6 class/subclass map"),
	Record_coverage = c(mean(is.finite(food0$PROCESSED_CODE)), mean(is.finite(food0$PREP_METHOD)), mean(!is.na(food0$CLASS2))),
	Gram_coverage = c(sum(food0$GRAMS[is.finite(food0$PROCESSED_CODE)], na.rm = TRUE) / sum(food0$GRAMS, na.rm = TRUE),
				sum(food0$GRAMS[is.finite(food0$PREP_METHOD)], na.rm = TRUE) / sum(food0$GRAMS, na.rm = TRUE), class_cov)
)
proc_vars <- c("processed_food_g", "upf_class_proxy_g", "deepfried_g", "sugary_drinks_week", "carbonated_drink_g", "salt_msg_proxy_g", "fat_proxy_ratio", "alcohol_beverage_g")
proc_labs <- c(processed_food_g = "V39C processed food (g/day)", upf_class_proxy_g = "CFCT processed/refined proxy (g/day)", deepfried_g = "Deep-fried food (g/day)",
			sugary_drinks_week = "Questionnaire sugary drinks/week", carbonated_drink_g = "CFCT carbonated drinks (g/day)", salt_msg_proxy_g = "Salt/MSG/other seasoning proxy (g/day)",
			fat_proxy_ratio = "Plant/(plant+animal) fat-source ratio", alcohol_beverage_g = "Alcoholic beverages (g/day)")
proc_profile <- dat %>% filter(is.finite(.data[[active_maha]])) %>% mutate(decile = ntile(.data[[active_maha]], 10)) %>%
	group_by(decile) %>% summarise(across(all_of(proc_vars), ~ mean(.x, na.rm = TRUE)), .groups = "drop") %>%
	pivot_longer( - decile, names_to = "Variable", values_to = "Value") %>% mutate(Metric = unname(proc_labs[Variable]))
proc_corr <- tibble(Metric = unname(proc_labs[proc_vars]), Variable = proc_vars,
			Spearman_rho = vapply(proc_vars, function(v)safe_cor(dat[[active_maha]], dat[[v]]), numeric(1)))
p4a <- proc_cov %>% pivot_longer(c(Record_coverage, Gram_coverage), names_to = "Coverage", values_to = "Value") %>%
	ggplot(aes(Variable, Value, fill = Coverage)) + geom_col(position = "dodge") + scale_y_continuous(labels = scales::percent) +
	labs(title = "a. Independent field/map coverage", x = NULL, y = "Coverage", fill = NULL) + theme_classic(10) + theme(axis.text.x = element_text(angle = 20, hjust = 1))
p4b <- ggplot(proc_profile, aes(decile, Value)) + geom_line() + geom_point() + facet_wrap( ~ Metric, scales = "free_y", ncol = 4) +
	labs(title = "b. External behavioral gradient across MAHA-CN deciles", x = paste(active_label, "decile"), y = NULL) + theme_classic(9)
p4c <- ggplot(proc_corr, aes(Spearman_rho, reorder(Metric, Spearman_rho))) + geom_vline(xintercept = 0, linetype = 2) + geom_segment(aes(x = 0, xend = Spearman_rho, yend = reorder(Metric, Spearman_rho))) + geom_point() +
	coord_cartesian(xlim = c( - 1, 1)) + labs(title = "c. Correlation with independent intake/behavior proxies", x = "Spearman rho", y = NULL) + theme_classic(10)
save_plot((p4a | p4c) / p4b + plot_annotation(title = "CHNS: processing, cooking and food-class construct checks"),
			"FigS4.processing_preparation_validation.png", 15, 10.5, 450)
write_xlsx(list(field_coverage = proc_cov, decile_profile = proc_profile, correlations = proc_corr), "FigS4.processing_preparation_validation.out.xlsx")


# 🚩 11) FigS5 socioeconomic / geographic gradients
step_header("CHNS FigS5: socioeconomic and geographic gradients")
socio_long <- bind_rows(
	if("log_income" %in% names(dat)) tibble(Factor = "Household income", rho = safe_cor(dat[[active_maha]], dat$log_income)) else NULL,
	if("education_years" %in% names(dat)) tibble(Factor = "Estimated schooling, y", rho = safe_cor(dat[[active_maha]], dat$education_years)) else NULL,
	tibble(Factor = "Urban residence (1/2 code)", rho = safe_cor(dat[[active_maha]], num(dat$urban)))
)
prov_score <- dat %>% group_by(province) %>% summarise(N = n(), mean_score = mean(.data[[active_maha]], na.rm = TRUE),
			mean_chdi = mean(diet.chdi.sum, na.rm = TRUE), .groups = "drop") %>% filter(N >= 30)
p5a <- ggplot(socio_long, aes(rho, reorder(Factor, rho))) + geom_vline(xintercept = 0, linetype = 2) + geom_segment(aes(x = 0, xend = rho, yend = reorder(Factor, rho))) + geom_point() +
	coord_cartesian(xlim = c( - 1, 1)) + labs(title = "a. Socioeconomic gradient", x = "Spearman rho", y = NULL) + theme_classic(11)
p5b <- ggplot(prov_score, aes(mean_score, reorder(as.character(province), mean_score))) + geom_point(aes(size = N)) +
	labs(title = "b. Provincial variation", x = paste("Mean", active_label), y = "Province code", size = "N") + theme_classic(11)
save_plot(p5a | p5b, "FigS5.socioeconomic_geographic.png", 12, 7, 420)
write_xlsx(list(socioeconomic = socio_long, province = prov_score), "FigS5.socioeconomic_geographic.out.xlsx")


# 🚩 12) FigS6 cardiometabolic validation
step_header("CHNS FigS6: cardiometabolic validation")
outcomes <- c(hypertension = "Hypertension", t2dm_self = "Diabetes", mi_self = "Myocardial infarction",
			stroke_self = "Stroke", cancer_self = "Cancer", obesity = "Obesity")
validation_scores <- validation_score_set
phen_all <- map_dfr(names(outcomes), function(y){
	imap_dfr(validation_scores, function(v, lab){
		fit_bin_score(dat, y, dat[[v]], covs = c("age", "female", "province", "urban", "log_income", "education_years", "current_smoke", "kcal")) %>%
			mutate(Outcome = outcomes[[y]], Score = lab, Analysis = "Cross-sectional at baseline")
	})
})

incident_outcomes <- c(inc_hypertension = "Incident hypertension", inc_obesity = "Incident obesity", inc_t2dm = "Incident diabetes",
			inc_mi = "Incident myocardial infarction", inc_stroke = "Incident stroke", inc_cancer = "Incident cancer")
inc_all <- map_dfr(names(incident_outcomes), function(y){
	imap_dfr(validation_scores, function(v, lab){
		fit_bin_score(dat, y, dat[[v]], covs = c("age", "female", "province", "urban", "log_income", "education_years", "current_smoke", "bmi", "kcal")) %>%
			mutate(Outcome = incident_outcomes[[y]], Score = lab, Analysis = paste0("Prospective to ", followup_end))
	})
})

change_outcomes <- c(delta_bmi = "Change in BMI", delta_sbp = "Change in SBP (mm Hg)", delta_dbp = "Change in DBP (mm Hg)")
change_baseline <- c(delta_bmi = "bmi", delta_sbp = "sbp", delta_dbp = "dbp")
change_all <- map_dfr(names(change_outcomes), function(y){
	imap_dfr(validation_scores, function(v, lab){
		cv <- c("age", "female", "province", "urban", "log_income", "education_years", "current_smoke", "kcal", change_baseline[[y]])
		fit_linear_score(dat, y, dat[[v]], covs = cv) %>%
			mutate(Outcome = change_outcomes[[y]], Score = lab, Analysis = paste0("Change to ", followup_end))
	})
})

p6a <- phen_all %>% filter(is.finite(OR)) %>% mutate(Outcome = factor(Outcome, levels = rev(unique(Outcome)))) %>%
	ggplot(aes(OR, Outcome, shape = Score)) + geom_vline(xintercept = 1, linetype = 2) +
	geom_segment(aes(x = CI_low, xend = CI_high, yend = Outcome), position = position_dodge(width = .35)) +
	geom_point(position = position_dodge(width = .35), size = 2.1) +
	labs(title = "a. Baseline cardiometabolic face validity", x = "Adjusted OR per 1 SD healthier score", y = NULL, shape = NULL) +
	theme_classic(9) + theme(legend.position = "bottom")
if(any(is.finite(inc_all$OR))) {
	p6b <- inc_all %>% filter(is.finite(OR)) %>% mutate(Outcome = factor(Outcome, levels = rev(unique(Outcome)))) %>%
		ggplot(aes(OR, Outcome, shape = Score)) + geom_vline(xintercept = 1, linetype = 2) +
		geom_segment(aes(x = CI_low, xend = CI_high, yend = Outcome), position = position_dodge(width = .35)) +
		geom_point(position = position_dodge(width = .35), size = 2.1) +
		labs(title = paste0("b. Prospective disease validation to ", followup_end), x = "Adjusted OR per 1 SD healthier score", y = NULL, shape = NULL) +
		theme_classic(9) + theme(legend.position = "bottom")
} else {
	p6b <- ggplot() + annotate("text", x = 0, y = 0, label = paste0("No usable prospective PEXAM disease outcomes at ", followup_end), size = 4) +
		xlim( - 1, 1) + ylim( - 1, 1) + labs(title = "b. Prospective disease validation") + theme_void()
}
if(any(is.finite(change_all$beta))) {
	p6c <- change_all %>% filter(is.finite(beta)) %>% mutate(Outcome = factor(Outcome, levels = rev(unique(Outcome)))) %>%
		ggplot(aes(beta, Outcome, shape = Score)) + geom_vline(xintercept = 0, linetype = 2) +
		geom_segment(aes(x = CI_low, xend = CI_high, yend = Outcome), position = position_dodge(width = .35)) +
		geom_point(position = position_dodge(width = .35), size = 2.1) +
		labs(title = paste0("c. Prospective continuous change to ", followup_end), x = "Adjusted change per 1 SD healthier score", y = NULL, shape = NULL) +
		theme_classic(9) + theme(legend.position = "bottom")
} else {
	p6c <- ggplot() + annotate("text", x = 0, y = 0, label = "No usable paired BMI/BP measurements", size = 4) + xlim( - 1, 1) + ylim( - 1, 1) + theme_void()
}
save_plot((p6a | p6b) / p6c, "FigS6.cardiometabolic_validation.png", 14, 11, 420)
write_xlsx(list(cross_sectional = phen_all, prospective_incident = inc_all, prospective_change = change_all), "FigS6.cardiometabolic_validation.out.xlsx")


# 🚩 13) FigS7 distributions/concordance
step_header("CHNS FigS7: score distributions and concordance")
score_vars <- validation_score_set[names(validation_score_set) != "CHDI"]
score_df <- dat %>% dplyr::select(all_of(unname(score_vars)))
names(score_df) <- names(score_vars)
score_long <- score_df %>% mutate(across(everything(), score100)) %>% pivot_longer(everything(), names_to = "Diet", values_to = "Score")
density_curve <- score_long %>% filter(is.finite(Score)) %>% group_by(Diet) %>% group_modify( ~ {
	dn <- density(.x$Score, from = 0, to = 100, n = 256, na.rm = TRUE)
	tibble(score = dn$x, density = dn$y)
}) %>% ungroup()
cm <- cor(score_df, use = "pairwise.complete.obs", method = "spearman")
p7a <- ggplot(score_long, aes(Score)) + geom_density() + facet_wrap( ~ Diet, scales = "free_y", ncol = 1) +
	labs(title = "a. Score distributions", x = "Rescaled 0-100", y = "Density") + theme_classic(10)
cd <- as.data.frame(as.table(cm))
p7b <- ggplot(cd, aes(Var2, Var1, fill = Freq)) + geom_tile() + geom_text(aes(label = sprintf("%.2f", Freq))) +
	scale_fill_gradient2(midpoint = 0, limits = c( - 1, 1)) + labs(title = "b. Spearman concordance", x = NULL, y = NULL, fill = NULL) + theme_minimal(10) +
	theme(axis.text.x = element_text(angle = 35, hjust = 1), panel.grid = element_blank())
save_plot(p7a | p7b, "FigS7.score_distribution_concordance.png", 11, 8, 420)
write_xlsx(list(summary = score_long %>% group_by(Diet) %>% summarise(N = sum(is.finite(Score)), mean = mean(Score, na.rm = TRUE), sd = sd(Score, na.rm = TRUE), .groups = "drop"), density_curve = density_curve,
				spearman = as.data.frame(cm)), "FigS7.score_distribution_concordance.out.xlsx")


# 🚩 12) FigS8 food-map QC / unresolved code-level components
step_header("CHNS FigS8: FOODCODE map QC")
class_prev <- food0 %>% group_by(CLASS2, CLASS_NAME_CN) %>% summarise(grams = sum(GRAMS, na.rm = TRUE), records = n(), .groups = "drop") %>%
	mutate(gram_fraction = grams / sum(grams, na.rm = TRUE)) %>% arrange(desc(gram_fraction))
unmatched <- partial_map %>% filter(is.na(CLASS2))
missing_code <- tibble(Field = required_code_fields) %>% mutate(
	gram_coverage = ifelse(Field %in% map_field_coverage$field, map_field_coverage$gram_coverage[match(Field, map_field_coverage$field)], 0),
	Status = ifelse(gram_coverage >= full_map_min_cov, "coverage meets threshold", "below threshold / unavailable")
)
p8a <- class_prev %>% slice_head(n = 15) %>% mutate(Group = factor(paste0(CLASS2, " ", CLASS_NAME_CN), levels = rev(paste0(CLASS2, " ", CLASS_NAME_CN)))) %>%
	ggplot(aes(gram_fraction, Group)) + geom_col() + scale_x_continuous(labels = scales::percent) + labs(title = "a. Baseline food grams by CFCT class", x = "Fraction of food grams", y = NULL) + theme_classic(10)
p8b <- map_field_coverage
if(!nrow(p8b)) p8b <- missing_code %>% transmute(field = Field, gram_coverage = gram_coverage)
p8bplot <- p8b %>% mutate(field = factor(field, levels = rev(field))) %>% ggplot(aes(gram_coverage, field)) + geom_col() + scale_x_continuous(labels = scales::percent, limits = c(0, 1)) +
	labs(title = "b. Code-level nutrient/fine-group coverage", subtitle = paste0("Threshold = ", scales::percent(full_map_min_cov), "; ", analysis_mode), x = "Gram-weighted coverage", y = NULL) + theme_classic(9)
save_plot(p8a | p8bplot, "FigS8.food_mapping_qc.png", 13, 7.5, 420)
write_xlsx(list(foodcode_audit = partial_map, class_prevalence = class_prev, unmatched_codes = unmatched,
				foodcode_field_coverage = map_field_coverage, unresolved_foodcode_fields = missing_code,
				metadata = tibble(map_mode = analysis_mode, class_gram_coverage = class_cov, foodcode_map_coverage = full_map_coverage)),
			"FigS8.food_mapping_qc.out.xlsx")


# 🚩 15) Main mortality comparison
step_header("CHNS main mortality comparison")

dm <- dat %>% filter(followup_known == 1) %>% drop_na(death)
mort_scores <- validation_score_set

mort <- imap_dfr(mort_scores, function(v, lab) {
	fit_primary_cox(dm, dm[[v]]) %>%
		mutate(Diet = lab, Effect = HR, Effect_measure = "HR", Method = "Cox", .before = 1)
})

# Reproducibility audit only.  Event/N changes are now expected when the user
# changes wave, age, or follow-up definition, so no result-direction guard is
# allowed to stop the analysis.
.guard_maha <- mort %>% filter(Diet == "MAHA") %>% slice(1)
.guard_cmp <- mort %>% filter(Diet %in% c("DASH", "MIND", "MEDI"))
primary_reproducibility_audit <- tibble(
	Expected_configuration = "2011 baseline, age >=40, end-wave complete case",
	Expected_N = 7387L, Expected_deaths = 219L,
	Observed_N = if(nrow(.guard_maha)).guard_maha$N[1]else NA_integer_,
	Observed_deaths = if(nrow(.guard_maha)).guard_maha$Events[1]else NA_integer_,
	Exact_match = isTRUE(nrow(.guard_maha) == 1L && .guard_maha$N[1] == 7387L && .guard_maha$Events[1] == 219L),
	Established_all_HR_below_1 = isTRUE(nrow(.guard_cmp) == 3L && all(is.finite(.guard_cmp$HR)) && all(.guard_cmp$HR < 1)),
	Note = "Audit only; neither N/events nor association direction is used to select a model."
)

quartiles <- imap_dfr(mort_scores, function(v, lab) {
	fit_primary_quartile(dm, dm[[v]]) %>%
		mutate(Diet = lab, Effect_measure = "HR", Method = "Cox", .before = 1)
})

pub_precision <- tibble::tibble(status = "Run diet.sh bayes chns to calculate precision and BF01")
pub_head2head <- bind_rows(
	fit_primary_head2head(dm, "diet.dash.sum", "DASH"),
	fit_primary_head2head(dm, "diet.mind.sum", "MIND"),
	fit_primary_head2head(dm, "diet.medi.sum", "MEDI")
)

reference_event_n <- if(nrow(mort) && any(is.finite(mort$Events))) max(mort$Events, na.rm = TRUE) else NA_integer_
pM <- mort %>% mutate(Diet = factor(Diet, levels = rev(names(mort_scores)))) %>%
	ggplot(aes(HR, Diet)) + geom_vline(xintercept = 1, linetype = 2) +
	geom_segment(aes(x = CI_low, xend = CI_high, yend = Diet)) + geom_point(size = 2.7) +
	labs(title = "CHNS all-cause mortality — Model 1 (default)",
			 subtitle = paste0("Fixed 2011 baseline; reference algorithms; ", reference_event_n, " deaths in the complete-case model"),
			 x = "Adjusted HR per 1 SD healthier score", y = NULL) + theme_classic(12)
save_plot(pM, "Fig6.mortality_validation.png", 8.5, 6.5, 420)
write_xlsx(list(mortality = mort, quartiles = quartiles, precision = pub_precision, head2head = pub_head2head,
				model_definition = tibble(Model = "Model1_default_reference", Baseline = 2011L,
					Design = "Fixed baseline", Energy_covariate = TRUE,
					Note = "Exact submitted reference benchmark; Models 4-5 are prespecified high-event cumulative-average sensitivities"),
				reproducibility_audit = primary_reproducibility_audit),
			"Fig6.mortality_validation.out.xlsx")

sensitivity_specs <- tribble(
 ~ Analysis, ~ require_three_days, ~ no_bmi, ~ exclude_first_year, ~ exclude_baseline_disease,
	"Primary", FALSE, FALSE, FALSE, FALSE,
	"All 3 diet-record days", TRUE, FALSE, FALSE, FALSE,
	"Exclude deaths in first year", FALSE, FALSE, TRUE, FALSE,
	"Exclude baseline disease", FALSE, FALSE, FALSE, TRUE
)

mort_sens <- pmap_dfr(sensitivity_specs, function(Analysis, require_three_days, no_bmi, exclude_first_year, exclude_baseline_disease) {
	imap_dfr(mort_scores, function(v, lab) {
		fit_primary_cox(dm, dm[[v]], require_three_days = require_three_days, no_bmi = no_bmi,
							exclude_first_year = exclude_first_year, exclude_baseline_disease = exclude_baseline_disease) %>%
			mutate(Diet = lab, Analysis = Analysis, .before = 1)
	})
})

# Evaluate the prespecified mortality model grid without outcome-based selection.
grid_scores <- c(MAHA = "diet.maha.sum", DASH = "diet.dash.sum", MIND = "diet.mind.sum", MEDI = "diet.medi.sum")
# Focused mortality set: the submitted fixed-wave reference plus the selected
# high-event socioeconomic model and its BMI sensitivity.  No other alternative
# model is fitted or written to the final workbooks.
adjustment_specs <- list(
	Model4_socioeconomic = c("age", "female", "province", "urban", "log_income", "education_years"),
	Model5_bmi = c("age", "female", "province", "urban", "log_income", "education_years", "bmi")
)
final_model_roles <- c(
	Model4_socioeconomic = "Primary high-event sensitivity",
	Model5_bmi = "BMI / potential-overadjustment sensitivity"
)
final_model_notes <- c(
	Model4_socioeconomic = "Age >=40 cumulative-average exposure across 2004-2011; adjusts for province, urban residence, household income and education",
	Model5_bmi = "Adds BMI to the socioeconomic model; interpreted as a possible mediator/overadjustment sensitivity"
)
mortality_definitions <- tibble(
	Mortality_definition = c("end_wave", "last_contact"),
	Definition_label = c(
		"Observed at end wave or death",
		"Censor at last observed contact"
	),
	Publication_eligible = c(TRUE, TRUE),
	Definition_note = c(
		"Conservative complete follow-up definition used in the submitted draft",
		"Retains participants lost before the end wave and censors them at their last observed wave"
	)
)
sample_specs <- "four_score_common"

mortality_flow <- pmap_dfr(mortality_definitions, function(Mortality_definition, Definition_label, Publication_eligible, Definition_note) {
	ff <- add_primary_survival_fields(dat, Mortality_definition) %>% filter(is.finite(surv_time), surv_event %in% c(0, 1))
	tibble(
		Baseline_wave = baseline_wave, Minimum_age = min_age,
		Mortality_definition = Mortality_definition, Definition_label = Definition_label,
		Publication_eligible = Publication_eligible, Definition_note = Definition_note,
		N_before_covariates = nrow(ff), Deaths_before_covariates = sum(ff$surv_event == 1),
		Person_years_before_covariates = sum(ff$surv_time),
		Events_ge_target = sum(ff$surv_event == 1) >= mortality_event_target,
		Event_target = mortality_event_target
	)
})

mortality_grid <- pmap_dfr(mortality_definitions, function(Mortality_definition, Definition_label, Publication_eligible, Definition_note) {
	imap_dfr(adjustment_specs, function(covs, adj_lab) {
		map_dfr(sample_specs, function(sample_lab) {
			dd0 <- dat
			if(identical(sample_lab, "four_score_common")) {
				dd0 <- dd0 %>% filter(if_all(all_of(unname(grid_scores)), ~ is.finite(num(.x))))
			}
			imap_dfr(grid_scores, function(v, diet_lab) {
				fit_primary_cox(dd0, dd0[[v]], covs = covs, mortality_definition = Mortality_definition) %>%
					mutate(
						Baseline_wave = baseline_wave, Minimum_age = min_age, Diet = diet_lab,
						Mortality_definition = Mortality_definition, Definition_label = Definition_label,
						Publication_eligible = Publication_eligible, Definition_note = Definition_note,
						Adjustment = adj_lab, Sample_strategy = sample_lab,
						Event_target = mortality_event_target, Events_ge_target = Events >= mortality_event_target,
						.before = 1
					)
			})
		})
	})
})

mortality_pattern <- mortality_grid %>%
	select(Baseline_wave, Minimum_age, Mortality_definition, Definition_label, Publication_eligible,
			Adjustment, Sample_strategy, Diet, HR, P, N, Events, Person_years) %>%
	pivot_wider(names_from = Diet, values_from = c(HR, P, N, Events, Person_years)) %>%
	mutate(
		Established_all_protective = is.finite(HR_DASH) & is.finite(HR_MIND) & is.finite(HR_MEDI) &
			HR_DASH < 1 & HR_MIND < 1 & HR_MEDI < 1,
		Established_all_significant_protective = Established_all_protective &
			P_DASH < .05 & P_MIND < .05 & P_MEDI < .05,
		MAHA_attenuated_or_null_vs_all = is.finite(HR_MAHA) & is.finite(HR_DASH) & is.finite(HR_MIND) & is.finite(HR_MEDI) &
			HR_MAHA >= pmax(HR_DASH, HR_MIND, HR_MEDI),
		MAHA_null_or_risky = is.finite(HR_MAHA) & is.finite(P_MAHA) & (HR_MAHA >= 1 | P_MAHA >= .05),
		Descriptive_expected_pattern = Established_all_significant_protective & MAHA_null_or_risky & MAHA_attenuated_or_null_vs_all,
		Events_ge_target = Events_MAHA >= mortality_event_target,
		Selection_rule = "Diagnostic only: HR direction is never used to choose a primary model"
	)

p_grid <- mortality_grid %>%
	filter(Adjustment == "Model5_bmi", Sample_strategy == "four_score_common", is.finite(HR)) %>%
	mutate(Definition_label = factor(Definition_label, levels = rev(mortality_definitions$Definition_label)),
			Diet = factor(Diet, levels = c("MAHA", "DASH", "MIND", "MEDI"))) %>%
	ggplot(aes(HR, Definition_label, color = Diet)) +
	geom_vline(xintercept = 1, linetype = 2, color = "grey55") +
	geom_segment(aes(x = CI_low, xend = CI_high, yend = Definition_label), position = position_dodge(width = .62), linewidth = .65) +
	geom_point(position = position_dodge(width = .62), size = 2.5) +
	labs(title = paste0("CHNS mortality-definition grid: baseline ", baseline_wave, ", age >=", min_age),
			 subtitle = "Model 5 BMI sensitivity (no energy covariate); common sample across all four scores",
			 x = "HR per 1 SD healthier score", y = NULL, color = NULL) +
	theme_classic(11) + theme(legend.position = "bottom")
save_plot(p_grid, "mortality_model_grid.png", 11, 7.5, 420)
write_xlsx(list(flow = mortality_flow, results = mortality_grid, pattern_diagnostic = mortality_pattern),
			"mortality_model_grid.out.xlsx")

analysis_metadata <- tibble(
	schooling_definition = "Estimated schooling, y: A11 0=0; 11-16=1-6; 21-26=7-12; 27-29=10-12; 31-36=13-18; other codes missing. EDUYRS retained as recorded. A11-derived years are estimated, not individually measured.",
	primary_model = analysis_spec,
	baseline_wave = baseline_wave,
	followup_end = followup_end,
	min_age = min_age,
	min_diet_days = min_diet_days,
	association_method = "Cox",
	association_measure = "HR",
	primary_covariates = paste(primary_covariates, collapse = " + "),
	primary_complete_case_N = mort$N[mort$Diet == "MAHA"][1],
	primary_deaths = mort$Events[mort$Diet == "MAHA"][1]
)

benchmark <- tibble(
	N = analysis_metadata$primary_complete_case_N,
	Deaths = analysis_metadata$primary_deaths,
	Baseline = baseline_wave,
	Followup_end = followup_end
)

write_xlsx(
	list(
		mortality = mort,
		quartiles = quartiles,
		precision = pub_precision,
		head2head = pub_head2head,
		sensitivity = mort_sens,
		benchmark = benchmark,
		reproducibility_audit = primary_reproducibility_audit,
		analysis_metadata = analysis_metadata
	),
	"publication_inputs.out.xlsx"
)


# 🚩 16) Follow-up QC
fu <- dat %>% summarise(
	Baseline_N = n(),
	Known_followup = sum(followup_known == 1, na.rm = TRUE),
	Observed_in_end_wave_roster = sum(observed_at_followup == 1, na.rm = TRUE),
	Observed_after_baseline = sum(observed_after_baseline == 1, na.rm = TRUE),
	Deaths = sum(death_reported == 1, na.rm = TRUE),
	Unknown_at_end_wave = sum(followup_known != 1 | is.na(followup_known), na.rm = TRUE),
	Note = "Roster presence is contact evidence, not a direct alive-status indicator"
)
fu_prov <- dat %>% group_by(province) %>% summarise(
	N = n(), known = sum(followup_known == 1, na.rm = TRUE), known_pct = known / N,
	deaths = sum(death_reported == 1, na.rm = TRUE), .groups = "drop"
)
p9 <- ggplot(fu_prov, aes(known_pct, reorder(as.character(province), known_pct))) +
	geom_col() + scale_x_continuous(labels = percent) +
	labs(title = "Follow-up ascertainment by province", x = "Known vital status", y = "Province code") + theme_classic(11)
save_plot(p9, "FigS9.followup_qc.png", 8, 6.5, 420)
write_xlsx(list(overall = fu, by_province = fu_prov), "FigS9.followup_qc.out.xlsx")


# 🚩 17) Mortality sensitivity figure
p10 <- mort_sens %>%
	filter(Diet %in% c("MAHA", "DASH", "MIND", "MEDI"), is.finite(HR)) %>%
	mutate(Analysis = factor(Analysis, levels = rev(sensitivity_specs$Analysis))) %>%
	ggplot(aes(HR, Analysis, color = Diet)) +
	geom_vline(xintercept = 1, linetype = 2, color = "grey55") +
	geom_segment(aes(x = CI_low, xend = CI_high, yend = Analysis), position = position_dodge(width = .55), linewidth = .7) +
	geom_point(position = position_dodge(width = .55), size = 2.5) +
	labs(title = "CHNS all-cause mortality sensitivity analyses", x = "Adjusted HR per 1 SD healthier score", y = NULL, color = NULL) +
	theme_classic(11) + theme(legend.position = "bottom")
save_plot(p10, "FigS10.mortality_sensitivity.png", 10.5, 7.5, 420)
write_xlsx(list(results = mort_sens, metadata = analysis_metadata), "FigS10.mortality_sensitivity.out.xlsx")


# 🚩 18) Table 1, metadata, reusable dataset
active_3c <- paste0(sub("\\.sum$", "", active_maha), ".3c")
tab1 <- dat %>% mutate(Diet_group = .data[[active_3c]]) %>% group_by(Diet_group) %>% summarise(
	N = n(), age_mean = mean(age, na.rm = TRUE), female_pct = mean(female, na.rm = TRUE), BMI_mean = mean(bmi, na.rm = TRUE),
	kcal_mean = mean(kcal, na.rm = TRUE), income_mean = mean(income_pc, na.rm = TRUE), education_mean = mean(education_years, na.rm = TRUE),
	known_followup_pct = mean(followup_known == 1, na.rm = TRUE), mortality_pct = mean(death, na.rm = TRUE), .groups = "drop"
)
meta <- tibble(
	schooling_definition = analysis_metadata$schooling_definition,
	primary_model = analysis_spec,
	baseline_wave = baseline_wave, followup_end = followup_end, followup_gap = followup_end - baseline_wave,
	min_age = min_age, min_diet_days = min_diet_days, association_method = "Cox", active_score = active_label,
	nutr3 = src$NUTR3[["path"]], c12diet = src$C12DIET[["path"]], master = src$MASTER[["path"]],
	roster = src$ROSTER[["path"]], pexam = src$PEXAM[["path"]], pact = src$PACT[["path"]],
	food_map = ifelse(is.na(food_map_file), "", food_map_file),
	nutrient_map = ifelse(is.na(nutrient_map_file), "", nutrient_map_file),
	classification_mapping_coverage = class_cov, foodcode_mapping_coverage = full_map_coverage,
	code_map_available = code_map_available, enhanced_map_mode = enhanced_map_mode, full_map_mode = full_map_mode,
	foodcode_field_threshold = full_map_min_cov
)
write_xlsx(list(Table1 = tab1, metadata = meta, input_sources = source_manifest, benchmark = benchmark), "Table1.out.xlsx")
saveRDS(dat, cohort_file("analysis_dataset.rds"))


# 🚩 19) Focused high-event cumulative-average sensitivity
# Only age >=40 Models 4 and 5 are retained.  The earlier exploratory age,
# first-entry, landmark, minimal-adjustment and all-at-end grids are not fitted.
if(run_mortality_sweep == 1L && !mortality_sweep_child) {
	step_header("CHNS focused high-event mortality sensitivity")
	stale_focused_outputs <- cohort_file(c(
		"mortality_design_sweep.out.xlsx", "mortality_design_sweep.png",
		"mortality_competing_models.out.xlsx", "mortality_competing_all4.png",
		"mortality_competing_pairwise.png"
	))
	unlink(stale_focused_outputs[file.exists(stale_focused_outputs)], force = TRUE)

	parse_int_list <- function(x) {
		z <- suppressWarnings(as.integer(trimws(strsplit(x, ",", fixed = TRUE)[[1]])))
		sort(unique(z[is.finite(z)]))
	}
	eligible_sweep_waves <- common_baseline[common_baseline < followup_end]
	if(!identical(tolower(mortality_sweep_wave_req), "auto")) {
		eligible_sweep_waves <- intersect(eligible_sweep_waves, parse_int_list(mortality_sweep_wave_req))
	}
	sweep_ages <- mortality_model_age
	if(!length(eligible_sweep_waves)) {
		warning("No eligible common waves for CHNS mortality sweep.")
	} else {
		sweep_root <- file.path(maha_state_dir, "waves")
		dir.create(sweep_root, recursive = TRUE, showWarnings = FALSE)
		rscript_bin <- Sys.getenv("RSCRIPT_BIN", unset = Sys.which("Rscript"))
		if(!nzchar(rscript_bin)) stop("Cannot locate Rscript for mortality sweep.", call. = FALSE)

		sweep_runs <- list()
		run_manifest <- list()
		md5_or_blank <- function(x) if(length(x) == 1L && !is.na(x) && file.exists(x))unname(tools::md5sum(x))else""
		kk <- 0L
		for(aa in sweep_ages) for(ww in eligible_sweep_waves) {
			kk <- kk + 1L
			run_id <- paste0("wave", ww, "_age", aa)
			child_root <- file.path(sweep_root, run_id)
			child_aux <- file.path(child_root, "chns")
			child_grid_file <- file.path(child_aux, "mortality_model_grid.out.xlsx")
			child_data_file <- file.path(child_aux, "analysis_dataset.rds")
			child_log <- file.path(sweep_root, paste0(run_id, ".log"))
			signature_file <- file.path(child_root, "run_signature.txt")
			run_signature <- paste(md5_or_blank(.this_file), ww, aa, followup_end, min_diet_days, full_map_min_cov,
						require_nutrient_map, require_enhanced_map,
						md5_or_blank(food_map_file), md5_or_blank(nutrient_map_file), sep = "|")
			cache_ok <- !tolower(Sys.getenv("DIET_REPLACE", "TRUE")) %in% c("true", "1", "yes") &&
				file.exists(child_grid_file) && file.exists(child_data_file) && file.exists(signature_file) &&
				identical(readLines(signature_file, warn = FALSE, n = 1L), run_signature)

			use_parent <- identical(as.integer(ww), as.integer(baseline_wave)) && identical(as.integer(aa), as.integer(min_age))
			status <- 0L
			if(use_parent) {
				child_grid_file <- cohort_file("mortality_model_grid.out.xlsx")
				child_data_file <- cohort_file("analysis_dataset.rds")
			} else if(!cache_ok) {
				dir.create(child_root, recursive = TRUE, showWarnings = FALSE)
				env_child <- c(
					paste0("MAHA_OUTDIR=", shQuote(child_root)),
					paste0("MAHA_AUXDIR=", shQuote(child_aux)),
                    paste0("DIET_MAHA_STATE=", shQuote(file.path(child_root, "state"))),
                    paste0("DIET_DISCO_PART=", shQuote(file.path(child_root, "scores.rds"))),
					paste0("CHNS_BASELINE_WAVE=", ww),
					paste0("CHNS_FOLLOWUP_END=", followup_end),
					paste0("CHNS_MIN_AGE=", aa),
					paste0("CHNS_MIN_DIET_DAYS=", min_diet_days),
					paste0("CHNS_MIN_MAP_COVERAGE=", full_map_min_cov),
					paste0("CHNS_REQUIRE_NUTRIENT_MAP=", require_nutrient_map),
					paste0("CHNS_REQUIRE_ENHANCED_MAP=", require_enhanced_map),
					paste0("CHNS_PROJECT_DIR=", shQuote(project_dir)),
					paste0("CHNS_FOOD_MAP=", shQuote(food_map_file)),
					paste0("CHNS_NUTRIENT_MAP=", shQuote(nutrient_map_file)),
					paste0("RSCRIPT_BIN=", shQuote(rscript_bin)),
					"CHNS_MORTALITY_SWEEP_CHILD=1",
					"CHNS_RUN_MORTALITY_SWEEP=0",
					"MAHA_NULL_MAX=1"
				)
				status <- tryCatch(
					system2(rscript_bin, args = c(shQuote(maha_runner), "chns", "all"), env = env_child, stdout = child_log, stderr = child_log),
					error = function(e) 999L
				)
				if(identical(as.integer(status), 0L) && file.exists(child_grid_file) && file.exists(child_data_file))
					writeLines(run_signature, signature_file, useBytes = TRUE)
			}
			ok <- identical(as.integer(status), 0L) && file.exists(child_grid_file) && file.exists(child_data_file)
			failure_tail <- if(!ok && file.exists(child_log)) {
				paste(utils::tail(readLines(child_log, warn = FALSE), 20L), collapse = " | ")
			} else ""
			run_manifest[[kk]] <- tibble(
				Run = run_id, Baseline_wave = ww, Minimum_age = aa, Followup_end = followup_end,
				Status = ifelse(ok, "complete", paste0("failed_status_", status)),
				Grid_file = child_grid_file, Dataset_file = child_data_file, Log_file = child_log,
				Failure_tail = failure_tail
			)
			if(ok) sweep_runs[[run_id]] <- list(wave = ww, age = aa, grid = child_grid_file, data = child_data_file)
			cat("[MORTALITY SWEEP] ", run_id, " : ", ifelse(ok, "complete", "failed"), "\n", sep = "")
			if(!ok && nzchar(failure_tail)) cat("[MORTALITY SWEEP FAILURE] ", run_id, " : ", failure_tail, "\n", sep = "")
		}
		run_manifest <- bind_rows(run_manifest)

		# Build the single prespecified age >=40 repeated-exposure panel.
		panel_by_age <- set_names(sweep_ages, as.character(sweep_ages)) %>% map(function(panel_age_i) {
			panel_runs_i <- sweep_runs[vapply(sweep_runs, function(x)x$age == panel_age_i, logical(1))]
			imap_dfr(panel_runs_i, function(rr, run_id) {
				readRDS(rr$data) %>% mutate(exposure_wave = as.integer(rr$wave), sweep_run = run_id)
			}) %>% arrange(IDIND, exposure_wave) %>% distinct(IDIND, exposure_wave, .keep_all = TRUE)
		})
		panel_by_age <- panel_by_age[vapply(panel_by_age, nrow, integer(1)) > 0L]
		if(!length(panel_by_age)) stop("No successful datasets were available for a repeated-wave panel.", call. = FALSE)
		available_panel_ages <- as.integer(names(panel_by_age))
		if(!mortality_model_age %in% available_panel_ages) {
			replacement_age <- available_panel_ages[which.min(abs(available_panel_ages - mortality_model_age))]
			warning("Requested named-model age ", mortality_model_age, " unavailable; using age ", replacement_age, ".")
			panel_age <- replacement_age
		} else panel_age <- mortality_model_age
		panel <- panel_by_age[[as.character(panel_age)]]

		cummean_known <- function(x) {
			x <- num(x); n <- cumsum(is.finite(x)); s <- cumsum(ifelse(is.finite(x), x, 0));
			out <- s / pmax(n, 1); out[n == 0] <- NA_real_; out
		}
		focused_algorithm_score_set <- algorithm_score_set[c(
			"MAHA", "MAHA-density", "MAHA-residual", "DASH-reference", "MIND", "MEDI"
		)]
		panel_cum_by_age <- map(panel_by_age, function(panel_i) {
			sweep_score_vars_i <- unique(c(unname(grid_scores), unname(focused_algorithm_score_set)))
			sweep_score_vars_i <- sweep_score_vars_i[sweep_score_vars_i %in% names(panel_i)]
			panel_i %>% group_by(IDIND) %>% arrange(exposure_wave, .by_group = TRUE) %>%
				mutate(across(all_of(sweep_score_vars_i), cummean_known, .names = "{.col}_cum")) %>% ungroup()
		})
		panel_cum <- panel_cum_by_age[[as.character(panel_age)]]

		build_time_updated <- cache_identical_inputs(function(d, mortality_definition) {
			d %>% group_by(IDIND) %>% arrange(exposure_wave, .by_group = TRUE) %>% group_modify( ~ {
				g <- .x
				dd <- c(num(g$DOD_Y), num(g$DOD_RPT)); dd <- dd[is.finite(dd) & dd > min(g$exposure_wave) & dd <= followup_end]
				dy <- if(length(dd))min(dd)else NA_real_
				lc0 <- num(g$last_contact_year); lc <- if(any(is.finite(lc0)))max(lc0, na.rm = TRUE)else NA_real_
				if(identical(mortality_definition, "end_wave")) {
					eligible <- any(g$observed_at_followup == 1, na.rm = TRUE) | is.finite(dy); analysis_end <- followup_end
				} else if(identical(mortality_definition, "last_contact")) {
					eligible <- is.finite(dy) | (is.finite(lc) & lc > min(g$exposure_wave)); analysis_end <- if(is.finite(dy))dy else lc
				} else stop("Unknown mortality definition: ", mortality_definition, call. = FALSE)
				if(!eligible || !is.finite(analysis_end)) return(g[0, , drop = FALSE])
				if(is.finite(dy)) analysis_end <- min(analysis_end, dy)
				g <- g %>% filter(exposure_wave < analysis_end)
				if(!nrow(g)) return(g)
				nxt <- dplyr::lead(g$exposure_wave, default = analysis_end)
				g$tstart <- g$exposure_wave
				g$tstop <- pmin(nxt, analysis_end)
				g$surv_event <- as.integer(is.finite(dy) & g$tstart < dy & g$tstop == dy)
				g %>% filter(tstop > tstart)
			}) %>% ungroup()
		})

		fit_time_updated_score <- function(d, score_var, covs, mortality_definition) {
			dd <- build_time_updated(d, mortality_definition)
			dd$score_z <- zstd(dd[[paste0(score_var, "_cum")]])
			covs_use <- good_covariates(dd, intersect(covs, names(dd)))
			dd <- dd %>% select(all_of(c("IDIND", "tstart", "tstop", "surv_event", "score_z", covs_use))) %>% drop_na()
			empty <- function() tibble(N = n_distinct(dd$IDIND), Rows = nrow(dd), Events = sum(dd$surv_event == 1),
				Person_years = sum(dd$tstop - dd$tstart), beta = NA_real_, se = NA_real_, HR = NA_real_, CI_low = NA_real_, CI_high = NA_real_, P = NA_real_,
				Covariates = paste(covs_use, collapse = " + "))
			if(nrow(dd) < 100 || n_distinct(dd$IDIND) < 100 || sum(dd$surv_event == 1) < 10 || length(unique(dd$score_z)) < 2)return(empty())
			fm <- as.formula(paste0("Surv(tstart,tstop,surv_event) ~ score_z",
				if(length(covs_use))paste0(" + ", paste(covs_use, collapse = " + "))else"", " + cluster(IDIND)"))
			fit <- tryCatch(coxph(fm, data = dd, ties = "efron"), error = function(e)NULL)
			if(is.null(fit))return(empty())
			tt <- broom::tidy(fit) %>% filter(term == "score_z")
			if(!nrow(tt))return(empty())
			se_beta <- if ("robust.se" %in% names(tt)) {
				tt$robust.se[1]
			} else {
				tt$std.error[1]
			}
			tibble(N = n_distinct(dd$IDIND), Rows = nrow(dd), Events = sum(dd$surv_event == 1), Person_years = sum(dd$tstop - dd$tstart),
				beta = tt$estimate[1], se = se_beta, HR = exp(tt$estimate[1]),
				CI_low = exp(tt$estimate[1] - 1.96 * se_beta), CI_high = exp(tt$estimate[1] + 1.96 * se_beta),
				P = tt$p.value[1], Covariates = paste(covs_use, collapse = " + "))
		}

		fit_time_updated_head2head <- function(d, maha_var, trad_var, trad_lab, covs, mortality_definition = "end_wave") {
			dd <- build_time_updated(d, mortality_definition)
			maha_cum <- paste0(maha_var, "_cum"); trad_cum <- paste0(trad_var, "_cum")
			if(!all(c(maha_cum, trad_cum) %in% names(dd))) return(tibble())
			dd$maha_z <- zstd(dd[[maha_cum]]); dd$trad_z <- zstd(dd[[trad_cum]])
			covs_use <- good_covariates(dd, intersect(covs, names(dd)))
			dd <- dd %>% select(all_of(c("IDIND", "tstart", "tstop", "surv_event", "maha_z", "trad_z", covs_use))) %>% drop_na()
			if(nrow(dd) < 100 || n_distinct(dd$IDIND) < 100 || sum(dd$surv_event == 1) < 10)return(tibble())
			fm <- as.formula(paste0("Surv(tstart,tstop,surv_event) ~ maha_z + trad_z",
				if(length(covs_use))paste0(" + ", paste(covs_use, collapse = " + "))else"", " + cluster(IDIND)"))
			fit <- tryCatch(coxph(fm, data = dd, ties = "efron"), error = function(e)NULL)
			if(is.null(fit))return(tibble())
			cf <- coef(fit); vc <- vcov(fit)
			if(!all(c("maha_z", "trad_z") %in% names(cf)))return(tibble())
			b1 <- unname(cf[["maha_z"]]); b2 <- unname(cf[["trad_z"]]); dff <- b1 - b2
			se_diff <- sqrt(vc["maha_z", "maha_z"] + vc["trad_z", "trad_z"] - 2 * vc["maha_z", "trad_z"])
			tibble(Comparison = paste0("MAHA vs ", trad_lab), Comparator = trad_lab,
				N = n_distinct(dd$IDIND), Events = sum(dd$surv_event == 1), HR_MAHA_mutual = exp(b1), HR_comparator_mutual = exp(b2),
				HR_ratio_MAHA_vs_comparator = exp(dff), CI_low_ratio = exp(dff - 1.96 * se_diff), CI_high_ratio = exp(dff + 1.96 * se_diff),
				P_equal_coefficients = 2 * pnorm( - abs(dff / se_diff)), Covariates = paste(covs_use, collapse = " + "))
		}

		time_updated_results <- imap_dfr(panel_cum_by_age, function(panel_cum_i, age_lab) {
			pmap_dfr(mortality_definitions, function(Mortality_definition, Definition_label, Publication_eligible, Definition_note) {
				imap_dfr(adjustment_specs, function(covs, adj_lab) {
					map_dfr(sample_specs, function(sample_lab) {
						dd0 <- panel_cum_i
						if(identical(sample_lab, "four_score_common")) dd0 <- dd0 %>%
							filter(if_all(all_of(paste0(unname(grid_scores), "_cum")), ~ is.finite(num(.x))))
						imap_dfr(grid_scores, function(v, diet_lab) {
							fit_time_updated_score(dd0, v, covs, Mortality_definition) %>% mutate(
								Model_id = paste0("time_updated_age", age_lab, "_", Mortality_definition, "_", adj_lab, "_", sample_lab),
								Design = "Chen-like time-updated cumulative average", Exposure_waves = paste(sort(unique(panel_cum_i$exposure_wave)), collapse = ","),
								Baseline_wave = NA_integer_, Minimum_age = as.integer(age_lab), Diet = diet_lab,
								Mortality_definition = Mortality_definition, Definition_label = Definition_label,
								Publication_eligible = Publication_eligible, Definition_note = Definition_note,
								Adjustment = adj_lab, Sample_strategy = sample_lab, Event_target = mortality_event_target,
								Events_ge_target = Events >= mortality_event_target,
								Method_note = "Start-stop Cox; cumulative mean updated only with current/prior waves; robust SE clustered by participant", .before = 1)
						})
					})
				})
			})
		})

		# Model 1 is the submitted 2011 analysis; Models 4-5 are the only retained
		# high-event cumulative-average sensitivities.
		model0_name <- "Model1_default_reference"
		model_definitions <- bind_rows(
			tibble(Model = model0_name, Design = "Fixed 2011 baseline", Adjustment_order = 1L,
				Minimum_age = min_age, Exposure_waves = as.character(baseline_wave),
				Covariates = paste(primary_covariates, collapse = " + "), Energy_covariate = TRUE,
				Role = "Default / submitted primary model",
				Note = paste0("Exact reproduction of Fig6; reference score algorithms; ", reference_event_n, "-death benchmark")),
			imap_dfr(adjustment_specs, function(covs, lab) tibble(
				Model = lab, Design = "Chen-like time-updated cumulative average", Adjustment_order = match(lab, names(adjustment_specs)) + 1L,
				Minimum_age = panel_age, Exposure_waves = paste(sort(unique(panel_cum$exposure_wave)), collapse = ","),
				Covariates = paste(covs, collapse = " + "), Energy_covariate = FALSE,
				Role = unname(final_model_roles[[lab]]), Note = unname(final_model_notes[[lab]])
			))
		) %>% mutate(
			Selection_rule = "Retained before final reporting for design coverage and interpretability; never selected from HR direction or P value")
		if(nrow(model_definitions) > 6L)
			stop("Focused CHNS analysis may retain no more than six named models.", call. = FALSE)

		model0_results <- mort %>% filter(Diet %in% c("MAHA", "DASH", "MIND", "MEDI")) %>%
			transmute(Model = model0_name, Design = "Fixed 2011 baseline",
				Score_algorithm = ifelse(Diet == "DASH", "DASH-reference", Diet), Diet, N, Events, Person_years,
				beta, se, HR, CI_low, CI_high, P, Covariates, Rows = NA_integer_)

		# Anchor alternatives to one reference-comparable starting cohort.
		chen_core <- c(MAHA = "diet.maha.sum", DASH = "diet.dash.sum", MIND = "diet.mind.sum", MEDI = "diet.medi.sum")
		chen_core_cum <- paste0(unname(chen_core), "_cum")
		chen_common <- panel_cum %>% filter(if_all(all_of(chen_core_cum), ~ is.finite(num(.x))))
		chen_model_results <- imap_dfr(adjustment_specs, function(covs, model_lab) {
			imap_dfr(focused_algorithm_score_set, function(v, score_lab) {
				fit_time_updated_score(chen_common, v, covs, "end_wave") %>% mutate(
					Model = model_lab, Design = "Chen-like time-updated cumulative average", Score_algorithm = score_lab,
					Diet = case_when(startsWith(score_lab, "MAHA") ~ "MAHA", startsWith(score_lab, "DASH") ~ "DASH", TRUE ~ score_lab),
					.before = 1)
			})
		})
		model_results_0N <- bind_rows(model0_results, chen_model_results)

		model_primary_results <- model_results_0N %>% filter(
			(Model == model0_name & Score_algorithm %in% c("MAHA", "DASH-reference", "MIND", "MEDI")) |
			(Model != model0_name & Score_algorithm %in% c("MAHA", "DASH-reference", "MIND", "MEDI"))
		)
		maha_algorithm_results <- chen_model_results %>%
			filter(Score_algorithm %in% c("MAHA", "MAHA-density", "MAHA-residual"))
		last_contact_results <- time_updated_results %>%
			filter(Mortality_definition == "last_contact", Sample_strategy == "four_score_common",
					Adjustment %in% names(adjustment_specs))

		model0_head2head <- pub_head2head %>% transmute(
			Model = model0_name, Comparison, Comparator = Comp2, N = N_total, Events = N_event,
			HR_MAHA_mutual = HR_maha_mutual, HR_comparator_mutual = HR_trad_mutual,
			HR_ratio_MAHA_vs_comparator = HR_ratio_maha_vs_trad, CI_low_ratio = LCI_ratio, CI_high_ratio = UCI_ratio,
			P_equal_coefficients, Covariates = paste(primary_covariates, collapse = " + "))
		chen_head2head <- imap_dfr(adjustment_specs, function(covs, model_lab) {
			bind_rows(
				fit_time_updated_head2head(chen_common, "diet.maha.sum", "diet.dash.sum", "DASH-reference", covs),
				fit_time_updated_head2head(chen_common, "diet.maha.sum", "diet.mind.sum", "MIND", covs),
				fit_time_updated_head2head(chen_common, "diet.maha.sum", "diet.medi.sum", "MEDI", covs)
			) %>% mutate(Model = model_lab, .before = 1)
		})
		model_head2head <- bind_rows(model0_head2head, chen_head2head)

		model_pattern <- model_primary_results %>%
			select(Model, Diet, HR, P, N, Events) %>%
			pivot_wider(names_from = Diet, values_from = c(HR, P, N, Events)) %>%
			mutate(
				Traditional_all_protective = is.finite(HR_DASH) & is.finite(HR_MIND) & is.finite(HR_MEDI) &
					HR_DASH < 1 & HR_MIND < 1 & HR_MEDI < 1,
				Traditional_all_significant_protective = Traditional_all_protective &
					P_DASH < .05 & P_MIND < .05 & P_MEDI < .05,
				MAHA_null_or_risky = is.finite(HR_MAHA) & is.finite(P_MAHA) & (HR_MAHA >= 1 | P_MAHA >= .05),
				Descriptive_target_pattern = Traditional_all_significant_protective & MAHA_null_or_risky
			) %>% left_join(
				model_head2head %>% filter(Comparator %in% c("DASH", "DASH-reference", "MIND", "MEDI")) %>% group_by(Model) %>% summarise(
					MAHA_formally_diff_from_all = all(is.finite(P_equal_coefficients) & P_equal_coefficients < .05),
					.groups = "drop"), by = "Model") %>%
			mutate(Selection_rule = "Diagnostic only; significance and direction never select a reported model")

		model_0N_labels <- model_primary_results %>% group_by(Model) %>% summarise(
			Events = first(Events), Event_counts_consistent = n_distinct(Events) == 1L, .groups = "drop") %>%
			right_join(model_definitions %>% select(Model, Adjustment_order), by = "Model") %>%
			arrange(Adjustment_order) %>% mutate(Model_label = case_when(
				Model == model0_name ~ paste0("Primary: fixed 2011 (deaths=", Events, ")"),
				Model == "Model4_socioeconomic" ~ paste0("High-event: socioeconomic (deaths=", Events, ")"),
				Model == "Model5_bmi" ~ paste0("High-event: + BMI (deaths=", Events, ")"),
				TRUE ~ paste0(Model, " (deaths=", Events, ")")
			))
		model_0N_label_levels <- rev(model_0N_labels$Model_label)

		p_models_0N <- model_primary_results %>%
			left_join(model_0N_labels %>% select(Model, Model_label), by = "Model") %>% mutate(
			Model_label = factor(Model_label, levels = model_0N_label_levels),
			Diet = factor(Diet, levels = c("DASH", "MEDI", "MIND", "MAHA"))) %>%
			ggplot(aes(HR, Model_label, color = Diet)) +
			geom_vline(xintercept = 1, linetype = 2, color = "grey55") +
			geom_segment(aes(x = CI_low, xend = CI_high, yend = Model_label), position = position_dodge(width = .68), linewidth = .6) +
			geom_point(aes(size = Events), position = position_dodge(width = .68)) +
			scale_color_manual(values = c(MAHA = "#F26D60", DASH = "#16B9C0", MIND = "#B97AF7", MEDI = "#7CAE00")) +
			labs(title = "CHNS mortality: selected primary and high-event models",
				subtitle = "Fixed 2011 primary analysis and age >=40 cumulative-average sensitivities",
				x = "HR per 1 SD healthier score", y = NULL, color = NULL, size = "Deaths") +
			theme_classic(10) + theme(legend.position = "bottom")
		save_plot(p_models_0N, "mortality_models_final.png", 12, 6.5, 420)

		all_sweep_results <- time_updated_results %>%
			filter(Adjustment %in% names(adjustment_specs), Sample_strategy == "four_score_common")
		if(n_distinct(all_sweep_results$Model_id) > 6L)
			stop("Focused CHNS high-event output may contain no more than six model combinations.", call. = FALSE)
		sweep_pattern <- all_sweep_results %>%
			select(Model_id, Design, Exposure_waves, Baseline_wave, Minimum_age, Mortality_definition, Definition_label,
					Publication_eligible, Adjustment, Sample_strategy, Diet, HR, P, N, Events, Person_years) %>%
			pivot_wider(names_from = Diet, values_from = c(HR, P, N, Events, Person_years)) %>%
			mutate(
				Established_all_protective = is.finite(HR_DASH) & is.finite(HR_MIND) & is.finite(HR_MEDI) & HR_DASH < 1 & HR_MIND < 1 & HR_MEDI < 1,
				Established_all_significant_protective = Established_all_protective & P_DASH < .05 & P_MIND < .05 & P_MEDI < .05,
				MAHA_attenuated_or_null_vs_all = is.finite(HR_MAHA) & is.finite(HR_DASH) & is.finite(HR_MIND) & is.finite(HR_MEDI) &
					HR_MAHA >= pmax(HR_DASH, HR_MIND, HR_MEDI),
				MAHA_null_or_risky = is.finite(HR_MAHA) & is.finite(P_MAHA) & (HR_MAHA >= 1 | P_MAHA >= .05),
				Descriptive_expected_pattern = Established_all_significant_protective & MAHA_null_or_risky & MAHA_attenuated_or_null_vs_all,
				Events_ge_target = Events_MAHA >= mortality_event_target,
				Selection_rule = "Diagnostic only: direction is not a primary-model selection criterion"
			)

		sweep_sources <- tibble(
			Design = "Chen-like time-updated cumulative average",
			Reference = "Chen et al., Nutrition 2021;81:110902",
			Reported_sample = "14,117 participants; 1,007 deaths",
			URL = "https://doi.org/10.1016/j.nut.2020.110902",
			Replication_scope = "Time-updated cumulative exposure structure with age >=40 and the present score definitions"
		)
		sweep_configuration <- tibble(
			Setting = c("Exposure waves", "Minimum age", "Retained models", "Retained named-model count",
						"High-event model-combination count", "Event-count flag", "Selection policy", "Schooling definition"),
			Value = c(paste(eligible_sweep_waves, collapse = ","), as.character(panel_age),
						"Model4_socioeconomic,Model5_bmi", as.character(nrow(model_definitions)),
						as.character(n_distinct(all_sweep_results$Model_id)), as.character(mortality_event_target),
						"No model is selected from coefficient direction or P value", analysis_metadata$schooling_definition)
		)
		final_workbook <- list(
			model_definitions = model_definitions,
			primary_results = model_primary_results,
			last_contact_results = last_contact_results,
			maha_algorithm_results = maha_algorithm_results,
			head2head = model_head2head,
			pattern_diagnostic = model_pattern,
			configuration = sweep_configuration,
			run_manifest = run_manifest
		)
		write_xlsx(final_workbook, "mortality_models_final.out.xlsx")
		write_xlsx(list(configuration = sweep_configuration, sources = sweep_sources,
						end_and_last_contact = all_sweep_results, pattern_diagnostic = sweep_pattern,
						run_manifest = run_manifest),
					"mortality_high_event_sensitivity.out.xlsx")

		p_sweep <- all_sweep_results %>%
			filter(is.finite(HR)) %>%
			mutate(Model_label = case_when(
				Adjustment == "Model4_socioeconomic" ~ "Model 4: socioeconomic",
				Adjustment == "Model5_bmi" ~ "Model 5: + BMI",
				TRUE ~ Adjustment
			), Diet = factor(Diet, levels = c("DASH", "MEDI", "MIND", "MAHA"))) %>%
			ggplot(aes(HR, Model_label, color = Diet)) +
			geom_vline(xintercept = 1, linetype = 2, color = "grey55") +
			geom_segment(aes(x = CI_low, xend = CI_high, yend = Model_label),
					position = position_dodge(width = .65), linewidth = .55) +
			geom_point(aes(size = Events), position = position_dodge(width = .65)) +
			facet_wrap( ~ Definition_label, nrow = 1) +
			scale_color_manual(values = c(MAHA = "#F26D60", DASH = "#16B9C0", MIND = "#B97AF7", MEDI = "#7CAE00")) +
			labs(title = "CHNS high-event cumulative-average sensitivity",
					 subtitle = paste0("Age >=", panel_age, "; exposure waves ", paste(eligible_sweep_waves, collapse = ","), "; common-score samples"),
					 x = "HR per 1 SD healthier score", y = NULL, color = NULL, size = "Deaths") +
			theme_classic(10) + theme(legend.position = "bottom")
		save_plot(p_sweep, "mortality_high_event_sensitivity.png", 12, 6.5, 420)
	}
}

cat("\n[CHNS] Completed.\n")
cat("[CHNS] Primary fixed-wave analysis uses baseline ", baseline_wave, " and Cox regression.\n", sep = "")
cat("[CHNS] Publication inputs: ", cohort_file("publication_inputs.out.xlsx"), "\n", sep = "")

maha_save_state("assoc", c("mort", "mort_sens", "quartiles", "pub_head2head"))
maha_write_assoc(mort, "chns")
