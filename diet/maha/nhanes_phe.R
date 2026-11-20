# MAHA manuscript engine, reorganized from jielab/scripts 971a61b44b84500ab46e2f5592478505a0eae50e.
maha_load_state("score")
dat <- maha_join_scores(dat_observed, "SEQN")
dat <- nhanes_refresh_scales(dat)
if (use_mortality) {
	mort <- read_mortality()
	if (!is.null(mort)) dat <- left_join(dat, mort, by = "SEQN")
}

dat <- dat %>%
	mutate(
		nonhdl_mgdl = tc_mgdl - hdl_mgdl,
		egfr = egfr_2021(creat_mgdl, age, female),
		t2dm = as.integer(diabetes_q == 1 | hba1c >= 6.5 | glu_mgdl >= 126),
		cad = as.integer(cad_q == 1 | mi_q == 1),
		mi = as.integer(mi_q == 1),
		stroke = as.integer(stroke_q == 1),
		heart_failure = as.integer(hf_q == 1),
		ckd = as.integer(ckd_q == 1 | egfr < 60 | uacr_mgg >= 30),
		hypertension = as.integer(htn_q == 1 | sbp >= 130 | dbp >= 80),
		copd = as.integer(copd_q == 1),
		asthma = as.integer(asthma_q == 1),
		depression = as.integer(phq9 >= 10),
		obesity = as.integer(bmi >= 30),
		mort_eligible = if ("eligstat" %in% names(.)) as.integer(eligstat == 1) else NA_integer_,
		death = if ("mortstat" %in% names(.)) ifelse(eligstat == 1, as.integer(mortstat == 1), NA_integer_) else NA_integer_,
		death_time_y = if ("permth_exm" %in% names(.)) as.numeric(permth_exm) / 12 else NA_real_,
		heart_death = if ("ucod_leading" %in% names(.)) ifelse(eligstat == 1, as.integer(mortstat == 1 & ucod_leading == 1), NA_integer_) else NA_integer_,
		cvd_death = if ("ucod_leading" %in% names(.)) ifelse(eligstat == 1, as.integer(mortstat == 1 & ucod_leading %in% c(1, 5)), NA_integer_) else NA_integer_,
		cancer_death = if ("ucod_leading" %in% names(.)) ifelse(eligstat == 1, as.integer(mortstat == 1 & ucod_leading == 2), NA_integer_) else NA_integer_,
		diabetes_death = if ("ucod_leading" %in% names(.)) ifelse(eligstat == 1, as.integer(mortstat == 1 & ucod_leading == 7), NA_integer_) else NA_integer_
	)

dat_observed <- dat
dat1_obs <- prep_analysis_dataset(dat_observed, main_cycles, require_mortality = FALSE)
dat_mort_obs <- prep_analysis_dataset(dat_observed, mort_cycles, require_mortality = TRUE)

dat1 <- dat1_obs
dat_mort1 <- dat_mort_obs

diet.sum.cols <- paste0("diet.", names(diet.lst), ".sum")
diet.sum.cols <- diet.sum.cols[diet.sum.cols %in% names(dat1_obs)]
diet.inc.pts <- paste0("diet.", diet.inc, ".pts")
diet.inc.s100 <- paste0("diet.", diet.inc, ".s100")
diet.inc.q5 <- paste0("diet.", diet.inc, ".q5")
diet.inc.3c <- paste0("diet.", diet.inc, ".3c")

message("N after QC, observed main analysis: ", nrow(dat1_obs))
message("N after QC, mortality analysis: ", nrow(dat_mort_obs))
message("Main cycles retained: ", paste(sort(unique(dat1_obs$cycle)), collapse = ", "))
message("Mortality cycles retained: ", paste(sort(unique(dat_mort_obs$cycle)), collapse = ", "))

covs.base <- c("age", "female", "race", "edu", "pir", "cycle", "kcal")
covs.add <- c("bmi", "smoke_ever", "smoke_now")
covs <- c(covs.base, covs.add)

covs.mort <- c("age", "female", "race", "edu", "pir", "marital", "health_insurance", "cycle", "kcal", "bmi", "smoke_ever", "smoke_now", "hypertension", "cad", "stroke", "heart_failure", "t2dm", "ckd")
Y.cross <- c("cad", "stroke", "heart_failure", "t2dm", "ckd")
Y.mort.plot <- c("death", "heart_death", "cvd_death", "cancer_death", "diabetes_death")


maha_save_state("phe", c("dat_observed", "dat1_obs", "dat_mort_obs", "dat1", "dat_mort1", "food_vars", "nutrient_vars", "diet.sum.cols", "diet.inc.pts", "diet.inc.s100", "diet.inc.q5", "diet.inc.3c", "covs.base", "covs.add", "covs", "covs.mort", "Y.cross", "Y.mort.plot"))
maha_write_phenotypes(dat1_obs, "SEQN", mortality = dat_mort_obs)
