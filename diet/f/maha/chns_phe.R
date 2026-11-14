# MAHA manuscript engine, reorganized from jielab/scripts 971a61b44b84500ab46e2f5592478505a0eae50e.
maha_load_state("score")
if (!identical(Sys.getenv("CHNS_MORTALITY_SWEEP_CHILD"), "1")) dat <- maha_join_scores(dat, "IDIND")
dat <- chns_refresh_scales(dat, all_score_vars)
dat <- dat %>% mutate(
	htn_measured = ifelse(is.finite(sbp) | is.finite(dbp), as.integer((sbp >= 140) | (dbp >= 90)), NA_integer_),
	hypertension = ifelse(htn_self == 1 | htn_measured == 1, 1L, ifelse(htn_self == 0 & htn_measured == 0, 0L, NA_integer_)),
	obesity = ifelse(is.finite(bmi), as.integer(bmi >= 28), NA_integer_),
	bmi_fu = ifelse(is.finite(weight_kg_fu) & is.finite(height_cm_fu) & height_cm_fu > 100, weight_kg_fu / (height_cm_fu / 100) ^ 2, NA_real_),
	htn_measured_fu = ifelse(is.finite(sbp_fu) | is.finite(dbp_fu), as.integer((sbp_fu >= 140) | (dbp_fu >= 90)), NA_integer_),
	hypertension_fu = ifelse(htn_self_fu == 1 | htn_measured_fu == 1, 1L, ifelse(htn_self_fu == 0 & htn_measured_fu == 0, 0L, NA_integer_)),
	obesity_fu = ifelse(is.finite(bmi_fu), as.integer(bmi_fu >= 28), NA_integer_),
	inc_hypertension = ifelse(hypertension == 0 & hypertension_fu %in% 0:1, as.integer(hypertension_fu == 1), NA_integer_),
	inc_obesity = ifelse(obesity == 0 & obesity_fu %in% 0:1, as.integer(obesity_fu == 1), NA_integer_),
	inc_t2dm = ifelse(t2dm_self == 0 & t2dm_fu %in% 0:1, as.integer(t2dm_fu == 1), NA_integer_),
	inc_mi = ifelse(mi_self == 0 & mi_fu %in% 0:1, as.integer(mi_fu == 1), NA_integer_),
	inc_stroke = ifelse(stroke_self == 0 & stroke_fu %in% 0:1, as.integer(stroke_fu == 1), NA_integer_),
	inc_cancer = ifelse(cancer_self == 0 & cancer_fu %in% 0:1, as.integer(cancer_fu == 1), NA_integer_),
	delta_bmi = ifelse(is.finite(bmi) & is.finite(bmi_fu), bmi_fu - bmi, NA_real_),
	delta_sbp = ifelse(is.finite(sbp) & is.finite(sbp_fu), sbp_fu - sbp, NA_real_),
	delta_dbp = ifelse(is.finite(dbp) & is.finite(dbp_fu), dbp_fu - dbp, NA_real_),
	observed_at_followup = replace_na(as.integer(observed_at_followup), 0L),
	observed_after_baseline = replace_na(as.integer(observed_after_baseline), 0L),
	death_by_year = is.finite(DOD_Y) & DOD_Y > baseline_wave & DOD_Y <= followup_end,
	death_by_report = is.finite(DOD_RPT) & DOD_RPT > baseline_wave & DOD_RPT <= followup_end,
	death_year = case_when(
		death_by_year & death_by_report ~ pmin(DOD_Y, DOD_RPT),
		death_by_year ~ DOD_Y,
		death_by_report ~ DOD_RPT,
		TRUE ~ NA_real_
	),
	death_reported = as.integer(death_by_year | death_by_report),
	followup_known = as.integer(observed_at_followup == 1 | death_reported == 1),
	followup_known_last_contact = as.integer(observed_after_baseline == 1 | death_reported == 1),
	death = ifelse(followup_known == 1, death_reported, NA_integer_)
)
dat <- add_primary_survival_fields(dat)
maha_save_required_state("phe", c("chns_assoc.R", "chns_bayes.R"))
maha_write_phenotypes(dat, "IDIND")
