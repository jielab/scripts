# MAHA manuscript engine, reorganized from jielab/scripts 971a61b44b84500ab46e2f5592478505a0eae50e.
nhanes_check_inputs(score_only = TRUE)

dat <- pmap_dfr(cycle_info, read_cycle)

if (nrow(dat) == 0) stop("No NHANES data were read. Check raw_dir: ", raw_dir)

for (v in c("cad_q", "copd_q")) dat[[v]][is.infinite(dat[[v]])] <- NA_real_

food_vars <- c("fruit", "berry", "vegetable", "green_leafy", "allium", "legumes", "nuts", "dairy", "lowfat_dairy",
	"whole_grain", "refined_grain", "fish", "poultry", "red_processed_meat", "processed_meat_proxy", "sweets_pastries", "ssb", "fried_fast", "coffee_tea")
nutrient_vars <- c("kcal", "protein_g", "carb_g", "sugar_g", "fiber_g", "fat_g", "sfat_g", "mufa_g", "pufa_g", "chol_mg", "sodium_mg", "alcohol_g")
for (v in food_vars) if (!v %in% names(dat)) dat[[v]] <- NA_real_
for (v in nutrient_vars) if (!v %in% names(dat)) dat[[v]] <- NA_real_

step_header("Diet score construction: observed repeated-recall mean")
dat_observed <- construct_diet_scores(dat, score_source = "observed_repeated_recall_mean")


maha_save_state("score", c("dat_observed", "food_vars", "nutrient_vars"))
maha_export_scores(dat_observed, "SEQN")
