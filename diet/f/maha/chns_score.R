# MAHA manuscript engine, reorganized from jielab/scripts 971a61b44b84500ab46e2f5592478505a0eae50e.
# Input-file QC and source provenance
step_header("CHNS input-file QC: single primary cohort")

req_path <- function(env_name, default_path) {
	x <- Sys.getenv(env_name, unset = default_path)
	normalizePath(x, winslash = "/", mustWork = FALSE)
}
opt_path <- function(env_name, default_path) {
	x <- Sys.getenv(env_name, unset = default_path)
	if (file.exists(x)) normalizePath(x, winslash = "/", mustWork = FALSE) else NA_character_
}

src <- list(
	NUTR3 = c(path = req_path("CHNS_NUTR3", file.path(raw_dir, "diet", "nutr3_00.sas7bdat")), source = "CHNS_public_use_plus_user_download"),
	C12DIET = c(path = req_path("CHNS_C12DIET", file.path(raw_dir, "diet", "c12diet.sas7bdat")), source = "CHNS_public_use"),
	MASTER = c(path = req_path("CHNS_MASTER", file.path(id2019_dir, "mast_pub_12.sas7bdat")), source = "CHNS_2019_archive"),
	ROSTER = c(path = req_path("CHNS_ROSTER", file.path(id2019_dir, "rst_12.sas7bdat")), source = "CHNS_2019_archive"),
	PEXAM = c(path = req_path("CHNS_PEXAM", pexam_public_default), source = "CHNS_public_use_extra_copy_user_approved"),
	PACT = c(path = req_path("CHNS_PACT", file.path(raw_dir, "individual", "pact_12.sas7bdat")), source = "CHNS_public_use"),
	EDUC = c(path = opt_path("CHNS_EDUC", file.path(raw_dir, "individual", "educ_12.sas7bdat")), source = "CHNS_public_use_optional"),
	HHINC = c(path = opt_path("CHNS_HHINC", file.path(raw_dir, "constructed", "hhinc_10.sas7bdat")), source = "CHNS_public_use_optional"),
	JOBS = c(path = opt_path("CHNS_JOBS", file.path(raw_dir, "individual", "jobs_12.sas7bdat")), source = "CHNS_public_use_optional"),
	SURVEY = c(path = opt_path("CHNS_SURVEY", file.path(id2019_dir, "surveys_pub_12.sas7bdat")), source = "CHNS_2019_archive_optional"),
	HEALTH = c(path = opt_path("CHNS_HEALTH", file.path(raw_dir, "individual", "hlth_12.sas7bdat")), source = "CHNS_public_use_optional")
)

food_map_file <- Sys.getenv("CHNS_FOOD_MAP", unset = food_map_default)
nutrient_map_file <- Sys.getenv("CHNS_NUTRIENT_MAP", unset = nutrient_map_default)
if (!file.exists(food_map_file)) food_map_file <- NA_character_
if (!file.exists(nutrient_map_file)) nutrient_map_file <- NA_character_

cat("[CHNS INPUT CHECK] project directory: ", project_dir, "\n", sep = "")
cat("[CHNS INPUT CHECK] raw directory: ", raw_dir, "\n", sep = "")
required_names <- c("NUTR3", "C12DIET", "MASTER", "ROSTER", "PEXAM", "PACT")
required_ok <- logical()
for (nm in names(src)) {
	ok <- print_source(nm, src[[nm]], required = nm %in% required_names)
	if (nm %in% required_names) required_ok <- c(required_ok, ok)
}
cat("  FOOD_MAP     ", if(!is.na(food_map_file))"OK      classification : " else "MISSING required       : ",
		ifelse(is.na(food_map_file), food_map_default, food_map_file), ifelse(is.na(food_map_file), " 文件不存在", ""), "\n", sep = "")
cat("  NUTRIENT_MAP ", if(!is.na(nutrient_map_file))"OK      code-level map       : " else if(require_nutrient_map == 1L)"MISSING required       : " else "OPTIONAL missing       : ",
		ifelse(is.na(nutrient_map_file), nutrient_map_default, nutrient_map_file), "\n", sep = "")

if (!all(required_ok) || is.na(food_map_file)) {
	stop("CHNS required raw input/classification map missing. See ", log_file, call. = FALSE)
}
if (require_nutrient_map == 1L && is.na(nutrient_map_file)) {
	stop("CHNS_NUTRIENT_MAP is required in the publication pipeline but was not found: ", nutrient_map_default, call. = FALSE)
}

# Source provenance is maintained beside the cohort data, not scoring outputs.
# Paths recorded on another machine are informational; scoring uses src above.
manifest_file <- file.path(project_dir, "chns.manifest.tsv")
if (file.exists(manifest_file)) {
    cat("[CHNS MANIFEST] ", manifest_file, " (provenance only)\n", sep = "")
} else {
    cat("[CHNS MANIFEST MISSING] ", manifest_file,
        "; sync the cohort provenance file to the data directory. Scoring uses configured raw inputs.\n", sep = "")
}


# 🚩 2) Read modules and resolve waves
step_header("CHNS variable and wave QC")

# nutr3_00 is large: only keep fields required here.
nutr3 <- read_any(src$NUTR3[["path"]], c("IDIND", "COMMID", "WAVE", "FOODCODE", "VD", "T1", "T2", "V39", "V39C", "V41", "V42"))
c12 <- read_any(src$C12DIET[["path"]])
master <- read_any(src$MASTER[["path"]])
roster <- read_any(src$ROSTER[["path"]])
pexam <- read_any(src$PEXAM[["path"]])
pexam_long <- pexam
pact <- read_any(src$PACT[["path"]])
educ <- if (!is.na(src$EDUC[["path"]])) read_any(src$EDUC[["path"]]) else NULL
hhinc <- if (!is.na(src$HHINC[["path"]])) read_any(src$HHINC[["path"]]) else NULL
jobs_module <- if (!is.na(src$JOBS[["path"]])) read_any(src$JOBS[["path"]]) else NULL
survey <- if (!is.na(src$SURVEY[["path"]])) read_any(src$SURVEY[["path"]]) else NULL
health <- if (!is.na(src$HEALTH[["path"]])) read_any(src$HEALTH[["path"]]) else NULL


# Official longitudinal pact_12 is a large person-year raw module. Reject small/derived substitutes.
if (nrow(pact) < 100000 || ncol(pact) < 250) {
	stop("PACT raw-data sanity check failed: expected official longitudinal pact_12 (~125k rows, 276 variables), got ",
			 nrow(pact), " rows and ", ncol(pact), " variables from ", src$PACT[["path"]],
			 ". Do not use pact12.sav or pact12 - 1/2/3/4.sav author-analysis files.", call. = FALSE)
}
cat("[CHNS RAW SANITY] PACT rows=", nrow(pact), "; variables=", ncol(pact), " [PASS]\n", sep = "")

idn <- pick(nutr3, "IDIND"); wn <- pick(nutr3, "WAVE"); fn <- pick(nutr3, "FOODCODE"); an <- pick(nutr3, "V39"); dn <- pick(nutr3, "VD")
proc_n <- pick(nutr3, "V39C", FALSE); prep_n <- pick(nutr3, "V42", FALSE); loc_n <- pick(nutr3, "V41", FALSE)
comm_n <- pick(nutr3, "COMMID", FALSE); prov_n <- pick(nutr3, "T1", FALSE); urban_n <- pick(nutr3, "T2", FALSE)

idc <- pick(c12, "IDIND"); wc <- pick(c12, "WAVE")
idm <- pick(master, "IDIND")
idr <- pick(roster, "IDIND"); wr <- pick(roster, "WAVE")
idp <- pick(pexam, "IDIND"); wp <- pick(pexam, c("WAVE", "YEAR"))
idpa <- pick(pact, "IDIND"); wpa <- pick(pact, "WAVE")

nutr_waves <- valid_years(nutr3[[wn]])
c12_waves <- valid_years(c12[[wc]])
pex_waves <- valid_years(pexam[[wp]])
idpl <- idp; wpl <- wp
pex_long_waves <- valid_years(pexam_long[[wpl]])
pact_waves <- valid_years(pact[[wpa]])
roster_waves <- valid_years(roster[[wr]])
common_baseline <- Reduce(intersect, list(nutr_waves, c12_waves, pex_waves, pact_waves))

cat("  NUTR3 waves : ", paste(nutr_waves, collapse = ", "), "\n", sep = "")
cat("  C12DIET waves: ", paste(c12_waves, collapse = ", "), "\n", sep = "")
cat("  PEXAMPUB waves: ", paste(pex_waves, collapse = ", "), "\n", sep = "")
cat("  PACT waves  : ", paste(pact_waves, collapse = ", "), "\n", sep = "")
cat("  ROSTER waves: ", paste(roster_waves, collapse = ", "), "\n", sep = "")
cat("  Common baseline waves: ", paste(common_baseline, collapse = ", "), "\n", sep = "")

baseline_wave <- parse_year_request(baseline_req, "CHNS_BASELINE_WAVE")
followup_end <- parse_year_request(followup_req, "CHNS_FOLLOWUP_END")

if (is.na(baseline_wave)) {
	eligible <- common_baseline[vapply(common_baseline, function(w) any(roster_waves >= w + min_followup_years), logical(1))]
	if (!length(eligible)) eligible <- common_baseline[vapply(common_baseline, function(w) any(roster_waves > w), logical(1))]
	if (!length(eligible)) stop("No common dietary/exam/preference wave has later roster follow-up.", call. = FALSE)
	baseline_wave <- max(eligible)
}
if (!baseline_wave %in% common_baseline) {
	stop("Requested baseline ", baseline_wave, " is unavailable jointly in NUTR3/C12DIET/PEXAM/PACT.", call. = FALSE)
}
if (is.na(followup_end)) {
	later <- roster_waves[roster_waves > baseline_wave]
	if (!length(later)) stop("No roster follow-up later than baseline ", baseline_wave, call. = FALSE)
	followup_end <- max(later)
}
if (!followup_end %in% roster_waves) stop("Requested follow-up ", followup_end, " not present in ROSTER.", call. = FALSE)

cat("[CHNS WAVE RESOLUTION] baseline=", baseline_wave, "; follow-up=", followup_end,
		"; gap=", followup_end - baseline_wave, " years\n", sep = "")


# 🚩 3) CFCT6 classification map + optional code-level nutrient map
step_header("CHNS FOODCODE classification and map QC")

class_map <- read_any(food_map_file)
required_flags <- c("CEREAL", "POTATO_STARCH", "LEGUME", "VEGETABLE", "LEAF_STEM_FLOWER_VEG_PROXY",
			"FRUIT", "BERRY", "NUT_SEED", "RED_MEAT", "POULTRY", "DAIRY", "CHEESE", "EGG",
			"FISH_SHELLFISH", "SWEET_PASTRY", "FAST_FOOD", "CONVENIENCE_FOOD", "SNACK_FOOD",
			"BEVERAGE", "CARBONATED_DRINK", "ALCOHOL_BEVERAGE", "SUGAR_CANDY_PRESERVE",
			"ANIMAL_FAT", "PLANT_OIL", "CONDIMENT", "SALT_MSG_OTHER", "UPF_CATEGORY_PROXY",
			"SSB_PROXY", "SWEET_SNACK_PROXY")
required_class <- c("CFCT6_PREFIX3", "CLASS2_TEXT", "SUBCLASS1_TEXT", "CLASS_NAME_CN", "SUBCLASS_NAME_CN", required_flags)
miss_class <- setdiff(required_class, names(class_map))
if(length(miss_class)) stop("Classification FOOD_MAP schema incomplete. Missing: ", paste(miss_class, collapse = ", "), call. = FALSE)
class_map <- class_map %>%
	mutate(CODE_PREFIX3 = sprintf("%03d", as.integer(num(CFCT6_PREFIX3))),
			CLASS2 = as.integer(num(CLASS2_TEXT)), SUBCLASS1 = as.integer(num(SUBCLASS1_TEXT))) %>%
	distinct(CODE_PREFIX3, .keep_all = TRUE)
for(v in required_flags) class_map[[v]] <- num(class_map[[v]])
required_join <- c("CODE_PREFIX3", "CLASS2", "SUBCLASS1", "CLASS_NAME_CN", "SUBCLASS_NAME_CN", required_flags)

food0 <- nutr3 %>%
	transmute(
		IDIND = num(.data[[idn]]), WAVE = num(.data[[wn]]), FOODCODE = num(.data[[fn]]),
		GRAMS = num(.data[[an]]), DAY = num(.data[[dn]]),
		COMMID = colnum(nutr3, comm_n), PROVINCE = colnum(nutr3, prov_n), URBAN = colnum(nutr3, urban_n),
		PROCESSED_CODE = colnum(nutr3, proc_n), PREP_METHOD = colnum(nutr3, prep_n), MEAL_LOCATION = colnum(nutr3, loc_n)
	) %>%
	filter(WAVE == baseline_wave, is.finite(IDIND), is.finite(FOODCODE), is.finite(GRAMS), GRAMS >= 0) %>%
	mutate(FOODCODE6 = sprintf("%06d", as.integer(FOODCODE)), CODE_PREFIX3 = substr(FOODCODE6, 1, 3)) %>%
	left_join(class_map %>% select(all_of(required_join)), by = "CODE_PREFIX3")
if (!nrow(food0)) stop("No nutr3_00 records at baseline wave ", baseline_wave, call. = FALSE)

class_cov <- sum(food0$GRAMS[!is.na(food0$CLASS2)], na.rm = TRUE) / sum(food0$GRAMS, na.rm = TRUE)
cat(sprintf("[CFCT6 classification] gram-weighted coverage %.1f%%\n", 100 * class_cov))
if(class_cov < .95) warning(sprintf("CFCT6 prefix classification covers only %.1f%% of baseline food grams; inspect unmatched FOODCODEs.", 100 * class_cov))

# Native CHNS processing/preparation indicators and class-derived proxies.
food0 <- food0 %>% mutate(
	OFFICIAL_PROCESSED = ifelse(is.finite(PROCESSED_CODE), as.integer(PROCESSED_CODE == 3), NA_integer_),
	RESTAURANT_MADE = ifelse(is.finite(PROCESSED_CODE), as.integer(PROCESSED_CODE == 2), NA_integer_),
	DEEPFRIED = ifelse(is.finite(PREP_METHOD), as.integer(PREP_METHOD == 3), NA_integer_)
)

# Code-level audit for transparent manual completion later.
partial_map <- food0 %>% group_by(FOODCODE, FOODCODE6, CODE_PREFIX3) %>% summarise(
	CLASS2 = first(na.omit(CLASS2), default = NA),
	SUBCLASS1 = first(na.omit(SUBCLASS1), default = NA),
	CLASS_NAME_CN = first(na.omit(CLASS_NAME_CN), default = NA_character_),
	SUBCLASS_NAME_CN = first(na.omit(SUBCLASS_NAME_CN), default = NA_character_),
	grams = sum(GRAMS, na.rm = TRUE), records = n(),
	V39C_processed_fraction = ifelse(any(is.finite(PROCESSED_CODE)), mean(PROCESSED_CODE == 3, na.rm = TRUE), NA_real_),
	V42_deepfried_fraction = ifelse(any(is.finite(PREP_METHOD)), mean(PREP_METHOD == 3, na.rm = TRUE), NA_real_),
	.groups = "drop"
) %>% arrange(desc(grams))
readr::write_tsv(partial_map, cohort_file("food_map.code_audit.tsv"))

# Prefer CFCT 2002/2004 numeric FOODCODE rows; evaluate coverage by field.
full_map <- NULL
code_map_available <- FALSE
full_map_mode <- FALSE
enhanced_map_mode <- FALSE
full_map_coverage <- NA_real_
map_field_coverage <- tibble(field = character(), gram_coverage = double())
map_field_ok <- setNames(logical(), character())
required_code_fields <- c("WHOLE_GRAIN", "REFINED_GRAIN", "LOWFAT_DAIRY", "GREEN_LEAFY",
			"PROCESSED_MEAT", "SSB", "SODIUM_MG_100G", "SFAT_G_100G", "MUFA_G_100G",
			"PUFA_G_100G", "ALCOHOL_G_100G")
fine_group_fields <- c("WHOLE_GRAIN", "REFINED_GRAIN", "LOWFAT_DAIRY", "GREEN_LEAFY", "PROCESSED_MEAT", "SSB")
fat_fields <- c("SFAT_G_100G", "MUFA_G_100G", "PUFA_G_100G")

if(!is.na(nutrient_map_file)) {
	tmp <- read_any(nutrient_map_file)
	miss <- setdiff(required_code_fields, names(tmp))
	if(length(miss) || (!"FOODCODE6" %in% names(tmp) && !"FOODCODE" %in% names(tmp))) {
		warning("NUTRIENT_MAP schema incomplete. Missing: ",
					paste(c(miss, if(!"FOODCODE6" %in% names(tmp) && !"FOODCODE" %in% names(tmp))"FOODCODE or FOODCODE6" else character()), collapse = ", "))
	} else {
		if("FOODCODE6" %in% names(tmp)) {
			tmp <- tmp %>% mutate(FOODCODE6 = sprintf("%06d", as.integer(num(FOODCODE6))))
		} else {
			tmp <- tmp %>% mutate(FOODCODE6 = sprintf("%06d", as.integer(num(FOODCODE))))
		}
		tmp <- tmp %>% filter(grepl("^[0-9]{6}$", FOODCODE6)) %>% distinct(FOODCODE6, .keep_all = TRUE)
		for(v in required_code_fields) tmp[[v]] <- num(tmp[[v]])
		gram_by_code <- food0 %>% group_by(FOODCODE6) %>% summarise(grams = sum(GRAMS, na.rm = TRUE), .groups = "drop")
		aud <- gram_by_code %>% left_join(tmp, by = "FOODCODE6")
		den <- sum(aud$grams, na.rm = TRUE)
		if(is.finite(den) && den > 0) {
			full_map_coverage <- sum(aud$grams[!is.na(aud$WHOLE_GRAIN)], na.rm = TRUE) / den
			map_field_coverage <- tibble(field = required_code_fields) %>% rowwise() %>%
				mutate(gram_coverage = sum(aud$grams[is.finite(num(aud[[field]]))], na.rm = TRUE) / den) %>% ungroup()
			map_field_ok <- setNames(map_field_coverage$gram_coverage >= full_map_min_cov, map_field_coverage$field)
			full_map <- tmp
			code_map_available <- TRUE
			fine_ok <- all(map_field_ok[fine_group_fields], na.rm = TRUE)
			sodium_ok <- isTRUE(map_field_ok[["SODIUM_MG_100G"]])
			fat_ok <- all(map_field_ok[fat_fields], na.rm = TRUE)
			alcohol_ok <- isTRUE(map_field_ok[["ALCOHOL_G_100G"]])
			enhanced_map_mode <- isTRUE(fine_ok)
			full_map_mode <- isTRUE(fine_ok && sodium_ok && fat_ok && alcohol_ok)
			cat("[CODE-LEVEL FOODCODE COVERAGE] ",
					paste0(map_field_coverage$field, "=", sprintf("%.1f%%", 100 * map_field_coverage$gram_coverage), collapse = "; "), "\n", sep = "")
			if(!enhanced_map_mode)
				warning(sprintf("Code-level fine-group coverage is below %.1f%% for at least one core field; proxy mode retained.", 100 * full_map_min_cov))
			else if(!full_map_mode)
				warning("Code-level fine groups are sufficiently covered, but one or more nutrient domains (sodium/fat/alcohol) are below threshold; foodcode-enhanced mode will use transparent proxy fallback for those domains.")
		}
	}
}
if(!exists("sodium_ok")) sodium_ok <- FALSE
if(!exists("fat_ok")) fat_ok <- FALSE
if(!exists("alcohol_ok")) alcohol_ok <- FALSE
if (require_enhanced_map == 1L && !enhanced_map_mode) {
	stop(sprintf("Code-level FOODCODE fine-group map coverage is insufficient for publication scoring (< %.0f%% in at least one required fine-group field). Inspect %s and %s.",
					100 * full_map_min_cov, food_map_file, nutrient_map_file), call. = FALSE)
}

analysis_mode <- if(full_map_mode) {
	"FOODCODE-complete mode: code-level map with all required nutrient domains above threshold"
} else if(enhanced_map_mode) {
	"FOODCODE-enhanced mode: code-level fine groups + field-specific nutrient/proxy fallback"
} else {
	"CFCT6 classification mode: MAHA-CN/DASH-CN/MIND-CN/MEDI-CN proxies"
}
cat("[CHNS MAP MODE] ", analysis_mode, "\n", sep = "")


# 🚩 4) Aggregate dietary intake
step_header("CHNS dietary harmonization")

food <- food0
if(code_map_available) food <- food %>% left_join(full_map %>% select(FOODCODE6, all_of(required_code_fields)), by = "FOODCODE6")
for(v in required_code_fields) if(!v %in% names(food)) food[[v]] <- NA_real_
for(v in required_code_fields) food[[v]] <- num(food[[v]])

food_day <- food %>% group_by(IDIND, DAY) %>% summarise(
	total_food_g = sum(GRAMS, na.rm = TRUE), n_food_types = n_distinct(FOODCODE),
	cereal = weighted_sum_known(GRAMS, CEREAL), potato_starch = weighted_sum_known(GRAMS, POTATO_STARCH),
	legumes = weighted_sum_known(GRAMS, LEGUME), vegetable = weighted_sum_known(GRAMS, VEGETABLE),
	leaf_stem_flower_veg = weighted_sum_known(GRAMS, LEAF_STEM_FLOWER_VEG_PROXY),
	fruit = weighted_sum_known(GRAMS, FRUIT), berry = weighted_sum_known(GRAMS, BERRY),
	nuts = weighted_sum_known(GRAMS, NUT_SEED), red_meat = weighted_sum_known(GRAMS, RED_MEAT),
	poultry = weighted_sum_known(GRAMS, POULTRY), dairy = weighted_sum_known(GRAMS, DAIRY),
	cheese = weighted_sum_known(GRAMS, CHEESE), eggs = weighted_sum_known(GRAMS, EGG),
	fish = weighted_sum_known(GRAMS, FISH_SHELLFISH),
	sweet_pastry_proxy = weighted_sum_known(GRAMS, SWEET_PASTRY),
	fast_food_proxy = weighted_sum_known(GRAMS, FAST_FOOD), convenience_food_proxy = weighted_sum_known(GRAMS, CONVENIENCE_FOOD),
	snack_food_proxy = weighted_sum_known(GRAMS, SNACK_FOOD), beverage_g = weighted_sum_known(GRAMS, BEVERAGE),
	carbonated_drink_g = weighted_sum_known(GRAMS, CARBONATED_DRINK), alcohol_beverage_g = weighted_sum_known(GRAMS, ALCOHOL_BEVERAGE),
	sugar_candy_preserve_g = weighted_sum_known(GRAMS, SUGAR_CANDY_PRESERVE),
	animal_fat_g = weighted_sum_known(GRAMS, ANIMAL_FAT), plant_oil_g = weighted_sum_known(GRAMS, PLANT_OIL),
	condiment_g = weighted_sum_known(GRAMS, CONDIMENT), salt_msg_proxy_g = weighted_sum_known(GRAMS, SALT_MSG_OTHER),
	upf_class_proxy_g = weighted_sum_known(GRAMS, UPF_CATEGORY_PROXY), ssb_proxy_g = weighted_sum_known(GRAMS, SSB_PROXY),
	sweet_snack_proxy_g = weighted_sum_known(GRAMS, SWEET_SNACK_PROXY),
	processed_food_g = weighted_sum_known(GRAMS, OFFICIAL_PROCESSED), restaurant_made_g = weighted_sum_known(GRAMS, RESTAURANT_MADE),
	deepfried_g = weighted_sum_known(GRAMS, DEEPFRIED),
	whole_grain = weighted_sum_known(GRAMS, WHOLE_GRAIN), refined_grain = weighted_sum_known(GRAMS, REFINED_GRAIN),
	lowfat_dairy = weighted_sum_known(GRAMS, LOWFAT_DAIRY), green_leafy = weighted_sum_known(GRAMS, GREEN_LEAFY),
	processed_meat = weighted_sum_known(GRAMS, PROCESSED_MEAT), ssb_code = weighted_sum_known(GRAMS, SSB),
	sodium_mg = sum_known(GRAMS / 100 * SODIUM_MG_100G), sfat_g = sum_known(GRAMS / 100 * SFAT_G_100G),
	mufa_g = sum_known(GRAMS / 100 * MUFA_G_100G), pufa_g = sum_known(GRAMS / 100 * PUFA_G_100G),
	alcohol_g_food = sum_known(GRAMS / 100 * ALCOHOL_G_100G), .groups = "drop"
)
food_person <- food_day %>% group_by(IDIND) %>% summarise(
	n_diet_days = n_distinct(DAY), across( - DAY, ~ mean(.x, na.rm = TRUE)), .groups = "drop"
) %>% filter(n_diet_days >= min_diet_days)

geo_person <- food0 %>% group_by(IDIND) %>% summarise(
	COMMID = first_num(COMMID), province = factor(first_num(PROVINCE)), urban = factor(first_num(URBAN)), .groups = "drop")

# c12diet
kcal_c <- pick(c12, c("D3KCAL", "KCAL"))
protein_c <- pick(c12, c("D3PROTN", "D3PROT", "PROTEIN"))
fat_c <- pick(c12, c("D3FAT", "FAT"), FALSE)
carb_c <- pick(c12, c("D3CARBO", "CARBO"), FALSE)
c12_base <- c12 %>% filter(num(.data[[wc]]) == baseline_wave)
c12b <- c12_base %>%
	transmute(IDIND = num(.data[[idc]]), kcal = num(.data[[kcal_c]]), protein_g = num(.data[[protein_c]]),
			fat_total_g = colnum(c12_base, fat_c), carb_g = colnum(c12_base, carb_c)) %>%
	group_by(IDIND) %>% summarise(across(everything(), ~ mean(.x, na.rm = TRUE)), .groups = "drop")

# pexampub12: anthropometrics, alcohol, smoking, disease indicators, BP
weight_p <- pick(pexam, c("WEIGHT", "WT"), FALSE); height_p <- pick(pexam, c("HEIGHT", "HT"), FALSE)
u25 <- pick(pexam, "U25", FALSE); u27 <- pick(pexam, "U27", FALSE); u40 <- pick(pexam, "U40", FALSE); u41 <- pick(pexam, "U41", FALSE)
u232 <- pick(pexam, "U232", FALSE); u234 <- pick(pexam, "U234", FALSE); u236 <- pick(pexam, "U236", FALSE)
u22 <- pick(pexam, "U22", FALSE); u24a <- pick(pexam, "U24A", FALSE); u24j <- pick(pexam, "U24J", FALSE)
u24l <- pick(pexam, "U24L", FALSE); u24w <- pick(pexam, "U24W", FALSE)
sys <- intersect(c("SYSTOL1", "SYSTOL2", "SYSTOL3"), names(pexam))
dia <- intersect(c("DIASTOL1", "DIASTOL2", "DIASTOL3"), names(pexam))

pex_base <- pexam %>% filter(num(.data[[wp]]) == baseline_wave)
sugar_mat <- cbind(colnum(pex_base, u232), colnum(pex_base, u234), colnum(pex_base, u236))
sugar_week <- rowSums(sugar_mat, na.rm = TRUE); sugar_week[rowSums(is.finite(sugar_mat)) == 0] <- NA_real_
sbp_base <- if(length(sys)) rowMeans(as.data.frame(lapply(pex_base[sys], num)), na.rm = TRUE) else rep(NA_real_, nrow(pex_base))
dbp_base <- if(length(dia)) rowMeans(as.data.frame(lapply(pex_base[dia], num)), na.rm = TRUE) else rep(NA_real_, nrow(pex_base))
sbp_base[is.nan(sbp_base)] <- NA_real_; dbp_base[is.nan(dbp_base)] <- NA_real_
pexb <- pex_base %>% transmute(
		IDIND = num(.data[[idp]]),
		weight_kg = colnum(pex_base, weight_p), height_cm = colnum(pex_base, height_p),
		ever_smoke = if(!is.na(u25))yes1(.data[[u25]])else NA_integer_,
		current_smoke = if(!is.na(u27))yes1(.data[[u27]])else NA_integer_,
		alcohol_drinker = if(!is.na(u40))yes1(.data[[u40]])else NA_integer_,
		alcohol_frequency = colnum(pex_base, u41),
		sugary_drinks_week = sugar_week,
		htn_self = if(!is.na(u22))yes1(.data[[u22]])else NA_integer_,
		t2dm_self = if(!is.na(u24a))yes1(.data[[u24a]])else NA_integer_,
		mi_self = if(!is.na(u24j))yes1(.data[[u24j]])else NA_integer_,
		stroke_self = if(!is.na(u24l))yes1(.data[[u24l]])else NA_integer_,
		cancer_self = if(!is.na(u24w))yes1(.data[[u24w]])else NA_integer_,
		sbp = sbp_base, dbp = dbp_base
	) %>% group_by(IDIND) %>% summarise(across(everything(), ~ first_num(.x)), .groups = "drop")


# The same longitudinal pexampub12 public-use copy supplies baseline examination data
# and later-wave prospective outcomes (including 2015 when present).
pexfu <- NULL
if(followup_end %in% pex_long_waves) {
	fu22 <- pick(pexam_long, "U22", FALSE); fu24a <- pick(pexam_long, "U24A", FALSE)
	fu24j <- pick(pexam_long, "U24J", FALSE); fu24l <- pick(pexam_long, "U24L", FALSE); fu24w <- pick(pexam_long, "U24W", FALSE)
	fu_weight <- pick(pexam_long, c("WEIGHT", "WT"), FALSE); fu_height <- pick(pexam_long, c("HEIGHT", "HT"), FALSE)
	fu_sys <- intersect(c("SYSTOL1", "SYSTOL2", "SYSTOL3"), names(pexam_long))
	fu_dia <- intersect(c("DIASTOL1", "DIASTOL2", "DIASTOL3"), names(pexam_long))
	pex_fu_base <- pexam_long %>% filter(num(.data[[wpl]]) == followup_end)
	sbp_fu_vec <- if(length(fu_sys))rowMeans(as.data.frame(lapply(pex_fu_base[fu_sys], num)), na.rm = TRUE)else rep(NA_real_, nrow(pex_fu_base))
	dbp_fu_vec <- if(length(fu_dia))rowMeans(as.data.frame(lapply(pex_fu_base[fu_dia], num)), na.rm = TRUE)else rep(NA_real_, nrow(pex_fu_base))
	sbp_fu_vec[is.nan(sbp_fu_vec)] <- NA_real_; dbp_fu_vec[is.nan(dbp_fu_vec)] <- NA_real_
	pexfu <- pex_fu_base %>% transmute(
		IDIND = num(.data[[idpl]]),
		weight_kg_fu = colnum(pex_fu_base, fu_weight), height_cm_fu = colnum(pex_fu_base, fu_height),
		sbp_fu = sbp_fu_vec, dbp_fu = dbp_fu_vec,
		htn_self_fu = if(!is.na(fu22))yes1(.data[[fu22]])else NA_integer_,
		t2dm_fu = if(!is.na(fu24a))yes1(.data[[fu24a]])else NA_integer_,
		mi_fu = if(!is.na(fu24j))yes1(.data[[fu24j]])else NA_integer_,
		stroke_fu = if(!is.na(fu24l))yes1(.data[[fu24l]])else NA_integer_,
		cancer_fu = if(!is.na(fu24w))yes1(.data[[fu24w]])else NA_integer_
	) %>% group_by(IDIND) %>% summarise(across(everything(), ~ first_num(.x)), .groups = "drop")
}

# PACT food preferences + sleep
pref_vars <- c(
	fast = pick(pact, c("U389", "U389_FF"), FALSE),
	salty = pick(pact, c("U390", "U390_SSF"), FALSE),
	fruit = pick(pact, c("U391", "U391_F"), FALSE),
	veg = pick(pact, c("U392", "U392_V"), FALSE),
	ssb = pick(pact, c("U393", "U393_SD"), FALSE)
)
sleep_var <- pick(pact, c("U324", "TIMEBED"), FALSE)

pact_base <- pact %>% filter(num(.data[[wpa]]) == baseline_wave)
pactb <- pact_base %>%
	transmute(
		IDIND = num(.data[[idpa]]),
		pref_fast = if(!is.na(pref_vars["fast"]))pref15(.data[[pref_vars["fast"]]])else NA_real_,
		pref_salty = if(!is.na(pref_vars["salty"]))pref15(.data[[pref_vars["salty"]]])else NA_real_,
		pref_fruit = if(!is.na(pref_vars["fruit"]))pref15(.data[[pref_vars["fruit"]]])else NA_real_,
		pref_veg = if(!is.na(pref_vars["veg"]))pref15(.data[[pref_vars["veg"]]])else NA_real_,
		pref_ssb = if(!is.na(pref_vars["ssb"]))pref15(.data[[pref_vars["ssb"]]])else NA_real_,
		sleep_hours = colnum(pact_base, sleep_var)
	) %>% group_by(IDIND) %>% summarise(across(everything(), ~ first_num(.x)), .groups = "drop") %>%
	mutate(
		pref_healthy = rowMeans(cbind(pref_fruit, pref_veg), na.rm = TRUE),
		pref_unhealthy = rowMeans(cbind(pref_fast, pref_salty, pref_ssb), na.rm = TRUE),
		pref_health_index = rowMeans(cbind(zstd(pref_fruit), zstd(pref_veg), - zstd(pref_fast), - zstd(pref_salty), - zstd(pref_ssb)), na.rm = TRUE)
	)
pactb$pref_healthy[!is.finite(pactb$pref_healthy)] <- NA_real_
pactb$pref_unhealthy[!is.finite(pactb$pref_unhealthy)] <- NA_real_
pactb$pref_health_index[!is.finite(pactb$pref_health_index)] <- NA_real_

# Education
educb <- NULL
if(!is.null(educ)) {
	ide <- pick(educ, "IDIND"); we <- pick(educ, "WAVE")
	a11 <- pick(educ, c("A11", "EDUYRS"), FALSE); a12 <- pick(educ, c("A12", "EDUCATION"), FALSE)
	educ_base <- educ %>% filter(num(.data[[we]]) == baseline_wave)
	schooling <- colnum(educ_base, a11)
	if (identical(a11, "A11")) {
		schooling <- recode_chns_school_years(schooling)
	}
	educb <- educ_base %>%
		transmute(IDIND = num(.data[[ide]]), education_years = schooling, education_level = colnum(educ_base, a12)) %>%
		group_by(IDIND) %>% summarise(across(everything(), ~ first_num(.x)), .groups = "drop")
}

# Household income
incb <- NULL
if(!is.null(hhinc)) {
	wh <- pick(hhinc, "WAVE"); hh <- pick(hhinc, "HHID")
	incv <- pick(hhinc, c("HHINCPC", "HHINC_PC", "HHINCPC_CPI", "HHINC_CPI"), FALSE)
	if(!is.na(incv)) {
		hhinc_base <- hhinc %>% filter(num(.data[[wh]]) == baseline_wave)
		incb <- hhinc_base %>%
			transmute(HHID = num(.data[[hh]]), income_pc = num(.data[[incv]])) %>%
			group_by(HHID) %>% summarise(income_pc = first_num(income_pc), .groups = "drop")
	}
}

# Master / mortality
gender_m <- pick(master, c("GENDER", "SEX"))
birth_m <- pick(master, c("WEST_DOB_Y", "DOB_Y", "BIRTH_Y", "MOON_DOB_Y"))
dody_m <- pick(master, c("DOD_Y", "DODY"), FALSE)
dodrpt_m <- pick(master, c("DOD_RPT", "DOD_REPORT"), FALSE)
masterb <- master %>% transmute(
	IDIND = num(.data[[idm]]), female = as.integer(num(.data[[gender_m]]) == 2),
	birth_year = num(.data[[birth_m]]), DOD_Y = colnum(master, dody_m), DOD_RPT = colnum(master, dodrpt_m)
) %>% group_by(IDIND) %>% summarise(across(everything(), ~ first_num(.x)), .groups = "drop")

# Retain both end-wave status and last observed contact year.
hh_r <- pick(roster, "HHID", FALSE)
roster_base0 <- roster %>% filter(num(.data[[wr]]) == baseline_wave)
roster_base <- roster_base0 %>%
	transmute(IDIND = num(.data[[idr]]), HHID = colnum(roster_base0, hh_r)) %>% distinct(IDIND, .keep_all = TRUE)
roster_follow <- roster %>%
	transmute(IDIND = num(.data[[idr]]), CONTACT_YEAR = num(.data[[wr]])) %>%
	filter(is.finite(IDIND), is.finite(CONTACT_YEAR), CONTACT_YEAR >= baseline_wave, CONTACT_YEAR <= followup_end) %>%
	group_by(IDIND) %>%
	summarise(
		last_contact_year = max(CONTACT_YEAR),
		observed_after_baseline = as.integer(any(CONTACT_YEAR > baseline_wave)),
		observed_at_followup = as.integer(any(CONTACT_YEAR == followup_end)),
		.groups = "drop"
	)

dat <- food_person %>%
	left_join(geo_person, by = "IDIND") %>%
	left_join(c12b, by = "IDIND") %>%
	left_join(pexb, by = "IDIND") %>%
	{if(!is.null(pexfu)) left_join(., pexfu, by = "IDIND") else .} %>%
	left_join(pactb, by = "IDIND") %>%
	left_join(masterb, by = "IDIND") %>%
	left_join(roster_base, by = "IDIND") %>%
	left_join(roster_follow, by = "IDIND")
if(!is.null(educb)) dat <- dat %>% left_join(educb, by = "IDIND")
if(!is.null(incb) && "HHID" %in% names(dat)) dat <- dat %>% left_join(incb, by = "HHID")
# Optional modules must never make the core pipeline fail simply because a covariate
# was unavailable in a particular download/format release.
for(v in c("income_pc", "education_years", "education_level", "weight_kg_fu", "height_cm_fu", "sbp_fu", "dbp_fu",
				"htn_self_fu", "t2dm_fu", "mi_fu", "stroke_fu", "cancer_fu")) {
	if(!v %in% names(dat)) dat[[v]] <- NA_real_
}

dat <- dat %>% mutate(
	age = baseline_wave - birth_year,
	bmi = ifelse(is.finite(weight_kg) & is.finite(height_cm) & height_cm > 100, weight_kg / (height_cm / 100) ^ 2, NA_real_),
	log_income = log1p(ifelse(is.finite(income_pc) & income_pc >= 0, income_pc, NA_real_)),
	protein_g_kg = protein_g / weight_kg,
	nuts_legumes = nuts + legumes,
	cereal_potato = cereal + potato_starch,
	meat_poultry_eggs = red_meat + poultry + eggs,
	fat_proxy_ratio = ifelse(is.finite(plant_oil_g) | is.finite(animal_fat_g), plant_oil_g / (plant_oil_g + animal_fat_g + 0.1), NA_real_),
	healthy_fat_ratio = ifelse(is.finite(mufa_g) | is.finite(pufa_g) | is.finite(sfat_g), (mufa_g + pufa_g) / (sfat_g + 0.1), NA_real_),
	upf_native = processed_food_g,
	upf_combined_proxy = rowMeans(cbind(zstd(processed_food_g), zstd(upf_class_proxy_g)), na.rm = TRUE),
	alcohol_proxy = alcohol_drinker
) %>%
	filter(age >= min_age, age <= 100, is.finite(kcal), kcal >= 500, kcal <= 5000)

# 🚩 5) CHDI + MAHA-CN / established-score proxy construction
step_header("CHNS diet-score construction")

# CHDI present-component version (0-70), using directly observed food-group intake.
dat <- dat %>% mutate(
	food_types_per_day = n_food_types,
	chdi_types = linear_types_score(food_types_per_day),
	chdi_cereal = linear_adequacy(cereal_potato / kcal * 1000, 140, 10),
	chdi_veg = linear_adequacy(vegetable / kcal * 1000, 270, 10),
	chdi_fruit = linear_adequacy(fruit / kcal * 1000, 110, 10),
	chdi_dairy = linear_adequacy(dairy / kcal * 1000, 100, 10),
	chdi_legnuts = linear_adequacy(nuts_legumes / kcal * 1000, 10, 10),
	chdi_meategg = linear_adequacy(meat_poultry_eggs / kcal * 1000, 50, 5),
	chdi_aquatic = linear_adequacy(fish / kcal * 1000, 30, 5),
	diet.chdi.sum = rowSums(cbind(chdi_types, chdi_cereal, chdi_veg, chdi_fruit, chdi_dairy, chdi_legnuts, chdi_meategg, chdi_aquatic), na.rm = FALSE)
)

# Build publication scores from fine-group maps when coverage passes the threshold.

# Harmonized food quantities used in several scores.
dat <- dat %>% mutate(
	red_processed_meat_g = ifelse(is.finite(red_meat) | is.finite(processed_meat),
				rowSums(cbind(red_meat, processed_meat), na.rm = TRUE), NA_real_),
	fried_fast_g = ifelse(is.finite(deepfried_g) | is.finite(fast_food_proxy),
				rowSums(cbind(deepfried_g, fast_food_proxy), na.rm = TRUE), NA_real_),
	sweets_pastries_g = ifelse(is.finite(sweet_snack_proxy_g) | is.finite(sugar_candy_preserve_g),
				rowSums(cbind(sweet_snack_proxy_g, sugar_candy_preserve_g), na.rm = TRUE), NA_real_),
	upf_pub_g = ifelse(
		is.finite(refined_grain) | is.finite(sweets_pastries_g) | is.finite(ssb_code) |
			is.finite(processed_meat) | is.finite(deepfried_g) | is.finite(fast_food_proxy) |
			is.finite(convenience_food_proxy),
		rowSums(cbind(refined_grain, sweets_pastries_g, ssb_code, processed_meat,
						deepfried_g, fast_food_proxy, convenience_food_proxy), na.rm = TRUE),
		NA_real_
	)
)

# Build prespecified sex-specific residual and density score inputs.
dat <- dat %>% mutate(
	maha_den_protein = energy_density(protein_g, kcal),
	maha_den_dairy = energy_density(dairy, kcal),
	maha_den_veg = energy_density(vegetable, kcal),
	maha_den_fruit = energy_density(fruit, kcal),
	maha_den_wholegrain = energy_density(whole_grain, kcal),
	maha_den_upf = energy_density(upf_pub_g, kcal),
	maha_den_alcohol = energy_density(alcohol_g_food, kcal),
	maha_den_sodium = energy_density(sodium_mg, kcal),
	maha_res_protein = energy_residual_grouped(protein_g, kcal, female),
	maha_res_dairy = energy_residual_grouped(dairy, kcal, female),
	maha_res_veg = energy_residual_grouped(vegetable, kcal, female),
	maha_res_fruit = energy_residual_grouped(fruit, kcal, female),
	maha_res_wholegrain = energy_residual_grouped(whole_grain, kcal, female),
	maha_res_upf = energy_residual_grouped(upf_pub_g, kcal, female),
	maha_res_alcohol = energy_residual_grouped(alcohol_g_food, kcal, female),
	maha_res_sodium = energy_residual_grouped(sodium_mg, kcal, female),
	dash_res_fruit = energy_residual_grouped(fruit, kcal, female),
	dash_res_veg = energy_residual_grouped(vegetable, kcal, female),
	dash_res_wholegrain = energy_residual_grouped(whole_grain, kcal, female),
	dash_res_lowfatdairy = energy_residual_grouped(lowfat_dairy, kcal, female),
	dash_res_nutslegumes = energy_residual_grouped(nuts_legumes, kcal, female),
	dash_res_sodium = energy_residual_grouped(if (sodium_ok) sodium_mg else salt_msg_proxy_g, kcal, female),
	dash_res_redmeat = energy_residual_grouped(red_processed_meat_g, kcal, female),
	dash_res_ssb = energy_residual_grouped(ssb_code, kcal, female),
	dash_den_fruit = energy_density(fruit, kcal),
	dash_den_veg = energy_density(vegetable, kcal),
	dash_den_wholegrain = energy_density(whole_grain, kcal),
	dash_den_lowfatdairy = energy_density(lowfat_dairy, kcal),
	dash_den_nutslegumes = energy_density(nuts_legumes, kcal),
	dash_den_sodium = energy_density(if (sodium_ok) sodium_mg else salt_msg_proxy_g, kcal),
	dash_den_redmeat = energy_density(red_processed_meat_g, kcal),
	dash_den_ssb = energy_density(ssb_code, kcal)
)

# Canonical MAHA components: same conceptual nine domains used in UKB/NHANES.
dat <- dat %>% mutate(
	maha_c_protein = qscore(protein_g_kg),
	maha_c_dairy = qscore(dairy),
	maha_c_veg = qscore(vegetable),
	maha_c_fruit = qscore(fruit),
	maha_c_wholegrain = qscore(whole_grain),
	maha_c_fat = if (fat_ok) qscore(healthy_fat_ratio) else qscore(fat_proxy_ratio),
	maha_c_upf = qscore(upf_pub_g, TRUE),
	maha_c_alcohol = if (alcohol_ok) qscore(alcohol_g_food, TRUE) else ifelse(alcohol_proxy == 0, 10, ifelse(alcohol_proxy == 1, 0, NA_real_)),
	maha_c_sodium = if (sodium_ok) qscore(sodium_mg, TRUE) else qscore(salt_msg_proxy_g, TRUE),
	maha_c_redmeat = qscore(red_processed_meat_g, TRUE),

	# Primary MAHA and sensitivity specifications aligned to UKB/NHANES naming.
	diet.maha.sum = rowmean_vecs(maha_c_protein, maha_c_dairy, maha_c_veg, maha_c_fruit,
				maha_c_wholegrain, maha_c_fat, maha_c_upf, maha_c_alcohol,
				maha_c_sodium, min_n = 6),
	# "balanced" removes the most definition-sensitive alcohol domain.
	diet.maha_bal.sum = rowmean_vecs(maha_c_protein, maha_c_dairy, maha_c_veg, maha_c_fruit,
				maha_c_wholegrain, maha_c_fat, maha_c_upf, maha_c_sodium,
				min_n = 5),
	# "strict" removes protein/dairy and adds low red/processed meat.
	diet.maha_strict.sum = rowmean_vecs(maha_c_veg, maha_c_fruit, maha_c_wholegrain, maha_c_fat,
				maha_c_upf, maha_c_alcohol, maha_c_sodium, maha_c_redmeat,
				min_n = 5),
	diet.maha_nodairy.sum = rowmean_vecs(maha_c_protein, maha_c_veg, maha_c_fruit, maha_c_wholegrain,
				maha_c_fat, maha_c_upf, maha_c_alcohol, maha_c_sodium,
				min_n = 5),
	diet.maha_noprotein.sum = rowmean_vecs(maha_c_dairy, maha_c_veg, maha_c_fruit, maha_c_wholegrain,
				maha_c_fat, maha_c_upf, maha_c_alcohol, maha_c_sodium,
				min_n = 5),

	# Energy-standardized MAHA variants. The original MAHA is retained unchanged
	# as Model 1/default reference. These alternatives test whether total intake drives the
	# association; neither weights nor cut points use mortality outcomes.
	maha_den_c_protein = qscore5_grouped(maha_den_protein, female),
	maha_den_c_dairy = qscore5_grouped(maha_den_dairy, female),
	maha_den_c_veg = qscore5_grouped(maha_den_veg, female),
	maha_den_c_fruit = qscore5_grouped(maha_den_fruit, female),
	maha_den_c_wholegrain = qscore5_grouped(maha_den_wholegrain, female),
	maha_den_c_fat = qscore5_grouped(if (fat_ok) healthy_fat_ratio else fat_proxy_ratio, female),
	maha_den_c_upf = qscore5_grouped(maha_den_upf, female, TRUE),
	maha_den_c_alcohol = if (alcohol_ok) qscore5_grouped(maha_den_alcohol, female, TRUE) else maha_c_alcohol / 2.5 + 1,
	maha_den_c_sodium = if (sodium_ok) qscore5_grouped(maha_den_sodium, female, TRUE) else qscore5_grouped(energy_density(salt_msg_proxy_g, kcal), female, TRUE),
	diet.maha_density.sum = rowmean_vecs(maha_den_c_protein, maha_den_c_dairy, maha_den_c_veg,
				maha_den_c_fruit, maha_den_c_wholegrain, maha_den_c_fat,
				maha_den_c_upf, maha_den_c_alcohol, maha_den_c_sodium,
				min_n = 6),
	maha_res_c_protein = qscore5_grouped(maha_res_protein, female),
	maha_res_c_dairy = qscore5_grouped(maha_res_dairy, female),
	maha_res_c_veg = qscore5_grouped(maha_res_veg, female),
	maha_res_c_fruit = qscore5_grouped(maha_res_fruit, female),
	maha_res_c_wholegrain = qscore5_grouped(maha_res_wholegrain, female),
	maha_res_c_fat = qscore5_grouped(if (fat_ok) healthy_fat_ratio else fat_proxy_ratio, female),
	maha_res_c_upf = qscore5_grouped(maha_res_upf, female, TRUE),
	maha_res_c_alcohol = if (alcohol_ok) qscore5_grouped(maha_res_alcohol, female, TRUE) else maha_c_alcohol / 2.5 + 1,
	maha_res_c_sodium = if (sodium_ok) qscore5_grouped(maha_res_sodium, female, TRUE) else qscore5_grouped(energy_residual_grouped(salt_msg_proxy_g, kcal, female), female, TRUE),
	diet.maha_residual.sum = rowmean_vecs(maha_res_c_protein, maha_res_c_dairy, maha_res_c_veg,
				maha_res_c_fruit, maha_res_c_wholegrain, maha_res_c_fat,
				maha_res_c_upf, maha_res_c_alcohol, maha_res_c_sodium,
				min_n = 6),

	# DASH: code-level whole grains / low-fat dairy / SSB / processed meat when mapped;
	# sodium falls back only if code-level nutrient coverage is insufficient.
	dash_c_fruit = qscore(fruit),
	dash_c_veg = qscore(vegetable),
	dash_c_wholegrain = qscore(whole_grain),
	dash_c_lowfatdairy = qscore(lowfat_dairy),
	dash_c_nutslegumes = qscore(nuts_legumes),
	dash_c_sodium = if (sodium_ok) qscore(sodium_mg, TRUE) else qscore(salt_msg_proxy_g, TRUE),
	dash_c_redmeat = qscore(red_processed_meat_g, TRUE),
	dash_c_sweets = qscore(rowSums(cbind(sweets_pastries_g, ssb_code), na.rm = TRUE), TRUE),
	diet.dash.sum = rowmean_vecs(dash_c_fruit, dash_c_veg, dash_c_wholegrain, dash_c_lowfatdairy,
				dash_c_nutslegumes, dash_c_sodium, dash_c_redmeat, dash_c_sweets,
				min_n = 5),

	# Literature-faithful Fung DASH: eight components, sex-specific quintiles,
	# energy-residual inputs, and sweetened beverages (not pastries) as the
	# negative beverage component. Mean and sum are linearly equivalent per SD.
	dash_fung_c_fruit = qscore5_grouped(dash_res_fruit, female),
	dash_fung_c_veg = qscore5_grouped(dash_res_veg, female),
	dash_fung_c_wholegrain = qscore5_grouped(dash_res_wholegrain, female),
	dash_fung_c_lowfatdairy = qscore5_grouped(dash_res_lowfatdairy, female),
	dash_fung_c_nutslegumes = qscore5_grouped(dash_res_nutslegumes, female),
	dash_fung_c_sodium = qscore5_grouped(dash_res_sodium, female, TRUE),
	dash_fung_c_redmeat = qscore5_grouped(dash_res_redmeat, female, TRUE),
	dash_fung_c_ssb = qscore5_grouped(dash_res_ssb, female, TRUE),
	diet.dash_fung.sum = rowmean_vecs(dash_fung_c_fruit, dash_fung_c_veg, dash_fung_c_wholegrain,
				dash_fung_c_lowfatdairy, dash_fung_c_nutslegumes,
				dash_fung_c_sodium, dash_fung_c_redmeat, dash_fung_c_ssb,
				min_n = 8),
	dash_den_c_fruit = qscore5_grouped(dash_den_fruit, female),
	dash_den_c_veg = qscore5_grouped(dash_den_veg, female),
	dash_den_c_wholegrain = qscore5_grouped(dash_den_wholegrain, female),
	dash_den_c_lowfatdairy = qscore5_grouped(dash_den_lowfatdairy, female),
	dash_den_c_nutslegumes = qscore5_grouped(dash_den_nutslegumes, female),
	dash_den_c_sodium = qscore5_grouped(dash_den_sodium, female, TRUE),
	dash_den_c_redmeat = qscore5_grouped(dash_den_redmeat, female, TRUE),
	dash_den_c_ssb = qscore5_grouped(dash_den_ssb, female, TRUE),
	diet.dash_density.sum = rowmean_vecs(dash_den_c_fruit, dash_den_c_veg, dash_den_c_wholegrain,
				dash_den_c_lowfatdairy, dash_den_c_nutslegumes,
				dash_den_c_sodium, dash_den_c_redmeat, dash_den_c_ssb,
				min_n = 8),

	# Use the harmonized 13-component MIND implementation.
	mind_c_green = qscore(green_leafy),
	mind_c_otherveg = qscore(pmax(vegetable - green_leafy, 0)),
	mind_c_berry = qscore(berry),
	mind_c_nuts = qscore(nuts),
	mind_c_wholegrain = qscore(whole_grain),
	mind_c_fish = qscore(fish),
	mind_c_poultry = qscore(poultry),
	mind_c_beans = qscore(legumes),
	mind_c_fat = if (fat_ok) qscore(healthy_fat_ratio) else qscore(fat_proxy_ratio),
	mind_c_meat = qscore(red_processed_meat_g, TRUE),
	mind_c_fried = qscore(fried_fast_g, TRUE),
	mind_c_sweets = qscore(sweets_pastries_g, TRUE),
	mind_c_cheese = qscore(cheese, TRUE),
	diet.mind.sum = rowmean_vecs(mind_c_green, mind_c_otherveg, mind_c_berry, mind_c_nuts,
				mind_c_wholegrain, mind_c_fish, mind_c_poultry, mind_c_beans,
				mind_c_fat, mind_c_meat, mind_c_fried, mind_c_sweets, mind_c_cheese,
				min_n = 8),

	# Mediterranean comparator, with code-level whole grain/fat/alcohol information when available.
	medi_c_fruit = qscore(fruit),
	medi_c_veg = qscore(vegetable),
	medi_c_legumes = qscore(legumes),
	medi_c_wholegrain = qscore(whole_grain),
	medi_c_nuts = qscore(nuts),
	medi_c_fish = qscore(fish),
	medi_c_fat = if (fat_ok) qscore(healthy_fat_ratio) else qscore(fat_proxy_ratio),
	medi_c_alcohol = if (alcohol_ok) moderate_alcohol_score(alcohol_g_food, female) else NA_real_,
	medi_c_meat = qscore(red_processed_meat_g, TRUE),
	medi_c_dairy = qscore(dairy, TRUE),
	diet.medi.sum = rowmean_vecs(medi_c_fruit, medi_c_veg, medi_c_legumes, medi_c_wholegrain,
				medi_c_nuts, medi_c_fish, medi_c_fat, medi_c_alcohol,
				medi_c_meat, medi_c_dairy, min_n = 6)
)

active_maha <- "diet.maha.sum"
active_label <- "MAHA"
component_keys <- c("protein", "dairy", "veg", "fruit", "wholegrain", "fat", "upf", "alcohol", "sodium")
component_cols <- paste0("maha_c_", component_keys)
component_labels <- c(
	protein = "Protein", dairy = "Dairy", veg = "Vegetables", fruit = "Fruit",
	wholegrain = "Whole grains", fat = if(fat_ok)"Healthy fats" else "Plant-vs-animal fat proxy",
	upf = "Low refined/processed foods", alcohol = if(alcohol_ok)"Low alcohol" else "No past-year alcohol proxy",
	sodium = if(sodium_ok)"Low sodium" else "Low salt/MSG seasoning proxy"
)
comparison_cols <- c("diet.dash.sum", "diet.mind.sum", "diet.medi.sum", "diet.chdi.sum")
comparison_labels <- c("DASH", "MIND", "MEDI", "CHDI")
validation_score_set <- c(
	MAHA = "diet.maha.sum",
	"MAHA-balanced" = "diet.maha_bal.sum",
	"MAHA-strict" = "diet.maha_strict.sum",
	"MAHA-no dairy" = "diet.maha_nodairy.sum",
	"MAHA-no protein" = "diet.maha_noprotein.sum",
	DASH = "diet.dash.sum",
	MIND = "diet.mind.sum",
	MEDI = "diet.medi.sum",
	CHDI = "diet.chdi.sum"
)
algorithm_score_set <- c(
	MAHA = "diet.maha.sum",
	"MAHA-density" = "diet.maha_density.sum",
	"MAHA-residual" = "diet.maha_residual.sum",
	DASH = "diet.dash_fung.sum",
	"DASH-reference" = "diet.dash.sum",
	"DASH-density" = "diet.dash_density.sum",
	MIND = "diet.mind.sum",
	MEDI = "diet.medi.sum"
)
score_algorithm_manifest <- tibble(
	Score = names(algorithm_score_set), Variable = unname(algorithm_score_set),
	Definition = c(
		"Original nine-domain MAHA; retained unchanged",
		"MAHA components per 1000 kcal with sex-specific quintiles",
		"MAHA residual energy adjustment within sex followed by sex-specific quintiles",
		"Fung eight-component DASH; residual energy adjustment, sex-specific quintiles, SSB-only negative domain",
		"Previous CHNS DASH proxy retained for exact Model 1 default reproduction",
		"Fung eight-component DASH using intake per 1000 kcal and sex-specific quintiles",
		"Existing harmonized 13-component MIND",
		"Existing harmonized Mediterranean score"
	),
	Reference = c(
		"Prespecified original MAHA implementation",
		"Prespecified MAHA energy-density sensitivity",
		"Prespecified MAHA residual-energy sensitivity",
		"Fung et al. Arch Intern Med. 2008;168:713-720. doi:10.1001/archinte.168.7.713",
		"Reference CHNS implementation retained for exact reproduction",
		"Fung domains with a commonly used per-1000-kcal sensitivity implementation",
		"Existing harmonized MIND implementation",
		"Existing harmonized Mediterranean implementation"
	),
	Outcome_informed = FALSE
)
all_score_vars <- unique(c(unname(validation_score_set), unname(algorithm_score_set)))

score_qc <- tibble(
	variable = unique(c(component_cols, active_maha, unname(algorithm_score_set))),
	finite_n = vapply(unique(c(component_cols, active_maha, unname(algorithm_score_set))),
				function(v)sum(is.finite(num(dat[[v]]))), integer(1))
)
cat("[CHNS SCORE QC] ", paste0(score_qc$variable, "=", score_qc$finite_n, collapse = "; "), "\n", sep = "")
anchor_cols <- comparison_cols

# Standardized/rescaled versions used by plots and sensitivity analyses.
all_score_vars <- unique(c(unname(validation_score_set), unname(algorithm_score_set)))
for(v in all_score_vars) if(v %in% names(dat)) {
	stem <- sub("\\.sum$", "", v); dat[[paste0(stem, ".s100")]] <- score100(dat[[v]]); dat[[paste0(stem, ".3c")]] <- f3c(dat[[v]])
}



maha_save_required_state("score", c("chns_phe.R", "chns_assoc.R", "chns_bayes.R"))
maha_export_scores(dat, "IDIND", wave = baseline_wave)
