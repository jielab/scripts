pacman::p_load(tidyverse, here, fs, config, foreach, doParallel)
scriptdir = "/work/sph-huangj/scripts/avia"


# 🚩 Simulation functions
cfig <- config::get(file = paste0(scriptdir, "/02.config.yml"), config = "dev")
n_cores <- cfig$n_cores
adherence_level <- seq(0, 1, .2)
summarize_results_column <- function(all_results, col) {
	all_results %>% dplyr::select({{col}}) %>%
	mutate(n_total = n(), n_infinite = sum(is.infinite({ {col }})), n_nan = sum(is.nan({ { col } })), n_missing = sum(is.na({ { col } })) ) %>% dplyr::filter(is.finite({ { col } })) %>%
	dplyr::summarize( metric = rlang::as_label(dplyr::enquo(col)), n = dplyr::n(), mean = mean({ { col } }, na.rm = TRUE), sd = stats::sd({ { col } }, na.rm = TRUE), median = stats::median({ { col } }, na.rm = TRUE), p025 = stats::quantile({ { col } }, na.rm = TRUE, .025), p100 = stats::quantile({ { col } }, na.rm = TRUE, .10), p250 = stats::quantile({ { col } }, na.rm = TRUE, .25), p750 = stats::quantile({ { col } }, na.rm = TRUE, .75), p900 = stats::quantile({ { col } }, na.rm = TRUE, .90), p975 = stats::quantile({ { col } }, na.rm = TRUE, .975), min = min({ { col } }, na.rm = TRUE), max = max({ { col } }, na.rm = TRUE), n_missing = mean(n_missing), n_infinite = mean(n_infinite), n_nan = mean(n_nan), n_total = mean(n_total) ) %>% dplyr::ungroup()
}
testing_scenarios <- basename(fs::dir_ls(
	here::here("intermediate_files"), type = "directory", regexp = "pcr|rapid|no_testing|perfect")
)
null_results <- readRDS(here::here("data_raw", "raw_simulations_no_testing.RDS"))
no_symp <- null_results %>% filter(symptom_screening == FALSE) %>%
	mutate(null_n_active_infection_clin = n_active_infection - n_active_infection_subclin) %>%
	dplyr::select(
		prob_inf, prop_subclin, risk_multiplier, if_threshold, time_step, round, rep,
		null_n_infected_all = n_infected_all,
		null_n_infected_new = n_infected_new,
		null_n_infection_daily = n_infection_daily,
		null_w_infection_daily = w_infection_daily,
		null_n_active_infection = n_active_infection,
		null_n_active_infection_clin,
		null_n_active_infection_subclin = n_active_infection_subclin,
		null_cume_n_infection_daily = cume_n_infection_daily,
		null_cume_w_infection_daily = cume_w_infection_daily,
		null_cume_n_infected_new = cume_n_infected_new
	)
categorize_prop_subclin <- function(all_results) {
	all_results %>% dplyr::mutate(prop_subclin_cat = factor(prop_subclin, levels = c(.3, .4), labels = c("30%", "40%"), ordered = TRUE))
}
categorize_sens_type <- function(all_results) {
	all_results %>% dplyr::mutate(sens_cat = factor(sens_type, levels = c("upper", "median"), labels = c("Upper Bound", "Median Value"), ordered = TRUE))
}
categorize_testing_types <- function(all_results) {
	all_results %>%
	mutate(testing_cat = factor(testing_type,
		levels = c( "no_testing", "no_testing_no_screening", "pcr_two_days_before", "pcr_three_days_before", "pcr_three_days_before_5_day_quarantine_pcr", "pcr_three_days_before_7_day_quarantine_pcr", "pcr_three_days_before_14_day_quarantine_pcr", "rapid_test_same_day", "rapid_same_day_5_day_quarantine_pcr", "pcr_five_days_after", "pcr_five_days_before", "pcr_seven_days_before", "pcr_two_days_before_5_day_quarantine_pcr", "pcr_five_days_before_5_day_quarantine_pcr", "pcr_seven_days_before_5_day_quarantine_pcr", "rapid_same_day_7_day_quarantine_pcr", "rapid_same_day_14_day_quarantine_pcr", "perfect_testing"),
		labels = c( "No testing", "No testing, no screening", "PCR 2 days before", "PCR 3 days before", "PCR 3 days before + 5-day quarantine", "PCR 3 days before + 7-day quarantine", "PCR 3 days before + 14-day quarantine", "Same-day Rapid Test", "Same-day Rapid Test + 5-day quarantine", "PCR 5 days after", "PCR 5 days before", "PCR 7 days before", "PCR 2 days before + 5-day quarantine", "PCR 5 days before + 5-day quarantine", "PCR 7 days before + 5-day quarantine", "Same-day Rapid Test + 7-day quarantine", "Same-day Rapid Test + 14-day quarantine", "Daily perfect testing"),
		ordered = TRUE)
	)
}
categorize_metric <- function(summarized_results) {
	summarized_results %>%
	mutate(metric_cat = factor(metric,
		levels = c("n_infection_daily", "w_infection_daily", "n_infected_all", "n_infected_new", "n_active_infection", "n_active_infection_subclin", "n_active_infection_clin", "cume_n_infection_daily", "cume_w_infection_daily", "cume_n_infected_new", "abs_cume_n_infected_new", "rel_n_infected_all", "abs_n_infection_daily", "rel_n_infection_daily", "ratio_n_infection_daily", "abs_w_infection_daily", "rel_w_infection_daily", "ratio_w_infection_daily", "abs_cume_n_infection_daily", "rel_cume_n_infection_daily", "ratio_cume_n_infection_daily", "abs_cume_w_infection_daily", "rel_cume_w_infection_daily", "ratio_cume_w_infection_daily", "frac_detected", "frac_active_detected", "any_positive_test", "n_test_true_pos", "n_test_false_pos", "ratio_false_true", "ratio_true_false", "ppv", "n_infected_day_of_flight", "n_active_infected_day_of_flight", "abs_n_active_infected_day_of_flight", "rel_n_active_infected_day_of_flight", "n_total_infections", "abs_n_infected_all", "abs_n_infected_new", "abs_n_active_infection", "abs_n_active_infection_subclin", "abs_n_active_infection_clin", "rel_n_active_infection", "rel_n_active_infection_subclin", "rel_n_active_infection_clin"),
		labels = c("Infectious days", "Weighted infections", "All infections", "New infections", "Active infections", "Active subclinical infections", "Active clinical infections", "Cumulative infectious days", "Cumulative weighted infections", "Cumulative new infections", "Abs. reduction in cumulative new infections", "Rel. reduction in total new infections", "Abs. difference in infectious days", "Rel. difference in infectious days", "Ratio of infectious days", "Abs. difference in weighted infectiousness", "Rel. difference in weighted infectiousness", "Ratio of weighted infectiousness", "Abs. difference in cumulative infection days", "Rel. difference in cumulative infection days", "Ratio of cumulative infection days", "Abs. difference in cumulative infections", "Rel. difference in cumulative infections", "Ratio of cumulative infections", "Fraction of infected detected", "Fraction of active infections detected", "Any positive test", "True positives", "False positives", "False/true positive results", "True/false positive results", "Positive Predictive Value", "Number infected on day of flight", "Number active infections on day of flight", "Abs. reduction in active infections on day of flight", "Rel. reduction in active infections on day of flight", "Total infections during observation", "Abs. reduction in total infections", "Abs. reduction in new infections", "Abs. reduction in active infections", "Abs. reduction in active subclinical infections", "Abs. reduction in active clinical infections", "Rel. reduction in active infections", "Rel. reduction in active subclinical infections", "Rel. reduction in active clinical infections"),
		ordered = TRUE)
	)
}
categorize_prob_inf <- function(all_results) {
	all_results %>% dplyr::mutate(prob_inf_cat = factor(prob_inf, levels = c(50, 100, 200, 500, 1000, 1500, 2500, 5000) / 1000000, labels = paste(c(5, 10, 20, 50, 100, 150, 250, 500), "per 100,000"), ordered = TRUE))
}
categorize_symptom_screening <- function(all_results) {
	all_results %>% dplyr::mutate(symptom_cat = factor(symptom_screening, levels = c(FALSE, TRUE), labels = c("No symptom screening", "With symptom screening"), ordered = TRUE))
}
categorize_rt_multiplier <- function(all_results) {
	all_results %>% mutate(rapid_test_cat = factor(rapid_test_multiplier, levels = c(.6, .75, .9, 1), labels = sprintf("%i%%", round(c(.6, .75, .9, 1) * 100)), ordered = TRUE))
}
categorize_risk_multiplier <- function(all_results) {
	all_results %>% mutate(risk_multi_cat = factor(risk_multiplier, levels = c(1, 2, 4, 10), labels = paste0(c(1, 2, 4, 10), "x"), ordered = TRUE))
}
shift_time_steps <- function(all_results, day_of_flight = 70) {
	all_results %>% dplyr::mutate(relative_time = time_step - day_of_flight)
}
add_cume_infections <- function(all_results) {
	all_results %>% dplyr::arrange(testing_type, symptom_screening, risk_multiplier, rapid_test_multiplier, sens_type, prob_inf, prop_subclin, if_threshold, round, rep, time_step ) %>% dplyr::group_by( testing_type, symptom_screening, risk_multiplier, rapid_test_multiplier, sens_type, prob_inf, prop_subclin, if_threshold, round, rep ) %>% dplyr::mutate( cume_n_infection_daily = cumsum(n_infection_daily), cume_w_infection_daily = cumsum(w_infection_daily), cume_n_infected_new = cumsum(n_infected_new) ) %>% dplyr::ungroup()
}
collect_results <- function(scenario_name) {
	all_files <- fs::dir_ls(paste0(outdir, "/intermediate_files/", scenario_name), recurse = TRUE, glob = "*.RDS")
	purrr::map_df(.x = all_files, .f =  ~ readRDS(.x))
}


# 🚩 Simulation workflow

outdir = "/work/sph-huangj/analysis/avia"


# 🚩 munge_intermediate_files
fs::dir_create(paste0(outdir, "/data_raw"))
testing_scenarios <- basename(fs::dir_ls(paste0(outdir, "/intermediate_files"), type = "directory", regexp = "pcr|rapid|no_testing|perfect"))
doParallel::registerDoParallel()
all_results <- foreach::foreach(s = testing_scenarios) %dopar% {
	temp_x <- collect_results(s) %>% filter(time_step >= 67, !is.na(if_threshold)) %>%
		shift_time_steps(day_of_flight = 70) %>% add_cume_infections() %>%
		categorize_testing_types() %>% categorize_prob_inf() %>% categorize_sens_type() %>% categorize_prop_subclin() %>%
		categorize_symptom_screening() %>% categorize_rt_multiplier() %>% categorize_risk_multiplier() %>%
		mutate(sim_id = sprintf("%03d.%02d", round, rep)) %>%
		dplyr::select(testing_type, testing_cat, prob_inf, prob_inf_cat, sens_type, sens_cat, prop_subclin, prop_subclin_cat, symptom_screening, symptom_cat, risk_multiplier, risk_multi_cat, rapid_test_multiplier, rapid_test_cat, if_threshold, time_step, round, rep, sim_id, dplyr::everything()) %>%
		arrange(testing_type, testing_cat, prob_inf, prob_inf_cat, prop_subclin, prop_subclin_cat, sens_type, symptom_screening, symptom_cat, risk_multiplier, risk_multi_cat, rapid_test_multiplier, rapid_test_cat, if_threshold, round, rep, time_step)
	if (!tibble::has_name(temp_x, "n_test_false_pos") & !tibble::has_name(temp_x, "n_test_true_pos")) {
		temp_x <- temp_x %>% dplyr::mutate(n_test_false_pos = NA, n_test_true_pos = NA)
	}
	s_name <- paste0(outdir, "/data_raw/", sprintf("raw_simulations_%s.RDS", s))
	saveRDS(temp_x, s_name, compress = "xz")
}
doParallel::stopImplicitCluster()


# 🚩 summarize_infection_quantities
doParallel::registerDoParallel()
temp_holder <- foreach::foreach(s = testing_scenarios) %dopar% {
	temp_x <- readRDS(here::here("data_raw", sprintf("raw_simulations_%s.RDS", s)))
	## Calculate number of active infections (clinical)
	temp_x <- temp_x %>% dplyr::mutate(n_active_infection_clin = n_active_infection - n_active_infection_subclin)
	## Join with null model (no testing, no symptom screening)
	temp_x <- temp_x %>% dplyr::left_join(no_symp)
	## Absolute differences
	temp_x <- temp_x %>% dplyr::mutate(
		abs_n_infected_all = n_infected_all - null_n_infected_all, abs_n_infected_new = n_infected_new - null_n_infected_new,
		abs_n_active_infection = n_active_infection - null_n_active_infection, abs_n_active_infection_subclin = n_active_infection_subclin - null_n_active_infection_subclin,
		abs_n_active_infection_clin = n_active_infection_clin - null_n_active_infection_clin, abs_n_infection_daily = n_infection_daily - null_n_infection_daily,
		abs_w_infection_daily = w_infection_daily - null_w_infection_daily,
		abs_cume_n_infected_new = cume_n_infected_new - null_cume_n_infected_new, abs_cume_n_infection_daily = cume_n_infection_daily - null_cume_n_infection_daily,
		abs_cume_w_infection_daily = cume_w_infection_daily - null_cume_w_infection_daily)
	## Relative differences
	temp_x <- temp_x %>% dplyr::mutate(
	rel_n_infected_all = (null_n_infected_all - n_infected_all) / null_n_infected_all,
	rel_n_infected_new = (null_n_infected_new - n_infected_new) / null_n_infected_new,
	rel_n_active_infection = (null_n_active_infection - n_active_infection) / null_n_active_infection,
	rel_n_active_infection_subclin = ( null_n_active_infection_subclin - n_active_infection_subclin ) / null_n_active_infection_subclin,
	rel_n_active_infection_clin = (null_n_active_infection_clin - n_active_infection_clin) / null_n_active_infection_clin,
	rel_n_infection_daily = (null_n_infection_daily - n_infection_daily) / null_n_infection_daily,
	rel_w_infection_daily = (null_w_infection_daily - w_infection_daily) / null_w_infection_daily,
	rel_cume_n_infected_new = (null_cume_n_infected_new - cume_n_infected_new) / null_cume_n_infected_new,
	rel_cume_n_infection_daily = (null_cume_n_infection_daily - cume_n_infection_daily) / null_cume_n_infection_daily,
	rel_cume_w_infection_daily = (null_cume_w_infection_daily - cume_w_infection_daily) / null_cume_w_infection_daily)
	## Pivot wider
	temp_x_wide <- temp_x %>% pivot_wider(
		id_cols = c(testing_type, testing_cat, prob_inf, prob_inf_cat, sens_type, sens_cat, prop_subclin, prop_subclin_cat, risk_multi_cat, risk_multiplier, rapid_test_multiplier, rapid_test_cat, if_threshold, time_step, round, rep, sim_id ),
		names_from = symptom_screening,
		values_from = c(n_infected_all, n_infected_new, n_active_infection, n_active_infection_clin, n_infection_daily, w_infection_daily, cume_n_infected_new, cume_n_infection_daily, cume_w_infection_daily, abs_n_infected_all:rel_cume_w_infection_daily )
	)
	temp_x_list <- vector("list", length = NROW(adherence_level))
	for (i in 1:NROW(adherence_level)) {
		a <- adherence_level[i]
		temp_x_list[[i]] <- temp_x_wide %>%
		dplyr::transmute(
		testing_type, testing_cat, prob_inf, prob_inf_cat, sens_type, sens_cat,
		prop_subclin, prop_subclin_cat, risk_multiplier, risk_multi_cat,
		rapid_test_multiplier, rapid_test_cat, if_threshold, time_step, round, rep, sim_id,
		symptom_adherence = a,
		n_infected_all = (n_infected_all_TRUE * a) + (n_infected_all_FALSE * (1 - a)),
		n_infected_new = (n_infected_new_TRUE * a) + (n_infected_new_FALSE * (1 - a)),
		n_active_infection = (n_active_infection_TRUE * a) + (n_active_infection_FALSE * (1 - a)),
		n_active_infection_clin = (n_active_infection_clin_TRUE * a) + (n_active_infection_clin_FALSE * (1 - a)),
		n_infection_daily = (n_infection_daily_TRUE * a) + (n_infection_daily_FALSE * (1 - a)),
		w_infection_daily = (w_infection_daily_TRUE * a) + (w_infection_daily_FALSE * (1 - a)),
		cume_n_infected_new = (cume_n_infected_new_TRUE * a) + (cume_n_infected_new_FALSE * (1 - a)),
		cume_n_infection_daily = (cume_n_infection_daily_TRUE * a) + (cume_n_infection_daily_FALSE * (1 - a)),
		cume_w_infection_daily = (cume_w_infection_daily_TRUE * a) + (cume_w_infection_daily_FALSE * (1 - a)),
		abs_n_infected_all = (abs_n_infected_all_TRUE * a) + (abs_n_infected_all_FALSE * (1 - a)),
		abs_n_infected_new = (abs_n_infected_new_TRUE * a) + (abs_n_infected_new_FALSE * (1 - a)),
		abs_n_active_infection = (abs_n_active_infection_TRUE * a) + (abs_n_active_infection_FALSE * (1 - a)),
		abs_n_active_infection_subclin = (abs_n_active_infection_subclin_TRUE * a) + (abs_n_active_infection_subclin_FALSE * (1 - a)),
		abs_n_active_infection_clin = (abs_n_active_infection_clin_TRUE * a) + (abs_n_active_infection_clin_FALSE * (1 - a)),
		abs_n_infection_daily = (abs_n_infection_daily_TRUE * a) + (abs_n_infection_daily_FALSE * (1 - a)),
		abs_w_infection_daily = (abs_w_infection_daily_TRUE * a) + (abs_w_infection_daily_FALSE * (1 - a)),
		abs_cume_n_infected_new = (abs_cume_n_infected_new_TRUE * a) + (abs_cume_n_infected_new_FALSE * (1 - a)),
		abs_cume_n_infection_daily = (abs_cume_n_infection_daily_TRUE * a) + (abs_cume_n_infection_daily_FALSE * (1 - a)),
		abs_cume_w_infection_daily = (abs_cume_w_infection_daily_TRUE * a) + (abs_cume_w_infection_daily_FALSE * (1 - a)),
		rel_n_infected_all = (rel_n_infected_all_TRUE * a) + (rel_n_infected_all_FALSE * (1 - a)),
		rel_n_infected_new = (rel_n_infected_new_TRUE * a) + (rel_n_infected_new_FALSE * (1 - a)),
		rel_n_active_infection = (rel_n_active_infection_TRUE * a) + (rel_n_active_infection_FALSE * (1 - a)),
		rel_n_active_infection_subclin = (rel_n_active_infection_subclin_TRUE * a) + (rel_n_active_infection_subclin_FALSE * (1 - a)),
		rel_n_active_infection_clin = (rel_n_active_infection_clin_TRUE * a) + (rel_n_active_infection_clin_FALSE * (1 - a)),
		rel_n_infection_daily = (rel_n_infection_daily_TRUE * a) + (rel_n_infection_daily_FALSE * (1 - a)),
		rel_w_infection_daily = (rel_w_infection_daily_TRUE * a) + (rel_w_infection_daily_FALSE * (1 - a)),
		rel_cume_n_infected_new = (rel_cume_n_infected_new_TRUE * a) + (rel_cume_n_infected_new_FALSE * (1 - a)),
		rel_cume_n_infection_daily = (rel_cume_n_infection_daily_TRUE * a) + (rel_cume_n_infection_daily_FALSE * (1 - a)),
		rel_cume_w_infection_daily = (rel_cume_w_infection_daily_TRUE * a) + (rel_cume_w_infection_daily_FALSE * (1 - a)) )
	}
	saveRDS(dplyr::bind_rows(temp_x_list), here::here("data_raw", sprintf("processed_simulations_wide_%s.RDS", s) ), compress = "xz")
}
doParallel::stopImplicitCluster()
# Summarize all results that do not incorporate quarantine
doParallel::registerDoParallel()
results_no_quarantine <-
	foreach::foreach(s = testing_scenarios[!grepl("quarantine", testing_scenarios)]) %dopar% {
		temp_x <- readRDS(here::here( "data_raw", sprintf("processed_simulations_wide_%s.RDS", s) ))
		temp_x <- temp_x %>% dplyr::group_by( testing_type, testing_cat, prob_inf, prob_inf_cat, sens_type, sens_cat, prop_subclin, prop_subclin_cat, symptom_adherence, risk_multiplier, risk_multi_cat, rapid_test_multiplier, rapid_test_cat, if_threshold, time_step )
		dplyr::bind_rows(
		temp_x %>% summarize_results_column(n_infected_all),
		temp_x %>% summarize_results_column(n_infected_new),
		temp_x %>% summarize_results_column(n_active_infection),
		temp_x %>% summarize_results_column(n_active_infection_clin),
		temp_x %>% summarize_results_column(n_infection_daily),
		temp_x %>% summarize_results_column(w_infection_daily),
		temp_x %>% summarize_results_column(cume_n_infected_new),
		temp_x %>% summarize_results_column(cume_n_infection_daily),
		temp_x %>% summarize_results_column(cume_w_infection_daily),
		temp_x %>% summarize_results_column(abs_n_infected_all),
		temp_x %>% summarize_results_column(abs_n_infected_new),
		temp_x %>% summarize_results_column(abs_n_active_infection),
		temp_x %>% summarize_results_column(abs_n_active_infection_subclin),
		temp_x %>% summarize_results_column(abs_n_active_infection_clin),
		temp_x %>% summarize_results_column(abs_n_infection_daily),
		temp_x %>% summarize_results_column(abs_w_infection_daily),
		temp_x %>% summarize_results_column(abs_cume_n_infected_new),
		temp_x %>% summarize_results_column(abs_cume_n_infection_daily),
		temp_x %>% summarize_results_column(abs_cume_w_infection_daily),
		temp_x %>% summarize_results_column(rel_n_infected_all),
		temp_x %>% summarize_results_column(rel_n_infected_new),
		temp_x %>% summarize_results_column(rel_n_active_infection),
		temp_x %>% summarize_results_column(rel_n_active_infection_subclin),
		temp_x %>% summarize_results_column(rel_n_active_infection_clin),
		temp_x %>% summarize_results_column(rel_n_infection_daily),
		temp_x %>% summarize_results_column(rel_w_infection_daily),
		temp_x %>% summarize_results_column(rel_cume_n_infected_new),
		temp_x %>% summarize_results_column(rel_cume_n_infection_daily),
		temp_x %>% summarize_results_column(rel_cume_w_infection_daily)
		) %>% dplyr::ungroup() %>% dplyr::mutate(quarantine_adherence = 0)
	}
doParallel::stopImplicitCluster()
## Now like above, we want to take a weighted average using different weights
## to estimate different levels of adherence to quarantining.
quarantine_comparisons <- dplyr::bind_rows(
	expand.grid(base_case = "rapid_test_same_day", comparison_case = c( "rapid_same_day_5_day_quarantine_pcr", "rapid_same_day_7_day_quarantine_pcr", "rapid_same_day_14_day_quarantine_pcr"), stringsAsFactors = FALSE),
	expand.grid(base_case = "pcr_three_days_before", comparison_case = c( "pcr_three_days_before_5_day_quarantine_pcr", "pcr_three_days_before_7_day_quarantine_pcr", "pcr_three_days_before_14_day_quarantine_pcr"), stringsAsFactors = FALSE),
) %>%
	tibble::add_case(base_case = "pcr_five_days_before", comparison_case = "pcr_five_days_before_5_day_quarantine_pcr") %>%
	tibble::add_case(base_case = "pcr_seven_days_before", comparison_case = "pcr_seven_days_before_5_day_quarantine_pcr") %>%
	tibble::add_case(base_case = "pcr_two_days_before", comparison_case = "pcr_two_days_before_5_day_quarantine_pcr") %>%
	tibble::add_case(base_case = "pcr_five_days_after", comparison_case = "5_day_quarantine_pcr_five_days_after")
doParallel::registerDoParallel()
results_quarantine <- foreach::foreach(i = 1:NROW(quarantine_comparisons)) %dopar% {
	b <- quarantine_comparisons$base_case[i]
	comp <- quarantine_comparisons$comparison_case[i]
	base_x <- readRDS(here::here("data_raw", sprintf("processed_simulations_wide_%s.RDS", b) ))
	comp_x <- readRDS(here::here("data_raw", sprintf("processed_simulations_wide_%s.RDS", comp) ))
	joined_x <- dplyr::left_join( comp_x, base_x %>% dplyr::select( - testing_type, - testing_cat), by = c( "prob_inf", "prob_inf_cat", "sens_type", "sens_cat", "prop_subclin", "prop_subclin_cat", "risk_multiplier", "risk_multi_cat", "rapid_test_multiplier", "rapid_test_cat", "if_threshold", "time_step", "round", "rep", "sim_id", "symptom_adherence" ))
	temp_x_list <- vector("list", length = NROW(adherence_level))
	for (i in 1:NROW(adherence_level)) {
		a <- adherence_level[i]
		temp_x_list[[i]] <- joined_x %>%
		dplyr::transmute(testing_type, testing_cat, prob_inf, prob_inf_cat, sens_type, sens_cat, prop_subclin, prop_subclin_cat, risk_multiplier, risk_multi_cat, rapid_test_multiplier, rapid_test_cat, if_threshold, time_step,
		round, rep, sim_id, symptom_adherence,
		quarantine_adherence = a,
		n_infected_all = (n_infected_all.x * a) + (n_infected_all.y * (1 - a)),
		n_infected_new = (n_infected_new.x * a) + (n_infected_new.y * (1 - a)),
		n_active_infection = (n_active_infection.x * a) + (n_active_infection.y * (1 - a)),
		n_active_infection_clin = (n_active_infection_clin.x * a) + (n_active_infection_clin.y * (1 - a)),
		n_infection_daily = (n_infection_daily.x * a) + (n_infection_daily.y * (1 - a)),
		w_infection_daily = (w_infection_daily.x * a) + (w_infection_daily.y * (1 - a)),
		cume_n_infected_new = (cume_n_infected_new.x * a) + (cume_n_infected_new.y * (1 - a)),
		cume_n_infection_daily = (cume_n_infection_daily.x * a) + (cume_n_infection_daily.y * (1 - a)),
		cume_w_infection_daily = (cume_w_infection_daily.x * a) + (cume_w_infection_daily.y * (1 - a)),
		abs_n_infected_all = (abs_n_infected_all.x * a) + (abs_n_infected_all.y * (1 - a)),
		abs_n_infected_new = (abs_n_infected_new.x * a) + (abs_n_infected_new.y * (1 - a)),
		abs_n_active_infection = (abs_n_active_infection.x * a) + (abs_n_active_infection.y * (1 - a)),
		abs_n_active_infection_subclin = (abs_n_active_infection_subclin.x * a) + (abs_n_active_infection_subclin.y * (1 - a)),
		abs_n_active_infection_clin = (abs_n_active_infection_clin.x * a) + (abs_n_active_infection_clin.y * (1 - a)),
		abs_n_infection_daily = (abs_n_infection_daily.x * a) + (abs_n_infection_daily.y * (1 - a)),
		abs_w_infection_daily = (abs_w_infection_daily.x * a) + (abs_w_infection_daily.y * (1 - a)),
		abs_cume_n_infected_new = (abs_cume_n_infected_new.x * a) + (abs_cume_n_infected_new.y * (1 - a)),
		abs_cume_n_infection_daily = (abs_cume_n_infection_daily.x * a) + (abs_cume_n_infection_daily.y * (1 - a)),
		abs_cume_w_infection_daily = (abs_cume_w_infection_daily.x * a) + (abs_cume_w_infection_daily.y * (1 - a)),
		rel_n_infected_all = (rel_n_infected_all.x * a) + (rel_n_infected_all.y * (1 - a)),
		rel_n_infected_new = (rel_n_infected_new.x * a) + (rel_n_infected_new.y * (1 - a)),
		rel_n_active_infection = (rel_n_active_infection.x * a) + (rel_n_active_infection.y * (1 - a)),
		rel_n_active_infection_subclin = (rel_n_active_infection_subclin.x * a) + (rel_n_active_infection_subclin.y * (1 - a)),
		rel_n_active_infection_clin = (rel_n_active_infection_clin.x * a) + (rel_n_active_infection_clin.y * (1 - a)),
		rel_n_infection_daily = (rel_n_infection_daily.x * a) + (rel_n_infection_daily.y * (1 - a)),
		rel_w_infection_daily = (rel_w_infection_daily.x * a) + (rel_w_infection_daily.y * (1 - a)),
		rel_cume_n_infected_new = (rel_cume_n_infected_new.x * a) + (rel_cume_n_infected_new.y * (1 - a)),
		rel_cume_n_infection_daily = (rel_cume_n_infection_daily.x * a) + (rel_cume_n_infection_daily.y * (1 - a)),
		rel_cume_w_infection_daily = (rel_cume_w_infection_daily.x * a) + (rel_cume_w_infection_daily.y * (1 - a)))
	}
	rm(base_x, comp_x); rm(joined_x); invisible(gc(FALSE))
	## Group these up and then use summarize function
	temp_x <- temp_x_list %>% dplyr::bind_rows() %>% dplyr::group_by( testing_type, testing_cat, prob_inf, prob_inf_cat, sens_type, sens_cat, prop_subclin, prop_subclin_cat, symptom_adherence, quarantine_adherence, risk_multiplier, risk_multi_cat, rapid_test_multiplier, rapid_test_cat, if_threshold, time_step )
	## summarize_results_column() takes a single column and returns
	## descriptive summary stats for that column (by the grouping variables)
	dplyr::bind_rows(
		temp_x %>% summarize_results_column(n_infected_all),
		temp_x %>% summarize_results_column(n_infected_new),
		temp_x %>% summarize_results_column(n_active_infection),
		temp_x %>% summarize_results_column(n_active_infection_clin),
		temp_x %>% summarize_results_column(n_infection_daily),
		temp_x %>% summarize_results_column(w_infection_daily),
		temp_x %>% summarize_results_column(cume_n_infected_new),
		temp_x %>% summarize_results_column(cume_n_infection_daily),
		temp_x %>% summarize_results_column(cume_w_infection_daily),
		temp_x %>% summarize_results_column(abs_n_infected_all),
		temp_x %>% summarize_results_column(abs_n_infected_new),
		temp_x %>% summarize_results_column(abs_n_active_infection),
		temp_x %>% summarize_results_column(abs_n_active_infection_subclin),
		temp_x %>% summarize_results_column(abs_n_active_infection_clin),
		temp_x %>% summarize_results_column(abs_n_infection_daily),
		temp_x %>% summarize_results_column(abs_w_infection_daily),
		temp_x %>% summarize_results_column(abs_cume_n_infected_new),
		temp_x %>% summarize_results_column(abs_cume_n_infection_daily),
		temp_x %>% summarize_results_column(abs_cume_w_infection_daily),
		temp_x %>% summarize_results_column(rel_n_infected_all),
		temp_x %>% summarize_results_column(rel_n_infected_new),
		temp_x %>% summarize_results_column(rel_n_active_infection),
		temp_x %>% summarize_results_column(rel_n_active_infection_subclin),
		temp_x %>% summarize_results_column(rel_n_active_infection_clin),
		temp_x %>% summarize_results_column(rel_n_infection_daily),
		temp_x %>% summarize_results_column(rel_w_infection_daily),
		temp_x %>% summarize_results_column(rel_cume_n_infected_new),
		temp_x %>% summarize_results_column(rel_cume_n_infection_daily),
		temp_x %>% summarize_results_column(rel_cume_w_infection_daily)
	) %>% dplyr::ungroup()
}
doParallel::stopImplicitCluster()
summarized_results <- dplyr::bind_rows(results_no_quarantine, results_quarantine)
saveRDS(summarized_results, here::here("data", "summarized_results.RDS"), compress = "xz")


# 🚩 summarize_testing_quantites
adherence_level <- seq(0, 1, .2) ## Different levels of adherence to symptom screening
testing_dict <- list(
	no_testing = 70,
	pcr_two_days_before = 68,
	pcr_two_days_before_5_day_quarantine_pcr = c(70, 75),
	pcr_three_days_before = 68,
	pcr_three_days_before_5_day_quarantine_pcr = c(70, 75),
	pcr_three_days_before_7_day_quarantine_pcr = c(70, 77),
	pcr_three_days_before_14_day_quarantine_pcr = c(70, 84),
	pcr_five_days_before = 66,
	pcr_five_days_before_5_day_quarantine_pcr = c(70, 75),
	pcr_seven_days_before = 64,
	pcr_seven_days_before_5_day_quarantine_pcr = c(70, 75),
	rapid_test_same_day = 70,
	rapid_same_day_5_day_quarantine_pcr = c(70, 75),
	rapid_same_day_7_day_quarantine_pcr = c(70, 77),
	rapid_same_day_14_day_quarantine_pcr = c(70, 84),
	pcr_five_days_after = 75,
	"5_day_quarantine_pcr_five_days_after" = 75
)
## Calculate testing quantities ----
doParallel::registerDoParallel(cores = n_cores)
testing_results <- foreach::foreach(i = 1:NROW(testing_dict)) %dopar% {
	testing_type <- names(testing_dict)[i]
	time_steps <- testing_dict[[testing_type]]
	temp_x <- readRDS(here::here("data_raw", sprintf("raw_simulations_%s.RDS", testing_type))) %>%
		dplyr::group_by(
			testing_type, testing_cat, prob_inf, prob_inf_cat, sens_type, sens_cat, prop_subclin,
			prop_subclin_cat, symptom_screening, risk_multiplier, risk_multi_cat, rapid_test_multiplier, rapid_test_cat, if_threshold, round, rep, sim_id
		)
	## Get day of flight infections and total infections
	temp_x_infs <- dplyr::left_join(
		temp_x %>% dplyr::filter(time_step == 70) %>% dplyr::select( n_infected_day_of_flight = n_infected_all, n_active_infected_day_of_flight = n_active_infection, n_active_infected_subclin_day_of_flight = n_active_infection_subclin ),
		temp_x %>% dplyr::summarize(n_total_infections = n_susceptible[time_step == 84] - n_susceptible[time_step == 67] + n_infected_all[time_step == 67]) %>% dplyr::ungroup() %>% dplyr::distinct()
	)
	## Get testing statistics based on the scenario dictionary because some tests are *before* t-3 (when we started counting cumulative infections), we need to get the raw files instead of processed files. 
	temp_x_test <- collect_results(testing_type) %>%
		dplyr::group_by( testing_type, prob_inf, sens_type, prop_subclin, symptom_screening, risk_multiplier, rapid_test_multiplier, if_threshold, round, rep) %>%
		dplyr::filter(time_step %in% time_steps)
	if (!tibble::has_name(temp_x_test, "n_test_false_pos")) {
		temp_x_test <- temp_x_test %>% dplyr::mutate(n_test_false_pos = NA_integer_, n_test_true_pos = NA_integer_)
	}
	temp_x_test <- temp_x_test %>%
		dplyr::select(time_step, n_test_false_pos, n_test_true_pos) %>%
		dplyr::summarize(
			n_test_false_pos_first = dplyr::case_when(
				testing_type == "no_testing" ~ NA_integer_,
				dplyr::n_distinct(time_steps) == 1 ~ as.integer(max(n_test_false_pos)),
				dplyr::n_distinct(time_steps) == 2 ~ as.integer(n_test_false_pos[time_step == min(time_steps)]) ), n_test_false_pos_second = dplyr::case_when( testing_type == "no_testing" ~ NA_integer_, dplyr::n_distinct(time_steps) == 1 ~ NA_integer_, dplyr::n_distinct(time_steps) == 2 ~ as.integer(n_test_false_pos[time_step == max(time_steps)]) ),
			n_test_true_pos_first = dplyr::case_when(
				testing_type == "no_testing" ~ NA_integer_,
				dplyr::n_distinct(time_steps) == 1 ~ as.integer(max(n_test_true_pos)),
				dplyr::n_distinct(time_steps) == 2 ~ as.integer(n_test_true_pos[time_step == min(time_steps)])
			),
			n_test_true_pos_second = dplyr::case_when(
				testing_type == "no_testing" ~ NA_integer_,
				dplyr::n_distinct(time_steps) == 1 ~ NA_integer_,
				dplyr::n_distinct(time_steps) == 2 ~ as.integer(n_test_true_pos[time_step == max(time_steps)])
			)
		) %>% ungroup() %>% dplyr::distinct() %>% dplyr::rowwise() %>%
		mutate(
			n_test_false_pos = dplyr::case_when(testing_type == "no_testing" ~ NA_integer_, TRUE ~ sum(n_test_false_pos_first, n_test_false_pos_second, na.rm = TRUE) ),
			n_test_true_pos = dplyr::case_when(testing_type == "no_testing" ~ NA_integer_, TRUE ~ sum(n_test_true_pos_first, n_test_true_pos_second, na.rm = TRUE) )
		) %>% ungroup()
	rm(temp_x); invisible(gc(FALSE))
	## Pivot wider to we can get "adherence" to symptom screening. This shouldn't change anything but reviewers want it. 
	temp_x_wide <- temp_x_infs %>% dplyr::left_join(temp_x_test) %>%
	tidyverse::pivot_wider(
		id_cols = c(testing_type, testing_cat, prob_inf, prob_inf_cat, sens_type, sens_cat, prop_subclin, prop_subclin_cat, risk_multi_cat, risk_multiplier, rapid_test_multiplier, rapid_test_cat, if_threshold, round, rep, sim_id ), names_from = symptom_screening, values_from = n_infected_day_of_flight:n_test_true_pos
	)
	## Loop through different levels of adherence to symptom screening. 
	temp_x_list <- vector("list", length = NROW(adherence_level))
	for (i in 1:NROW(adherence_level)) {
		a <- adherence_level[i]
		temp_x_list[[i]] <- temp_x_wide %>% dplyr::transmute(
		testing_type, testing_cat, prob_inf, prob_inf_cat, sens_type, sens_cat, prop_subclin, prop_subclin_cat, risk_multiplier, risk_multi_cat, rapid_test_multiplier, rapid_test_cat, if_threshold, round, rep, sim_id, symptom_adherence = a,
		n_infected_day_of_flight = round((n_infected_day_of_flight_TRUE * a) + (n_infected_day_of_flight_FALSE * (1 - a)) ),
		n_active_infected_day_of_flight = round((n_active_infected_day_of_flight_TRUE * a) + (n_active_infected_day_of_flight_FALSE * (1 - a)) ),
		n_active_infected_subclin_day_of_flight = round((n_active_infected_subclin_day_of_flight_TRUE * a) + (n_active_infected_subclin_day_of_flight_FALSE * (1 - a)) ),
		n_total_infections = round((n_total_infections_TRUE * a) + (n_total_infections_FALSE * (1 - a))),
		n_test_false_pos_first = round((n_test_false_pos_first_TRUE * a) + (n_test_false_pos_first_FALSE * (1 - a)) ),
		n_test_false_pos_second = round((n_test_false_pos_second_TRUE * a) + (n_test_false_pos_second_FALSE * (1 - a)) ),
		n_test_true_pos_first = round((n_test_true_pos_first_TRUE * a) + (n_test_true_pos_first_FALSE * (1 - a)) ),
		n_test_true_pos_second = round((n_test_true_pos_second_TRUE * a) + (n_test_true_pos_second_FALSE * (1 - a)) ),
		n_test_false_pos = round((n_test_false_pos_TRUE * a) + (n_test_false_pos_FALSE * (1 - a))),
		n_test_true_pos = round((n_test_true_pos_TRUE * a) + (n_test_true_pos_FALSE * (1 - a)))
		)
	}
	## Calculate testing quantities we are interested in but don't summarize
	temp_x_list %>% dplyr::bind_rows() %>% dplyr::mutate( any_positive_test = n_test_false_pos + n_test_true_pos, ratio_false_true = n_test_false_pos / n_test_true_pos, ratio_true_false = n_test_true_pos / n_test_false_pos, ppv = n_test_true_pos / (n_test_true_pos + n_test_false_pos), frac_detected = n_test_true_pos / n_infected_day_of_flight, frac_active_detected = n_test_true_pos / n_active_infected_day_of_flight, time_step = 84)
}
doParallel::stopImplicitCluster()
closeAllConnections()
testing_results <- dplyr::bind_rows(testing_results)
saveRDS(testing_results, here::here("data", "all_testing_results.RDS"), compress = "xz")
# Summarize these results
testing_results <- testing_results %>%
	dplyr::group_by(
		testing_type, testing_cat, prob_inf, prob_inf_cat, sens_type, sens_cat, prop_subclin, prop_subclin_cat,
		symptom_adherence, risk_multiplier, risk_multi_cat, rapid_test_multiplier, rapid_test_cat, if_threshold, time_step
	)
testing_summary <- dplyr::bind_rows(
	testing_results %>% summarize_results_column(any_positive_test),
	testing_results %>% summarize_results_column(n_test_true_pos),
	testing_results %>% summarize_results_column(n_test_false_pos),
	testing_results %>% summarize_results_column(n_test_true_pos_first),
	testing_results %>% summarize_results_column(n_test_false_pos_first),
	testing_results %>% summarize_results_column(n_test_true_pos_second),
	testing_results %>% summarize_results_column(n_test_false_pos_second),
	testing_results %>% summarize_results_column(n_active_infected_day_of_flight),
	testing_results %>% summarize_results_column(n_active_infected_subclin_day_of_flight),
	testing_results %>% summarize_results_column(n_infected_day_of_flight),
	testing_results %>% summarize_results_column(ratio_false_true),
	testing_results %>% summarize_results_column(ratio_true_false),
	testing_results %>% summarize_results_column(ppv),
	testing_results %>% summarize_results_column(frac_detected),
	testing_results %>% summarize_results_column(frac_active_detected),
	testing_results %>% summarize_results_column(n_total_infections)
) %>% dplyr::ungroup()
saveRDS(testing_summary, here::here("data", "summarized_testing_results.RDS"), compress = "xz")

