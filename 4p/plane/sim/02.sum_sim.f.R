cfig <- config::get(file=paste0(scriptdir,"/02.config.yml"), config="dev")
n_cores <- cfig$n_cores
adherence_level <- seq(0, 1, .2)
summarize_results_column <- function(all_results, col) {
	all_results %>% dplyr::select({{col}}) %>% 
	mutate(n_total=n(), n_infinite=sum(is.infinite({ {col }})), n_nan=sum(is.nan({ { col } })), n_missing=sum(is.na({ { col } })) ) %>% dplyr::filter(is.finite({ { col } })) %>% 
	dplyr::summarize( metric=rlang::as_label(dplyr::enquo(col)), n=dplyr::n(), mean=mean({ { col } }, na.rm=TRUE), sd=stats::sd({ { col } }, na.rm=TRUE), median=stats::median({ { col } }, na.rm=TRUE), p025=stats::quantile({ { col } }, na.rm=TRUE, .025), p100=stats::quantile({ { col } }, na.rm=TRUE, .10), p250=stats::quantile({ { col } }, na.rm=TRUE, .25), p750=stats::quantile({ { col } }, na.rm=TRUE, .75), p900=stats::quantile({ { col } }, na.rm=TRUE, .90), p975=stats::quantile({ { col } }, na.rm=TRUE, .975), min=min({ { col } }, na.rm=TRUE), max=max({ { col } }, na.rm=TRUE), n_missing=mean(n_missing), n_infinite=mean(n_infinite), n_nan=mean(n_nan), n_total=mean(n_total) ) %>% dplyr::ungroup()
}
testing_scenarios <- basename(fs::dir_ls(
	here::here("intermediate_files"),type="directory",regexp="pcr|rapid|no_testing|perfect")
)
null_results <- readRDS(here::here("data_raw", "raw_simulations_no_testing.RDS"))
no_symp <- null_results %>% filter(symptom_screening == FALSE) %>%
	mutate(null_n_active_infection_clin=n_active_infection - n_active_infection_subclin) %>%
	dplyr::select(
		prob_inf,prop_subclin,risk_multiplier,if_threshold,time_step, round,rep,
		null_n_infected_all=n_infected_all,
		null_n_infected_new=n_infected_new,
		null_n_infection_daily=n_infection_daily,
		null_w_infection_daily=w_infection_daily,
		null_n_active_infection=n_active_infection,
		null_n_active_infection_clin,
		null_n_active_infection_subclin=n_active_infection_subclin, 
		null_cume_n_infection_daily=cume_n_infection_daily, 
		null_cume_w_infection_daily=cume_w_infection_daily, 
		null_cume_n_infected_new=cume_n_infected_new
	)
categorize_prop_subclin <- function(all_results) {
	all_results %>% dplyr::mutate(prop_subclin_cat=factor(prop_subclin,levels=c(.3,.4),labels=c("30%", "40%"),ordered=TRUE))
}
categorize_sens_type <- function(all_results) {
	all_results %>% dplyr::mutate(sens_cat=factor(sens_type,levels=c("upper", "median"),labels=c("Upper Bound","Median Value"),ordered=TRUE))
}
categorize_testing_types <- function(all_results) {
	all_results %>%
	mutate(testing_cat=factor(testing_type,
		levels=c( "no_testing", "no_testing_no_screening", "pcr_two_days_before", "pcr_three_days_before", "pcr_three_days_before_5_day_quarantine_pcr", "pcr_three_days_before_7_day_quarantine_pcr", "pcr_three_days_before_14_day_quarantine_pcr", "rapid_test_same_day", "rapid_same_day_5_day_quarantine_pcr", "pcr_five_days_after", "pcr_five_days_before", "pcr_seven_days_before", "pcr_two_days_before_5_day_quarantine_pcr", "pcr_five_days_before_5_day_quarantine_pcr", "pcr_seven_days_before_5_day_quarantine_pcr", "rapid_same_day_7_day_quarantine_pcr", "rapid_same_day_14_day_quarantine_pcr", "perfect_testing"),
		labels=c( "No testing", "No testing, no screening", "PCR 2 days before", "PCR 3 days before", "PCR 3 days before + 5-day quarantine", "PCR 3 days before + 7-day quarantine", "PCR 3 days before + 14-day quarantine", "Same-day Rapid Test", "Same-day Rapid Test + 5-day quarantine", "PCR 5 days after", "PCR 5 days before", "PCR 7 days before", "PCR 2 days before + 5-day quarantine", "PCR 5 days before + 5-day quarantine", "PCR 7 days before + 5-day quarantine", "Same-day Rapid Test + 7-day quarantine", "Same-day Rapid Test + 14-day quarantine", "Daily perfect testing"),
		ordered=TRUE)
	)
}
categorize_metric <- function(summarized_results) {
	summarized_results %>%
	mutate(metric_cat =factor(metric,
		levels=c("n_infection_daily", "w_infection_daily", "n_infected_all", "n_infected_new", "n_active_infection", "n_active_infection_subclin", "n_active_infection_clin", "cume_n_infection_daily", "cume_w_infection_daily", "cume_n_infected_new", "abs_cume_n_infected_new", "rel_n_infected_all", "abs_n_infection_daily", "rel_n_infection_daily", "ratio_n_infection_daily", "abs_w_infection_daily", "rel_w_infection_daily", "ratio_w_infection_daily", "abs_cume_n_infection_daily", "rel_cume_n_infection_daily", "ratio_cume_n_infection_daily", "abs_cume_w_infection_daily", "rel_cume_w_infection_daily", "ratio_cume_w_infection_daily", "frac_detected", "frac_active_detected", "any_positive_test", "n_test_true_pos", "n_test_false_pos", "ratio_false_true", "ratio_true_false", "ppv", "n_infected_day_of_flight", "n_active_infected_day_of_flight", "abs_n_active_infected_day_of_flight", "rel_n_active_infected_day_of_flight", "n_total_infections", "abs_n_infected_all", "abs_n_infected_new", "abs_n_active_infection", "abs_n_active_infection_subclin", "abs_n_active_infection_clin", "rel_n_active_infection", "rel_n_active_infection_subclin", "rel_n_active_infection_clin"),
		labels=c("Infectious days", "Weighted infections", "All infections", "New infections", "Active infections", "Active subclinical infections", "Active clinical infections", "Cumulative infectious days", "Cumulative weighted infections", "Cumulative new infections", "Abs. reduction in cumulative new infections", "Rel. reduction in total new infections", "Abs. difference in infectious days", "Rel. difference in infectious days", "Ratio of infectious days", "Abs. difference in weighted infectiousness", "Rel. difference in weighted infectiousness", "Ratio of weighted infectiousness", "Abs. difference in cumulative infection days", "Rel. difference in cumulative infection days", "Ratio of cumulative infection days", "Abs. difference in cumulative infections", "Rel. difference in cumulative infections", "Ratio of cumulative infections", "Fraction of infected detected", "Fraction of active infections detected", "Any positive test", "True positives", "False positives", "False/true positive results", "True/false positive results", "Positive Predictive Value", "Number infected on day of flight", "Number active infections on day of flight", "Abs. reduction in active infections on day of flight", "Rel. reduction in active infections on day of flight", "Total infections during observation", "Abs. reduction in total infections", "Abs. reduction in new infections", "Abs. reduction in active infections", "Abs. reduction in active subclinical infections", "Abs. reduction in active clinical infections", "Rel. reduction in active infections", "Rel. reduction in active subclinical infections", "Rel. reduction in active clinical infections"), 
		ordered=TRUE)
	)
}
categorize_prob_inf <- function(all_results) {
	all_results %>% dplyr::mutate(prob_inf_cat =factor(prob_inf,levels=c(50, 100, 200, 500, 1000, 1500, 2500, 5000) / 1000000,labels=paste(c(5, 10, 20, 50, 100, 150, 250, 500),"per 100,000"),ordered=TRUE))
}
categorize_symptom_screening <- function(all_results) {
	all_results %>% dplyr::mutate(symptom_cat=factor(symptom_screening,levels=c(FALSE, TRUE),labels=c("No symptom screening","With symptom screening"),ordered=TRUE))
}
categorize_rt_multiplier <- function(all_results) {
	all_results %>% mutate(rapid_test_cat=factor(rapid_test_multiplier,levels=c(.6, .75, .9, 1), labels=sprintf("%i%%", round(c(.6, .75, .9, 1) * 100)),ordered=TRUE))
}
categorize_risk_multiplier <- function(all_results) {
	all_results %>% mutate(risk_multi_cat =factor(risk_multiplier,levels=c(1, 2, 4, 10),labels=paste0(c(1, 2, 4, 10), "x"),ordered=TRUE))
}
shift_time_steps <- function(all_results, day_of_flight=70) {
	all_results %>% dplyr::mutate(relative_time=time_step - day_of_flight)
}
add_cume_infections <- function(all_results) {
	all_results %>% dplyr::arrange(testing_type, symptom_screening, risk_multiplier, rapid_test_multiplier, sens_type, prob_inf, prop_subclin, if_threshold, round, rep, time_step ) %>% dplyr::group_by( testing_type, symptom_screening, risk_multiplier, rapid_test_multiplier, sens_type, prob_inf, prop_subclin, if_threshold, round, rep ) %>% dplyr::mutate( cume_n_infection_daily=cumsum(n_infection_daily), cume_w_infection_daily=cumsum(w_infection_daily), cume_n_infected_new=cumsum(n_infected_new) ) %>% dplyr::ungroup()
}
collect_results <- function(scenario_name) {
	all_files <- fs::dir_ls(paste0(outdir,"/intermediate_files/", scenario_name), recurse=TRUE, glob="*.RDS")
	purrr::map_df(.x=all_files, .f=~ readRDS(.x))
}
