pacman::p_load(tidyverse, here, fs, config, foreach, doParallel)
scriptdir="/work/sph-huangj/scripts/avia"
source(paste0(scriptdir,"/02.sum_sim.f.R"))
outdir="/work/sph-huangj/analysis/avia"

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# munge_intermediate_files
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
fs::dir_create(paste0(outdir,"/data_raw"))
testing_scenarios <- basename(fs::dir_ls(paste0(outdir,"/intermediate_files"), type="directory", regexp="pcr|rapid|no_testing|perfect"))
doParallel::registerDoParallel()
all_results <- foreach::foreach(s=testing_scenarios) %dopar% {
	temp_x <- collect_results(s) %>% filter(time_step >= 67,!is.na(if_threshold)) %>%
		shift_time_steps(day_of_flight=70) %>% add_cume_infections() %>%
		categorize_testing_types() %>% categorize_prob_inf() %>% categorize_sens_type() %>% categorize_prop_subclin() %>% 
		categorize_symptom_screening() %>% categorize_rt_multiplier() %>% categorize_risk_multiplier() %>%
		mutate(sim_id=sprintf("%03d.%02d", round, rep)) %>%
		dplyr::select(testing_type,testing_cat,prob_inf,prob_inf_cat,sens_type,sens_cat,prop_subclin,prop_subclin_cat,symptom_screening,symptom_cat,risk_multiplier,risk_multi_cat,rapid_test_multiplier,rapid_test_cat,if_threshold,time_step,round,rep,sim_id,dplyr::everything()) %>%
		arrange(testing_type,testing_cat,prob_inf,prob_inf_cat, prop_subclin, prop_subclin_cat,sens_type,symptom_screening,symptom_cat,risk_multiplier,risk_multi_cat,rapid_test_multiplier,rapid_test_cat,if_threshold,round, rep,time_step)
	if (!tibble::has_name(temp_x, "n_test_false_pos") & !tibble::has_name(temp_x, "n_test_true_pos")) {
		temp_x <- temp_x %>% dplyr::mutate(n_test_false_pos=NA,n_test_true_pos=NA)
	}
	s_name <- paste0(outdir,"/data_raw/", sprintf("raw_simulations_%s.RDS", s))
	saveRDS(temp_x, s_name, compress="xz")
}
doParallel::stopImplicitCluster()


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# summarize_infection_quantities
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
doParallel::registerDoParallel()
temp_holder <- foreach::foreach(s=testing_scenarios) %dopar% {
	temp_x <- readRDS(here::here("data_raw", sprintf("raw_simulations_%s.RDS", s)))
	## Calculate number of active infections (clinical)
	temp_x <- temp_x %>% dplyr::mutate(n_active_infection_clin=n_active_infection - n_active_infection_subclin)
	## Join with null model (no testing, no symptom screening)
	temp_x <- temp_x %>% dplyr::left_join(no_symp)
	## Absolute differences
	temp_x <- temp_x %>% dplyr::mutate( 
		abs_n_infected_all=n_infected_all - null_n_infected_all, abs_n_infected_new=n_infected_new - null_n_infected_new, 
		abs_n_active_infection=n_active_infection - null_n_active_infection, abs_n_active_infection_subclin=n_active_infection_subclin - null_n_active_infection_subclin, 
		abs_n_active_infection_clin=n_active_infection_clin - null_n_active_infection_clin, abs_n_infection_daily=n_infection_daily - null_n_infection_daily, 
		abs_w_infection_daily=w_infection_daily - null_w_infection_daily, 
		abs_cume_n_infected_new=cume_n_infected_new - null_cume_n_infected_new, abs_cume_n_infection_daily=cume_n_infection_daily - null_cume_n_infection_daily, 
		abs_cume_w_infection_daily=cume_w_infection_daily - null_cume_w_infection_daily)
	## Relative differences
	temp_x <- temp_x %>% dplyr::mutate( 
	rel_n_infected_all=(null_n_infected_all - n_infected_all) / null_n_infected_all, 
	rel_n_infected_new=(null_n_infected_new - n_infected_new) / null_n_infected_new, 
	rel_n_active_infection=(null_n_active_infection - n_active_infection) / null_n_active_infection, 
	rel_n_active_infection_subclin=( null_n_active_infection_subclin - n_active_infection_subclin ) / null_n_active_infection_subclin, 
	rel_n_active_infection_clin=(null_n_active_infection_clin - n_active_infection_clin) / null_n_active_infection_clin, 
	rel_n_infection_daily=(null_n_infection_daily - n_infection_daily) / null_n_infection_daily, 
	rel_w_infection_daily=(null_w_infection_daily - w_infection_daily) / null_w_infection_daily, 
	rel_cume_n_infected_new=(null_cume_n_infected_new - cume_n_infected_new) / null_cume_n_infected_new, 
	rel_cume_n_infection_daily=(null_cume_n_infection_daily - cume_n_infection_daily) / null_cume_n_infection_daily, 
	rel_cume_w_infection_daily=(null_cume_w_infection_daily - cume_w_infection_daily) / null_cume_w_infection_daily)
	## Pivot wider
	temp_x_wide <- temp_x %>% pivot_wider( 
		id_cols=c(testing_type,testing_cat,prob_inf,prob_inf_cat,sens_type,sens_cat,prop_subclin,prop_subclin_cat,risk_multi_cat,risk_multiplier,rapid_test_multiplier,rapid_test_cat,if_threshold,time_step,round,rep,sim_id ), 
		names_from=symptom_screening, 
		values_from=c(n_infected_all,n_infected_new,n_active_infection,n_active_infection_clin,n_infection_daily,w_infection_daily,cume_n_infected_new,cume_n_infection_daily,cume_w_infection_daily,abs_n_infected_all:rel_cume_w_infection_daily )
	)
	temp_x_list <- vector("list", length=NROW(adherence_level))
	for (i in 1:NROW(adherence_level)) {
		a <- adherence_level[i]
		temp_x_list[[i]] <- temp_x_wide %>% 
		dplyr::transmute(
		testing_type,testing_cat,prob_inf,prob_inf_cat,sens_type,sens_cat,
		prop_subclin,prop_subclin_cat,risk_multiplier,risk_multi_cat,
		rapid_test_multiplier,rapid_test_cat,if_threshold,time_step,round,rep,sim_id,
		symptom_adherence=a,
		n_infected_all=(n_infected_all_TRUE * a) + (n_infected_all_FALSE * (1 - a)),
		n_infected_new=(n_infected_new_TRUE * a) + (n_infected_new_FALSE * (1 - a)),
		n_active_infection=(n_active_infection_TRUE * a) + (n_active_infection_FALSE * (1 - a)),
		n_active_infection_clin=(n_active_infection_clin_TRUE * a) + (n_active_infection_clin_FALSE * (1 - a)),
		n_infection_daily=(n_infection_daily_TRUE * a) + (n_infection_daily_FALSE * (1 - a)),
		w_infection_daily=(w_infection_daily_TRUE * a) + (w_infection_daily_FALSE * (1 - a)),
		cume_n_infected_new=(cume_n_infected_new_TRUE * a) + (cume_n_infected_new_FALSE * (1 - a)),
		cume_n_infection_daily=(cume_n_infection_daily_TRUE * a) + (cume_n_infection_daily_FALSE * (1 - a)),
		cume_w_infection_daily=(cume_w_infection_daily_TRUE * a) + (cume_w_infection_daily_FALSE * (1 - a)),
		abs_n_infected_all=(abs_n_infected_all_TRUE * a) + (abs_n_infected_all_FALSE * (1 - a)),
		abs_n_infected_new=(abs_n_infected_new_TRUE * a) + (abs_n_infected_new_FALSE * (1 - a)),
		abs_n_active_infection=(abs_n_active_infection_TRUE * a) + (abs_n_active_infection_FALSE * (1 - a)),
		abs_n_active_infection_subclin=(abs_n_active_infection_subclin_TRUE * a) + (abs_n_active_infection_subclin_FALSE * (1 - a)),
		abs_n_active_infection_clin=(abs_n_active_infection_clin_TRUE * a) + (abs_n_active_infection_clin_FALSE * (1 - a)),
		abs_n_infection_daily=(abs_n_infection_daily_TRUE * a) + (abs_n_infection_daily_FALSE * (1 - a)),
		abs_w_infection_daily=(abs_w_infection_daily_TRUE * a) + (abs_w_infection_daily_FALSE * (1 - a)),
		abs_cume_n_infected_new=(abs_cume_n_infected_new_TRUE * a) + (abs_cume_n_infected_new_FALSE * (1 - a)),
		abs_cume_n_infection_daily=(abs_cume_n_infection_daily_TRUE * a) + (abs_cume_n_infection_daily_FALSE * (1 - a)),
		abs_cume_w_infection_daily=(abs_cume_w_infection_daily_TRUE * a) + (abs_cume_w_infection_daily_FALSE * (1 - a)),
		rel_n_infected_all=(rel_n_infected_all_TRUE * a) + (rel_n_infected_all_FALSE * (1 - a)),
		rel_n_infected_new=(rel_n_infected_new_TRUE * a) + (rel_n_infected_new_FALSE * (1 - a)),
		rel_n_active_infection=(rel_n_active_infection_TRUE * a) + (rel_n_active_infection_FALSE * (1 - a)),
		rel_n_active_infection_subclin=(rel_n_active_infection_subclin_TRUE * a) + (rel_n_active_infection_subclin_FALSE * (1 - a)),
		rel_n_active_infection_clin=(rel_n_active_infection_clin_TRUE * a) + (rel_n_active_infection_clin_FALSE * (1 - a)),
		rel_n_infection_daily=(rel_n_infection_daily_TRUE * a) + (rel_n_infection_daily_FALSE * (1 - a)),
		rel_w_infection_daily=(rel_w_infection_daily_TRUE * a) + (rel_w_infection_daily_FALSE * (1 - a)),
		rel_cume_n_infected_new=(rel_cume_n_infected_new_TRUE * a) + (rel_cume_n_infected_new_FALSE * (1 - a)),
		rel_cume_n_infection_daily=(rel_cume_n_infection_daily_TRUE * a) + (rel_cume_n_infection_daily_FALSE * (1 - a)),
		rel_cume_w_infection_daily=(rel_cume_w_infection_daily_TRUE * a) + (rel_cume_w_infection_daily_FALSE * (1 - a)) )
	}
	saveRDS(dplyr::bind_rows(temp_x_list), here::here("data_raw", sprintf("processed_simulations_wide_%s.RDS", s) ), compress="xz")
}
doParallel::stopImplicitCluster()
# Summarize all results that do not incorporate quarantine
doParallel::registerDoParallel()
results_no_quarantine <-
	foreach::foreach(s=testing_scenarios[!grepl("quarantine", testing_scenarios)]) %dopar% {
		temp_x <- readRDS(here::here( "data_raw",sprintf("processed_simulations_wide_%s.RDS", s) ))
		temp_x <- temp_x %>% dplyr::group_by( testing_type,testing_cat,prob_inf,prob_inf_cat,sens_type,sens_cat,prop_subclin,prop_subclin_cat,symptom_adherence,risk_multiplier,risk_multi_cat,rapid_test_multiplier,rapid_test_cat,if_threshold,time_step )
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
		) %>% dplyr::ungroup() %>% dplyr::mutate(quarantine_adherence=0)
	}
doParallel::stopImplicitCluster()
## Now like above, we want to take a weighted average using different weights
## to estimate different levels of adherence to quarantining.
quarantine_comparisons <- dplyr::bind_rows(
	expand.grid(base_case="rapid_test_same_day",comparison_case=c( "rapid_same_day_5_day_quarantine_pcr", "rapid_same_day_7_day_quarantine_pcr", "rapid_same_day_14_day_quarantine_pcr"),stringsAsFactors=FALSE),
	expand.grid(base_case="pcr_three_days_before",comparison_case=c( "pcr_three_days_before_5_day_quarantine_pcr", "pcr_three_days_before_7_day_quarantine_pcr", "pcr_three_days_before_14_day_quarantine_pcr"),stringsAsFactors=FALSE),
) %>%
	tibble::add_case(base_case="pcr_five_days_before", comparison_case="pcr_five_days_before_5_day_quarantine_pcr") %>%
	tibble::add_case(base_case="pcr_seven_days_before", comparison_case="pcr_seven_days_before_5_day_quarantine_pcr") %>%
	tibble::add_case(base_case="pcr_two_days_before", comparison_case="pcr_two_days_before_5_day_quarantine_pcr") %>%
	tibble::add_case(base_case="pcr_five_days_after",comparison_case="5_day_quarantine_pcr_five_days_after")
doParallel::registerDoParallel()
results_quarantine <- foreach::foreach(i=1:NROW(quarantine_comparisons)) %dopar% {
	b <- quarantine_comparisons$base_case[i]
	comp <- quarantine_comparisons$comparison_case[i]
	base_x <- readRDS(here::here("data_raw",sprintf("processed_simulations_wide_%s.RDS", b) ))
	comp_x <- readRDS(here::here("data_raw",sprintf("processed_simulations_wide_%s.RDS", comp) ))
	joined_x <- dplyr::left_join( comp_x, base_x %>%  dplyr::select(-testing_type, -testing_cat), by=c( "prob_inf","prob_inf_cat","sens_type","sens_cat","prop_subclin","prop_subclin_cat","risk_multiplier","risk_multi_cat","rapid_test_multiplier","rapid_test_cat","if_threshold","time_step","round","rep","sim_id","symptom_adherence" ))
	temp_x_list <- vector("list", length=NROW(adherence_level))
	for (i in 1:NROW(adherence_level)) { 
		a <- adherence_level[i]
		temp_x_list[[i]] <- joined_x %>% 
		dplyr::transmute(testing_type,testing_cat,prob_inf,prob_inf_cat,sens_type,sens_cat,prop_subclin,prop_subclin_cat,risk_multiplier,risk_multi_cat,rapid_test_multiplier,rapid_test_cat,if_threshold,time_step,
		round,rep,sim_id,symptom_adherence,
		quarantine_adherence=a,
		n_infected_all=(n_infected_all.x * a) + (n_infected_all.y * (1 - a)),
		n_infected_new=(n_infected_new.x * a) + (n_infected_new.y * (1 - a)),
		n_active_infection=(n_active_infection.x * a) + (n_active_infection.y * (1 - a)),
		n_active_infection_clin=(n_active_infection_clin.x * a) + (n_active_infection_clin.y * (1 - a)),
		n_infection_daily=(n_infection_daily.x * a) + (n_infection_daily.y * (1 - a)),
		w_infection_daily=(w_infection_daily.x * a) + (w_infection_daily.y * (1 - a)),
		cume_n_infected_new=(cume_n_infected_new.x * a) + (cume_n_infected_new.y * (1 - a)),
		cume_n_infection_daily=(cume_n_infection_daily.x * a) + (cume_n_infection_daily.y * (1 - a)),
		cume_w_infection_daily=(cume_w_infection_daily.x * a) + (cume_w_infection_daily.y * (1 - a)),
		abs_n_infected_all=(abs_n_infected_all.x * a) + (abs_n_infected_all.y * (1 - a)),
		abs_n_infected_new=(abs_n_infected_new.x * a) + (abs_n_infected_new.y * (1 - a)),
		abs_n_active_infection=(abs_n_active_infection.x * a) + (abs_n_active_infection.y * (1 - a)),
		abs_n_active_infection_subclin=(abs_n_active_infection_subclin.x * a) + (abs_n_active_infection_subclin.y * (1 - a)),
		abs_n_active_infection_clin=(abs_n_active_infection_clin.x * a) + (abs_n_active_infection_clin.y * (1 - a)),
		abs_n_infection_daily=(abs_n_infection_daily.x * a) + (abs_n_infection_daily.y * (1 - a)),
		abs_w_infection_daily=(abs_w_infection_daily.x * a) + (abs_w_infection_daily.y * (1 - a)),
		abs_cume_n_infected_new=(abs_cume_n_infected_new.x * a) + (abs_cume_n_infected_new.y * (1 - a)),
		abs_cume_n_infection_daily=(abs_cume_n_infection_daily.x * a) + (abs_cume_n_infection_daily.y * (1 - a)),
		abs_cume_w_infection_daily=(abs_cume_w_infection_daily.x * a) + (abs_cume_w_infection_daily.y * (1 - a)),
		rel_n_infected_all=(rel_n_infected_all.x * a) + (rel_n_infected_all.y * (1 - a)),
		rel_n_infected_new=(rel_n_infected_new.x * a) + (rel_n_infected_new.y * (1 - a)),
		rel_n_active_infection=(rel_n_active_infection.x * a) + (rel_n_active_infection.y * (1 - a)),
		rel_n_active_infection_subclin=(rel_n_active_infection_subclin.x * a) + (rel_n_active_infection_subclin.y * (1 - a)),
		rel_n_active_infection_clin=(rel_n_active_infection_clin.x * a) + (rel_n_active_infection_clin.y * (1 - a)),
		rel_n_infection_daily=(rel_n_infection_daily.x * a) + (rel_n_infection_daily.y * (1 - a)),
		rel_w_infection_daily=(rel_w_infection_daily.x * a) + (rel_w_infection_daily.y * (1 - a)),
		rel_cume_n_infected_new=(rel_cume_n_infected_new.x * a) + (rel_cume_n_infected_new.y * (1 - a)),
		rel_cume_n_infection_daily=(rel_cume_n_infection_daily.x * a) + (rel_cume_n_infection_daily.y * (1 - a)),
		rel_cume_w_infection_daily=(rel_cume_w_infection_daily.x * a) + (rel_cume_w_infection_daily.y * (1 - a))) 
	}
	rm(base_x, comp_x); rm(joined_x); invisible(gc(FALSE))
	## Group these up and then use summarize function
	temp_x <- temp_x_list %>% dplyr::bind_rows() %>% dplyr::group_by( testing_type,testing_cat,prob_inf,prob_inf_cat,sens_type,sens_cat,prop_subclin,prop_subclin_cat,symptom_adherence,quarantine_adherence,risk_multiplier,risk_multi_cat,rapid_test_multiplier,rapid_test_cat,if_threshold,time_step )
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
saveRDS(summarized_results,here::here("data", "summarized_results.RDS"),compress="xz")


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# summarize_testing_quantites
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
adherence_level <- seq(0, 1, .2) ## Different levels of adherence to symptom screening
testing_dict <- list(
	no_testing=70,
	pcr_two_days_before=68,
	pcr_two_days_before_5_day_quarantine_pcr=c(70, 75),
	pcr_three_days_before=68,
	pcr_three_days_before_5_day_quarantine_pcr=c(70, 75),
	pcr_three_days_before_7_day_quarantine_pcr=c(70, 77),
	pcr_three_days_before_14_day_quarantine_pcr=c(70, 84),
	pcr_five_days_before=66, 
	pcr_five_days_before_5_day_quarantine_pcr=c(70, 75),
	pcr_seven_days_before=64, 
	pcr_seven_days_before_5_day_quarantine_pcr=c(70, 75),
	rapid_test_same_day=70,
	rapid_same_day_5_day_quarantine_pcr=c(70, 75),
	rapid_same_day_7_day_quarantine_pcr=c(70, 77),
	rapid_same_day_14_day_quarantine_pcr=c(70, 84),
	pcr_five_days_after=75,
	"5_day_quarantine_pcr_five_days_after"=75
)
## Calculate testing quantities ----
doParallel::registerDoParallel(cores=n_cores)
testing_results <- foreach::foreach(i=1:NROW(testing_dict)) %dopar% {
	testing_type <- names(testing_dict)[i]
	time_steps <- testing_dict[[testing_type]]
	temp_x <- readRDS(here::here("data_raw", sprintf("raw_simulations_%s.RDS",testing_type))) %>% 
		dplyr::group_by(
			testing_type, testing_cat, prob_inf, prob_inf_cat, sens_type, sens_cat, prop_subclin, 
			prop_subclin_cat, symptom_screening, risk_multiplier, risk_multi_cat, rapid_test_multiplier, rapid_test_cat, if_threshold, round, rep, sim_id
		) 
	## Get day of flight infections and total infections
	temp_x_infs <- dplyr::left_join(
		temp_x %>% dplyr::filter(time_step == 70) %>% dplyr::select( n_infected_day_of_flight=n_infected_all,n_active_infected_day_of_flight=n_active_infection,n_active_infected_subclin_day_of_flight=n_active_infection_subclin ),
		temp_x  %>% dplyr::summarize(n_total_infections=n_susceptible[time_step == 84] -				  n_susceptible[time_step == 67] +				  n_infected_all[time_step == 67]) %>% dplyr::ungroup() %>% dplyr::distinct()
	)
	## Get testing statistics based on the scenario dictionary because some tests are *before* t-3 (when we started counting cumulative infections), we need to get the raw files instead of processed files. 
	temp_x_test <- collect_results(testing_type) %>% 
		dplyr::group_by( testing_type, prob_inf, sens_type, prop_subclin, symptom_screening, risk_multiplier, rapid_test_multiplier, if_threshold, round, rep) %>% 
		dplyr::filter(time_step %in% time_steps) 
	if (!tibble::has_name(temp_x_test, "n_test_false_pos")) {
		temp_x_test <- temp_x_test %>%  dplyr::mutate(n_test_false_pos=NA_integer_, n_test_true_pos=NA_integer_)
	}
	temp_x_test <- temp_x_test %>% 
		dplyr::select(time_step, n_test_false_pos, n_test_true_pos) %>% 
		dplyr::summarize( 
			n_test_false_pos_first=dplyr::case_when( 
				testing_type == "no_testing"  ~ NA_integer_,
				dplyr::n_distinct(time_steps) == 1 ~ as.integer(max(n_test_false_pos)),
				dplyr::n_distinct(time_steps) == 2 ~ as.integer(n_test_false_pos[time_step == min(time_steps)]) ), n_test_false_pos_second=dplyr::case_when( testing_type == "no_testing"  ~ NA_integer_,dplyr::n_distinct(time_steps) == 1 ~ NA_integer_,dplyr::n_distinct(time_steps) == 2 ~ as.integer(n_test_false_pos[time_step == max(time_steps)]) ), 
			n_test_true_pos_first=dplyr::case_when(
				testing_type == "no_testing"  ~ NA_integer_,
				dplyr::n_distinct(time_steps) == 1 ~ as.integer(max(n_test_true_pos)),
				dplyr::n_distinct(time_steps) == 2 ~ as.integer(n_test_true_pos[time_step == min(time_steps)]) 
			), 
			n_test_true_pos_second=dplyr::case_when(
				testing_type == "no_testing"  ~ NA_integer_,
				dplyr::n_distinct(time_steps) == 1 ~ NA_integer_,
				dplyr::n_distinct(time_steps) == 2 ~ as.integer(n_test_true_pos[time_step == max(time_steps)]) 
			)
		) %>% ungroup() %>% dplyr::distinct() %>% dplyr::rowwise() %>% 
		mutate(
			n_test_false_pos=dplyr::case_when(testing_type == "no_testing"  ~ NA_integer_,TRUE ~ sum(n_test_false_pos_first, n_test_false_pos_second, na.rm=TRUE) ), 
			n_test_true_pos =dplyr::case_when(testing_type == "no_testing"  ~ NA_integer_,TRUE ~ sum(n_test_true_pos_first, n_test_true_pos_second, na.rm=TRUE) )
		) %>% ungroup()
	rm(temp_x); invisible(gc(FALSE))
	## Pivot wider to we can get "adherence" to symptom screening. This shouldn't change anything but reviewers want it. 
	temp_x_wide <- temp_x_infs %>% dplyr::left_join(temp_x_test) %>% 
	tidyverse::pivot_wider( 
		id_cols=c(testing_type,testing_cat,prob_inf,prob_inf_cat,sens_type,sens_cat,prop_subclin,prop_subclin_cat,risk_multi_cat,risk_multiplier,rapid_test_multiplier,rapid_test_cat,if_threshold,round,rep,sim_id ), names_from=symptom_screening, values_from=n_infected_day_of_flight:n_test_true_pos
	)
	## Loop through different levels of adherence to symptom screening. 
	temp_x_list <- vector("list", length=NROW(adherence_level))
	for (i in 1:NROW(adherence_level)) {
		a <- adherence_level[i]
		temp_x_list[[i]] <- temp_x_wide %>% dplyr::transmute(
		testing_type,testing_cat,prob_inf,prob_inf_cat,sens_type,sens_cat,prop_subclin,prop_subclin_cat,risk_multiplier,risk_multi_cat,rapid_test_multiplier,rapid_test_cat,if_threshold,round,rep,sim_id,symptom_adherence=a,
		n_infected_day_of_flight=round((n_infected_day_of_flight_TRUE * a) + (n_infected_day_of_flight_FALSE * (1 - a)) ),
		n_active_infected_day_of_flight=round((n_active_infected_day_of_flight_TRUE * a) + (n_active_infected_day_of_flight_FALSE * (1 - a)) ),
		n_active_infected_subclin_day_of_flight=round((n_active_infected_subclin_day_of_flight_TRUE * a) + (n_active_infected_subclin_day_of_flight_FALSE * (1 - a)) ),
		n_total_infections=round((n_total_infections_TRUE * a) + (n_total_infections_FALSE * (1 - a))),
		n_test_false_pos_first=round((n_test_false_pos_first_TRUE * a) + (n_test_false_pos_first_FALSE * (1 - a)) ),
		n_test_false_pos_second=round((n_test_false_pos_second_TRUE * a) + (n_test_false_pos_second_FALSE * (1 - a)) ),
		n_test_true_pos_first=round((n_test_true_pos_first_TRUE * a) + (n_test_true_pos_first_FALSE * (1 - a)) ),
		n_test_true_pos_second=round((n_test_true_pos_second_TRUE * a) + (n_test_true_pos_second_FALSE * (1 - a)) ),
		n_test_false_pos=round((n_test_false_pos_TRUE * a) + (n_test_false_pos_FALSE * (1 - a))),
		n_test_true_pos=round((n_test_true_pos_TRUE * a) +	(n_test_true_pos_FALSE * (1 - a))) 
		)
	}
	## Calculate testing quantities we are interested in but don't summarize
	temp_x_list %>% dplyr::bind_rows() %>% dplyr::mutate( any_positive_test=n_test_false_pos + n_test_true_pos, ratio_false_true=n_test_false_pos / n_test_true_pos, ratio_true_false=n_test_true_pos / n_test_false_pos, ppv=n_test_true_pos / (n_test_true_pos + n_test_false_pos), frac_detected=n_test_true_pos / n_infected_day_of_flight, frac_active_detected=n_test_true_pos / n_active_infected_day_of_flight, time_step=84)
}
doParallel::stopImplicitCluster()
closeAllConnections()
testing_results <- dplyr::bind_rows(testing_results)
saveRDS(testing_results, here::here("data", "all_testing_results.RDS"), compress="xz")
# Summarize these results
testing_results <- testing_results %>%
	dplyr::group_by(
		testing_type,testing_cat,prob_inf, prob_inf_cat,sens_type,sens_cat,prop_subclin,prop_subclin_cat,
		symptom_adherence,risk_multiplier,risk_multi_cat,rapid_test_multiplier,rapid_test_cat,if_threshold,time_step
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
saveRDS(testing_summary,here::here("data", "summarized_testing_results.RDS"),compress="xz")

