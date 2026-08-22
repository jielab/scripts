pacman::p_load(tidyverse, here, fs, config, foreach, doParallel)
scriptdir="/work/sph-huangj/scripts/avia"
source(paste0(scriptdir,"/02.run_sim.f.R"))

outdir="/work/sph-huangj/analysis/avia"
cfig <- config::get(file=paste0(scriptdir,"/02.config.yml"), config="dev")
n_reps <- cfig$n_reps_per_round
n_rounds <- cfig$n_rounds_of_sims
n_pop <- cfig$n_passengers
prob_infs <- cfig$prob_inf / 1000000
prop_subclin <- cfig$prop_subclin
subclin_infectious <- cfig$subclin_infectious
n_cores <- cfig$n_cores
n_burnin <- cfig$n_burnin_days
sens_type <- cfig$sens_type
days_quarantine <- cfig$days_quarantine
risk_multipliers <- cfig$risk_multiplier
rapid_test_multipliers <- cfig$rt_sens_multiplier


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# main analyses
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
day_of_flight <- cfig$day_of_flight
n_outcome_day <- cfig$n_outcome_day
days_quarantine <- cfig$days_quarantine
fs::dir_create(paste0(outdir,"/intermediate_files/", "logs"))
param_grid <- dplyr::bind_rows(
	expand.grid(testing_type="no_testing", prob_inf=prob_infs,round=1:n_rounds,sens_type=NA,risk_multiplier=risk_multipliers,rapid_test_multiplier=NA,symptom_screening=c(TRUE, FALSE), prop_subclin=prop_subclin, stringsAsFactors=FALSE),
	expand.grid(testing_type="perfect_testing", prob_inf=prob_infs,round=1:20,sens_type=NA,risk_multiplier=risk_multipliers,rapid_test_multiplier=NA,symptom_screening=c(TRUE, FALSE),prop_subclin=prop_subclin, stringsAsFactors=FALSE),
	expand.grid(testing_type=c("rapid_test_same_day","rapid_same_day_5_day_quarantine_pcr"),prob_inf=prob_infs,sens_type=sens_type,risk_multiplier=risk_multipliers,rapid_test_multiplier=rapid_test_multipliers,symptom_screening=c(TRUE, FALSE),round=1:n_rounds,prop_subclin=prop_subclin, stringsAsFactors=FALSE),
	expand.grid(testing_type=c("pcr_five_days_after","pcr_three_days_before","pcr_three_days_before_5_day_quarantine_pcr"),prob_inf=prob_infs,sens_type=sens_type,risk_multiplier=risk_multipliers,rapid_test_multiplier=NA,symptom_screening=c(TRUE, FALSE),round=1:n_rounds,prop_subclin=prop_subclin, stringsAsFactors=FALSE)
)
param_grid <- param_grid[with(param_grid,!file.exists(return_sim_file_name(scenario_name=testing_type,symptom_screening=symptom_screening,prob_inf=prob_inf,prop_subclin=prop_subclin, risk_multiplier=risk_multiplier,sens_type=sens_type, rapid_test_multiplier=rapid_test_multiplier,round=round))),]
doParallel::registerDoParallel(cores=n_cores)
foreach::foreach(i=sample(1:nrow(param_grid))) %dopar% {
	print(paste("DO:", i))
	run_sim_pop( # 🏮
		param_grid$testing_type[i], param_grid$prob_inf[i], param_grid$sens_type[i], param_grid$risk_multiplier[i],
		param_grid$rapid_test_multiplier[i], param_grid$symptom_screening[i], param_grid$round[i], n_reps,
		param_grid$prop_subclin[i], subclin_infectious, n_burnin, day_of_flight, n_outcome_day, days_quarantine
	)
	print(paste("DONE:", i))
}
doParallel::stopImplicitCluster()
closeAllConnections()


