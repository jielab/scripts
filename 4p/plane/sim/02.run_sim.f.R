get_seed_alpha <- function(x) {
	hexval <- paste0("0x", digest::digest(x, "crc32"))
	intval <- utils::type.convert(hexval) %% .Machine$integer.max
	return(intval)
}
sens_by_day_data <- function() {
	structure( list(lower=c( 0.38478627031159, 0.57516731208943, 0.604228529934369, 0.615523437684205, 0.634889417483735, 0.658734818554505, 0.683467990118063, 0.704796690799587, 0.715626316838783, 0.708161673878986, 0.677331542210053, 0.628960600707915, 0.57159750289502, 0.513790902293819, 0.46408945242676, 0.427816321164893, 0.397392733773658, 0.362014429867096, 0.318503776983481, 0.294189654358031, 0.268771607253158, 0.243947269839712, 0.219122932426265, 0.194298595012819, 0.169474257599372), median=c( 0.507191235810699, 0.661649603931173, 0.734764923759696, 0.774741063381784, 0.79627353568264, 0.805680764747389, 0.809281174661157, 0.811418403300466, 0.808536943707438, 0.795106502715592, 0.76729880640491, 0.72809365784123, 0.682172879336852, 0.634218293204075, 0.588911721755199, 0.549417948006275, 0.512833597788347, 0.474738257636212, 0.434289189502867, 0.404954357014122, 0.377078552659996, 0.349468986554785, 0.321859420449574, 0.294249854344363, 0.266640288239152 ), upper=c( 0.630862363616549, 0.746928293064444, 0.836882533155805, 0.888289249824968, 0.907958695940319, 0.907512646103758, 0.898572874917185, 0.89014329958142, 0.880756407692975, 0.866326829447286, 0.84380027533739, 0.814246777046741, 0.779767446556394, 0.742463395847405, 0.704435736900831, 0.667378678469655, 0.631358816394583, 0.596035843288247, 0.561457791759062, 0.52922605439859, 0.499840993181371, 0.47051893746343, 0.441196881745488, 0.411874826027547, 0.382552770309606) ), row.names=c(NA, -25L), class="data.frame")
}
return_sim_state_file <- function(round, rep, prob_inf, prop_subclin) {
	sprintf("%s/round_%02d/round_%02d-rep_%03d-prob_inf_%03d-asx_%02d.rds", paste0(outdir,"/intermediate_files/", "simulation_states"),round, round,rep,prob_inf * 1000000,prop_subclin * 100)
}
return_sim_file_name <- function(scenario_name,symptom_screening,prob_inf,prop_subclin,risk_multiplier,rapid_test_multiplier,sens_type,round) {
	sprintf("%s/%s/%03d-%s-prob_inf_%03d-asx_%02d-risk_multi_%i-sens_type_%s-rt_multi_%i-results.RDS",
	here::here("intermediate_files",scenario_name), ifelse(symptom_screening, "with_screening", "no_screening"),
	round,scenario_name,prob_inf * 1000000,prop_subclin * 100,risk_multiplier,sens_type,rapid_test_multiplier * 100)
}
return_test_sensitivity <- function(sens_type) {
	sens_by_day <- sens_by_day_data()
	if (is.na(sens_type)) {sens <- NA
	} else if (sens_type %in% names(sens_by_day)) {sens <- sens_by_day %>% dplyr::pull(sens_type)
	} else if (sens_type == "perfect") {sens <- rep(1, NROW(sens_by_day))
	} else if (sens_type == "random") {sens <- apply(sens_by_day, 1, function(x) {truncnorm::rtruncnorm( n=1,a=x[["lower"]],mean=x[["median"]],b=x[["upper"]],sd=(x[["upper"]] - x[["median"]]) / 2)})
	} else {stop(paste(c("Not a valid sens_type. Must be one of:", names(sens_by_day), "perfect", "random"),collapse=" "))}
	sens
}
draw_incubation <- function(n, prop_pre_sx=1/2) {
	ptri <- function(x, a, b, c) {ifelse(x > a & x < b,ifelse(x <= c,((x-a)^2)/((b-a)*(c-a)), 1-((b-x)^2)/((b-a)*(b-c))),ifelse(x <= a,0,1))}
	qtri <- function(p, a, b, c) {ifelse(p > 0 & p < 1, ifelse(p <= ptri(c, a, b, c),a+sqrt((a^2-a*b-(a-b)*c)*p),b-sqrt(b^2-a*b+(a-b)*c+(a*b-b^2-(a-b)*c)*p)), NA)}
	rtri <- function(n, a, b, c) {qtri(stats::runif(n, min=0, max=1), a, b, c)}
	## This is a triangle distribution with median/mode at 3 (mean ~3.3) and bounds at [2, 5] after rounding.
	pre_sx_not_infectious <- round((rtri(n, 1.5, 5.49, 3)))
	## This is a truncated normal with some proportion as pre-symptomatic
	pre_sx_infectious <- round(draw_duration_of_infectiousness(n) * prop_pre_sx)
	## We return the sum of both days and keep the number of infectious pre-symptomatic days as the names
	res <- pre_sx_infectious + pre_sx_not_infectious
	names(res) <- pre_sx_infectious
	res
}
draw_symptomatic <- function(n, days_test_pos_noninf=14, prop_pre_sx=1/2) {
	sx_infectious <- round(draw_duration_of_infectiousness(n) * prop_pre_sx)
	res <- sx_infectious + days_test_pos_noninf
	names(res) <- sx_infectious
	res
}
draw_duration_of_infectiousness <- function(n, mean=5, a=1, b=Inf, sd=1) {
	round(truncnorm::rtruncnorm( n=n, a=a, b=b, mean=mean, sd=sd))
}
draw_infectiousness_weights <- function(days_incubation, days_symptomatic) {
	n_days_pre_sx_inf <- as.numeric(names(days_incubation))
	n_days_no_sx_no_inf <- days_incubation - n_days_pre_sx_inf
	n_days_sx_inf <- as.numeric(names(days_symptomatic))
	n_days_sx_no_inf <- days_symptomatic - n_days_sx_inf
	d <- c( rep(0, n_days_no_sx_no_inf), rep(1, n_days_pre_sx_inf), rep(1, n_days_sx_inf), rep(0, n_days_sx_no_inf))
	names(d) <- c(paste0("e", 1:days_incubation), paste0("l", 1:days_symptomatic))
	return(d)
}
draw_specificity <- function(n=1,lower=.995,upper=1,mean=.998,sd=.001) {truncnorm::rtruncnorm( n=n, a=lower, b=upper, mean=mean, sd=sd)}
create_sim_pop <- function(n_pop, prop_subclin, days_incubation, days_symptomatic) {
	res <- tibble::tibble( ID=as.character(1:n_pop), InfectionType=stats::rbinom(n_pop, 1, 1 - prop_subclin), State=0, DaysIncubation=days_incubation, DaysPreSxInf=as.numeric(names(days_incubation)), DaysSxInf=as.numeric(names(days_symptomatic)), DaysSymptomatic=days_symptomatic, NewInfection=0, DayOfInfection=0, DayInObs=0, DayOfDetection=0)
	res <- data.table::as.data.table(res)
	data.table::setkey(res, ID)
	rownames(res) <- res$ID
	res
}
burnin_sim <- function(sim_pop, n_iter=50, ...) {
	for (i in 1:n_iter) { sim_pop <- increment_sim(sim_pop=sim_pop, t_step=i, ...)}
	return(sim_pop)
}
day_in_infection <- function(x, t) {
    cols <- c("State", "DayOfInfection")
    x <- as.numeric(x[cols])
    names(x) <- cols
    if (x["State"] %in% c(1, 2)) { return(t - x["DayOfInfection"] + 1) }
    return(NA)
}
test_pos <- function(sim_pop_pos, t, p_sens, days_incubation, days_symptomatic, n_delay=0, min_days=5) {
	p_i <- ceiling(apply(sim_pop_pos, 1, function(y) day_in_infection(y, t)))
	# can only test positive X days prior to symptom onset
	toffset <- days_incubation - min_days
	p_i <- p_i - toffset 
	sTested_tp <- sapply(p_i, function(y) ifelse(y > 0 & y <= length(p_sens), stats::rbinom(1, 1, p_sens[y]), 0))	
	sim_pop_pos[sTested_tp == 1, DayOfDetection := t + n_delay]
	pos_ids <- sim_pop_pos[sTested_tp == 1, ID]
	return(list(sim_pop_pos, pos_ids))
}

increment_sim <- function(sim_pop, t_step, if_weights, days_quarantine, prob_inf) {
	sim_pop[, DayInObs := ifelse(DayInObs <= days_quarantine &  DayInObs > 0, DayInObs + 1, 0)]
	sim_pop[DayOfDetection == t_step, DayInObs := 1]
	sim_pop[, NewInfection := 0]
	# S --> E; Get IDs of susceptible workers
	sWorkers <- sim_pop[State == 0 & DayInObs == 0, ID]
	## Passengers have a static per-day risk of infection
	sim_pop[sWorkers, NewInfection := stats::rbinom(length(sWorkers), 1, prob_inf)]
	## Update states
	sim_pop <- sim_pop %>% dplyr::mutate(State=dplyr::case_when( State == 1 & t_step >= (DaysIncubation + DayOfInfection) ~ 2, State == 2 & t_step >= (DaysIncubation + DaysSymptomatic + DayOfInfection) ~ 3, TRUE ~ State ))
	# Update infectious timelines
	sim_pop[NewInfection == 1, c("State", "DayOfInfection") := list(1, t_step + 1)]
	# Add time step for reference
	sim_pop$time_step <- t_step
	return(sim_pop)
}
unpack_simulation_state <- function(f_name) {
	x <- readRDS(f_name)
	assign("sim_pop", data.table::copy(data.table::setDT(x$sim_pop)), envir=.GlobalEnv)
	assign("p_sens", x$p_sens, envir=.GlobalEnv)
	assign("p_spec", x$p_spec, envir=.GlobalEnv)
	assign("days_incubation", x$days_incubation, envir=.GlobalEnv)
	assign("days_symptomatic", x$days_symptomatic, envir=.GlobalEnv)
	assign("if_weights", x$if_weights, envir=.GlobalEnv)
	assign(".Random.seed", x$rseed, envir=.GlobalEnv)
}
rapid_antigen_testing <- function(sim_pop, t, p_sens, p_spec, rapid_antigen_multiplier, days_incubation, days_symptomatic, n_delay=0, min_days=5) {
	x <- data.table::copy(sim_pop) 
	## Remove everybody who is in quarantine
	wids <- x[DayInObs == 0 & DayOfDetection < t, ID]
	## Test the rest
	tested_pop <- test_pop( x[wids,], t=t, p_sens=p_sens * rapid_antigen_multiplier, p_spec=p_spec, days_incubation=days_incubation, days_symptomatic=days_symptomatic, n_delay=n_delay, min_days=min_days
	)
	## Put new results back
	x[wids,] <- tested_pop[["sim_pop"]]
	## Add columns for the false/true positive flags
	x %>% dplyr::mutate( sTested_tp=ifelse(ID %in% tested_pop[["sTested_tp"]], 1, 0), sTested_fp=ifelse(ID %in% tested_pop[["sTested_fp"]], 1, 0) )
}
pcr_test_subset_simpop <- function(sim_pop, k, n_subsets=2, subset_to_test=1, p_sens, p_spec, days_incubation, days_symptomatic, n_delay=0, min_days=5) {
	test_ids <- return_portion_of_ids(sim_pop, n_subsets, subset_to_test)
	data.table::setDF(sim_pop)
	sim_pop_to_test <- sim_pop %>% dplyr::filter(ID %in% test_ids)
	sim_pop_no_test <- sim_pop %>% dplyr::filter(!(ID %in% test_ids))
	## NOTE: YOU HAVE TO REMOVE THE EXTRA COLUMNS OR DATA.TABLE
	## WON'T LET YOU RUN THE TESTING CODE.
	if (tibble::has_name(sim_pop_to_test, "sTested_tp") | tibble::has_name(sim_pop_to_test, "sTested_fp")) { if (any(!is.na( c( sim_pop_to_test$sTested_tp, sim_pop_to_test$sTested_fp ) ))) { stop("This population has already been tested before.") } else { sim_pop_to_test <- sim_pop_to_test %>% dplyr::select(-sTested_tp, -sTested_fp) }
	}
	## Test the ones we need to test
	sim_pop_to_test <- set_id_as_key(sim_pop_to_test)
	sim_pop_to_test <- pcr_testing(sim_pop_to_test, k, p_sens, p_spec, days_incubation, days_symptomatic)
	data.table::setDF(sim_pop_to_test)
	## Put it back together again and reset the key
	sim_pop <- dplyr::bind_rows(sim_pop_to_test, sim_pop_no_test)
	sim_pop <- set_id_as_key(sim_pop)
	sim_pop
}
clear_old_tests <- function(sim_pop) {
	sim_pop[, c("sTested_fp", "sTested_tp") := NULL]
}
print_sim_start <- function(scenario_name, symptom_screening, prob_inf, sens_type, risk_multiplier, rapid_test_multiplier, prop_subclin, round) {
	print( sprintf( paste0( "Starting (%s): ", "test: %s, ", "round: %s, ", "sens: %s, ", "symptom_testing: %s, ", "prob_inf: %s, ", "prop_subclin: %s, ", "risk_multi: %s, ", "rt_multi: %s" ), Sys.time(), scenario_name, round, sens_type, symptom_screening, prob_inf, prop_subclin, risk_multiplier, rapid_test_multiplier )
	)
}
print_sim_rep <- function(rep, days_incubation, days_symptomatic, prob_inf, prop_subclin) {
	print( sprintf( "Rep %s (%s): days_inc: %s, days_sympt: %s, prob_inf: %s, prop_subclin: %s", rep, Sys.time(), days_incubation, days_symptomatic, prob_inf, prop_subclin )
	)
}
set_id_as_key <- function(res) {
	res <- data.table::as.data.table(res)
	data.table::setkey(res, ID)
	rownames(res) <- res$ID
	res
}
return_portion_of_ids <- function(sim_pop, n_split=2, return_chunk=1) {
	x <- sim_pop$ID
	split(x, sort(1:length(x) %% n_split))[[return_chunk]]
}
summarize_sim_state <- function(sim_pop) {
	sim_pop %>% dplyr::group_by(time_step) %>% dplyr::summarize( n_susceptible=sum(State == 0), n_incubation_sub=sum(State == 1 & InfectionType == 0), n_infected_sub=sum(State == 2 & InfectionType == 0), n_incubation_clin=sum(State == 1 & InfectionType == 1), n_infected_clin=sum(State == 2 & InfectionType == 1), n_infected_all=sum(State %in% 1:2), n_infected_new=sum(NewInfection), n_recovered=sum(State == 3), n_quarantined=sum(DayInObs > 0), cume_n_detected=sum(DayOfDetection > 0), min_day_inf=min(DayOfInfection), max_day_inf=max(DayOfInfection), min_day_detect=min(DayOfDetection), max_day_detect=max(DayOfDetection) ) %>% dplyr::ungroup()
}
summarize_infectiousness <- function(sim_pop, if_weights, day_of_flight, subclin_infectious, min_if_weights=0) {
	res <- tibble::tibble()
	for (w in min_if_weights) {
		x2 <- calculate_daily_infectiousness(sim_pop, if_weights, subclin_infectious, fast_mode=TRUE, min_if_weight=w )  
		if (NROW(x2) > 0) {
			res <- dplyr::bind_rows( res, tibble::tibble(time_step=unique(x2$time_step), if_threshold=w, n_infection_daily=sum(x2$n_infection_daily, na.rm=TRUE), w_infection_daily=sum(x2$w_infection_daily, na.rm=TRUE), n_active_infection=sum(x2$n_infection_daily > 0, na.rm=TRUE), n_active_infection_subclin=sum(x2$n_infection_daily == .5, na.rm=TRUE) ) ) 
		} else { 
			res <- dplyr::bind_rows( res, tibble::tibble(time_step=unique(x2$time_step), if_threshold=w, n_infection_daily=0, w_infection_daily=0, n_active_infection=0, n_active_infection_subclin=0 ) ) 
		} 
	}
	res
}
pcr_testing <- function(sim_pop, t, p_sens, p_spec, days_incubation, days_symptomatic, n_delay=0, min_days=5) {
	## The test_pop function doesn't work if we introduce new columns so
	## we need to check for previous testing columns.
	x <- data.table::copy(sim_pop)
	if (tibble::has_name(x, "sTested_tp") | tibble::has_name(x, "sTested_fp")) { if (any(!is.na(c(x$sTested_tp, x$sTested_fp)))) { stop("This population has already been tested before.") } else { stop("Remove sTested_tp and sTested_fp.") }}
	## Remove everybody who is in quarantine
	wids <- x[DayInObs == 0 & DayOfDetection < t, ID]
	## Test the rest
	tested_pop <- test_pop( x[wids,], t=t, p_sens=p_sens, p_spec=p_spec, days_incubation=days_incubation, days_symptomatic=days_symptomatic, n_delay=n_delay, min_days=min_days
	)
	## Put new results back
	x[wids,] <- tested_pop[["sim_pop"]]
	## Add columns for the false/true positive flags
	x %>% dplyr::mutate( sTested_tp=ifelse(ID %in% tested_pop[["sTested_tp"]], 1, 0), sTested_fp=ifelse(ID %in% tested_pop[["sTested_fp"]], 1, 0) )
}
calculate_daily_infectiousness <- function(sim_pop, if_weights, subclin_infectious, fast_mode=FALSE, min_if_weight=0, keep_cols=FALSE) {
	## Not an infectious day unless meets minimum threshold
	if_binary <- (if_weights > min_if_weight) + 0
	## Reweigh infectiousness limited to just infectious days (at threshold)
	if_weights_x <- ifelse(if_weights <= min_if_weight, 0, if_weights)
	if_weights_x <- if_weights_x / sum(if_weights_x)
	## Prep the dataframe by adding number where they are in infectious
	## trajectory and number of symptomatic days
	x <- sim_pop %>% add_day_since_infection()
	if (fast_mode) { x <- x %>% dplyr::filter(infection_day > 0, State != 3)
	}
	## On this day, are they infectious (0, .5, 1) and if so, how infectious?
	x <- x %>% dplyr::rowwise() %>% 
		mutate(n_infection_daily=ifelse( infection_day > 0 & infection_day <= (DaysIncubation + DaysSymptomatic) & DayInObs == 0 & !is.na(infection_day), if_binary[infection_day], 0 ), w_infection_daily=ifelse( infection_day > 0 & infection_day <= (DaysIncubation + DaysSymptomatic) & DayInObs == 0 & !is.na(infection_day), if_weights_x[infection_day], 0 ) ) %>% dplyr::ungroup() %>% 
		mutate(w_infection_daily=ifelse( InfectionType == 0, subclin_infectious * w_infection_daily, w_infection_daily ) ) %>% 
		mutate(n_infection_daily=ifelse( InfectionType == 0, subclin_infectious * n_infection_daily, n_infection_daily ) )
	if (keep_cols) { x }
	x %>% dplyr::select(-infection_day)
}
add_day_since_infection <-  function(sim_pop) {
	sim_pop %>% dplyr::mutate(infection_day=time_step - ifelse(DayOfInfection > 0, DayOfInfection, NA))
}
test_pop <- function(sim_pop, t, p_sens, p_spec, days_incubation, days_symptomatic, n_delay=0, min_days=5) {
	## Test the infected. Subset to infected IDs
	pids <- sim_pop[State %in% 1:2 & DayOfDetection < DayOfInfection, ID]
	test_pos <- test_pos( sim_pop[pids,], t, p_sens=p_sens, days_incubation=days_incubation, days_symptomatic=days_symptomatic, n_delay=n_delay, min_days=min_days
	)
	## Replace infected with their rest results and save true positives
	sim_pop[pids,] <- test_pos[[1]]
	## Test susceptible population
	test_neg <- test_neg(sim_pop[State == 0,], t, p_spec=p_spec, n_delay=0)
	sim_pop[State == 0,] <- test_neg[[1]]
	## Test recovered
	rids <- sim_pop[State == 3 & DayOfDetection < DayOfInfection, ID]
	test_neg_recovered <- test_neg(sim_pop[rids,], t, p_spec=p_spec, n_delay=0)
	sim_pop[rids,] <- test_neg_recovered[[1]]
	return(list( sim_pop=sim_pop, sTested_tp=test_pos[[2]], sTested_fp=unique(c(test_neg[[2]], test_neg_recovered[[2]]))
	))
}
test_neg <- function(sim_pop_neg, t, p_spec, n_delay=0) {
	sTested_fp <- stats::rbinom(nrow(sim_pop_neg), 1, 1 - p_spec)
	sim_pop_neg[sTested_fp == 1, DayOfDetection := t + n_delay]
	neg_ids <- sim_pop_neg[sTested_fp == 1, ID]
	return(list(sim_pop_neg, neg_ids))
}


run_sim_pop <- function(
	testing_type,prob_inf,sens_type,risk_multiplier, rapid_test_multiplier, symptom_screening,
	round,n_reps, prop_subclin, subclin_infectious,n_burnin,day_of_flight,n_outcome_day, days_quarantine
	) {
sim_file <- return_sim_file_name(scenario_name=testing_type,symptom_screening=symptom_screening,prob_inf=prob_inf,prop_subclin=prop_subclin, risk_multiplier=risk_multiplier,rapid_test_multiplier=rapid_test_multiplier,sens_type=sens_type, round=round)
fs::dir_create(dirname(sim_file))
sink(sprintf("%s/log_%s.txt", paste0(outdir,"/intermediate_files/", "logs"), Sys.getpid()), append=TRUE)
if (!file_exists(sim_file)) {
print_sim_start(testing_type,symptom_screening, prob_inf,sens_type,risk_multiplier,rapid_test_multiplier,prop_subclin, round)
holder <- vector("list", n_reps)
for (r in 1:n_reps) {
	unpack_simulation_state(return_sim_state_file(round, r, prob_inf, prop_subclin))
	p_sens <- return_test_sensitivity(sens_type)
	print_sim_rep(r, days_incubation, days_symptomatic, prob_inf, prop_subclin)
	if (!(testing_type %in% c("no_testing", "perfect_testing", "pcr_three_days_before", "pcr_three_days_before_5_day_quarantine_pcr","rapid_test_same_day", "rapid_same_day_5_day_quarantine_pcr", "pcr_five_days_after"))) { stop("Not a valid testing_type")}
	if (testing_type == "no_testing") {
		## 1. Pre-flight iterations ----
		res <- tibble()
		for (k in n_burnin:(day_of_flight - 1)) {
			sim_pop <- increment_sim(sim_pop=sim_pop, t_step=k,prob_inf=prob_inf,if_weights=if_weights,days_quarantine=days_quarantine)
			res <- bind_rows(res,left_join(summarize_sim_state(sim_pop),summarize_infectiousness(sim_pop,if_weights,day_of_flight,subclin_infectious),by="time_step"))
		}
		## 2. Day of flight iteration ### Day of flight symptom screening (if applicable)
		if (symptom_screening) {sim_pop[DayInObs == 0 & InfectionType == 1 &State == 2, DayOfDetection := day_of_flight]}
		## Increment on the day of flight. Keep this separate because day-of-flight has a risk multiplier
		sim_pop <- increment_sim(sim_pop=sim_pop,t_step=day_of_flight,prob_inf=prob_inf * risk_multiplier, if_weights=if_weights,days_quarantine=days_quarantine)
		res <- bind_rows(res,left_join(summarize_sim_state(sim_pop),summarize_infectiousness(sim_pop,if_weights,day_of_flight,subclin_infectious),by="time_step"))
		## 3. Post-flight iterations ----
		for (k in (day_of_flight + 1):(day_of_flight + n_outcome_day)) {
			## We only want to count infections we could have potentially averted so after our t+3 (last possible testing day across all our scenarios), we stop infecting the population.
			if (k > (day_of_flight + 3)) {
				sim_pop <- increment_sim(sim_pop=sim_pop,t_step=k,prob_inf=0,if_weights=if_weights,days_quarantine=days_quarantine)
			} else {
				sim_pop <- increment_sim(sim_pop=sim_pop, t_step=k,prob_inf=prob_inf,if_weights=if_weights,days_quarantine=days_quarantine)
			}
			res <- bind_rows(res,left_join(summarize_sim_state(sim_pop),summarize_infectiousness(sim_pop,if_weights,day_of_flight,subclin_infectious),by="time_step"))
		}
	}
	if (testing_type == "perfect_testing") {
		res <- tibble()
		for (k in n_burnin:(day_of_flight - 1)) {
			sim_pop <- pcr_test_subset_simpop(sim_pop,k,n_subsets=1, subset_to_test=1, p_sens=return_test_sensitivity("perfect"), p_spec=1,days_incubation,days_symptomatic)
			sim_pop <- increment_sim(sim_pop=sim_pop,t_step=k, prob_inf=prob_inf,if_weights=if_weights,days_quarantine=days_quarantine)
			res <- bind_rows(res, left_join(summarize_sim_state(sim_pop) %>% mutate(n_test_false_pos=sum(sim_pop$sTested_fp, na.rm=TRUE),n_test_true_pos=sum(sim_pop$sTested_tp, na.rm=TRUE)),summarize_infectiousness(sim_pop,if_weights,day_of_flight,subclin_infectious),by="time_step") )
			sim_pop <- clear_old_tests(sim_pop)
		}
		if (symptom_screening) { sim_pop[DayInObs == 0 &InfectionType == 1 &State == 2, DayOfDetection := day_of_flight] }
		sim_pop <- pcr_test_subset_simpop(sim_pop,k,n_subsets=1,subset_to_test=1,p_sens=return_test_sensitivity("perfect"),p_spec=1,days_incubation, days_symptomatic)
		sim_pop <- increment_sim( sim_pop=sim_pop,t_step=day_of_flight, prob_inf=prob_inf * risk_multiplier, if_weights=if_weights,days_quarantine=days_quarantine)
		res <- bind_rows(res,left_join(summarize_sim_state(sim_pop) %>%mutate(n_test_false_pos=sum(sim_pop$sTested_fp, na.rm=TRUE),n_test_true_pos=sum(sim_pop$sTested_tp, na.rm=TRUE)),summarize_infectiousness(sim_pop,if_weights,day_of_flight,subclin_infectious),by="time_step"))
		sim_pop <- clear_old_tests(sim_pop)
		for (k in (day_of_flight + 1):(day_of_flight + n_outcome_day)) {
			sim_pop <- pcr_test_subset_simpop(sim_pop,k,n_subsets=1,subset_to_test=1, p_sens=return_test_sensitivity("perfect"), p_spec=1,days_incubation,days_symptomatic)
			sim_pop <- increment_sim(sim_pop=sim_pop,t_step=k,prob_inf=prob_inf,if_weights=if_weights,days_quarantine=days_quarantine)
			res <- bind_rows(res,left_join(summarize_sim_state(sim_pop) %>% mutate(n_test_false_pos=sum(sim_pop$sTested_fp, na.rm=TRUE),n_test_true_pos=sum(sim_pop$sTested_tp, na.rm=TRUE)),summarize_infectiousness(sim_pop,if_weights,day_of_flight,subclin_infectious),by="time_step"))
			sim_pop <- clear_old_tests(sim_pop)
		}
	}
	if (testing_type == "pcr_three_days_before") {
		test_days <- day_of_flight - 2:3
		res <- tibble()
		for (k in n_burnin:(day_of_flight - 1)) {
			if (k %in% test_days) {sim_pop <- pcr_test_subset_simpop(sim_pop,k,n_subsets=NROW(test_days),subset_to_test=which(test_days == k),p_sens,p_spec,days_incubation,days_symptomatic)}
			sim_pop <- increment_sim(sim_pop=sim_pop,t_step=k,prob_inf=prob_inf,if_weights=if_weights,days_quarantine=days_quarantine)
			if (k %in% test_days) {
				res <- bind_rows(res,left_join(summarize_sim_state(sim_pop) %>% mutate(n_test_false_pos=sum(sim_pop$sTested_fp, na.rm=TRUE),n_test_true_pos=sum(sim_pop$sTested_tp, na.rm=TRUE)),summarize_infectiousness(sim_pop,if_weights,day_of_flight,subclin_infectious),by="time_step"))
			} else {
				res <- bind_rows(res,left_join(summarize_sim_state(sim_pop),summarize_infectiousness(sim_pop,if_weights,day_of_flight,subclin_infectious),by="time_step"))
			}
		}
		if (symptom_screening) {sim_pop[DayInObs == 0 & InfectionType == 1 & State == 2, DayOfDetection := day_of_flight]}
		sim_pop <- increment_sim(sim_pop=sim_pop,t_step=day_of_flight,prob_inf=prob_inf * risk_multiplier,if_weights=if_weights,days_quarantine=days_quarantine)		
		res <- bind_rows(res,left_join(summarize_sim_state(sim_pop),summarize_infectiousness(sim_pop,if_weights,day_of_flight,subclin_infectious),by="time_step"))
		for (k in (day_of_flight + 1):(day_of_flight + n_outcome_day)) {
			if (k > (day_of_flight + 3)) {
				sim_pop <- increment_sim(sim_pop=sim_pop,t_step=k,prob_inf=0,if_weights=if_weights,days_quarantine=days_quarantine)
			} else {
				sim_pop <- increment_sim(sim_pop=sim_pop,t_step=k,prob_inf=prob_inf,if_weights=if_weights,days_quarantine=days_quarantine)
			}
			res <- bind_rows(res,left_join(summarize_sim_state(sim_pop),summarize_infectiousness(sim_pop,if_weights,day_of_flight,subclin_infectious),by="time_step"))
		}
	}
	if (testing_type == "pcr_three_days_before_5_day_quarantine_pcr") {
		test_days <- day_of_flight - 2:3
		res <- tibble()
		for (k in n_burnin:(day_of_flight - 1)) {
			if (k %in% test_days) {
				sim_pop <- pcr_test_subset_simpop(sim_pop,k,n_subsets=NROW(test_days),subset_to_test=which(test_days == k),p_sens,p_spec,days_incubation,days_symptomatic)
			}
			sim_pop <- increment_sim(sim_pop=sim_pop,t_step=k,prob_inf=prob_inf,if_weights=if_weights,days_quarantine=days_quarantine)
			if (k %in% test_days) {
				res <- bind_rows(res,left_join( summarize_sim_state(sim_pop) %>% mutate(n_test_false_pos=sum(sim_pop$sTested_fp, na.rm=TRUE),n_test_true_pos=sum(sim_pop$sTested_tp, na.rm=TRUE)),summarize_infectiousness(sim_pop,if_weights,day_of_flight,subclin_infectious),by="time_step"))
			} else {
				res <- bind_rows(res,left_join(summarize_sim_state(sim_pop),summarize_infectiousness(sim_pop,if_weights,day_of_flight,subclin_infectious),by="time_step"))
			}
		}
		if (symptom_screening) {sim_pop[DayInObs == 0 & InfectionType == 1 & State == 2, DayOfDetection := day_of_flight]}
		sim_pop <- increment_sim(sim_pop=sim_pop,t_step=day_of_flight,prob_inf=prob_inf * risk_multiplier,if_weights=if_weights,days_quarantine=days_quarantine)
		res <- bind_rows(res,left_join(summarize_sim_state(sim_pop) %>% mutate(n_test_false_pos=sum(sim_pop$sTested_fp),n_test_true_pos=sum(sim_pop$sTested_tp)),summarize_infectiousness(sim_pop,if_weights,day_of_flight,subclin_infectious),by="time_step"))
		test_days <- day_of_flight + 5
		sim_pop <- clear_old_tests(sim_pop)
		sim_pop$DayInObs <- 1
		for (k in (day_of_flight + 1):(day_of_flight + n_outcome_day)) {
			if (k %in% test_days) {
				sim_pop <- pcr_test_subset_simpop(sim_pop,k,n_subsets=NROW(test_days),subset_to_test=which(test_days == k),p_sens,p_spec,days_incubation,days_symptomatic)
			}
			sim_pop <- increment_sim(sim_pop=sim_pop,t_step=k,prob_inf=0,if_weights=if_weights,days_quarantine=days_quarantine)
			if (k %in% test_days) {
				res <- bind_rows(res,left_join(summarize_sim_state(sim_pop) %>% mutate(n_test_false_pos=sum(sim_pop$sTested_fp, na.rm=TRUE),n_test_true_pos=sum(sim_pop$sTested_tp, na.rm=TRUE)),summarize_infectiousness(sim_pop,if_weights,day_of_flight,subclin_infectious),by="time_step"))
			} else {
				res <- bind_rows(res,left_join(summarize_sim_state(sim_pop),summarize_infectiousness(sim_pop,if_weights,day_of_flight,subclin_infectious),by="time_step"))
			}
		}
	}
	if (testing_type == "rapid_test_same_day") {
		res <- tibble()
		for (k in n_burnin:(day_of_flight - 1)) {
			sim_pop <- increment_sim(sim_pop=sim_pop,t_step=k,prob_inf=prob_inf,if_weights=if_weights,days_quarantine=days_quarantine)
			res <- bind_rows(res,left_join(summarize_sim_state(sim_pop),summarize_infectiousness(sim_pop,if_weights,day_of_flight,subclin_infectious),by="time_step"))
		}
		if (symptom_screening) {sim_pop[DayInObs == 0 & InfectionType == 1 & State == 2, DayOfDetection := day_of_flight]}
		sim_pop <- rapid_antigen_testing(sim_pop, t=day_of_flight,p_sens,p_spec,rapid_test_multiplier,days_incubation,days_symptomatic)
		sim_pop <- increment_sim(sim_pop=sim_pop,t_step=day_of_flight,prob_inf=prob_inf * risk_multiplier,if_weights=if_weights, days_quarantine=days_quarantine)
		res <- bind_rows(res,left_join(summarize_sim_state(sim_pop) %>% mutate(n_test_false_pos=sum(sim_pop$sTested_fp),n_test_true_pos=sum(sim_pop$sTested_tp)),summarize_infectiousness(sim_pop,if_weights,day_of_flight,subclin_infectious),by="time_step"))
		for (k in (day_of_flight + 1):(day_of_flight + n_outcome_day)) {
			if (k > (day_of_flight + 3)) {
				sim_pop <- increment_sim(sim_pop=sim_pop,t_step=k,prob_inf=0,if_weights=if_weights,days_quarantine=days_quarantine)
			} else {
				sim_pop <- increment_sim(sim_pop=sim_pop,t_step=k,prob_inf=prob_inf,if_weights=if_weights,days_quarantine=days_quarantine)
			}
			res <- bind_rows(res,left_join(summarize_sim_state(sim_pop),summarize_infectiousness(sim_pop,if_weights,day_of_flight,subclin_infectious),by="time_step"))
		}
	}		
	if (testing_type == "rapid_same_day_5_day_quarantine_pcr") {
		res <- tibble()
		for (k in n_burnin:(day_of_flight - 1)) {
			sim_pop <- increment_sim(sim_pop=sim_pop,t_step=k,prob_inf=prob_inf,if_weights=if_weights,days_quarantine=days_quarantine)
			res <- bind_rows(res,left_join(summarize_sim_state(sim_pop),summarize_infectiousness(sim_pop,if_weights,day_of_flight,subclin_infectious),by="time_step"))
		}
		if (symptom_screening) {sim_pop[DayInObs == 0 & InfectionType == 1 & State == 2, DayOfDetection := day_of_flight]}
		sim_pop <- rapid_antigen_testing(sim_pop,t=day_of_flight,p_sens, p_spec,rapid_test_multiplier,days_incubation,days_symptomatic)
		sim_pop <- increment_sim(sim_pop=sim_pop,t_step=day_of_flight, prob_inf=prob_inf * risk_multiplier,if_weights=if_weights,days_quarantine=days_quarantine)
		res <- bind_rows(res,left_join(summarize_sim_state(sim_pop) %>% mutate(n_test_false_pos=sum(sim_pop$sTested_fp),n_test_true_pos=sum(sim_pop$sTested_tp)),summarize_infectiousness(sim_pop,if_weights,day_of_flight,subclin_infectious),by="time_step"))
		test_days <- day_of_flight + 5
		sim_pop <- clear_old_tests(sim_pop)
		sim_pop$DayInObs <- 1
		for (k in (day_of_flight + 1):(day_of_flight + n_outcome_day)) {
			if (k %in% test_days) {
				sim_pop <- pcr_test_subset_simpop(sim_pop,k,n_subsets=NROW(test_days),subset_to_test=which(test_days == k),p_sens,p_spec,days_incubation,days_symptomatic)
			}
			sim_pop <- increment_sim(sim_pop=sim_pop,t_step=k,prob_inf=0,if_weights=if_weights,days_quarantine=days_quarantine)
			if (k %in% test_days) {
				res <- bind_rows(res,left_join(summarize_sim_state(sim_pop) %>% mutate(n_test_false_pos=sum(sim_pop$sTested_fp, na.rm=TRUE),n_test_true_pos=sum(sim_pop$sTested_tp, na.rm=TRUE)),summarize_infectiousness(sim_pop,if_weights,day_of_flight,subclin_infectious),by="time_step"))
			} else {
				res <- bind_rows(res,left_join(summarize_sim_state(sim_pop),summarize_infectiousness(sim_pop,if_weights,day_of_flight,subclin_infectious),by="time_step"))
			}
		}
	}
	if (testing_type == "pcr_five_days_after") {
		res <- tibble()
		for (k in n_burnin:(day_of_flight - 1)) {
			sim_pop <- increment_sim(sim_pop=sim_pop,t_step=k,prob_inf=prob_inf,if_weights=if_weights,days_quarantine=days_quarantine)
			res <- bind_rows(res,left_join(summarize_sim_state(sim_pop),summarize_infectiousness(sim_pop,if_weights,day_of_flight,subclin_infectious),by="time_step"))
		}
		if (symptom_screening) {
			sim_pop[DayInObs == 0 & InfectionType == 1 & State == 2, DayOfDetection := day_of_flight]
		}
		sim_pop <- increment_sim(sim_pop=sim_pop, t_step=day_of_flight,prob_inf=prob_inf * risk_multiplier,if_weights=if_weights,days_quarantine=days_quarantine)
		res <- bind_rows(res,left_join(summarize_sim_state(sim_pop),summarize_infectiousness(sim_pop,if_weights,day_of_flight,subclin_infectious),by="time_step")) 
		test_days <- day_of_flight + 5
		for (k in (day_of_flight + 1):(day_of_flight + n_outcome_day)) {
			if (k %in% test_days) {sim_pop <- pcr_test_subset_simpop(sim_pop,k,n_subsets=NROW(test_days),subset_to_test=which(test_days == k),p_sens,p_spec,days_incubation,days_symptomatic)}  
			if (k > (day_of_flight + 3)) {
				sim_pop <- increment_sim(sim_pop=sim_pop,t_step=k,prob_inf=0,if_weights=if_weights,days_quarantine=days_quarantine)
			} else {
				sim_pop <- increment_sim(sim_pop=sim_pop,t_step=k,prob_inf=prob_inf,if_weights=if_weights,days_quarantine=days_quarantine)
			}
			if (k %in% test_days) {
				res <- bind_rows(res,left_join(summarize_sim_state(sim_pop) %>% mutate(n_test_false_pos=sum(sim_pop$sTested_fp, na.rm=TRUE),n_test_true_pos=sum(sim_pop$sTested_tp, na.rm=TRUE)),summarize_infectiousness(sim_pop,if_weights,day_of_flight,subclin_infectious),by="time_step")) 
			} else {
				res <- bind_rows(res,left_join(summarize_sim_state(sim_pop),summarize_infectiousness(sim_pop,if_weights,day_of_flight,subclin_infectious),by="time_step"))
			}
		}
	}
	### Close out simulations. Store summarized simulation results with parameter data ----
	rm(sim_pop); invisible(gc(FALSE))	
	holder[[r]] <- res %>%
	mutate(testing_type=testing_type,symptom_screening=symptom_screening,sens_type=sens_type,round=round,rep=r,prob_inf=prob_inf,prop_subclin=prop_subclin, risk_multiplier,rapid_test_multiplier=rapid_test_multiplier)
	}
	saveRDS(bind_rows(holder),file=sim_file,compress="xz")
	print(sprintf("Finished (%s): %s", Sys.time(), sim_file))
} else {
	print(sprintf("Skipping (%s): %s", Sys.time(), sim_file))
}
sink()
}
