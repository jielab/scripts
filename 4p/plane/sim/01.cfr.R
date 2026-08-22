# https://github.com/thimotei/CFR_calculation 🏮在 LINUX系统里面运行！
pacman::p_load(data.table, dplyr, tidyverse, reshape2, lubridate, countrycode)# 必须先在sh里面 conda create -n r-reticulate python
library(reticulate); use_condaenv('r-reticulate', required=TRUE)# 必须先在R里面 install_greta_deps(timeout=50)
library(greta); library(greta.gp) # 必须先在sh里面 pip install tensorflow
#greta::install_tensorflow(method="conda", extra_packages="tensorflow-probability") # 安装一次即可

code <- "?"

dir='/work/sph-huangj'
dat.jhu <- paste0(dir,'/data/avia/jhu.rds')
cfr_baseline <- 1.4; cfr_range <- c(1.2, 1.7); cfr_trend=NULL
mean <- 13; median <- 9.1; n_inducing=5; 
mu_cfr <- log(median); sigma_cfr <- sqrt(2*(log(mean) - mu_cfr))
hospital_to_death <- function(x) { plnorm(x + 1, mu_cfr, sigma_cfr) - plnorm(x, mu_cfr, sigma_cfr)}


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 读入数据
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
import_jhu_data <- function (dat_file) {
	dat <- data.table::fread(dat_file) %>% 
    dplyr::select(-c("Province/State", "Lat", "Long")) %>% rename("country"="Country/Region") %>% 
	melt(., id.vars="country", variable.name="date", value.name="new_cases") %>% mutate(date := mdy(date)) %>%
	group_by(date, country) %>% summarize(new_cases=sum(as.numeric(new_cases))) %>%
	group_by(country) %>% mutate(new_cases = new_cases - lag(new_cases)) %>% as.data.table() %>%
	mutate(new_cases =ifelse(new_cases < 0, 0, new_cases)) %>% 
	mutate(iso3c := countrycode(country, 'country.name', 'iso3c', custom_match=c("Micronesia"="FSM", "Kosovo"="XXK"))) %>% 
	filter(!is.na(new_cases), !is.na(iso3c)) %>% setcolorder("iso3c") %>% .[order(country, date)] %>% unique()
	return(dat)
}
if (file.exists(dat.jhu)) {
	dat0 <- readRDS(dat.jhu)
} else {
	dat.case <- import_jhu_data(paste0(dir,'/data/avia/covid19_confirmed_global.csv')) 
	dat.death <- import_jhu_data(paste0(dir,'/data/avia/covid19_deaths_global.csv')) %>% rename(new_deaths=new_cases) 
	dat0 <- merge.data.table(dat.case, dat.death, by=c("country", "iso3c", "date"), all=FALSE) 
	saveRDS(dat0, dat.jhu)
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 计算CFR
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
dat <- dat0 %>% filter(iso3c==code) %>% dplyr::mutate(date=date - mean)
	if (nrow(dat)==0) break
	death_threshold_date <- dat %>% mutate(death_cum_sum=cumsum(new_deaths)) %>% filter(death_cum_sum >= 10) %>% pull(date) %>% min() 
	reporting_spike <- dat %>% filter(new_deaths < 0) %>% pull(date) %>% min() # date where negative death reporting spikes occur
	case_incidence <- dat$new_cases; death_incidence <- dat$new_deaths
	cumulative_known_t <- NULL # cumulative cases with known outcome at time tt
	for(ii in 1:nrow(dat)){ # Sum over cases up to time tt
		known_i <- 0 # number of cases with known outcome at time ii
		for(jj in 0:(ii - 1)){
			known_jj <- (case_incidence[ii - jj]*hospital_to_death(jj))
			known_i <- known_i + known_jj
		}
		cumulative_known_t <- c(cumulative_known_t,known_i) # Tally cumulative known
	}
	b_tt <- sum(death_incidence)/sum(case_incidence) # naive CFR value
	p_tt <- (death_incidence/cumulative_known_t) %>% pmin(.,1) # corrected CFR estimator
	cfr <- data.frame(nCFR=b_tt, ccfr=p_tt,total_deaths=sum(death_incidence), deaths =death_incidence,cum_known_t=round(cumulative_known_t), total_cases=sum(case_incidence)) %>% dplyr::as_tibble() %>%
	mutate(reporting_estimate=cfr_baseline/ccfr) %>% mutate(reporting_estimate=pmin(reporting_estimate, 1), country=dat$country, date=dat$date, date_num=as.numeric(dat$date), deaths=dat$new_deaths, cases_known=cum_known_t) %>% 
		filter(date >= death_threshold_date, date < reporting_spike) %>% 
		dplyr::select(country, date, date_num, reporting_estimate, deaths, cases_known) # %>% filter(date > "2020-09-01")
	if (nrow(cfr)==0) break


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 贝叶斯
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
cfr <- cfr %>% mutate(cases_known=ifelse(cases_known==0, 1, cases_known)) ## !! 🏮🐕
	if (nrow(cfr)==0) break
	write.table(cfr, file=paste0(code,'.cfr.txt'), append=FALSE, quote=FALSE, sep='\t', row.names=FALSE, col.names=TRUE)
	# %>% complete(date = seq.Date(min(date), max(date), by="day")) %>% fill(col1, col2)
	n <- nrow(cfr)
	times <- seq(min(cfr$date_num), max(cfr$date_num)) ## times <- as.integer(cfr$date_num)
	if (n != length(times)) { print(paste('ERROR: cfr N and Dates different', n, length(times))); break }
	lengthscale <- greta::lognormal(4, 0.5)
	sigma <- greta::lognormal(-1, 1)
	temporal <- greta.gp::rbf(lengthscales=lengthscale, variance=sigma ^ 2)
	intercept <- greta.gp::bias(1)
	reporting_kernel <- temporal + intercept
	# IID noise kernel for observation overdispersion (clumped death reports)
	sigma_obs <- greta::normal(0, 0.5, truncation=c(0, Inf))
	observation_kernel <- greta.gp::white(sigma_obs ^ 2)
	# combined kernel (marginalises a bunch of parameters for easier sampling)
	kernel <- reporting_kernel + observation_kernel
	inducing_points <- seq(min(times), max(times), length.out=n_inducing + 1)[-1] 
	z <- greta.gp::gp(times, inducing=inducing_points, kernel) 
	reporting_rate <- greta::iprobit(z)
	true_cfr_mean <- cfr_trend
	if (is.null(cfr_trend)) {true_cfr_mean <- cfr_baseline}
	true_cfr_sigma <- mean(abs(cfr_range - cfr_baseline)) / 1.96
	cfr_error <- greta::normal(0, true_cfr_sigma)
	baseline_cfr_perc <- greta::normal(true_cfr_mean, true_cfr_sigma, truncation=c(0, 100))
	log_expected_deaths <- log(baseline_cfr_perc / 100) + log(cfr$cases_known) - log(reporting_rate)
	expected_deaths <- exp(log_expected_deaths)
	greta::distribution(cfr$deaths) <- greta::poisson(expected_deaths)
	# construct the model
	m <- greta::model(reporting_rate)
	n_chains <- 200
	inits <- replicate(n_chains, greta::initials(lengthscale=rlnorm(1, 4, 0.5),sigma=abs(rnorm(1, 0, 0.5)),sigma_obs=abs(rnorm(1, 0, 0.5)),
		baseline_cfr_perc=max(0.001, min(99.999,rnorm(1, true_cfr_mean, true_cfr_sigma)))),simplify=FALSE
	)
	country <- cfr$country[1]; message("running model for ", country)
	draws <- greta::mcmc( m,sampler=greta::hmc(Lmin=15, Lmax=20),chains=n_chains, warmup=1000, n_samples=1000,initial_values=inits,one_by_one=TRUE,verbose=TRUE)
	# extend the number of chains until convergence (or give up)
	for (i in 1:5) {
		r_hats <- coda::gelman.diag(draws, autoburnin=FALSE, multivariate=FALSE)$psrf[, 1]
		n_eff <- coda::effectiveSize(draws)
		decent_samples <- max(r_hats) <= 1.1 & min(n_eff) > 1000 
		if (!decent_samples) {
			message("maximum R-hat: ", max(r_hats), "\nminimum n-eff: ", min(n_eff))
			draws <- greta::extra_samples(draws, 2000, one_by_one=TRUE, verbose =TRUE)
		}
	}
	# predict without IID noise (true reporting rate, without clumped death reporting)
	z_smooth <- greta.gp::project(z, times, kernel=reporting_kernel)
	reporting_rate_smooth <- greta::iprobit(z_smooth)
	draws_pred <- greta::calculate(reporting_rate_smooth, values=draws)
	draws_pred_mat <- as.matrix(draws_pred)
	# compute posterior mean and 95% credible interval and return
	pred <- tibble::tibble(date=cfr$date,estimate=colMeans(draws_pred_mat),lower=apply(draws_pred_mat, 2, quantile, 0.025),upper=apply(draws_pred_mat, 2, quantile, 0.975))
	write.table(pred, file=paste0(code,'.pred.txt'), append=FALSE, quote=FALSE, sep='\t', row.names=FALSE, col.names=TRUE)
#	ci_poly <- tibble::tibble(x=c(dat$date, rev(dat$date)), y=c(pred$upper, rev(pred$lower)))


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 分析汇总
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
dir="D:/data/avia"
airport <- read.csv(paste0(dir,"/iata-icao.csv"),header=TRUE) %>% 
	mutate(iso3c=countrycode(country_code, "iso2c", "iso3c")) %>% dplyr::select(-c(country_code,region_name,icao))
guests <- read.csv(paste0(dir,"/schedule/2019_01.csv"), header=T) %>% 
	filter(Time.series >= as.Date("2019-01-21") & Time.series <= as.Date("2019-01-27")) %>%
	merge(airport, by.x="Origin", by.y="iata") %>% merge(airport, by.x="Destination", by.y="iata") %>%
	filter(iso3c.x !="CHN" & iso3c.y=="CHN") %>% rename(iso3c=iso3c.x) %>% 
	dplyr::select(iso3c, Seats) %>% group_by(iso3c) %>% summarize(Seats =sum(Seats))
pop <- subset(world_bank_pop, indicator=="SP.POP.TOTL", select=c("country", "2017")) %>% rename(iso3c=country, pop="2017")
unrate <- read.table(paste0(dir,"/cfr.txt"),header=TRUE) %>% mutate(estimate=as.numeric(estimate), date=as.Date(unrate$date))
jhu.raw <- data.table::fread(paste0(dir,"/covid19_confirmed_global.csv")) %>% 
    .[, c("Province/State", "Lat", "Long") := NULL] %>% melt(., id.vars="Country/Region", variable.name="date", value.name="new_cases") %>% 
    setnames(., old="Country/Region", new="country") %>% .[, date := mdy(date)] %>%
	.[, new_cases := sum(as.numeric(new_cases)), by=c("date", "country")] %>% unique() %>% 
	.[, iso3c := countrycode(country, 'country.name', 'iso3c', custom_match=c("Micronesia"="FSM", "Kosovo"="XXK"))] %>% 
	merge(unrate, by=c("iso3c","date")) %>% merge(pop, by="iso3c") %>% 
	mutate(truerate=(new_cases/estimate)/pop, week=ceiling(as.numeric(date-as.Date("2020-01-19"))/7)) #把1月19日礼拜天作为起点，1月22日开始有JHU感染数据。
jhu <- jhu.raw %>% filter(week ==10) %>% # 2020年1月19日礼拜天到2023年1月8日礼拜天，一共155周。week==10🏮，2020年3月26日，是第10周开始，民航开始“五个一”，至2023年1月8日结束。 
	dplyr::select(iso3c, truerate, week) %>% group_by(iso3c,week) %>% summarize(truerate =mean(truerate,na.rm=TRUE)) 
	sPDF <- joinCountryData2Map(jhu, joinCode = "ISO3", nameJoinColumn = "iso3c")
	mapPlot <- mapCountryData(sPDF, nameColumnToPlot = "truerate", catMethod = "quantiles", colourPalette = "heat", addLegend = TRUE)
dat <- merge(jhu, guests, by="iso3c") %>% mutate(cases=round(Seats * truerate)); table(dat$iso3c); sum(dat$cases)


