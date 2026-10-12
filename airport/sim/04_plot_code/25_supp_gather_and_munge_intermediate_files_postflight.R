# 🚩 25_supp_gather_and_munge_intermediate_files_postflight
## 25_supp_gather_and_munge_intermediate_files_postflight.R ----
## Collect scenario outputs and compare runs with identical random draws.

## Imports ----
library(tidyverse)
library(here)
library(fs)
library(foreach)
library(doParallel)
source(here::here("code", "utils.R"))
fs::dir_create(here::here("data_raw"))

testing_scenarios <- basename(
	fs::dir_ls(here::here("intermediate_files"),
			   type = "directory",
			   regexp = "pcr|rapid|no_testing|perfect")
)

## Collect and save intermediate files in raw form ----
doParallel::registerDoParallel()
all_results <- foreach::foreach(s = testing_scenarios) %dopar% {
	s_name <- here::here("data_raw", sprintf("raw_simulations_%s_postflight.RDS", s))

	temp_x <- collect_and_munge_post_flight_simulations(s)

	### Save raw simulations ----
	saveRDS(temp_x,
			s_name,
			compress = "xz")
}
doParallel::stopImplicitCluster()
