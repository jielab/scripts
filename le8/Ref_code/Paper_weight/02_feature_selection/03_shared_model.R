###################################################
####         Shared model performance          ####
###################################################

# This script trains, tests the shared model

rm(list=ls())
args <- commandArgs(TRUE)
options(stringsAsFactors = F)
options(scipen = 1)
setwd("<path>")

## --> packages needed <-- ##
require(data.table)
require(tidyverse)
require(doMC)
require(glmnet)
require(dplyr)
require(magrittr)
require(caret)
require(ROSE)
require(rms)
require(survival)
require(Rfast)
require(influential)
require(sna)
require(ggsci)

## load data
ukb.phe <-  fread("<path>")

## load labels
lab.phe <- fread("<path>")

## load LOFO results
final_leave_one_out_results <- fread("<path>")

## define outcomes
outc_list <- c("<outcomes>")

## define follow-up
cdat_list <- c("<followup>")

## set aside outcomes
ukb.outcfol <- ukb.phe %>%
  select("f.eid","<outcomes>", "<followup>")

lab.phe %<>% filter(!short_name %in% "centre")

intersect <- intersect(colnames(ukb.phe), lab.phe$short_name_new)
eid_split_col <- c("f.eid", "split")
intersect <- c(eid_split_col,intersect)
ukb.phe <- ukb.phe[, ..intersect]
ukb.phe <-  na.omit(ukb.phe)

## norm
binary_cols <- sapply(ukb.phe, function(col) all(col %in% c(0, 1)))
continuous_cols <- !binary_cols & names(ukb.phe) != "f.eid" & names(ukb.phe) != "split"
continuous_cols <- names(ukb.phe)[continuous_cols]  # Get actual column names
ukb.phe[, (continuous_cols) := lapply(.SD, scale), .SDcols = continuous_cols]

## add ten random variables
for(j in 1:10){ukb.phe[, paste0("rand_", j)] <- rbinom(nrow(ukb.phe), size=1, p=1/j)}
new_rows <- data.frame(short_name_new = paste0("rand_", 1:10))
lab.phe <- as.data.table(lab.phe)
new_rows <- as.data.table(new_rows)
## append the new rows
lab.phe <- rbindlist(list(lab.phe, new_rows), fill = TRUE)

## add outcomes and follow-up
ukb.phe.original <- left_join(ukb.phe, ukb.outcfol, by="f.eid")

## clear some space
rm(new_rows,ukb.phe, ukb.outcfol); gc(reset=T)

########################################
####           train-test           ####
########################################

## --> load function <-- ##
source("scripts/coxnet_ridge_optim_boot.R")

## parallel
registerDoMC(5)

## seed
set.seed(123)

## get mean across outcomes for each feature and take top 20
## testing of other benchmark models
## can be done by integrating the relevant features
top <- final_leave_one_out_results %>% filter(ExcludedFeature != "None") %>%
  group_by(ExcludedFeature) %>%
  summarise(mean_Cindex_diff = mean(Cindex.diff, na.rm = TRUE)) %>%
  arrange(desc(-mean_Cindex_diff)) %>% slice(1:20) %>% pull(ExcludedFeature)

## list to store results
all_cindex_data <- list()

## loop across outcomes
for(i in seq_along(outc_list)){
  outc <- outc_list[i]
  cdat <- cdat_list[i]
  
  ## set seed
  set.seed(123)
  
  ## remove participants with event during first 6 months
  ukb.phe <- ukb.phe.original %>% filter(!(((eval(as.name(cdat)) < 0.5)) & (eval(as.name(outc)) == 1)))
  
  ukb.phe <- as.data.frame(ukb.phe)
  ##--  create data splits  --##
  # training
  t.idx <- ukb.phe$f.eid[ukb.phe$split=="t.idx"]
  # optimization
  o.idx <- ukb.phe$f.eid[ukb.phe$split=="o.idx"]
  # remaining samples for validation
  v.idx <- ukb.phe$f.eid[ukb.phe$split=="v.idx"]
  # convert to data frame
  ukb.phe  <- as.data.frame(ukb.phe)
    
  cat("1: Performing opt & test with top selected features\n")
    
  ## train test
  model <- coxnet.optim.r(ukb.phe[ ukb.phe$f.eid %in% o.idx, ], top, 
                          Surv(ukb.phe[ ukb.phe$f.eid %in% o.idx, cdat],
                               ukb.phe[ ukb.phe$f.eid %in% o.idx, outc]),
                          1000, ukb.phe[ ukb.phe$f.eid %in% v.idx, ],
                          Surv(ukb.phe[ ukb.phe$f.eid %in% v.idx, cdat],
                               ukb.phe[ ukb.phe$f.eid %in% v.idx, outc]))
  
  feature_top_name <- paste0("feature.LOFO.mean", outc)
  assign(feature_top_name, model)
  
  ## save model
  save(list=feature_top_name, file=paste0("<path>", outc, ".RData"))
  
  ## df performances
  cindex_data <- data.frame(
    Outcome = outc,
    Cindex.mean = model$Cindex.mean,
    ci.low = model$ci.low,
    ci.upp = model$ci.upp,
    Nfeatures = "20",
    Model = "LOFO20_Mean")
  
  all_cindex_data <- append(all_cindex_data, list(cindex_data))
}

## combine all
all_res <- do.call(rbind, all_cindex_data)

## save performances of shared model
fwrite(all_res, "<path>", sep="\t", row.names=FALSE, na=NA)