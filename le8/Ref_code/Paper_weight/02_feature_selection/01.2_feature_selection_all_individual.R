####################################################################
####         Individual domain top 20 feature selection         ####
####################################################################
#!/usr/bin/env Rscript

# This script runs feature selection, optimisation and testing
# for a set of predefined obesity complications using multimodal data

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
require(readxl)
require(dplyr)
require(magrittr)
require(caret)
require(ROSE)
require(rms)
require(survival)
require(Rfast)
require(influential)
require(sna)

## get outcome and follow-up from command line
outc <- args[1]
cdat <- args[2]

## import imputed UKB dataset
ukb.phe <-  fread("<path>")

## import labels for columns 
lab.phe <- fread("<path>")

## outcomes and follow-up
ukb.outcfol <- ukb.phe %>%
  select("f.eid","<outcomes>", "<followup>")

## match
lab.phe %<>% filter(!short_name %in% "centre")
intersect <- intersect(colnames(ukb.phe), lab.phe$short_name_new)
eid_split_col <- c("f.eid", "split")
intersect <- c(eid_split_col,intersect)
ukb.phe <- ukb.phe[, ..intersect]

## clear some space
gc(reset=T)

#--------------------------------------#
##--    prune variables for LASSO   --##
#--------------------------------------#
set.seed(123)

## sample a fraction of the dataset for comp. feasibility
ukb.phe.sample <- ukb.phe %>%
  sample_frac(size = 0.03)

## compute simple correlation
cont.cor                                   <- lab.phe$short_name_new
## compute simple correlation
cont.cor                                   <- Rfast::cora(ukb.phe.sample[, cont.cor, with=F])
## drop redundancies
cont.cor[ lower.tri(cont.cor, diag = T)]   <- NA
## convert into table
cont.cor                                   <- reshape2::melt(cont.cor, na.rm=T)
cont.cor                                   <- as.data.table(cont.cor)
## convert to character
cont.cor[, Var1 := as.character(Var1)]
cont.cor[, Var2 := as.character(Var2)]

## map to network
cont.net <- graph_from_data_frame(cont.cor[ value >= .9])
## get all separate components
cont.net <- igraph::components(cont.net)$membership
## convert to data frame
cont.net <- data.table(short_name=names(cont.net), id.group=cont.net)
## add missingness
cont.net <- merge(cont.net, lab.phe)
## order accordingly
cont.net <- cont.net[ order(id.group, miss.per)]
## add indicator
cont.net[, ind := 1:.N, by="id.group"]

## manually change highly correlated due to clinical availability
## chol for ApoB and hdl_chol for ApoA1 and body fat % for leg fat $
## keep 
cont.net$ind[ cont.net$short_name == "chol"]          <- 1
## drop
cont.net$ind[ cont.net$short_name == "ApoB"]    <- 2

## keep
cont.net$ind[ cont.net$short_name == "Body_fat_percentage"]         <- 1
## drop
cont.net$ind[ cont.net$short_name == "leg_fat_percentage__avg_"]    <- 2

## keep
cont.net$ind[ cont.net$short_name == "hdl_chol"]      <- 1
## drop
cont.net$ind[ cont.net$short_name == "ApoA1"]   <- 2

#--------------------------------------#
##--   final feature set for pred.  --##
#--------------------------------------#
## drop by indicator, i.e. lowest missingness
lab.final <- lab.phe[ !(short_name %in% c(cont.net[ ind != 1]$short_name))]
lab.phe <- lab.final

intersect <- intersect(colnames(ukb.phe), lab.phe$short_name_new)
eid_split_col <- c("f.eid", "split")
intersect <- c(eid_split_col,intersect)
ukb.phe <- ukb.phe[, ..intersect]

## get binary cols
binary_cols <- sapply(ukb.phe, function(col) all(col %in% c(0, 1)))

## get cont cols, except f.eid
continuous_cols <- !binary_cols & names(ukb.phe) != "f.eid" & names(ukb.phe) != "split"
continuous_cols <- names(ukb.phe)[continuous_cols]

## standardise cont cols
ukb.phe[ , (continuous_cols) := lapply(.SD, scale), .SDcols = continuous_cols]

## add ten random variables
for(j in 1:10){
  ukb.phe[, paste0("rand_", j)] <- rbinom(nrow(ukb.phe), size=1, p=1/j)
}

## create a new data frame with the additional rows
new_rows <- data.frame(short_name_new = paste0("rand_", 1:10))

lab.phe <- as.data.table(lab.phe)
new_rows <- as.data.table(new_rows)

## append the new rows
lab.phe <- rbindlist(list(lab.phe, new_rows), fill = TRUE)

## add outcomes and follow-up
ukb.phe.original <- left_join(ukb.phe, ukb.outcfol, by="f.eid")

## clear some space
rm(new_rows,ukb.phe, ukb.outcfol); gc(reset=T)

#--------------------------------------#
##--     categories of predictors   --##
#--------------------------------------#
lab.phe <- lab.phe %>%
  mutate(
    category_group = case_when(
      category %in% c("Ancestry", "Basic demographics", "Health", "Family history", "Diet", "Socioeconomic") ~ "Basic parameters",
      category %in% c("Biomarker", "Blood cell counts") ~ "Blood markers",
      category %in% c("Cardiovascular", "Pulmonary") ~ "Cardiopulmonary parameters",
      category %in% c("Body composition") ~ "Body composition",
      category %in% c("Diseases", "Drugs") ~ "Diseases and drug intake",
      category %in% c("NMR") ~ "NMR",
      category %in% c("PGS") ~ "PGS",
      TRUE ~ "Random"
    )
  )

## define order of individual data domains
## data-domain specific outcome models
category_groups <- list(
  "Diseases & Drugs" = c("Diseases and drug intake", "Random"),
  "Biomarker" = c("Blood markers", "Random"),
  "Cardiopulmonary" = c("Cardiopulmonary parameters", "Random"),
  "Body composition" = c("Body composition", "Random"),
  "NMR" = c("NMR","Random"),
  "PGS" = c("PGS","Random")
)

########################################
####        feature selection       ####
########################################

## --> load function <-- ##
source("<path>/coxnet_ridge_optim_boot.R")

## run parallel
registerDoMC(5)

## for each data modality
for (group_name in names(category_groups)) {
  cat("\nProcessing category group:", group_name, "\n")
  
  ## get the variables in each category
  var.phe_group <- lab.phe[category_group %in% category_groups[[group_name]]]$short_name_new
  
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
  
  ## seed
  set.seed(123)
  ## draw 250 random subsets of 40% 
  ll <- lapply(1:250, function(x) sample(t.idx, ceiling(length(t.idx)*.4)))
  
  cat("\nProcessing outcome:", outc, "with follow-up time:", cdat, "\n")
  
  ## --> time start <-- ##
  start_time <- proc.time()
  
  ## run feature selection
  feat <- mclapply(1:250, function(i) {
    message("Processing iteration: ", i)
    x <- ll[[i]]
    
    var.phe <- c(var.phe_group, paste0("rand_", 1:10))  # include random variables
    
    tmp <- train(ukb.phe[ukb.phe$f.eid %in% x, c(var.phe)], 
                 as.factor(ukb.phe[ukb.phe$f.eid %in% x, outc]), 
                 method="glmnet", family="binomial", metric="Kappa",
                 tuneGrid=as.data.frame(expand_grid(alpha=1, lambda=10^-seq(6,.25,-.5))),
                 trControl= trainControl(method = "repeatedcv", number = 3, repeats = 3,
                                         sampling = "rose", allowParallel = FALSE))
    
    return(list(mat = as.matrix(coef(tmp$finalModel, s = tmp$finalModel$lambdaOpt)), 
                lasso = tmp$results))
  }, mc.cores=5)
  
  end_time <- proc.time()
  ## --> time stop <-- ##
  
  time_taken <- end_time - start_time
  print(paste(time_taken[3]/60, "minutes")) ## print time needed for feature selection
    
  ## store feature selection matrix
  feat.sel <- do.call(cbind, lapply(feat, function(x) x$mat))
  feat.sel <- feat.sel[-1,]  # drop intercept
  feat.sel <- data.table(id=rownames(feat.sel), feat.sel)
  feat.sel[, src := apply(feat.sel[,-1], 1, function(x) median(abs(x)))]
  feat.sel[, sel := apply(feat.sel[, -c(1, ncol(feat.sel)), with=FALSE], 1, function(x) sum(x != 0))]
    
  ## get performance measures from each run
  per.mat <- unlist(lapply(feat, function(x) return(max(x$lasso$Kappa, na.rm=TRUE))))
  feat.sel[, per := apply(feat.sel[, -c(1, ncol(feat.sel)-c(1,0)), with=FALSE], 1, function(x){
    x <- ifelse(x != 0, per.mat, 0)
    return(mean(x))
  })]
  
  ## store feature selection results
  fwrite(feat.sel, paste0("<path>", group_name, ".", outc, ".txt"), 
         sep="\t", row.names=FALSE, na=NA)
    
  ## select features based on:
  ## 1) times selected across subsamples and 
  ## 2) median absolute coefficient values across subsamples 
  feat.sel <- feat.sel[order(-sel,-src)]
  jj <- grep("rand", feat.sel$id)
  var.sel <- feat.sel[1:(jj[1]-1)]
  
  ## report number of overall non-random features
  cat(nrow(var.sel), "non-random features for", outc, "\n")
  
  ## select up to top 20 features from non-random features
  top.sel <- feat.sel[1:min(20, nrow(var.sel))]
  
  cat("1: Performing opt & test with top selected features\n")
  feature_top_name <- paste0("feature.top.", group_name)
  
  ## run optimisation and testing in the held out opt & test set using regularised Cox regression
  assign(feature_top_name, coxnet.optim.r(ukb.phe[ ukb.phe$f.eid %in% o.idx, ], top.sel$id, 
                                          Surv(ukb.phe[ ukb.phe$f.eid %in% o.idx, cdat],
                                               ukb.phe[ ukb.phe$f.eid %in% o.idx, outc]),
                                          1000, ukb.phe[ ukb.phe$f.eid %in% v.idx, ],
                                          Surv(ukb.phe[ ukb.phe$f.eid %in% v.idx, cdat],
                                               ukb.phe[ ukb.phe$f.eid %in% v.idx, outc])))
  
  ## save model results from optimisation and testing
  save(list=feature_top_name, file=paste0("<path>", group_name, ".", outc, ".RData"))
}
