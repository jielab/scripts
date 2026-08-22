###################################################
####       Leave one feature out (LOFO)        ####
####    framework for obesity complications    ####
###################################################

# This script runs the LOFO framework across outcomes and features:
# for each outcome, a baseline performance is derived using the 
# full pool of features of top 20 across outcomes
# leaving out each feature for each outcome, the lost performance is calculated

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

## load data
ukb.phe <-  fread("<path>")

## define outcomes
outc_list <- c("<outcomes>")

## define follow-up
cdat_list <- c("<followup>")

## get labels
lab.phe <- fread("<path>")

#----------------------------#
##--     top features     --##
#----------------------------#
## list to store
var_selection_simple <- list()

## get top 20 features across outcomes
for (outc in outc_list) {
  
  ## load clinical features (i.e. up to body composition)
  feat.sel.simple <- fread(paste0("<path>", outc, ".txt"))
  
  ## order
  feat.sel.simple <- data.frame(feat.sel.simple) %>% arrange(-sel, -src)
  jj <- grep("rand", feat.sel.simple$id)
  var.sel.simple <- feat.sel.simple[1:(jj[1]-1),]
  var.sel.simple <- as.data.table(var.sel.simple)
  
  ## get top 20 and merge with labels
  var.sel.simple <- left_join(var.sel.simple[1:20,c(1,253)],lab.phe[,c(4,6,9)], by=c("id"="short_name_new"))
  var.sel.simple <- var.sel.simple[,Outcome := outc]
  var_selection_simple[[outc]] <- var.sel.simple
}

## combine everything
var_selection_simple <- bind_rows(var_selection_simple)

## pull top features
top <- var_selection_simple %>%  distinct(id) %>% pull(id)

## set aside outcomes
ukb.outcfol <- ukb.phe %>%
  select("f.eid","<outcomes>", "<followup>")

## keep labs needed
lab.phe %<>% filter(short_name_new %in% top)

intersect <- intersect(colnames(ukb.phe), lab.phe$short_name_new)
eid_split_col <- c("f.eid", "split")
intersect <- c(eid_split_col,intersect)
ukb.phe <- ukb.phe[, ..intersect]
ukb.phe <-  na.omit(ukb.phe)

## clear some space
gc(reset=T)

## norm
binary_cols <- sapply(ukb.phe, function(col) all(col %in% c(0, 1)))
continuous_cols <- !binary_cols & names(ukb.phe) != "f.eid" & names(ukb.phe) != "split"
continuous_cols <- names(ukb.phe)[continuous_cols]  ## get actual column names
ukb.phe[, (continuous_cols) := lapply(.SD, scale), .SDcols = continuous_cols]

#----------------------------#
##--    prune variables   --##
#----------------------------#
set.seed(123)

## sample a fraction of the dataset
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
cont.net <- graph_from_data_frame(cont.cor[ value >= .6])
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

## manually keep whtr and body fat percentage remove waist
## keep
cont.net$ind[ cont.net$short_name == "whtr"]   <- 1
cont.net$ind[ cont.net$short_name == "Body_fat_percentage"]   <- 1
## drop
cont.net$ind[ cont.net$short_name == "waist"] <- 2

## manually keep HbA1c and T2D instead of Metformin
## keep
cont.net$ind[ cont.net$short_name == "hba1c"]          <- 1
cont.net$ind[ cont.net$short_name == "bin_250.2"]          <- 1
## drop
cont.net$ind[ cont.net$short_name == "A10BA02"]   <- 2

#--------------------------------------#
##--   final feature set for pred.  --##
#--------------------------------------#
lab.final <- lab.phe[ !(short_name %in% c(cont.net[ ind != 1]$short_name))]

lab.phe <- lab.final

top <- lab.phe$short_name_new

intersect <- intersect(colnames(ukb.phe), lab.phe$short_name_new)
eid_split_col <- c("f.eid", "split")
intersect <- c(eid_split_col,intersect)
ukb.phe <- ukb.phe[, ..intersect]
ukb.phe <-  na.omit(ukb.phe)

## add outcomes and follow-up
ukb.phe.original <- left_join(ukb.phe, ukb.outcfol, by="f.eid")

## clear some space
rm(ukb.phe, ukb.outcfol, ukb.phe.sample); gc(reset=T)

########################################
####               LOFO             ####
########################################

## --> load function <-- ##
source("<path>/coxnet_ridge_optim_boot.R")

## list to store
all_leave_one_out_results <- list()

## parallel
registerDoMC(5)
## set seed
set.seed(123)

## loop through outcomes
for (i in seq_along(outc_list)) {
  outc <- outc_list[i]
  cdat <- cdat_list[i]
  
  cat(paste0("Analysing outcome: ", outc, "\n"))
  
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
  
  ## --> baseline performance <-- ##
  baseline_model <- coxnet.optim.r(ukb.phe[ukb.phe$f.eid %in% t.idx, ], top, 
                                   Surv(ukb.phe[ukb.phe$f.eid %in% t.idx, cdat],
                                        ukb.phe[ukb.phe$f.eid %in% t.idx, outc]),
                                   1000, ukb.phe[ukb.phe$f.eid %in% o.idx, ],
                                   Surv(ukb.phe[ukb.phe$f.eid %in% o.idx, cdat],
                                        ukb.phe[ukb.phe$f.eid %in% o.idx, outc]))
  ## df baseline performance
  baseline_cindex <- data.frame(
    Outcome = outc,
    ExcludedFeature = "None",
    Cindex.mean = baseline_model$Cindex.mean,
    ci.low = baseline_model$ci.low,
    ci.upp = baseline_model$ci.upp
  )
  
  ## list for LOFO results
  leave_one_out_results <- list()
  
  ## loop through features to exclude
  for (feature_to_exclude in top) {
    cat(paste0("Excluding feature: ", feature_to_exclude, "\n"))
    
    ## features after left out
    reduced_features <- setdiff(top, feature_to_exclude)
    
    ## --> LOFO peformance <-- ##
    model <- coxnet.optim.r(ukb.phe[ukb.phe$f.eid %in% t.idx, ], reduced_features, 
                            Surv(ukb.phe[ukb.phe$f.eid %in% t.idx, cdat],
                                 ukb.phe[ukb.phe$f.eid %in% t.idx, outc]),
                            1000, ukb.phe[ukb.phe$f.eid %in% o.idx, ],
                            Surv(ukb.phe[ukb.phe$f.eid %in% o.idx, cdat],
                                 ukb.phe[ukb.phe$f.eid %in% o.idx, outc]))
    
    ##df LOFO performance
    cindex_diff <- data.frame(
      Outcome = outc,
      ExcludedFeature = feature_to_exclude,
      Cindex.diff = model$Cindex.mean - baseline_model$Cindex.mean,
      ci.low.diff = model$ci.low - baseline_model$ci.low,
      ci.upp.diff = model$ci.upp - baseline_model$ci.upp
    )
    
    ## store result for each feature exclusion
    leave_one_out_results[[feature_to_exclude]] <- cindex_diff
  }
  
  ## combine results for the current outcome
  all_leave_one_out_results[[outc]] <- do.call(rbind, leave_one_out_results) %>%
    bind_rows(baseline_cindex, .)
}

## combine results for across outcomes
final_leave_one_out_results <- do.call(rbind, all_leave_one_out_results)

## write results
fwrite(final_leave_one_out_results, "<path>", sep="\t", row.names=FALSE, na=NA)
