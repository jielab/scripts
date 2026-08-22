
#Author: Linsey Jackson
#Date: 12/9/2025
#Script name: standardize & compute score
#Purpose: standardize predictors and compute OBSCORE risk scores 

library(dplyr)
library(data.table)
library(purrr)
library(haven)

#set a working directory
setwd("") 

#read in individual level predictors as cleaneddf_complete
#merging, formatting, and cleaning of ADaM tables not shown
#this step must be completed for each time point

#read in summary stats for predictors from UKB cohort
stats_df <- read.csv("summary_stats.csv")
means <- setNames(as.list(stats_df$mean), stats_df$variable)
sds <- setNames(as.list(stats_df$sd), stats_df$variable)

# Standardize predictors using the specified means and SDs from UKB cohort
df1 <- cleaneddf_complete %>%
  mutate(across(all_of(names(means)), 
                ~ (. - means[[cur_column()]]) / sds[[cur_column()]]))

#read in coefficients 
df2<-read.csv("OBSCORE_Light_Coefficients_23102025.csv")
rownames(df2) <- df2$Variable

#calculate risk scores 
predict_outcomes <- function(df_features, df_coefs, id_col = "ID") {
  # 1) Extract IDs and drop from features
  ids <- df_features[[id_col]]
  X <- df_features[, setdiff(names(df_features), id_col), drop = FALSE]
  
  # 2) Check that features match coefficients
  if (!all(rownames(df_coefs) %in% names(X))) {
    stop("Some features in coefficients do not match columns in the features data frame.")
  }
  # Reorder X to match coefficients
  X <- X[, rownames(df_coefs), drop = FALSE]
  
  # 3) Matrix multiply: subjects × features  %*%  features × outcomes
  scores <- as.matrix(X) %*% as.matrix(df_coefs)
  
  # 4) Return as data frame with IDs
  result <- data.frame(ID = ids, scores, check.names = FALSE)
  return(result)
}

predicted_scores <- predict_outcomes(df1, df2, id_col = "ID")