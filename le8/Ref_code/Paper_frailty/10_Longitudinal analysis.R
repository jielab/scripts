# =====================================================================================
# Longitudinal Analyses
# Xueqing Jia, 2025
# Description: Longitudinal analyses, including the calculation of PFS/FI change rates, 
#   and evaluation of their associations with baseline age, baseline frailty severity, 
#   and accumulative disease counts. The code used for developing the feature-reduced
#   PFS model see "5_Development of PFS model"
# =====================================================================================

##### Load R packages
library(ggplot2)
library(impute)
library(glmnet)
library(caret)
library(dplyr)
library(plyr)
library(dplyr)
library(fdrci)

##### 1 Calculate PFS for the third and fourth follow-up visits using the baseline LASSO model for which the participant was included in the baseline testing (hold-out) set ######
### Find intersection of samples across three waves
load(file = "Protein_matrix_knn_re2.Rdata")
sample_re2<-rownames(Protein_matrix_knn)

load(file = "Protein_matrix_knn_re3.Rdata")
sample_re3<-rownames(Protein_matrix_knn)

load(file = "Analysis sample_split_longi.Rdata")
sample_ba<-temp_data_split$n_eid

sample_inter <- intersect(intersect(sample_ba, sample_re2), sample_re3)

### Load baseline PFS results
pred_FI_ba <- read.csv(file = "./Results/results_nested_cv_longi/pred_vs_true_allfolds.csv")
pred_FI_ba <- pred_FI_ba[pred_FI_ba$n_eid %in% sample_inter, c("n_eid", "pred_exp")]
colnames(pred_FI_ba)[2] <- "pred_FI_ba"

### Define function to calculate PFS for follow-ups
PFS_calculator_fun <- function(data, EN_list, mm_sd_list){
  
  pred_list <- vector("list", 10)
  
  for (i in c(1:10)) {
    
    cat("Processing fold", i, "\n")
    
    # Extract testing set data
    data_tem <- data[data$model == paste0("model_",i), ]
    
    # Extract protein column names
    protein_names <- colnames(data_tem)[37:1495]
    
    # Extract and order standardization parameters (named vectors)
    mm <- setNames(mm_sd_list[[i]]$mean$mean, rownames(mm_sd_list[[i]]$mean))
    sd <- setNames(mm_sd_list[[i]]$sd$sd, rownames(mm_sd_list[[i]]$sd))
    
    mm <- mm[protein_names]
    sd <- sd[protein_names]
    
    if (any(names(mm) != protein_names) || any(names(sd) != protein_names)) {
      stop("Standardization parameters do not match data columns (fold ", i, ")")
    }
    
    # standardization
    data_scaled  <- as.data.frame(scale(as.matrix(data_tem[, protein_names]),  center = mm, scale = sd))
    
    # Extract model coefficients
    var_list <- EN_list[[i]][EN_list[[i]] != 0]
    intercept <- var_list[names(var_list) == "(Intercept)"]
    coef_vec <- var_list[names(var_list) != "(Intercept)"]
    
    cat("Number of non-zero proteins (excluding intercept):", length(coef_vec), "\n")
    
    #  Calculate predicted values: Xβ + intercept
    data_matched <- data_scaled[, names(coef_vec), drop = FALSE]
    pred_scores <- as.matrix(data_matched) %*% coef_vec + intercept
    
    #  Build results data frame
    result_df <- data.frame(
      n_eid = rownames(data_matched),
      pred_FI = as.numeric(pred_scores),
      model_fold = paste0("model_",i),
      stringsAsFactors = FALSE
    )
    
    pred_list[[i]] <- result_df
  }
  final_pred <- do.call(rbind, pred_list)
  final_pred <- final_pred[, c("n_eid", "pred_FI", "model_fold")]
  final_pred$pred_FI <- exp(final_pred$pred_FI)-0.1
  
  return(final_pred)
}

## calculate PFS for the third follow-up
load(file = "./Results/results_nested_cv_longi/EN_coeffs_list.Rdata")
load(file = "./Results/results_nested_cv_longi/mm_sd_list.Rdata")

load(file = "./Datasets/Analysis sample_split_longi.Rdata")
protein_names <- colnames(temp_data_split)[37:1495]

load(file = "Protein_matrix_knn_re2.Rdata")
colnames(Protein_matrix_knn)<-gsub("-","_",colnames(Protein_matrix_knn))

data<-Protein_matrix_knn[sample_inter,protein_names]
data$n_eid<-rownames(data)
data<-merge(temp_data_split[,c(1:36)],data, by="n_eid")
rownames(data)<-data$n_eid

pred_FI_re2<-PFS_calculator_fun(data, EN_list, mm_sd_list)
colnames(pred_FI_re2)[2] <- "pred_FI_re2"

## calculate PFS for the fourth follow-up
load(file = "Protein_matrix_knn_re3.Rdata")
colnames(Protein_matrix_knn)<-gsub("-","_",colnames(Protein_matrix_knn))

data<-Protein_matrix_knn[sample_inter,protein_names]
data$n_eid<-rownames(data)
data<-merge(temp_data_split[,c(1:36)],data, by="n_eid")
rownames(data)<-data$n_eid

pred_FI_re3<-PFS_calculator_fun(data, EN_list, mm_sd_list)
colnames(pred_FI_re3)[2] <- "pred_FI_re3"

## Merge PFS results from all three waves
final_pred <- merge(pred_FI_ba, merge(pred_FI_re2[,-3], pred_FI_re3, by = "n_eid"), by = "n_eid")

load(file = "Analysis sample_split_longi.Rdata")
final_pred <- merge(temp_data_split[,1:36],final_pred,by="n_eid")

save(final_pred, file = "./Results/results_nested_cv_longi/final_pred_longi.Rdata")


##### 2 Method 1: Calculation of the change rate of PFS and FI using (later wave - previous wave) / follow-up time #######
load(file = "./Results/results_nested_cv_longi/final_pred_longi.Rdata")
followup_date<-read.csv(file = "ukb_field53_XJ.csv")

temp_data<-merge(final_pred, followup_date,by="n_eid")

temp_data$followup_time_ba_re2<-(temp_data$s_53_2_0-temp_data$s_53_0_0)/30/12
temp_data$followup_time_ba_re3<-(temp_data$s_53_3_0-temp_data$s_53_0_0)/30/12
temp_data$followup_time_re2_re3<-(temp_data$s_53_3_0-temp_data$s_53_2_0)/30/12

### Calculate PFS change rate
temp_data$pred_FI_rate_ba_re2<-(temp_data$pred_FI_re2-temp_data$pred_FI_ba)/temp_data$followup_time_ba_re2
temp_data$pred_FI_rate_ba_re3<-(temp_data$pred_FI_re3-temp_data$pred_FI_ba)/temp_data$followup_time_ba_re3
temp_data$pred_FI_rate_re2_re3<-(temp_data$pred_FI_re3-temp_data$pred_FI_re2)/temp_data$followup_time_re2_re3

### Calculate FI change rate
FI_data<-read.csv(file = "UKB_FI_XC.csv",header = T)
FI_data<-FI_data[,c("n_eid","FIsum_2","FIfrail2_2","FIfrail3_2","FI_Nmiss_2","FIsum_3","FIfrail2_3","FIfrail3_3","FI_Nmiss_3")]
temp_data<-merge(FI_data,temp_data,by="n_eid")

temp_data$FI_rate_ba_re2<-(temp_data$FIsum_2-temp_data$FIsum_0)/temp_data$followup_time_ba_re2
temp_data$FI_rate_ba_re3<-(temp_data$FIsum_3-temp_data$FIsum_0)/temp_data$followup_time_ba_re3
temp_data$FI_rate_re2_re3<-(temp_data$FIsum_3-temp_data$FIsum_2)/temp_data$followup_time_re2_re3

### Keep only complete cases
temp_data<-temp_data[complete.cases(temp_data[,c("FI_rate_ba_re2","FI_rate_ba_re3")]),]
temp_data$FI_rate_ba_re2_resid<-resid(lm(temp_data$FI_rate_ba_re2~temp_data$FIsum_0))
temp_data$FI_rate_ba_re3_resid<-resid(lm(temp_data$FI_rate_ba_re3~temp_data$FIsum_0))


##### 3 Method 2: Calculation of the change rate of PFS and FI using mixed linear model ############################
pred_FI_ba<-temp_data[,c("n_eid","sex","pred_FI_ba","FIsum_0")]
pred_FI_ba$dif_Age<-0
pred_FI_re2<-temp_data[,c("n_eid","sex","pred_FI_re2","FIsum_2")]
pred_FI_re2$dif_Age<-(temp_data$s_53_2_0-temp_data$s_53_0_0)/30/12
pred_FI_re3<-temp_data[,c("n_eid","sex","pred_FI_re3","FIsum_3")]
pred_FI_re3$dif_Age<-(temp_data$s_53_3_0-temp_data$s_53_0_0)/30/12

colnames(pred_FI_ba)[3]<-"pred_FI";colnames(pred_FI_ba)[4]<-"FIsum"
colnames(pred_FI_re2)[3]<-"pred_FI";colnames(pred_FI_re2)[4]<-"FIsum"
colnames(pred_FI_re3)[3]<-"pred_FI";colnames(pred_FI_re3)[4]<-"FIsum"

merged_data<-rbind(pred_FI_ba,pred_FI_re2,pred_FI_re3)

### Random intercept + random slope model
library(blme)
data<-merged_data

## PFS mixed model
fit_slope <- blmer(pred_FI~1+dif_Age+sex+(1+dif_Age|n_eid), control = lmerControl(optimizer = "Nelder_Mead"),data = data)
fit_summary <- summary(fit_slope)
out <- data.frame(t(fit_summary$coefficients["dif_Age",]))
out$P <- parameters::p_value(fit_slope)[2,2]
out_slope <- cbind(out,as.data.frame(performance::icc(fit_slope))[,1:2])

# Extract random effects
random_effects <- ranef(fit_slope)
random_effects_pred_FI<- as.data.frame(random_effects$n_eid)
random_effects_pred_FI$longi_slope_pred_FI<-random_effects_pred_FI$dif_Age+fit_summary$coefficients["dif_Age",1]
random_effects_pred_FI$n_eid<-rownames(random_effects_pred_FI)

## FI mixed model
data$n_eid<-as.character(data$n_eid)
fit_slope <- blmer(FIsum~1+dif_Age+sex+(1+dif_Age|n_eid), control = lmerControl(optimizer = "Nelder_Mead"),data = data)
fit_summary <- summary(fit_slope)
out <- data.frame(t(fit_summary$coefficients["dif_Age",]))
out$P <- parameters::p_value(fit_slope)[2,2]
out_slope <- cbind(out,as.data.frame(performance::icc(fit_slope))[,1:2])

# Extract random effects
random_effects <- ranef(fit_slope)
random_effects_FIsum<- as.data.frame(random_effects$n_eid)
random_effects_FIsum$longi_slope_FI<-random_effects_FIsum$dif_Age+fit_summary$coefficients["dif_Age",1]
random_effects_FIsum$n_eid<-rownames(random_effects_FIsum)

temp_data<-merge(temp_data,merge(random_effects_pred_FI[,c("n_eid","longi_slope_pred_FI")],random_effects_FIsum[,c("n_eid","longi_slope_FI")],by="n_eid"),by="n_eid")

save(temp_data,file = "./Results/results_nested_cv_longi/Results_longitudial.Rdata")

