# ==========================================================================
# Development of Proteomic Frailty Score (PFS)
# Xueqing Jia, 2025
# Description:
#   Develop Proteomic Frailty Score (PFS) using a ten-fold iteration scheme 
#   of LASSO algorithm with ten-fold cross-validation.
# ==========================================================================

### Load R packages
library(ggplot2)
library(impute)
library(glmnet)
library(caret)
library(dplyr)
library(dplyr)
library(fdrci)
library(tableone)

########## Split the cohort to 10 sets #########
load(file = "Analysis sample_part1.Rdata")
set.seed(123)
group_sizes <- c(5050, 5050, 5050, 5050, 5050, 5050, 5050, 5050, 5050, 5056)
set_indices <- rep(1:10, times = group_sizes)
set_indices <- sample(set_indices)
sets <- split(temp_data, set_indices)
lapply(sets, head)

for (i in 1:10) {
  sets[[i]]$model<-paste0("model_",i)
}
temp_data_split<-do.call("rbind",sets)
save(temp_data_split,file = "Analysis sample_split.Rdata")

### Check independency among the sets
mBL<- CreateTableOne(vars<-c("Age_XC","sex","ethnic_v1","TD_index","Edu_XC","alco_0_XL","smok_0","hediet_0","exerc_0_XL","BMI_0","FIsum_0","FIfrail2_0","FIfrail3_0"),
                     data = temp_data_split,
                     strata = 'model',
                     factorVars=c("ethnic_v1","sex","Edu_XC","smok_0","alco_0_XL","hediet_0","exerc_0_XL","FIfrail2_0","FIfrail3_0"),
                     addOverall = TRUE)
print(mBL,showAllLevels = TRUE,smd=F)
tb_base_Mat1 <- print(mBL,showAllLevels = TRUE,quote = FALSE, noSpaces = TRUE, printToggle = FALSE)
write.csv(tb_base_Mat1, file = "baseline_characteristic_split.csv")


########## LASSO with cross-validation #########
# nested CV: imputation/standardization + feature selection (PWAS thresholding or other) + model fitting within each training fold, with evaluation strictly on held-out data.

### Define imputation function
library(impute)
impute_olink <- function(olink_data){
  olink.values <- t(olink_data)
  values.imputed <- impute::impute.knn(olink.values, k = 10, rowmax = 0.5, colmax = 0.8, maxp = 1500, rng.seed = 1714933057)$data
  df_imputed <- as.data.frame(t(values.imputed))
  return(df_imputed)
}

### Define linear regression feature selection function (training set only)
linear_reg <- function(x,reg_data){
  PCs <- paste0("n_22009_0_", 1:20, collapse = "+")
  FML <- as.formula(paste0("FIsum_0~",x,'+Age_XC+sex+ethnic_v1+Edu_XC+TD_index+BMI_0+smok_0+alco_0_XL+hediet_0+exerc_0_XL+',PCs))
  fit <- lm(FML,data = reg_data)
  GSum <- summary(fit)
  res_conti = as.data.frame(t(GSum$coefficients[x,]))
  res = res_conti |>
    mutate(protein = x)
  return(res)
}

### Define function for nested CV
nested_cv_fold <- function(trainset, testset, feature_idx, save_dir, model = i) {
  
  # Split protein and covariate data
  train_protein <- trainset[, feature_idx]
  test_protein <- testset[, feature_idx]
  train_cov_FI <- trainset[, !(names(trainset) %in% feature_idx)]
  test_cov_FI <- testset[, !(names(testset) %in% feature_idx)]
  
  # imputation
  train_imp <- impute_olink(train_protein)
  test_imp <- impute_olink(test_protein)
  
  # standardization (fit on train, apply to test)
  mm <- apply(train_imp, 2, mean)
  sd <- apply(train_imp, 2, sd)
  train_scaled <- as.data.frame(scale(train_imp, center = mm, scale = sd))
  test_scaled  <- as.data.frame(scale(test_imp,  center = mm, scale = sd))
  
  train_scaled$n_eid<-rownames(train_scaled)
  test_scaled$n_eid<-rownames(test_scaled)
  
  # Merge covariates and features (train only)
  train_final<-merge(train_cov_FI,train_scaled,by="n_eid")
  test_final<-merge(test_cov_FI,test_scaled,by="n_eid")
  
  # feature selection (train only)
  reg_data <- train_final %>%
    mutate(
      FIsum_0 = scale(FIsum_0),
      Edu_XC = as.factor(Edu_XC),
      smok_0 = as.factor(smok_0),
      alco_0_XL = as.factor(alco_0_XL)
    )
  
  UniVar<-lapply(feature_idx, linear_reg, reg_data = reg_data)
  lm_results <- bind_rows(UniVar) %>%
    select(protein, Estimate, `Std. Error`, statistic = `t value`, p_value = `Pr(>|t|)`, FI) %>%
    mutate(FDR = p.adjust(p_value, method = "fdr"),
           P_bon = p.adjust(p_value, method = "bonferroni")
    )
  
  n_sig_fdr <- nrow(filter(lm_results, FDR <= 0.05))
  n_sig_bon <- nrow(filter(lm_results, P_bon <= 0.05))
  cat("Number of proteins with FDR <= 0.05:", n_sig_fdr, "\n")
  cat("Number of proteins with Bonferroni <= 0.05:", n_sig_bon, "\n")
  
  selected_feature  <- lm_results$protein[lm_results$P_bon <= 0.05]
  
  # Save results
  save(train_final,file = paste0(save_dir, "train_imputed_model_", i, ".Rdata"))
  save(test_final,file = paste0(save_dir, "test_imputed_model_", i, ".Rdata"))
  write.csv(lm_results,file = paste0(save_dir, "LinearReg_train_model_", i, ".csv"))
  
  return(list(train = train_final, test = test_final, selected_feature = selected_feature, lm_results = lm_results))
}

### Nested CV
load(file = "Analysis sample_split.Rdata")
rownames(temp_data_split)<-temp_data_split$n_eid
colnames(temp_data_split)<-gsub("-","_",colnames(temp_data_split))

bestTune<-data.frame(lambda_1se = rep(NA, 10))
EN_list <- vector("list", 10)
model_check<-data.frame(matrix(NA, nrow = 5, ncol = 10))
rownames(model_check)<-c("R","R2","RMSE","MAE","MMSE")
pred_FI <- data.frame(n_eid = numeric(0), pred = numeric(0), FI = numeric(0), pred_exp = numeric(0))

for (i in 1:10) {
  cat("Starting fold", i, "\n")
  
  trainset<-subset(temp_data_split,temp_data_split$model!=paste0("model_",i))
  testset<-subset(temp_data_split,temp_data_split$model==paste0("model_",i))
  
  # Nested CV: preprocessing + feature selection
  nested_cv<-nested_cv_fold(
    trainset, testset, 
    feature_idx = colnames(temp_data_split)[16:2926], 
    model = i,
    save_dir = "./results_nested_cv/")
  
  train_final<-nested_cv$train
  test_final<-nested_cv$test
  selected_feature<-nested_cv$selected_feature
  
  train_final$FIsum_0 <- log(train_final$FIsum_0+0.1)
  test_final$FIsum_0 <- log(test_final$FIsum_0+0.1)
  
  # Train LASSO model
  y<-as.matrix(train_final[,c("FIsum_0")])
  x<-as.matrix(train_final[,selected_feature])
  
  set.seed(123 + i)
  mod_cv <- cv.glmnet(x, y, family = "gaussian", nfolds = 10, alpha = 1, intercept = TRUE)
  best_lambda <- mod_cv$lambda.1se
  bestTune$lambda_1se[i] <- best_lambda
  
  best_model <- glmnet(x, y, lambda = best_lambda, family = "gaussian", alpha = 1, intercept = TRUE)
  
  # Extract coefficients
  coef_vec <- as.numeric(coef(best_model, s = best_lambda))
  names(coef_vec) <- rownames(coef(best_model, s = best_lambda))
  EN_list[[i]] <- coef_vec
  
  # Evaluate performance on test set
  newx <- as.matrix(test_final[, selected_feature, drop = FALSE])
  pred <- as.numeric(predict(best_model, newx = newx, s = best_lambda))
  FI <- as.numeric(test_final$FIsum_0)
  
  model_check[,i]<-c(cor(pred, FI),R2(pred, FI),RMSE(pred, FI),MAE(pred, FI),RMSE(pred, FI)/mean(FI))
  
  fold_df <- data.frame(n_eid = test_final$n_eid, pred = pred, FI = FI, pred_exp = exp(pred) - 0.1)
  pred_FI <- rbind(pred_FI, fold_df)
}

# Save coefficient list and results
save_dir = "./results_nested_cv/"
save(EN_list, file = file.path(save_dir, "EN_coeffs_list.Rdata"))
write.csv(bestTune, file = file.path(save_dir, "bestTune.csv"), row.names = FALSE)
write.csv(model_check, file = file.path(save_dir, "model_check.csv"))
write.csv(pred_FI, file = file.path(save_dir, "pred_vs_true_allfolds.csv"), row.names = FALSE)


########## Visualization of results ############
load(file = "dataforplot.Rdata")

### Scatter plot: FI vs. PFS
cor_result <- cor.test(data$FIsum_0, data$FI_sum_0_pre)
pearson_r <- round(cor_result$estimate, 3)
p_value <- format.pval(cor_result$p.value)

p <- ggplot(data, aes(x = FIsum_0, y = FI_sum_0_pre)) +
  geom_point(alpha = 0.6, color = "#478ecc") +
  geom_smooth(method = "lm", se = TRUE, color = "gray20", fill = "gray20", alpha = 0.3) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray") +
  labs(title = "Actual FI vs. Predicted FI",
       x = "FI",
       y = "Proteomic FI") +
  theme_minimal() +
  annotate("text", x = 0.1, y = 0.35,
           label = paste("Pearson's r =", pearson_r, "\n",
                         "P value ", p_value),
           size = 5, vjust = -1) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 15),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 12),
    axis.line = element_line(color = "black"),
    axis.ticks.length = unit(0.1, "cm"),
    axis.ticks.x = element_line(),
    axis.ticks.y = element_line())

p

### Violin plot by sex and age
# Stratify by sex
data$sex<-factor(data$sex,levels=c(0,1),labels=c("Male","Female"))

library("ggpubr")
P<- ggplot(data, aes(x = sex, y = FI_sum_0_pre, fill=sex)) + 
  rotate_x_text(angle = 45) +
  # scale_x_discrete(labels = c("0" = "Female","1" = "Male")) +
  geom_violin(trim=FALSE,color="white",alpha=0.4) + 
  geom_boxplot(width=0.2,position=position_dodge(0.9),alpha=0.6, outlier.shape = NA)+
  scale_fill_manual(values = c("#404080", "#d56516")) +
  theme_bw()+ theme(legend.position = "top") +
  theme(axis.text.x=element_text(colour="black",family="Times",size=10),
        axis.text.y=element_text(family="Times",size=10,face="plain"),
        axis.title.x=element_text(family="Times",size = 12,face="plain"),
        axis.title.y=element_text(family="Times",size = 12,face="plain"),
        panel.border = element_blank(),axis.line = element_line(colour = "black",size=1), 
        title=element_text(family="Times",size = 14,face="plain")
  )+ 
  ylab("PFS")+xlab("sex")
P + stat_compare_means(method = "t.test",
                       aes(label = "p.format"),
                       size = 3)

# Stratify by age
data$age<-as.factor(ifelse(data$Age_XC<60,"Young","Old"))
data$age<-factor(data$age,levels=c("Young","Old"))

P<- ggplot(data, aes(x = age, y = FI_sum_0_pre, fill=age)) + 
  rotate_x_text(angle = 45) +
  # scale_x_discrete(labels = c("0" = "Female","1" = "Male")) +
  geom_violin(trim=FALSE,color="white",alpha=0.4) + 
  geom_boxplot(width=0.2,position=position_dodge(0.9),alpha=0.6, outlier.shape = NA)+
  scale_fill_manual(values = c("#404080", "#d56516")) +
  theme_bw()+ theme(legend.position = "top") +
  theme(axis.text.x=element_text(colour="black",family="Times",size=10),
        axis.text.y=element_text(family="Times",size=10,face="plain"),
        axis.title.x=element_text(family="Times",size = 12,face="plain"),
        axis.title.y=element_text(family="Times",size = 12,face="plain"),
        panel.border = element_blank(),axis.line = element_line(colour = "black",size=1), 
        title=element_text(family="Times",size = 14,face="plain")
  )+ 
  ylab("PFS")+xlab("Age")
P + stat_compare_means(method = "t.test",
                       aes(label = "p.format"),
                       size = 3)


########## Extract mean and standard deviation for each fold #######
### Define imputation function
impute_olink <- function(olink_data){
  olink.values <- t(olink_data)
  values.imputed <- impute::impute.knn(olink.values, k = 10, rowmax = 0.5, colmax = 0.8, maxp = 1500, rng.seed = 1714933057)$data
  df_imputed <- as.data.frame(t(values.imputed))
  return(df_imputed)
}

### Define function to compute mean and SD in each training set
nested_cv_mean_sd <- function(trainset, feature_idx) {
  
  # extract protein data
  train_protein <- trainset[, feature_idx]
  
  # imputation
  train_imp <- impute_olink(train_protein)
  
  # compute mean and standard deviation
  means <- apply(train_imp, 2, mean)
  sds   <- apply(train_imp, 2, sd)
  
  mm <- data.frame(mean = means, row.names = names(means))
  sd <- data.frame(sd = sds, row.names = names(sds))
  
  return(list(mean = mm, sd = sd))
}

### Main loop
load(file = "Analysis sample_split.Rdata")
rownames(temp_data_split)<-temp_data_split$n_eid
colnames(temp_data_split)<-gsub("-","_",colnames(temp_data_split))

feature_idx = colnames(temp_data_split)[16:2926]

mm_sd_list <- vector("list", 10)

for (i in 1:10) {
  cat("Starting fold", i, "\n")
  trainset <- subset(temp_data_split, model != paste0("model_", i))
  mm_sd_list[[i]] <- nested_cv_mean_sd(trainset, feature_idx)
}

# Save mean and SD results
save(mm_sd_list, file = "mm_sd_list.Rdata")
