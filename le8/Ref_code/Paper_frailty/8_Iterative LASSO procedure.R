# ==========================================================================
# Iterative LASSO procedure: Assess LASSO selection stability
# Xueqing Jia, 2025
# ==========================================================================

### Load R packages
library(parallel)
library(doParallel)
library(foreach)
library(glmnet)
library(caret)
library(tidyverse)
library(ggplot2)

######## 1 Iterative LASSO procedure ########

### Define LASSO procedure function ###
compute_cv_R2_with_glmnet <- function(to_remove, save_dir, iter, n_cores = 14){
  
  # Parallel cores
  cl <- makeCluster(n_cores)
  registerDoParallel(cl)
  
  # Parallel computation
  res_list <- foreach(i = 1:10, .packages = c("glmnet","caret")) %dopar% {
    
    cat("Starting fold", i, "\n")
    
    # Load preprocessed nested CV fold data and feature selection results
    load(file = paste0(save_dir, "train_imputed_model_", i, ".Rdata"))
    load(file = paste0(save_dir, "test_imputed_model_", i, ".Rdata"))
    
    lm_results <- read.csv(file = paste0(save_dir, "LinearReg_train_model_", i, ".csv"))
    selected_feature <- lm_results$protein[lm_results$P_bon <= 0.05]
    
    cat("selected_feature", selected_feature, "\n")
    
    # Remove top proteins
    selected_feature <- selected_feature[! selected_feature %in% to_remove]
    
    cat("selected_feature_after_remove", selected_feature, "\n")
    
    # Preprocess FIsum_0
    train_final$FIsum_0 <- log(train_final$FIsum_0+0.1)
    test_final$FIsum_0 <- log(test_final$FIsum_0+0.1)
    
    # Train LASSO model
    y<-as.matrix(train_final[,c("FIsum_0")])
    x<-as.matrix(train_final[,selected_feature])
    
    set.seed(123 + i)
    mod_cv <- cv.glmnet(x, y, family = "gaussian", nfolds = 10, alpha = 1, intercept = TRUE)
    best_lambda <- mod_cv$lambda.1se
    
    best_model <- glmnet(x, y, lambda = best_lambda, family = "gaussian", alpha = 1, intercept = TRUE)
    
    # Extract coefficients
    coef_vec <- as.numeric(coef(best_model, s = best_lambda))
    names(coef_vec) <- rownames(coef(best_model, s = best_lambda))
    
    # Evaluate performance on test set
    newx <- as.matrix(test_final[, selected_feature, drop = FALSE])
    pred <- as.numeric(predict(best_model, newx = newx, s = best_lambda))
    FI <- as.numeric(test_final$FIsum_0)
    
    model_metrics <- c(cor(pred, FI),R2(pred, FI),RMSE(pred, FI),MAE(pred, FI),RMSE(pred, FI)/mean(FI))
    
    list(
      coef_vec = coef_vec,
      best_lambda = best_lambda,
      model_metrics = model_metrics ,
      pred_df = data.frame(
        n_eid = test_final$n_eid,
        pred = pred,
        FI = FI,
        pred_exp = exp(pred) - 0.1
      )
    )
  }
  
  stopCluster(cl)
  
  # output results
  EN_list <- lapply(res_list, `[[`, "coef_vec")
  bestTune <- data.frame(lambda_1se = sapply(res_list, `[[`, "best_lambda"))
  model_check <- do.call(cbind, lapply(res_list, `[[`, "model_metrics"))
  rownames(model_check) <- c("R","R2","RMSE","MAE","MMSE")
  pred_FI <- do.call(rbind, lapply(res_list, `[[`, "pred_df"))
  
  results <- list(EN_list = EN_list, model_check = model_check, bestTune = bestTune, pred_FI = pred_FI)
  
  save(results, file = file.path(save_dir,paste0("results_iter_recur_", iter, ".Rdata")))
  
  return(results)
}


### Iterative removal procedure ###

results_iter <- list()
cur_removed <- character(0)
save_dir = "./Results/"

for (iter in 1:500) {
  
  cat("=== Iteration", iter, "===\n")
  
  # 1) retrain model
  cv_results <- compute_cv_R2_with_glmnet(to_remove = cur_removed, save_dir = save_dir, iter = iter)
  
  # 2) Align coefficients from cv_results$EN_list and calculate mean(|beta|)
  EN_list <- cv_results$EN_list
  all_feats <- unique(unlist(lapply(EN_list, names)))
  
  EN_aligned <- lapply(EN_list, function(cf) {
    cf_full <- setNames(rep(NA, length(all_feats)), all_feats)
    cf_full[names(cf)] <- cf
    cf_full
  })
  
  coef_mat <- do.call(cbind, EN_aligned)
  rownames(coef_mat) <- all_feats
  
  # 3) Decide which feature to remove (i.e., remove protein appearing in all folds with largest mean_abs_beta)
  stable_nonzero <- all_feats[apply(!is.na(coef_mat) & coef_mat != 0, 1, all)]
  
  if(length(stable_nonzero) == 0){
    cat("No stable non-zero protein left!\n")
    break
  }
  
  coef_mat_stable <- coef_mat[stable_nonzero, , drop = FALSE]
  mean_abs_beta <- rowMeans(abs(coef_mat_stable), na.rm = TRUE)
  mean_abs_beta <- mean_abs_beta[names(mean_abs_beta) != "(Intercept)"]
  
  to_remove_now <- names(which.max(mean_abs_beta))
  cat("Removing:", to_remove_now, "\n")
  
  # 4) Record performance for this iteration (mean R2 + CI)
  # mean
  R_folds <- as.numeric(cv_results$model_check[1,])
  R2_folds <- as.numeric(cv_results$model_check[2,])
  
  R_mean <- mean(R_folds, na.rm = TRUE)
  R2_mean <- mean(R2_folds, na.rm = TRUE)
  
  # 95% CI
  ci <- function(x) {
    n <- length(x)
    se <- sd(x, na.rm = TRUE)/sqrt(n)
    t_val <- qt(0.975, df=n-1)
    c(mean=mean(x, na.rm=TRUE), lower=mean(x, na.rm=TRUE)-t_val*se, upper=mean(x, na.rm=TRUE)+t_val*se)
  }
  R_ci <- ci(R_folds)
  R2_ci <- ci(R2_folds)
  
  # Summarize results
  results_iter[[iter]] <- list(
    iter = iter,
    feature_removed = to_remove_now,
    R_mean = R_ci["mean"], R_lower = R_ci["lower"], R_upper = R_ci["upper"],
    R2_mean = R2_ci["mean"], R2_lower = R2_ci["lower"], R2_upper = R2_ci["upper"]
  )
  
  # 5) Update removed feature list
  cur_removed <- c(cur_removed, to_remove_now)
}

save(results_iter, file = file.path(save_dir,"iterative_removal_summary.Rdata"))


######## 2 Organize results and visualize ############

### Data preparation ###
load(file = "iterative_removal_summary.Rdata")

results_df <- results_iter %>%
  map_dfr(~ data.frame(
    iter = .x$iter,
    feature_removed = .x$feature_removed,
    R_mean = .x$R_mean,
    R_lower = .x$R_lower,
    R_upper = .x$R_upper,
    R2_mean = .x$R2_mean,
    R2_lower = .x$R2_lower,
    R2_upper = .x$R2_upper,
    stringsAsFactors = FALSE
  ))
head(results_df)

### Visualization  ###
# Step 1: Create a new column marking whether the protein is among core features retained in all folds
core_features <- read.csv(file = "core_features.csv")
core_features <-core_features$feature[-1]
results_df$feature_removed_flag <- ifelse(results_df$feature_removed %in% core_features, 1, 0)

# Step 2: Plot
y_breaks <- seq(0, 0.3, by = 0.05)
y_labels <- paste0(round(y_breaks, 2))
x_breaks <- seq(0, 500, by = 50)
p <- ggplot(results_df, aes(x = iter, y = R2_mean)) +
  geom_line(color = "black", size = 1) +
  
  geom_ribbon(aes(ymin = R2_lower, ymax = R2_upper), 
              fill = "gray80", alpha = 0.5) +
  
  geom_bar(data = subset(results_df, feature_removed_flag == 1),
           aes(x = iter), stat = "identity", width = 0.7,
           fill = "#989cc8", alpha = 0.3) +
  labs(x = "Iteration number", y = "Out-of-sample R²") +
  
  scale_x_continuous(breaks = x_breaks, limits = c(0, 500), labels = x_breaks, expand = c(0, 0)) +
  scale_y_continuous(breaks = y_breaks, limits = c(0, 0.3), labels = y_labels, expand = c(0, 0)) +
  theme_minimal() +
  theme(
    panel.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black"),
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 14),
    axis.ticks = element_line(color = "black", size = 0.5),
    axis.ticks.length = unit(4, "pt")
  )

print(p)

