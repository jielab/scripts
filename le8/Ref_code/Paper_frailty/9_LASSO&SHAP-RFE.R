# ==================================================================================
# LASSO + SHAP-RFE: Evaluate LASSO selection stability + select simplified PFS model
# Xueqing Jia, 2025
# ==================================================================================

### Load R packages
library(glmnet)
library(caret)
library(parallel)
library(doParallel)
library(foreach)
library(iml)
library(dplyr)
library(tidyverse)
library(ggplot2)

######## 1 LASSO + SHAP-RFE procedure ########

### Proteins retained in all 10 models ###
load("coef_stats.Rdata")
core_features <- coef_stats %>% filter(freq == 10 & feature != "(Intercept)") %>% pull(feature)

### Define compute_foldaware_shap_rfe function ###
compute_foldaware_shap_rfe <- function(to_remove, save_dir, iter, fold_n = 10, core_features,
                                       n_cores = 14, n_shap_samples){
  # Parallel cores
  cl <- makeCluster(n_cores)
  registerDoParallel(cl)
  
  # Parallel computation
  res_list <- foreach(f = 1:fold_n, .packages = c("glmnet","caret","iml","dplyr")) %dopar% {
    
    cat("  [iter", iter, "] processing fold", f, "\n")
    
    # Load preprocessed Nested CV fold datasets and feature selection results
    load(file = paste0("./Results/", "train_imputed_model_", f, ".Rdata"))
    load(file = paste0("./Results/", "test_imputed_model_", f, ".Rdata"))
    
    # Remove already deleted proteins
    selected_feature <- setdiff(core_features, to_remove)
    
    cat("Selected features after remove:", length(selected_feature), "features\n")
    
    # Preprocess FIsum_0
    train_final$FIsum_0 <- log(train_final$FIsum_0+0.1)
    test_final$FIsum_0 <- log(test_final$FIsum_0+0.1)
    
    y_train <- as.matrix(train_final$FIsum_0)
    x_train <- as.matrix(train_final[, selected_feature, drop = FALSE])
    y_test  <- as.numeric(test_final$FIsum_0)
    x_test  <- as.matrix(test_final[, selected_feature, drop = FALSE])
    
    # LASSO CV
    set.seed(123 + f)
    mod_cv <- cv.glmnet(x_train, y_train, family = "gaussian", nfolds = 10, alpha = 1, intercept = TRUE)
    best_lambda <- mod_cv$lambda.1se
    
    best_model <- glmnet(x_train, y_train, lambda = best_lambda, family = "gaussian", alpha = 1, intercept = TRUE)
    
    # Extract coefficients
    coef_vec <- as.numeric(coef(best_model, s = best_lambda))
    names(coef_vec) <- rownames(coef(best_model, s = best_lambda))
    
    # Test set prediction R²
    pred <- as.numeric(predict(best_model, newx = x_test, s = best_lambda))
    R2_val <- R2(pred, y_test)
    
    model_metrics <- c(cor(pred, y_test),R2(pred, y_test),RMSE(pred, y_test),MAE(pred, y_test),RMSE(pred, y_test)/mean(y_test))
    
    # SHAP calculation (iml)
    predict_fun <- function(model, newdata){
      as.numeric(predict(model, newx = as.matrix(newdata), s = best_lambda))
    }
    
    pred_obj <- Predictor$new(
      model = best_model,
      data = as.data.frame(x_train),
      y = y_train,
      predict.function = predict_fun,
      type = "regression"
    )
    
    n_samp <- min(n_shap_samples, nrow(x_train))
    set.seed(2000 + iter + f)
    samp_idx <- sample(seq_len(nrow(x_train)), n_samp)
    
    shap_mat <- sapply(samp_idx, function(ii){
      shap_i <- Shapley$new(pred_obj, x.interest = as.data.frame(x_train[ii, , drop = FALSE]))
      phi_vec <- shap_i$results$phi[match(selected_feature, shap_i$results$feature)]
      phi_vec[is.na(phi_vec)] <- 0
      return(phi_vec)
    })
    
    if(is.null(dim(shap_mat))) 
      shap_mat <- matrix(shap_mat, nrow = length(selected_feature), ncol = 1)
    rownames(shap_mat) <- selected_feature
    
    shap_mean <- rowMeans(abs(shap_mat), na.rm = TRUE)
    
    return(list(
      fold = f,
      shap_mean = shap_mean,
      R2 = R2_val,
      coef_vec = coef_vec,
      best_lambda = best_lambda,
      model_metrics = model_metrics ,
      pred_df = data.frame(
        n_eid = test_final$n_eid,
        pred = pred,
        FI = y_test,
        pred_exp = exp(pred) - 0.1
      )
    ))
  }
  
  stopCluster(cl)
  
  # Organize results
  EN_list <- lapply(res_list, `[[`, "coef_vec")
  bestTune <- data.frame(lambda_1se = sapply(res_list, `[[`, "best_lambda"))
  model_check <- do.call(cbind, lapply(res_list, `[[`, "model_metrics"))
  rownames(model_check) <- c("R","R2","RMSE","MAE","MMSE")
  pred_FI <- do.call(rbind, lapply(res_list, `[[`, "pred_df"))
  
  results <- list(EN_list = EN_list, model_check = model_check, bestTune = bestTune, pred_FI = pred_FI)
  save(results, file = file.path(save_dir,paste0("results_iter_recur_shap_", iter, ".Rdata")))
  
  
  # Summarize SHAP results across folds
  fold_shap_list <- lapply(res_list, function(x) x$shap_mean)
  all_feats <- sort(unique(unlist(lapply(fold_shap_list, names))))
  shap_mat_aligned <- sapply(fold_shap_list, function(v){
    tmp <- rep(NA, length(all_feats))
    names(tmp) <- all_feats
    tmp[names(v)] <- v
    tmp
  })
  rownames(shap_mat_aligned) <- all_feats
  
  # Step 1: Global average SHAP, find minimum protein
  shap_mean_global <- rowMeans(shap_mat_aligned, na.rm = TRUE)
  min_val <- min(shap_mean_global, na.rm = TRUE)
  candidate_feats <- names(shap_mean_global)[shap_mean_global == min_val]
  
  if(length(candidate_feats) == 1){
    feature_to_remove <- candidate_feats
  } else {
    # Step 2: minimum vote across folds
    min_per_fold <- apply(shap_mat_aligned, 2, function(col){
      col2 <- col
      col2[is.na(col2)] <- Inf
      names(which.min(col2))
    })
    vote_tab <- table(min_per_fold)
    vote_tab <- vote_tab[candidate_feats]
    max_vote <- max(vote_tab, na.rm = TRUE)
    top_candidates <- names(vote_tab)[vote_tab == max_vote]
    
    if(length(top_candidates) == 1){
      feature_to_remove <- top_candidates
    } else {
      # Step 3: tie-break → random selection
      set.seed(1000 + iter)
      feature_to_remove <- sample(top_candidates, 1)
    }
  }
  
  # Average R²
  R2_vals <- sapply(res_list, function(x) x$R2)
  R2_mean <- mean(R2_vals, na.rm = TRUE)
  
  ci <- function(x) {
    n <- length(x)
    se <- sd(x, na.rm = TRUE)/sqrt(n)
    t_val <- qt(0.975, df=n-1)
    c(lower=mean(x, na.rm=TRUE)-t_val*se, upper=mean(x, na.rm=TRUE)+t_val*se)
  }
  R2_ci <- ci(R2_vals)
  
  return(list(feature_to_remove = feature_to_remove,
              R2_mean = R2_mean, R2_lower = R2_ci["lower"], R2_upper = R2_ci["upper"],
              shap_mat_aligned = shap_mat_aligned,
              res_list = res_list))
}


### Iterative removal ###
results_iter <- list()
cur_removed <- character(0)
save_dir <- "./Results/"

for(iter in 231:270){
  remaining_feats <- setdiff(core_features, cur_removed)
  
  if(length(remaining_feats) <= 5){
    cat("Reached 5 proteins, stopping iteration.\n")
    break
  }
  
  cat("=== Outer iteration", iter, "===\n")
  
  out <- compute_foldaware_shap_rfe(to_remove = cur_removed, save_dir = save_dir,
                                    iter = iter, fold_n = 10, core_features = core_features,
                                    n_shap_samples = n_shap_samples, n_cores = 14)
  
  feature_to_remove <- out$feature_to_remove
  cat("Feature to remove:", feature_to_remove, "\n")
  
  cur_removed <- c(cur_removed, feature_to_remove)
  results_iter[[as.character(iter)]] <- list(
    iter = iter,
    removed = out$feature_to_remove,
    R2_mean = out$R2_mean,
    R2_lower = out$R2_lower, R2_upper = out$R2_upper,
    shap_mat = out$shap_mat_aligned
  )
}

save(results_iter, file = file.path(save_dir,"iterative_shap_removal_summary.Rdata"))


######## 2 Organize results and visualize ############

### Data preparation ###
load(file = "D:/Projects_ZJU/UKB_proteomic&frailty/修回/Results/iterative_shap_removal_summary.Rdata")

results_df <- results_iter %>%
  map_dfr(~ data.frame(
    iter = .x$iter,
    protein_number = nrow(.x$shap_mat),
    feature_removed = .x$removed,
    R2_mean = .x$R2_mean,
    R2_lower = as.numeric(.x$R2_lower),
    R2_upper = as.numeric(.x$R2_upper),
    stringsAsFactors = FALSE
  ))

head(results_df)

### Visualization ###
y_breaks <- seq(0.1, 0.3, by = 0.025)
y_labels <- paste0(round(y_breaks, 3))
x_breaks <- seq(0, 274, by = 20)

p <- ggplot(results_df, aes(x = protein_number, y = R2_mean)) +
  geom_line(color = "black", size = 1) +
  
  geom_ribbon(aes(ymin = R2_lower, ymax = R2_upper), 
              fill = "gray80", alpha = 0.5) +
  
  labs(x = "Number of proteins", y = "Out-of-sample R²") +

  scale_x_continuous(breaks = x_breaks, limits = c(0, 274), labels = x_breaks, expand = c(0, 0)) +
  scale_y_continuous(breaks = y_breaks, limits = c(0.15, 0.28), labels = y_labels) +
  
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

