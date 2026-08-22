# ==========================================================================
# Function used for calculating PFS in TwinGene cohort
# Xueqing Jia, 2025
# ==========================================================================


### Function used for calculating PFS in TwinGene

PFS_calculator_fun <- function(data_demo, olink_proteins, lasso_list, mm_sd_list){
  
  # data_demo: external dataset (rows = individuals, columns = proteins)
  # olink_proteins: vector of protein names used in the model
  # lasso_list: list of 10 trained lasso models (each a named coefficient vector)
  # mm_sd_list: list of corresponding mean and SD values for standardization
  
  if (is.null(rownames(data_demo))) stop("data_demo must have rownames.")
  
  pred_list <- vector("list", length(lasso_list)) # store predictions from each model
  
  for (i in seq_along(lasso_list)) {
    cat("Processing fold", i, "\n")
    
    # Extract protein data for the current model
    data_tem <- data_demo[, olink_proteins, drop = FALSE]
    
    # Retrieve mean and SD for standardization
    mm <- setNames(mm_sd_list[[i]]$mean$mean, rownames(mm_sd_list[[i]]$mean))
    sd <- setNames(mm_sd_list[[i]]$sd$sd, rownames(mm_sd_list[[i]]$sd))
    mm <- mm[olink_proteins]
    sd <- sd[olink_proteins]
    
    # Standardize external data using training parameters
    data_scaled  <- as.data.frame(scale(as.matrix(data_tem[, olink_proteins]),  center = mm, scale = sd))
    rownames(data_scaled) <- rownames(data_tem)
    
    # Extract non-zero coefficients from the lasso model
    var_list <- lasso_list[[i]][lasso_list[[i]] != 0]
    intercept <- if ("(Intercept)" %in% names(var_list)) var_list["(Intercept)"] else 0
    coef_vec <- var_list[names(var_list) != "(Intercept)"]
    
    cat("Number of non-zero proteins (excluding intercept):", length(coef_vec), "\n")
    
    # Select the corresponding protein columns for prediction
    data_matched <- data_scaled[, names(coef_vec), drop = FALSE]
    
    # Linear predictor: Xβ + intercept
    pred_scores <- as.matrix(data_matched) %*% as.numeric(coef_vec) + as.numeric(intercept)
    
    # Back-transform the predictions (exp(·) − 0.1, as in the original training step)
    result_df <- data.frame(
      n_eid = rownames(data_matched),
      PFS = as.numeric(exp(pred_scores) - 0.1),
      stringsAsFactors = FALSE
    )
    colnames(result_df)[2]<-paste0("PFS_model_",i)
    
    pred_list[[i]] <- result_df
  }
  
  # Merge predictions from all 10 models by individual ID
  final_pred <- Reduce(function(x, y) merge(x, y, by = "n_eid"), pred_list)
  
  # Average PFS predictions across the 10 models
  final_pred$PFS_output <- rowMeans(final_pred[, -1])
  final_pred <- final_pred[, c("n_eid", "PFS_output")]
  return(final_pred)
}

#
results <- PFS_calculator_fun(data_demo, olink_proteins, lasso_list, mm_sd_list)
