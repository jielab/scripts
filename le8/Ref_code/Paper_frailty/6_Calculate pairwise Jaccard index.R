# ==========================================================================
# Calculate Pairwise Jaccard Index
# Xueqing Jia, 2025
# ==========================================================================

### Load R packages
library(dplyr)
library(tidyr)
library(proxy)
library(corrplot)


### Calculate Jaccard index
load("EN_coeffs_list.Rdata")

# Get proteins that appear in at least one model
all_features <- unique(unlist(lapply(EN_list, names)))

# Align coefficient vectors
EN_aligned <- lapply(EN_list, function(cf) {
  cf_full <- setNames(rep(NA, length(all_features)), all_features)
  cf_full[names(cf)] <- cf
  return(cf_full)
})

# Build matrix: whether retained in model (0/1/NA) + calculate frequency of each protein across 10 models
sel_mat <- do.call(cbind, lapply(EN_aligned, function(cf) as.numeric(cf != 0)))
rownames(sel_mat) <- all_features
colnames(sel_mat) <- paste0("model_", seq_along(EN_list))
sel_mat<-sel_mat[-1,]

# Compute Jaccard similarity between folds
jaccard_sim <- proxy::simil(sel_mat, method = "Jaccard", by_rows = FALSE)
jaccard_matrix <- as.matrix(jaccard_sim)
print(jaccard_matrix)

# Calculate mean Jaccard similarity
jaccard_values <- jaccard_matrix[lower.tri(jaccard_matrix)]
mean_jaccard <- mean(jaccard_values)
mean_jaccard


### Visualization-Heatmap ######
palette_5_blue <- c("#C9D6E0", "#C4D9EC", "#8FB8DE", "#6CA3D4")
morandi_colors <- colorRampPalette(palette_5_blue)(100)

corrplot(
  jaccard_matrix,
  method = "circle", 
  type = "lower",
  tl.col = "black",
  tl.cex = 0.8,
  tl.srt = 45,
  diag = FALSE,
  is.corr = FALSE,
  addCoef.col = "grey20",
  number.cex = 0.7,
  
  col = morandi_colors,
  col.lim = c(0.7,0.8),
)


