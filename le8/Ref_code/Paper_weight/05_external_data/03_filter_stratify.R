
#Author: Linsey Jackson
#Date: 12/9/2025
#Script name: filter & stratify
#Purpose: filter patients by those with complete data (scores at baseline and week; change in endpoint at week 72)
#then stratify patients by baseline risk scores 

library(dplyr)   
library(tidyr)   
library(tibble)  

#final_endpoints - change in endpoints 
#final_endscores - risk scores at week 72
#predicted_scores - risk scores at baseline

#filter for patients with complete data 
final_endpoints<-read.csv("final_endpoints.csv")
final_endpoints <- final_endpoints %>% 
  rename(ID = USUBJID)
final_endscores<-read.csv("scores_wk72.csv")
predicted_scores<-read.csv("")

#find intersection of IDs in both week 72 and endpoints
common_ids <- intersect(final_endpoints$ID, final_endscores$ID)
#filter baseline scores to keep only complete data
predicted_scores_filtered<- predicted_scores %>% 
  filter(ID %in% common_ids)


#stratify patients by baseline risk score
library(dplyr)

# Identify score columns explicitly
score_cols <- setdiff(names(predicted_scores_filtered), "ID")

# --- Parameters ---
n_groups <- 4  # number of strata/quantiles desired

# --- Helper: compute robust global breaks (quantiles) per outcome ---
compute_breaks <- function(x, n_groups = 4) {
  # Quantiles across ALL subjects (all arms combined)
  qs <- quantile(x, probs = seq(0, 1, length.out = n_groups + 1),
                 na.rm = TRUE, type = 7)
  br <- unique(as.numeric(qs))  # remove duplicates to avoid cut() errors
  
  # If everything is identical, you'll only get one value.
  # In that case, cut() will produce NA bins (no variation).
  # We'll return the unique vector and let labeling adapt to the available bins.
  br
}

# Compute breaks for each score column globally
breaks_list <- lapply(merged_df[score_cols], compute_breaks, n_groups = n_groups)
names(breaks_list) <- score_cols

# --- Build a tidy table of cutoffs used
library(tidyr)
library(tibble)
cutoffs_long <- do.call(rbind, lapply(names(breaks_list), function(nm) {
  br <- breaks_list[[nm]]
  k <- length(br)
  tibble(
    outcome = nm,
    boundary_index = seq_len(k),
    prob = if (k > 1) seq(0, 1, length.out = k) else 0,
    value = br,
    # For convenience: the upper bound label of each bin (Q1, Q2, ...)
    bin_upper_label = c(NA, paste0("Q", seq_len(max(1, k - 1))))
  )
}))


# --- Assign strata globally using precomputed breaks ---
result_stratified <- predicted_scores_filtered
for (nm in score_cols) {
  br <- breaks_list[[nm]]
  # Labels must match the number of bins (length(br) - 1). If no variation, br length==1.
  labs <- if (length(br) >= 2) paste0("Q", seq_len(length(br) - 1)) else NA
  result_stratified[[paste0(nm, "_group")]] <- cut(
    result_stratified[[nm]],
    breaks = br,
    include.lowest = TRUE,
    labels = labs,
    right = TRUE
  )
}

#join endpoints with stratified baseline scores 
df_joined <- left_join(results_stratified, final_endpoints, by = c("ID" = "USUBJID"))