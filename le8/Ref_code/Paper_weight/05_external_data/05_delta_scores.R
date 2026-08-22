#Author: Linsey Jackson
#Date: 12/9/2025
#Script name: delta scores
#Purpose: calculate change in risk scores


baseline_score <- predicted_scores_filtered #from filter step
final_endscores #from filter step

merged_scores_final <- baseline_score %>%
  left_join(final_endscores, by = "ID", suffix = c(".x", ".y")) %>%
  # Calculate differences for all paired columns
  mutate(across(
    ends_with(".y"),
    ~ . - get(sub(".y$", ".x", cur_column())),
    .names = "{sub('.y$', '', {.col})}"
  )) %>%
  # Keep only ID and the difference columns (no .x or .y suffix)
  select(ID, !ends_with(c(".x", ".y"))) %>%
  drop_na() #remove patients with missing scores