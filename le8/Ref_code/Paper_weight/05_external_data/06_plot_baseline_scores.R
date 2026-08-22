
#Author: Linsey Jackson
#Date: 12/9/2025
#Script name: plot baseline scores


library(tidyverse)
library(broom)

#predicted_scores_filtered df is from filtering step


# Reshape baseline data to long format
baseline_long <- predicted_scores_filtered %>%
  pivot_longer(
    cols = 2:19,  # adjust to match your outcome columns
    names_to = "Outcome",
    values_to = "Score"
  ) %>%
  mutate(
    Outcome = gsub("\\.", " ", Outcome)
  )

# Summarise: mean and standard error
baseline_summary <- baseline_long %>%
  group_by(Outcome, ARM) %>%
  summarise(
    Mean_Score = mean(Score, na.rm = TRUE),
    SE = sd(Score, na.rm = TRUE) / sqrt(sum(!is.na(Score))),
    N = sum(!is.na(Score)),
    .groups = "drop"
  )

# Define order of arms
baseline_summary <- baseline_summary %>%
  mutate(
    ARM = factor(ARM, levels = c("Placebo", "TZP 5mg", "TZP 10mg", "TZP 15mg"))
  )

write.csv(baseline_summary, "baseline_summary_stats.csv")

arm_colors <- c(
  "Placebo"     = "darkgray",
  "TZP 5mg"     = "#FBCFC8",  
  "TZP 10mg"    = "#F58E7D", 
  "TZP 15mg"    = "#E1251B"
)

# Plot baseline scores
p <- ggplot(baseline_summary, aes(x = ARM, y = Mean_Score, fill = ARM)) +
  geom_col(width = 0.6) +
  geom_errorbar(aes(ymin = Mean_Score - SE, ymax = Mean_Score + SE),
                width = 0.2, linewidth = 0.6) +
  facet_wrap(~Outcome, nrow = 3, scales = "fixed") +
  labs(
    title = "Baseline Risk Scores by Treatment Arm",
    x = "Treatment Arm",
    y = "Mean Risk Score"
  ) +
  
  theme(
    axis.text.x = element_text(angle = 30, hjust = 1),
    strip.text  = element_text(face = "bold"),
    legend.position = "bottom",
    legend.text   = element_text(size = 10),
    legend.title  = element_text(size = 12)
  ) +
  scale_fill_manual(values = arm_colors) +
  theme_bw() + 
  theme( axis.text.x = element_text(angle = 30, hjust = 1), strip.text = element_text(face = "bold"), legend.position = "bottom", legend.text = element_text(size = legend_text_size), legend.title = element_text(size = legend_title_size), legend.key.size = unit(0.5, "lines"), plot.caption = element_text(hjust = 0.5, size = 9, margin = margin(t = 10)) ) 