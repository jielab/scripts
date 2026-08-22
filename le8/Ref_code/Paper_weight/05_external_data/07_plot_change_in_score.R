
#Author: Linsey Jackson
#Date: 12/9/2025
#Script name: plot change in score
#Purpose: calculate one-way ANOVA across treatment arms
#unadjusted pairwise contrasts of the group means
#plot results


library(dplyr)
library(tidyr)
library(ggplot2)
library(tidyverse)
library(broom)
library(emmeans)


df<-merged_scores_final #from delta score step


# Reshape from wide to long
long_df <- df %>%
  pivot_longer(
    cols = 3:20,  # adjust to match your outcome columns
    names_to = "Outcome",
    values_to = "Change"
  ) %>%
  mutate(
    Outcome = gsub("\\.", " ", Outcome)  # replace periods with spaces
  )


# Perform pairwise comparisons: Placebo vs each treatment
pairwise_results <- long_df %>%
  group_by(Outcome) %>%
  do({
    model <- aov(Change ~ ARM, data = .)
    emm <- emmeans(model, ~ ARM)
    pairs(emm, adjust = "none") %>%
      as.data.frame() %>%
      filter(grepl("Placebo", contrast))
  }) %>%
  ungroup() %>%
  mutate(
    # Extract treatment arm from contrast
    treatment_arm = case_when(
      grepl("TZP 5mg", contrast) ~ "TZP 5mg",
      grepl("TZP 10mg", contrast) ~ "TZP 10mg",
      grepl("TZP 15mg", contrast) ~ "TZP 15mg"
    ),
    # Assign significance symbols
    significance = case_when(
      p.value < 0.001 ~ "***",
      p.value < 0.01  ~ "**",
      p.value < 0.05  ~ "*",
      TRUE ~ ""
    )
  )


########### all comparisons 

pairwise_results <- long_df %>%
  group_by(Outcome) %>%
  do({
    # Fit ANOVA for each outcome
    model <- aov(Change ~ ARM, data = .)
    
    # Get estimated marginal means
    emm <- emmeans(model, ~ ARM)
    
    
    pairs(emm, adjust = "none") %>%
      as.data.frame() %>%
      # Split the contrast into two arms
      mutate(
        arm1 = sub(" -.*", "", contrast),
        arm2 = sub(".*- ", "", contrast)
      )
  }) %>%
  ungroup()



# Summarise: mean, standard error, and N
summary_df <- long_df %>%
  group_by(Outcome, ARM) %>%
  summarise(
    Mean_Change = mean(Change, na.rm = TRUE),
    SE = sd(Change, na.rm = TRUE) / sqrt(sum(!is.na(Change))),
    N = sum(!is.na(Change)),
    .groups = "drop"
  )


write.csv(summary_df, "changeinscore_summary_stats.csv")
write.csv(pairwise_results, "deltascore_pairwise_unadj.csv")


# Define order of arms
summary_df <- summary_df %>%
  mutate(
    ARM = factor(ARM, levels = c("Placebo", "TZP 5mg", "TZP 10mg", "TZP 15mg"))
  )

arm_colors <- c(
  "Placebo"     = "darkgray",
  "TZP 5mg"     = "#FBCFC8",  
  "TZP 10mg"    = "#F58E7D", 
  "TZP 15mg"    = "#E1251B"
)

# Legend sizes
legend_text_size  <- 10
legend_title_size <- 12




#####   final plot

library(ggplot2)
library(dplyr)

# Step 1: Precompute label positions
summary_df <- summary_df %>%
  mutate(
    # Position numeric labels just outside the error bars
    label_y = ifelse(
      Mean_Change >= 0,
      Mean_Change + SE + 0.05 * max(abs(Mean_Change + SE), na.rm = TRUE),
      Mean_Change - SE - 0.05 * max(abs(Mean_Change + SE), na.rm = TRUE)
    ),
    # Position significance asterisks slightly above/below the numeric labels
    asterisk_y = ifelse(
      Mean_Change >= 0,
      label_y + 0.05 * max(abs(Mean_Change + SE), na.rm = TRUE),
      label_y - 0.05 * max(abs(Mean_Change + SE), na.rm = TRUE)
    )
  )

# Calculate the minimum (most negative) change across all 4 arms for each outcome
outcome_order <- summary_df %>%
  group_by(Outcome) %>%
  summarise(min_change = min(Mean_Change)) %>%  # finds the most negative value among the 4 arms
  arrange(min_change) %>%  # sorts outcomes by this minimum value
  pull(Outcome)

# Convert Outcome to factor with this ordering
summary_df <- summary_df %>%
  mutate(Outcome = factor(Outcome, levels = outcome_order))

# Step 2: Create plot
p <- ggplot(summary_df, aes(x = ARM, y = Mean_Change, fill = ARM)) +
  geom_col(width = 0.6) +
  geom_errorbar(aes(ymin = Mean_Change - SE, ymax = Mean_Change + SE),
                width = 0.2, linewidth = 0.6) +
  
  # Numeric labels
  geom_text(aes(y = label_y, label = sprintf("%.2f", Mean_Change)),
            size = 2.5, fontface = "bold", color = "black",
            vjust = ifelse(summary_df$Mean_Change >= 0, 0, 1)) +
  
  
  facet_wrap(~Outcome, nrow = 3, scales = "fixed") +
  labs(
    title = "Average Change in Risk Score by Treatment Arm",
    x = "Treatment Arm",
    y = "Mean Change"
    
  ) +
  theme(
    axis.text.x = element_text(angle = 30, hjust = 1),
    strip.text  = element_text(face = "bold"),
    legend.position = "bottom",
    legend.text   = element_text(size = legend_text_size),
    legend.title  = element_text(size = legend_title_size),
    legend.key.size = unit(0.5, "lines"),
    plot.caption = element_text(hjust = 0.5, size = 9, margin = margin(t = 10))
  ) +
  scale_fill_manual(values = arm_colors) + 
  theme_bw() + 
  theme(axis.text.x = element_text(size = 8, angle = 45, hjust = 1))

# Step 3: Save
ggsave("mean_change_plot.pdf", plot = p, width = 12, height = 10)