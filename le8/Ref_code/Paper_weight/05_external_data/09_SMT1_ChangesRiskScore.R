###################################################
####            RCT changes in risk            ####
###################################################

rm(list = ls())
options(stringsAsFactors = F)
options(scipen = 1)
setwd("<path>")

## --> packages needed <-- ##
require(ggplot2)
require(dplyr)
require(tidyr)
require(forcats)

## read estimates in change of risk score
est.df <- read.csv("<path>/changeinscore_summary_stats.csv")
est.df$X <- NULL

## read pvalues for comparison between arms
diff.df <- read.csv("<path>/deltascore_pairwise_unadj.csv")
diff.df <- diff.df %>% filter(contrast == "Placebo - TZP 15mg")
diff.df$X <- NULL

## order outcomes by mean change in 15 mg arm
outcome_order <- est.df %>%
  filter(ARM == "TZP 15mg") %>%
  arrange(-desc(Mean_Change)) %>%
  pull(Outcome)

## factor levels
est_plot <- est.df %>%
  mutate(
    Outcome = factor(Outcome, levels = outcome_order),
    ARM = factor(ARM,
                 levels = c("Placebo", "TZP 5mg", "TZP 10mg", "TZP 15mg")),
    Mean_round = round(Mean_Change, 2)
  )

## get pvalues, assign positions for plotting
pvals.df <- est_plot %>%
  group_by(Outcome) %>%
  summarise(y_bracket = min(Mean_Change - SE) - 0.35) %>%
  left_join(diff.df %>% select(Outcome, p.value), by = "Outcome") %>%
  mutate(
    Outcome = factor(Outcome, levels = outcome_order),
    x_start = 1,
    x_end   = 4,
    x_mid   = (x_start + x_end) / 2,
    p_label = paste0("p = ", format(
      p.value, digits = 2, scientific = TRUE
    ))
  )

## plot
pdf("<path>/ChangesRiskScores.pdf",
    width = 7,
    height = 6)
ggplot(est_plot, aes(x = ARM, y = Mean_Change, fill = ARM)) +
  facet_wrap(~ Outcome, ncol = 6) +
  geom_col(width = 0.8) +
  geom_errorbar(aes(ymin = Mean_Change - SE, ymax = Mean_Change + SE),
                width = 0.1) +
  ## mean estimates
  geom_text(
    aes(y = Mean_Change - SE, label = Mean_round),
    vjust = 1.5,
    size = 1.8,
    color = "black"
  ) +
  ## pvalues
  geom_text(
    data = pvals.df,
    aes(x = x_mid, y = y_bracket, label = p_label),
    vjust = 1.3,
    size = 1.8,
    inherit.aes = FALSE
  ) +
  ## horizontal line for pvalue bracket
  geom_segment(
    data = pvals.df,
    aes(
      x = x_start,
      xend = x_end,
      y = y_bracket,
      yend = y_bracket
    ),
    inherit.aes = FALSE,
    size = .35,
    color = "grey"
  ) +
  ## left vertical line for pvalue bracket
  geom_segment(
    data = pvals.df,
    aes(
      x = x_start,
      xend = x_start,
      y = y_bracket,
      yend = y_bracket + 0.1
    ),
    inherit.aes = FALSE,
    size = .35,
    color = "grey"
  ) +
  ## right vertical line for pvalue bracket
  geom_segment(
    data = pvals.df,
    aes(
      x = x_end,
      xend = x_end,
      y = y_bracket,
      yend = y_bracket + 0.1
    ),
    inherit.aes = FALSE,
    size = .35,
    color = "grey"
  ) +
  scale_fill_manual(
    values = c(
      "Placebo"  =  "darkgrey",
      "TZP 5mg"  = "#a4cbce",
      "TZP 10mg" = "#6eb7b7",
      "TZP 15mg" = "#2d6274"
    )
  ) +
  labs(x = NULL,
       y = "Mean change after 72 week intervention (± SE)",
       fill = NULL) +
  ylim(-2.7, 0.1) +
  theme_bw() +
  theme(
    legend.position = "none",
    legend.key.size = unit(0.4, "lines"),
    text = element_text(size = 8),
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    strip.text = element_text(color = "white"),
    strip.background = element_rect(
      color = "black",
      fill = "black",
      linetype = "solid"
    ),
    panel.spacing = unit(0.1, "lines")
  )
dev.off()