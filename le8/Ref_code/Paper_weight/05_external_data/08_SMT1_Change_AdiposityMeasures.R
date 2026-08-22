###################################################
####         Plot changes in adiposity         ####
###################################################

rm(list = ls())
options(stringsAsFactors = F)
options(scipen = 1)
setwd("<path>")

## --> packages needed <-- ##
require(ggplot2)
require(dplyr)
require(tidyr)

## read estimates in change of adiposity measures according to risk score
est.df <- read.csv("<path>/emmeans_cells.csv")
est.df <-
  est.df %>% filter(Risk_Score_Strata == "Type.2.Diabetes_group")

pval.df <- read.csv("<path>/interaction_vs_placebo.csv")
pval.df <-
  pval.df %>% filter(Endpoint %in% c("CHG_Weight (kg)",
                                     "PCHG_Weight (kg)",
                                     "WHtR_Change"))

pval.df <-
  pval.df %>% filter(Risk_Score_Strata == "Type.2.Diabetes_group")

## to rename endpoints
ylabs <- c(
  "PCHG_Weight (kg)" = "Change in weight (%)",
  "CHG_Weight (kg)"  = "Change in weight (kg)",
  "WHtR_Change"      = "Change in waist-height ratio"
)

## prep data with estimates
est_plot <- est.df %>%
  mutate(
    ARM = factor(ARM, levels = c("Placebo", "TZP 5mg", "TZP 10mg", "TZP 15mg")),
    Stratum = factor(Stratum),
    Endpoint = factor(Endpoint, levels = names(ylabs)),
    Endpoint = factor(
      Endpoint,
      levels = c("PCHG_Weight (kg)",
                 "CHG_Weight (kg)",
                 "WHtR_Change")
    )
  )

## prep data with pvalues
pval_annot <- est_plot %>%
  group_by(Endpoint, ARM) %>%
  summarise(
    x_pos  = mean(as.numeric(Stratum)),
    y_max  = max(emmean + SE),
    y_min  = min(emmean - SE),
    .groups = "drop"
  ) %>%
  filter(ARM != "Placebo") %>%
  left_join(pval.df %>% rename(ARM = Treatment_ARM),
            by = c("Endpoint", "ARM")) %>%
  mutate(
    y_pos   = y_max + 0.05 * (y_max - y_min),
    p_label = paste0("p interaction = ",
                     format(Interaction_P_Value, digits = 1)),
    Endpoint = factor(
      Endpoint,
      levels = c("PCHG_Weight (kg)",
                 "CHG_Weight (kg)",
                 "WHtR_Change")
    )
  )

## plot
pdf("<path>/ChangesAdiposityMeasures.pdf",
    width = 7,
    height = 2)
ggplot(est_plot, aes(
  x = Stratum,
  y = emmean,
  color = ARM,
  group = ARM
)) +
  facet_wrap(
    ~ Endpoint,
    scales = "free_y",
    labeller = labeller(Endpoint = as_labeller(ylabs)),
    ncol = 4
  ) +
  geom_line(linewidth = 0.7) +
  geom_errorbar(aes(ymin = emmean - SE, ymax = emmean + SE),
                width = 0.2) +
  geom_point(size = 1.5) +
  geom_text(
    data = pval_annot,
    aes(
      x = 1,
      y = y_pos,
      label = p_label,
      color = ARM
    ),
    size = 1.8,
    hjust = 0,
    inherit.aes = FALSE
  ) +
  scale_color_manual(
    values = c(
      "Placebo"  = "darkgrey",
      "TZP 5mg"  = "#a4cbce",
      "TZP 10mg" = "#6eb7b7",
      "TZP 15mg" = "#2d6274"
    )
  ) +
  labs(x = "Predicted T2D risk quartiles",
       y = "Mean change (± SE)") +
  theme_bw() +
  theme(
    legend.position = "right",
    legend.key.size = unit(0.4, "lines"),
    text = element_text(size = 8),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    strip.text = element_text(color = "white"),
    strip.background = element_rect(
      color = "black",
      fill = "black",
      linetype = "solid"
    ),
    panel.spacing = unit(0.2, "lines")
  )
dev.off()