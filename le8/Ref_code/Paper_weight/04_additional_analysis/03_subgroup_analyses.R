#############################################################
####                Subgroup performances                ####
#############################################################

# This script generates
## SFig 9: Subgroup analysis EUR vs. non-EUR and high and low Townsend index

rm(list = ls())

## --> packages needed <-- ##
require(data.table)
require(tidyverse)
require(dplyr)
require(magrittr)
require(ggrepel)
require(patchwork)
require(tidytext)
require(ggsci)

## set working directory
setwd("<path>")

## define outcomes
outc_list <- c("<outcomes")

#############################################
####    Subgroup performances (SFig 9)   ####
#############################################

## list to store results
cindex_data_list <- list()
cindex_test_list <- list()

## loop through each outcome
for (outc in outc_list) {
  ## load models
  m.EUR <-
    get(load(file = paste0("<path>", outc, ".EUR", ".RData")))
  m.nonEUR <-
    get(load(file = paste0("<path>", outc, ".nonEUR", ".RData")))
  m.dep <-
    get(load(file = paste0("<path>", outc, ".Deprived", ".RData")))
  m.aff <-
    get(load(file = paste0("<path>", outc, ".Affluent", ".RData")))
  
  ## data frame with performances
  cindex_data_list[[outc]] <- data.frame(
    Outcome = outc,
    FeatureSet = factor(
      c("European", "Non-European", "More deprived", "Less deprived"),
      levels = c("Less deprived", "More deprived", "Non-European", "European")
    ),
    Cindex.mean = c(
      m.EUR$Cindex.mean,
      m.nonEUR$Cindex.mean,
      m.dep$Cindex.mean,
      m.aff$Cindex.mean
    ),
    ci.low = c(m.EUR$ci.low,
               m.nonEUR$ci.low,
               m.dep$ci.low,
               m.aff$ci.low),
    ci.upp = c(m.EUR$ci.upp,
               m.nonEUR$ci.upp,
               m.dep$ci.upp,
               m.aff$ci.upp)
  )
  
  ## bootstrap test across subgroups
  delta.anc <- m.EUR$Cindex.vec - m.nonEUR$Cindex.vec
  ci.anc    <- quantile(delta.anc, c(0.025, 0.975))
  mean_diff.anc <- mean(delta.anc)
  pval.anc <- 2 * min(mean(delta.anc <= 0), mean(delta.anc >= 0))
  
  delta.dep <- m.dep$Cindex.vec - m.aff$Cindex.vec
  ci.dep    <- quantile(delta.dep, c(0.025, 0.975))
  mean_diff.dep <- mean(delta.dep)
  pval.dep <- 2 * min(mean(delta.dep <= 0), mean(delta.dep >= 0))
  
  ## data frame with tests
  cindex_test_list[[outc]] <- data.frame(
    Outcome = outc,
    mean.delta.anc = mean_diff.anc,
    ci.low.anc   = ci.anc[1],
    ci.upp.anc   = ci.anc[2],
    pval.anc     = pval.anc,
    mean.delta.dep = mean_diff.dep,
    ci.low.dep   = ci.dep[1],
    ci.upp.dep   = ci.dep[2],
    pval.dep     = pval.dep
  )
}

## combine all performances
cindex_data_all <- as.data.table(bind_rows(cindex_data_list))
## combine all comparisons
cindex_test_all <- as.data.table(bind_rows(cindex_test_list))

## define outcome labels
new_labels <- c("<nicenames>")

## replace outcome names
cindex_data_all$Outcome <- new_labels[cindex_data_all$Outcome]
cindex_test_all$Outcome <- new_labels[cindex_test_all$Outcome]

## reshape for plot
labels_long <- cindex_test_all %>%
  transmute(
    Outcome,
    FeatureSet_anc = "Non-European",
    label_anc = sprintf(
      "%.3f\n(%.3f to %.3f)\np = %s",
      mean.delta.anc,
      ci.low.anc,
      ci.upp.anc,
      round(pval.anc, digits = 2)
    )
  ) %>%
  rename(FeatureSet = FeatureSet_anc, label = label_anc) %>%
  bind_rows(
    cindex_test_all %>%
      transmute(
        Outcome,
        FeatureSet_dep = "More deprived",
        label_dep = sprintf(
          "%.3f\n(%.3f to %.3f)\np = %s",
          mean.delta.dep,
          ci.low.dep,
          ci.upp.dep,
          round(pval.dep, digits = 2)
        )
      ) %>%
      rename(FeatureSet = FeatureSet_dep, label = label_dep)
  )

labels_long <- labels_long %>%
  left_join(
    cindex_data_all %>% select(Outcome, FeatureSet, Cindex.mean, ci.upp),
    by = c("Outcome", "FeatureSet")
  )

## plot
pdf("<path>/subgroup_performances.pdf",
    width = 7,
    height = 7)
ggplot(cindex_data_all, aes(x = FeatureSet, y = Cindex.mean)) +
  coord_flip(clip = "off") +
  facet_wrap(~ Outcome, ncol = 3) +
  geom_rect(
    aes(
      xmin = 3.5,
      xmax = 4.5,
      ymin = 0.4,
      ymax = 1,
      group = Outcome
    ),
    fill = "grey80",
    alpha = 0.25,
    inherit.aes = FALSE
  ) +
  ## add bracket for p-value
  geom_segment(
    x = 3,
    xend = 4,
    y = 1.025,
    yend = 1.025,
    color = "grey50",
    size = .2
  ) +
  geom_segment(
    x = 4,
    xend = 4,
    y = 1.01,
    yend = 1.025,
    color = "grey50",
    size = .2
  ) +
  geom_segment(
    x = 3,
    xend = 3,
    y = 1.01,
    yend = 1.025,
    color = "grey50",
    size = .2
  ) +
  geom_segment(
    x = 1,
    xend = 2,
    y = 1.025,
    yend = 1.025,
    color = "grey50",
    size = .2
  ) +
  geom_segment(
    x = 2,
    xend = 2,
    y = 1.01,
    yend = 1.025,
    color = "grey50",
    size = .2
  ) +
  geom_segment(
    x = 1,
    xend = 1,
    y = 1.01,
    yend = 1.025,
    color = "grey50",
    size = .2
  ) +
  geom_rect(
    aes(
      xmin = 1.5,
      xmax = 2.5,
      ymin = 0.4,
      ymax = 1,
      group = Outcome
    ),
    fill = "grey80",
    alpha = 0.25,
    inherit.aes = FALSE
  ) +
  geom_hline(yintercept = 0.5,
             linetype = "dotted",
             color = "grey50") +
  geom_hline(
    yintercept = 0.4,
    linetype = "solid",
    color = "grey75",
    size = .2
  ) +
  geom_hline(
    yintercept = 0.6,
    linetype = "solid",
    color = "grey75",
    size = .2
  ) +
  geom_hline(
    yintercept = 0.8,
    linetype = "solid",
    color = "grey75",
    size = .2
  ) +
  geom_hline(
    yintercept = 1.0,
    linetype = "solid",
    color = "grey75",
    size = .2
  ) +
  geom_pointrange(
    aes(ymin = ci.low, ymax = ci.upp, fill = FeatureSet),
    linewidth = 0.25,
    shape = 21,
    size = 0.35
  ) +
  geom_text(
    data = labels_long[FeatureSet == "Non-European"],
    aes(x = 3.5, y = 1.15, label = label),
    hjust = 0.5,
    vjust = 0.5,
    size = 1.5,
    angle = 0
  ) +
  geom_text(
    data = labels_long[FeatureSet == "More deprived"],
    aes(x = 1.5, y = 1.15, label = label),
    hjust = 0.5,
    vjust = 0.5,
    size = 1.5,
    angle = 0
  ) +
  scale_fill_nejm() +
  labs(x = "Subgroups", y = "C-index (95% CI)") +
  theme_light() + ylim(0.4, 1.25) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 0, hjust = 0.5),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    strip.text = element_text(color = "white"),
    strip.background = element_rect(
      color = "black",
      fill = "black",
      linetype = "solid"
    ),
    text = element_text(size = 10),
    plot.margin = margin(5.5, 30, 5.5, 5.5)
  )
dev.off()
