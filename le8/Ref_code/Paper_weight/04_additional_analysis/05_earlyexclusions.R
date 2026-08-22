#############################################################
####              Exclusion of early events              ####
#############################################################
# This script generates:
## SFig 8: Performances excluding early events up to 2 years after baseline

rm(list = ls())
setwd("<path>")

## --> packages needed <-- ##
require(data.table)
require(tidyverse)
require(dplyr)
require(magrittr)
require(ggrepel)
require(patchwork)
require(tidytext)
require(ggsci)

## define outcomes
outc_list <- c("<outcomes>")

## list to store results
cindex_data_list <- list()

## loop through each outcome
for (outc in outc_list) {
  # load models
  m.0.5 <- 
    get(load(file = paste0("<path>", outc, ".RData")))
  m.1   <-
    get(load(file = paste0("<path>", outc, ".1year.excl.RData")))
  m.1.5 <-
    get(load(file = paste0("<path>", outc, ".1.5year.excl.RData")))
  m.2   <-
    get(load(file = paste0("<path>", outc, ".2year.excl.RData")))
  
  # data frame with performances
  cindex_data <- data.frame(
    Outcome = outc,
    FeatureSet = factor(
      c("0.5 years", "1 years", "1.5 years", "2 years"),
      levels = c("2 years", "1.5 years", "1 years", "0.5 years")
    ),
    Cindex.mean = c(
      m.0.5$Cindex.mean,
      m.1$Cindex.mean,
      m.1.5$Cindex.mean,
      m.2$Cindex.mean
    ),
    ci.low = c(m.0.5$ci.low,
               m.1$ci.low,
               m.1.5$ci.low,
               m.2$ci.low),
    ci.upp = c(m.0.5$ci.upp,
               m.1$ci.upp,
               m.1.5$ci.upp,
               m.2$ci.upp)
  )
  
  ## append to the list
  cindex_data_list[[outc]] <- cindex_data
  
}

## combine all C-index data
cindex_data_all <- as.data.table(bind_rows(cindex_data_list))

## define outcome labels
new_labels <- c("<nicenames>")

## replace outcome names
cindex_data_all$Outcome <- new_labels[cindex_data_all$Outcome]

## plot
cindex_data_all <- cindex_data_all %>%
  mutate(FeatureSet = factor(
    as.character(FeatureSet),
    levels = c("2 years", "1.5 years", "1 years", "0.5 years")
  ))

base <- cindex_data_all %>%
  filter(FeatureSet == "0.5 years") %>%
  select(Outcome, base_mean = Cindex.mean)

ann <- cindex_data_all %>%
  filter(FeatureSet %in% c("1 years", "1.5 years", "2 years")) %>%
  left_join(base, by = "Outcome") %>%
  mutate(
    delta = Cindex.mean - base_mean,
    label = number(delta, accuracy = 0.001, signed = TRUE),
    y = pmin(Cindex.mean + 0.1, 0.99)
  )

pdf("graphics/earlyexclusion_performances.pdf",
    width = 7,
    height = 6)
ggplot(cindex_data_all, aes(x = FeatureSet, y = Cindex.mean)) +
  coord_flip() +
  facet_wrap( ~ Outcome, ncol = 3) +
  geom_rect(
    aes(
      xmin = 3.5,
      xmax = 4.5,
      ymin = 0.4,
      ymax = 1,
      group = Outcome
    ),
    fill = "grey80",
    alpha = 0.25
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
    alpha = 0.25
  ) +
  geom_pointrange(
    aes(ymin = ci.low, ymax = ci.upp, fill = FeatureSet),
    linewidth = 0.25,
    shape = 21,
    size = 0.35
  ) +
  geom_label(
    data = ann,
    aes(x = FeatureSet, y = y, label = label),
    inherit.aes = FALSE,
    fill = NA,
    color = "black",
    size = 2,
    label.size = 0
  ) +
  scale_fill_nejm() +
  labs(x = "Exclusion of early cases", y = "C-index (95% CI)") +
  theme_light() + ylim(0.4, 1) +
  geom_hline(yintercept = 0.5,
             linetype = "dotted",
             color = "grey50") +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 0, hjust = 0.5),
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_blank(),
    strip.text = element_text(color = "white"),
    strip.background = element_rect(
      color = "black",
      fill = "black",
      linetype = "solid"
    ),
    text = element_text(size = 10)
  )
dev.off()
