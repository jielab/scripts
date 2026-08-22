#############################################################
####                Comparison WHR vs. WHtR              ####
#############################################################

# This script generates:
## Reviewer Fig 1: performances comparing OBSCORE with WHt-R vs. WH-R

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
cindex_test_list <- list()

## loop through each outcome
for (outc in outc_list) {
  # load models
  m.whr <-
    get(load(
      file = paste0(
        "<path>",
        outc,
        ".RData"
      )
    ))
  m.WHtR <-
    get(load(
      file = paste0(
        "<path>",
        outc,
        "WHR.RData"
      )
    ))
  
  # c
  cindex_data <- data.frame(
    Outcome = outc,
    FeatureSet = factor(c("WHtR", "WHR"),
                        levels = c("WHR", "WHtR")),
    Cindex.mean = c(m.whr$Cindex.mean,
                    m.WHtR$Cindex.mean),
    ci.low = c(m.whr$ci.low,
               m.WHtR$ci.low),
    ci.upp = c(m.whr$ci.upp,
               m.WHtR$ci.upp)
  )
  
  # --- bootstrap test: whr vs whtr
  delta <- m.whr$Cindex.vec - m.WHtR$Cindex.vec
  ci    <- quantile(delta, c(0.025, 0.975))
  mean_diff <- mean(delta)
  pval <- 2 * min(mean(delta <= 0), mean(delta >= 0))
  
  cindex_test_list[[outc]] <- data.frame(
    Outcome = outc,
    mean.delta = mean_diff,
    ci.low   = ci[1],
    ci.upp   = ci[2],
    pval     = pval
  )
  
  ## append to the list
  cindex_data_list[[outc]] <- cindex_data
  
}

## combine all C-index data
cindex_data_all <- as.data.table(bind_rows(cindex_data_list))
cindex_test_all <- as.data.table(bind_rows(cindex_test_list))

## define outcome labels
new_labels <- c("<nicenames>")

## replace outcome names
cindex_data_all$Outcome <- new_labels[cindex_data_all$Outcome]

pdf("<path>/WHR_WHtR.pdf", width = 7, height = 7)
cindex_data_all[FeatureSet == "WHR"] %>%
  left_join(., cindex_data_all[FeatureSet == "WHtR"], by = "Outcome") %>%
  ggplot(aes(Cindex.mean.y, Cindex.mean.x)) +
  geom_abline(
    slope = 1,
    intercept = 0,
    linetype = "dashed",
    color = "black"
  ) +
  geom_errorbarh(
    aes(xmin = ci.low.y, xmax = ci.upp.y, fill = Outcome),
    linewidth = 0.25,
    shape = 21,
    size = 0.35
  ) +
  geom_errorbar(
    aes(ymin = ci.low.x, ymax = ci.upp.x, fill = Outcome),
    linewidth = 0.25,
    shape = 21,
    size = 0.35
  ) +
  geom_point(
    shape = 21,
    fill = "skyblue",
    size = 2,
    alpha = .9
  ) +
  geom_text_repel(
    aes(label = Outcome),
    size = 2.5,
    hjust = -1.5,
    vjust = 1,
    max.overlaps = Inf,
    segment.size = 0.3,
    segment.color = "gray70"
  ) +
  labs(y = "C-index (95% CI) OBSCORE w/ waist-to-hip ratio", x = "C-index (95% CI) OBSCORE w/ waist-to-height ratio") +
  theme_light() +
  theme(
    panel.grid.minor = element_blank(),
    legend.title = element_blank(),
    legend.position = "none"
  )
dev.off()
